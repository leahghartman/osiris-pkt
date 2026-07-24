# 1 "tiles/os-nconf-tiles.f03"
# 1 "<built-in>" 1
# 1 "<built-in>" 3
# 467 "<built-in>" 3
# 1 "<command line>" 1
# 1 "<built-in>" 2
# 1 "tiles/os-nconf-tiles.f03" 2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! node configuration class for tiles
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# 1 "./os-config.h" 1
! Configuration file for osiris

! ----------------------------------------------------------------------------------------
! Algorithm options
! ----------------------------------------------------------------------------------------

! ----------------------------------------------------------------------------------------
! System options
! ----------------------------------------------------------------------------------------

! MPI supports MPI_IN_PLACE in global operations
!#define 1

! Use OpenMP
!#define 1

! Use SIMD optimized code
!#define SIMD
!#define SIMD_SSE
!#define SIMD_AVX
!#define SIMD_BGQ
!#define SIMD_MIC

! Use SION for checkpointing (not available on all systems)
!#define __RST_IO__ = __RST_SION__

! Use log files
!#define __USE_LOG__

! Use MPE for logging and profiling
!#define __USE_MPE__

! Compiler does not fully support the sizeof intrinsic which is a Fortran 2003 feature
!#define __NO_SIZEOF__

! Use the PAPI library for profiling
!#define __USE_PAPI__


! ----------------------------------------------------------------------------------------
! Distribution Options
! ----------------------------------------------------------------------------------------
!
! Optional modules to be removed for distribution
! These are all turned on by default unless the __DISTRO__ preprocessor macro is defined
!
! ----------------------------------------------------------------------------------------



! Include Ionization module


! Include binary collisions module


! Include particle tracking module


! Include perfectly matched layers boundary conditions for EMF


! Include spin advance module which is disabled by default
! #define __HAS_SPIN__

! Include a debug flag for cylindrical modes simulations
! #define __CYL_MODES_DEBUG__

! Compile functions for reporting B fields in the RaDiO module
! #define __HAS_RAD_BFLD__
# 6 "tiles/os-nconf-tiles.f03" 2
# 1 "./os-preprocess.fpp" 1
!
! File: os-preprocess.fpp
!
! A set of preprocessing macros for osiris, and the fpp/cpp preprocessor
!
!



! Macros for the IBM CPP/GNU preprocessor




! Assertion macro



! Debug functions
# 33 "./os-preprocess.fpp"
! SCR hostname functions







! LOG functions




! ERROR functions






! WARNING functions




! functions for restart io
# 7 "tiles/os-nconf-tiles.f03" 2

module m_node_conf_tiles

# 1 "./memory/memory.h" 1
! Include file for the memory module
! Files must include this file using #include "memory/memory.h" rather than just using the module
! the #include must be placed where the module use statement would usually be i.e.
!
! module module2
!
! use module1
! #include "memory/memory.h"
!
# 19 "./memory/memory.h"
use memory
# 11 "tiles/os-nconf-tiles.f03" 2

use m_parameters
use m_system
use m_restart
use m_node_conf, only : t_node_conf, p_topology_mpi
use m_tile_conf
use m_input_file, only : t_input_file, get_namelist
use m_hilbert
use m_snake

implicit none

! restrict access to things explicitly declared public
private

! string to id restart data
character(len=*), parameter :: p_nconf_tiles_rst_id = "node conf tiles rst data - 0x0001"

! forcing heavy/light tiles
integer, parameter :: p_default = 0
integer, parameter :: p_heavy = 1
integer, parameter :: p_light = 2

! minimal block size to be used when growing particle buffers in a tile
! units are particle number
! moved here instead of os-spec-comm-tiles to be set in the input deck
integer, public :: p_spec_buf_block_t = 1024

type :: t_omp_patt

  ! force heavy, light, or default (calculate load)
  integer :: force_omp

  ! number of heavy/light tiles
  integer :: n_heavy, n_light

  ! arrays containing tile numbers of heavy/light tiles
  integer, dimension(:), pointer :: heavy_tils => null()
  integer, dimension(:), pointer :: light_tils => null()

end type

type, extends(t_node_conf) :: t_node_conf_tiles

  ! tile configuration class
  type( t_tile_conf ) :: ti_co

  logical :: group

  type( t_omp_patt ) :: omp_patt

contains

  procedure :: init => init_node_conf_tiles
  procedure :: read_input => read_input_node_conf_tiles
  procedure :: cleanup => cleanup_node_conf_tiles
  procedure :: write_checkpoint => write_checkpoint_node_conf_tiles
  procedure :: restart_read => restart_read_node_conf_tiles

  procedure :: setup_tiles => setup_node_conf_tiles
  procedure :: on_edge => on_edge_tiles
  procedure :: physical_boundary_local => physical_boundary_local_tiles
  procedure :: physical_boundary_node => physical_boundary_node_tiles
  procedure :: ngp => ngp_tiles
  procedure :: ngp_id => ngp_id_tiles
  procedure :: ifpr => relevant_ifpr_tiles

end type t_node_conf_tiles

interface node_of_til
  module procedure node_of_til
end interface node_of_til

interface node_of_til_old
  module procedure node_of_til_old
end interface node_of_til_old

interface set_ti_co_viktor
  module procedure set_ti_co_viktor
end interface

public :: p_default, p_heavy, p_light
public :: node_of_til, node_of_til_old
public :: set_ti_co_viktor, t_node_conf_tiles
public :: t_omp_patt

contains

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine read_input_node_conf_tiles( this, input_file, x_dim )

  implicit none

  class( t_node_conf_tiles ), intent(inout) :: this
  class( t_input_file ), intent(inout) :: input_file
  integer, intent(in) :: x_dim

  integer, dimension(p_x_dim) :: node_number
  logical, dimension(p_x_dim) :: if_periodic
  integer :: n_threads, grow_spec_buff

  character(len=16) :: topology, tile_topology, force_omp

  integer, dimension(p_x_dim) :: tile_number

  namelist /nl_node_conf/ node_number, if_periodic, n_threads, topology, tile_number, &
                          tile_topology, force_omp, grow_spec_buff

  integer :: i, ierr

  this%x_dim = x_dim

  node_number = 1
  if_periodic = .false.
  n_threads = 1
  tile_number = 0
  grow_spec_buff = p_spec_buf_block_t

  ! Default topology routines
  topology = "mpi"
  tile_topology = "snake"
  force_omp = "default"

  ! Get namelist text from input file
  call get_namelist( input_file, "nl_node_conf", ierr )

  if (ierr /= 0) then
    if ( mpi_node() == 0 ) then
      if (ierr < 0) then
        write(0,*) "Error reading node_conf parameters"
      else
        write(0,*) "Error: node_conf parameters missing"
      endif
      write(0,*) "aborting..."
    endif
    stop
  endif

  read (input_file%nml_text, nml = nl_node_conf, iostat = ierr)
  if (ierr /= 0) then
    if ( mpi_node() == 0 ) then
      write(0,*) "Error reading node_conf parameters"
      write(0,*) "aborting..."
    endif
    stop
  endif

  this%no_num = 1
  this%ti_co%g_til_num = 1
  do i=1, x_dim
    this%nx(i) = node_number(i)
    this%ifpr_(i) = if_periodic(i)
    this%no_num = this%no_num * this%nx(i)
    this%ti_co%g_til(i) = tile_number(i)
    this%ti_co%g_til_num = this%ti_co%g_til_num * this%ti_co%g_til(i)
  enddo

  ! process node topology type
  select case ( trim( topology ))
  case ( "mpi" )
    this%topology_type = p_topology_mpi
  case ( "old" )
    if ( mpi_node() == 0 ) then
      write(0,*) "(*warning*) Defaulting to MPI topology with tiling"
    endif
    this%topology_type = p_topology_mpi
# 187 "tiles/os-nconf-tiles.f03"
  case default
    if ( mpi_node() == 0 ) then
      write(0,*) 'Error reading node_conf parameters'
      write(0,*) 'invalid topology type: "', trim( topology ), '"'
      write(0,*) 'With tiling, only valid value is "mpi".'
      write(0,*) 'aborting...'
    endif
    stop
  end select

  ! process tile topology type
  select case ( trim( tile_topology ))
  case ( "hilbert" )
    this%ti_co%topology_type = p_topology_hilbert

  case ( "snake" )
    this%ti_co%topology_type = p_topology_snake

  case ( "rectangular", "viktor" )
    this%ti_co%topology_type = p_topology_viktor

  case default
    if ( mpi_node() == 0 ) then
      write(0,*) 'Error reading node_conf parameters'
      write(0,*) 'invalid tile_topology type: "', trim( tile_topology ), '"'
      write(0,*) 'Valid values are "hilbert", "snake" and "rectangular".'
      write(0,*) 'aborting...'
    endif
    stop
  end select

  ! process forcing heavy/light tiles
  select case( trim(force_omp) )
  case( "default" )
    this%omp_patt%force_omp = p_default
  case( "heavy" )
    this%omp_patt%force_omp = p_heavy
  case( "light" )
    this%omp_patt%force_omp = p_light
  case default
    if ( mpi_node() == 0 ) then
      write(0,*) 'Error reading node_conf parameters'
      write(0,*) 'invalid force_omp type: "', trim( force_omp ), '"'
      write(0,*) 'Valid values are "default", "heavy" and "light".'
      write(0,*) 'aborting...'
    endif
    stop
  end select

  if ( this%no_num < 1 ) then
    if ( mpi_node() == 0 ) then
      write(0,*) "Error reading node_conf parameters"
      write(0,*) 'Number of nodes (node_number) must be > 0 in all directions'
      write(0,*) "aborting..."
    endif
    stop
  endif

  ! set how many particles to grow species buffer by in update boundary/dlb routines
  if ( grow_spec_buff > 0 ) then
    p_spec_buf_block_t = grow_spec_buff
  else
    if ( mpi_node() == 0 ) then
      write(0,*) "Error reading node_conf parameters"
      write(0,*) "grow_spec_buff must be > 0"
      write(0,*) "aborting..."
    endif
    stop
  endif


  ! Multithread support
  this%n_threads = n_threads

  if ( n_threads < 1 ) then
    if ( mpi_node() == 0 ) then
      write(0,*) "Error reading node_conf parameters"
      write(0,*) "Invalid number of threads per node (n_threads)."
      write(0,*) "aborting..."
    endif
    stop
  endif
# 281 "tiles/os-nconf-tiles.f03"
  if ( this%ti_co%g_til_num < 1 ) then
    if ( mpi_node() == 0 ) then
      write(0,*) "Error reading node_conf parameters"
      write(0,*) "Number of tiles (tile_number) must be > 0 in all directions."
      write(0,*) "aborting..."
    endif
    stop
  endif

  ! This is necessary for diagnostics
  do i=1, x_dim
    if ( this%ti_co%g_til(i) < this%nx(i) ) then
      if ( mpi_node() == 0 ) then
        write(0,*) "Error reading node_conf parameters"
        write(0,*) "Number of tiles (tile_number) must be > number of nodes (node_number) in each dimension"
        write(0,*) "aborting..."
      endif
      stop
    endif
  enddo

end subroutine read_input_node_conf_tiles
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! sets up this data structure from the given information
!-----------------------------------------------------------------------------------------
subroutine init_node_conf_tiles( this )


  use omp_lib


  implicit none

  class( t_node_conf_tiles ), intent( inout ) :: this

  integer :: started, id, i, ierr, n_threads

  ! Belongs to sim_group
  this % group = .true.

  ! Initialize MPI if required
  if ( this%ti_co%g_til_num > 1 ) then

    call MPI_COMM_SIZE( mpi_comm_world, started, ierr )
    if ( ierr /= MPI_SUCCESS ) then
      write(err_buf__,*) "Unable to get MPI size.";call err__("tiles/os-nconf-tiles.f03",328)
      call abort_program( p_err_mpi )
    endif

    if ( started /= this%no_num ) then
      write(0,'(A,I0,A)') "(*error*) ", started, " nodes started up"
      write(0,'(A,I0,A)') "(*error*) ", this%no_num, " nodes should have started"
      call abort_program( p_err_invalid )
    endif

    call MPI_COMM_RANK( mpi_comm_world, id, ierr )
    if ( ierr /= MPI_SUCCESS ) then
      write(err_buf__,*) "Unable to get MPI rank.";call err__("tiles/os-nconf-tiles.f03",340)
      call abort_program( p_err_mpi )
    endif

    ! store using fortran ordering
    this%my_aid_ = id + 1

  else
    this%my_aid_ = 1
  endif




  ! Initilize OpenMP if required
  if ( this%n_threads > 1 ) then

    if ( this%n_threads > omp_get_max_threads() ) then
      write(0,*) '(*error*) The number of threads requested, ', this%n_threads , ' is more than the'
      write(0,*) '(*error*) maximum number of threads available, ', omp_get_max_threads()
      call abort_program()
    endif

  endif

  ! Set the number of threads to use
  call omp_set_num_threads( this%n_threads )

  ! Check the operation was successfull (on some systems you need to set the number of threads
  ! externally e.g. in the job script)

  !$omp parallel shared(n_threads)
  n_threads = omp_get_num_threads()
  !$omp end parallel

  if ( this%n_threads /= n_threads ) then
    write(0,*) '(*error*) Unable to set the number of threads requested'
    write(0,*) '(*error*) You may need to set the number of threads outside OSIRIS e.g. ', &
                          'in the job script'
    call abort_program()
  endif

  ! Disable dynamic variation of number of threads
  call omp_set_dynamic (.false.)

  ! Disable nested parallelism
  ! call omp_set_nested(.false.) ! No longer supported in OpenMP 6
  call omp_set_max_active_levels(1)

  ! Print parallel information
  if ( this%my_aid_ == 1 ) then
    print '(A,I0,A,I0,A)', "  - Running on ", this%no_num, " processes, ", this%n_threads, &
               " thread(s) per process"
  endif
# 404 "tiles/os-nconf-tiles.f03"
  ! setup comm, ngp and neighbor data
  if (mpi_node()==0) print *," setting up tile topology"
  if ( this%ti_co%g_til_num > 1 ) then

    select case ( this%ti_co%topology_type )
    case ( p_topology_hilbert )
      if ( this%my_aid_ == 1 ) then
        print '(A)', "  - Using MPI_Cart_create based topology"
      endif
      call create_topology_hilbert( this )
    case ( p_topology_snake )
      if ( this%my_aid_ == 1 ) then
        print '(A)', "  - Using MPI_Cart_create based topology"
      endif
      call create_topology_snake( this )
    case ( p_topology_viktor )
      if ( this%my_aid_ == 1 ) then
        print '(A)', "  - Using MPI_Cart_create based topology"
      endif
      call create_topology_viktor( this )
    end select

  else

    ! Serial run, no need to setup the topology
    this%my_ngp = 1
    do i = 1, this%x_dim
      if ( this%ifpr_(i) ) then
        this%neighbor(:,i) = 1
      else
        this%neighbor(:,i) = -1
      endif
    enddo

    ! Need to initialize the tile list parameters
    select case ( this%ti_co%topology_type )
    case ( p_topology_hilbert )
      call create_topology_hilbert( this )
    case ( p_topology_snake )
      call create_topology_snake( this )
    case ( p_topology_viktor )
      call create_topology_viktor( this )
    end select

  endif

end subroutine init_node_conf_tiles
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine cleanup_node_conf_tiles( this, barrier )
!-----------------------------------------------------------------------------------------
! Object destructor, cleans up all dynamically
! allocated memory
!-----------------------------------------------------------------------------------------

  implicit none

  class( t_node_conf_tiles ), intent( inout ) :: this
  logical, intent(in), optional :: barrier

  ! Synchronize nodes
  if (present(barrier)) then
    if (barrier) then
      call this % barrier()
    endif
  else
    call this % barrier()
  endif

  ! Only cleanup if we're on sim_group
  if (this%free_comm) call this%ti_co%cleanup()

  ! Cleanup omp_patt
  if ( associated(this%omp_patt%heavy_tils) ) call freemem(this%omp_patt%heavy_tils,"tiles/os-nconf-tiles.f03",478)
  if ( associated(this%omp_patt%light_tils) ) call freemem(this%omp_patt%light_tils,"tiles/os-nconf-tiles.f03",479)

  ! Call superclass cleanup
  if (present(barrier)) then
    call this % t_node_conf % cleanup(barrier=barrier)
  else
    call this % t_node_conf % cleanup()
  endif

end subroutine cleanup_node_conf_tiles
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! write object information into a restart file
!-----------------------------------------------------------------------------------------
subroutine write_checkpoint_node_conf_tiles( this, restart_handle )

  implicit none

  class( t_node_conf_tiles ), intent(in) :: this
  type( t_restart_handle ), intent(inout) :: restart_handle

  character(len=*), parameter :: err_msg = 'error writing restart data for node_conf_tiles object.'
  integer :: i, ierr

  call restart_io_write("p_nconf_tiles_rst_id", p_nconf_tiles_rst_id, restart_handle, ierr)
  call check_error(ierr,err_msg,p_err_rstwrt,"tiles/os-nconf-tiles.f03",505)


  ! this data is not needed for restart, this is just used to insure
  ! that the input deck has not changed
  call restart_io_write("this%nx", this%nx, restart_handle, ierr)
  call check_error(ierr,err_msg,p_err_rstwrt,"tiles/os-nconf-tiles.f03",511)

  call restart_io_write("this%ifpr_", this%ifpr_, restart_handle, ierr)
  call check_error(ierr,err_msg,p_err_rstwrt,"tiles/os-nconf-tiles.f03",514)

  call restart_io_write("this%my_ngp", this%my_ngp, restart_handle, ierr)
  call check_error(ierr,err_msg,p_err_rstwrt,"tiles/os-nconf-tiles.f03",517)

  call restart_io_write("this%ti_co%g_til", this%ti_co%g_til, restart_handle, ierr)
  call check_error(ierr,err_msg,p_err_rstwrt,"tiles/os-nconf-tiles.f03",520)

  call restart_io_write("this%ti_co%topology_type", this%ti_co%topology_type, restart_handle, ierr)
  call check_error(ierr,err_msg,p_err_rstwrt,"tiles/os-nconf-tiles.f03",523)

  ! Write other things that will be copied over
  call restart_io_write("this%ti_co%h", this%ti_co%h, restart_handle, ierr)
  call check_error(ierr,err_msg,p_err_rstwrt,"tiles/os-nconf-tiles.f03",527)

  select case( this%x_dim )
  case(1)
    call restart_io_write("this%ti_co%hi1", this%ti_co%hi1, restart_handle, ierr)
    call check_error(ierr,err_msg,p_err_rstwrt,"tiles/os-nconf-tiles.f03",532)
  case(2)
    call restart_io_write("this%ti_co%hi2", this%ti_co%hi2, restart_handle, ierr)
    call check_error(ierr,err_msg,p_err_rstwrt,"tiles/os-nconf-tiles.f03",535)
  case(3)
    call restart_io_write("this%ti_co%hi3", this%ti_co%hi3, restart_handle, ierr)
    call check_error(ierr,err_msg,p_err_rstwrt,"tiles/os-nconf-tiles.f03",538)
  end select

  call restart_io_write("this%ti_co%g_til_aid_min_max", this%ti_co%g_til_aid_min_max, restart_handle, ierr)
  call check_error(ierr,err_msg,p_err_rstwrt,"tiles/os-nconf-tiles.f03",542)

  if ( this%ti_co%topology_type == p_topology_viktor ) then

    do i = 1, this%x_dim
      call restart_io_write("this%ti_co%cart_til_mm(i)%ii", this%ti_co%cart_til_mm(i)%ii, restart_handle, ierr)
      call check_error(ierr,err_msg,p_err_rstwrt,"tiles/os-nconf-tiles.f03",548)
    enddo

  endif

end subroutine write_checkpoint_node_conf_tiles
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! read object information from a restart file
!-----------------------------------------------------------------------------------------
subroutine restart_read_node_conf_tiles( this, restart_handle )

  implicit none

  class( t_node_conf_tiles ), intent(inout) :: this
  type( t_restart_handle ), intent(in) :: restart_handle

  character(len=len(p_nconf_tiles_rst_id)) :: rst_id
  character(len=*), parameter :: err_msg = 'error reading restart data for node_conf_tiles object.'
  integer, dimension(p_max_dim) :: nx, my_ngp
  logical, dimension(p_max_dim) :: ifpr
  integer, dimension(p_x_dim) :: g_til
  integer :: tile_topology, i, ierr

  call restart_io_read(rst_id, restart_handle, ierr)
  call check_error(ierr,err_msg,p_err_rstrd,"tiles/os-nconf-tiles.f03",574)

  ! check if restart file is compatible
  if ( rst_id /= p_nconf_tiles_rst_id ) then
    write(err_buf__,*) 'Corrupted restart file, or restart file ';call err__("tiles/os-nconf-tiles.f03",578)
    write(err_buf__,*) 'from incompatible binary (node_conf_tiles)';call err__("tiles/os-nconf-tiles.f03",579)
    write(err_buf__,*) 'rst_id = ', rst_id;call err__("tiles/os-nconf-tiles.f03",580)
    call abort_program(p_err_rstrd)
  endif

  call restart_io_read(nx, restart_handle, ierr)
  call check_error(ierr,err_msg,p_err_rstrd,"tiles/os-nconf-tiles.f03",585)

  call restart_io_read(ifpr, restart_handle, ierr)
  call check_error(ierr,err_msg,p_err_rstrd,"tiles/os-nconf-tiles.f03",588)

  call restart_io_read(my_ngp, restart_handle, ierr)
  call check_error(ierr,err_msg,p_err_rstrd,"tiles/os-nconf-tiles.f03",591)

  call restart_io_read(g_til, restart_handle, ierr)
  call check_error(ierr,err_msg,p_err_rstrd,"tiles/os-nconf-tiles.f03",594)

  call restart_io_read(tile_topology, restart_handle, ierr)
  call check_error(ierr,err_msg,p_err_rstrd,"tiles/os-nconf-tiles.f03",597)

  ! check if the user messed up the input file
  do i = 1, p_x_dim
    if ( nx(i) /= this%nx(i) ) then
      write(err_buf__,*) 'The node_number specified on the input deck does';call err__("tiles/os-nconf-tiles.f03",602)
      write(err_buf__,*) 'not match the restart data';call err__("tiles/os-nconf-tiles.f03",603)
      call abort_program(p_err_rstrd)
    endif

    if ( ifpr(i) .neqv. this%ifpr_(i) ) then
      write(err_buf__,*) 'The if_periodic specified on the input deck does';call err__("tiles/os-nconf-tiles.f03",608)
      write(err_buf__,*) 'not match the restart data';call err__("tiles/os-nconf-tiles.f03",609)
      call abort_program(p_err_rstrd)
    endif

    if ( my_ngp(i) /= this%my_ngp(i) ) then
      write(err_buf__,*) 'NGP mismatch, reading data for different node';call err__("tiles/os-nconf-tiles.f03",614)
      call abort_program(p_err_rstrd)
    endif

    if ( g_til(i) /= this%ti_co%g_til(i) ) then
      write(err_buf__,*) 'The tile_number specified on the input deck does';call err__("tiles/os-nconf-tiles.f03",619)
      write(err_buf__,*) 'not match the restart data';call err__("tiles/os-nconf-tiles.f03",620)
      call abort_program(p_err_rstrd)
    endif

  enddo

  if ( tile_topology /= this%ti_co%topology_type ) then
    write(err_buf__,*) 'The tile_topology specified on the input deck does';call err__("tiles/os-nconf-tiles.f03",627)
    write(err_buf__,*) 'not match the restart data';call err__("tiles/os-nconf-tiles.f03",628)
    call abort_program(p_err_rstrd)
  endif

  ! Read other things that will be copied over
  call restart_io_read(this%ti_co%h, restart_handle, ierr)
  call check_error(ierr,err_msg,p_err_rstrd,"tiles/os-nconf-tiles.f03",634)

  select case( this%x_dim )
  case(1)
    call restart_io_read(this%ti_co%hi1, restart_handle, ierr)
    call check_error(ierr,err_msg,p_err_rstrd,"tiles/os-nconf-tiles.f03",639)
  case(2)
    call restart_io_read(this%ti_co%hi2, restart_handle, ierr)
    call check_error(ierr,err_msg,p_err_rstrd,"tiles/os-nconf-tiles.f03",642)
  case(3)
    call restart_io_read(this%ti_co%hi3, restart_handle, ierr)
    call check_error(ierr,err_msg,p_err_rstrd,"tiles/os-nconf-tiles.f03",645)
  end select

  call restart_io_read(this%ti_co%g_til_aid_min_max, restart_handle, ierr)
  call check_error(ierr,err_msg,p_err_rstrd,"tiles/os-nconf-tiles.f03",649)

  if ( this%ti_co%topology_type == p_topology_viktor ) then

    do i = 1, this%x_dim
      call restart_io_read(this%ti_co%cart_til_mm(i)%ii, restart_handle, ierr)
      call check_error(ierr,err_msg,p_err_rstrd,"tiles/os-nconf-tiles.f03",655)
    enddo

  endif

end subroutine restart_read_node_conf_tiles
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Returns true if tile boundary is on the edge of the simulation box
!-----------------------------------------------------------------------------------------
function on_edge_tiles( this, dim, bnd )

  implicit none

  class( t_node_conf_tiles ), intent(in) :: this
  integer, intent(in) :: dim, bnd

  logical :: on_edge_tiles

  if ( bnd == p_lower ) then
    on_edge_tiles = ( this%ti_co%g_tgp(dim) == 1 )
  else
    on_edge_tiles = ( this%ti_co%g_tgp(dim) == this%ti_co%g_til(dim) )
  endif

end function on_edge_tiles
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Returns true if local node boundary is a physical boundary
! (non communication / periodic )
!-----------------------------------------------------------------------------------------
function physical_boundary_local_tiles( this, dim, bnd )

  implicit none

  class( t_node_conf_tiles ), intent(in) :: this
  integer, intent(in) :: dim, bnd

  logical :: physical_boundary_local_tiles

  if ( this%ifpr_( dim ) ) then
    ! Periodic boundaries are not physical boundaries
    physical_boundary_local_tiles = .false.
  else
    ! To be a physical boundary the boundary must be on the edge of the simulation box
    if ( bnd == p_lower ) then
      physical_boundary_local_tiles = ( this%ti_co%g_tgp(dim) == 1 )
    else
      physical_boundary_local_tiles = ( this%ti_co%g_tgp(dim) == this%ti_co%g_til(dim) )
    endif
  endif

end function physical_boundary_local_tiles
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Returns true if node boundary is a physical boundary (non communication / periodic )
!-----------------------------------------------------------------------------------------
function physical_boundary_node_tiles( this, node, dim, bnd )

  implicit none

  class( t_node_conf_tiles ), intent(in) :: this
  integer, intent(in) :: node, dim, bnd

  logical :: physical_boundary_node_tiles

  if ( this%ifpr_( dim ) ) then
    ! Periodic boundaries are not physical boundaries
    physical_boundary_node_tiles = .false.
  else
    print *, "physical_boundary_node_tiles not currently supported in tiling"
    print *, "aborting..."
    stop
    ! To be a physical boundary the boundary must be on the edge of the simulation box
    if ( bnd == p_lower ) then
      physical_boundary_node_tiles = ( this%ngp( node, dim ) == 1 )
    else
      physical_boundary_node_tiles = ( this%ngp( node, dim ) == this%nx(dim) )
    endif
  endif

end function physical_boundary_node_tiles
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Returns the position of the requested node on the parallel partition
!-----------------------------------------------------------------------------------------
function ngp_tiles( this, node, dim )

  implicit none

  class( t_node_conf_tiles ), intent( in ) :: this
  integer, intent( in ) :: node
  integer, intent(in) :: dim

  integer :: ngp_tiles

  integer, dimension( p_max_dim ) :: coords
  integer :: ierr

  if ( this%no_num > 1 ) then
    select case ( this%ti_co%topology_type )
    case ( p_topology_hilbert, p_topology_snake, p_topology_viktor )
      call mpi_cart_coords( this%comm, node - 1, this%x_dim, coords, ierr )
      if ( ierr /= MPI_SUCCESS ) then
        write(err_buf__,*) "mpi_cart_coords failed";call err__("tiles/os-nconf-tiles.f03",763)
        call abort_program( p_err_mpi )
      endif
      ngp_tiles = coords(dim) + 1
    case default
      write(err_buf__,*) "unknown tile topology";call err__("tiles/os-nconf-tiles.f03",768)
      call abort_program( p_err_invalid )
      ngp_tiles = -1
    end select
  else
    ngp_tiles = 1
  endif


end function ngp_tiles
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Returns unique id integer based on the tile grid position instead of node grid position
! Calls the parent class if this belongs to sim_group, else returns essentially tgp_id
! We use this rather than g_til_aid since it should remain constant even after dlb
! The special case should only be used when assigning tags for initialized particles
!-----------------------------------------------------------------------------------------
pure function ngp_id_tiles(this)

  implicit none

  class( t_node_conf_tiles ), intent(in) :: this
  integer :: i,st,ngp_id_tiles

  if (this%group) then

    ngp_id_tiles = this % t_node_conf % ngp_id()

  else

    st = 1
    ngp_id_tiles = 1

    do i = 1, this%x_dim
      ngp_id_tiles = ngp_id_tiles + ( this%ti_co%g_tgp(i)-1 ) * st
      st = st * this%ti_co%g_til(i)
    enddo

  endif

end function ngp_id_tiles
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
function relevant_ifpr_tiles( this )
!-----------------------------------------------------------------------------------------
! gives flags for the periodic node configuration WHERE IT IS RELEVANT for the rest of
! the code. That means in directions which are set for periodicity and that have only one
! node if there is more then one node in that directions the rest of the code does not
! notice the periodicity in that direction since it just changes which nodes data are
! shipped to
!-----------------------------------------------------------------------------------------

  implicit none

  class( t_node_conf_tiles ), intent(in) :: this
  logical, dimension(this%x_dim) :: relevant_ifpr_tiles

  integer :: i

  do i=1, this%x_dim
    if ( this%ifpr_(i) .and. ( this%ti_co%g_til(i) .eq. 1 ) ) then
      relevant_ifpr_tiles(i) = .true.
    else
      relevant_ifpr_tiles(i) = .false.
    endif
  enddo

end function relevant_ifpr_tiles
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine setup_node_conf_tiles( n_a, n_b )
!-----------------------------------------------------------------------------------------
! sets up this data structure from the given information
!-----------------------------------------------------------------------------------------

  implicit none

  class( t_node_conf_tiles ), intent(inout) :: n_a
  class( t_node_conf ), intent(in) :: n_b

  ! Does not belong to sim_group
  n_a%group = .false.

  select type(n_b); class is(t_node_conf_tiles)

  n_a%x_dim = n_b%x_dim

  n_a%comm = n_b%comm

  n_a%free_comm = n_b%free_comm

  n_a%no_num = n_b%no_num

  n_a%nx = n_b%nx

  n_a%my_aid_ = n_b%my_aid_

  n_a%my_ngp = n_b%my_ngp

  n_a%neighbor = n_b%neighbor

  n_a%ifpr_ = n_b%ifpr_

  n_a%n_threads = n_b%n_threads

  n_a%topology_type = n_b%topology_type

  n_a%ti_co%g_til = n_b%ti_co%g_til

  n_a%ti_co%g_til_num = n_b%ti_co%g_til_num

  n_a%ti_co%h => n_b%ti_co%h

  n_a%ti_co%g_til_aid_min_max => n_b%ti_co%g_til_aid_min_max

  n_a%ti_co%g_til_aid_min_max_old => n_b%ti_co%g_til_aid_min_max_old

  n_a%ti_co%topology_type = n_b%ti_co%topology_type

  select case( n_a%x_dim )
  case(1)
    n_a%ti_co%hi1 => n_b%ti_co%hi1
  case(2)
    n_a%ti_co%hi2 => n_b%ti_co%hi2
  case(3)
    n_a%ti_co%hi3 => n_b%ti_co%hi3
  end select

  call setup_tile_topology( n_a, n_b )

  end select

end subroutine setup_node_conf_tiles
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Create simulation topology for tiles
! Nodes will be distributed by x3 first, then x2, then x1.
! This is the same ordering as mpi_cart_create
!-----------------------------------------------------------------------------------------
subroutine setup_tile_topology( n_a, n_b )
!-----------------------------------------------------------------------------------------
  implicit none

  class(t_node_conf_tiles), intent(inout) :: n_a
  class(t_node_conf_tiles), intent(in) :: n_b

  integer :: i
  integer, dimension(p_max_dim) :: l, u, ind

  n_a%ti_co%g_tgp(1:n_a%x_dim) = n_b%ti_co%h(:,n_a%ti_co%g_til_aid)

  do i = 1, n_a%x_dim

    ! lower tile id
    ind = n_a%ti_co%g_tgp
    if (ind(i)==1) then
      if (n_a%ifpr_(i)) then
        ind(i) = n_a%ti_co%g_til(i)
        select case( n_a%x_dim )
        case(1)
          l(i) = n_b%ti_co%hi1( ind(1) )
        case(2)
          l(i) = n_b%ti_co%hi2( ind(1), ind(2) )
        case(3)
          l(i) = n_b%ti_co%hi3( ind(1), ind(2), ind(3) )
        end select
      else
        l(i) = -1
      endif
    else
      ind(i) = ind(i) - 1
      select case( n_a%x_dim )
      case(1)
        l(i) = n_b%ti_co%hi1( ind(1) )
      case(2)
        l(i) = n_b%ti_co%hi2( ind(1), ind(2) )
      case(3)
        l(i) = n_b%ti_co%hi3( ind(1), ind(2), ind(3) )
      end select
    endif

    ! lower node id
    if (l(i)==-1) then
      n_a%neighbor(p_lower,i) = -1
    else
      n_a%neighbor(p_lower,i) = node_of_til( n_a, l(i), n_b )
    endif

    ! upper tile id
    ind = n_a%ti_co%g_tgp
    if (ind(i)==n_a%ti_co%g_til(i)) then
      if (n_a%ifpr_(i)) then
        ind(i) = 1
        select case( n_a%x_dim )
        case(1)
          u(i) = n_b%ti_co%hi1( ind(1) )
        case(2)
          u(i) = n_b%ti_co%hi2( ind(1), ind(2) )
        case(3)
          u(i) = n_b%ti_co%hi3( ind(1), ind(2), ind(3) )
        end select
      else
        u(i) = -1
      endif
    else
      ind(i) = ind(i) + 1
      select case( n_a%x_dim )
      case(1)
        u(i) = n_b%ti_co%hi1( ind(1) )
      case(2)
        u(i) = n_b%ti_co%hi2( ind(1), ind(2) )
      case(3)
        u(i) = n_b%ti_co%hi3( ind(1), ind(2), ind(3) )
      end select
    endif

    ! upper node id
    if (u(i)==-1) then
      n_a%neighbor(p_upper,i) = -1
    else
      n_a%neighbor(p_upper,i) = node_of_til( n_a, u(i), n_b )
    endif

  enddo

  n_a%ti_co%neighbor_g_til(p_lower,1:n_a%x_dim) = l(1:n_a%x_dim)
  n_a%ti_co%neighbor_g_til(p_upper,1:n_a%x_dim) = u(1:n_a%x_dim)

  ! print *, "[",mpi_node(),"] ",n_a%neighbor(:,1:n_a%x_dim)

  ! (* debug *) Test topology
  ! call test_topology( n_a )

end subroutine setup_tile_topology

!-----------------------------------------------------------------------------------------
! Create simulation topology using hilbert routines
!-----------------------------------------------------------------------------------------
subroutine create_topology_hilbert( this )

  implicit none

  class(t_node_conf_tiles), intent(inout) :: this

  integer :: ndim, dims(p_max_dim), coords(p_max_dim), myid, i, ierr
  integer :: source, dest, j, k, id
  logical :: isperiodic(p_max_dim+1)
  integer, dimension(p_max_dim) :: lb, ub
  integer :: t

  ! For now we can leave this all the same as create_topology_mpi()
  !TODO: Either make 1d array or use graph topology

  this%free_comm = .true.
  ndim = this%x_dim

  ! Get a cartesian topology through MPI
  dims(1:ndim) = this%nx(1:ndim)
  isperiodic(1:ndim) = this%ifpr_(1:ndim)

  if ( this%ti_co%g_til_num > 1 ) then

    call mpi_cart_create( mpi_comm_world, ndim, dims, isperiodic, &
             .true., this%comm, ierr )
    if ( ierr /= MPI_SUCCESS ) then
      write(err_buf__,*) "mpi_cart_create failed";call err__("tiles/os-nconf-tiles.f03",1037)
      call abort_program( p_err_mpi )
    endif

    ! this will probably reorder the nodes so we need to get our new node id
    call mpi_comm_rank( this%comm, myid, ierr )
    if ( ierr /= MPI_SUCCESS ) then
      write(err_buf__,*) "mpi_comm_rank failed";call err__("tiles/os-nconf-tiles.f03",1044)
      call abort_program( p_err_mpi )
    endif

    this%my_aid_= myid + 1 ! we use 1 for the 1st node in osiris (i.e. fortran numbering)

    ! get local ngp
    call mpi_cart_coords( this%comm, myid, ndim, coords, ierr )
    if ( ierr /= MPI_SUCCESS ) then
      write(err_buf__,*) "mpi_cart_coords failed";call err__("tiles/os-nconf-tiles.f03",1053)
      call abort_program( p_err_mpi )
    endif
    this%my_ngp(1:this%x_dim) = coords(1:this%x_dim) + 1


    ! get neighbors; direction is 0 -> x1, 1 -> x2, 2 -> x3
    do i = 1, this%x_dim
      call mpi_cart_shift( this%comm, i-1, 1, source, dest, ierr)
      if ( ierr /= MPI_SUCCESS ) then
        write(err_buf__,*) "mpi_cart_shift failed";call err__("tiles/os-nconf-tiles.f03",1063)
        call abort_program( p_err_mpi )
      endif

      this%neighbor( p_lower, i ) = source + 1
      this%neighbor( p_upper, i ) = dest + 1
    enddo

  endif

  ! (* debug *) Test topology
  !call test_topology( this )

  ! Allocate objects for tile topology
  lb = 1
  do i = 1, p_x_dim
    ub(i) = this%ti_co%g_til(i)
  enddo

  select case ( this%x_dim )
  case(1)
    call alloc(this%ti_co%hi1, lb, ub,"tiles/os-nconf-tiles.f03",1084)
  case(2)
    call alloc(this%ti_co%hi2, lb, ub,"tiles/os-nconf-tiles.f03",1086)
  case(3)
    call alloc(this%ti_co%hi3, lb, ub,"tiles/os-nconf-tiles.f03",1088)
  end select

  ub(1) = this%x_dim
  ub(2) = this%ti_co%g_til_num

  call alloc(this%ti_co%h, lb, ub,"tiles/os-nconf-tiles.f03",1094)

  ub(1) = 3
  ub(2) = this%no_num

  call alloc(this%ti_co%g_til_aid_min_max, lb, ub,"tiles/os-nconf-tiles.f03",1099)
  call alloc(this%ti_co%g_til_aid_min_max_old, lb, ub,"tiles/os-nconf-tiles.f03",1100)

  ! Fill in hilbert and inverse hilbert mappings
  select case ( this%x_dim )
  case(1)
    do i = 1, this%ti_co%g_til_num
      this%ti_co%h(1,i) = i
      this%ti_co%hi1(i) = i
    enddo
  case(2)
    ! Populate inverse hilbert mapping
    t = 1
    call fill_hi2_hilbert( this%ti_co%hi2, this%ti_co%g_til, t )
    ! Populate hilbert mapping from inverse hilbert
    do j = 1, this%ti_co%g_til(2)
      do i = 1, this%ti_co%g_til(1)
        id = this%ti_co%hi2(i,j)
        this%ti_co%h(:,id) = (/i,j/)
      enddo
    enddo
  case(3)
    ! Populate inverse hilbert mapping
    call fill_hi3_hilbert( this%ti_co%hi3, this%ti_co%g_til )
    ! Populate hilbert mapping from inverse hilbert
    do k = 1, this%ti_co%g_til(3)
      do j = 1, this%ti_co%g_til(2)
        do i = 1, this%ti_co%g_til(1)
          id = this%ti_co%hi3(i,j,k)
          this%ti_co%h(:,id) = (/i,j,k/)
        enddo
      enddo
    enddo
  end select

  ! ! (*debug*)
  ! ! print tile topology
  ! if ( root(this) ) then
  ! call vis_hid(this)
  ! endif

end subroutine create_topology_hilbert
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Create simulation topology using snake routines
!-----------------------------------------------------------------------------------------
subroutine create_topology_snake( this )

  implicit none

  class(t_node_conf_tiles), intent(inout) :: this

  integer :: ndim, dims(p_max_dim), coords(p_max_dim), myid, i, ierr
  integer :: source, dest, j, k, id
  logical :: isperiodic(p_max_dim+1)
  integer, dimension(p_max_dim) :: lb, ub
  integer :: t

  ! For now we can leave this all the same as create_topology_mpi()
  !TODO: Either make 1d array or use graph topology

  ! -- Node topology --
  this%free_comm = .true.
  ndim = this%x_dim

  ! Get a cartesian topology through MPI
  dims(1:ndim) = this%nx(1:ndim)
  isperiodic(1:ndim) = this%ifpr_(1:ndim)

  if ( this%ti_co%g_til_num > 1 ) then

    call mpi_cart_create( mpi_comm_world, ndim, dims, isperiodic, &
             .true., this%comm, ierr )
    if ( ierr /= MPI_SUCCESS ) then
      write(err_buf__,*) "mpi_cart_create failed";call err__("tiles/os-nconf-tiles.f03",1174)
      call abort_program( p_err_mpi )
    endif

    ! this will probably reorder the nodes so we need to get our new node id
    call mpi_comm_rank( this%comm, myid, ierr )
    if ( ierr /= MPI_SUCCESS ) then
      write(err_buf__,*) "mpi_comm_rank failed";call err__("tiles/os-nconf-tiles.f03",1181)
      call abort_program( p_err_mpi )
    endif

    this%my_aid_= myid + 1 ! we use 1 for the 1st node in osiris (i.e. fortran numbering)

    ! get local ngp
    call mpi_cart_coords( this%comm, myid, ndim, coords, ierr )
    if ( ierr /= MPI_SUCCESS ) then
      write(err_buf__,*) "mpi_cart_coords failed";call err__("tiles/os-nconf-tiles.f03",1190)
      call abort_program( p_err_mpi )
    endif
    this%my_ngp(1:this%x_dim) = coords(1:this%x_dim) + 1


    ! get neighbors; direction is 0 -> x1, 1 -> x2, 2 -> x3
    do i = 1, this%x_dim
      call mpi_cart_shift( this%comm, i-1, 1, source, dest, ierr)
      if ( ierr /= MPI_SUCCESS ) then
        write(err_buf__,*) "mpi_cart_shift failed";call err__("tiles/os-nconf-tiles.f03",1200)
        call abort_program( p_err_mpi )
      endif

      this%neighbor( p_lower, i ) = source + 1
      this%neighbor( p_upper, i ) = dest + 1
    enddo

  endif

  ! (* debug *) Test topology
  !call test_topology( this )

  ! Tile topology
  lb = 1
  do i = 1, p_x_dim
    ub(i) = this%ti_co%g_til(i)
  enddo

  select case ( this%x_dim )
  case(1)
    call alloc(this%ti_co%hi1, lb, ub,"tiles/os-nconf-tiles.f03",1221)
  case(2)
    call alloc(this%ti_co%hi2, lb, ub,"tiles/os-nconf-tiles.f03",1223)
  case(3)
    call alloc(this%ti_co%hi3, lb, ub,"tiles/os-nconf-tiles.f03",1225)
  end select

  ub(1) = this%x_dim
  ub(2) = this%ti_co%g_til_num

  call alloc(this%ti_co%h, lb, ub,"tiles/os-nconf-tiles.f03",1231)

  ub(1) = 3
  ub(2) = this%no_num

  call alloc(this%ti_co%g_til_aid_min_max, lb, ub,"tiles/os-nconf-tiles.f03",1236)
  call alloc(this%ti_co%g_til_aid_min_max_old, lb, ub,"tiles/os-nconf-tiles.f03",1237)

  ! Fill in hilbert and inverse hilbert mappings
  select case ( this%x_dim )
    case(1)
      do i = 1, this%ti_co%g_til_num
        this%ti_co%h(1,i) = i
        this%ti_co%hi1(i) = i
      enddo
    case(2)
      ! Populate inverse hilbert mapping
      t = 1
      call fill_hi2_snake( this%ti_co%hi2, this%ti_co%g_til, t )
      ! Populate hilbert mapping from inverse hilbert
      do j = 1, this%ti_co%g_til(2)
        do i = 1, this%ti_co%g_til(1)
          id = this%ti_co%hi2(i,j)
          this%ti_co%h(:,id) = (/i,j/)
        enddo
      enddo
    case(3)
      ! Populate inverse hilbert mapping
      t = 1
      call fill_hi3_snake( this%ti_co%hi3, this%ti_co%g_til, t )
      ! Populate hilbert mapping from inverse hilbert
      do k = 1, this%ti_co%g_til(3)
        do j = 1, this%ti_co%g_til(2)
          do i = 1, this%ti_co%g_til(1)
            id = this%ti_co%hi3(i,j,k)
            this%ti_co%h(:,id) = (/i,j,k/)
          enddo
        enddo
      enddo
  end select

  ! ! (*debug*)
  ! ! print tile topology
  ! if ( root(this) ) then
  ! call vis_hid(this)
  ! endif


end subroutine create_topology_snake
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Create simulation topology using viktor routines
!-----------------------------------------------------------------------------------------
subroutine create_topology_viktor( this )

  implicit none

  class(t_node_conf_tiles), intent(inout) :: this

  integer :: ndim, dims(p_max_dim), coords(p_max_dim), myid, i, ierr
  integer :: source, dest, j, k, id, t
  logical :: isperiodic(p_max_dim+1)
  integer, dimension(p_max_dim) :: lb, ub

  ! For now we can leave this all the same as create_topology_mpi()
  !TODO: Either make 1d array or use graph topology

  ! -- Node topology --
  this%free_comm = .true.
  ndim = this%x_dim

  ! Get a cartesian topology through MPI
  dims(1:ndim) = this%nx(1:ndim)
  isperiodic(1:ndim) = this%ifpr_(1:ndim)

  if ( this%ti_co%g_til_num > 1 ) then

    call mpi_cart_create( mpi_comm_world, ndim, dims, isperiodic, &
             .true., this%comm, ierr )
    if ( ierr /= MPI_SUCCESS ) then
      write(err_buf__,*) "mpi_cart_create failed";call err__("tiles/os-nconf-tiles.f03",1312)
      call abort_program( p_err_mpi )
    endif

    ! this will probably reorder the nodes so we need to get our new node id
    call mpi_comm_rank( this%comm, myid, ierr )
    if ( ierr /= MPI_SUCCESS ) then
      write(err_buf__,*) "mpi_comm_rank failed";call err__("tiles/os-nconf-tiles.f03",1319)
      call abort_program( p_err_mpi )
    endif

    this%my_aid_= myid + 1 ! we use 1 for the 1st node in osiris (i.e. fortran numbering)

    ! get local ngp
    call mpi_cart_coords( this%comm, myid, ndim, coords, ierr )
    if ( ierr /= MPI_SUCCESS ) then
      write(err_buf__,*) "mpi_cart_coords failed";call err__("tiles/os-nconf-tiles.f03",1328)
      call abort_program( p_err_mpi )
    endif
    this%my_ngp(1:this%x_dim) = coords(1:this%x_dim) + 1


    ! get neighbors; direction is 0 -> x1, 1 -> x2, 2 -> x3
    do i = 1, this%x_dim
      call mpi_cart_shift( this%comm, i-1, 1, source, dest, ierr)
      if ( ierr /= MPI_SUCCESS ) then
        write(err_buf__,*) "mpi_cart_shift failed";call err__("tiles/os-nconf-tiles.f03",1338)
        call abort_program( p_err_mpi )
      endif

      this%neighbor( p_lower, i ) = source + 1
      this%neighbor( p_upper, i ) = dest + 1
    enddo

  endif

  ! (* debug *) Test topology
  !call test_topology( this )

  ! Tile topology
  lb = 1
  do i = 1, p_x_dim
    ub(i) = this%ti_co%g_til(i)
  enddo

  select case ( this%x_dim )
  case(1)
    call alloc(this%ti_co%hi1, lb, ub,"tiles/os-nconf-tiles.f03",1359)
    call alloc(this%ti_co%hi1_old, lb, ub,"tiles/os-nconf-tiles.f03",1360)
  case(2)
    call alloc(this%ti_co%hi2, lb, ub,"tiles/os-nconf-tiles.f03",1362)
    call alloc(this%ti_co%hi2_old, lb, ub,"tiles/os-nconf-tiles.f03",1363)
  case(3)
    call alloc(this%ti_co%hi3, lb, ub,"tiles/os-nconf-tiles.f03",1365)
    call alloc(this%ti_co%hi3_old, lb, ub,"tiles/os-nconf-tiles.f03",1366)
  end select

  ub(1) = this%x_dim
  ub(2) = this%ti_co%g_til_num

  call alloc(this%ti_co%h, lb, ub,"tiles/os-nconf-tiles.f03",1372)
  call alloc(this%ti_co%h_old, lb, ub,"tiles/os-nconf-tiles.f03",1373)

  ub(1) = 3
  ub(2) = this%no_num

  call alloc(this%ti_co%g_til_aid_min_max, lb, ub,"tiles/os-nconf-tiles.f03",1378)
  call alloc(this%ti_co%g_til_aid_min_max_old, lb, ub,"tiles/os-nconf-tiles.f03",1379)

  ! allocate cart_til_mm/cart_til_mm_old
  select case (p_x_dim)
  case (1)
    allocate( this%ti_co%cart_til_mm(this%x_dim) )
    allocate( this%ti_co%cart_til_mm(1)%ii(3,this%nx(1),1,1) )
    allocate( this%ti_co%cart_til_mm_old(this%x_dim) )
    allocate( this%ti_co%cart_til_mm_old(1)%ii(3,this%nx(1),1,1) )
  case (2)
    allocate( this%ti_co%cart_til_mm(this%x_dim) )
    allocate( this%ti_co%cart_til_mm(1)%ii(3,this%nx(1),1,1) )
    allocate( this%ti_co%cart_til_mm(2)%ii(3,this%nx(1),this%nx(2),1) )
    allocate( this%ti_co%cart_til_mm_old(this%x_dim) )
    allocate( this%ti_co%cart_til_mm_old(1)%ii(3,this%nx(1),1,1) )
    allocate( this%ti_co%cart_til_mm_old(2)%ii(3,this%nx(1),this%nx(2),1) )
  case (3)
    allocate( this%ti_co%cart_til_mm(this%x_dim) )
    allocate( this%ti_co%cart_til_mm(1)%ii(3,this%nx(1),1,1) )
    allocate( this%ti_co%cart_til_mm(2)%ii(3,this%nx(1),this%nx(2),1) )
    allocate( this%ti_co%cart_til_mm(3)%ii(3,this%nx(1),this%nx(2),this%nx(3)) )
    allocate( this%ti_co%cart_til_mm_old(this%x_dim) )
    allocate( this%ti_co%cart_til_mm_old(1)%ii(3,this%nx(1),1,1) )
    allocate( this%ti_co%cart_til_mm_old(2)%ii(3,this%nx(1),this%nx(2),1) )
    allocate( this%ti_co%cart_til_mm_old(3)%ii(3,this%nx(1),this%nx(2),this%nx(3)) )
  end select

  ! Fill in hilbert and inverse hilbert mappings in case doing static load balance
  select case ( this%x_dim )
  case(1)

    do i = 1, this%ti_co%g_til_num
      this%ti_co%h(1,i) = i
      this%ti_co%hi1(i) = i
    enddo

  case(2)

    ! Populate inverse hilbert mapping
    t = 1
    call fill_hi2_snake( this%ti_co%hi2, this%ti_co%g_til, t )

    ! ! print tile topology
    ! if ( root(this) ) then
    ! call vis_hid(this)
    ! endif

    ! Populate hilbert mapping from inverse hilbert
    do j = 1, this%ti_co%g_til(2)
      do i = 1, this%ti_co%g_til(1)

        id = this%ti_co%hi2(i,j)
        this%ti_co%h(:,id) = (/i,j/)

      enddo
    enddo

  case(3)

    ! Populate inverse hilbert mapping
    t = 1
    call fill_hi3_snake( this%ti_co%hi3, this%ti_co%g_til, t )

    ! print tile topology
    ! if ( root(this) ) then
      ! call vis_hid(this)
    ! endif

    ! Populate hilbert mapping from inverse hilbert
    do k = 1, this%ti_co%g_til(3)
      do j = 1, this%ti_co%g_til(2)
        do i = 1, this%ti_co%g_til(1)

          id = this%ti_co%hi3(i,j,k)
          this%ti_co%h(:,id) = (/i,j,k/)

        enddo
      enddo
    enddo

  end select

end subroutine create_topology_viktor
!-----------------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Find which node tile with g_til_aid=til belongs to:
! Need to search for where til fits into g_til_aid_min_max (which node it
! belongs to)
! Start by getting flipped list of g_til_aid_min for each node,
! n_b%ti_co%g_til_aid_min_max(p_lower,n_a%no_num:1:-1)
! Subtract 1 (so indexing begins at 0), then divide by til
! The correct node that til belongs to will be the index of the first 0
! Subtract this from no_num and add 1 to get neighbor node
!-------------------------------------------------------------------------------
pure function node_of_til( n_a, til, n_b )

  implicit none

  class( t_node_conf_tiles ), intent(in) :: n_a
  integer, intent(in) :: til
  class( t_node_conf_tiles ), intent(in), optional :: n_b

  integer :: node_of_til

  if (present(n_b)) then
    node_of_til = n_a%no_num - minloc( &
                  (n_b%ti_co%g_til_aid_min_max(p_lower,n_a%no_num:1:-1)-1)/til,&
                   dim=1) + 1
  else
    node_of_til = n_a%no_num - minloc( &
                  (n_a%ti_co%g_til_aid_min_max(p_lower,n_a%no_num:1:-1)-1)/til,&
                   dim=1) + 1
  endif

end function node_of_til
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Find which node tile with g_til_aid=til belongs to:
! Need to search for where til fits into g_til_aid_min_max (which node it
! belongs to)
! Start by getting flipped list of g_til_aid_min for each node,
! n_b%ti_co%g_til_aid_min_max(p_lower,n_a%no_num:1:-1)
! Subtract 1 (so indexing begins at 0), then divide by til
! The correct node that til belongs to will be the index of the first 0
! Subtract this from no_num and add 1 to get neighbor node
!-------------------------------------------------------------------------------
pure function node_of_til_old( n_a, til, n_b )

  implicit none

  class( t_node_conf_tiles ), intent(in) :: n_a
  integer, intent(in) :: til
  class( t_node_conf_tiles ), intent(in), optional :: n_b

  integer :: node_of_til_old

  if (present(n_b)) then
    node_of_til_old = n_a%no_num - minloc( &
                  (n_b%ti_co%g_til_aid_min_max_old(p_lower,n_a%no_num:1:-1)-1)/til,&
                   dim=1) + 1
  else
    node_of_til_old = n_a%no_num - minloc( &
                  (n_a%ti_co%g_til_aid_min_max_old(p_lower,n_a%no_num:1:-1)-1)/til,&
                   dim=1) + 1
  endif

end function node_of_til_old
!-------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Store tile configuration data for viktor topology
!-----------------------------------------------------------------------------------------
subroutine set_ti_co_viktor( this, min_til )

  implicit none

  class( t_node_conf_tiles ), intent(inout) :: this
  integer, intent(in) :: min_til

  integer :: i,j,k
  integer :: t, node
  integer, dimension(this%x_dim) :: lb, ub

  ! Set g_til_aid_min_max
  select case( this%x_dim )
  case (2)
    t = 0
    do node = 1, this%no_num
      j = this%ngp(node,1)
      k = this%ngp(node,2)
      this%ti_co%g_til_aid_min_max(1,node) = t + 1
      this%ti_co%g_til_aid_min_max(2,node) = this%ti_co%cart_til_mm(1)%ii(3,j,1,1) &
                                       *this%ti_co%cart_til_mm(2)%ii(3,j,k,1) + t
      this%ti_co%g_til_aid_min_max(3,node) = this%ti_co%g_til_aid_min_max(2,node) &
                                       - this%ti_co%g_til_aid_min_max(1,node) + 1
      t = t + this%ti_co%g_til_aid_min_max(3,node)
    enddo
  case (3)
    t = 0
    do node = 1, this%no_num
      j = this%ngp(node,1)
      k = this%ngp(node,2)
      i = this%ngp(node,3)
      this%ti_co%g_til_aid_min_max(1,node) = t + 1
      this%ti_co%g_til_aid_min_max(2,node) = this%ti_co%cart_til_mm(1)%ii(3,j,1,1) &
                                       *this%ti_co%cart_til_mm(2)%ii(3,j,k,1) &
                                       *this%ti_co%cart_til_mm(3)%ii(3,j,k,i) + t
      this%ti_co%g_til_aid_min_max(3,node) = this%ti_co%g_til_aid_min_max(2,node) &
                                       - this%ti_co%g_til_aid_min_max(1,node) + 1
      t = t + this%ti_co%g_til_aid_min_max(3,node)
    enddo
  end select

  ! Set inverse space filling operator
  t = 1
  select case( this%x_dim )
  case (2)
    do node = 1, this%no_num
      i = this%ngp(node,1)
      lb(1) = this%ti_co%cart_til_mm(1)%ii(1,i,1,1)
      ub(1) = this%ti_co%cart_til_mm(1)%ii(2,i,1,1)
      j = this%ngp(node,2)
      lb(2) = this%ti_co%cart_til_mm(2)%ii(1,i,j,1)
      ub(2) = this%ti_co%cart_til_mm(2)%ii(2,i,j,1)
      call fill_hi2_snake( this%ti_co%hi2(lb(1):ub(1),lb(2):ub(2)), ub-lb+1, t )
    enddo
  case (3)
    do node = 1, this%no_num
      i = this%ngp(node,1)
      lb(1) = this%ti_co%cart_til_mm(1)%ii(1,i,1,1)
      ub(1) = this%ti_co%cart_til_mm(1)%ii(2,i,1,1)
      j = this%ngp(node,2)
      lb(2) = this%ti_co%cart_til_mm(2)%ii(1,i,j,1)
      ub(2) = this%ti_co%cart_til_mm(2)%ii(2,i,j,1)
      k = this%ngp(node,3)
      lb(3) = this%ti_co%cart_til_mm(3)%ii(1,i,j,k)
      ub(3) = this%ti_co%cart_til_mm(3)%ii(2,i,j,k)
      call fill_hi3_snake( this%ti_co%hi3(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3)), &
        ub-lb+1, t )
    enddo
  end select

  ! Set space filling operator
  select case( this%x_dim )
    case(2)
      do j = 1, this%ti_co%g_til(2)
        do i = 1, this%ti_co%g_til(1)
          t = this%ti_co%hi2(i,j)
          this%ti_co%h(:,t) = (/i,j/)
        enddo
      enddo
    case(3)
      do k = 1, this%ti_co%g_til(3)
        do j = 1, this%ti_co%g_til(2)
          do i = 1, this%ti_co%g_til(1)
            t = this%ti_co%hi3(i,j,k)
            this%ti_co%h(:,t) = (/i,j,k/)
          enddo
        enddo
      enddo
  end select

  ! ! print tile topology
  ! if ( root(this) ) then
  ! call vis_hid(this)
  ! endif

  ! (* debug *) Validate the partition
  ! Check partition sizes
  do j = 1, this%no_num
    if ( this%ti_co%g_til_aid_min_max(3,j) < min_til ) then
      write(err_buf__,*) 'Partition too small';call err__("tiles/os-nconf-tiles.f03",1632)
      call abort_program()
    endif
  enddo

  ! Check overlap
  do j = 2, this%no_num
    if ( this%ti_co%g_til_aid_min_max(1,j) /= &
      this%ti_co%g_til_aid_min_max(2,j-1) + 1 ) then
      write(err_buf__,*) 'Partition overlap';call err__("tiles/os-nconf-tiles.f03",1641)
      call abort_program()
    endif
  enddo

end subroutine set_ti_co_viktor
!-----------------------------------------------------------------------------------------

end module m_node_conf_tiles
