# 1 "gr/os-species-gr.f03"
# 1 "<built-in>" 1
# 1 "<built-in>" 3
# 467 "<built-in>" 3
# 1 "<command line>" 1
# 1 "<built-in>" 2
# 1 "gr/os-species-gr.f03" 2
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
# 2 "gr/os-species-gr.f03" 2
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
# 3 "gr/os-species-gr.f03" 2

module m_species_gr

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
# 7 "gr/os-species-gr.f03" 2

use m_system
use m_parameters

use m_geometry_gr
use m_emf_define_gr

end module m_species_gr

!-------------------------------------------------------------------------------
! Read information from input file
!-------------------------------------------------------------------------------
subroutine read_input_species_gr( this, input_file, def_name, periodic, &
  if_move, grid, dt, read_prof, sim_options )

  use m_system
  use m_parameters
  use m_input_file

  use m_species_define
  use m_species_diagnostics
  use m_species_boundary
  use m_species_udist
  use m_grid_define
  use m_psource_gr

  use m_species_define_gr

  implicit none

  class( t_species_gr ), intent(inout) :: this
  class( t_input_file ), intent(inout) :: input_file
  character(len = *), intent(in) :: def_name
  logical, dimension(:), intent(in) :: periodic, if_move
  class( t_grid ), intent(in) :: grid
  type(t_options), intent(in) :: sim_options

  ! local variables
  real(p_double), intent(in) :: dt
  logical, intent(in) :: read_prof

  character(len = p_max_spname_len) :: name
  character(len = 10 ) :: push_type

  integer :: num_par_max
  integer :: n_sort
  real(p_k_part) :: rqm

  real(p_k_part) :: q_real

  integer, dimension(p_x_dim) :: num_par_x
  integer(p_int64), dimension(p_x_dim) :: tot_par_x

  logical :: add_tag
  logical :: free_stream

  real(p_double) :: push_start_time

  logical :: init_fields

  character(len = 16) :: init_type

  namelist /nl_species/ name, num_par_max, n_sort, rqm, q_real, &
    num_par_x, tot_par_x, push_type, push_start_time, add_tag, free_stream, &
    init_fields, init_type

  integer :: i, ierr

  ! copy in relevant simulation parameters
  this%dt = dt
  this%g_nx( 1:p_x_dim ) = grid%g_nx( 1:p_x_dim )
  this%coordinates = grid%coordinates

  push_type = "standard"

  ! Initialization type
  init_type = "standard"

  num_par_max = 0
  rqm = 0
  q_real = 0
  num_par_x = -1
  tot_par_x = -1
  n_sort = 25

  push_start_time = -1.0_p_double
  name = trim(adjustl(def_name))

  add_tag = .false.
  free_stream = .false.
  init_fields = .false.

  ! Get namelist text from input file
  call get_namelist( input_file, "nl_species", ierr )
  if (ierr /= 0) then
    if ( mpi_node() == 0 ) then
      if (ierr < 0) then
        write(0,*) "Error reading species parameters"
      else
        write(0,*) "Error: species parameters missing"
      endif
      write(0,*) "aborting..."
    endif
    stop
  endif


  read (input_file%nml_text, nml = nl_species, iostat = ierr)
  if (ierr /= 0) then
    if ( mpi_node() == 0 ) then
      write(0,*) "   Error reading species parameters"
      write(0,*) "   aborting..."
    endif
    stop
  endif

  if (disp_out(input_file)) then
  if (mpi_node()==0) print *,"   Species name : ", trim(name)
  endif

  ! -------------- Allocate any sub-objects needed --------------------
  call this%allocate_objs()
  this%energy => null()

  this%name = name
  this%free_stream = free_stream

  this%num_par_max = num_par_max
  this%rqm = rqm

  ! Validate num_par_x and tot_par_x
  do i = 1, p_x_dim
    if ( (num_par_x(i) <= 0) .and. (tot_par_x(i) <= 0) ) then
      if ( mpi_node() == 0 ) then
        write(0,*) "   Error reading species parameters"
        write(0,*) "   Number of particles per cell must be >= 1 in all directions:"
        write(0,"(A,I0,A,I0,A,I0,A,I0)") "    -> num_par_x(",i,") = ", &
        num_par_x(i), ", tot_par_x(",i,") = ", tot_par_x(i)
        write(0,*) "   aborting..."
      endif
      stop
    endif

    ! num_par_x overrides tot_par_x
    if ( num_par_x(i) <= 0 ) then
      this%num_par_x(i) = tot_par_x(i) / this % g_nx(i)
      this%tot_par_x(i) = tot_par_x(i)
    else
      this%num_par_x(i) = num_par_x(i)
      this%tot_par_x(i) = num_par_x(i) * this % g_nx(i)
    endif
  enddo

  this%n_sort = n_sort

  ! Pusher type
  select case ( trim( push_type ) )
  case('standard')
    this%push_type = p_std
  case('vay')
    this%push_type = p_vay
  case('rk4')
    this%push_type = p_rk4
  case('rk6')
    this%push_type = p_rk6
  end select

  if ( free_stream ) then
    this%push_type = p_std
  endif

  this%push_start_time = push_start_time
  this%add_tag = add_tag

  ! this is also used by some people even without collisions
  ! q_real and rqm having different signs is inconsistent, give rqm precedence
  this%q_real = sign(q_real, rqm)

  this%if_collide = .false.
  this%if_like_collide = .false.

  ! read momentum distribution
  call read_nml( this%udist, input_file )

  ! Initialize fields from initial density / momentum profile
  this%init_fields = init_fields

  ! read profile definition unless read_prof was set to false
  if ( read_prof ) then
    this % source => new_source( init_type, input_file, this%coordinates )
    if ( .not. associated( this%source )) then
      if ( mpi_node() == 0) then
        write(0,*) "   Error reading species parameters"
        write(0,*) "   Invalid 'init_type' value (", trim(init_type), ")"
        write(0,*) "   aborting..."
      endif
      stop
    endif
  else
    this % source => null()
  endif

  ! read boundary condtions
  call read_nml( this%bnd_con, input_file, periodic, if_move, this%coordinates )

  this%num_pistons = 0

  ! set emission diagnostic flag to false by default
  this%track_nemit = .false.

  ! read diagnostics
  call this%diag%read_input( input_file, sim_options%gamma )

end subroutine read_input_species_gr
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine init_species_gr( this, sp_id, interpolation, grid_center, grid, g_space, &
  emf, jay, no_co, send_vdf, recv_vdf, bnd_cross, node_cross, send_spec, &
  recv_spec, ndump_fac, restart, restart_handle, sim_options, tstep, tmin, tmax )

  use m_system
  use m_species_define_gr
  use m_species_define
  use m_grid_define, only : t_grid
  use m_space
  use m_emf_define, only : t_emf
  use m_current_define, only : t_current
  use m_vdf_define
  use m_vdf_memory, only : alloc
  use m_vdf_comm, only : t_vdf_msg
  use m_node_conf
  use m_restart
  use m_time_step, only : t_time_step

  use m_parameters

  !use m_species, only : restart_read
  use m_species_diagnostics
  use m_species_boundary
  use m_species_initfields
  use m_species_push_gr
  use m_species_current_gr

  implicit none

  class( t_species_gr ), intent(inout) :: this

  integer, intent(in) :: sp_id
  integer, intent(in) :: interpolation
  logical, intent(in) :: grid_center
  class( t_grid ), intent(in) :: grid
  type( t_space ), intent(in) :: g_space
  class( t_emf ), intent(inout) :: emf
  class( t_current ), intent(inout) :: jay
  class( t_node_conf ), intent(in) :: no_co
  type(t_vdf_msg), dimension(2), intent(inout) :: send_vdf, recv_vdf
  type( t_part_idx ), dimension(2), intent(inout) :: bnd_cross
  type( t_part_idx ), intent(inout) :: node_cross
  type( t_spec_msg ), dimension(2), intent(inout) :: send_spec, recv_spec
  integer, intent(in) :: ndump_fac
  logical, intent(in) :: restart
  type( t_restart_handle ), intent(in) :: restart_handle
  type(t_options), intent(in) :: sim_options
  type( t_time_step ), intent(in) :: tstep
  real(p_double), intent(in) :: tmin, tmax

  integer :: i

  ! Interpolation level is now defined by the global particles object
  this%interpolation = interpolation

  ! Field values are grid centered
  this%grid_center = grid_center

  if ( this%grid_center .and. this%push_type == p_simd ) then
    if (mpi_node()==0) then
      write(0,*) '(*warning*) Grid centered field values are not supported by SIMD code'
      write(0,*) '(*warning*) defaulting to the standard fortran pusher.'
    endif
    this%push_type = p_std
  endif

  ! initialize diagnostics
  ! dianostics must be setup before inject_area because of particle tracking
  call this%diag%init(this%name,this% get_n_x_dims(), ndump_fac, interpolation, restart, &
                      restart_handle )

  ! allocate aux. array for openMP particle energy calculations
  ! (if OpenMP support was not compiled the function n_threads always returns 1)
  call alloc(this%energy, (/ n_threads(no_co) /),"gr/os-species-gr.f03",297)

  ! initialize boundary conditions
  call setup( this%bnd_con, g_space, no_co, grid%coordinates, restart, restart_handle )

  ! Get unique id based on node grid position (independent of mpi rank)
  ! This needs to happen before calling inject_area
  this%ngp_id = no_co % ngp_id()


  ! Ensure that the particles buffer size, if set in the input file, is a multiple of 8 because
  ! of vector code
  this%num_par_max = ((this%num_par_max + 7) / 8 ) *8


  if ( restart ) then
  call this % restart_read( restart_handle )
  else
  this%sp_id = sp_id

  this%coordinates = grid%coordinates

  ! process position type (only linear interpolation is available)
  this%pos_type = p_cell_low

  ! Global box boundaries are shifted from the simulation box boundaries by + 0.5*dx to simplify
  ! global particle position calculations, which are only used for injection (density profile)
  ! and diagnostics
  do i = 1, p_x_dim
    this%dx( i ) = emf%e%dx(i)
    this%g_box( p_lower, i ) = xmin( g_space, i ) + 0.5_p_double * this%dx(i)
    this%g_box( p_upper, i ) = xmax( g_space, i ) + 0.5_p_double * this%dx(i)
  enddo

  this%g_nx( 1:p_x_dim ) = grid%g_nx( 1:p_x_dim )

  do i = 1, p_x_dim
    this%total_xmoved( i ) = total_xmoved( g_space, i, p_lower )
  enddo

  this%my_nx_p(:,1:p_x_dim) = grid%my_nx(:,1:p_x_dim)

  ! allocate needed buffers
  if ( this%num_par_max > 0 .and. this%init_part ) then
    ! Check if a superclass did not already allocate the memory
    if (.not. associated(this%x)) call this%init_buffer( this%num_par_max )
  endif

  ! actual initialization of particles of this species
  this%num_par = 0
  this%num_created = 0

  call this%inject( this%my_nx_p, jay, no_co, bnd_cross, node_cross, send_spec, recv_spec )

  ! Initialize fields associated with the initial distribution / momentum
  if ( this%init_fields ) then
    call init_fields( this, emf, no_co, grid, send_vdf, recv_vdf )
  endif

  endif

  ! Setup data that is not in the restart file

  ! Set pointer to dudt function
  select case ( this % push_type )
  case( p_std )
    this % dudt => dudt_boris_gr
  case( p_vay )
    ! this % dudt => dudt_vay_gr
  case( p_rk4 )
    !this % dudt => dudt_rk_gr
    ! gets chosen anyway via advance func
  case( p_rk6 )
    !this % dudt => dudt_rk_gr
    ! gets chosen anyway via advance func

  case default
    write(err_buf__,*) "Invalid pusher type";call err__("gr/os-species-gr.f03",374)
    call abort_program( p_err_invalid )
  end select

  ! Set pointer to current deposition functions
  select case ( this%interpolation )
  case( p_linear )
    this % dep_current_1d => null()
    this % dep_current_2d => null() !dep_current_2d_gr
    this % dep_current_3d => null()
  case default
    write(err_buf__,*) "Invalid interpolation type";call err__("gr/os-species-gr.f03",385)
    call abort_program( p_err_invalid )
  end select

  ! calculate m_real (this is currently only used by collisions)
  this%m_real = this%q_real * this%rqm
  if( this%q_real /= 0.0_p_k_part ) this%rq_real = 1.0_p_k_part / this%q_real

end subroutine init_species_gr
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine push_species_gr( this, emf, jay, t, tstep, tid, n_threads, &
  options )

  use m_system
  use m_parameters

  use m_species_define_gr, only : t_species_gr
  use m_emf_define, only : t_emf
  use m_vdf_define, only : t_vdf
  use m_species_push_gr, only : advance_deposit_gr

  use m_time_step
  ! use m_species_push_gr
  use m_species_radcool

  implicit none

  class( t_species_gr ), intent(inout) :: this
  class( t_emf ), intent(inout) :: emf
  type( t_vdf ), intent(inout) :: jay

  real(p_double), intent(in) :: t
  type( t_time_step ), intent(in) :: tstep
  type( t_options ), intent(in) :: options

  integer, intent(in) :: tid ! local thread id
  integer, intent(in) :: n_threads ! total number of threads

  integer :: push_type

  ! if before push start time return silently
  if ( t >= this%push_start_time ) then

    push_type = this%push_type
    select case ( push_type )
      case default
        call advance_deposit_gr( this, emf, jay, t, tstep, tid, n_threads )
    end select

  endif

end subroutine push_species_gr
!-----------------------------------------------------------------------------------------
