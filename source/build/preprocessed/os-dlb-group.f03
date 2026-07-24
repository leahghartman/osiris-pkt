# 1 "tiles/os-dlb-group.f03"
# 1 "<built-in>" 1
# 1 "<built-in>" 3
# 467 "<built-in>" 3
# 1 "<command line>" 1
# 1 "<built-in>" 2
# 1 "tiles/os-dlb-group.f03" 2
module m_dynamic_loadbalance_group

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
# 4 "tiles/os-dlb-group.f03" 2
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
# 5 "tiles/os-dlb-group.f03" 2
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
# 6 "tiles/os-dlb-group.f03" 2

use m_system
use m_parameters
use m_node_conf, only : root, p_max, no_num, reduce_array, comm, nx
use m_node_conf_tiles, only : t_node_conf_tiles, node_of_til
use m_tile_conf, only : p_topology_viktor, g_tgp
use m_simulation_group, only : t_simulation_group
use m_species_define, only : t_species
use m_diagnostic_utilities, only : t_diag_file, p_diag_create, p_diag_grid, get_filename
use m_diagnostic_utilities, only : create_diag_file
use m_diagfile, only : p_time_length
use m_diag_utilities_tiles, only : write_vdf
use m_vdf_define, only : t_vdf, t_vdf_arr, t_vdf_report
use m_vdf_memory
use stringutil, only : idx_string
use m_species_loadbalance, only : deposit_cell_load

implicit none

private

interface report_global_load
  module procedure report_global_load_group
end interface

interface report_node_load
  module procedure report_node_load_group
end interface

interface report_grid_load
  module procedure report_grid_load_group
end interface

interface report_tile_node
  module procedure report_tile_node_group
end interface

interface report_heavy_light_tils
  module procedure report_heavy_light_tils_group
end interface

interface report_tile_load
  module procedure report_tile_load_group
end interface

public :: report_global_load, report_node_load, report_grid_load
public :: report_tile_node, report_tile_load, report_heavy_light_tils

contains

!---------------------------------------------------------------------------------------------------
! Reports global particle load information (total, min, max and avg. number of particles
! per node)
!---------------------------------------------------------------------------------------------------
subroutine report_global_load_group( sim, n )

  implicit none

  class(t_simulation_group), intent(in) :: sim
  integer, intent(in) :: n

  integer( p_int64 ), dimension(2) :: npart, tmp
  integer( p_int64 ), dimension(1) :: parts

  character(80) :: path, full_name
  class( t_species ), pointer :: species
  integer :: i, ierr

  parts(1) = 0

  ! get total particles on local node
  do i = sim%til_id(1), sim%til_id(2)
    species => sim%tiles(i)%t%part%species

    do
      if (.not. associated(species)) exit
      parts(1) = parts(1) + species%num_par
      species => species % next
    enddo
  enddo

  if ( no_num(sim%no_co) > 1 ) then

    ! get maximum and minimum number of particles on all nodes
    npart(1) = parts(1)
    npart(2) = -parts(1)

    call MPI_REDUCE( npart, tmp, 2, MPI_INTEGER8, MPI_MAX, 0, comm(sim%no_co), ierr )

    npart(1) = tmp(1)
    npart(2) = -tmp(2)

    ! get total number of particles (this needs to be done in 64 bits because the total number
    ! of particles can be larger than the 32 bit signed limit, 2^32-1 = 2.147e9)
    call MPI_REDUCE( parts, tmp, 1, MPI_INTEGER8, MPI_SUM, 0, comm(sim%no_co), ierr )
    parts(1) = tmp(1)

  else

    npart(1) = parts(1)
    npart(2) = parts(1)

  endif

  if ( root( sim%no_co ) ) then

    path = trim(path_hist) // 'LOAD' // p_dir_sep // 'GLOBAL'
    full_name = trim(path) // p_dir_sep // 'global_particle_load'

    if ( n == 0 ) then
      ! create directory
      call mkdir( path, ierr )
      ! create file
      open (unit=file_id_part_load, file=full_name, status = 'replace' , &
        form='formatted')
      ! print header
      write( file_id_part_load, '(A6, 4(1X,A18))' ) &
      'Iter', 'Total Parts.', 'Min', 'Max', 'Avg'
      write( file_id_part_load, '(A)' ) &
      '----------------------------------------------------------------------------------'
    else
      open (unit=file_id_part_load, file=full_name, position = 'append', &
        form='formatted')
    endif

    ! write particle load information to disk
    write( file_id_part_load, '(I6, 4(1X,I18))' ) &
      n, parts(1), npart(2), npart(1), parts(1) / no_num( sim%no_co )

    close(file_id_part_load)
  endif

end subroutine report_global_load_group
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! Reports number of particles per node for all nodes (small grid file)
!---------------------------------------------------------------------------------------------------
subroutine report_node_load_group( sim, n, ndump, t )

  implicit none

  class(t_simulation_group), intent(in) :: sim
  integer, intent(in) :: n, ndump
  real( p_double ), intent(in) :: t

  real( p_double ), dimension(:), pointer :: node_load
  real( p_double ), dimension(:), pointer :: f1
  real( p_double ), dimension(:,:), pointer :: f2
  real( p_double ), dimension(:,:,:), pointer :: f3

  class( t_diag_file ), allocatable :: diagFile

  integer, dimension(p_x_dim) :: lnx

  class( t_species ), pointer :: species
  integer :: i, npart, nnodes, ierr
  real( p_double ) :: npart_dbl

  node_load => null()
  f1 => null()
  f2 => null()
  f3 => null()

  ! get total number of particles on node
  npart = 0
  do i = sim%til_id(1), sim%til_id(2)
    species => sim%tiles(i)%t%part%species

    do
      if (.not. associated(species)) exit
      npart = npart + species%num_par
      species => species % next
    enddo
  enddo
  npart_dbl = npart

  ! gather data
  nnodes = no_num( sim%no_co )
  call alloc(node_load, (/nnodes/),"tiles/os-dlb-group.f03",185)

  if ( nnodes > 1 ) then
    call mpi_gather( npart_dbl, 1, MPI_DOUBLE_PRECISION, &
      node_load, 1, MPI_DOUBLE_PRECISION, &
      0, comm(sim%no_co), ierr )
    if ( ierr /= MPI_SUCCESS ) then
      write(err_buf__,*) 'MPI Error';call err__("tiles/os-dlb-group.f03",192)
      call abort_program( p_err_mpi )
    endif
  else
    node_load = npart_dbl
  endif

  ! save data to disk
  if ( root( sim%no_co ) ) then

    select type(no_co => sim%no_co); class is( t_node_conf_tiles )

      select case( no_co%ti_co%topology_type )
      case( p_topology_viktor )

        lnx(1:p_x_dim) = nx( no_co )

        ! Open the output file
        call create_diag_file( diagFile )
        diagFile%ftype = p_diag_grid

        diagFile%filepath = trim(path_mass) // 'LOAD' // p_dir_sep // 'NODE' // p_dir_sep
        diagFile%filename = trim(get_filename( n/ndump, 'node_load'))
        diagFile%name = 'Particles per node'
        diagFile%iter%n = n
        diagFile%iter%t = t
        diagFile%iter%time_units = '1 / \omega_p'

        diagFile%grid % ndims = p_x_dim
        diagFile%grid % name = "load"
        diagFile%grid % label = "particles per node"
        diagFile%grid % units = "particles"

        do i = 1, p_x_dim
          diagFile % grid % axis(i) % min = 0
          diagFile % grid % axis(i) % max = lnx(i)
          write( diagFile % grid % axis(i) % name, '(A,I0)' ) 'x', i
          write( diagFile % grid % axis(i) % label, '(A,I0)' ) 'x_', i
          diagFile % grid % axis(i) % units = 'nodes'
        enddo

        call diagFile % open( p_diag_create )

        ! Move node load values to the proper position and write node load values
        select case (p_x_dim)
        case(1)
          call alloc(f1, lnx,"tiles/os-dlb-group.f03",238)
          f1 = -1.0
          do i = 1, nnodes
            f1( no_co%ngp( i, 1 ) ) = node_load(i)
          enddo
          call diagFile % add_dataset( "load", f1 )
          call freemem(f1,"tiles/os-dlb-group.f03",244)

        case(2)
          call alloc(f2, lnx,"tiles/os-dlb-group.f03",247)
          f2 = -1.0
          do i = 1, nnodes
            f2( no_co%ngp( i, 1 ), no_co%ngp( i, 2 ) ) = node_load(i)
          enddo
          call diagFile % add_dataset( "load", f2 )
          call freemem(f2,"tiles/os-dlb-group.f03",253)

        case(3)
          call alloc(f3, lnx,"tiles/os-dlb-group.f03",256)
          f3 = -1.0
          do i = 1, nnodes
            f3( no_co%ngp( i, 1 ), no_co%ngp( i, 2 ), no_co%ngp( i, 3 ) ) = node_load(i)
          enddo
          call diagFile % add_dataset( "load", f3 )
          call freemem(f3,"tiles/os-dlb-group.f03",262)
        end select

      case default

        ! Since the topology is hilbert or snake, do a 1d node curve
        lnx(1) = nnodes

        ! Open the output file
        ! Use tile grid object with correct par_nx data
        call create_diag_file( diagFile )
        diagFile%ftype = p_diag_grid

        diagFile%filepath = trim(path_mass) // 'LOAD' // p_dir_sep // 'NODE' // p_dir_sep
        diagFile%filename = trim(get_filename( n/ndump, 'node_load'))
        diagFile%name = 'Particles per node'
        diagFile%iter%n = n
        diagFile%iter%t = t
        diagFile%iter%time_units = '1 / \omega_p'

        diagFile%grid % ndims = 1
        diagFile%grid % name = "load"
        diagFile%grid % label = "particles per node"
        diagFile%grid % units = "particles"

        diagFile % grid % axis(1) % min = 0
        diagFile % grid % axis(1) % max = lnx(1)
        write( diagFile % grid % axis(1) % name, '(A,I0)' ) 'x', 1
        write( diagFile % grid % axis(1) % label, '(A,I0)' ) 'x_', 1
        diagFile % grid % axis(1) % units = 'nodes'

        call diagFile % open( p_diag_create )

        ! Move node load values to the proper position and write node load values
        call alloc(f1, (/ lnx(1) /),"tiles/os-dlb-group.f03",296)
        f1 = -1.0
        do i = 1, nnodes
          f1( i ) = node_load(i)
        enddo
        call diagFile % add_dataset( "load", f1 )
        call freemem(f1,"tiles/os-dlb-group.f03",302)

      end select

    end select

    ! close output file
    call diagFile % close( )
    deallocate( diagFile )

  endif

  ! free load array
  call freemem(node_load,"tiles/os-dlb-group.f03",315)

end subroutine report_node_load_group
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! Reports number of particles per node for each grid cell
!---------------------------------------------------------------------------------------------------
subroutine report_grid_load_group( sim, n, ndump, t )

  implicit none

  class(t_simulation_group), intent(inout) :: sim
  integer, intent(in) :: n, ndump
  real( p_double ), intent(in) :: t

  ! local variables
  type( t_vdf_arr ), dimension(:), pointer :: load
  type( t_vdf_report ) :: load_report
  class( t_species ), pointer :: species

  integer, dimension(2,3) :: gc_num
  real( p_double ), dimension(3) :: dx
  integer :: i

  allocate( load(sim%til_id(1):sim%til_id(2)) )

  ! for the old version of p_lower_cell positions we may have particles at nx+1 if a physical
  ! boundary exists at that edge. This is just a simple hack to aacount for that situation
  gc_num( p_lower, : ) = 0
  gc_num( p_upper, : ) = 1

  ! create vdfs to hold the load
  dx = 1.0
  do i = sim%til_id(1), sim%til_id(2)
    call alloc(load(i)%v,"tiles/os-dlb-group.f03",350)
    call load(i) % v % new( p_x_dim, 1, sim%tiles(i)%t%grid%my_nx( 3, : ), gc_num, dx, .true. )
  enddo

  do i = sim%til_id(1), sim%til_id(2)
    ! get the load on each cell
    species => sim%tiles(i)%t % part % species
    do
      if (.not. associated(species)) exit
      call deposit_cell_load( species, load(i)%v )
      species => species % next
    enddo
  enddo

  ! write the data
  load_report%name = 'cell_load'
  load_report%ndump = ndump

  load_report%xname = (/'x1', 'x2', 'x3'/)
  load_report%xlabel = (/'x_1', 'x_2', 'x_3'/)
  load_report%xunits = (/'c / \omega_p', 'c / \omega_p', 'c / \omega_p'/)

  load_report%time_units = '1 / \omega_p'
  load_report%dt = 1.0

  load_report%fileLabel = ''
  load_report%path = trim(path_mass) // 'LOAD' // p_dir_sep // 'CELL'

  load_report%label = 'Particles per cell'
  load_report%units = 'particles'

  load_report%n = n
  load_report%t = t

  load_report%path = trim(path_mass) // 'LOAD' // p_dir_sep // 'GRID' // p_dir_sep
  load_report%filename = 'cell_load-' // idx_string( n/ndump, p_time_length )

  call write_vdf( load, load_report, sim%til_id(1), 1, sim%g_space, sim%grid, sim%no_co, &
                  sim%diag_patt )

  ! cleanup the vdf
  do i = sim%til_id(1), sim%til_id(2)
    call load(i)%v % cleanup()
    call freemem(load(i)%v,"tiles/os-dlb-group.f03",393)
  enddo
  deallocate(load)

end subroutine report_grid_load_group
!---------------------------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Reports which tiles are heavy (denoted by 1 in array) and which are light (denoted by 0)
!-------------------------------------------------------------------------------
subroutine report_heavy_light_tils_group( sim, n, ndump, t )

  implicit none

  class(t_simulation_group), intent(in) :: sim
  integer, intent(in) :: n
  integer, intent(in) :: ndump
  real( p_double ), intent(in) :: t

  class( t_diag_file ), allocatable :: diagFile

  real( p_single ), dimension(:), pointer :: hi1_node
  real( p_single ), dimension(:,:), pointer :: hi2_node
  real( p_single ), dimension(:,:,:), pointer :: hi3_node

  integer, dimension(p_max_dim) :: lb, ub
  integer :: i, til

  select type(no_co => sim%no_co); class is(t_node_conf_tiles)

    do i = 1, p_x_dim
      lb(i) = 1
      ub(i) = no_co%ti_co%g_til( i )
    enddo

    select case ( p_x_dim )
    case (1)

      call alloc(hi1_node, lb, ub,"tiles/os-dlb-group.f03",431)

      hi1_node = 0.0

      ! fill in heavy tiles
      do i = 1, no_co%omp_patt%n_heavy
        ! til id of heavy tile
        til = no_co%omp_patt%heavy_tils(i)
        ! fill heavy into proper position in array
        hi1_node( no_co%ti_co%h(1,til) ) = 1
      enddo

      ! fill in light tiles
      do i = 1, no_co%omp_patt%n_light
        ! til id of light tile
        til = no_co%omp_patt%light_tils(i)
        ! fill light into proper position in array
        hi1_node( no_co%ti_co%h(1,til) ) = 0
      enddo

      ! gather data
      if ( no_num( no_co ) > 1 ) then
        call reduce_array( no_co, hi1_node, operation = p_max )
      endif

    case (2)

      call alloc(hi2_node, lb, ub,"tiles/os-dlb-group.f03",458)
      hi2_node = 0.0

      ! fill in heavy tiles
      do i = 1, no_co%omp_patt%n_heavy
        ! til id of heavy tile
        til = no_co%omp_patt%heavy_tils(i)
        ! fill heavy into proper position in array
        hi2_node( no_co%ti_co%h(1,til), no_co%ti_co%h(2,til) ) = 1
      enddo

      ! fill in light tiles
      do i = 1, no_co%omp_patt%n_light
        ! til id of light tile
        til = no_co%omp_patt%light_tils(i)
        ! fill light into proper position in array
        hi2_node( no_co%ti_co%h(1,til), no_co%ti_co%h(2,til) ) = 0
      enddo

      ! gather data
      if ( no_num( no_co ) > 1 ) then
        call reduce_array( no_co, hi2_node, operation = p_max )
      endif

    case (3)

      call alloc(hi3_node, lb, ub,"tiles/os-dlb-group.f03",484)
      hi3_node = 0.0

      ! fill in heavy tiles
      do i = 1, no_co%omp_patt%n_heavy
        ! til id of heavy tile
        til = no_co%omp_patt%heavy_tils(i)
        ! fill heavy into proper position in array
        hi3_node( no_co%ti_co%h(1,til), no_co%ti_co%h(2,til), no_co%ti_co%h(3,til) ) = 1
      enddo

      ! fill in light tiles
      do i = 1, no_co%omp_patt%n_light
        ! til id of light tile
        til = no_co%omp_patt%light_tils(i)
        ! fill light into proper position in array
        hi3_node( no_co%ti_co%h(1,til), no_co%ti_co%h(2,til), no_co%ti_co%h(3,til) ) = 0
      enddo

      ! gather data
      if ( no_num( no_co ) > 1 ) then
        call reduce_array( no_co, hi3_node, operation = p_max )
      endif

    end select

    ! save data to disk
    if ( root( no_co ) ) then

      ! Use tile grid object with correct par_nx data
      call create_diag_file( diagFile )
      diagFile%ftype = p_diag_grid

      diagFile%filepath = trim(path_mass) // 'LOAD' // p_dir_sep // 'HEAVY-LIGHT' // p_dir_sep
      diagFile%filename = trim(get_filename( n/ndump, 'heavy_light'))
      diagFile%name = 'Heavy-light tile topology'
      diagFile%iter%n = n
      diagFile%iter%t = t
      diagFile%iter%time_units = '1 / \omega_p'

      diagFile%grid % ndims = p_x_dim
      diagFile%grid % name = "topology"
      diagFile%grid % label = "heavy-light tile topology"
      diagFile%grid % units = "light=0 heavy=1"

      do i = 1, p_x_dim
        diagFile % grid % axis(i) % min = 0
        diagFile % grid % axis(i) % max = ub(i)
        write( diagFile % grid % axis(i) % name, '(A,I0)' ) 'x', i
        write( diagFile % grid % axis(i) % label, '(A,I0)' ) 'x_', i
        diagFile % grid % axis(i) % units = 'tiles'
      enddo

      call diagFile % open( p_diag_create )

      ! move node load values to the proper position and write node load values
      select case (p_x_dim)
      case(1)
        call diagFile % add_dataset( "topology", hi1_node )
        call freemem(hi1_node,"tiles/os-dlb-group.f03",543)
      case(2)
        call diagFile % add_dataset( "topology", hi2_node )
        call freemem(hi2_node,"tiles/os-dlb-group.f03",546)
      case(3)
        call diagFile % add_dataset( "topology", hi3_node )
        call freemem(hi3_node,"tiles/os-dlb-group.f03",549)
      end select

      ! close output file
      call diagFile % close( )
      deallocate( diagFile )

    endif

    ! all nodes cleanup
    select case (p_x_dim)
    case(1)
      call freemem(hi1_node,"tiles/os-dlb-group.f03",561)
    case(2)
      call freemem(hi2_node,"tiles/os-dlb-group.f03",563)
    case(3)
      call freemem(hi3_node,"tiles/os-dlb-group.f03",565)
    end select

  end select

end subroutine report_heavy_light_tils_group
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Reports which tiles are on which node
!-------------------------------------------------------------------------------
subroutine report_tile_node_group( this, n, ndump, t )

  implicit none

  class( t_simulation_group ), intent(in) :: this
  integer, intent(in) :: n, ndump
  real( p_double ), intent(in) :: t

  class( t_diag_file ), allocatable :: diagFile

  integer :: node
  integer :: i
  integer, dimension(p_max_dim) :: lb, ub

  ! inverse space filling operator for nodes, coordinate maps to node
  real( p_single ), dimension(:), pointer :: hi1_node
  real( p_single ), dimension(:,:), pointer :: hi2_node
  real( p_single ), dimension(:,:,:), pointer :: hi3_node

  select type(no_co => this%no_co); class is(t_node_conf_tiles)

  ! -- Gather data --
  do i = 1, p_x_dim
    lb(i) = 1
    ub(i) = no_co%ti_co%g_til( i )
  enddo

  select case ( p_x_dim )
  case (1)
    call alloc(hi1_node, lb, ub,"tiles/os-dlb-group.f03",605)
    do i = 1, product(no_co%ti_co%g_til)
      ! find node that g_til_aid=i belongs to
      node = node_of_til( no_co, i )
      ! fill tile's element of inv. space filling op. with node of tile
      hi1_node( g_tgp(no_co%ti_co,i,1) ) = node
    enddo
  case (2)
    call alloc(hi2_node, lb, ub,"tiles/os-dlb-group.f03",613)
    do i = 1, product(no_co%ti_co%g_til)
      ! find node that g_til_aid=i belongs to
      node = node_of_til( no_co, i )
      ! fill tile's element of inv. space filling op. with node of tile
      hi2_node( g_tgp(no_co%ti_co,i,1), g_tgp(no_co%ti_co,i,2) ) = node
    enddo
  case (3)
    call alloc(hi3_node, lb, ub,"tiles/os-dlb-group.f03",621)
    do i = 1, product(no_co%ti_co%g_til)
      ! find node that g_til_aid=i belongs to
      node = node_of_til( no_co, i )
      ! fill tile's element of inv. space filling op. with node of tile
      hi3_node( g_tgp(no_co%ti_co,i,1), g_tgp(no_co%ti_co,i,2), g_tgp(no_co%ti_co,i,3) ) = node
    enddo
  end select

  end select

  ! -- Save data to disk --
  ! Use tile grid object with correct par_nx data
  call create_diag_file( diagFile )
  diagFile%ftype = p_diag_grid

  diagFile%filepath = trim(path_mass) // 'LOAD' // p_dir_sep // 'TILE-NODE' // p_dir_sep
  diagFile%filename = trim(get_filename( n/ndump, 'tile_node'))
  diagFile%name = 'Node topology'
  diagFile%iter%n = n
  diagFile%iter%t = t
  diagFile%iter%time_units = '1 / \omega_p'

  diagFile%grid % ndims = p_x_dim
  diagFile%grid % name = "topology"
  diagFile%grid % label = "node topology"
  diagFile%grid % units = "node id"

  do i = 1, p_x_dim
    diagFile % grid % axis(i) % min = 0
    diagFile % grid % axis(i) % max = ub(i)
    write( diagFile % grid % axis(i) % name, '(A,I0)' ) 'x', i
    write( diagFile % grid % axis(i) % label, '(A,I0)' ) 'x_', i
    diagFile % grid % axis(i) % units = 'tiles'
  enddo

  call diagFile % open( p_diag_create )

  ! move node load values to the proper position and write node load values
  select case (p_x_dim)
  case(1)
    call diagFile % add_dataset( "topology", hi1_node )
    call freemem(hi1_node,"tiles/os-dlb-group.f03",663)
  case(2)
    call diagFile % add_dataset( "topology", hi2_node )
    call freemem(hi2_node,"tiles/os-dlb-group.f03",666)
  case(3)
    call diagFile % add_dataset( "topology", hi3_node )
    call freemem(hi3_node,"tiles/os-dlb-group.f03",669)
  end select

  ! close output file
  call diagFile % close( )
  deallocate( diagFile )

end subroutine report_tile_node_group
!-------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! Reports number of particles per tile for all tiles (small grid file)
!---------------------------------------------------------------------------------------------------
subroutine report_tile_load_group( sim, n, ndump, t )

  implicit none

  class(t_simulation_group), intent(in) :: sim
  integer, intent(in) :: n, ndump
  real( p_double ), intent(in) :: t

  real( p_double ), dimension(:), pointer :: tile_load
  real( p_double ), dimension(:), pointer :: f1
  real( p_double ), dimension(:,:), pointer :: f2
  real( p_double ), dimension(:,:,:), pointer :: f3

  class( t_diag_file ), allocatable :: diagFile

  integer, dimension(p_x_dim) :: lnx

  class( t_species ), pointer :: species
  integer :: i, npart, nnodes, g_til_num

  ! allocate arrays to hold particles per tile
  select type(no_co => sim%no_co); class is(t_node_conf_tiles)
    g_til_num = no_co%ti_co%g_til_num
    do i = 1, p_x_dim
      lnx(i) = no_co%ti_co%g_til(i)
    enddo
  end select
  call alloc(tile_load, (/g_til_num/),"tiles/os-dlb-group.f03",709)
  tile_load = 0.0

  ! get total number of particles on each tile
  do i = sim%til_id(1), sim%til_id(2)
    npart = 0
    species => sim%tiles(i)%t%part%species

    do
      if (.not. associated(species)) exit
      npart = npart + species%num_par
      species => species % next
    enddo
    tile_load(i) = npart
  enddo

  ! gather data
  nnodes = no_num( sim%no_co )
  if ( nnodes > 1 ) then
    call reduce_array( sim%no_co, tile_load, operation = p_max )
  endif

  ! save data to disk
  if ( root( sim%no_co ) ) then

    ! Open the output file
    call create_diag_file( diagFile )
    diagFile%ftype = p_diag_grid

    diagFile%filepath = trim(path_mass) // 'LOAD' // p_dir_sep // 'TILE' // p_dir_sep
    diagFile%filename = trim(get_filename( n/ndump, 'tile_load'))
    diagFile%name = 'Particles per tile'
    diagFile%iter%n = n
    diagFile%iter%t = t
    diagFile%iter%time_units = '1 / \omega_p'

    diagFile%grid % ndims = p_x_dim
    diagFile%grid % name = "load"
    diagFile%grid % label = "particles per tile"
    diagFile%grid % units = "particles"

    do i = 1, p_x_dim
      diagFile % grid % axis(i) % min = 0
      diagFile % grid % axis(i) % max = lnx(i)
      write( diagFile % grid % axis(i) % name, '(A,I0)' ) 'x', i
      write( diagFile % grid % axis(i) % label, '(A,I0)' ) 'x_', i
      diagFile % grid % axis(i) % units = 'tiles'
    enddo

    call diagFile % open( p_diag_create )

    ! move tile load values to the proper position and write tile load values
    select type(no_co => sim%no_co); class is(t_node_conf_tiles)
      select case (p_x_dim)
      case(1)
        call alloc(f1, lnx,"tiles/os-dlb-group.f03",764)
        f1 = -1.0
        do i = 1, g_til_num
          f1( no_co%ti_co%h(1,i) ) = tile_load(i)
        enddo
        call diagFile % add_dataset( "load", f1 )
        call freemem(f1,"tiles/os-dlb-group.f03",770)
      case(2)
        call alloc(f2, lnx,"tiles/os-dlb-group.f03",772)
        f2 = -1.0
        do i = 1, g_til_num
          f2( no_co%ti_co%h(1,i), no_co%ti_co%h(2,i) ) = tile_load(i)
        enddo
        call diagFile % add_dataset( "load", f2 )
        call freemem(f2,"tiles/os-dlb-group.f03",778)
      case(3)
        call alloc(f3, lnx,"tiles/os-dlb-group.f03",780)
        f3 = -1.0
        do i = 1, g_til_num
          f3( no_co%ti_co%h(1,i), no_co%ti_co%h(2,i), no_co%ti_co%h(3,i) ) = tile_load(i)
        enddo
        call diagFile % add_dataset( "load", f3 )
        call freemem(f3,"tiles/os-dlb-group.f03",786)
      end select
    end select

    ! close output file
    call diagFile % close( )
    deallocate( diagFile )

  endif

  ! free load array
  call freemem(tile_load,"tiles/os-dlb-group.f03",797)

end subroutine report_tile_load_group
!---------------------------------------------------------------------------------------------------

end module m_dynamic_loadbalance_group

!-------------------------------------------------------------------------------
! Report on the load of the simulation for tiles algorithm
!-------------------------------------------------------------------------------
subroutine report_load_group( this )

  use m_dynamic_loadbalance_group
  use m_simulation_group, only : t_simulation_group
  use m_simulation, only : report_load_ev
  use m_grid_define, only : p_global, p_node, p_grid
  use m_grid_tiles, only : p_tile_node, p_tile_load, p_hl_tils
  use m_time_step, only : n
  use m_time, only : t
  use m_node_conf, only : root
  use m_logprof, only : begin_event, end_event

  implicit none

  class( t_simulation_group ), intent(inout) :: this

  call begin_event(report_load_ev)

  ! global load : total, max, min, avg parts per node
  if ( this%grid%if_report( n(this%tstep), p_global ) ) then
     call report_global_load( this, n(this%tstep) )
  endif

  ! node load : number of particles per node for all nodes
  if ( this%grid%if_report( n(this%tstep), p_node ) ) then
     call report_node_load( this, n(this%tstep), this%grid%ndump( p_node ), t(this%time) )
  endif

  ! grid load : number of particles per cell for all grid points
  if ( this%grid%if_report( n(this%tstep), p_grid ) ) then
     call report_grid_load( this, n(this%tstep), this%grid%ndump( p_grid ), t(this%time) )
  endif

  ! root node reports which tiles are on which node
  if ( this%grid%if_report( n(this%tstep), p_tile_node) .and. root(this%no_co) ) then
    call report_tile_node( this, n(this%tstep), this%grid%ndump( p_tile_node ), t(this%time) )
  endif

  ! root node reports distribution of heavy and light tiles accross simulation box
  if ( this%grid%if_report( n(this%tstep), p_hl_tils) ) then
    call report_heavy_light_tils( this, n(this%tstep), this%grid%ndump( p_hl_tils ), t(this%time) )
  endif

  ! tile load : number of particles per tile for all nodes
  if ( this%grid%if_report( n(this%tstep), p_tile_load) ) then
    call report_tile_load( this, n(this%tstep), this%grid%ndump( p_tile_load ), t(this%time) )
  endif

  call end_event(report_load_ev)

end subroutine report_load_group
!-------------------------------------------------------------------------------
