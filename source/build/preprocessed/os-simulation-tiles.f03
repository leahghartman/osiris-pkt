# 1 "tiles/os-simulation-tiles.f03"
# 1 "<built-in>" 1
# 1 "<built-in>" 3
# 467 "<built-in>" 3
# 1 "<command line>" 1
# 1 "<built-in>" 2
# 1 "tiles/os-simulation-tiles.f03" 2
module m_simulation_tiles

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
# 4 "tiles/os-simulation-tiles.f03" 2
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
# 5 "tiles/os-simulation-tiles.f03" 2

use m_system
use m_parameters
use m_simulation, only : t_simulation
use m_node_conf, only : t_node_conf, periodic
use m_node_conf_tiles, only : t_node_conf_tiles
use m_grid_tiles, only : t_grid_tiles
use m_restart, only : read_nml
use m_space, only : read_nml, if_move, xmin, xmax
use m_time, only : read_nml
use m_time_step, only : read_nml, dt
use m_input_file, only : t_input_file

implicit none

private

type, extends( t_simulation ) :: t_simulation_tiles

contains

  procedure :: cleanup => cleanup_sim_tiles
  procedure :: allocate_objs => allocate_objs_tiles
  procedure :: read_copy_sim

end type t_simulation_tiles

public :: t_simulation_tiles

contains

!-----------------------------------------------------------------------------------------
! Cleanup simulation objects
!-----------------------------------------------------------------------------------------
subroutine cleanup_sim_tiles( sim )

  implicit none

  class( t_simulation_tiles ), intent(inout) :: sim

  call sim % grid % cleanup()

  ! cleanup particles
  call sim % part % cleanup( )

  ! cleanup current object
  call sim%jay%cleanup()

  ! cleanup emf object
  call sim%emf%cleanup()

  ! cleanup laser pulses
  call sim % zpulse_list % cleanup( )

  ! cleanup antennas
  call sim % antenna_array % cleanup( )

  ! cleanup vdf and spec comm buffers
  call sim % bnd % cleanup()

  sim%no_co%free_comm = .false.
  call sim % no_co % cleanup( barrier=.false. )

  ! deallocate polymorphic objects
  deallocate( sim%no_co )
  deallocate( sim%grid )
  deallocate( sim%emf )
  deallocate( sim%part )
  deallocate( sim%jay )
  deallocate( sim%zpulse_list )
  deallocate( sim%bnd )

end subroutine cleanup_sim_tiles

!-----------------------------------------------------------------------------------------
! Initialize simulation objects
!-----------------------------------------------------------------------------------------
subroutine allocate_objs_tiles( sim )

  implicit none

  class( t_simulation_tiles ), intent(inout) :: sim

  if ( .not. associated( sim%no_co ) ) then
    allocate( t_node_conf_tiles :: sim%no_co )
  endif

  if ( .not. associated( sim%grid ) ) then
    allocate( t_grid_tiles :: sim%grid )
  endif

  call sim % t_simulation % allocate_objs()

end subroutine allocate_objs_tiles

!-----------------------------------------------------------------------------------------
! Read in tile data after having read in group data (and copy options)
!-----------------------------------------------------------------------------------------
subroutine read_copy_sim( this, sim, input_file, i )

  implicit none

  class( t_simulation_tiles ), intent(inout) :: this
  class( t_simulation ), intent(in) :: sim
  class( t_input_file ), intent(inout) :: input_file
  integer, intent(in) :: i

  class( t_node_conf ), pointer :: nc_til
  real( p_double ), dimension( p_x_dim ) :: ldx

  this % options = sim % options

  ! We don't read no_co here because we set it up later and copy over information

  ! set aid for tile
  nc_til => this%no_co; select type(nc_til); class is(t_node_conf_tiles)
    nc_til % ti_co % g_til_aid = i
  end select

  call this % grid % read_input( input_file, sim%no_co%nx )

  call read_nml( this%tstep, input_file )

  call read_nml( this%restart, input_file, this%options%restart )

  ! We only have one g_space to avoid roundoff error in moving window
  ! But we'll read it in here for namelist reasons and for test_input function
  call read_nml( this%g_space, input_file, periodic(sim%no_co), this%grid%coordinates )

  call read_nml( this%time, input_file )

  ! determine cell size (use sim objects for this)
  ldx = (xmax( sim%g_space ) - xmin( sim%g_space ))/ sim%grid%g_nx( 1:p_x_dim )
  call this % emf % read_input( input_file, periodic(sim%no_co), if_move(sim%g_space), &
                                this%grid, ldx, dt(this%tstep), this%options%gamma )

  call this%part%read_input( input_file, periodic(sim%no_co), &
                  if_move(sim%g_space), this%grid, dt(this%tstep), &
                  this%options )

  call this % zpulse_list % read_input( input_file, sim%g_space, this%emf%bnd_con, &
                 periodic(sim%no_co), this%grid, this%options )

  call this % jay % read_input( input_file )

  call this % antenna_array % read_input( input_file )

end subroutine read_copy_sim
!-------------------------------------------------------------------------------


end module m_simulation_tiles
