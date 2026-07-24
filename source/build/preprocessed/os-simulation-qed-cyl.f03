# 1 "qed-cyl/os-simulation-qed-cyl.f03"
# 1 "<built-in>" 1
# 1 "<built-in>" 3
# 467 "<built-in>" 3
# 1 "<command line>" 1
# 1 "<built-in>" 2
# 1 "qed-cyl/os-simulation-qed-cyl.f03" 2
module m_simulation_qedcyl

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
# 4 "qed-cyl/os-simulation-qed-cyl.f03" 2
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
# 5 "qed-cyl/os-simulation-qed-cyl.f03" 2

use m_system
use m_simulation_cyl_modes, only: t_simulation_cyl_modes
use m_particles_qedcyl, only: t_particles_qedcyl

implicit none

private

type, extends( t_simulation_cyl_modes ) :: t_simulation_qedcyl

contains
  procedure :: allocate_objs => allocate_objs_qedcyl
  procedure :: report_load => report_load_sim_qedcyl

end type t_simulation_qedcyl

public :: t_simulation_qedcyl

contains

!-----------------------------------------------------------------------------------------
! Initialize simulation objects
!-----------------------------------------------------------------------------------------
subroutine allocate_objs_qedcyl( sim )

  implicit none

  class( t_simulation_qedcyl ), intent(inout) :: sim

  if (mpi_node()==0) print *,'Allocating QEDCYL objects.'

  allocate( t_particles_qedcyl :: sim%part )

  ! allocate any remaining objects in superclass
  call sim % t_simulation_cyl_modes % allocate_objs()

end subroutine allocate_objs_qedcyl


!-----------------------------------------------------------------------------------------
! Report on the load of the simulation, including particles + photons.
!-----------------------------------------------------------------------------------------
subroutine report_load_sim_qedcyl( this )

  use m_grid_define
  use m_particles
  use m_dynamic_loadbalance
  use m_simulation
  use m_time_step
  use m_time
  use m_logprof
  use m_particles_qedcyl

  implicit none

  class ( t_simulation_qedcyl ), intent(inout) :: this

  call begin_event(report_load_ev)

  ! global load : total, max, min, avg parts per node
  if ( this%grid%if_report( n(this%tstep), p_global ) ) then
      call report_global_load( this%part, n(this%tstep), this%no_co )
  endif

  ! node load : number of particles per node for all nodes
  if ( this%grid%if_report( n(this%tstep), p_node ) ) then
      call report_node_load_part_qedcyl( this%part, n(this%tstep), this%grid%ndump( p_node ), t(this%time), &
      this%no_co )
  endif

  ! grid load : number of particles per cell for all grid points
  if ( this%grid%if_report( n(this%tstep), p_grid ) ) then
      call report_grid_load( this%part, n(this%tstep), this%grid%ndump( p_grid ), t(this%time), &
      this%g_space, this%grid, this%no_co )
  endif

  call end_event(report_load_ev)

end subroutine report_load_sim_qedcyl
!-----------------------------------------------------------------------------------------

end module m_simulation_qedcyl
