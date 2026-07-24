# 1 "radio/os-radio-simulation.f03"
# 1 "<built-in>" 1
# 1 "<built-in>" 3
# 467 "<built-in>" 3
# 1 "<command line>" 1
# 1 "<built-in>" 2
# 1 "radio/os-radio-simulation.f03" 2
module m_simulation_rad

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
# 4 "radio/os-radio-simulation.f03" 2
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
# 5 "radio/os-radio-simulation.f03" 2

use m_simulation, only : t_simulation
use m_system

private

type, extends( t_simulation ) :: t_simulation_rad

contains

  procedure :: allocate_objs => allocate_objs_rad
  procedure :: list_algorithm => list_algorithm_rad

end type t_simulation_rad

public :: t_simulation_rad

contains

!-----------------------------------------------------------------------------------------
! Initialize simulation objects
!-----------------------------------------------------------------------------------------
subroutine allocate_objs_rad( sim )

  use m_particles_rad

  implicit none

  class( t_simulation_rad ), intent(inout) :: sim

  if (mpi_node()==0) print *,'Allocating Radiative particle objects.'

  allocate( t_particles_rad :: sim%part )

  ! allocate any remaining objects in superclass
  call sim % t_simulation % allocate_objs()

end subroutine allocate_objs_rad
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine list_algorithm_rad( sim )

  implicit none

  class( t_simulation_rad ), intent(in) :: sim

  print *, ""
  print *, "Simulation mode: RaDiO"
  call sim % t_simulation % list_algorithm()

end subroutine list_algorithm_rad
!-----------------------------------------------------------------------------------------

end module m_simulation_rad
