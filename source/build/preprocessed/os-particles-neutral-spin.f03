# 1 "neutral_spin/os-particles-neutral-spin.f03"
# 1 "<built-in>" 1
# 1 "<built-in>" 3
# 467 "<built-in>" 3
# 1 "<command line>" 1
# 1 "<built-in>" 2
# 1 "neutral_spin/os-particles-neutral-spin.f03" 2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Particles module for TDSE model
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
# 6 "neutral_spin/os-particles-neutral-spin.f03" 2
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
# 7 "neutral_spin/os-particles-neutral-spin.f03" 2

module m_particles_neutral_spin

use m_system
use m_parameters
use m_particles_define, only : t_particles
use m_species_define, only : t_species
use m_neutral_spin, only : t_neutral_spin

implicit none

private

type, extends( t_particles ) :: t_particles_neutral_spin

contains

  procedure :: allocate_objs => allocate_objs_part_neutral_spin
  ! procedure :: init => init_part_neutral_spin
  ! procedure :: advance_deposit => advance_deposit_part_neutral_spin

end type t_particles_neutral_spin

public :: t_particles_neutral_spin


contains

!-----------------------------------------------------------------------------------------
! Allocate particle species, cathodes, etc.
!-----------------------------------------------------------------------------------------
subroutine allocate_objs_part_neutral_spin( this )

  implicit none

  class( t_particles_neutral_spin ), intent(inout) :: this

  integer :: i
  class(t_species), pointer :: spec, tail


  if ( this%num_species > 0 ) then

    allocate( spec )
    this % species => spec
    do i = 2, this%num_species
      tail => spec
      allocate( spec )
      tail % next => spec
    enddo

    if ( this%num_cathode > 0 ) then
      allocate( this%cathode( this%num_cathode ) )
    endif


    if ( this%num_neutral > 0 ) then
      allocate( t_neutral_spin :: this%neutral( this%num_neutral ) )
    endif


  endif

end subroutine allocate_objs_part_neutral_spin
!-----------------------------------------------------------------------------------------

end module m_particles_neutral_spin
