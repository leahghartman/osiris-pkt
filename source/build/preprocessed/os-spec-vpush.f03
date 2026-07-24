# 1 "spec/os-spec-vpush.f03"
# 1 "<built-in>" 1
# 1 "<built-in>" 3
# 467 "<built-in>" 3
# 1 "<command line>" 1
# 1 "<built-in>" 2
# 1 "spec/os-spec-vpush.f03" 2
!-----------------------------------------------------------------------------------------
! Particle pusher using hardware (vector) acceleration
!-----------------------------------------------------------------------------------------

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
# 6 "spec/os-spec-vpush.f03" 2
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
# 7 "spec/os-spec-vpush.f03" 2

module m_species_vpush

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
# 11 "spec/os-spec-vpush.f03" 2

use m_parameters
use m_species_define
use m_emf_define
use m_vdf_define

implicit none

private

interface vadvance_deposit
    module procedure vadvance_deposit
end interface

public :: vadvance_deposit

contains

!---------------------------------------------------------------------------------------------------
! Interfaces to hardware optimized pushers
!---------------------------------------------------------------------------------------------------
# 438 "spec/os-spec-vpush.f03"
!---------------------------------------------------------------------------------------------------
! Push particles and deposit electric current using vector code
!---------------------------------------------------------------------------------------------------
subroutine vadvance_deposit( this, emf, jay, t, tstep, tid, n_threads )

    use m_time_step

    implicit none

    class( t_species ), intent(inout) :: this
    class( t_emf ), intent( in ) :: emf
    type( t_vdf ), dimension(:), intent(inout) :: jay

    real(p_double), intent(in) :: t
    type( t_time_step ) :: tstep

    integer, intent(in) :: tid ! local thread id
    integer, intent(in) :: n_threads ! total number of threads

    write(wrn_buf__,*) 'Vector code not available';call wrn__("spec/os-spec-vpush.f03",457)

end subroutine vadvance_deposit



end module m_species_vpush
