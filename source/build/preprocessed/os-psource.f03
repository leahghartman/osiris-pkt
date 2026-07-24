# 1 "spec/psource/os-psource.f03"
# 1 "<built-in>" 1
# 1 "<built-in>" 3
# 467 "<built-in>" 3
# 1 "<command line>" 1
# 1 "<built-in>" 2
# 1 "spec/psource/os-psource.f03" 2
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
# 2 "spec/psource/os-psource.f03" 2
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
# 3 "spec/psource/os-psource.f03" 2

module m_profile

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
# 7 "spec/psource/os-psource.f03" 2

use m_psource_std
use m_psource_std_transposed
use m_psource_constq
use m_psource_beam
use m_psource_btfel
use m_psource_file
use m_psource_python

use m_species_define

implicit none

private

interface new_source
    module procedure new_source
end interface new_source

public :: new_source

contains

function new_source( type, input_file, coordinates )

  use m_input_file

  implicit none

  class( t_psource ), pointer :: new_source
  character(len=*), intent(in) :: type
  class( t_input_file ), intent(inout) :: input_file
  integer, intent(in) :: coordinates
  ! create t_profile object of the selected kind
  select case ( trim(type) )
  case('standard','profile')
    allocate( t_psource_std :: new_source )
  case('transposed')
    allocate( t_psource_std_transposed :: new_source )
  case('constq')
    allocate( t_psource_constq :: new_source )
  case('beamfocus')
    allocate( t_psource_beam :: new_source )
  case('betatronfel')
    allocate( t_psource_betatronfel :: new_source )
  case('file')
    allocate( t_psource_file :: new_source )
  case('python')
    allocate( t_psource_python :: new_source )
  case default
    new_source => null()
  end select

  ! If successfull read input parameters and finish initalization
  if ( associated(new_source) ) then
    call new_source % read_input( input_file, coordinates )
  endif

end function new_source

end module m_profile
