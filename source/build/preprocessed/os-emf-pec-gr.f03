# 1 "gr/os-emf-pec-gr.f03"
# 1 "<built-in>" 1
# 1 "<built-in>" 3
# 467 "<built-in>" 3
# 1 "<command line>" 1
# 1 "<built-in>" 2
# 1 "gr/os-emf-pec-gr.f03" 2
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
# 2 "gr/os-emf-pec-gr.f03" 2
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
# 3 "gr/os-emf-pec-gr.f03" 2

module m_emf_pec_gr

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
# 7 "gr/os-emf-pec-gr.f03" 2

use m_parameters
use m_geometry_gr

implicit none

private

interface update_b_2d_pec
  module procedure update_b_2d_pec
end interface

interface update_e_2d_pec
  module procedure update_e_2d_pec
end interface

public :: update_b_2d_pec, update_e_2d_pec

contains

!-------------------------------------------------------------------------------
subroutine update_b_2d_pec( geometry, b_vdf, bnd )

  use m_vdf_define, only : t_vdf

  implicit none

  ! dummy variables
  class( t_geometry_gr ), intent(inout) :: geometry
  class ( t_vdf ), intent( inout ) :: b_vdf
  integer, intent( in ) :: bnd

  integer :: i1, i2
  real( p_k_fld ), pointer :: b(:,:,:)

  b => b_vdf%f2

  select case (bnd)

  case (p_lower)
    ! inner radial boundary
    do i2 = 0, b_vdf%nx_(p_tc_dim) + 1
      ! br (at the surface):
      b( 1, 1, i2 ) = 0.

      ! inside the conductor:
      b( 1, 0, i2 ) = 0.
      b( 2, 0, i2 ) = 0.
      b( 3, 0, i2 ) = 0.
    enddo
  case (p_upper)
    ! outer radial boundary
    i1 = b_vdf%nx_(p_rc_dim)

    do i2 = 0, b_vdf%nx_(p_tc_dim)
      ! br (at the surface):
      b( 1, i1+1, i2) = 0.

      ! inside the conductor:
      b( 2, i1+1, i2) = 0.
      b( 3, i1+1, i2) = 0.
    enddo
  end select

  ! clear pointers
  nullify(b)

end subroutine

!-------------------------------------------------------------------------------
subroutine update_e_2d_pec( geometry, e_vdf, bnd )

  use m_vdf_define, only : t_vdf

  implicit none

  ! dummy variables
  class( t_geometry_gr ), intent(inout) :: geometry
  class ( t_vdf ), intent( inout ) :: e_vdf
  integer, intent( in ) :: bnd

  integer :: i1, i2
  real( p_k_fld ), pointer :: e(:,:,:)

  e => e_vdf%f2

  select case (bnd)

  case (p_lower)
    ! inner radial boundary
    do i2 = 0, e_vdf%nx_(p_tc_dim) + 1
      ! et (at the surface)
      e( 2, 1, i2 ) = 0.

      ! ep
      e( 3, 1, i2 ) = 0.

      ! inside the conductor:
      e( 1, 0, i2 ) = 0.
      e( 2, 0, i2 ) = 0.
      e( 3, 0, i2 ) = 0.
    enddo

  case (p_upper)
    ! outer radial boundary
    i1 = e_vdf%nx_(p_rc_dim)

    do i2 = 0, e_vdf%nx_(p_tc_dim) + 1
      ! et (at the surface):
      e( 2, i1+1, i2 ) = 0.

      ! ep
      e( 3, i1+1, i2 ) = 0.

      ! inside the conductor:
      e( 1, i1+1, i2 ) = 0.
    enddo
  end select

  ! clear pointers
  nullify(e)

end subroutine

end module m_emf_pec_gr
