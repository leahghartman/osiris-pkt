# 1 "gr/os-emf-star-gr.f03"
# 1 "<built-in>" 1
# 1 "<built-in>" 3
# 467 "<built-in>" 3
# 1 "<command line>" 1
# 1 "<built-in>" 2
# 1 "gr/os-emf-star-gr.f03" 2
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
# 2 "gr/os-emf-star-gr.f03" 2
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
# 3 "gr/os-emf-star-gr.f03" 2

module m_emf_star_gr

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
# 7 "gr/os-emf-star-gr.f03" 2

use m_parameters
use m_geometry_gr

implicit none

private

interface update_b_2d_star
  module procedure update_b_2d_star
end interface

interface update_e_2d_star
  module procedure update_e_2d_star
end interface

public :: update_b_2d_star, update_e_2d_star

contains

!-------------------------------------------------------------------------------
subroutine update_b_2d_star( geometry, b_vdf )

  use m_vdf_define, only : t_vdf

  implicit none

  ! dummy variables
  class( t_geometry_gr ), intent(inout) :: geometry
  class ( t_vdf ), intent( inout ) :: b_vdf

  integer :: i2
  real( p_k_fld ), pointer :: r(:,:), st(:,:), ct(:,:)
  real( p_k_fld ), pointer :: b(:,:,:)
  real( p_k_fld ) :: rs, bs

  b => b_vdf%f2
  r => geometry%r%f1
  st => geometry%st%f1
  ct => geometry%ct%f1
  rs = geometry%rs
  bs = geometry%bs

  if ( geometry%metric == p_geometry_minkowski ) then
    select case (geometry%bprofile)
    case (p_bprofile_monopole)
      do i2 = 0, b_vdf%nx_(p_tc_dim) + 1
        b( 1, 1, i2 ) = bs

        ! inside the conductor:
        b( 1, 0, i2 ) = bs / r(1, 0)**2
        b( 2, 0, i2 ) = 0._p_double
        b( 3, 0, i2 ) = 0._p_double
      enddo
    case (p_bprofile_dipole)
      do i2 = 0, b_vdf%nx_(p_tc_dim) + 1
        b( 1, 1, i2 ) = bs * ct(2, i2)

        ! inside the conductor:
        b( 1, 0, i2 ) = bs * ct(2, i2) / r(1, 0)**3
        b( 2, 0, i2 ) = 0.5_p_double * bs * st(1, i2) / r(2, 0)**3
        b( 3, 0, i2 ) = 0._p_double
      enddo
    case default
    end select
  else
    do i2 = 0, b_vdf%nx_(p_tc_dim) + 1
      b( 1, 1, i2 ) = - 3._p_double * bs * (log(1._p_double - rs) +&
                       rs * (1._p_double + rs / 2._p_double)) * ct(2, i2) / rs**3

      ! inside the conductor:
      b( 1, 0, i2 ) = - 3._p_double * bs * (log(1._p_double - rs/r(1, 0)) +&
                       rs * (1._p_double + rs / (2._p_double*r(1, 0)))/r(1, 0)) * ct(2, i2) / rs**3
      !b( 2, 0, i2 ) = 1.5_p_double * bs * sqrt(1._p_double - rs) *&
      ! (log(1 - rs) + rs * (1._p_double + rs / 2._p_double)) * st(1, i2) / rs**3
      b( 2, 0, i2 ) = 1.5_p_double * bs * sqrt(1._p_double - rs/r(2, 0)) *&
                      (2._p_double * r(2, 0) * log(1 - rs/r(2, 0)) / rs + 1._p_double +&
                       1._p_double / (1._p_double - rs/r(2, 0))) * st(1, i2) / rs**2 / r(2, 0)
      b( 3, 0, i2 ) = 0._p_double
    enddo
  endif

  ! clear pointers
  nullify(b)
  nullify(r)
  nullify(st)
  nullify(ct)


end subroutine

!-------------------------------------------------------------------------------
subroutine update_e_2d_star( geometry, e_vdf )

  use m_vdf_define, only : t_vdf

  implicit none

  ! dummy variables
  class( t_geometry_gr ), intent(inout) :: geometry
  class ( t_vdf ), intent( inout ) :: e_vdf

  integer :: i2
  real( p_k_fld ), pointer :: r(:,:), st(:,:), ct(:,:)
  real( p_k_fld ), pointer :: e(:,:,:)
  real( p_k_fld ) :: rs, bs, omega, beta0

  e => e_vdf%f2
  r => geometry%r%f1
  st => geometry%st%f1
  ct => geometry%ct%f1
  rs = geometry%rs
  bs = geometry%bs
  omega = geometry%omega
  beta0 = geometry%beta0

  if ( geometry%metric == p_geometry_minkowski ) then
    select case (geometry%bprofile)
    case (p_bprofile_monopole)
      do i2 = 0, e_vdf%nx_(p_tc_dim) + 1
        e( 2, 1, i2 ) = - bs * omega * st(2, i2)
        e( 3, 1, i2 ) = 0._p_double

        ! inside the star:
        e( 1, 0, i2 ) = 0._p_double
        e( 2, 0, i2 ) = - bs * omega * st(2, i2) / r(1, 0)
        e( 3, 0, i2 ) = 0._p_double
      end do
    case (p_bprofile_dipole)
      do i2 = 0, e_vdf%nx_(p_tc_dim) + 1
        e( 2, 1, i2 ) = - bs * omega * st(2, i2) * ct(2, i2)
        e( 3, 1, i2 ) = 0._p_double

        ! inside the star:
        e( 1, 0, i2 ) = 0.5_p_double * bs * omega * st(1, i2)**2 / r(2, 0)**2
        e( 2, 0, i2 ) = - bs * omega * st(2, i2) * ct(2, i2) / r(1, 0)**2
        e( 3, 0, i2 ) = 0._p_double
      end do
    case default
    end select
  else
    do i2 = 0, e_vdf%nx_(p_tc_dim) + 1
      e( 2, 1, i2 ) = st(2, i2) * ct(2, i2) * bs / ( 1._p_double - rs)**(0.5) * &
       ( &
          3._p_double * omega / rs**3 * (log(1._p_double - rs) + rs * (1._p_double + rs/2._p_double)) - &
          90._p_double * beta0 / rs**5 * ( rs**2 / 30._p_double * (log(1._p_double - rs) + rs) + rs**4 / 60._p_double ) &
        )

      e( 3, 1, i2 ) = 0._p_double

      ! inside the star:
      !e( 1, 0, i2 ) = 3 * (omega-(beta0/r(2, 0)**3)) * r(2, 0) * st(1, i2)**2 * bs * &
      ! (log(1 - rs) + rs * (1 + rs/2.)) / rs**3
      !e( 2, 0, i2 ) = 3 * (omega-(beta0/r(1, 0)**3)) * r(1, 0) * st(2, i2) * ct(2, i2) * bs * &
      ! (log(1 - rs) + rs * (1 + rs/2.)) / rs**3 / sqrt(1 - rs/r(1, 0))

      !e( 1, 0, i2 ) = 0._p_double
      !e( 2, 0, i2 ) = 0._p_double
      !e( 3, 0, i2 ) = 0._p_double

      e( 1, 0, i2 ) = 1.5_p_double * bs * st(1, i2)**2 * (omega - beta0/r(2, 0)**3) *&
                     (2._p_double * r(2, 0) * log(1 - rs / r(2, 0)) / rs +&
                      1._p_double + 1._p_double / (1._p_double - rs/r(2, 0))) / rs**2
      e( 2, 0, i2 ) = 3._p_double * bs * st(2, i2) * ct(2, i2) * (omega - beta0/r(1, 0)**3) *&
                      (r(1, 0) * log(1._p_double - rs / r(1, 0)) + rs * (1._p_double + rs / (2._p_double*r(1, 0)))) /&
                      rs**3 / sqrt(1._p_double - rs/r(1, 0))
      e( 3, 0, i2 ) = 0._p_double

    end do
  endif

  ! clear pointers
  nullify(e)
  nullify(r)
  nullify(st)
  nullify(ct)

end subroutine

end module m_emf_star_gr
