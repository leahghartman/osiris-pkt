# 1 "gr/os-species-current-gr.f03"
# 1 "<built-in>" 1
# 1 "<built-in>" 3
# 467 "<built-in>" 3
# 1 "<command line>" 1
# 1 "<built-in>" 2
# 1 "gr/os-species-current-gr.f03" 2
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
# 2 "gr/os-species-current-gr.f03" 2
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
# 3 "gr/os-species-current-gr.f03" 2

module m_species_current_gr

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
# 7 "gr/os-species-current-gr.f03" 2

use m_system
use m_parameters

use m_geometry_gr
use m_math, only : pi_2

interface dep_current_2d_srect_mink_gr
   module procedure dep_current_2d_srect_mink_gr
end interface

interface dep_current_2d_sint_mink_gr
   module procedure dep_current_2d_sint_mink_gr
end interface

interface dep_current_2d_srect_curved_gr
   module procedure dep_current_2d_srect_curved_gr
end interface

interface dep_current_2d_sint_curved_gr
   module procedure dep_current_2d_sint_curved_gr
end interface

contains

!-------------------------------------------------------------------------------
! Charge conserving current deposit - Rectangular rule - Minkowski
!-------------------------------------------------------------------------------
subroutine dep_current_2d_srect_mink_gr( this, jay, dxi, xnew, ixold, xold, q, &
  rgamma, u, np, dt )

  use m_species_current
  use m_vdf_define, only : t_vdf

  use m_species_define_gr, only : t_species_gr

  implicit none

  ! dummy variables

  integer, parameter :: rank = 2
  class( t_species_gr ), intent(inout) :: this
  type( t_vdf ), intent(inout) :: jay

  integer, dimension(:,:), intent(inout) :: dxi, ixold
  real(p_k_part), dimension(:,:), intent(inout) :: xnew, xold ,u
  real(p_k_part), dimension( : ), intent(inout) :: q,rgamma

  integer, intent(in) :: np
  real(p_double), intent(in) :: dt

  ! local variables
  integer :: l, nsplit, ix, jx, k1, k2

  real(p_k_fld) :: x0, x1, y0, y1, vp
  real(p_k_fld) :: jnorm, qn

  type(t_vp2D), dimension( 3*p_cache_size ) :: vpbuf2D

  real(p_k_fld), dimension(-1:1) :: S0r, S1r
  real(p_k_fld), dimension(0:1) :: S0t, S1t, Sat

  real(p_k_fld), dimension(-1:1,0:1) :: divJr, divJt, Jr
  real(p_k_fld), dimension(-1:1) :: Jt
  real(p_k_fld) :: aux
  class( t_geometry_gr ), pointer :: g

  g => this%geometry

  ! executable statements

  ! split particles
  call split_2d_gr( dxi, xnew, ixold, xold, q, rgamma, u, np, vpbuf2D, nsplit )

  ! now accumulate jay looping through all virtual particles
  ! and shifting grid indexes
  ! flat case - rect. rule
  jnorm = real( 0.0795774715459477 * g%A / dt, p_k_fld ) ! 1/4pi

  do l=1,nsplit

    x0 = vpbuf2D(l)%x0
    y0 = vpbuf2D(l)%y0
    x1 = vpbuf2D(l)%x1
    y1 = vpbuf2D(l)%y1
    ix = vpbuf2D(l)%i
    jx = vpbuf2D(l)%j
    vp = vpbuf2D(l)%vz

    ! Normalize charge
    qn = real( vpbuf2D(l)%q, p_k_fld ) * jnorm

    ! get spline weights for r
    ! the spline_r is not modified by selecting dif. spline rules
    call g % spline_r_rectangular( ix, x0, S0r )
    call g % spline_r_rectangular( ix, x1, S1r )

    ! get spline weights for t
    call g % spline_t( jx, y0, S0t )
    call g % spline_t( jx, y1, S1t )
    call g % spline_t( jx, 0.5_p_k_fld*(y0+y1), Sat )

    ! get charge density time derivatives
    aux = 0.6666666666666667_p_double * g%vt%f1(1, jx) / g%st%f1(2, jx)

    do k1 = -1, 1
      do k2 = 0, 1
        divJr(k1,k2) = - qn * (S1r(k1) - S0r(k1)) * Sat(k2)
        divJt(k1,k2) = - qn * (S1r(k1) * S1t(k2) - S0r(k1) * S0t(k2)) - divJr(k1,k2)
      enddo

      ! get current polar components
      Jt(k1) = g%vr%f1(1,ix+k1) / (g%r%f1(2, ix+k1)**2 - g%r%f1(2, ix-1+k1)**2) * aux * divJt(k1,0)
    enddo

    ! get current radial components
    do k2 = 0, 1
      Jr(-1,k2) = 0.3333333333333333_p_double * g%vr%f1(1,ix-1) / g%r%f1(2, ix-1)**2 * divJr(-1,k2)
      Jr( 0,k2) = 0.3333333333333333_p_double * g%vr%f1(1,ix) / g%r%f1(2, ix)**2 * divJr( 0,k2) + g%r%f1(2, ix-1)**2/g%r%f1(2, ix)**2 * Jr(-1,k2)
    enddo

    ! deposit in-plane current components
    jay%f2(1,ix-1,jx ) = jay%f2(1,ix-1,jx ) + Jr(-1,0)
    jay%f2(2,ix-1,jx ) = jay%f2(2,ix-1,jx ) + Jt(-1)
    jay%f2(1,ix-1,jx+1) = jay%f2(1,ix-1,jx+1) + Jr(-1,1)
    jay%f2(1,ix ,jx ) = jay%f2(1,ix ,jx ) + Jr( 0,0)
    jay%f2(2,ix ,jx ) = jay%f2(2,ix ,jx ) + Jt( 0)
    jay%f2(1,ix ,jx+1) = jay%f2(1,ix ,jx+1) + Jr( 0,1)
    jay%f2(2,ix+1,jx ) = jay%f2(2,ix+1,jx ) + Jt( 1)

    ! deposit out-of-plane current component
    do k1 = -1, 1
      do k2 = 0, 1
        jay%f2(3,ix+k1,jx+k2) = jay%f2(3,ix+k1,jx+k2) + vp * qn * 0.5_p_double * dt * (S1r(k1) * S1t(k2) + S0r(k1) * S0t(k2))
      enddo
    enddo

  enddo

end subroutine dep_current_2d_srect_mink_gr
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Charge conserving current deposit - Integral rule - Minkowski
!-------------------------------------------------------------------------------
subroutine dep_current_2d_sint_mink_gr( this, jay, dxi, xnew, ixold, xold, q, &
  rgamma, u, np, dt )

  use m_species_current
  use m_vdf_define, only : t_vdf

  use m_species_define_gr, only : t_species_gr

  implicit none

  ! dummy variables

  integer, parameter :: rank = 2
  class( t_species_gr ), intent(inout) :: this
  type( t_vdf ), intent(inout) :: jay

  integer, dimension(:,:), intent(inout) :: dxi, ixold
  real(p_k_part), dimension(:,:), intent(inout) :: xnew, xold ,u
  real(p_k_part), dimension( : ), intent(inout) :: q,rgamma

  integer, intent(in) :: np
  real(p_double), intent(in) :: dt

  ! local variables
  integer :: l, nsplit, ix, jx, k1, k2

  real(p_k_fld) :: x0, x1, y0, y1, vp
  real(p_k_fld) :: jnorm, qn

  type(t_vp2D), dimension( 3*p_cache_size ) :: vpbuf2D

  real(p_k_fld), dimension(-1:1) :: S0r, S1r
  real(p_k_fld), dimension(0:1) :: S0t, S1t, Sat

  real(p_k_fld), dimension(-1:1,0:1) :: divJr, divJt, Jr
  real(p_k_fld), dimension(-1:1) :: Jt
  real(p_k_fld) :: aux
  class( t_geometry_gr ), pointer :: g

  g => this%geometry

  ! executable statements

  ! split particles
  call split_2d_gr( dxi, xnew, ixold, xold, q, rgamma, u, np, vpbuf2D, nsplit )

  ! now accumulate jay looping through all virtual particles
  ! and shifting grid indexes
  ! flat case - integral rule
  jnorm = real( 0.0596831036594608 / dt * (1._p_double + (1._p_double / g%delta))**3 / (1._p_double - (1._p_double / g%delta**3) ), p_k_fld ) ! 3/16pi

  do l=1,nsplit

    x0 = vpbuf2D(l)%x0
    y0 = vpbuf2D(l)%y0
    x1 = vpbuf2D(l)%x1
    y1 = vpbuf2D(l)%y1
    ix = vpbuf2D(l)%i
    jx = vpbuf2D(l)%j
    vp = vpbuf2D(l)%vz

    ! Normalize charge
    qn = real( vpbuf2D(l)%q, p_k_fld ) * jnorm

    ! get spline weights for r
    ! the spline_r is not modified by selecting dif. spline rules
    call g % spline_r_integral( ix, x0, S0r )
    call g % spline_r_integral( ix, x1, S1r )

    ! get spline weights for t
    call g % spline_t( jx, y0, S0t )
    call g % spline_t( jx, y1, S1t )
    call g % spline_t( jx, 0.5_p_k_fld*(y0+y1), Sat )

    ! get charge density time derivatives
    aux = 0.6666666666666667_p_double * g%vt%f1(1, jx) / g%st%f1(2, jx)

    do k1 = -1, 1
      do k2 = 0, 1
        divJr(k1,k2) = - qn * (S1r(k1) - S0r(k1)) * Sat(k2)
        divJt(k1,k2) = - qn * (S1r(k1) * S1t(k2) - S0r(k1) * S0t(k2)) - divJr(k1,k2)
      enddo

      ! get current polar components
      Jt(k1) = g%vr%f1(1,ix+k1) / (g%r%f1(2, ix+k1)**2 - g%r%f1(2, ix-1+k1)**2) * aux * divJt(k1,0)
    enddo

    ! get current radial components
    do k2 = 0, 1
      Jr(-1,k2) = 0.3333333333333333_p_double * g%vr%f1(1,ix-1) / g%r%f1(2, ix-1)**2 * divJr(-1,k2)
      Jr( 0,k2) = 0.3333333333333333_p_double * g%vr%f1(1,ix) / g%r%f1(2, ix)**2 * divJr( 0,k2) + g%r%f1(2, ix-1)**2/g%r%f1(2, ix)**2 * Jr(-1,k2)
    enddo

    ! deposit in-plane current components
    jay%f2(1,ix-1,jx ) = jay%f2(1,ix-1,jx ) + Jr(-1,0)
    jay%f2(2,ix-1,jx ) = jay%f2(2,ix-1,jx ) + Jt(-1)
    jay%f2(1,ix-1,jx+1) = jay%f2(1,ix-1,jx+1) + Jr(-1,1)
    jay%f2(1,ix ,jx ) = jay%f2(1,ix ,jx ) + Jr( 0,0)
    jay%f2(2,ix ,jx ) = jay%f2(2,ix ,jx ) + Jt( 0)
    jay%f2(1,ix ,jx+1) = jay%f2(1,ix ,jx+1) + Jr( 0,1)
    jay%f2(2,ix+1,jx ) = jay%f2(2,ix+1,jx ) + Jt( 1)

    ! deposit out-of-plane current component
    do k1 = -1, 1
      do k2 = 0, 1
        jay%f2(3,ix+k1,jx+k2) = jay%f2(3,ix+k1,jx+k2) + vp * qn * 0.5_p_double * dt * (S1r(k1) * S1t(k2) + S0r(k1) * S0t(k2))
      enddo
    enddo

  enddo

end subroutine dep_current_2d_sint_mink_gr
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Charge conserving current deposit - Rectangular rule - Schwarzschild
!-------------------------------------------------------------------------------
subroutine dep_current_2d_srect_curved_gr(this, jay, dxi, xnew, xposnew, ixold, xold, xposold, q, &
  rgamma, u, uold, np, dt )

  use m_species_current
  use m_vdf_define, only : t_vdf

  use m_species_define_gr, only : t_species_gr

  implicit none

  ! dummy variables

  integer, parameter :: rank = 2
  class( t_species_gr ), intent(inout) :: this
  type( t_vdf ), intent(inout) :: jay

  integer, dimension(:,:), intent(inout) :: dxi, ixold
  real(p_k_part), dimension(:,:), intent(inout) :: xnew, xposnew, xposold, xold ,u, uold
  real(p_k_part), dimension( : ), intent(inout) :: q,rgamma

  integer, intent(in) :: np
  real(p_double), intent(in) :: dt

  ! local variables
  integer :: l, nsplit, ix, jx, k1, k2

  real(p_k_fld) :: x0, x1, y0, y1, vp
  real(p_k_fld) :: jnorm, qn

  type(t_vp2D), dimension( 3*p_cache_size ) :: vpbuf2D

  real(p_k_fld), dimension(-1:1) :: S0r, S1r
  real(p_k_fld), dimension(0:1) :: S0t, S1t, Sat

  real(p_k_fld), dimension(-1:1,0:1) :: divJr, divJt, Jr
  real(p_k_fld), dimension(-1:1) :: Jt
  class( t_geometry_gr ), pointer :: g

  g => this%geometry

  ! executable statements

  ! split particles
  call split_2d_curved_gr( this, dxi, xnew, xposnew, ixold, xold, xposold, q, rgamma, u, uold, np, vpbuf2D, nsplit )

  ! now accumulate jay looping through all virtual particles
  ! and shifting grid indexes
  ! curved case - rect. rule
  jnorm = real( 0.0795774715459477 * g%A / dt, p_k_fld ) ! 1/4pi

  do l=1,nsplit

    x0 = vpbuf2D(l)%x0
    y0 = vpbuf2D(l)%y0
    x1 = vpbuf2D(l)%x1
    y1 = vpbuf2D(l)%y1
    ix = vpbuf2D(l)%i
    jx = vpbuf2D(l)%j
    vp = vpbuf2D(l)%vz

    ! Normalize charge
    qn = real( vpbuf2D(l)%q, p_k_fld ) * jnorm

    ! get spline weights for r
    call g % spline_r_rectangular_gr( ix, x0, S0r )
    call g % spline_r_rectangular_gr( ix, x1, S1r )


    ! get spline weights for t
    call g % spline_t( jx, y0, S0t )
    call g % spline_t( jx, y1, S1t )
    call g % spline_t( jx, 0.5_p_k_fld*(y0+y1), Sat )

    ! get charge density time derivatives
    do k1 = -1, 1
      do k2 = 0, 1
        divJr(k1,k2) = - qn * (S1r(k1) - S0r(k1)) * Sat(k2)
        divJt(k1,k2) = - qn * (S1r(k1) * S1t(k2) - S0r(k1) * S0t(k2)) - divJr(k1,k2)
      enddo

      ! get current polar components
      Jt(k1) = divJt(k1,0) * g%vr%f1(1,ix+k1) * g%vt%f1(1, jx) / g%at%f2(2, ix+k1, jx)
    enddo

    ! get current radial components
    do k2 = 0, 1
      Jr(-1,k2) = divJr(-1,k2) * g%vr%f1(1,ix-1) * g%vt%f1(1,jx) / g%ar%f2(1, ix-1, jx)
      Jr( 0,k2) = divJr( 0,k2) * g%vr%f1(1,ix) * g%vt%f1(1, jx) / g%ar%f2(1, ix, jx) + g%ar%f2(1, ix-1, jx)/g%ar%f2(1, ix, jx) * Jr(-1,k2)

    enddo

    ! deposit in-plane current components
    jay%f2(1,ix-1,jx ) = jay%f2(1,ix-1,jx ) + Jr(-1,0)
    jay%f2(2,ix-1,jx ) = jay%f2(2,ix-1,jx ) + Jt(-1)
    jay%f2(1,ix-1,jx+1) = jay%f2(1,ix-1,jx+1) + Jr(-1,1)
    jay%f2(1,ix ,jx ) = jay%f2(1,ix ,jx ) + Jr( 0,0)
    jay%f2(2,ix ,jx ) = jay%f2(2,ix ,jx ) + Jt( 0)
    jay%f2(1,ix ,jx+1) = jay%f2(1,ix ,jx+1) + Jr( 0,1)
    jay%f2(2,ix+1,jx ) = jay%f2(2,ix+1,jx ) + Jt( 1)

    !! deposit out-of-plane current component
    do k1 = -1, 1
      do k2 = 0, 1
                                                                    ! this dt appears to cancel the one in qn
                                                        ! alpha * jay - rho * beta
                                                                        ! this is rho
                                                        ! this vp is alpha*w_phi/gamma - beta_phi
        jay%f2(3,ix+k1,jx+k2) = jay%f2(3,ix+k1,jx+k2) + (vp) * qn * dt * 0.5_p_double * (S1r(k1) * S1t(k2) + S0r(k1) * S0t(k2))
      enddo
    enddo

  enddo

end subroutine dep_current_2d_srect_curved_gr
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Charge conserving current deposit - Integral rule - Schwarzschild
!-------------------------------------------------------------------------------
subroutine dep_current_2d_sint_curved_gr( this, jay, dxi, xnew, xposnew, ixold, xold, xposold, q, &
  rgamma, u, uold, np, dt )

  use m_species_current
  use m_vdf_define, only : t_vdf

  use m_species_define_gr, only : t_species_gr

  implicit none

  ! dummy variables

  integer, parameter :: rank = 2
  class( t_species_gr ), intent(inout) :: this
  type( t_vdf ), intent(inout) :: jay

  integer, dimension(:,:), intent(inout) :: dxi, ixold
  real(p_k_part), dimension(:,:), intent(inout) :: xnew, xposnew, xposold, xold ,u, uold
  real(p_k_part), dimension( : ), intent(inout) :: q,rgamma

  integer, intent(in) :: np
  real(p_double), intent(in) :: dt

  ! local variables
  integer :: l, nsplit, ix, jx, k1, k2

  real(p_k_fld) :: x0, x1, y0, y1, vp
  real(p_k_fld) :: jnorm, qn

  type(t_vp2D), dimension( 3*p_cache_size ) :: vpbuf2D

  real(p_k_fld), dimension(-1:1) :: S0r, S1r
  real(p_k_fld), dimension(0:1) :: S0t, S1t, Sat

  real(p_k_fld), dimension(-1:1,0:1) :: divJr, divJt, Jr
  real(p_k_fld), dimension(-1:1) :: Jt
  class( t_geometry_gr ), pointer :: g

  g => this%geometry

  ! executable statements
  ! split particles
  call split_2d_curved_gr( this, dxi, xnew, xposnew, ixold, xold, xposold, q, rgamma, u, uold, np, vpbuf2D, nsplit )

  ! now accumulate jay looping through all virtual particles
  ! and shifting grid indexes
  ! curved case - integral rule
  jnorm = real( 0.1591549430918953 / dt , p_k_fld ) ! 1/2pi

  do l=1,nsplit

    x0 = vpbuf2D(l)%x0
    y0 = vpbuf2D(l)%y0
    x1 = vpbuf2D(l)%x1
    y1 = vpbuf2D(l)%y1
    ix = vpbuf2D(l)%i
    jx = vpbuf2D(l)%j
    vp = vpbuf2D(l)%vz

    ! Normalize charge
    qn = real( vpbuf2D(l)%q, p_k_fld ) * jnorm

    ! get spline weights for r
    call g % spline_r_integral_gr( ix, x0, S0r )
    call g % spline_r_integral_gr( ix, x1, S1r )


    ! get spline weights for t
    call g % spline_t( jx, y0, S0t )
    call g % spline_t( jx, y1, S1t )
    call g % spline_t( jx, 0.5_p_k_fld*(y0+y1), Sat )

    ! get charge density time derivatives
    do k1 = -1, 1
      do k2 = 0, 1
        divJr(k1,k2) = - qn * (S1r(k1) - S0r(k1)) * Sat(k2)
        divJt(k1,k2) = - qn * (S1r(k1) * S1t(k2) - S0r(k1) * S0t(k2)) - divJr(k1,k2)
      enddo

      ! get current polar components
      Jt(k1) = divJt(k1,0) * g%vr%f1(1,ix+k1) * g%vt%f1(1, jx) / g%at%f2(2, ix+k1, jx)
    enddo

    ! get current radial components
    do k2 = 0, 1
      Jr(-1,k2) = divJr(-1,k2) * g%vr%f1(1,ix-1) * g%vt%f1(1,jx) / g%ar%f2(1, ix-1, jx)
      Jr( 0,k2) = divJr( 0,k2) * g%vr%f1(1,ix) * g%vt%f1(1, jx) / g%ar%f2(1, ix, jx) + g%ar%f2(1, ix-1, jx)/g%ar%f2(1, ix, jx) * Jr(-1,k2)

    enddo

    ! deposit in-plane current components
    jay%f2(1,ix-1,jx ) = jay%f2(1,ix-1,jx ) + Jr(-1,0)
    jay%f2(2,ix-1,jx ) = jay%f2(2,ix-1,jx ) + Jt(-1)
    jay%f2(1,ix-1,jx+1) = jay%f2(1,ix-1,jx+1) + Jr(-1,1)
    jay%f2(1,ix ,jx ) = jay%f2(1,ix ,jx ) + Jr( 0,0)
    jay%f2(2,ix ,jx ) = jay%f2(2,ix ,jx ) + Jt( 0)
    jay%f2(1,ix ,jx+1) = jay%f2(1,ix ,jx+1) + Jr( 0,1)
    jay%f2(2,ix+1,jx ) = jay%f2(2,ix+1,jx ) + Jt( 1)

    ! deposit out-of-plane current component
    do k1 = -1, 1
      do k2 = 0, 1
                                                                    ! this dt appears to cancel the one in qn
                                                        ! alpha * jay - rho * beta
                                                                        ! this is rho
                                                        ! this vp is alpha*w_phi/gamma - beta_phi
        jay%f2(3,ix+k1,jx+k2) = jay%f2(3,ix+k1,jx+k2) + (vp) * qn * dt * 0.5_p_double * (S1r(k1) * S1t(k2) + S0r(k1) * S0t(k2))
      enddo
    enddo

  enddo

end subroutine dep_current_2d_sint_curved_gr
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine split_2d_gr( dxi, xnew, ixold, xold, q, rgamma, u, np, &
  vpbuf2D, nsplit )

  use m_species_current, only : t_vp2D

  implicit none

  integer, dimension(:,:), intent(in) :: dxi, ixold
  real(p_k_part), dimension(:,:), intent(in) :: xnew, xold ,u
  real(p_k_part), dimension( : ), intent(in) :: q,rgamma

  integer, intent(in) :: np ! number of particles to process

  type(t_vp2D), dimension(:), intent(out) :: vpbuf2D
  integer, intent(out) :: nsplit ! number of virtual particles created

  ! local variables
  real(p_k_part), dimension( np ) :: vptemp
  real(p_k_part) :: xint, yint, xint2, yint2, delta
  real(p_k_part) :: vz, vzint
  real(p_k_part) :: xmid, ymid, rmid
  integer :: k,l

  integer :: cross

  k=0

  ! create virtual particles that correspond to a motion that starts and ends
  ! in the same grid cell
  do l=1,np
    ! positions in integer times and velocities in semi
    xmid = 0.5_p_k_fld * (xold(3,l) + xnew(3,l))
    ymid = 0.5_p_k_fld * (xold(4,l) + xnew(4,l))
    rmid = sqrt(xmid**2 + ymid**2)

    ! vphi must be evaluated in semi
    vptemp(l) = ( -ymid * u(1,l) + xmid * u(2,l) ) * rgamma(l) / rmid
  enddo

  do l=1,np

    cross = abs(dxi(1,l)) + 2 * abs(dxi(2,l))

    select case( cross )
      case(0) ! no cross
          k=k+1
          vpbuf2D(k)%x0 = xold(1,l)
          vpbuf2D(k)%y0 = xold(2,l)
          vpbuf2D(k)%x1 = xnew(1,l)
          vpbuf2D(k)%y1 = xnew(2,l)
          vpbuf2D(k)%q = q(l)

          vpbuf2D(k)%vz = vptemp(l)

          vpbuf2D(k)%i = ixold(1,l)
          vpbuf2D(k)%j = ixold(2,l)


      case(1) ! x cross only

          xint = 0.5_p_k_fld * dxi(1,l)
          delta = ( xint - xold(1,l) ) / ( xnew(1,l) - xold(1,l) )
          yint = xold(2,l) + (xnew(2,l) - xold(2,l)) * delta

          vz = vptemp(l)

          k=k+1

          vpbuf2D(k)%x0 = xold(1,l)
          vpbuf2D(k)%y0 = xold(2,l)

          vpbuf2D(k)%x1 = xint
          vpbuf2D(k)%y1 = yint

          vpbuf2D(k)%q = q(l)
          vpbuf2D(k)%vz = vz * delta

          vpbuf2D(k)%i = ixold(1,l)
          vpbuf2D(k)%j = ixold(2,l)

          k=k+1
          vpbuf2D(k)%x0 = -xint
          vpbuf2D(k)%y0 = yint

          vpbuf2D(k)%x1 = xnew(1,l)
          vpbuf2D(k)%y1 = xnew(2,l)

          vpbuf2D(k)%q = q(l)
          vpbuf2D(k)%vz = vz * (1._p_double-delta)

          vpbuf2D(k)%i = ixold(1,l) + dxi(1,l)
          vpbuf2D(k)%j = ixold(2,l)

      case(2) ! y cross only

          yint = 0.5_p_k_fld * dxi(2,l)
          delta = ( yint - xold(2,l) ) / ( xnew(2,l) - xold(2,l))
          xint = xold(1,l) + (xnew(1,l) - xold(1,l)) * delta

          vz = vptemp(l)

          k=k+1
          vpbuf2D(k)%x0 = xold(1,l)
          vpbuf2D(k)%y0 = xold(2,l)

          vpbuf2D(k)%x1 = xint
          vpbuf2D(k)%y1 = yint
          vpbuf2D(k)%q = q(l)
          vpbuf2D(k)%vz = vz * delta

          vpbuf2D(k)%i = ixold(1,l)
          vpbuf2D(k)%j = ixold(2,l)

          k=k+1
          vpbuf2D(k)%x0 = xint
          vpbuf2D(k)%y0 = -yint

          vpbuf2D(k)%x1 = xnew(1,l)
          vpbuf2D(k)%y1 = xnew(2,l)
          vpbuf2D(k)%q = q(l)
          vpbuf2D(k)%vz = vz * (1._p_double-delta)

          vpbuf2D(k)%i = ixold(1,l)
          vpbuf2D(k)%j = ixold(2,l) + dxi(2,l)

      case(3) ! x,y cross

          ! split in x direction first
          xint = 0.5_p_k_fld * dxi(1,l)
          delta = ( xint - xold(1,l) ) / ( xnew(1,l) - xold(1,l))
          yint = xold(2,l) + ( xnew(2,l) - xold(2,l)) * delta

          vz = vptemp(l)

          ! check if y intersection occured for 1st or 2nd split
          if ((yint >= -0.5_p_k_fld) .and. ( yint < 0.5_p_k_fld )) then

            ! no y cross on 1st vp
            k=k+1
            vpbuf2D(k)%x0 = xold(1,l)
            vpbuf2D(k)%y0 = xold(2,l)

            vpbuf2D(k)%x1 = xint
            vpbuf2D(k)%y1 = yint

            vpbuf2D(k)%q = q(l)
            vpbuf2D(k)%vz = vz * delta

            vzint = vz*(1-delta)

            vpbuf2D(k)%i = ixold(1,l)
            vpbuf2D(k)%j = ixold(2,l)

            ! y split 2nd vp
            k=k+1

            yint2 = 0.5_p_k_fld * dxi(2,l)
            delta = ( yint2 - yint ) / ( xnew(2,l) - yint )
            xint2 = -xint + ( xnew(1,l) - xint ) * delta


            vpbuf2D(k)%x0 = -xint
            vpbuf2D(k)%y0 = yint

            vpbuf2D(k)%x1 = xint2
            vpbuf2D(k)%y1 = yint2

            vpbuf2D(k)%q = q(l)
            vpbuf2D(k)%vz = vzint * delta

            vpbuf2D(k)%i = ixold(1,l) + dxi(1,l)
            vpbuf2D(k)%j = ixold(2,l)

            k=k+1

            vpbuf2D(k)%x0 = xint2
            vpbuf2D(k)%y0 = -yint2

            vpbuf2D(k)%x1 = xnew(1,l)
            vpbuf2D(k)%y1 = xnew(2,l)
            vpbuf2D(k)%q = q(l)
            vpbuf2D(k)%vz = vzint * (1._p_double-delta)

            vpbuf2D(k)%i = ixold(1,l) + dxi(1,l)
            vpbuf2D(k)%j = ixold(2,l) + dxi(2,l)

          else

            vzint = vz * delta

            ! y split 1st vp
            yint2 = 0.5_p_k_fld * dxi(2,l)
            delta = ( yint2 - xold(2,l) ) / ( yint - xold(2,l))
            xint2 = xold(1,l) + (xint - xold(1,l)) * delta

            k=k+1
            vpbuf2D(k)%x0 = xold(1,l)
            vpbuf2D(k)%y0 = xold(2,l)

            vpbuf2D(k)%x1 = xint2
            vpbuf2D(k)%y1 = yint2

            vpbuf2D(k)%q = q(l)
            vpbuf2D(k)%vz = vzint*delta

            vpbuf2D(k)%i = ixold(1,l)
            vpbuf2D(k)%j = ixold(2,l)


            k=k+1
            vpbuf2D(k)%x0 = xint2
            vpbuf2D(k)%y0 = -yint2

            vpbuf2D(k)%x1 = xint
            vpbuf2D(k)%y1 = yint - dxi(2,l)

            vpbuf2D(k)%q = q(l)
            vpbuf2D(k)%vz = vzint*(1._p_double-delta)

            vpbuf2D(k)%i = ixold(1,l)
            vpbuf2D(k)%j = ixold(2,l) + dxi(2,l)

            ! no y cross on second vp
            k=k+1
            vpbuf2D(k)%x0 = - xint
            vpbuf2D(k)%y0 = yint - dxi(2,l)
            vpbuf2D(k)%x1 = xnew(1,l)
            vpbuf2D(k)%y1 = xnew(2,l)
            vpbuf2D(k)%q = q(l)
            vpbuf2D(k)%vz = vz - vzint

            vpbuf2D(k)%i = ixold(1,l) + dxi(1,l)
            vpbuf2D(k)%j = ixold(2,l) + dxi(2,l)

          endif

    end select

  enddo

  nsplit = k

end subroutine split_2d_gr
!-----------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine split_2d_curved_gr( this, dxi, xnew, xposnew, ixold, xold, xposold, q, rgamma, u, uold, np, &
  vpbuf2D, nsplit )

  use m_species_current, only : t_vp2D
  use m_species_define_gr, only : t_species_gr

  implicit none

  class( t_species_gr ), intent(inout) :: this
  integer, dimension(:,:), intent(in) :: dxi, ixold
  real(p_k_part), dimension(:,:), intent(in) :: xnew, xold ,u, uold
  real(p_k_part), dimension( : ), intent(in) :: q,rgamma
  real(p_k_part), dimension(:,:), intent(in) :: xposnew, xposold

  integer, intent(in) :: np ! number of particles to process

  type(t_vp2D), dimension(:), intent(out) :: vpbuf2D
  integer, intent(out) :: nsplit ! number of virtual particles created

  ! local variables
  real(p_k_part), dimension( np ) :: vptemp
  real(p_k_part) :: xint, yint, xint2, yint2, delta
  real(p_k_part) :: vz, vzint
  real(p_k_part) :: alpgold, alpgnew, betamed
  integer :: k,l
  real(p_double) :: shift, tc, tc2, tc4, tc6, tc8, tc10
  real(p_double) :: stold, stnew

  integer :: cross
  class( t_geometry_gr ), pointer :: g
  g => this%geometry


  k=0

  ! used for the polynomial approx of trig functions
  shift = pi_2

  ! create virtual particles that correspond to a motion that starts and ends
  ! in the same grid cell
  do l=1,np
    ! positions and velocities in integer times
    ! vphi must be evaluated in semi-integer times

    ! evaluate the old sin of theta around pi/2
    tc = xposold(2,l) - shift
    tc2 = tc * tc
    tc4 = tc2 * tc2
    tc6 = tc4 * tc2
    tc8 = tc6 * tc2
    tc10 = tc8 * tc2

    stold = 1.0_p_double - 0.500000000000000_p_double * tc2 +&
                           0.041666666666667_p_double * tc4 - 0.001388888888889_p_double * tc6 +&
                           0.000024801587302_p_double * tc8 - 0.000000275573192_p_double * tc10

    ! evaluate the new sin of theta around pi/2
    tc = xposnew(2,l) - shift
    tc2 = tc * tc
    tc4 = tc2 * tc2
    tc6 = tc4 * tc2
    tc8 = tc6 * tc2
    tc10 = tc8 * tc2

    stnew = 1.0_p_double - 0.500000000000000_p_double * tc2 +&
                           0.041666666666667_p_double * tc4 - 0.001388888888889_p_double * tc6 +&
                           0.000024801587302_p_double * tc8 - 0.000000275573192_p_double * tc10

    ! time averaged shift vector
    betamed = - g%beta0 * 0.5_p_double * (stold / xposold(1,l)**2 + stnew / xposnew(1,l)**2)

    ! lapse function divided by the particle's gamma for the old and new timesteps
    alpgold = sqrt(1.0_p_double - g%rs/xposold(1,l)) / sqrt(1 + uold(1,l)**2 + uold(2,l)**2 + uold(3,l)**2)
    alpgnew = sqrt(1.0_p_double - g%rs/xposnew(1,l)) / sqrt(1 + u(1,l)**2 + u(2,l)**2 + u(3,l)**2)

    ! time averaged orthonormalized phi velocity
    vptemp(l) = 0.5_p_double * (uold(3,l)*alpgold + u(3,l)*alpgnew) - betamed

  enddo

  do l=1,np

    cross = abs(dxi(1,l)) + 2 * abs(dxi(2,l))

    select case( cross )
      case(0) ! no cross
          k=k+1
          vpbuf2D(k)%x0 = xold(1,l)
          vpbuf2D(k)%y0 = xold(2,l)
          vpbuf2D(k)%x1 = xnew(1,l)
          vpbuf2D(k)%y1 = xnew(2,l)
          vpbuf2D(k)%q = q(l)

          vpbuf2D(k)%vz = vptemp(l)

          vpbuf2D(k)%i = ixold(1,l)
          vpbuf2D(k)%j = ixold(2,l)

      case(1) ! x cross only

          xint = 0.5_p_k_fld * dxi(1,l)
          delta = ( xint - xold(1,l) ) / ( xnew(1,l) - xold(1,l) )
          yint = xold(2,l) + (xnew(2,l) - xold(2,l)) * delta

          vz = vptemp(l)

          k=k+1

          vpbuf2D(k)%x0 = xold(1,l)
          vpbuf2D(k)%y0 = xold(2,l)

          vpbuf2D(k)%x1 = xint
          vpbuf2D(k)%y1 = yint

          vpbuf2D(k)%q = q(l)
          vpbuf2D(k)%vz = vz * delta

          vpbuf2D(k)%i = ixold(1,l)
          vpbuf2D(k)%j = ixold(2,l)

          k=k+1
          vpbuf2D(k)%x0 = -xint
          vpbuf2D(k)%y0 = yint

          vpbuf2D(k)%x1 = xnew(1,l)
          vpbuf2D(k)%y1 = xnew(2,l)

          vpbuf2D(k)%q = q(l)
          vpbuf2D(k)%vz = vz * (1._p_double-delta)

          vpbuf2D(k)%i = ixold(1,l) + dxi(1,l)
          vpbuf2D(k)%j = ixold(2,l)

      case(2) ! y cross only

          yint = 0.5_p_k_fld * dxi(2,l)
          delta = ( yint - xold(2,l) ) / ( xnew(2,l) - xold(2,l))
          xint = xold(1,l) + (xnew(1,l) - xold(1,l)) * delta

          vz = vptemp(l)

          k=k+1
          vpbuf2D(k)%x0 = xold(1,l)
          vpbuf2D(k)%y0 = xold(2,l)

          vpbuf2D(k)%x1 = xint
          vpbuf2D(k)%y1 = yint
          vpbuf2D(k)%q = q(l)
          vpbuf2D(k)%vz = vz * delta

          vpbuf2D(k)%i = ixold(1,l)
          vpbuf2D(k)%j = ixold(2,l)

          k=k+1
          vpbuf2D(k)%x0 = xint
          vpbuf2D(k)%y0 = -yint

          vpbuf2D(k)%x1 = xnew(1,l)
          vpbuf2D(k)%y1 = xnew(2,l)
          vpbuf2D(k)%q = q(l)
          vpbuf2D(k)%vz = vz * (1._p_double-delta)

          vpbuf2D(k)%i = ixold(1,l)
          vpbuf2D(k)%j = ixold(2,l) + dxi(2,l)

      case(3) ! x,y cross

          ! split in x direction first
          xint = 0.5_p_k_fld * dxi(1,l)
          delta = ( xint - xold(1,l) ) / ( xnew(1,l) - xold(1,l))
          yint = xold(2,l) + ( xnew(2,l) - xold(2,l)) * delta

          vz = vptemp(l)

          ! check if y intersection occured for 1st or 2nd split
          if ((yint >= -0.5_p_k_fld) .and. ( yint < 0.5_p_k_fld )) then

            ! no y cross on 1st vp
            k=k+1
            vpbuf2D(k)%x0 = xold(1,l)
            vpbuf2D(k)%y0 = xold(2,l)

            vpbuf2D(k)%x1 = xint
            vpbuf2D(k)%y1 = yint

            vpbuf2D(k)%q = q(l)
            vpbuf2D(k)%vz = vz * delta

            vzint = vz*(1._p_double-delta)

            vpbuf2D(k)%i = ixold(1,l)
            vpbuf2D(k)%j = ixold(2,l)

            ! y split 2nd vp
            k=k+1

            yint2 = 0.5_p_k_fld * dxi(2,l)
            delta = ( yint2 - yint ) / ( xnew(2,l) - yint )
            xint2 = -xint + ( xnew(1,l) - xint ) * delta


            vpbuf2D(k)%x0 = -xint
            vpbuf2D(k)%y0 = yint

            vpbuf2D(k)%x1 = xint2
            vpbuf2D(k)%y1 = yint2

            vpbuf2D(k)%q = q(l)
            vpbuf2D(k)%vz = vzint * delta

            vpbuf2D(k)%i = ixold(1,l) + dxi(1,l)
            vpbuf2D(k)%j = ixold(2,l)

            k=k+1

            vpbuf2D(k)%x0 = xint2
            vpbuf2D(k)%y0 = -yint2

            vpbuf2D(k)%x1 = xnew(1,l)
            vpbuf2D(k)%y1 = xnew(2,l)
            vpbuf2D(k)%q = q(l)
            vpbuf2D(k)%vz = vzint * (1-delta)

            vpbuf2D(k)%i = ixold(1,l) + dxi(1,l)
            vpbuf2D(k)%j = ixold(2,l) + dxi(2,l)

          else

            vzint = vz * delta

            ! y split 1st vp
            yint2 = 0.5_p_k_fld * dxi(2,l)
            delta = ( yint2 - xold(2,l) ) / ( yint - xold(2,l))
            xint2 = xold(1,l) + (xint - xold(1,l)) * delta

            k=k+1
            vpbuf2D(k)%x0 = xold(1,l)
            vpbuf2D(k)%y0 = xold(2,l)

            vpbuf2D(k)%x1 = xint2
            vpbuf2D(k)%y1 = yint2

            vpbuf2D(k)%q = q(l)
            vpbuf2D(k)%vz = vzint*delta

            vpbuf2D(k)%i = ixold(1,l)
            vpbuf2D(k)%j = ixold(2,l)

            k=k+1
            vpbuf2D(k)%x0 = xint2
            vpbuf2D(k)%y0 = -yint2

            vpbuf2D(k)%x1 = xint
            vpbuf2D(k)%y1 = yint - dxi(2,l)

            vpbuf2D(k)%q = q(l)
            vpbuf2D(k)%vz = vzint*(1._p_double-delta)

            vpbuf2D(k)%i = ixold(1,l)
            vpbuf2D(k)%j = ixold(2,l) + dxi(2,l)

            ! no y cross on second vp
            k=k+1
            vpbuf2D(k)%x0 = - xint
            vpbuf2D(k)%y0 = yint - dxi(2,l)
            vpbuf2D(k)%x1 = xnew(1,l)
            vpbuf2D(k)%y1 = xnew(2,l)
            vpbuf2D(k)%q = q(l)
            vpbuf2D(k)%vz = vz - vzint

            vpbuf2D(k)%i = ixold(1,l) + dxi(1,l)
            vpbuf2D(k)%j = ixold(2,l) + dxi(2,l)

          endif

    end select

  enddo

  nsplit = k

end subroutine split_2d_curved_gr
!-----------------------------------------------------------------------------

end module m_species_current_gr
