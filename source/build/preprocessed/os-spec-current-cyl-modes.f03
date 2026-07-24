# 1 "cyl_modes/os-spec-current-cyl-modes.f03"
# 1 "<built-in>" 1
# 1 "<built-in>" 3
# 467 "<built-in>" 3
# 1 "<command line>" 1
# 1 "<built-in>" 2
# 1 "cyl_modes/os-spec-current-cyl-modes.f03" 2
!-----------------------------------------------------------------------------------------
! quasi-3D specific current deposition schemes
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
# 5 "cyl_modes/os-spec-current-cyl-modes.f03" 2
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
# 6 "cyl_modes/os-spec-current-cyl-modes.f03" 2

module m_species_current_cyl_modes

use m_parameters
use m_species_cyl_modes_define
use m_species_define, only : p_cell_low, p_cell_near
use m_current_define
use m_current_cyl_modes
use m_emf_define
use m_emf_interpolate
use m_emf_cyl_modes
use m_vdf_define
use m_time_step


private

type :: t_vp_cyl_m
  ! sequence
  real(p_k_part) :: x0, y0
  real(p_k_part) :: x1, y1
  real(p_k_part) :: q
  real(p_k_part) :: vz
  integer :: i, j
end type t_vp_cyl_m


interface get_coef
  module procedure get_coef
end interface

interface split_cyl_m
  module procedure split_cyl_m
end interface

interface getjr_cyl_m_s1
   module procedure getjr_cyl_m_s1
end interface

interface getjr_cyl_m_s2
   module procedure getjr_cyl_m_s2
end interface

interface getjr_cyl_m_s3
   module procedure getjr_cyl_m_s3
end interface

public :: get_coef, getjr_cyl_m_s1, getjr_cyl_m_s2, getjr_cyl_m_s3

contains

! ----------------------------------------------------------------------------------------
subroutine get_coef( re, im, yz, np, mode )

  implicit none

  real( p_k_part ), dimension(:), intent(inout) :: re, im
  real( p_k_part ), dimension(:,:), intent(in) :: yz
  integer, intent(in) :: np, mode

  real( p_k_part ) :: r, r2, r3, re_tmp, im_tmp
  integer :: i, m

  select case (mode)

  case (0) ! for debugging purposes
    do i = 1, np
      re(i) = 1.0_p_k_part
      im(i) = 0.0_p_k_part
    enddo

  case (1)
    do i = 1, np
      r = sqrt( yz( 1, i )**2 + yz( 2, i )**2 )
      ! the if statement shouldn't be necessary once I fix the particle injection scheme
      ! if ( r /= 0.0 ) then
      re(i) = yz( 1, i ) / r
      im(i) = yz( 2, i ) / r
      ! endif
    enddo

  case (2)
    do i = 1, np
      r2 = yz( 1, i )**2 + yz( 2, i )**2
      re(i) = ( yz( 1, i )**2 - yz( 2, i )**2 ) / r2
      im(i) = 2 * yz( 1, i ) * yz( 2, i ) / r2
    enddo

  case (3)
    do i = 1, np
      r2 = yz( 1, i )**2 + yz( 2, i )**2
      r = sqrt( r2 )
      r3 = r*r2
      re(i) = +4*( yz( 1, i )**3 / r3 ) - 3*yz( 1, i ) / r
      im(i) = -4*( yz( 2, i )**3 / r3 ) + 3*yz( 2, i ) / r
    enddo

  case (4)
    do i = 1, np
      r2 = yz( 1, i )**2 + yz( 2, i )**2
      re(i) = ( (yz( 1, i )**2 - yz( 2, i )**2) / r2 )**2 - 4 * yz( 1, i )**2 * yz( 2, i )**2 / r2**2
      im(i) = 4 * ( yz( 1, i )*yz( 2, i )/r2 ) * ( (yz( 1, i )**2 - yz( 2, i )**2) / r2 )
    enddo

  ! use an algorithm to calculate the coefficients for an arbitrary mode number
  ! by the way this is very unoptimized. I just wanted to see if it would work.
  case default

    do i = 1, np

      re_tmp = 1.0_p_k_part
      im_tmp = 0.0_p_k_part
      r = sqrt( yz( 1, i )**2 + yz( 2, i )**2 )

      ! here I use the logic that e^{i m \phi} = e^{i \phi} x e^{i (m - 1) \phi}
      ! and therefore may be solved recursively
      do m = 1, mode

        re(i) = re_tmp*(yz( 1, i ) / r) - im_tmp*(yz( 2, i ) / r)
        im(i) = re_tmp*(yz( 2, i ) / r) + im_tmp*(yz( 1, i ) / r)

        re_tmp = re(i)
        im_tmp = im(i)

      enddo ! m
    enddo ! i

  end select


end subroutine get_coef
! ----------------------------------------------------------------------------------------

!-------------------------------------------------------------------------------------------------
subroutine split_cyl_m( dxi, xnew, ixold, xold, q, rgamma, u, np, vpbuf2D, nsplit, gshift_i2 , pos_type, dr, cart_xold, cart_xnew)
!-------------------------------------------------------------------------------------------------
! Splits particle trajectories so that all virtual particles have a motion starting and ending
! in the same cell. This routines also calculates vz for each virtual particle.
!
! The result is stored in the module variable vpbuf2D
!-------------------------------------------------------------------------------------------------
  implicit none

  integer, dimension(:,:), intent(in) :: dxi, ixold
  real(p_k_part), dimension(:,:), intent(in) :: xnew, xold ,u
  real(p_k_part), dimension( : ), intent(in) :: q,rgamma

  integer, intent(in) :: np ! number of particles to process

  type(t_vp_cyl_m), dimension(:), intent(out) :: vpbuf2D
  integer, intent(out) :: nsplit ! number of virtual particles created
  integer, intent(in) :: gshift_i2
  integer, intent(in) :: pos_type
  real(p_double), intent(in) :: dr
  real(p_k_part), dimension(:,:), intent(out) :: cart_xold, cart_xnew

  real(p_k_part) :: xint, yint, xint2, yint2, delta, cart_xint, cart_yint, cart_xint2, cart_yint2
  real(p_k_part) :: vz, vzint
  integer :: k,l

  integer :: cross

  k=0

  ! create virtual particles that correspond to a motion that starts and ends
  ! in the same grid cell

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

      vpbuf2D(k)%vz = u(3,l)*rgamma(l)

      vpbuf2D(k)%i = ixold(1,l)
      vpbuf2D(k)%j = ixold(2,l)

      cart_xold(1:2,k) = xold(3:4,l)
      cart_xnew(1:2,k) = xnew(3:4,l)

    case(1) ! x cross only

      xint = 0.5_p_k_part * dxi(1,l) ! + 0.5 if crosses on the right, -0.5 if crosses on the left
      delta = ( xint - xold(1,l) ) / ( xnew(1,l) - xold(1,l) )
      yint = xold(2,l) + (xnew(2,l) - xold(2,l)) * delta
      cart_xint = xold(3,l) + (xnew(3,l) - xold(3,l)) * delta
      cart_yint = xold(4,l) + (xnew(4,l) - xold(4,l)) * delta

      vz = u(3,l)*rgamma(l)

      k=k+1

      vpbuf2D(k)%x0 = xold(1,l) ! initial pos first particle is simply initial pos
      vpbuf2D(k)%y0 = xold(2,l)

      vpbuf2D(k)%x1 = xint ! left wall point (-0.5) if crosses left, right wall point (0.5) if crosses right
      vpbuf2D(k)%y1 = yint ! calc from slope with respect to r

      vpbuf2D(k)%q = q(l)
      vpbuf2D(k)%vz = vz * delta

      vpbuf2D(k)%i = ixold(1,l) ! initial position of first part is simply initial pos
      vpbuf2D(k)%j = ixold(2,l)

      cart_xold(1:2,k) = xold(3:4,l)
      cart_xnew(1,k) = cart_xint
      cart_xnew(2,k) = cart_yint

      k=k+1
      vpbuf2D(k)%x0 = -xint ! whatever the cell position the first particle ended, for the second particle it is opposite
      vpbuf2D(k)%y0 = yint

      vpbuf2D(k)%x1 = xnew(1,l) - dxi(1,l) ! subtract out the amount attributed by cell index
      vpbuf2D(k)%y1 = xnew(2,l)

      vpbuf2D(k)%q = q(l)
      vpbuf2D(k)%vz = vz * (1-delta)

      vpbuf2D(k)%i = ixold(1,l) + dxi(1,l) ! add back in the subtracted amount into cell index
      vpbuf2D(k)%j = ixold(2,l)

      cart_xold(1,k) = cart_xint
      cart_xold(2,k) = cart_yint
      cart_xnew(1:2,k) = xnew(3:4,l)


    case(2) ! y cross only

      yint = 0.5_p_k_part * dxi(2,l)
      delta = ( yint - xold(2,l) ) / ( xnew(2,l) - xold(2,l))
      xint = xold(1,l) + (xnew(1,l) - xold(1,l)) * delta
      cart_xint = xold(3,l) + (xnew(3,l) - xold(3,l)) * delta
      cart_yint = xold(4,l) + (xnew(4,l) - xold(4,l)) * delta

      vz = u(3,l)*rgamma(l)

      k=k+1
      vpbuf2D(k)%x0 = xold(1,l) ! beginning point of first particle - initial position
      vpbuf2D(k)%y0 = xold(2,l)

      vpbuf2D(k)%x1 = xint
      vpbuf2D(k)%y1 = yint
      vpbuf2D(k)%q = q(l)
      vpbuf2D(k)%vz = vz * delta

      vpbuf2D(k)%i = ixold(1,l)
      vpbuf2D(k)%j = ixold(2,l)

      cart_xold(1:2,k) = xold(3:4,l)
      cart_xnew(1,k) = cart_xint
      cart_xnew(2,k) = cart_yint

      k=k+1
      vpbuf2D(k)%x0 = xint
      vpbuf2D(k)%y0 = -yint

      vpbuf2D(k)%x1 = xnew(1,l)
      vpbuf2D(k)%y1 = xnew(2,l) - dxi(2,l)
      vpbuf2D(k)%q = q(l)
      vpbuf2D(k)%vz = vz * (1-delta)

      vpbuf2D(k)%i = ixold(1,l)
      vpbuf2D(k)%j = ixold(2,l) + dxi(2,l)

      cart_xold(1,k) = cart_xint
      cart_xold(2,k) = cart_yint
      cart_xnew(1:2,k) = xnew(3:4,l)

    case(3) ! x,y cross

      ! split in x direction first
      xint = 0.5_p_k_part * dxi(1,l)
      delta = ( xint - xold(1,l) ) / ( xnew(1,l) - xold(1,l))
      yint = xold(2,l) + ( xnew(2,l) - xold(2,l)) * delta
      cart_xint = xold(3,l) + (xnew(3,l) - xold(3,l)) * delta
      cart_yint = xold(4,l) + (xnew(4,l) - xold(4,l)) * delta

      vz = u(3,l)*rgamma(l)

      ! check if y intersection occured for 1st or 2nd split
      if ((yint >= -0.5_p_k_part) .and. ( yint < 0.5_p_k_part )) then

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

        cart_xold(1:2,k) = xold(3:4,l)
        cart_xnew(1,k) = cart_xint
        cart_xnew(2,k) = cart_yint

        ! y split 2nd vp
        k=k+1

        yint2 = 0.5_p_k_part * dxi(2,l)
        delta = ( yint2 - yint ) / ( xnew(2,l) - yint )
        xint2 = -xint + ( xnew(1,l) - xint ) * delta
        cart_xint2 = cart_xint + (xnew(3,l) - cart_xint) * delta
        cart_yint2 = cart_yint + (xnew(4,l) - cart_yint) * delta

        vpbuf2D(k)%x0 = -xint
        vpbuf2D(k)%y0 = yint

        vpbuf2D(k)%x1 = xint2
        vpbuf2D(k)%y1 = yint2

        vpbuf2D(k)%q = q(l)
        vpbuf2D(k)%vz = vzint * delta

        vpbuf2D(k)%i = ixold(1,l) + dxi(1,l)
        vpbuf2D(k)%j = ixold(2,l)

        cart_xold(1,k) = cart_xint
        cart_xold(2,k) = cart_yint
        cart_xnew(1,k) = cart_xint2
        cart_xnew(2,k) = cart_yint2

        k=k+1

        vpbuf2D(k)%x0 = xint2
        vpbuf2D(k)%y0 = -yint2

        vpbuf2D(k)%x1 = xnew(1,l) - dxi(1,l)
        vpbuf2D(k)%y1 = xnew(2,l) - dxi(2,l)
        vpbuf2D(k)%q = q(l)
        vpbuf2D(k)%vz = vzint * (1-delta)

        vpbuf2D(k)%i = ixold(1,l) + dxi(1,l)
        vpbuf2D(k)%j = ixold(2,l) + dxi(2,l)

        cart_xold(1,k) = cart_xint2
        cart_xold(2,k) = cart_yint2
        cart_xnew(1:2,k) = xnew(3:4,l)


      else

        vzint = vz * delta

        ! y split 1st vp
        yint2 = 0.5_p_k_part * dxi(2,l)
        delta = ( yint2 - xold(2,l) ) / ( yint - xold(2,l))
        xint2 = xold(1,l) + (xint - xold(1,l)) * delta
        cart_xint2 = xold(3,l) + (cart_xint - xold(3,l)) * delta
        cart_yint2 = xold(4,l) + (cart_yint - xold(4,l)) * delta

        k=k+1
        vpbuf2D(k)%x0 = xold(1,l)
        vpbuf2D(k)%y0 = xold(2,l)

        vpbuf2D(k)%x1 = xint2
        vpbuf2D(k)%y1 = yint2

        vpbuf2D(k)%q = q(l)
        vpbuf2D(k)%vz = vzint*delta

        vpbuf2D(k)%i = ixold(1,l)
        vpbuf2D(k)%j = ixold(2,l)

        cart_xold(1:2,k) = xold(3:4,l)
        cart_xnew(1,k) = cart_xint2
        cart_xnew(2,k) = cart_yint2


        k=k+1
        vpbuf2D(k)%x0 = xint2
        vpbuf2D(k)%y0 = -yint2

        vpbuf2D(k)%x1 = xint
        vpbuf2D(k)%y1 = yint - dxi(2,l)

        vpbuf2D(k)%q = q(l)
        vpbuf2D(k)%vz = vzint*(1-delta)

        vpbuf2D(k)%i = ixold(1,l)
        vpbuf2D(k)%j = ixold(2,l) + dxi(2,l)

        cart_xold(1,k) = cart_xint2
        cart_xold(2,k) = cart_yint2
        cart_xnew(1,k) = cart_xint
        cart_xnew(2,k) = cart_yint

        ! no y cross on second vp
        k=k+1
        vpbuf2D(k)%x0 = - xint
        vpbuf2D(k)%y0 = yint - dxi(2,l)
        vpbuf2D(k)%x1 = xnew(1,l) - dxi(1,l)
        vpbuf2D(k)%y1 = xnew(2,l) - dxi(2,l)
        vpbuf2D(k)%q = q(l)
        vpbuf2D(k)%vz = vz - vzint

        vpbuf2D(k)%i = ixold(1,l) + dxi(1,l)
        vpbuf2D(k)%j = ixold(2,l) + dxi(2,l)

        cart_xold(1,k) = cart_xint
        cart_xold(2,k) = cart_yint
        cart_xnew(1:2,k) = xnew(3:4,l)

      endif

    end select

  enddo

  nsplit = k

end subroutine split_cyl_m
!-------------------------------------------------------------------------------------------------

!-------------------------------------------------------------------------------------------------
subroutine getjr_cyl_m_s1( jay_re, jay_im, dxi, xnew, ixold, xold, q, rgamma, u, np, dt , gshift_i2, mode)
!---------------------------------------------------
! Linear interpolation with cell indexed positions
!---------------------------------------------------
  implicit none

  ! dummy variables
  integer, parameter :: rank = 2
  type( t_vdf ), intent(inout) :: jay_re
  type( t_vdf ), intent(inout) :: jay_im

  integer, dimension(:,:), intent(in) :: dxi, ixold
  real(p_k_part), dimension(:,:), intent(in) :: xnew, xold , u
  real(p_k_part), dimension( : ), intent(in) :: q,rgamma

  integer, intent(in) :: np
  real(p_double), intent(in) :: dt
  integer, intent(in) :: gshift_i2
  integer, intent(in) :: mode

  ! local variables
  integer :: l, nsplit, ix, jx
  ! integer :: source_p

  real(p_k_fld), dimension(0:1) :: S0x, S1x, S0y, S1y
  real(p_k_fld), dimension(0:1) :: wp1, wp2
  real(p_k_fld) :: wl1, wl2

  real(p_k_fld) :: x0,x1,y0,y1, rmid
  real(p_k_fld) :: qnx, qny, qvz
  real(p_k_fld) :: jnorm1, jnorm2, jnorm3

  real(p_k_part), dimension(rank,3*p_cache_size) :: x_center, x_new, x_old
  real(p_k_part), dimension(3*p_cache_size) :: coeff_re, coeff_im, coeff_re_p, coeff_im_p, coeff_re_m, coeff_im_m
  real(p_k_part) :: factor_1, factor_2

  type(t_vp_cyl_m), dimension( 3*p_cache_size ) :: vpbuf2D

  ! executable statements
  ! split particles
  call split_cyl_m( dxi, xnew, ixold, xold, q, rgamma, u, np, vpbuf2D, nsplit, gshift_i2, p_cell_low, jay_re%dx_(2), x_old, x_new)

  jnorm1 = real( jay_re%dx_(1) / dt / 2, p_k_fld )
  jnorm2 = real( jay_re%dx_(2) / dt / 2, p_k_fld )
  jnorm3 = real( jay_re%dx_(2) / (mode*dt), p_k_fld ) ! dr / m dt

  x_center(:,1:nsplit) = (x_new(:,1:nsplit) + x_old(:,1:nsplit))/2.0_p_k_fld

  call get_coef(coeff_re, coeff_im, x_center(:,:), nsplit, mode) ! e^(i m th_av )
  call get_coef(coeff_re_p, coeff_im_p, x_new(:,:), nsplit, mode) ! e^(i m th_p )
  call get_coef(coeff_re_m, coeff_im_m, x_old(:,:), nsplit, mode) ! e^(i m th_m )

  ! now accumulate jay looping through all virtual particles

  do l=1,nsplit ! is nplit the total number of particles in the split particle buffer? (guess so)
    ! order 1 charge conserving current deposition
    ! generated automatically by z-2.1

    x0 = vpbuf2D(l)%x0
    x1 = vpbuf2D(l)%x1
    y0 = vpbuf2D(l)%y0
    y1 = vpbuf2D(l)%y1
    ix = vpbuf2D(l)%i
    jx = vpbuf2D(l)%j
    ! integer pointing to the source particle for theta dependent coefficients
    ! source_p = vpbuf2D(l)%sp

    rmid = jx + gshift_i2 - 0.5_p_k_fld ! note this is the r of center cell
    ! rmid = jx - 0.5_p_k_fld ! note this is the r of center cell -- the gshift_i2 is redundant since jx is already shifted

    ! Normalize charge
    qnx = real( vpbuf2D(l)%q, p_k_fld ) * jnorm1
    qny = real( vpbuf2D(l)%q, p_k_fld ) * jnorm2
    qvz = real( vpbuf2D(l)%q, p_k_fld ) * jnorm3
    ! qvz = real( vpbuf2D(l)%q * vpbuf2D(l)%vz, p_k_fld )/3

    ! I am wondering whether the current deposited at the axis should be doubled to take into account a
    ! "mirror" current
    ! if (jx == 1) then
    ! qnx = 2.0*qnx
    ! qny = 2.0*qny
    ! qvz = 2.0*qvz
    ! endif

    ! if (jx == 2) then
    ! mult_factor = 2.0
    ! else
    ! mult_factor = 1.0
    ! endif

    ! get spline weitghts for x and y
    S0x(0) = 0.5 - x0
    S0x(1) = 0.5 + x0

    S1x(0) = 0.5 - x1
    S1x(1) = 0.5 + x1

    S0y(0) = 0.5 - y0
    S0y(1) = 0.5 + y0

    S1y(0) = 0.5 - y1
    S1y(1) = 0.5 + y1


    ! get longitudinal motion weights
    wl1 = qnx*(-x0 + x1)

    wl2 = qny*(-y0 + y1)

    ! get perpendicular motion weights
    wp1(0) = S0y(0) + S1y(0)
    wp1(1) = S0y(1) + S1y(1)

    wp2(0) = S0x(0) + S1x(0)
    wp2(1) = S0x(1) + S1x(1)

    ! accumulate j1_re
    jay_re%f2(1,ix ,jx ) = jay_re%f2(1,ix ,jx ) + wl1 * wp1(0) * coeff_re(l)
    jay_re%f2(1,ix ,jx+1) = jay_re%f2(1,ix ,jx+1) + wl1 * wp1(1) * coeff_re(l)

    ! accumulate j2_re
    jay_re%f2(2,ix ,jx ) = jay_re%f2(2,ix ,jx ) + wl2 * wp2(0) * coeff_re(l)
    jay_re%f2(2,ix+1,jx ) = jay_re%f2(2,ix+1,jx ) + wl2 * wp2(1) * coeff_re(l)

    ! accumulate j1_im
    jay_im%f2(1,ix ,jx ) = jay_im%f2(1,ix ,jx ) + wl1 * wp1(0) * coeff_im(l)
    jay_im%f2(1,ix ,jx+1) = jay_im%f2(1,ix ,jx+1) + wl1 * wp1(1) * coeff_im(l)

    ! accumulate j2_im
    jay_im%f2(2,ix ,jx ) = jay_im%f2(2,ix ,jx ) + wl2 * wp2(0) * coeff_im(l)
    jay_im%f2(2,ix+1,jx ) = jay_im%f2(2,ix+1,jx ) + wl2 * wp2(1) * coeff_im(l)


    ! accumulate j3
    ! charge conserving cylindrical mode current deposition algorithm
    ! real part
    factor_1 = coeff_im_p(l) - coeff_im(l) !-(coeff_im_p(l) - coeff_im(l))
    factor_2 = coeff_im_m(l) - coeff_im(l) !-(coeff_im_m(l) - coeff_im(l))

    jay_re%f2(3,ix , jx ) = jay_re%f2(3,ix ,jx ) + qvz * ABS((rmid )) * (S1x(0)*S1y(0)*factor_1 - S0x(0)*S0y(0)*factor_2)
    jay_re%f2(3,ix+1, jx ) = jay_re%f2(3,ix+1,jx ) + qvz * ABS((rmid )) * (S1x(1)*S1y(0)*factor_1 - S0x(1)*S0y(0)*factor_2)
    jay_re%f2(3,ix , jx+1) = jay_re%f2(3,ix ,jx+1) + qvz * ABS((rmid + 1 )) * (S1x(0)*S1y(1)*factor_1 - S0x(0)*S0y(1)*factor_2)
    jay_re%f2(3,ix+1, jx+1) = jay_re%f2(3,ix+1,jx+1) + qvz * ABS((rmid + 1 )) * (S1x(1)*S1y(1)*factor_1 - S0x(1)*S0y(1)*factor_2)


    ! imaginary part
    factor_1 = -(coeff_re_p(l) - coeff_re(l)) !coeff_re_p(l) - coeff_re(l)
    factor_2 = -(coeff_re_m(l) - coeff_re(l)) !coeff_re_m(l) - coeff_re(l)

    jay_im%f2(3,ix , jx ) = jay_im%f2(3,ix ,jx ) + qvz * ABS((rmid )) * (S1x(0)*S1y(0)*factor_1 - S0x(0)*S0y(0)*factor_2)
    jay_im%f2(3,ix+1, jx ) = jay_im%f2(3,ix+1,jx ) + qvz * ABS((rmid )) * (S1x(1)*S1y(0)*factor_1 - S0x(1)*S0y(0)*factor_2)
    jay_im%f2(3,ix , jx+1) = jay_im%f2(3,ix ,jx+1) + qvz * ABS((rmid + 1 )) * (S1x(0)*S1y(1)*factor_1 - S0x(0)*S0y(1)*factor_2)
    jay_im%f2(3,ix+1, jx+1) = jay_im%f2(3,ix+1,jx+1) + qvz * ABS((rmid + 1 )) * (S1x(1)*S1y(1)*factor_1 - S0x(1)*S0y(1)*factor_2)

  enddo


end subroutine getjr_cyl_m_s1
!-------------------------------------------------------------------------------------------------


!-------------------------------------------------------------------------------------------------
subroutine getjr_cyl_m_s2( jay_re, jay_im, dxi, xnew, ixold, xold, q, rgamma, u, np, dt, gshift_i2, mode)
!-------------------------------------------------------------------------------------------------
! Quadratic (n=2) interpolation with cell indexed positions
!-------------------------------------------------------------------------------------------------

  implicit none

  ! dummy variables
  integer, parameter :: rank = 2
  type( t_vdf ), intent(inout) :: jay_re
  type( t_vdf ), intent(inout) :: jay_im

  integer, dimension(:,:), intent(in) :: dxi, ixold
  real(p_k_part), dimension(:,:), intent(in) :: xnew, xold ,u
  real(p_k_part), dimension( : ), intent(in) :: q,rgamma

  integer, intent(in) :: np
  real(p_double), intent(in) :: dt
  integer, intent(in) :: gshift_i2
  integer, intent(in) :: mode

  ! local variables
  integer :: l, nsplit, ix, jx
  ! integer :: source_p

  real(p_k_fld), dimension(-1:1) :: S0x, S1x, S0y, S1y
  real(p_k_fld), dimension(-1:1) :: wp1, wp2
  real(p_k_fld), dimension(-1:0) :: wl1, wl2

  real(p_k_fld) :: x0,x1,y0,y1, rmid
  real(p_k_fld) :: qnx, qny, qvz
  real(p_k_fld) :: jnorm1, jnorm2, jnorm3

  real(p_k_part), dimension(2,3*p_cache_size) :: x_center, x_new, x_old
  real(p_k_part), dimension(3*p_cache_size) :: coeff_re, coeff_im, coeff_re_p, coeff_im_p, coeff_re_m, coeff_im_m
  real(p_k_part) :: factor_1, factor_2

  type(t_vp_cyl_m), dimension( 3*p_cache_size ) :: vpbuf2D ! must change this to cyl_m

  ! executable statements

  ! split particles
  call split_cyl_m( dxi, xnew, ixold, xold, q, rgamma, u, np, vpbuf2D, nsplit, gshift_i2 , p_cell_near, jay_re%dx_(2), x_old, x_new)

  ! now accumulate jay looping through all virtual particles
  ! and shifting grid indexes
  jnorm1 = real( jay_re%dx_(1)/dt/4, p_k_fld )
  jnorm2 = real( jay_re%dx_(2)/dt/4, p_k_fld )
  jnorm3 = real( jay_re%dx_(2) / (mode*dt), p_k_fld )

  x_center(:,1:nsplit) = (x_new(:,1:nsplit) + x_old(:,1:nsplit))/2.0_p_k_fld

  call get_coef(coeff_re, coeff_im, x_center(:,:), nsplit, mode) ! e^(i m th_av )
  call get_coef(coeff_re_p, coeff_im_p, x_new(:,:), nsplit, mode) ! e^(i m th_p )
  call get_coef(coeff_re_m, coeff_im_m, x_old(:,:), nsplit, mode) ! e^(i m th_m )


  ! get e^{i m (delta th / 2.0} coefficients
  ! call get_coef(coeff_re_p, coeff_im_p, xnew(3:4,:), np, mode) ! e^(i m th_p )
  ! call get_coef(coeff_re_m, coeff_im_m, xold(3:4,:), np, mode) ! e^(i m th_m )
  ! coeff_im_m(:) = -coeff_im_m(:) ! flip sign of imaginary part to make into e^(-i m th_m)
  ! coeff_re_sq(:) = coeff_re_p(:)*coeff_re_m(:) - coeff_im_p(:)*coeff_im_m(:)
  ! coeff_im_sq(:) = coeff_im_p(:)*coeff_re_m(:) + coeff_re_p*coeff_im_m(:) ! e^{im(th_p - th_m)} = e^{i m th_p} e^{- i m th_m }
  ! coeff_re(:) = sqrt((coeff_re_sq(:) + sqrt(coeff_re_sq(:)**2 + coeff_im_sq(:)**2))/2.0) ! real part of complex sqrt()
  ! coeff_im(:) = sign(1.0, coeff_im_sq(:)) * sqrt( (-coeff_re_sq(:) + sqrt(coeff_re_sq(:)**2 + coeff_im_sq(:)**2))/2.0) !im part of complex sqrt()


  do l=1,nsplit

    ! order 2 charge conserving current deposition
    ! requires near point splitting with z motion split
    ! generated automatically by z-2.0

    ! manual optimizations :
    ! i) optimized dx1sdt, dx2sdt (moved a division by 4 out of the loop)
    ! ii) optimized longitudinal weights
    ! iii) optimized spline functions

    x0 = vpbuf2D(l)%x0
    y0 = vpbuf2D(l)%y0
    x1 = vpbuf2D(l)%x1
    y1 = vpbuf2D(l)%y1
    ix = vpbuf2D(l)%i
    jx = vpbuf2D(l)%j


    ! pointer to the source particle, for extracting theta-dependent coefficients
    ! source_p = vpbuf2D(l)%sp
    ! source_p = l

    rmid = jx + gshift_i2 - 0.5_p_k_fld ! note this is the r of center cell


    ! Normalize charge
    qnx = real( vpbuf2D(l)%q, p_k_fld ) * jnorm1
    qny = real( vpbuf2D(l)%q, p_k_fld ) * jnorm2
    qvz = real( vpbuf2D(l)%q, p_k_fld ) * jnorm3
    ! qvz = real( vpbuf2D(l)%q * vpbuf2D(l)%vz, p_k_fld ) / 3

    ! get spline weitghts for x and y
    S0x(-1) = 0.5_p_k_fld*(0.5_p_k_fld - x0)**2
    S0x(0) = 0.75 - x0**2
    S0x(1) = 0.5_p_k_fld*(0.5_p_k_fld + x0)**2

    S1x(-1) = 0.5_p_k_fld*(0.5_p_k_fld - x1)**2
    S1x(0) = 0.75 - x1**2
    S1x(1) = 0.5_p_k_fld*(0.5_p_k_fld + x1)**2

    S0y(-1) = 0.5_p_k_fld*(0.5_p_k_fld - y0)**2
    S0y(0) = 0.75 - y0**2
    S0y(1) = 0.5_p_k_fld*(0.5_p_k_fld + y0)**2

    S1y(-1) = 0.5_p_k_fld*(0.5_p_k_fld - y1)**2
    S1y(0) = 0.75 - y1**2
    S1y(1) = 0.5_p_k_fld*(0.5_p_k_fld + y1)**2

    ! get longitudinal motion weights
    wl1(-1) = (qnx*(x0 - x1)*(-1 + x0 + x1))
    wl1(0) = (qnx*(-(x0*(1 + x0)) + x1*(1 + x1)))

    wl2(-1) = (qny*(y0 - y1)*(-1 + y0 + y1))
    wl2(0) = (qny*(-(y0*(1 + y0)) + y1*(1 + y1)))

    ! get perpendicular motion weights
    wp1(-1) = (S0y(-1)+S1y(-1))
    wp1(0) = (S0y(0)+S1y(0))
    wp1(1) = (S0y(1)+S1y(1))

    wp2(-1) = (S0x(-1)+S1x(-1))
    wp2(0) = (S0x(0)+S1x(0))
    wp2(1) = (S0x(1)+S1x(1))

    ! accumulate j1_re
    jay_re%f2(1,ix-1,jx-1) = jay_re%f2(1,ix-1,jx-1) + wl1(-1) * wp1(-1) * coeff_re(l)
    jay_re%f2(1,ix ,jx-1) = jay_re%f2(1,ix ,jx-1) + wl1(0) * wp1(-1) * coeff_re(l)
    jay_re%f2(1,ix-1,jx ) = jay_re%f2(1,ix-1,jx ) + wl1(-1) * wp1(0) * coeff_re(l)
    jay_re%f2(1,ix ,jx ) = jay_re%f2(1,ix ,jx ) + wl1(0) * wp1(0) * coeff_re(l)
    jay_re%f2(1,ix-1,jx+1) = jay_re%f2(1,ix-1,jx+1) + wl1(-1) * wp1(1) * coeff_re(l)
    jay_re%f2(1,ix ,jx+1) = jay_re%f2(1,ix ,jx+1) + wl1(0) * wp1(1) * coeff_re(l)

    ! accumulate j2_re
    jay_re%f2(2,ix-1,jx-1) = jay_re%f2(2,ix-1,jx-1) + wl2(-1) * wp2(-1) * coeff_re(l)
    jay_re%f2(2,ix ,jx-1) = jay_re%f2(2,ix ,jx-1) + wl2(-1) * wp2(0) * coeff_re(l)
    jay_re%f2(2,ix+1,jx-1) = jay_re%f2(2,ix+1,jx-1) + wl2(-1) * wp2(1) * coeff_re(l)
    jay_re%f2(2,ix-1,jx ) = jay_re%f2(2,ix-1,jx ) + wl2(0) * wp2(-1) * coeff_re(l)
    jay_re%f2(2,ix ,jx ) = jay_re%f2(2,ix ,jx ) + wl2(0) * wp2(0) * coeff_re(l)
    jay_re%f2(2,ix+1,jx ) = jay_re%f2(2,ix+1,jx ) + wl2(0) * wp2(1) * coeff_re(l)

    ! accumulate j1_im
    jay_im%f2(1,ix-1,jx-1) = jay_im%f2(1,ix-1,jx-1) + wl1(-1) * wp1(-1) * coeff_im(l)
    jay_im%f2(1,ix ,jx-1) = jay_im%f2(1,ix ,jx-1) + wl1(0) * wp1(-1) * coeff_im(l)
    jay_im%f2(1,ix-1,jx ) = jay_im%f2(1,ix-1,jx ) + wl1(-1) * wp1(0) * coeff_im(l)
    jay_im%f2(1,ix ,jx ) = jay_im%f2(1,ix ,jx ) + wl1(0) * wp1(0) * coeff_im(l)
    jay_im%f2(1,ix-1,jx+1) = jay_im%f2(1,ix-1,jx+1) + wl1(-1) * wp1(1) * coeff_im(l)
    jay_im%f2(1,ix ,jx+1) = jay_im%f2(1,ix ,jx+1) + wl1(0) * wp1(1) * coeff_im(l)

    ! accumulate j2_im
    jay_im%f2(2,ix-1,jx-1) = jay_im%f2(2,ix-1,jx-1) + wl2(-1) * wp2(-1) * coeff_im(l)
    jay_im%f2(2,ix ,jx-1) = jay_im%f2(2,ix ,jx-1) + wl2(-1) * wp2(0) * coeff_im(l)
    jay_im%f2(2,ix+1,jx-1) = jay_im%f2(2,ix+1,jx-1) + wl2(-1) * wp2(1) * coeff_im(l)
    jay_im%f2(2,ix-1,jx ) = jay_im%f2(2,ix-1,jx ) + wl2(0) * wp2(-1) * coeff_im(l)
    jay_im%f2(2,ix ,jx ) = jay_im%f2(2,ix ,jx ) + wl2(0) * wp2(0) * coeff_im(l)
    jay_im%f2(2,ix+1,jx ) = jay_im%f2(2,ix+1,jx ) + wl2(0) * wp2(1) * coeff_im(l)

    ! accumulate j3
    ! charge conserving cylindrical mode current deposition algorithm
    ! real part
    factor_1 = coeff_im_p(l) - coeff_im(l) !-(coeff_im_p(l) - coeff_im(l))
    factor_2 = coeff_im_m(l) - coeff_im(l) !-(coeff_im_m(l) - coeff_im(l))

    jay_re%f2(3,ix-1, jx-1) = jay_re%f2(3,ix-1,jx-1) + qvz * ABS((rmid-1)) * (S1x(-1)*S1y(-1)*factor_1 - S0x(-1)*S0y(-1)*factor_2)
    jay_re%f2(3,ix , jx-1) = jay_re%f2(3,ix ,jx-1) + qvz * ABS((rmid-1)) * (S1x(0)*S1y(-1)*factor_1 - S0x(0)*S0y(-1)*factor_2)
    jay_re%f2(3,ix+1, jx-1) = jay_re%f2(3,ix+1,jx-1) + qvz * ABS((rmid-1)) * (S1x(1)*S1y(-1)*factor_1 - S0x(1)*S0y(-1)*factor_2)
    jay_re%f2(3,ix-1, jx ) = jay_re%f2(3,ix-1,jx ) + qvz * ABS((rmid )) * (S1x(-1)*S1y(0)*factor_1 - S0x(-1)*S0y(0)*factor_2)
    jay_re%f2(3,ix , jx ) = jay_re%f2(3,ix ,jx ) + qvz * ABS((rmid )) * (S1x(0)*S1y(0)*factor_1 - S0x(0)*S0y(0)*factor_2)
    jay_re%f2(3,ix+1, jx ) = jay_re%f2(3,ix+1,jx ) + qvz * ABS((rmid )) * (S1x(1)*S1y(0)*factor_1 - S0x(1)*S0y(0)*factor_2)
    jay_re%f2(3,ix-1, jx+1) = jay_re%f2(3,ix-1,jx+1) + qvz * ABS((rmid+1)) * (S1x(-1)*S1y(1)*factor_1 - S0x(-1)*S0y(1)*factor_2)
    jay_re%f2(3,ix , jx+1) = jay_re%f2(3,ix ,jx+1) + qvz * ABS((rmid+1)) * (S1x(0)*S1y(1)*factor_1 - S0x(0)*S0y(1)*factor_2)
    jay_re%f2(3,ix+1, jx+1) = jay_re%f2(3,ix+1,jx+1) + qvz * ABS((rmid+1)) * (S1x(1)*S1y(1)*factor_1 - S0x(1)*S0y(1)*factor_2)

    ! imaginary part
    factor_1 = -(coeff_re_p(l) - coeff_re(l)) !coeff_re_p(l) - coeff_re(l)
    factor_2 = -(coeff_re_m(l) - coeff_re(l)) !coeff_re_m(l) - coeff_re(l)


    jay_im%f2(3,ix-1, jx-1) = jay_im%f2(3,ix-1,jx-1) + qvz * ABS((rmid-1)) * (S1x(-1)*S1y(-1)*factor_1 - S0x(-1)*S0y(-1)*factor_2)
    jay_im%f2(3,ix , jx-1) = jay_im%f2(3,ix ,jx-1) + qvz * ABS((rmid-1)) * (S1x(0)*S1y(-1)*factor_1 - S0x(0)*S0y(-1)*factor_2)
    jay_im%f2(3,ix+1, jx-1) = jay_im%f2(3,ix+1,jx-1) + qvz * ABS((rmid-1)) * (S1x(1)*S1y(-1)*factor_1 - S0x(1)*S0y(-1)*factor_2)
    jay_im%f2(3,ix-1, jx ) = jay_im%f2(3,ix-1,jx ) + qvz * ABS((rmid )) * (S1x(-1)*S1y(0)*factor_1 - S0x(-1)*S0y(0)*factor_2)
    jay_im%f2(3,ix , jx ) = jay_im%f2(3,ix ,jx ) + qvz * ABS((rmid )) * (S1x(0)*S1y(0)*factor_1 - S0x(0)*S0y(0)*factor_2)
    jay_im%f2(3,ix+1, jx ) = jay_im%f2(3,ix+1,jx ) + qvz * ABS((rmid )) * (S1x(1)*S1y(0)*factor_1 - S0x(1)*S0y(0)*factor_2)
    jay_im%f2(3,ix-1, jx+1) = jay_im%f2(3,ix-1,jx+1) + qvz * ABS((rmid+1)) * (S1x(-1)*S1y(1)*factor_1 - S0x(-1)*S0y(1)*factor_2)
    jay_im%f2(3,ix , jx+1) = jay_im%f2(3,ix ,jx+1) + qvz * ABS((rmid+1)) * (S1x(0)*S1y(1)*factor_1 - S0x(0)*S0y(1)*factor_2)
    jay_im%f2(3,ix+1, jx+1) = jay_im%f2(3,ix+1,jx+1) + qvz * ABS((rmid+1)) * (S1x(1)*S1y(1)*factor_1 - S0x(1)*S0y(1)*factor_2)


  enddo

end subroutine getjr_cyl_m_s2
!-------------------------------------------------------------------------------------------------

!-------------------------------------------------------------------------------------------------
subroutine getjr_cyl_m_s3( jay_re, jay_im, dxi, xnew, ixold, xold, q, rgamma, u, np, dt, gshift_i2, mode)
!-------------------------------------------------------------------------------------------------
! Cubic (n=3) interpolation with cell indexed positions
!-------------------------------------------------------------------------------------------------

  implicit none


  ! dummy variables
  integer, parameter :: rank = 2
  type( t_vdf ), intent(inout) :: jay_re
  type( t_vdf ), intent(inout) :: jay_im

  integer, dimension(:,:), intent(in) :: dxi, ixold
  real(p_k_part), dimension(:,:), intent(in) :: xnew, xold ,u
  real(p_k_part), dimension( : ), intent(in) :: q,rgamma

  integer, intent(in) :: np
  real(p_double), intent(in) :: dt
  integer, intent(in) :: gshift_i2
  integer, intent(in) :: mode

  ! local variables
  integer :: l, nsplit, ix, jx
  ! integer :: source_p

  real(p_k_fld), dimension(-1:2) :: S0x, S1x, S0y, S1y
  real(p_k_fld), dimension(-1:2) :: wp1, wp2
  real(p_k_fld), dimension(-1:1) :: wl1, wl2

  real(p_k_fld) :: x0,x1,y0,y1, rmid
  real(p_k_fld) :: qnx, qny, qvz
  real(p_k_fld) :: jnorm1, jnorm2, jnorm3

  real(p_k_part), dimension(2,3*p_cache_size) :: x_center, x_old, x_new
  real(p_k_part), dimension(3*p_cache_size) :: coeff_re, coeff_im, coeff_re_p, coeff_im_p, coeff_re_m, coeff_im_m
  real(p_k_part) :: factor_1, factor_2

  type(t_vp_cyl_m), dimension( 3*p_cache_size ) :: vpbuf2D

  ! executable statements

  ! split particles
  call split_cyl_m( dxi, xnew, ixold, xold, q, rgamma, u, np, vpbuf2D, nsplit, gshift_i2 , p_cell_low, jay_re%dx_(2), x_old, x_new)

  x_center(:,1:nsplit) = (x_new(:,1:nsplit) + x_old(:,1:nsplit))/2.0_p_k_fld

  ! now accumulate jay looping through all virtual particles
  ! and shifting grid indexes
  jnorm1 = real( jay_re%dx_(1) / dt , p_k_fld )
  jnorm2 = real( jay_re%dx_(2) / dt , p_k_fld )
  jnorm3 = real( jay_re%dx_(2) / (mode*dt), p_k_fld )

  ! x_center(:,:) = (xnew(3:4,:) + xold(3:4,:))/2.0

  call get_coef(coeff_re, coeff_im, x_center(:,:), nsplit, mode) ! e^(i m th_av )
  call get_coef(coeff_re_p, coeff_im_p, x_new(:,:), nsplit, mode) ! e^(i m th_p )
  call get_coef(coeff_re_m, coeff_im_m, x_old(:,:), nsplit, mode) ! e^(i m th_m )

  do l=1,nsplit

    ! order 3 charge conserving current deposition
    ! generated automatically by z-2.1

    x0 = vpbuf2D(l)%x0
    x1 = vpbuf2D(l)%x1
    y0 = vpbuf2D(l)%y0
    y1 = vpbuf2D(l)%y1
    ix = vpbuf2D(l)%i
    jx = vpbuf2D(l)%j
    ! pointer to source particle for extracting theta-dependent coefficients
    ! source_p = vpbuf2D(l)%sp

    rmid = jx + gshift_i2 - 0.5_p_k_fld

    ! Normalize charge
    qnx = real( vpbuf2D(l)%q, p_k_fld ) * jnorm1
    qny = real( vpbuf2D(l)%q, p_k_fld ) * jnorm2
    ! qvz = real( vpbuf2D(l)%q * vpbuf2D(l)%vz, p_k_fld )/3
    qvz = real( vpbuf2D(l)%q, p_k_fld) * jnorm3

    ! get spline weitghts for x and y
    S0x(-1) = -(-0.5 + x0)**3/6.
    S0x(0) = (4 - 6*(0.5 + x0)**2 + 3*(0.5 + x0)**3)/6.
    S0x(1) = (23 + 30*x0 - 12*x0**2 - 24*x0**3)/48.
    S0x(2) = (0.5 + x0)**3/6.

    S1x(-1) = -(-0.5 + x1)**3/6.
    S1x(0) = (4 - 6*(0.5 + x1)**2 + 3*(0.5 + x1)**3)/6.
    S1x(1) = (23 + 30*x1 - 12*x1**2 - 24*x1**3)/48.
    S1x(2) = (0.5 + x1)**3/6.

    S0y(-1) = -(-0.5 + y0)**3/6.
    S0y(0) = (4 - 6*(0.5 + y0)**2 + 3*(0.5 + y0)**3)/6.
    S0y(1) = (23 + 30*y0 - 12*y0**2 - 24*y0**3)/48.
    S0y(2) = (0.5 + y0)**3/6.

    S1y(-1) = -(-0.5 + y1)**3/6.
    S1y(0) = (4 - 6*(0.5 + y1)**2 + 3*(0.5 + y1)**3)/6.
    S1y(1) = (23 + 30*y1 - 12*y1**2 - 24*y1**3)/48.
    S1y(2) = (0.5 + y1)**3/6.


    ! get longitudinal motion weights
    wl1(-1) = (qnx*(-(-0.5 + x0)**3 + (-0.5 + x1)**3))/6.
    wl1(0) = (qnx*(-9*x0 + 4*x0**3 + 9*x1 - 4*x1**3))/12.
    wl1(1) = (qnx*(-(x0*(3 + 6*x0 + 4*x0**2)) + x1*(3 + 6*x1 + 4*x1**2)))/24.

    wl2(-1) = (qny*(-(-0.5 + y0)**3 + (-0.5 + y1)**3))/6.
    wl2(0) = (qny*(-9*y0 + 4*y0**3 + 9*y1 - 4*y1**3))/12.
    wl2(1) = (qny*(-(y0*(3 + 6*y0 + 4*y0**2)) + y1*(3 + 6*y1 + 4*y1**2)))/24.

    ! get perpendicular motion weights
    wp1(-1) = 0.5_p_k_part*(S0y(-1)+S1y(-1))
    wp1(0) = 0.5_p_k_part*(S0y(0)+S1y(0))
    wp1(1) = 0.5_p_k_part*(S0y(1)+S1y(1))
    wp1(2) = 0.5_p_k_part*(S0y(2)+S1y(2))

    wp2(-1) = 0.5_p_k_part*(S0x(-1)+S1x(-1))
    wp2(0) = 0.5_p_k_part*(S0x(0)+S1x(0))
    wp2(1) = 0.5_p_k_part*(S0x(1)+S1x(1))
    wp2(2) = 0.5_p_k_part*(S0x(2)+S1x(2))

    ! accumulate j1_re
    jay_re%f2(1,ix-1,jx-1) = jay_re%f2(1,ix-1,jx-1) + wl1(-1) * wp1(-1) * coeff_re(l)
    jay_re%f2(1,ix ,jx-1) = jay_re%f2(1,ix ,jx-1) + wl1(0) * wp1(-1) * coeff_re(l)
    jay_re%f2(1,ix+1,jx-1) = jay_re%f2(1,ix+1,jx-1) + wl1(1) * wp1(-1) * coeff_re(l)
    jay_re%f2(1,ix-1,jx ) = jay_re%f2(1,ix-1,jx ) + wl1(-1) * wp1(0) * coeff_re(l)
    jay_re%f2(1,ix ,jx ) = jay_re%f2(1,ix ,jx ) + wl1(0) * wp1(0) * coeff_re(l)
    jay_re%f2(1,ix+1,jx ) = jay_re%f2(1,ix+1,jx ) + wl1(1) * wp1(0) * coeff_re(l)
    jay_re%f2(1,ix-1,jx+1) = jay_re%f2(1,ix-1,jx+1) + wl1(-1) * wp1(1) * coeff_re(l)
    jay_re%f2(1,ix ,jx+1) = jay_re%f2(1,ix ,jx+1) + wl1(0) * wp1(1) * coeff_re(l)
    jay_re%f2(1,ix+1,jx+1) = jay_re%f2(1,ix+1,jx+1) + wl1(1) * wp1(1) * coeff_re(l)
    jay_re%f2(1,ix-1,jx+2) = jay_re%f2(1,ix-1,jx+2) + wl1(-1) * wp1(2) * coeff_re(l)
    jay_re%f2(1,ix ,jx+2) = jay_re%f2(1,ix ,jx+2) + wl1(0) * wp1(2) * coeff_re(l)
    jay_re%f2(1,ix+1,jx+2) = jay_re%f2(1,ix+1,jx+2) + wl1(1) * wp1(2) * coeff_re(l)

    ! accumulate j2_re
    jay_re%f2(2,ix-1,jx-1) = jay_re%f2(2,ix-1,jx-1) + wl2(-1) * wp2(-1) * coeff_re(l)
    jay_re%f2(2,ix ,jx-1) = jay_re%f2(2,ix ,jx-1) + wl2(-1) * wp2(0) * coeff_re(l)
    jay_re%f2(2,ix+1,jx-1) = jay_re%f2(2,ix+1,jx-1) + wl2(-1) * wp2(1) * coeff_re(l)
    jay_re%f2(2,ix+2,jx-1) = jay_re%f2(2,ix+2,jx-1) + wl2(-1) * wp2(2) * coeff_re(l)
    jay_re%f2(2,ix-1,jx ) = jay_re%f2(2,ix-1,jx ) + wl2(0) * wp2(-1) * coeff_re(l)
    jay_re%f2(2,ix ,jx ) = jay_re%f2(2,ix ,jx ) + wl2(0) * wp2(0) * coeff_re(l)
    jay_re%f2(2,ix+1,jx ) = jay_re%f2(2,ix+1,jx ) + wl2(0) * wp2(1) * coeff_re(l)
    jay_re%f2(2,ix+2,jx ) = jay_re%f2(2,ix+2,jx ) + wl2(0) * wp2(2) * coeff_re(l)
    jay_re%f2(2,ix-1,jx+1) = jay_re%f2(2,ix-1,jx+1) + wl2(1) * wp2(-1) * coeff_re(l)
    jay_re%f2(2,ix ,jx+1) = jay_re%f2(2,ix ,jx+1) + wl2(1) * wp2(0) * coeff_re(l)
    jay_re%f2(2,ix+1,jx+1) = jay_re%f2(2,ix+1,jx+1) + wl2(1) * wp2(1) * coeff_re(l)
    jay_re%f2(2,ix+2,jx+1) = jay_re%f2(2,ix+2,jx+1) + wl2(1) * wp2(2) * coeff_re(l)

    ! accumulate j1_im
    jay_im%f2(1,ix-1,jx-1) = jay_im%f2(1,ix-1,jx-1) + wl1(-1) * wp1(-1) * coeff_im(l)
    jay_im%f2(1,ix ,jx-1) = jay_im%f2(1,ix ,jx-1) + wl1(0) * wp1(-1) * coeff_im(l)
    jay_im%f2(1,ix+1,jx-1) = jay_im%f2(1,ix+1,jx-1) + wl1(1) * wp1(-1) * coeff_im(l)
    jay_im%f2(1,ix-1,jx ) = jay_im%f2(1,ix-1,jx ) + wl1(-1) * wp1(0) * coeff_im(l)
    jay_im%f2(1,ix ,jx ) = jay_im%f2(1,ix ,jx ) + wl1(0) * wp1(0) * coeff_im(l)
    jay_im%f2(1,ix+1,jx ) = jay_im%f2(1,ix+1,jx ) + wl1(1) * wp1(0) * coeff_im(l)
    jay_im%f2(1,ix-1,jx+1) = jay_im%f2(1,ix-1,jx+1) + wl1(-1) * wp1(1) * coeff_im(l)
    jay_im%f2(1,ix ,jx+1) = jay_im%f2(1,ix ,jx+1) + wl1(0) * wp1(1) * coeff_im(l)
    jay_im%f2(1,ix+1,jx+1) = jay_im%f2(1,ix+1,jx+1) + wl1(1) * wp1(1) * coeff_im(l)
    jay_im%f2(1,ix-1,jx+2) = jay_im%f2(1,ix-1,jx+2) + wl1(-1) * wp1(2) * coeff_im(l)
    jay_im%f2(1,ix ,jx+2) = jay_im%f2(1,ix ,jx+2) + wl1(0) * wp1(2) * coeff_im(l)
    jay_im%f2(1,ix+1,jx+2) = jay_im%f2(1,ix+1,jx+2) + wl1(1) * wp1(2) * coeff_im(l)

    ! accumulate j2_im
    jay_im%f2(2,ix-1,jx-1) = jay_im%f2(2,ix-1,jx-1) + wl2(-1) * wp2(-1) * coeff_im(l)
    jay_im%f2(2,ix ,jx-1) = jay_im%f2(2,ix ,jx-1) + wl2(-1) * wp2(0) * coeff_im(l)
    jay_im%f2(2,ix+1,jx-1) = jay_im%f2(2,ix+1,jx-1) + wl2(-1) * wp2(1) * coeff_im(l)
    jay_im%f2(2,ix+2,jx-1) = jay_im%f2(2,ix+2,jx-1) + wl2(-1) * wp2(2) * coeff_im(l)
    jay_im%f2(2,ix-1,jx ) = jay_im%f2(2,ix-1,jx ) + wl2(0) * wp2(-1) * coeff_im(l)
    jay_im%f2(2,ix ,jx ) = jay_im%f2(2,ix ,jx ) + wl2(0) * wp2(0) * coeff_im(l)
    jay_im%f2(2,ix+1,jx ) = jay_im%f2(2,ix+1,jx ) + wl2(0) * wp2(1) * coeff_im(l)
    jay_im%f2(2,ix+2,jx ) = jay_im%f2(2,ix+2,jx ) + wl2(0) * wp2(2) * coeff_im(l)
    jay_im%f2(2,ix-1,jx+1) = jay_im%f2(2,ix-1,jx+1) + wl2(1) * wp2(-1) * coeff_im(l)
    jay_im%f2(2,ix ,jx+1) = jay_im%f2(2,ix ,jx+1) + wl2(1) * wp2(0) * coeff_im(l)
    jay_im%f2(2,ix+1,jx+1) = jay_im%f2(2,ix+1,jx+1) + wl2(1) * wp2(1) * coeff_im(l)
    jay_im%f2(2,ix+2,jx+1) = jay_im%f2(2,ix+2,jx+1) + wl2(1) * wp2(2) * coeff_im(l)


    ! accumulate j3
    ! charge conserving cylindrical mode current deposition algorithm
    ! real part
    factor_1 = coeff_im_p(l) - coeff_im(l) !-(coeff_im_p(l) - coeff_im(l))
    factor_2 = coeff_im_m(l) - coeff_im(l) !-(coeff_im_m(l) - coeff_im(l))

    jay_re%f2(3,ix-1,jx-1) = jay_re%f2(3,ix-1,jx-1) + qvz * ABS(rmid-1) * (S1x(-1)*S1y(-1) * factor_1 - S0x(-1)*S0y(-1) * factor_2)
    jay_re%f2(3,ix ,jx-1) = jay_re%f2(3,ix ,jx-1) + qvz * ABS(rmid-1) * (S1x(0)*S1y(-1) * factor_1 - S0x(0)*S0y(-1) * factor_2)
    jay_re%f2(3,ix+1,jx-1) = jay_re%f2(3,ix+1,jx-1) + qvz * ABS(rmid-1) * (S1x(1)*S1y(-1) * factor_1 - S0x(1)*S0y(-1) * factor_2)
    jay_re%f2(3,ix+2,jx-1) = jay_re%f2(3,ix+2,jx-1) + qvz * ABS(rmid-1) * (S1x(2)*S1y(-1) * factor_1 - S0x(2)*S0y(-1) * factor_2)
    jay_re%f2(3,ix-1,jx ) = jay_re%f2(3,ix-1,jx ) + qvz * ABS(rmid) * (S1x(-1)*S1y(0) * factor_1 - S0x(-1)*S0y(0) * factor_2)
    jay_re%f2(3,ix ,jx ) = jay_re%f2(3,ix ,jx ) + qvz * ABS(rmid) * (S1x(0)*S1y(0) * factor_1 - S0x(0)*S0y(0) * factor_2)
    jay_re%f2(3,ix+1,jx ) = jay_re%f2(3,ix+1,jx ) + qvz * ABS(rmid) * (S1x(1)*S1y(0) * factor_1 - S0x(1)*S0y(0) * factor_2)
    jay_re%f2(3,ix+2,jx ) = jay_re%f2(3,ix+2,jx ) + qvz * ABS(rmid) * (S1x(2)*S1y(0) * factor_1 - S0x(2)*S0y(0) * factor_2)
    jay_re%f2(3,ix-1,jx+1) = jay_re%f2(3,ix-1,jx+1) + qvz * ABS(rmid+1) * (S1x(-1)*S1y(1) * factor_1 - S0x(-1)*S0y(1) * factor_2)
    jay_re%f2(3,ix ,jx+1) = jay_re%f2(3,ix ,jx+1) + qvz * ABS(rmid+1) * (S1x(0)*S1y(1) * factor_1 - S0x(0)*S0y(1) * factor_2)
    jay_re%f2(3,ix+1,jx+1) = jay_re%f2(3,ix+1,jx+1) + qvz * ABS(rmid+1) * (S1x(1)*S1y(1) * factor_1 - S0x(1)*S0y(1) * factor_2)
    jay_re%f2(3,ix+2,jx+1) = jay_re%f2(3,ix+2,jx+1) + qvz * ABS(rmid+1) * (S1x(2)*S1y(1) * factor_1 - S0x(2)*S0y(1) * factor_2)
    jay_re%f2(3,ix-1,jx+2) = jay_re%f2(3,ix-1,jx+2) + qvz * ABS(rmid+2) * (S1x(-1)*S1y(2) * factor_1 - S0x(-1)*S0y(2) * factor_2)
    jay_re%f2(3,ix ,jx+2) = jay_re%f2(3,ix ,jx+2) + qvz * ABS(rmid+2) * (S1x(0)*S1y(2) * factor_1 - S0x(0)*S0y(2) * factor_2)
    jay_re%f2(3,ix+1,jx+2) = jay_re%f2(3,ix+1,jx+2) + qvz * ABS(rmid+2) * (S1x(1)*S1y(2) * factor_1 - S0x(1)*S0y(2) * factor_2)
    jay_re%f2(3,ix+2,jx+2) = jay_re%f2(3,ix+2,jx+2) + qvz * ABS(rmid+2) * (S1x(2)*S1y(2) * factor_1 - S0x(2)*S0y(2) * factor_2)

    ! imaginary part
    factor_1 = -(coeff_re_p(l) - coeff_re(l)) !coeff_re_p(l) - coeff_re(l)
    factor_2 = -(coeff_re_m(l) - coeff_re(l)) !coeff_re_m(l) - coeff_re(l)

    jay_im%f2(3,ix-1,jx-1) = jay_im%f2(3,ix-1,jx-1) + qvz * ABS(rmid-1) * (S1x(-1)*S1y(-1) * factor_1 - S0x(-1)*S0y(-1) * factor_2)
    jay_im%f2(3,ix ,jx-1) = jay_im%f2(3,ix ,jx-1) + qvz * ABS(rmid-1) * (S1x(0)*S1y(-1) * factor_1 - S0x(0)*S0y(-1) * factor_2)
    jay_im%f2(3,ix+1,jx-1) = jay_im%f2(3,ix+1,jx-1) + qvz * ABS(rmid-1) * (S1x(1)*S1y(-1) * factor_1 - S0x(1)*S0y(-1) * factor_2)
    jay_im%f2(3,ix+2,jx-1) = jay_im%f2(3,ix+2,jx-1) + qvz * ABS(rmid-1) * (S1x(2)*S1y(-1) * factor_1 - S0x(2)*S0y(-1) * factor_2)
    jay_im%f2(3,ix-1,jx ) = jay_im%f2(3,ix-1,jx ) + qvz * ABS(rmid) * (S1x(-1)*S1y(0) * factor_1 - S0x(-1)*S0y(0) * factor_2)
    jay_im%f2(3,ix ,jx ) = jay_im%f2(3,ix ,jx ) + qvz * ABS(rmid) * (S1x(0)*S1y(0) * factor_1 - S0x(0)*S0y(0) * factor_2)
    jay_im%f2(3,ix+1,jx ) = jay_im%f2(3,ix+1,jx ) + qvz * ABS(rmid) * (S1x(1)*S1y(0) * factor_1 - S0x(1)*S0y(0) * factor_2)
    jay_im%f2(3,ix+2,jx ) = jay_im%f2(3,ix+2,jx ) + qvz * ABS(rmid) * (S1x(2)*S1y(0) * factor_1 - S0x(2)*S0y(0) * factor_2)
    jay_im%f2(3,ix-1,jx+1) = jay_im%f2(3,ix-1,jx+1) + qvz * ABS(rmid+1) * (S1x(-1)*S1y(1) * factor_1 - S0x(-1)*S0y(1) * factor_2)
    jay_im%f2(3,ix ,jx+1) = jay_im%f2(3,ix ,jx+1) + qvz * ABS(rmid+1) * (S1x(0)*S1y(1) * factor_1 - S0x(0)*S0y(1) * factor_2)
    jay_im%f2(3,ix+1,jx+1) = jay_im%f2(3,ix+1,jx+1) + qvz * ABS(rmid+1) * (S1x(1)*S1y(1) * factor_1 - S0x(1)*S0y(1) * factor_2)
    jay_im%f2(3,ix+2,jx+1) = jay_im%f2(3,ix+2,jx+1) + qvz * ABS(rmid+1) * (S1x(2)*S1y(1) * factor_1 - S0x(2)*S0y(1) * factor_2)
    jay_im%f2(3,ix-1,jx+2) = jay_im%f2(3,ix-1,jx+2) + qvz * ABS(rmid+2) * (S1x(-1)*S1y(2) * factor_1 - S0x(-1)*S0y(2) * factor_2)
    jay_im%f2(3,ix ,jx+2) = jay_im%f2(3,ix ,jx+2) + qvz * ABS(rmid+2) * (S1x(0)*S1y(2) * factor_1 - S0x(0)*S0y(2) * factor_2)
    jay_im%f2(3,ix+1,jx+2) = jay_im%f2(3,ix+1,jx+2) + qvz * ABS(rmid+2) * (S1x(1)*S1y(2) * factor_1 - S0x(1)*S0y(2) * factor_2)
    jay_im%f2(3,ix+2,jx+2) = jay_im%f2(3,ix+2,jx+2) + qvz * ABS(rmid+2) * (S1x(2)*S1y(2) * factor_1 - S0x(2)*S0y(2) * factor_2)

    ! end of automatic code
  enddo

end subroutine getjr_cyl_m_s3
!-------------------------------------------------------------------------------------------------

end module
