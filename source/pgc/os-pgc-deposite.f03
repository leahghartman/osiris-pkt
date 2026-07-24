!*******************************************************************************
! Field interpolation routines for PGC§
!
!*******************************************************************************


! if __TEMPLATE__ is defined then read the template definition at the end of
! the file
#ifndef __TEMPLATE__


#include "os-preprocess.fpp"
#include "os-config.h"

module m_pgc_deposite

use m_vdf_define
use m_parameters

implicit none

private

interface deposite_pgc
   module procedure deposite_pgc_cell
end interface

public :: deposite_pgc


contains

!-------------------------------------------------------------------------------
! Deposite chi of particles onto cell for cell based particle positions
!-------------------------------------------------------------------------------
subroutine deposite_pgc_cell( chi_buffer, chi, ix, x, np, interpolation )

  implicit none

  real( p_k_part ), dimension(:),   intent(in)  :: chi_buffer

  type( t_vdf ),    intent(inout)               :: chi

  integer,        dimension(:,:),   intent(in)  :: ix
  real(p_k_part), dimension(:,:),   intent(in)  :: x

  integer, intent(in) :: np

  integer, intent(in) :: interpolation


  ! executable statements
	select case (interpolation)
  case(p_linear)
    select case ( p_x_dim )
    case (2)
      call set_pgc_2d_s1( chi_buffer, chi, ix, x, np )
    case (3)
      call set_pgc_3d_s1( chi_buffer, chi, ix, x, np )
    end select

  case(p_quadratic)
    select case ( p_x_dim )
    case (2)
      call set_pgc_2d_s2( chi_buffer, chi, ix, x, np )
    case (3)
      call set_pgc_3d_s2( chi_buffer, chi, ix, x, np )
    end select

  case(p_cubic)
    select case ( p_x_dim )
    case (2)
      call set_pgc_2d_s3( chi_buffer, chi, ix, x, np )
    case (3)
      call set_pgc_3d_s3( chi_buffer, chi, ix, x, np )
    end select

	case(p_quartic)
    select case ( p_x_dim )
    case (2)
      call set_pgc_2d_s4( chi_buffer, chi, ix, x, np )
    case (3)
      call set_pgc_3d_s4( chi_buffer, chi, ix, x, np )
    end select

  case default
    ERROR('Not implemented yet')
    call abort_program( p_err_notimplemented )

	end select


end subroutine deposite_pgc_cell
!-----------------------------------------------------------------------------------------

! ----------------------------------------------------------------------------------------
! Splines
! ----------------------------------------------------------------------------------------

! ----------------------------------------------------------------------------------------
! Linear
! ----------------------------------------------------------------------------------------
subroutine spline_s1( x, s )

  implicit none

  real( p_k_fld ), intent(in) :: x
  real( p_k_fld ), dimension(0:1), intent(out) :: s

  s(0) = 0.5 - x
  s(1) = 0.5 + x

end subroutine spline_s1

! ----------------------------------------------------------------------------------------
! Quadratic
! ----------------------------------------------------------------------------------------
subroutine spline_s2( x, s )

  implicit none

  real( p_k_fld ), intent(in) :: x
  real( p_k_fld ), dimension(-1:1), intent(out) :: s

  real( p_k_fld ) :: t0, t1

  t0 = 0.5 - x
  t1 = 0.5 + x

  s(-1) = 0.5 * t0**2
  s( 0) = 0.5 + t0*t1
  s( 1) = 0.5 * t1**2

end subroutine spline_s2

! ----------------------------------------------------------------------------------------
! Cubic
! ----------------------------------------------------------------------------------------

subroutine spline_s3( x, s )

  implicit none

  real( p_k_fld ), intent(in) :: x
  real( p_k_fld ), dimension(-1:2), intent(out) :: s

  real( p_k_fld ) :: t0, t1, t2, t3

  t0 = 0.5 - x
  t1 = 0.5 + x

  t2 = t0 * t0
  t3 = t1 * t1

  t0 = t0 * t2
  t1 = t1 * t3

  s(-1) = t0/6.
  s( 0) = (0.6666666666666666_p_k_fld - t3 ) + 0.5*t1
  s( 1) = (0.6666666666666666_p_k_fld - t2 ) + 0.5*t0
  s( 2) = t1/6.

end subroutine spline_s3

! ----------------------------------------------------------------------------------------
! Quartic
! ----------------------------------------------------------------------------------------

subroutine spline_s4( x, s )

  implicit none

  real( p_k_fld ), intent(in) :: x
  real( p_k_fld ), dimension(-2:2), intent(out) :: s

  real( p_k_fld ) :: t0, t1, t2

  t0 = 0.5 - x
  t1 = 0.5 + x
  t2 = t0 * t1

  s(-2) = t0**4/24.
  s(-1) = ( 0.25 + t0 * ( 1 + t0 * ( 1.5 + t2 ) ) ) / 6.
  s( 0) = 0.4583333333333333_p_k_fld + t2 * ( 0.5 + 0.25 * t2 )
  s( 1) = ( 0.25 + t1 * ( 1 + t1 * ( 1.5 + t2 ) ) ) / 6.
  s( 2) = t1**4/24.

end subroutine spline_s4

!-----------------------------------------------------------------------------------------
! Signbit function
!-----------------------------------------------------------------------------------------
function signbit(x)

  implicit none
  real(p_k_fld), intent(in) :: x

  integer :: signbit

  if ( x < 0 ) then
    signbit = 1
  else
    signbit = 0
  endif

end function signbit
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
!  Generate specific template functions for linear, quadratic, cubic and quartic
!  interpolation levels.
!-----------------------------------------------------------------------------------------

#define __TEMPLATE__


!********************************** Linear interpolation ********************************

#define SET_PGC_2D   set_pgc_2d_s1
#define SET_PGC_3D   set_pgc_3d_s1

#define SPLINE  spline_s1

! Lower point
#define LP 0
! Upper point
#define UP 1

#include __FILE__

!********************************** Quadratic interpolation ********************************

#define SET_PGC_2D   set_pgc_2d_s2
#define SET_PGC_3D   set_pgc_3d_s2

#define SPLINE  spline_s2

! Lower point
#define LP -1
! Upper point
#define UP 1

#include __FILE__

!********************************** Cubic interpolation **********************************

#define SET_PGC_2D   set_pgc_2d_s3
#define SET_PGC_3D   set_pgc_3d_s3

#define SPLINE  spline_s3

! Lower point
#define LP -1
! Upper point
#define UP 2

#include __FILE__

!********************************** Quartic interpolation ********************************

#define SET_PGC_2D   set_pgc_2d_s4
#define SET_PGC_3D   set_pgc_3d_s4

#define SPLINE  spline_s4

! Lower point
#define LP -2
! Upper point
#define UP 2

#include __FILE__

!****************************************************************************************

end module m_pgc_deposite

#else

!*****************************************************************************************
!
!  Template function definitions for field interpolation
!
!*****************************************************************************************

!-----------------------------------------------------------------------------------------
subroutine SET_PGC_2D( chi_buffer, chi, ix, x, np )
!-----------------------------------------------------------------------------------------

  implicit none

  integer, parameter :: rank = 2

  real( p_k_part ), dimension(:),   intent(in)  :: chi_buffer

  type( t_vdf ),    intent(inout)               :: chi

  integer,        dimension(:,:),   intent(in)  :: ix
  real(p_k_part), dimension(:,:),   intent(in)  :: x

  integer, intent(in) :: np

  ! local variables
  integer :: i1, i2, l, k1, k2
  real(p_k_fld) :: dx1, dx2
  real(p_k_fld), dimension(LP:UP) :: w1, w2

  real(p_k_fld) :: chil


  do l = 1, np

    i1 = ix(1,l)
    i2 = ix(2,l)
    dx1 = real( x(1,l), p_k_fld )
    dx2 = real( x(2,l), p_k_fld )
    chil = real( chi_buffer(l), p_k_fld )

    ! get spline weitghts for x and y
    call SPLINE( dx1, w1 )
    call SPLINE( dx2, w2 )

    ! Deposite chi_buffer on chi field
    do k2 = LP, UP
      do k1 = LP, UP
        chi%f2(1, i1+k1, i2+k2) = chi%f2(1, i1+k1, i2+k2) + chil*w1(k1)*w2(k2)
      enddo
    enddo

  enddo


end subroutine SET_PGC_2D
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine SET_PGC_3D( chi_buffer, chi, ix, x, np )
!-----------------------------------------------------------------------------------------

  implicit none

  ! dummy variables
  integer, parameter :: rank = 3

  real( p_k_part ), dimension(:),   intent(in)  :: chi_buffer

  type( t_vdf ),    intent(inout)               :: chi

  integer,        dimension(:,:),   intent(in)  :: ix
  real(p_k_part), dimension(:,:),   intent(in)  :: x

  integer, intent(in) :: np

  ! local variables
  integer :: i1, i2, i3, l, k1, k2, k3
  real(p_k_fld) :: dx1, dx2, dx3
  real(p_k_fld), dimension(LP:UP) :: w1, w2, w3

  real(p_k_fld) :: chil


  do l = 1, np

    i1 = ix(1,l)
    i2 = ix(2,l)
    i3 = ix(3,l)
    dx1 = real( x(1,l), p_k_fld )
    dx2 = real( x(2,l), p_k_fld )
    dx3 = real( x(3,l), p_k_fld )
    chil = real( chi_buffer(l), p_k_fld )

    ! get spline weights for x,y and z
    call SPLINE( dx1, w1 )
    call SPLINE( dx2, w2 )
    call SPLINE( dx3, w3 )

    ! Deposite chi_buffer on chi field
    do k3 = LP, UP
      do k2 = LP, UP
        do k1 = LP, UP
          chi%f3(1, i1+k1, i2+k2, i3+k3) = &
            chi%f3(1, i1+k1, i2+k2, i3+k3) + chil * w1(k1) * w2(k2) * w3(k3)
        enddo
      enddo
    enddo
  enddo

end subroutine SET_PGC_3D
!-----------------------------------------------------------------------------------------

! Clear template definitions

#undef SET_PGC_2D
#undef SET_PGC_3D

#undef SPLINE

#undef LP
#undef UP

#endif
