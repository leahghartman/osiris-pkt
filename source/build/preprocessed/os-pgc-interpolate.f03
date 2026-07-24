# 1 "pgc/os-pgc-interpolate.f03"
# 1 "<built-in>" 1
# 1 "<built-in>" 3
# 467 "<built-in>" 3
# 1 "<command line>" 1
# 1 "<built-in>" 2
# 1 "pgc/os-pgc-interpolate.f03" 2
!*******************************************************************************
! Field interpolation routines for PGC§
!
!*******************************************************************************


! if __TEMPLATE__ is defined then read the template definition at the end of
! the file



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
# 13 "pgc/os-pgc-interpolate.f03" 2
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
# 14 "pgc/os-pgc-interpolate.f03" 2

module m_pgc_interpolate

use m_vdf_define
use m_parameters

implicit none

private

interface get_pgc
   module procedure get_pgc_cell
end interface

public :: get_pgc


contains

!-------------------------------------------------------------------------------
! Interpolate fields (such as ponderomotive force 'f_p'
! and envelope squared 'a^2') at particle positions for cell based positions.
!-------------------------------------------------------------------------------
subroutine get_pgc_cell( f, a2, fp, a2p, ix, x, np, interpolation )

  implicit none

  type( t_vdf ), intent(in) :: f, a2

  real( p_k_part ), dimension(:,:), intent(out) :: fp
  real( p_k_part ), dimension(:), intent(out) :: a2p

  integer, dimension(:,:), intent(in) :: ix
  real(p_k_part), dimension(:,:), intent(in) :: x

  integer, intent(in) :: np

  integer, intent(in) :: interpolation


  ! executable statements
 select case (interpolation)
  case(p_linear)
    select case ( p_x_dim )
    case (2)
      call get_pgc_2d_s1( fp, a2p, f, a2, ix, x, np )
    case (3)
      call get_pgc_3d_s1( fp, a2p, f, a2, ix, x, np )
    end select

  case(p_quadratic)
    select case ( p_x_dim )
    case (2)
      call get_pgc_2d_s2( fp, a2p, f, a2, ix, x, np )
    case (3)
      call get_pgc_3d_s2( fp, a2p, f, a2, ix, x, np )
    end select

  case(p_cubic)
    select case ( p_x_dim )
    case (2)
      call get_pgc_2d_s3( fp, a2p, f, a2, ix, x, np )
    case (3)
      call get_pgc_3d_s3( fp, a2p, f, a2, ix, x, np )
    end select

 case(p_quartic)
    select case ( p_x_dim )
    case (2)
      call get_pgc_2d_s4( fp, a2p, f, a2, ix, x, np )
    case (3)
      call get_pgc_3d_s4( fp, a2p, f, a2, ix, x, np )
    end select

  case default
    write(err_buf__,*) 'Not implemented yet';call err__("pgc/os-pgc-interpolate.f03",89)
    call abort_program( p_err_notimplemented )

 end select


end subroutine get_pgc_cell
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

subroutine splineh_s1( x, h, s )

  implicit none

  real( p_k_fld ), intent(in) :: x
  integer, intent(in) :: h
  real( p_k_fld ), dimension(0:1), intent(out) :: s

  s(0) = ( 1 - h ) - x
  s(1) = ( h ) + x

end subroutine splineh_s1

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

subroutine splineh_s2( x, h, s )

  implicit none

  real( p_k_fld ), intent(in) :: x
  integer, intent(in) :: h
  real( p_k_fld ), dimension(-1:1), intent(out) :: s

  real( p_k_fld ) :: t0, t1

  t0 = ( 1 - h ) - x
  t1 = ( h ) + x

  s(-1) = 0.5 * t0**2
  s( 0) = 0.5 + t0*t1
  s( 1) = 0.5 * t1**2

end subroutine splineh_s2

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

subroutine splineh_s3( x, h, s )

  implicit none

  real( p_k_fld ), intent(in) :: x
  integer, intent(in) :: h
  real( p_k_fld ), dimension(-1:2), intent(out) :: s

  real( p_k_fld ) :: t0, t1, t2, t3

  t0 = ( 1 - h ) - x
  t1 = ( h ) + x

  t2 = t0 * t0
  t3 = t1 * t1

  t0 = t0 * t2
  t1 = t1 * t3

  s(-1) = t0/6.
  s( 0) = (0.6666666666666666_p_k_fld - t3 ) + 0.5*t1
  s( 1) = (0.6666666666666666_p_k_fld - t2 ) + 0.5*t0
  s( 2) = t1/6.

end subroutine splineh_s3

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

subroutine splineh_s4( x, h, s )

  implicit none

  real( p_k_fld ), intent(in) :: x
  integer, intent(in) :: h
  real( p_k_fld ), dimension(-2:2), intent(out) :: s

  real( p_k_fld ) :: t0, t1, t2

  t0 = ( 1 - h ) - x
  t1 = ( h ) + x
  t2 = t0 * t1

  s(-2) = t0**4/24.
  s(-1) = ( 0.25 + t0 * ( 1 + t0 * ( 1.5 + t2 ) ) ) / 6.
  s( 0) = 0.4583333333333333_p_k_fld + t2 * ( 0.5 + 0.25 * t2 )
  s( 1) = ( 0.25 + t1 * ( 1 + t1 * ( 1.5 + t2 ) ) ) / 6.
  s( 2) = t1**4/24.

end subroutine splineh_s4


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
! Generate specific template functions for linear, quadratic, cubic and quartic
! interpolation levels.
!-----------------------------------------------------------------------------------------




!********************************** Linear interpolation ********************************







! Lower point

! Upper point


# 1 "pgc/os-pgc-interpolate.f03" 1
!*******************************************************************************
! Field interpolation routines for PGC§
!
!*******************************************************************************


! if is defined then read the template definition at the end of
! the file
# 367 "pgc/os-pgc-interpolate.f03"
!*****************************************************************************************
!
! Template function definitions for field interpolation
!
!*****************************************************************************************

!-----------------------------------------------------------------------------------------
subroutine get_pgc_2d_s1( fp, a2p, f, a2, ix, x, np )
!-----------------------------------------------------------------------------------------

  implicit none

  integer, parameter :: rank = 2

  real(p_k_part), dimension(:,:), intent( out ) :: fp
  real(p_k_part), dimension(:), intent( out ) :: a2p

  type( t_vdf ), intent(in) :: f, a2

  integer, dimension(:,:), intent(in) :: ix
  real(p_k_part), dimension(:,:), intent(in) :: x
  integer, intent(in) :: np

  ! local variables
  integer :: i1, i2, i1h, i2h, l, k1, k2, h1, h2
  real(p_k_fld) :: dx1, dx2
  real(p_k_fld), dimension(0:1) :: w1, w1h, w2, w2h

  real(p_k_fld) :: f1, f2, f1line, f2line

  do l = 1, np

    i1 = ix(1,l)
    i2 = ix(2,l)
    dx1 = real( x(1,l), p_k_fld )
    dx2 = real( x(2,l), p_k_fld )

    h1 = signbit(dx1)
    h2 = signbit(dx2)

    i1h = i1 - h1
    i2h = i2 - h2

    ! get spline weitghts for x and y
    call spline_s1( dx1, w1 )
    call splineh_s1( dx1, h1, w1h )

    call spline_s1( dx2, w2 )
    call splineh_s1( dx2, h2, w2h )

    ! Interpolate ponderomotive force - vector field
    f1 = 0
    f2 = 0

    do k2 = 0, 1
      f1line = 0
      f2line = 0

      do k1 = 0, 1
        f1line = f1line + f%f2(1,i1h + k1, i2 + k2) * w1h(k1)
        f2line = f2line + f%f2(2,i1 + k1, i2h + k2) * w1(k1)
      enddo

      f1 = f1 + f1line * w2(k2)
      f2 = f2 + f2line * w2h(k2)
    enddo

    fp( 1, l ) = f1
    fp( 2, l ) = f2

    ! Interpolate a^2 - scalar field
    f1 = 0
    f2 = 0

    do k2 = 0, 1
      f1line = 0
      f2line = 0

      do k1 = 0, 1
        f1line = f1line + a2%f2(1,i1 + k1, i2 + k2) * w1(k1)
      enddo

      f1 = f1 + f1line * w2(k2)
    enddo

    a2p( l ) = f1

  enddo


end subroutine get_pgc_2d_s1
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine get_pgc_3d_s1( fp, a2p, f, a2, ix, x, np )
!-----------------------------------------------------------------------------------------

  implicit none

  ! dummy variables
  integer, parameter :: rank = 3

  real(p_k_part), dimension(:,:), intent( out ) :: fp
  real(p_k_part), dimension(:), intent( out ) :: a2p

  type( t_vdf ), intent(in) :: f, a2

  integer, dimension(:,:), intent(in) :: ix
  real(p_k_part), dimension(:,:), intent(in) :: x
  integer, intent(in) :: np

  ! local variables
  real(p_k_fld) :: dx1, dx2, dx3
  integer :: h1, h2, h3
  real(p_k_fld), dimension(0:1) :: w1, w1h, w2, w2h, w3, w3h
  real(p_k_fld) :: f1, f2, f3
  real(p_k_fld) :: f1line, f2line, f3line
  real(p_k_fld) :: f1plane, f2plane, f3plane

  integer :: i1, i2, i3, i1h, i2h, i3h, l
  integer :: k1, k2, k3

  do l = 1, np

    i1 = ix(1,l)
    i2 = ix(2,l)
    i3 = ix(3,l)
    dx1 = real( x(1,l), p_k_fld )
    dx2 = real( x(2,l), p_k_fld )
    dx3 = real( x(3,l), p_k_fld )

    h1 = signbit(dx1)
    h2 = signbit(dx2)
    h3 = signbit(dx3)

    i1h = i1 - h1
    i2h = i2 - h2
    i3h = i3 - h3

    ! get spline weights for x,y and z
    call spline_s1( dx1, w1 )
    call splineh_s1( dx1, h1, w1h )

    call spline_s1( dx2, w2 )
    call splineh_s1( dx2, h2, w2h )

    call spline_s1( dx3, w3 )
    call splineh_s1( dx3, h3, w3h )

    ! Interpolate ponderomotive force - vector field
    f1 = 0
    f2 = 0
    f3 = 0

    do k3 = 0, 1
      f1plane = 0
      f2plane = 0
      f3plane = 0

      do k2 = 0, 1
        f1line = 0
        f2line = 0
        f3line = 0

        do k1 = 0, 1
          f1line = f1line + f%f3(1,i1h + k1, i2 + k2, i3 + k3) * w1h(k1)
          f2line = f2line + f%f3(2,i1 + k1, i2h + k2, i3 + k3) * w1(k1)
          f3line = f3line + f%f3(3,i1 + k1, i2 + k2, i3h + k3) * w1(k1)
        enddo

        f1plane = f1plane + f1line * w2(k2)
        f2plane = f2plane + f2line * w2h(k2)
        f3plane = f3plane + f3line * w2(k2)
      enddo

      f1 = f1 + f1plane * w3(k3)
      f2 = f2 + f2plane * w3(k3)
      f3 = f3 + f3plane * w3h(k3)
    enddo

    fp( 1, l ) = f1
    fp( 2, l ) = f2
    fp( 3, l ) = f3

    ! Interpolate a^2 - scalar field
    f1 = 0

    do k3 = 0, 1
      f1plane = 0

      do k2 = 0, 1
        f1line = 0

        do k1 = 0, 1
          f1line = f1line + a2%f3(1,i1 + k1, i2 + k2, i3 + k3) * w1(k1)
        enddo

        f1plane = f1plane + f1line * w2(k2)
      enddo

      f1 = f1 + f1plane * w3(k3)
    enddo

    a2p( l ) = f1

  enddo

end subroutine get_pgc_3d_s1
!-----------------------------------------------------------------------------------------

! Clear template definitions
# 315 "pgc/os-pgc-interpolate.f03" 2

!********************************** Quadratic interpolation ********************************







! Lower point

! Upper point


# 1 "pgc/os-pgc-interpolate.f03" 1
!*******************************************************************************
! Field interpolation routines for PGC§
!
!*******************************************************************************


! if is defined then read the template definition at the end of
! the file
# 367 "pgc/os-pgc-interpolate.f03"
!*****************************************************************************************
!
! Template function definitions for field interpolation
!
!*****************************************************************************************

!-----------------------------------------------------------------------------------------
subroutine get_pgc_2d_s2( fp, a2p, f, a2, ix, x, np )
!-----------------------------------------------------------------------------------------

  implicit none

  integer, parameter :: rank = 2

  real(p_k_part), dimension(:,:), intent( out ) :: fp
  real(p_k_part), dimension(:), intent( out ) :: a2p

  type( t_vdf ), intent(in) :: f, a2

  integer, dimension(:,:), intent(in) :: ix
  real(p_k_part), dimension(:,:), intent(in) :: x
  integer, intent(in) :: np

  ! local variables
  integer :: i1, i2, i1h, i2h, l, k1, k2, h1, h2
  real(p_k_fld) :: dx1, dx2
  real(p_k_fld), dimension(-1:1) :: w1, w1h, w2, w2h

  real(p_k_fld) :: f1, f2, f1line, f2line

  do l = 1, np

    i1 = ix(1,l)
    i2 = ix(2,l)
    dx1 = real( x(1,l), p_k_fld )
    dx2 = real( x(2,l), p_k_fld )

    h1 = signbit(dx1)
    h2 = signbit(dx2)

    i1h = i1 - h1
    i2h = i2 - h2

    ! get spline weitghts for x and y
    call spline_s2( dx1, w1 )
    call splineh_s2( dx1, h1, w1h )

    call spline_s2( dx2, w2 )
    call splineh_s2( dx2, h2, w2h )

    ! Interpolate ponderomotive force - vector field
    f1 = 0
    f2 = 0

    do k2 = -1, 1
      f1line = 0
      f2line = 0

      do k1 = -1, 1
        f1line = f1line + f%f2(1,i1h + k1, i2 + k2) * w1h(k1)
        f2line = f2line + f%f2(2,i1 + k1, i2h + k2) * w1(k1)
      enddo

      f1 = f1 + f1line * w2(k2)
      f2 = f2 + f2line * w2h(k2)
    enddo

    fp( 1, l ) = f1
    fp( 2, l ) = f2

    ! Interpolate a^2 - scalar field
    f1 = 0
    f2 = 0

    do k2 = -1, 1
      f1line = 0
      f2line = 0

      do k1 = -1, 1
        f1line = f1line + a2%f2(1,i1 + k1, i2 + k2) * w1(k1)
      enddo

      f1 = f1 + f1line * w2(k2)
    enddo

    a2p( l ) = f1

  enddo


end subroutine get_pgc_2d_s2
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine get_pgc_3d_s2( fp, a2p, f, a2, ix, x, np )
!-----------------------------------------------------------------------------------------

  implicit none

  ! dummy variables
  integer, parameter :: rank = 3

  real(p_k_part), dimension(:,:), intent( out ) :: fp
  real(p_k_part), dimension(:), intent( out ) :: a2p

  type( t_vdf ), intent(in) :: f, a2

  integer, dimension(:,:), intent(in) :: ix
  real(p_k_part), dimension(:,:), intent(in) :: x
  integer, intent(in) :: np

  ! local variables
  real(p_k_fld) :: dx1, dx2, dx3
  integer :: h1, h2, h3
  real(p_k_fld), dimension(-1:1) :: w1, w1h, w2, w2h, w3, w3h
  real(p_k_fld) :: f1, f2, f3
  real(p_k_fld) :: f1line, f2line, f3line
  real(p_k_fld) :: f1plane, f2plane, f3plane

  integer :: i1, i2, i3, i1h, i2h, i3h, l
  integer :: k1, k2, k3

  do l = 1, np

    i1 = ix(1,l)
    i2 = ix(2,l)
    i3 = ix(3,l)
    dx1 = real( x(1,l), p_k_fld )
    dx2 = real( x(2,l), p_k_fld )
    dx3 = real( x(3,l), p_k_fld )

    h1 = signbit(dx1)
    h2 = signbit(dx2)
    h3 = signbit(dx3)

    i1h = i1 - h1
    i2h = i2 - h2
    i3h = i3 - h3

    ! get spline weights for x,y and z
    call spline_s2( dx1, w1 )
    call splineh_s2( dx1, h1, w1h )

    call spline_s2( dx2, w2 )
    call splineh_s2( dx2, h2, w2h )

    call spline_s2( dx3, w3 )
    call splineh_s2( dx3, h3, w3h )

    ! Interpolate ponderomotive force - vector field
    f1 = 0
    f2 = 0
    f3 = 0

    do k3 = -1, 1
      f1plane = 0
      f2plane = 0
      f3plane = 0

      do k2 = -1, 1
        f1line = 0
        f2line = 0
        f3line = 0

        do k1 = -1, 1
          f1line = f1line + f%f3(1,i1h + k1, i2 + k2, i3 + k3) * w1h(k1)
          f2line = f2line + f%f3(2,i1 + k1, i2h + k2, i3 + k3) * w1(k1)
          f3line = f3line + f%f3(3,i1 + k1, i2 + k2, i3h + k3) * w1(k1)
        enddo

        f1plane = f1plane + f1line * w2(k2)
        f2plane = f2plane + f2line * w2h(k2)
        f3plane = f3plane + f3line * w2(k2)
      enddo

      f1 = f1 + f1plane * w3(k3)
      f2 = f2 + f2plane * w3(k3)
      f3 = f3 + f3plane * w3h(k3)
    enddo

    fp( 1, l ) = f1
    fp( 2, l ) = f2
    fp( 3, l ) = f3

    ! Interpolate a^2 - scalar field
    f1 = 0

    do k3 = -1, 1
      f1plane = 0

      do k2 = -1, 1
        f1line = 0

        do k1 = -1, 1
          f1line = f1line + a2%f3(1,i1 + k1, i2 + k2, i3 + k3) * w1(k1)
        enddo

        f1plane = f1plane + f1line * w2(k2)
      enddo

      f1 = f1 + f1plane * w3(k3)
    enddo

    a2p( l ) = f1

  enddo

end subroutine get_pgc_3d_s2
!-----------------------------------------------------------------------------------------

! Clear template definitions
# 330 "pgc/os-pgc-interpolate.f03" 2

!********************************** Cubic interpolation **********************************







! Lower point

! Upper point


# 1 "pgc/os-pgc-interpolate.f03" 1
!*******************************************************************************
! Field interpolation routines for PGC§
!
!*******************************************************************************


! if is defined then read the template definition at the end of
! the file
# 367 "pgc/os-pgc-interpolate.f03"
!*****************************************************************************************
!
! Template function definitions for field interpolation
!
!*****************************************************************************************

!-----------------------------------------------------------------------------------------
subroutine get_pgc_2d_s3( fp, a2p, f, a2, ix, x, np )
!-----------------------------------------------------------------------------------------

  implicit none

  integer, parameter :: rank = 2

  real(p_k_part), dimension(:,:), intent( out ) :: fp
  real(p_k_part), dimension(:), intent( out ) :: a2p

  type( t_vdf ), intent(in) :: f, a2

  integer, dimension(:,:), intent(in) :: ix
  real(p_k_part), dimension(:,:), intent(in) :: x
  integer, intent(in) :: np

  ! local variables
  integer :: i1, i2, i1h, i2h, l, k1, k2, h1, h2
  real(p_k_fld) :: dx1, dx2
  real(p_k_fld), dimension(-1:2) :: w1, w1h, w2, w2h

  real(p_k_fld) :: f1, f2, f1line, f2line

  do l = 1, np

    i1 = ix(1,l)
    i2 = ix(2,l)
    dx1 = real( x(1,l), p_k_fld )
    dx2 = real( x(2,l), p_k_fld )

    h1 = signbit(dx1)
    h2 = signbit(dx2)

    i1h = i1 - h1
    i2h = i2 - h2

    ! get spline weitghts for x and y
    call spline_s3( dx1, w1 )
    call splineh_s3( dx1, h1, w1h )

    call spline_s3( dx2, w2 )
    call splineh_s3( dx2, h2, w2h )

    ! Interpolate ponderomotive force - vector field
    f1 = 0
    f2 = 0

    do k2 = -1, 2
      f1line = 0
      f2line = 0

      do k1 = -1, 2
        f1line = f1line + f%f2(1,i1h + k1, i2 + k2) * w1h(k1)
        f2line = f2line + f%f2(2,i1 + k1, i2h + k2) * w1(k1)
      enddo

      f1 = f1 + f1line * w2(k2)
      f2 = f2 + f2line * w2h(k2)
    enddo

    fp( 1, l ) = f1
    fp( 2, l ) = f2

    ! Interpolate a^2 - scalar field
    f1 = 0
    f2 = 0

    do k2 = -1, 2
      f1line = 0
      f2line = 0

      do k1 = -1, 2
        f1line = f1line + a2%f2(1,i1 + k1, i2 + k2) * w1(k1)
      enddo

      f1 = f1 + f1line * w2(k2)
    enddo

    a2p( l ) = f1

  enddo


end subroutine get_pgc_2d_s3
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine get_pgc_3d_s3( fp, a2p, f, a2, ix, x, np )
!-----------------------------------------------------------------------------------------

  implicit none

  ! dummy variables
  integer, parameter :: rank = 3

  real(p_k_part), dimension(:,:), intent( out ) :: fp
  real(p_k_part), dimension(:), intent( out ) :: a2p

  type( t_vdf ), intent(in) :: f, a2

  integer, dimension(:,:), intent(in) :: ix
  real(p_k_part), dimension(:,:), intent(in) :: x
  integer, intent(in) :: np

  ! local variables
  real(p_k_fld) :: dx1, dx2, dx3
  integer :: h1, h2, h3
  real(p_k_fld), dimension(-1:2) :: w1, w1h, w2, w2h, w3, w3h
  real(p_k_fld) :: f1, f2, f3
  real(p_k_fld) :: f1line, f2line, f3line
  real(p_k_fld) :: f1plane, f2plane, f3plane

  integer :: i1, i2, i3, i1h, i2h, i3h, l
  integer :: k1, k2, k3

  do l = 1, np

    i1 = ix(1,l)
    i2 = ix(2,l)
    i3 = ix(3,l)
    dx1 = real( x(1,l), p_k_fld )
    dx2 = real( x(2,l), p_k_fld )
    dx3 = real( x(3,l), p_k_fld )

    h1 = signbit(dx1)
    h2 = signbit(dx2)
    h3 = signbit(dx3)

    i1h = i1 - h1
    i2h = i2 - h2
    i3h = i3 - h3

    ! get spline weights for x,y and z
    call spline_s3( dx1, w1 )
    call splineh_s3( dx1, h1, w1h )

    call spline_s3( dx2, w2 )
    call splineh_s3( dx2, h2, w2h )

    call spline_s3( dx3, w3 )
    call splineh_s3( dx3, h3, w3h )

    ! Interpolate ponderomotive force - vector field
    f1 = 0
    f2 = 0
    f3 = 0

    do k3 = -1, 2
      f1plane = 0
      f2plane = 0
      f3plane = 0

      do k2 = -1, 2
        f1line = 0
        f2line = 0
        f3line = 0

        do k1 = -1, 2
          f1line = f1line + f%f3(1,i1h + k1, i2 + k2, i3 + k3) * w1h(k1)
          f2line = f2line + f%f3(2,i1 + k1, i2h + k2, i3 + k3) * w1(k1)
          f3line = f3line + f%f3(3,i1 + k1, i2 + k2, i3h + k3) * w1(k1)
        enddo

        f1plane = f1plane + f1line * w2(k2)
        f2plane = f2plane + f2line * w2h(k2)
        f3plane = f3plane + f3line * w2(k2)
      enddo

      f1 = f1 + f1plane * w3(k3)
      f2 = f2 + f2plane * w3(k3)
      f3 = f3 + f3plane * w3h(k3)
    enddo

    fp( 1, l ) = f1
    fp( 2, l ) = f2
    fp( 3, l ) = f3

    ! Interpolate a^2 - scalar field
    f1 = 0

    do k3 = -1, 2
      f1plane = 0

      do k2 = -1, 2
        f1line = 0

        do k1 = -1, 2
          f1line = f1line + a2%f3(1,i1 + k1, i2 + k2, i3 + k3) * w1(k1)
        enddo

        f1plane = f1plane + f1line * w2(k2)
      enddo

      f1 = f1 + f1plane * w3(k3)
    enddo

    a2p( l ) = f1

  enddo

end subroutine get_pgc_3d_s3
!-----------------------------------------------------------------------------------------

! Clear template definitions
# 345 "pgc/os-pgc-interpolate.f03" 2

!********************************** Quartic interpolation ********************************







! Lower point

! Upper point


# 1 "pgc/os-pgc-interpolate.f03" 1
!*******************************************************************************
! Field interpolation routines for PGC§
!
!*******************************************************************************


! if is defined then read the template definition at the end of
! the file
# 367 "pgc/os-pgc-interpolate.f03"
!*****************************************************************************************
!
! Template function definitions for field interpolation
!
!*****************************************************************************************

!-----------------------------------------------------------------------------------------
subroutine get_pgc_2d_s4( fp, a2p, f, a2, ix, x, np )
!-----------------------------------------------------------------------------------------

  implicit none

  integer, parameter :: rank = 2

  real(p_k_part), dimension(:,:), intent( out ) :: fp
  real(p_k_part), dimension(:), intent( out ) :: a2p

  type( t_vdf ), intent(in) :: f, a2

  integer, dimension(:,:), intent(in) :: ix
  real(p_k_part), dimension(:,:), intent(in) :: x
  integer, intent(in) :: np

  ! local variables
  integer :: i1, i2, i1h, i2h, l, k1, k2, h1, h2
  real(p_k_fld) :: dx1, dx2
  real(p_k_fld), dimension(-2:2) :: w1, w1h, w2, w2h

  real(p_k_fld) :: f1, f2, f1line, f2line

  do l = 1, np

    i1 = ix(1,l)
    i2 = ix(2,l)
    dx1 = real( x(1,l), p_k_fld )
    dx2 = real( x(2,l), p_k_fld )

    h1 = signbit(dx1)
    h2 = signbit(dx2)

    i1h = i1 - h1
    i2h = i2 - h2

    ! get spline weitghts for x and y
    call spline_s4( dx1, w1 )
    call splineh_s4( dx1, h1, w1h )

    call spline_s4( dx2, w2 )
    call splineh_s4( dx2, h2, w2h )

    ! Interpolate ponderomotive force - vector field
    f1 = 0
    f2 = 0

    do k2 = -2, 2
      f1line = 0
      f2line = 0

      do k1 = -2, 2
        f1line = f1line + f%f2(1,i1h + k1, i2 + k2) * w1h(k1)
        f2line = f2line + f%f2(2,i1 + k1, i2h + k2) * w1(k1)
      enddo

      f1 = f1 + f1line * w2(k2)
      f2 = f2 + f2line * w2h(k2)
    enddo

    fp( 1, l ) = f1
    fp( 2, l ) = f2

    ! Interpolate a^2 - scalar field
    f1 = 0
    f2 = 0

    do k2 = -2, 2
      f1line = 0
      f2line = 0

      do k1 = -2, 2
        f1line = f1line + a2%f2(1,i1 + k1, i2 + k2) * w1(k1)
      enddo

      f1 = f1 + f1line * w2(k2)
    enddo

    a2p( l ) = f1

  enddo


end subroutine get_pgc_2d_s4
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine get_pgc_3d_s4( fp, a2p, f, a2, ix, x, np )
!-----------------------------------------------------------------------------------------

  implicit none

  ! dummy variables
  integer, parameter :: rank = 3

  real(p_k_part), dimension(:,:), intent( out ) :: fp
  real(p_k_part), dimension(:), intent( out ) :: a2p

  type( t_vdf ), intent(in) :: f, a2

  integer, dimension(:,:), intent(in) :: ix
  real(p_k_part), dimension(:,:), intent(in) :: x
  integer, intent(in) :: np

  ! local variables
  real(p_k_fld) :: dx1, dx2, dx3
  integer :: h1, h2, h3
  real(p_k_fld), dimension(-2:2) :: w1, w1h, w2, w2h, w3, w3h
  real(p_k_fld) :: f1, f2, f3
  real(p_k_fld) :: f1line, f2line, f3line
  real(p_k_fld) :: f1plane, f2plane, f3plane

  integer :: i1, i2, i3, i1h, i2h, i3h, l
  integer :: k1, k2, k3

  do l = 1, np

    i1 = ix(1,l)
    i2 = ix(2,l)
    i3 = ix(3,l)
    dx1 = real( x(1,l), p_k_fld )
    dx2 = real( x(2,l), p_k_fld )
    dx3 = real( x(3,l), p_k_fld )

    h1 = signbit(dx1)
    h2 = signbit(dx2)
    h3 = signbit(dx3)

    i1h = i1 - h1
    i2h = i2 - h2
    i3h = i3 - h3

    ! get spline weights for x,y and z
    call spline_s4( dx1, w1 )
    call splineh_s4( dx1, h1, w1h )

    call spline_s4( dx2, w2 )
    call splineh_s4( dx2, h2, w2h )

    call spline_s4( dx3, w3 )
    call splineh_s4( dx3, h3, w3h )

    ! Interpolate ponderomotive force - vector field
    f1 = 0
    f2 = 0
    f3 = 0

    do k3 = -2, 2
      f1plane = 0
      f2plane = 0
      f3plane = 0

      do k2 = -2, 2
        f1line = 0
        f2line = 0
        f3line = 0

        do k1 = -2, 2
          f1line = f1line + f%f3(1,i1h + k1, i2 + k2, i3 + k3) * w1h(k1)
          f2line = f2line + f%f3(2,i1 + k1, i2h + k2, i3 + k3) * w1(k1)
          f3line = f3line + f%f3(3,i1 + k1, i2 + k2, i3h + k3) * w1(k1)
        enddo

        f1plane = f1plane + f1line * w2(k2)
        f2plane = f2plane + f2line * w2h(k2)
        f3plane = f3plane + f3line * w2(k2)
      enddo

      f1 = f1 + f1plane * w3(k3)
      f2 = f2 + f2plane * w3(k3)
      f3 = f3 + f3plane * w3h(k3)
    enddo

    fp( 1, l ) = f1
    fp( 2, l ) = f2
    fp( 3, l ) = f3

    ! Interpolate a^2 - scalar field
    f1 = 0

    do k3 = -2, 2
      f1plane = 0

      do k2 = -2, 2
        f1line = 0

        do k1 = -2, 2
          f1line = f1line + a2%f3(1,i1 + k1, i2 + k2, i3 + k3) * w1(k1)
        enddo

        f1plane = f1plane + f1line * w2(k2)
      enddo

      f1 = f1 + f1plane * w3(k3)
    enddo

    a2p( l ) = f1

  enddo

end subroutine get_pgc_3d_s4
!-----------------------------------------------------------------------------------------

! Clear template definitions
# 360 "pgc/os-pgc-interpolate.f03" 2

!****************************************************************************************

end module m_pgc_interpolate
