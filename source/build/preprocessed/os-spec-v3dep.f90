# 1 "spec/os-spec-v3dep.f90"
# 1 "<built-in>" 1
# 1 "<built-in>" 3
# 467 "<built-in>" 3
# 1 "<command line>" 1
# 1 "<built-in>" 2
# 1 "spec/os-spec-v3dep.f90" 2
!*****************************************************************************************
! m_species_v3dep module
!
! Handles deposition of a vector 3 quantity at the cell corner
!*****************************************************************************************

! if __TEMPLATE__ is defined then read the template definition at the end of the file



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
# 12 "spec/os-spec-v3dep.f90" 2
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
# 13 "spec/os-spec-v3dep.f90" 2

module m_species_v3dep

use m_species_define

use m_system
use m_parameters
use m_vdf_define

private

interface deposit_v3
  module procedure deposit_v3_spec
end interface

public :: deposit_v3

! -----------------------------------------------------------------------------
contains

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


!*****************************************************************************************
!* Interface *
!*****************************************************************************************

!-----------------------------------------------------------------------------------------
subroutine deposit_v3_spec( this, f, i1, i2, vq )

  implicit none

  ! dummy variables
  class( t_species ), intent(in) :: this
  type( t_vdf ), intent(inout) :: f
  integer, intent(in) :: i1, i2

  real(p_k_part), dimension(:,:) :: vq

  integer :: np

  ! number of particles to deposit
  np = i2 - i1 + 1

  ! deposit given density
  select case ( this%interpolation )

     case (p_linear)
        select case ( p_x_dim )
         case (1)
           call deposit_v3_1d_s1( f, this%ix(:,i1:i2), this%x(:,i1:i2), vq, np )
         case (2)
           call deposit_v3_2d_s1( f, this%ix(:,i1:i2), this%x(:,i1:i2), vq, np )
         case (3)
           call deposit_v3_3d_s1( f, this%ix(:,i1:i2), this%x(:,i1:i2), vq, np )
        end select

     case (p_quadratic)
        select case ( p_x_dim )
         case (1)
           call deposit_v3_1d_s2( f, this%ix(:,i1:i2), this%x(:,i1:i2), vq, np )
         case (2)
           call deposit_v3_2d_s2( f, this%ix(:,i1:i2), this%x(:,i1:i2), vq, np )
         case (3)
           call deposit_v3_3d_s2( f, this%ix(:,i1:i2), this%x(:,i1:i2), vq, np )
        end select

     case (p_cubic)
        select case ( p_x_dim )
         case (1)
           call deposit_v3_1d_s3( f, this%ix(:,i1:i2), this%x(:,i1:i2), vq, np )
         case (2)
           call deposit_v3_2d_s3( f, this%ix(:,i1:i2), this%x(:,i1:i2), vq, np )
         case (3)
           call deposit_v3_3d_s3( f, this%ix(:,i1:i2), this%x(:,i1:i2), vq, np )
        end select

     case (p_quartic)
        select case ( p_x_dim )
         case (1)
           call deposit_v3_1d_s4( f, this%ix(:,i1:i2), this%x(:,i1:i2), vq, np )
         case (2)
           call deposit_v3_2d_s4( f, this%ix(:,i1:i2), this%x(:,i1:i2), vq, np )
         case (3)
           call deposit_v3_3d_s4( f, this%ix(:,i1:i2), this%x(:,i1:i2), vq, np )
        end select

     case default
        write(err_buf__,*) 'Not implemented';call err__("spec/os-spec-v3dep.f90",193)
        call abort_program( p_err_notimplemented )

  end select



end subroutine deposit_v3_spec
!-------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
! Generate specific template functions for linear, quadratic, cubic and quartic
! interpolation levels.
!-----------------------------------------------------------------------------------------



!********************************** Linear interpolation ********************************







! Lower point

! Upper point


# 1 "spec/os-spec-v3dep.f90" 1
!*****************************************************************************************
! m_species_v3dep module
!
! Handles deposition of a vector 3 quantity at the cell corner
!*****************************************************************************************

! if is defined then read the template definition at the end of the file
# 279 "spec/os-spec-v3dep.f90"
!*****************************************************************************************
!
! Template function definitions for v3 deposition
!
!*****************************************************************************************


!*****************************************************************************************
!* Linear Interpolation *
!*****************************************************************************************

! ----------------------------------------------------------------------------------------
subroutine deposit_v3_1d_s1( f, ix, x, vq, np )

  implicit none

  integer, parameter :: rank = 1

  type( t_vdf ), intent(inout) :: f
  integer, dimension(:,:), intent(in) :: ix
  real(p_k_part), dimension(:,:), intent(in) :: x
  real(p_k_part), dimension(:,:), intent(in) :: vq
  integer, intent(in) :: np

  ! local variables
  integer :: l, i1, k1
  real( p_k_fld ) :: lq1, lq2, lq3, dx1
  real( p_k_fld ), dimension(0:1) :: w1

  ! quantity is defined in the lower corner of the cell

  do l = 1, np

     i1 = ix(1,l)
     dx1 = real( x(1,l), p_k_fld )
     lq1 = real( vq(1,l) , p_k_fld )
     lq2 = real( vq(2,l) , p_k_fld )
     lq3 = real( vq(3,l) , p_k_fld )

     ! get spline weitghts for x
     call spline_s1( dx1, w1 )

     ! Deposit Charge
     do k1 = 0, 1
       f%f1(1,i1+k1) = f%f1(1,i1+k1) + lq1 * w1(k1)
       f%f1(2,i1+k1) = f%f1(2,i1+k1) + lq2 * w1(k1)
       f%f1(3,i1+k1) = f%f1(3,i1+k1) + lq3 * w1(k1)
     enddo

  enddo

end subroutine deposit_v3_1d_s1
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine deposit_v3_2d_s1( f, ix, x, vq, np )

  implicit none

  integer, parameter :: rank = 2

  type( t_vdf ), intent(inout) :: f
  integer, dimension(:,:), intent(in) :: ix
  real(p_k_part), dimension(:,:), intent(in) :: x
  real(p_k_part), dimension(:,:), intent(in) :: vq
  integer, intent(in) :: np

  ! local variables
  integer :: l, i1, i2, k1, k2
  real(p_k_fld) :: dx1, dx2, lq1, lq2, lq3, w12
  real( p_k_fld ), dimension(0:1) :: w1, w2

  ! quantity is defined in the lower corner of the cell

  do l = 1, np

    i1 = ix(1,l)
    i2 = ix(2,l)
    dx1 = real( x(1,l), p_k_fld )
    dx2 = real( x(2,l), p_k_fld )
    lq1 = real( vq(1,l) , p_k_fld )
    lq2 = real( vq(2,l) , p_k_fld )
    lq3 = real( vq(3,l) , p_k_fld )

    ! get spline weitghts for x and y
    call spline_s1( dx1, w1 )
    call spline_s1( dx2, w2 )

    ! Deposit vector3 quantity
    do k2 = 0, 1
      do k1 = 0, 1
        w12 = w1(k1)* w2(k2)
        f%f2(1,i1+k1,i2+k2) = f%f2(1,i1+k1,i2+k2) + lq1 * w12
        f%f2(2,i1+k1,i2+k2) = f%f2(2,i1+k1,i2+k2) + lq2 * w12
        f%f2(3,i1+k1,i2+k2) = f%f2(3,i1+k1,i2+k2) + lq3 * w12
      enddo
    enddo

  enddo

end subroutine deposit_v3_2d_s1
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine deposit_v3_3d_s1( f, ix, x, vq, np )

  implicit none

  integer, parameter :: rank = 3

  type( t_vdf ), intent(inout) :: f
  integer, dimension(:,:), intent(in) :: ix
  real(p_k_part), dimension(:,:), intent(in) :: x
  real(p_k_part), dimension(:,:), intent(in) :: vq
  integer, intent(in) :: np

  ! local variables
  integer :: l, i1, i2, i3, k1, k2, k3
  real( p_k_fld ) :: dx1, dx2, dx3, lq1, lq2, lq3, w123
  real( p_k_fld ), dimension(0:1) :: w1, w2, w3

  ! quantity is defined in the lower corner of the cell

  do l = 1, np

    i1 = ix(1,l)
    i2 = ix(2,l)
    i3 = ix(3,l)
    dx1 = real( x(1,l), p_k_fld )
    dx2 = real( x(2,l), p_k_fld )
    dx3 = real( x(3,l), p_k_fld )
    lq1 = real( vq(1,l) , p_k_fld )
    lq2 = real( vq(2,l) , p_k_fld )
    lq3 = real( vq(3,l) , p_k_fld )

    ! get spline weitghts for x, y and z
    call spline_s1( dx1, w1 )
    call spline_s1( dx2, w2 )
    call spline_s1( dx3, w3 )

    ! Deposit vector3 quantity
    do k3 = 0, 1
      do k2 = 0, 1
        do k1 = 0, 1
          w123 = w1(k1)* w2(k2) * w3(k3)
          f%f3(1,i1+k1,i2+k2,i3+k3) = f%f3(1,i1+k1,i2+k2,i3+k3) + lq1 * w123
          f%f3(2,i1+k1,i2+k2,i3+k3) = f%f3(2,i1+k1,i2+k2,i3+k3) + lq2 * w123
          f%f3(3,i1+k1,i2+k2,i3+k3) = f%f3(3,i1+k1,i2+k2,i3+k3) + lq3 * w123
        enddo
      enddo
    enddo

  enddo

end subroutine deposit_v3_3d_s1
!-----------------------------------------------------------------------------------------

! Clear template definitions
# 225 "spec/os-spec-v3dep.f90" 2

!********************************** Quadratic interpolation ********************************







! Lower point

! Upper point


# 1 "spec/os-spec-v3dep.f90" 1
!*****************************************************************************************
! m_species_v3dep module
!
! Handles deposition of a vector 3 quantity at the cell corner
!*****************************************************************************************

! if is defined then read the template definition at the end of the file
# 279 "spec/os-spec-v3dep.f90"
!*****************************************************************************************
!
! Template function definitions for v3 deposition
!
!*****************************************************************************************


!*****************************************************************************************
!* Linear Interpolation *
!*****************************************************************************************

! ----------------------------------------------------------------------------------------
subroutine deposit_v3_1d_s2( f, ix, x, vq, np )

  implicit none

  integer, parameter :: rank = 1

  type( t_vdf ), intent(inout) :: f
  integer, dimension(:,:), intent(in) :: ix
  real(p_k_part), dimension(:,:), intent(in) :: x
  real(p_k_part), dimension(:,:), intent(in) :: vq
  integer, intent(in) :: np

  ! local variables
  integer :: l, i1, k1
  real( p_k_fld ) :: lq1, lq2, lq3, dx1
  real( p_k_fld ), dimension(-1:1) :: w1

  ! quantity is defined in the lower corner of the cell

  do l = 1, np

     i1 = ix(1,l)
     dx1 = real( x(1,l), p_k_fld )
     lq1 = real( vq(1,l) , p_k_fld )
     lq2 = real( vq(2,l) , p_k_fld )
     lq3 = real( vq(3,l) , p_k_fld )

     ! get spline weitghts for x
     call spline_s2( dx1, w1 )

     ! Deposit Charge
     do k1 = -1, 1
       f%f1(1,i1+k1) = f%f1(1,i1+k1) + lq1 * w1(k1)
       f%f1(2,i1+k1) = f%f1(2,i1+k1) + lq2 * w1(k1)
       f%f1(3,i1+k1) = f%f1(3,i1+k1) + lq3 * w1(k1)
     enddo

  enddo

end subroutine deposit_v3_1d_s2
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine deposit_v3_2d_s2( f, ix, x, vq, np )

  implicit none

  integer, parameter :: rank = 2

  type( t_vdf ), intent(inout) :: f
  integer, dimension(:,:), intent(in) :: ix
  real(p_k_part), dimension(:,:), intent(in) :: x
  real(p_k_part), dimension(:,:), intent(in) :: vq
  integer, intent(in) :: np

  ! local variables
  integer :: l, i1, i2, k1, k2
  real(p_k_fld) :: dx1, dx2, lq1, lq2, lq3, w12
  real( p_k_fld ), dimension(-1:1) :: w1, w2

  ! quantity is defined in the lower corner of the cell

  do l = 1, np

    i1 = ix(1,l)
    i2 = ix(2,l)
    dx1 = real( x(1,l), p_k_fld )
    dx2 = real( x(2,l), p_k_fld )
    lq1 = real( vq(1,l) , p_k_fld )
    lq2 = real( vq(2,l) , p_k_fld )
    lq3 = real( vq(3,l) , p_k_fld )

    ! get spline weitghts for x and y
    call spline_s2( dx1, w1 )
    call spline_s2( dx2, w2 )

    ! Deposit vector3 quantity
    do k2 = -1, 1
      do k1 = -1, 1
        w12 = w1(k1)* w2(k2)
        f%f2(1,i1+k1,i2+k2) = f%f2(1,i1+k1,i2+k2) + lq1 * w12
        f%f2(2,i1+k1,i2+k2) = f%f2(2,i1+k1,i2+k2) + lq2 * w12
        f%f2(3,i1+k1,i2+k2) = f%f2(3,i1+k1,i2+k2) + lq3 * w12
      enddo
    enddo

  enddo

end subroutine deposit_v3_2d_s2
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine deposit_v3_3d_s2( f, ix, x, vq, np )

  implicit none

  integer, parameter :: rank = 3

  type( t_vdf ), intent(inout) :: f
  integer, dimension(:,:), intent(in) :: ix
  real(p_k_part), dimension(:,:), intent(in) :: x
  real(p_k_part), dimension(:,:), intent(in) :: vq
  integer, intent(in) :: np

  ! local variables
  integer :: l, i1, i2, i3, k1, k2, k3
  real( p_k_fld ) :: dx1, dx2, dx3, lq1, lq2, lq3, w123
  real( p_k_fld ), dimension(-1:1) :: w1, w2, w3

  ! quantity is defined in the lower corner of the cell

  do l = 1, np

    i1 = ix(1,l)
    i2 = ix(2,l)
    i3 = ix(3,l)
    dx1 = real( x(1,l), p_k_fld )
    dx2 = real( x(2,l), p_k_fld )
    dx3 = real( x(3,l), p_k_fld )
    lq1 = real( vq(1,l) , p_k_fld )
    lq2 = real( vq(2,l) , p_k_fld )
    lq3 = real( vq(3,l) , p_k_fld )

    ! get spline weitghts for x, y and z
    call spline_s2( dx1, w1 )
    call spline_s2( dx2, w2 )
    call spline_s2( dx3, w3 )

    ! Deposit vector3 quantity
    do k3 = -1, 1
      do k2 = -1, 1
        do k1 = -1, 1
          w123 = w1(k1)* w2(k2) * w3(k3)
          f%f3(1,i1+k1,i2+k2,i3+k3) = f%f3(1,i1+k1,i2+k2,i3+k3) + lq1 * w123
          f%f3(2,i1+k1,i2+k2,i3+k3) = f%f3(2,i1+k1,i2+k2,i3+k3) + lq2 * w123
          f%f3(3,i1+k1,i2+k2,i3+k3) = f%f3(3,i1+k1,i2+k2,i3+k3) + lq3 * w123
        enddo
      enddo
    enddo

  enddo

end subroutine deposit_v3_3d_s2
!-----------------------------------------------------------------------------------------

! Clear template definitions
# 240 "spec/os-spec-v3dep.f90" 2

!********************************** Cubic interpolation ********************************







! Lower point

! Upper point


# 1 "spec/os-spec-v3dep.f90" 1
!*****************************************************************************************
! m_species_v3dep module
!
! Handles deposition of a vector 3 quantity at the cell corner
!*****************************************************************************************

! if is defined then read the template definition at the end of the file
# 279 "spec/os-spec-v3dep.f90"
!*****************************************************************************************
!
! Template function definitions for v3 deposition
!
!*****************************************************************************************


!*****************************************************************************************
!* Linear Interpolation *
!*****************************************************************************************

! ----------------------------------------------------------------------------------------
subroutine deposit_v3_1d_s3( f, ix, x, vq, np )

  implicit none

  integer, parameter :: rank = 1

  type( t_vdf ), intent(inout) :: f
  integer, dimension(:,:), intent(in) :: ix
  real(p_k_part), dimension(:,:), intent(in) :: x
  real(p_k_part), dimension(:,:), intent(in) :: vq
  integer, intent(in) :: np

  ! local variables
  integer :: l, i1, k1
  real( p_k_fld ) :: lq1, lq2, lq3, dx1
  real( p_k_fld ), dimension(-1:2) :: w1

  ! quantity is defined in the lower corner of the cell

  do l = 1, np

     i1 = ix(1,l)
     dx1 = real( x(1,l), p_k_fld )
     lq1 = real( vq(1,l) , p_k_fld )
     lq2 = real( vq(2,l) , p_k_fld )
     lq3 = real( vq(3,l) , p_k_fld )

     ! get spline weitghts for x
     call spline_s3( dx1, w1 )

     ! Deposit Charge
     do k1 = -1, 2
       f%f1(1,i1+k1) = f%f1(1,i1+k1) + lq1 * w1(k1)
       f%f1(2,i1+k1) = f%f1(2,i1+k1) + lq2 * w1(k1)
       f%f1(3,i1+k1) = f%f1(3,i1+k1) + lq3 * w1(k1)
     enddo

  enddo

end subroutine deposit_v3_1d_s3
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine deposit_v3_2d_s3( f, ix, x, vq, np )

  implicit none

  integer, parameter :: rank = 2

  type( t_vdf ), intent(inout) :: f
  integer, dimension(:,:), intent(in) :: ix
  real(p_k_part), dimension(:,:), intent(in) :: x
  real(p_k_part), dimension(:,:), intent(in) :: vq
  integer, intent(in) :: np

  ! local variables
  integer :: l, i1, i2, k1, k2
  real(p_k_fld) :: dx1, dx2, lq1, lq2, lq3, w12
  real( p_k_fld ), dimension(-1:2) :: w1, w2

  ! quantity is defined in the lower corner of the cell

  do l = 1, np

    i1 = ix(1,l)
    i2 = ix(2,l)
    dx1 = real( x(1,l), p_k_fld )
    dx2 = real( x(2,l), p_k_fld )
    lq1 = real( vq(1,l) , p_k_fld )
    lq2 = real( vq(2,l) , p_k_fld )
    lq3 = real( vq(3,l) , p_k_fld )

    ! get spline weitghts for x and y
    call spline_s3( dx1, w1 )
    call spline_s3( dx2, w2 )

    ! Deposit vector3 quantity
    do k2 = -1, 2
      do k1 = -1, 2
        w12 = w1(k1)* w2(k2)
        f%f2(1,i1+k1,i2+k2) = f%f2(1,i1+k1,i2+k2) + lq1 * w12
        f%f2(2,i1+k1,i2+k2) = f%f2(2,i1+k1,i2+k2) + lq2 * w12
        f%f2(3,i1+k1,i2+k2) = f%f2(3,i1+k1,i2+k2) + lq3 * w12
      enddo
    enddo

  enddo

end subroutine deposit_v3_2d_s3
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine deposit_v3_3d_s3( f, ix, x, vq, np )

  implicit none

  integer, parameter :: rank = 3

  type( t_vdf ), intent(inout) :: f
  integer, dimension(:,:), intent(in) :: ix
  real(p_k_part), dimension(:,:), intent(in) :: x
  real(p_k_part), dimension(:,:), intent(in) :: vq
  integer, intent(in) :: np

  ! local variables
  integer :: l, i1, i2, i3, k1, k2, k3
  real( p_k_fld ) :: dx1, dx2, dx3, lq1, lq2, lq3, w123
  real( p_k_fld ), dimension(-1:2) :: w1, w2, w3

  ! quantity is defined in the lower corner of the cell

  do l = 1, np

    i1 = ix(1,l)
    i2 = ix(2,l)
    i3 = ix(3,l)
    dx1 = real( x(1,l), p_k_fld )
    dx2 = real( x(2,l), p_k_fld )
    dx3 = real( x(3,l), p_k_fld )
    lq1 = real( vq(1,l) , p_k_fld )
    lq2 = real( vq(2,l) , p_k_fld )
    lq3 = real( vq(3,l) , p_k_fld )

    ! get spline weitghts for x, y and z
    call spline_s3( dx1, w1 )
    call spline_s3( dx2, w2 )
    call spline_s3( dx3, w3 )

    ! Deposit vector3 quantity
    do k3 = -1, 2
      do k2 = -1, 2
        do k1 = -1, 2
          w123 = w1(k1)* w2(k2) * w3(k3)
          f%f3(1,i1+k1,i2+k2,i3+k3) = f%f3(1,i1+k1,i2+k2,i3+k3) + lq1 * w123
          f%f3(2,i1+k1,i2+k2,i3+k3) = f%f3(2,i1+k1,i2+k2,i3+k3) + lq2 * w123
          f%f3(3,i1+k1,i2+k2,i3+k3) = f%f3(3,i1+k1,i2+k2,i3+k3) + lq3 * w123
        enddo
      enddo
    enddo

  enddo

end subroutine deposit_v3_3d_s3
!-----------------------------------------------------------------------------------------

! Clear template definitions
# 255 "spec/os-spec-v3dep.f90" 2

!********************************** Quartic interpolation ********************************







! Lower point

! Upper point


# 1 "spec/os-spec-v3dep.f90" 1
!*****************************************************************************************
! m_species_v3dep module
!
! Handles deposition of a vector 3 quantity at the cell corner
!*****************************************************************************************

! if is defined then read the template definition at the end of the file
# 279 "spec/os-spec-v3dep.f90"
!*****************************************************************************************
!
! Template function definitions for v3 deposition
!
!*****************************************************************************************


!*****************************************************************************************
!* Linear Interpolation *
!*****************************************************************************************

! ----------------------------------------------------------------------------------------
subroutine deposit_v3_1d_s4( f, ix, x, vq, np )

  implicit none

  integer, parameter :: rank = 1

  type( t_vdf ), intent(inout) :: f
  integer, dimension(:,:), intent(in) :: ix
  real(p_k_part), dimension(:,:), intent(in) :: x
  real(p_k_part), dimension(:,:), intent(in) :: vq
  integer, intent(in) :: np

  ! local variables
  integer :: l, i1, k1
  real( p_k_fld ) :: lq1, lq2, lq3, dx1
  real( p_k_fld ), dimension(-2:2) :: w1

  ! quantity is defined in the lower corner of the cell

  do l = 1, np

     i1 = ix(1,l)
     dx1 = real( x(1,l), p_k_fld )
     lq1 = real( vq(1,l) , p_k_fld )
     lq2 = real( vq(2,l) , p_k_fld )
     lq3 = real( vq(3,l) , p_k_fld )

     ! get spline weitghts for x
     call spline_s4( dx1, w1 )

     ! Deposit Charge
     do k1 = -2, 2
       f%f1(1,i1+k1) = f%f1(1,i1+k1) + lq1 * w1(k1)
       f%f1(2,i1+k1) = f%f1(2,i1+k1) + lq2 * w1(k1)
       f%f1(3,i1+k1) = f%f1(3,i1+k1) + lq3 * w1(k1)
     enddo

  enddo

end subroutine deposit_v3_1d_s4
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine deposit_v3_2d_s4( f, ix, x, vq, np )

  implicit none

  integer, parameter :: rank = 2

  type( t_vdf ), intent(inout) :: f
  integer, dimension(:,:), intent(in) :: ix
  real(p_k_part), dimension(:,:), intent(in) :: x
  real(p_k_part), dimension(:,:), intent(in) :: vq
  integer, intent(in) :: np

  ! local variables
  integer :: l, i1, i2, k1, k2
  real(p_k_fld) :: dx1, dx2, lq1, lq2, lq3, w12
  real( p_k_fld ), dimension(-2:2) :: w1, w2

  ! quantity is defined in the lower corner of the cell

  do l = 1, np

    i1 = ix(1,l)
    i2 = ix(2,l)
    dx1 = real( x(1,l), p_k_fld )
    dx2 = real( x(2,l), p_k_fld )
    lq1 = real( vq(1,l) , p_k_fld )
    lq2 = real( vq(2,l) , p_k_fld )
    lq3 = real( vq(3,l) , p_k_fld )

    ! get spline weitghts for x and y
    call spline_s4( dx1, w1 )
    call spline_s4( dx2, w2 )

    ! Deposit vector3 quantity
    do k2 = -2, 2
      do k1 = -2, 2
        w12 = w1(k1)* w2(k2)
        f%f2(1,i1+k1,i2+k2) = f%f2(1,i1+k1,i2+k2) + lq1 * w12
        f%f2(2,i1+k1,i2+k2) = f%f2(2,i1+k1,i2+k2) + lq2 * w12
        f%f2(3,i1+k1,i2+k2) = f%f2(3,i1+k1,i2+k2) + lq3 * w12
      enddo
    enddo

  enddo

end subroutine deposit_v3_2d_s4
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine deposit_v3_3d_s4( f, ix, x, vq, np )

  implicit none

  integer, parameter :: rank = 3

  type( t_vdf ), intent(inout) :: f
  integer, dimension(:,:), intent(in) :: ix
  real(p_k_part), dimension(:,:), intent(in) :: x
  real(p_k_part), dimension(:,:), intent(in) :: vq
  integer, intent(in) :: np

  ! local variables
  integer :: l, i1, i2, i3, k1, k2, k3
  real( p_k_fld ) :: dx1, dx2, dx3, lq1, lq2, lq3, w123
  real( p_k_fld ), dimension(-2:2) :: w1, w2, w3

  ! quantity is defined in the lower corner of the cell

  do l = 1, np

    i1 = ix(1,l)
    i2 = ix(2,l)
    i3 = ix(3,l)
    dx1 = real( x(1,l), p_k_fld )
    dx2 = real( x(2,l), p_k_fld )
    dx3 = real( x(3,l), p_k_fld )
    lq1 = real( vq(1,l) , p_k_fld )
    lq2 = real( vq(2,l) , p_k_fld )
    lq3 = real( vq(3,l) , p_k_fld )

    ! get spline weitghts for x, y and z
    call spline_s4( dx1, w1 )
    call spline_s4( dx2, w2 )
    call spline_s4( dx3, w3 )

    ! Deposit vector3 quantity
    do k3 = -2, 2
      do k2 = -2, 2
        do k1 = -2, 2
          w123 = w1(k1)* w2(k2) * w3(k3)
          f%f3(1,i1+k1,i2+k2,i3+k3) = f%f3(1,i1+k1,i2+k2,i3+k3) + lq1 * w123
          f%f3(2,i1+k1,i2+k2,i3+k3) = f%f3(2,i1+k1,i2+k2,i3+k3) + lq2 * w123
          f%f3(3,i1+k1,i2+k2,i3+k3) = f%f3(3,i1+k1,i2+k2,i3+k3) + lq3 * w123
        enddo
      enddo
    enddo

  enddo

end subroutine deposit_v3_3d_s4
!-----------------------------------------------------------------------------------------

! Clear template definitions
# 270 "spec/os-spec-v3dep.f90" 2

!*****************************************************************************************



end module m_species_v3dep
