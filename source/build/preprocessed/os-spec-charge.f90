# 1 "spec/os-spec-charge.f90"
# 1 "<built-in>" 1
# 1 "<built-in>" 3
# 467 "<built-in>" 3
# 1 "<command line>" 1
# 1 "<built-in>" 2
# 1 "spec/os-spec-charge.f90" 2
!*****************************************************************************************
! m_species_charge module
!
! Charge deposition routines
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
# 12 "spec/os-spec-charge.f90" 2
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
# 13 "spec/os-spec-charge.f90" 2

module m_species_charge

use m_species_define

use m_system
use m_parameters
use m_vdf_define

private

interface deposit_rho
  module procedure deposit_rho_spec
end interface

interface norm_charge_cyl
  module procedure norm_charge_cyl
end interface

public :: deposit_rho, norm_charge_cyl

public :: deposit_rho_1d_s1, deposit_rho_1d_s2, deposit_rho_1d_s3, deposit_rho_1d_s4
public :: deposit_rho_2d_s1, deposit_rho_2d_s2, deposit_rho_2d_s3, deposit_rho_2d_s4
public :: deposit_rho_3d_s1, deposit_rho_3d_s2, deposit_rho_3d_s3, deposit_rho_3d_s4

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


!-------------------------------------------------------------------------------
! Deposit species charge using correct interpolation scheme.
!-------------------------------------------------------------------------------
subroutine deposit_rho_spec( this, rho )

  implicit none

  ! dummy variables
  class( t_species ), intent(in) :: this
  type( t_vdf ), intent(inout) :: rho

  ! executable statements

  select case ( this%interpolation )

  case (p_linear)
    select case ( p_x_dim )
    case (1)
      call deposit_rho_1d_s1( rho, this%ix, this%x, this%q, this%num_par )
    case (2)
      call deposit_rho_2d_s1( rho, this%ix, this%x, this%q, this%num_par )
    case (3)
      call deposit_rho_3d_s1( rho, this%ix, this%x, this%q, this%num_par )
    end select

  case (p_quadratic)
    select case ( p_x_dim )
    case (1)
      call deposit_rho_1d_s2( rho, this%ix, this%x, this%q, this%num_par )
    case (2)
      call deposit_rho_2d_s2( rho, this%ix, this%x, this%q, this%num_par )
    case (3)
      call deposit_rho_3d_s2( rho, this%ix, this%x, this%q, this%num_par )
    end select

  case (p_cubic)
    select case ( p_x_dim )
    case (1)
      call deposit_rho_1d_s3( rho, this%ix, this%x, this%q, this%num_par )
    case (2)
      call deposit_rho_2d_s3( rho, this%ix, this%x, this%q, this%num_par )
    case (3)
      call deposit_rho_3d_s3( rho, this%ix, this%x, this%q, this%num_par )
    end select

  case (p_quartic)
    select case ( p_x_dim )
    case (1)
      call deposit_rho_1d_s4( rho, this%ix, this%x, this%q, this%num_par )
    case (2)
      call deposit_rho_2d_s4( rho, this%ix, this%x, this%q, this%num_par )
    case (3)
      call deposit_rho_3d_s4( rho, this%ix, this%x, this%q, this%num_par )
    end select

  case default
    write(err_buf__,*) 'Not implemented';call err__("spec/os-spec-charge.f90",192)
    call abort_program( p_err_notimplemented )

  end select

end subroutine deposit_rho_spec
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
subroutine deposit_density_spec( this, rho, i1, i2, q )

  implicit none

  ! dummy variables
  class( t_species ), intent(in) :: this
  type( t_vdf ), intent(inout) :: rho
  integer, intent(in) :: i1, i2

  real(p_k_part), dimension(:), intent(in) :: q

  integer :: np

  ! number of particles to deposit
  np = i2 - i1 + 1

  ! deposit given density
  select case ( this%interpolation )

  case (p_linear)
    select case ( p_x_dim )
    case (1)
      call deposit_rho_1d_s1( rho, this%ix(:,i1:i2), this%x(:,i1:i2), q, np )
    case (2)
      call deposit_rho_2d_s1( rho, this%ix(:,i1:i2), this%x(:,i1:i2), q, np )
    case (3)
      call deposit_rho_3d_s1( rho, this%ix(:,i1:i2), this%x(:,i1:i2), q, np )
    end select

  case (p_quadratic)
    select case ( p_x_dim )
    case (1)
      call deposit_rho_1d_s2( rho, this%ix(:,i1:i2), this%x(:,i1:i2), q, np )
    case (2)
      call deposit_rho_2d_s2( rho, this%ix(:,i1:i2), this%x(:,i1:i2), q, np )
    case (3)
      call deposit_rho_3d_s2( rho, this%ix(:,i1:i2), this%x(:,i1:i2), q, np )
    end select

  case (p_cubic)
    select case ( p_x_dim )
    case (1)
      call deposit_rho_1d_s3( rho, this%ix(:,i1:i2), this%x(:,i1:i2), q, np )
    case (2)
      call deposit_rho_2d_s3( rho, this%ix(:,i1:i2), this%x(:,i1:i2), q, np )
    case (3)
      call deposit_rho_3d_s3( rho, this%ix(:,i1:i2), this%x(:,i1:i2), q, np )
    end select

  case (p_quartic)
    select case ( p_x_dim )
    case (1)
      call deposit_rho_1d_s4( rho, this%ix(:,i1:i2), this%x(:,i1:i2), q, np )
    case (2)
      call deposit_rho_2d_s4( rho, this%ix(:,i1:i2), this%x(:,i1:i2), q, np )
    case (3)
      call deposit_rho_3d_s4( rho, this%ix(:,i1:i2), this%x(:,i1:i2), q, np )
    end select

  case default
    write(err_buf__,*) 'Not implemented';call err__("spec/os-spec-charge.f90",263)
    call abort_program( p_err_notimplemented )

  end select

end subroutine deposit_density_spec
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Normalize charge in cylindrical coordinates.
!-------------------------------------------------------------------------------
subroutine norm_charge_cyl( rho, gir_pos, dr )

  implicit none

  type( t_vdf ), intent(inout) :: rho
  integer, intent(in) :: gir_pos ! local grid position on global grid
  real( p_k_fld ), intent(in) :: dr

  integer :: ir
  real( p_k_fld ) :: r

  ! normalize for 'ring' particles
  ! Note that the lower radial spatial boundary of a cylindrical geometry simulation is always
  ! -dr/2, where dr is the radial cell size. Also note that charge beyond axial boundary is reversed
  ! since r is considered to be < 0.

  do ir = lbound( rho%f2, 3), ubound( rho%f2, 3)
     r = ( (ir+ gir_pos - 2) - 0.5_p_k_fld )* dr
     rho%f2(1,:,ir) = rho%f2(1,:,ir) / r
  enddo

  ! Fold axial guard cells back into simulation space
  if ( gir_pos == 1 ) then
    do ir = 0, (1 - lbound( rho%f2, 3))
      rho%f2(1,:,ir+2) = rho%f2(1,:,ir+2) + rho%f2(1,:,1-ir)
      rho%f2(1,:,1-ir) = rho%f2(1,:,ir+2)
    enddo
  endif


end subroutine norm_charge_cyl
!-------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Generate specific template functions for linear, quadratic, cubic and quartic
! interpolation levels.
!-----------------------------------------------------------------------------------------



!********************************** Linear interpolation ********************************







! Lower point

! Upper point


# 1 "spec/os-spec-charge.f90" 1
!*****************************************************************************************
! m_species_charge module
!
! Charge deposition routines
!*****************************************************************************************

! if is defined then read the template definition at the end of the file
# 457 "spec/os-spec-charge.f90"
!*****************************************************************************************
!
! Template function definitions for charge deposition
!
!*****************************************************************************************

subroutine deposit_rho_1d_s1( rho, ix, x, q, np )

  implicit none

  integer, parameter :: rank = 1

  type( t_vdf ), intent(inout) :: rho
  integer, dimension(:,:), intent(in) :: ix
  real(p_k_part), dimension(:,:), intent(in) :: x
  real(p_k_part), dimension( : ), intent(in) :: q
  integer, intent(in) :: np

  ! local variables
  integer :: l, i1, k1
  real( p_k_fld ) :: lq, dx1
  real( p_k_fld ), dimension(0:1) :: w1

! rho is defined in the lower corner of the cell

  do l = 1, np

     i1 = ix(1,l)
     dx1 = real( x(1,l), p_k_fld )
     lq = real( q(l) , p_k_fld )

     ! get spline weitghts for x
     call spline_s1( dx1, w1 )

     ! Deposit Charge
     do k1 = 0, 1
       rho%f1(1,i1 + k1 ) = rho%f1(1, i1 + k1 ) + lq * w1(k1)
     enddo
  enddo

end subroutine deposit_rho_1d_s1


subroutine deposit_rho_2d_s1( rho, ix, x, q, np )

  implicit none

  integer, parameter :: rank = 2

  type( t_vdf ), intent(inout) :: rho
  integer, dimension(:,:), intent(in) :: ix
  real(p_k_part), dimension(:,:), intent(in) :: x
  real(p_k_part), dimension( : ), intent(in) :: q
  integer, intent(in) :: np

  ! local variables
  integer :: l, i1, i2, k1, k2
  real(p_k_fld) :: dx1, dx2, lq
  real( p_k_fld ), dimension(0:1) :: w1, w2

! rho is defined in the lower corner of the cell

  do l = 1, np

    i1 = ix(1,l)
    i2 = ix(2,l)
    dx1 = real( x(1,l), p_k_fld )
    dx2 = real( x(2,l), p_k_fld )
    lq = real( q(l) , p_k_fld )

    ! get spline weitghts for x and y
    call spline_s1( dx1, w1 )
    call spline_s1( dx2, w2 )

    ! Deposit Charge
    do k2 = 0, 1
      do k1 = 0, 1
        rho%f2(1,i1 + k1,i2 + k2) = rho%f2(1,i1 + k1,i2 + k2) + lq * w1(k1)* w2(k2)
      enddo
    enddo

  enddo

end subroutine deposit_rho_2d_s1
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine deposit_rho_3d_s1( rho, ix, x, q, np )

  implicit none

  integer, parameter :: rank = 3

  type( t_vdf ), intent(inout) :: rho
  integer, dimension(:,:), intent(in) :: ix
  real(p_k_part), dimension(:,:), intent(in) :: x
  real(p_k_part), dimension( : ), intent(in) :: q
  integer, intent(in) :: np

  ! local variables
  integer :: l, i1, i2, i3, k1, k2, k3
  real( p_k_fld ) :: dx1, dx2, dx3, lq
  real( p_k_fld ), dimension(0:1) :: w1, w2, w3

  ! rho is defined in the lower corner of the cell

  do l = 1, np

    i1 = ix(1,l)
    i2 = ix(2,l)
    i3 = ix(3,l)
    dx1 = real( x(1,l), p_k_fld )
    dx2 = real( x(2,l), p_k_fld )
    dx3 = real( x(3,l), p_k_fld )
    lq = real( q(l) , p_k_fld )

    ! get spline weitghts for x, y and z
    call spline_s1( dx1, w1 )
    call spline_s1( dx2, w2 )
    call spline_s1( dx3, w3 )

    ! Deposit Charge
    do k3 = 0, 1
      do k2 = 0, 1
        do k1 = 0, 1
          rho%f3(1,i1+k1,i2+k2,i3+k3) = rho%f3(1,i1+k1,i2+k2,i3+k3) + &
                                        lq * w1(k1)* w2(k2)* w3(k3)
        enddo
      enddo
    enddo
  enddo

end subroutine deposit_rho_3d_s1
!-------------------------------------------------------------------------------

! Clear template definitions
# 328 "spec/os-spec-charge.f90" 2

!********************************** Quadratic interpolation ********************************







! Lower point

! Upper point


# 1 "spec/os-spec-charge.f90" 1
!*****************************************************************************************
! m_species_charge module
!
! Charge deposition routines
!*****************************************************************************************

! if is defined then read the template definition at the end of the file
# 457 "spec/os-spec-charge.f90"
!*****************************************************************************************
!
! Template function definitions for charge deposition
!
!*****************************************************************************************

subroutine deposit_rho_1d_s2( rho, ix, x, q, np )

  implicit none

  integer, parameter :: rank = 1

  type( t_vdf ), intent(inout) :: rho
  integer, dimension(:,:), intent(in) :: ix
  real(p_k_part), dimension(:,:), intent(in) :: x
  real(p_k_part), dimension( : ), intent(in) :: q
  integer, intent(in) :: np

  ! local variables
  integer :: l, i1, k1
  real( p_k_fld ) :: lq, dx1
  real( p_k_fld ), dimension(-1:1) :: w1

! rho is defined in the lower corner of the cell

  do l = 1, np

     i1 = ix(1,l)
     dx1 = real( x(1,l), p_k_fld )
     lq = real( q(l) , p_k_fld )

     ! get spline weitghts for x
     call spline_s2( dx1, w1 )

     ! Deposit Charge
     do k1 = -1, 1
       rho%f1(1,i1 + k1 ) = rho%f1(1, i1 + k1 ) + lq * w1(k1)
     enddo
  enddo

end subroutine deposit_rho_1d_s2


subroutine deposit_rho_2d_s2( rho, ix, x, q, np )

  implicit none

  integer, parameter :: rank = 2

  type( t_vdf ), intent(inout) :: rho
  integer, dimension(:,:), intent(in) :: ix
  real(p_k_part), dimension(:,:), intent(in) :: x
  real(p_k_part), dimension( : ), intent(in) :: q
  integer, intent(in) :: np

  ! local variables
  integer :: l, i1, i2, k1, k2
  real(p_k_fld) :: dx1, dx2, lq
  real( p_k_fld ), dimension(-1:1) :: w1, w2

! rho is defined in the lower corner of the cell

  do l = 1, np

    i1 = ix(1,l)
    i2 = ix(2,l)
    dx1 = real( x(1,l), p_k_fld )
    dx2 = real( x(2,l), p_k_fld )
    lq = real( q(l) , p_k_fld )

    ! get spline weitghts for x and y
    call spline_s2( dx1, w1 )
    call spline_s2( dx2, w2 )

    ! Deposit Charge
    do k2 = -1, 1
      do k1 = -1, 1
        rho%f2(1,i1 + k1,i2 + k2) = rho%f2(1,i1 + k1,i2 + k2) + lq * w1(k1)* w2(k2)
      enddo
    enddo

  enddo

end subroutine deposit_rho_2d_s2
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine deposit_rho_3d_s2( rho, ix, x, q, np )

  implicit none

  integer, parameter :: rank = 3

  type( t_vdf ), intent(inout) :: rho
  integer, dimension(:,:), intent(in) :: ix
  real(p_k_part), dimension(:,:), intent(in) :: x
  real(p_k_part), dimension( : ), intent(in) :: q
  integer, intent(in) :: np

  ! local variables
  integer :: l, i1, i2, i3, k1, k2, k3
  real( p_k_fld ) :: dx1, dx2, dx3, lq
  real( p_k_fld ), dimension(-1:1) :: w1, w2, w3

  ! rho is defined in the lower corner of the cell

  do l = 1, np

    i1 = ix(1,l)
    i2 = ix(2,l)
    i3 = ix(3,l)
    dx1 = real( x(1,l), p_k_fld )
    dx2 = real( x(2,l), p_k_fld )
    dx3 = real( x(3,l), p_k_fld )
    lq = real( q(l) , p_k_fld )

    ! get spline weitghts for x, y and z
    call spline_s2( dx1, w1 )
    call spline_s2( dx2, w2 )
    call spline_s2( dx3, w3 )

    ! Deposit Charge
    do k3 = -1, 1
      do k2 = -1, 1
        do k1 = -1, 1
          rho%f3(1,i1+k1,i2+k2,i3+k3) = rho%f3(1,i1+k1,i2+k2,i3+k3) + &
                                        lq * w1(k1)* w2(k2)* w3(k3)
        enddo
      enddo
    enddo
  enddo

end subroutine deposit_rho_3d_s2
!-------------------------------------------------------------------------------

! Clear template definitions
# 343 "spec/os-spec-charge.f90" 2

!********************************** Cubic interpolation ********************************







! Lower point

! Upper point


# 1 "spec/os-spec-charge.f90" 1
!*****************************************************************************************
! m_species_charge module
!
! Charge deposition routines
!*****************************************************************************************

! if is defined then read the template definition at the end of the file
# 457 "spec/os-spec-charge.f90"
!*****************************************************************************************
!
! Template function definitions for charge deposition
!
!*****************************************************************************************

subroutine deposit_rho_1d_s3( rho, ix, x, q, np )

  implicit none

  integer, parameter :: rank = 1

  type( t_vdf ), intent(inout) :: rho
  integer, dimension(:,:), intent(in) :: ix
  real(p_k_part), dimension(:,:), intent(in) :: x
  real(p_k_part), dimension( : ), intent(in) :: q
  integer, intent(in) :: np

  ! local variables
  integer :: l, i1, k1
  real( p_k_fld ) :: lq, dx1
  real( p_k_fld ), dimension(-1:2) :: w1

! rho is defined in the lower corner of the cell

  do l = 1, np

     i1 = ix(1,l)
     dx1 = real( x(1,l), p_k_fld )
     lq = real( q(l) , p_k_fld )

     ! get spline weitghts for x
     call spline_s3( dx1, w1 )

     ! Deposit Charge
     do k1 = -1, 2
       rho%f1(1,i1 + k1 ) = rho%f1(1, i1 + k1 ) + lq * w1(k1)
     enddo
  enddo

end subroutine deposit_rho_1d_s3


subroutine deposit_rho_2d_s3( rho, ix, x, q, np )

  implicit none

  integer, parameter :: rank = 2

  type( t_vdf ), intent(inout) :: rho
  integer, dimension(:,:), intent(in) :: ix
  real(p_k_part), dimension(:,:), intent(in) :: x
  real(p_k_part), dimension( : ), intent(in) :: q
  integer, intent(in) :: np

  ! local variables
  integer :: l, i1, i2, k1, k2
  real(p_k_fld) :: dx1, dx2, lq
  real( p_k_fld ), dimension(-1:2) :: w1, w2

! rho is defined in the lower corner of the cell

  do l = 1, np

    i1 = ix(1,l)
    i2 = ix(2,l)
    dx1 = real( x(1,l), p_k_fld )
    dx2 = real( x(2,l), p_k_fld )
    lq = real( q(l) , p_k_fld )

    ! get spline weitghts for x and y
    call spline_s3( dx1, w1 )
    call spline_s3( dx2, w2 )

    ! Deposit Charge
    do k2 = -1, 2
      do k1 = -1, 2
        rho%f2(1,i1 + k1,i2 + k2) = rho%f2(1,i1 + k1,i2 + k2) + lq * w1(k1)* w2(k2)
      enddo
    enddo

  enddo

end subroutine deposit_rho_2d_s3
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine deposit_rho_3d_s3( rho, ix, x, q, np )

  implicit none

  integer, parameter :: rank = 3

  type( t_vdf ), intent(inout) :: rho
  integer, dimension(:,:), intent(in) :: ix
  real(p_k_part), dimension(:,:), intent(in) :: x
  real(p_k_part), dimension( : ), intent(in) :: q
  integer, intent(in) :: np

  ! local variables
  integer :: l, i1, i2, i3, k1, k2, k3
  real( p_k_fld ) :: dx1, dx2, dx3, lq
  real( p_k_fld ), dimension(-1:2) :: w1, w2, w3

  ! rho is defined in the lower corner of the cell

  do l = 1, np

    i1 = ix(1,l)
    i2 = ix(2,l)
    i3 = ix(3,l)
    dx1 = real( x(1,l), p_k_fld )
    dx2 = real( x(2,l), p_k_fld )
    dx3 = real( x(3,l), p_k_fld )
    lq = real( q(l) , p_k_fld )

    ! get spline weitghts for x, y and z
    call spline_s3( dx1, w1 )
    call spline_s3( dx2, w2 )
    call spline_s3( dx3, w3 )

    ! Deposit Charge
    do k3 = -1, 2
      do k2 = -1, 2
        do k1 = -1, 2
          rho%f3(1,i1+k1,i2+k2,i3+k3) = rho%f3(1,i1+k1,i2+k2,i3+k3) + &
                                        lq * w1(k1)* w2(k2)* w3(k3)
        enddo
      enddo
    enddo
  enddo

end subroutine deposit_rho_3d_s3
!-------------------------------------------------------------------------------

! Clear template definitions
# 358 "spec/os-spec-charge.f90" 2

!********************************** Quartic interpolation ********************************







! Lower point

! Upper point


# 1 "spec/os-spec-charge.f90" 1
!*****************************************************************************************
! m_species_charge module
!
! Charge deposition routines
!*****************************************************************************************

! if is defined then read the template definition at the end of the file
# 457 "spec/os-spec-charge.f90"
!*****************************************************************************************
!
! Template function definitions for charge deposition
!
!*****************************************************************************************

subroutine deposit_rho_1d_s4( rho, ix, x, q, np )

  implicit none

  integer, parameter :: rank = 1

  type( t_vdf ), intent(inout) :: rho
  integer, dimension(:,:), intent(in) :: ix
  real(p_k_part), dimension(:,:), intent(in) :: x
  real(p_k_part), dimension( : ), intent(in) :: q
  integer, intent(in) :: np

  ! local variables
  integer :: l, i1, k1
  real( p_k_fld ) :: lq, dx1
  real( p_k_fld ), dimension(-2:2) :: w1

! rho is defined in the lower corner of the cell

  do l = 1, np

     i1 = ix(1,l)
     dx1 = real( x(1,l), p_k_fld )
     lq = real( q(l) , p_k_fld )

     ! get spline weitghts for x
     call spline_s4( dx1, w1 )

     ! Deposit Charge
     do k1 = -2, 2
       rho%f1(1,i1 + k1 ) = rho%f1(1, i1 + k1 ) + lq * w1(k1)
     enddo
  enddo

end subroutine deposit_rho_1d_s4


subroutine deposit_rho_2d_s4( rho, ix, x, q, np )

  implicit none

  integer, parameter :: rank = 2

  type( t_vdf ), intent(inout) :: rho
  integer, dimension(:,:), intent(in) :: ix
  real(p_k_part), dimension(:,:), intent(in) :: x
  real(p_k_part), dimension( : ), intent(in) :: q
  integer, intent(in) :: np

  ! local variables
  integer :: l, i1, i2, k1, k2
  real(p_k_fld) :: dx1, dx2, lq
  real( p_k_fld ), dimension(-2:2) :: w1, w2

! rho is defined in the lower corner of the cell

  do l = 1, np

    i1 = ix(1,l)
    i2 = ix(2,l)
    dx1 = real( x(1,l), p_k_fld )
    dx2 = real( x(2,l), p_k_fld )
    lq = real( q(l) , p_k_fld )

    ! get spline weitghts for x and y
    call spline_s4( dx1, w1 )
    call spline_s4( dx2, w2 )

    ! Deposit Charge
    do k2 = -2, 2
      do k1 = -2, 2
        rho%f2(1,i1 + k1,i2 + k2) = rho%f2(1,i1 + k1,i2 + k2) + lq * w1(k1)* w2(k2)
      enddo
    enddo

  enddo

end subroutine deposit_rho_2d_s4
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine deposit_rho_3d_s4( rho, ix, x, q, np )

  implicit none

  integer, parameter :: rank = 3

  type( t_vdf ), intent(inout) :: rho
  integer, dimension(:,:), intent(in) :: ix
  real(p_k_part), dimension(:,:), intent(in) :: x
  real(p_k_part), dimension( : ), intent(in) :: q
  integer, intent(in) :: np

  ! local variables
  integer :: l, i1, i2, i3, k1, k2, k3
  real( p_k_fld ) :: dx1, dx2, dx3, lq
  real( p_k_fld ), dimension(-2:2) :: w1, w2, w3

  ! rho is defined in the lower corner of the cell

  do l = 1, np

    i1 = ix(1,l)
    i2 = ix(2,l)
    i3 = ix(3,l)
    dx1 = real( x(1,l), p_k_fld )
    dx2 = real( x(2,l), p_k_fld )
    dx3 = real( x(3,l), p_k_fld )
    lq = real( q(l) , p_k_fld )

    ! get spline weitghts for x, y and z
    call spline_s4( dx1, w1 )
    call spline_s4( dx2, w2 )
    call spline_s4( dx3, w3 )

    ! Deposit Charge
    do k3 = -2, 2
      do k2 = -2, 2
        do k1 = -2, 2
          rho%f3(1,i1+k1,i2+k2,i3+k3) = rho%f3(1,i1+k1,i2+k2,i3+k3) + &
                                        lq * w1(k1)* w2(k2)* w3(k3)
        enddo
      enddo
    enddo
  enddo

end subroutine deposit_rho_3d_s4
!-------------------------------------------------------------------------------

! Clear template definitions
# 373 "spec/os-spec-charge.f90" 2

!*****************************************************************************************

end module m_species_charge

!-------------------------------------------------------------------------------
subroutine deposit_density_spec( this, rho, i1, i2, q )

  use m_system
  use m_parameters

  use m_species_define, only : t_species
  use m_vdf_define, only : t_vdf
  use m_species_charge

  implicit none

  ! dummy variables
  class( t_species ), intent(in) :: this
  type( t_vdf ), intent(inout) :: rho
  integer, intent(in) :: i1, i2

  real(p_k_part), dimension(:), intent(in) :: q

  integer :: np

  ! number of particles to deposit
  np = i2 - i1 + 1

  ! deposit given density
  select case ( this%interpolation )

     case (p_linear)
        select case ( p_x_dim )
         case (1)
           call deposit_rho_1d_s1( rho, this%ix(:,i1:i2), this%x(:,i1:i2), q, np )
         case (2)
           call deposit_rho_2d_s1( rho, this%ix(:,i1:i2), this%x(:,i1:i2), q, np )
         case (3)
           call deposit_rho_3d_s1( rho, this%ix(:,i1:i2), this%x(:,i1:i2), q, np )
        end select

     case (p_quadratic)
        select case ( p_x_dim )
         case (1)
           call deposit_rho_1d_s2( rho, this%ix(:,i1:i2), this%x(:,i1:i2), q, np )
         case (2)
           call deposit_rho_2d_s2( rho, this%ix(:,i1:i2), this%x(:,i1:i2), q, np )
         case (3)
           call deposit_rho_3d_s2( rho, this%ix(:,i1:i2), this%x(:,i1:i2), q, np )
        end select

     case (p_cubic)
        select case ( p_x_dim )
         case (1)
           call deposit_rho_1d_s3( rho, this%ix(:,i1:i2), this%x(:,i1:i2), q, np )
         case (2)
           call deposit_rho_2d_s3( rho, this%ix(:,i1:i2), this%x(:,i1:i2), q, np )
         case (3)
           call deposit_rho_3d_s3( rho, this%ix(:,i1:i2), this%x(:,i1:i2), q, np )
        end select

     case (p_quartic)
        select case ( p_x_dim )
         case (1)
           call deposit_rho_1d_s4( rho, this%ix(:,i1:i2), this%x(:,i1:i2), q, np )
         case (2)
           call deposit_rho_2d_s4( rho, this%ix(:,i1:i2), this%x(:,i1:i2), q, np )
         case (3)
           call deposit_rho_3d_s4( rho, this%ix(:,i1:i2), this%x(:,i1:i2), q, np )
        end select

     case default
        write(err_buf__,*) 'Not implemented';call err__("spec/os-spec-charge.f90",446)
        call abort_program( p_err_notimplemented )

  end select


end subroutine deposit_density_spec
!-------------------------------------------------------------------------------
