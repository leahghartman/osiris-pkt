# 1 "gr/os-species-charge-gr.f03"
# 1 "<built-in>" 1
# 1 "<built-in>" 3
# 467 "<built-in>" 3
# 1 "<command line>" 1
# 1 "<built-in>" 2
# 1 "gr/os-species-charge-gr.f03" 2
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
# 2 "gr/os-species-charge-gr.f03" 2
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
# 3 "gr/os-species-charge-gr.f03" 2

module m_species_charge_gr

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
# 7 "gr/os-species-charge-gr.f03" 2

use m_parameters
use m_geometry_gr

implicit none

private

interface deposit_rho_2d_srect_gr
  module procedure deposit_rho_2d_srect_gr
end interface

interface deposit_rho_2d_srect_curved_gr
  module procedure deposit_rho_2d_srect_curved_gr
end interface

interface deposit_rho_2d_sint_gr
  module procedure deposit_rho_2d_sint_gr
end interface

interface deposit_rho_2d_sint_curved_gr
  module procedure deposit_rho_2d_sint_curved_gr
end interface

interface deposit_rho_2d_shalf_gr
  module procedure deposit_rho_2d_shalf_gr
end interface

public deposit_rho_2d_srect_gr, deposit_rho_2d_srect_curved_gr
public deposit_rho_2d_sint_gr , deposit_rho_2d_sint_curved_gr
public deposit_rho_2d_shalf_gr

contains

!-------------------------------------------------------------------------------
! Deposit charge using the rectangular rule on the spline_r - Minkowski
!-------------------------------------------------------------------------------
subroutine deposit_rho_2d_srect_gr( rho, geometry, ix, x, q, np )

  use m_vdf_define, only : t_vdf

  implicit none

  integer, parameter :: rank = 2

  type( t_vdf ), intent(inout) :: rho
  type( t_geometry_gr ), intent(in), target :: geometry
  integer, dimension(:,:), intent(in) :: ix
  real(p_k_part), dimension(:,:), intent(in) :: x
  real(p_k_part), dimension( : ), intent(in) :: q
  integer, intent(in) :: np

  ! local variables
  integer :: l, i1, i2
  real(p_k_fld) :: x1, x2, lq
  real(p_k_fld) :: qnorm

  real(p_k_fld), dimension(-1:1) :: Sr
  real(p_k_fld), dimension(-1:1) :: St
  integer :: k1, k2

  ! flat and curved case - rectangular rule
  qnorm = real( 0.0795774715459477 * geometry%A, p_k_fld ) ! 1/4pi

  do l = 1, np

    i1 = ix(1,l)
    i2 = ix(2,l)
    x1 = real( x(1,l), p_k_fld )
    x2 = real( x(2,l), p_k_fld )
    lq = real( q(l) * qnorm, p_k_fld )

    ! get spline weights for r at i,j
    call geometry % spline_r_rectangular( i1, x1, Sr )

    ! get spline weights for t at i,j
    call geometry % spline_t( i2, x2, St )

    ! deposit charge
    do k1 = -1, 1
      do k2 = -1, 1
        rho%f2(1,i1+k1,i2+k2) = rho%f2(1,i1+k1,i2+k2) + lq * Sr(k1) * St(k2)
      enddo
    enddo
  enddo

end subroutine deposit_rho_2d_srect_gr
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Deposit charge using the integral rule on the spline_r - Minkowski
!-------------------------------------------------------------------------------
subroutine deposit_rho_2d_sint_gr( rho, geometry, ix, x, q, np )

  use m_vdf_define, only : t_vdf

  implicit none

  integer, parameter :: rank = 2

  type( t_vdf ), intent(inout) :: rho
  type( t_geometry_gr ), intent(in), target :: geometry
  integer, dimension(:,:), intent(in) :: ix
  real(p_k_part), dimension(:,:), intent(in) :: x
  real(p_k_part), dimension( : ), intent(in) :: q
  integer, intent(in) :: np

  ! local variables
  integer :: l, i1, i2
  real(p_k_fld) :: x1, x2, lq
  real(p_k_fld) :: qnorm

  real(p_k_fld), dimension(-1:1) :: Sr
  real(p_k_fld), dimension(-1:1) :: St
  integer :: k1, k2


  ! flat case - integral rule
  qnorm = real( 0.0596831036594608 * (1._p_double + (1._p_double / geometry%delta))**3 /&
               (1._p_double - (1._p_double / geometry%delta**3) ), p_k_fld ) ! 3/16pi


  do l = 1, np

    i1 = ix(1,l)
    i2 = ix(2,l)
    x1 = real( x(1,l), p_k_fld )
    x2 = real( x(2,l), p_k_fld )
    lq = real( q(l) * qnorm, p_k_fld )

    ! get spline weights for r at i,j
    call geometry % spline_r_integral( i1, x1, Sr )

    ! get spline weights for t at i,j
    call geometry % spline_t( i2, x2, St )

    ! deposit charge
    do k1 = -1, 1
      do k2 = -1, 1
        rho%f2(1,i1+k1,i2+k2) = rho%f2(1,i1+k1,i2+k2) + lq * Sr(k1) * St(k2)
      enddo
    enddo
  enddo

end subroutine deposit_rho_2d_sint_gr
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Deposit charge using the rectangular rule on the spline_r - Schwarzschild
!-------------------------------------------------------------------------------
subroutine deposit_rho_2d_srect_curved_gr( rho, geometry, ix, x, q, np )

  use m_vdf_define, only : t_vdf

  implicit none

  integer, parameter :: rank = 2

  type( t_vdf ), intent(inout) :: rho
  type( t_geometry_gr ), intent(in), target :: geometry
  integer, dimension(:,:), intent(in) :: ix
  real(p_k_part), dimension(:,:), intent(in) :: x
  real(p_k_part), dimension( : ), intent(in) :: q
  integer, intent(in) :: np

  ! local variables
  integer :: l, i1, i2
  real(p_k_fld) :: x1, x2, lq
  real(p_k_fld) :: qnorm

  real(p_k_fld), dimension(-1:1) :: Sr
  real(p_k_fld), dimension(-1:1) :: St
  integer :: k1, k2

  ! flat and curved case - rectangular rule
  qnorm = real( 0.0795774715459477 * geometry%A, p_k_fld ) ! 1/4pi

  do l = 1, np

    i1 = ix(1,l)
    i2 = ix(2,l)
    x1 = real( x(1,l), p_k_fld )
    x2 = real( x(2,l), p_k_fld )
    lq = real( q(l) * qnorm, p_k_fld )

    ! get spline weights for r at i,j
    call geometry % spline_r_rectangular_gr( i1, x1, Sr )

    ! get spline weights for t at i,j
    call geometry % spline_t( i2, x2, St )

    ! deposit charge
    do k1 = -1, 1
      do k2 = -1, 1
        rho%f2(1,i1+k1,i2+k2) = rho%f2(1,i1+k1,i2+k2) + lq * Sr(k1) * St(k2)
      enddo
    enddo
  enddo

end subroutine deposit_rho_2d_srect_curved_gr
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Deposit charge using the integral rule on the spline_r - Schwarzschild
!-------------------------------------------------------------------------------
subroutine deposit_rho_2d_sint_curved_gr( rho, geometry, ix, x, q, np )

  use m_vdf_define, only : t_vdf

  implicit none

  integer, parameter :: rank = 2

  type( t_vdf ), intent(inout) :: rho
  type( t_geometry_gr ), intent(in), target :: geometry
  integer, dimension(:,:), intent(in) :: ix
  real(p_k_part), dimension(:,:), intent(in) :: x
  real(p_k_part), dimension( : ), intent(in) :: q
  integer, intent(in) :: np

  ! local variables
  integer :: l, i1, i2
  real(p_k_fld) :: x1, x2, lq
  real(p_k_fld) :: qnorm

  real(p_k_fld), dimension(-1:1) :: Sr
  real(p_k_fld), dimension(-1:1) :: St
  integer :: k1, k2

  ! curved case - integral rule
  qnorm = real( 0.1591549430918953 , p_k_fld ) ! 1/2pi

  do l = 1, np

    i1 = ix(1,l)
    i2 = ix(2,l)
    x1 = real( x(1,l), p_k_fld )
    x2 = real( x(2,l), p_k_fld )
    lq = real( q(l) * qnorm, p_k_fld )

    ! get spline weights for r at i,j
    call geometry % spline_r_integral_gr( i1, x1, Sr )

    ! get spline weights for t at i,j
    call geometry % spline_t( i2, x2, St )

    ! deposit charge
    do k1 = -1, 1
      do k2 = -1, 1
        rho%f2(1,i1+k1,i2+k2) = rho%f2(1,i1+k1,i2+k2) + lq * Sr(k1) * St(k2)
      enddo
    enddo
  enddo

end subroutine deposit_rho_2d_sint_curved_gr
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Deposit charge at half grid positions - Both Minkowski and Schwarzschild
!-------------------------------------------------------------------------------
subroutine deposit_rho_2d_shalf_gr( rho, geometry, ix, x, q, np )

  use m_vdf_define, only : t_vdf

  implicit none

  integer, parameter :: rank = 2

  type( t_vdf ), intent(inout) :: rho
  type( t_geometry_gr ), intent(in), target :: geometry
  integer, dimension(:,:), intent(in) :: ix
  real(p_k_part), dimension(:,:), intent(in) :: x
  real(p_k_part), dimension( : ), intent(in) :: q
  integer, intent(in) :: np

  ! local variables
  integer :: l, i1, i2
  real(p_k_fld) :: x1, x2, lq
  real(p_k_fld) :: qnorm

  real(p_k_fld), dimension(-1:1) :: Sr
  real(p_k_fld), dimension(-1:1) :: St
  integer :: k1, k2

  ! flat and curved case - rectangular rule
  qnorm = real( 0.0795774715459477 * geometry%A, p_k_fld ) ! 1/4pi

  do l = 1, np

    i1 = ix(1,l)
    i2 = ix(2,l)
    x1 = real( x(1,l), p_k_fld )
    x2 = real( x(2,l), p_k_fld )
    lq = real( q(l) * qnorm, p_k_fld )

    ! get spline weights for r at i+1/2,j+1/2
    call geometry % spline_rh( i1, x1, Sr )

    ! get spline weights for t at i+1/2,j+1/2
    call geometry % spline_th( i2, x2, St )

    ! deposit charge
    do k1 = -1, 1
      do k2 = -1, 1
        rho%f2(1,i1+k1,i2+k2) = rho%f2(1,i1+k1,i2+k2) + lq * Sr(k1) * St(k2)
      enddo
    enddo
  enddo

end subroutine deposit_rho_2d_shalf_gr
!-------------------------------------------------------------------------------

end module m_species_charge_gr

!-------------------------------------------------------------------------------
subroutine deposit_density_species_gr( this, rho, i1, i2, q )

  use m_system
  use m_parameters

  use m_species_define_gr, only : t_species_gr
  use m_vdf_define, only : t_vdf
  use m_species_charge_gr
  use m_geometry_gr

  implicit none

  ! dummy variables
  class( t_species_gr ), intent(in) :: this
  type( t_vdf ), intent(inout) :: rho
  integer, intent(in) :: i1, i2
  real(p_k_part), dimension(:), intent(in) :: q

  ! local variables
  class( t_geometry_gr ), pointer :: g
  integer :: np
  g => this%geometry

  ! number of particles to deposit
  np = i2 - i1 + 1

  ! deposit given density
  select case ( this%interpolation )
  case (p_linear)

      select case ( p_x_dim )
      case (2)

        select case ( g%spline_rule )
        case ( p_spline_rectangular )

          select case ( g%metric )
          case ( p_geometry_minkowski )
            call deposit_rho_2d_srect_gr( rho, this%geometry, this%ix(:,i1:i2), this%x(:,i1:i2), q, np )
          case default
            call deposit_rho_2d_srect_curved_gr( rho, this%geometry, this%ix(:,i1:i2), this%x(:,i1:i2), q, np )
          end select

        case ( p_spline_integral )

          select case ( g%metric )
          case ( p_geometry_minkowski )
            call deposit_rho_2d_sint_gr( rho, this%geometry, this%ix(:,i1:i2), this%x(:,i1:i2), q, np )
          case default
            call deposit_rho_2d_sint_curved_gr( rho, this%geometry, this%ix(:,i1:i2), this%x(:,i1:i2), q, np )
          end select

        case default

          call deposit_rho_2d_shalf_gr( rho, this%geometry, this%ix(:,i1:i2), this%x(:,i1:i2), q, np )

        end select

      end select

  case default
      write(err_buf__,*) 'Not implemented';call err__("gr/os-species-charge-gr.f03",383)
      call abort_program( p_err_notimplemented )

  end select


end subroutine deposit_density_species_gr
!-------------------------------------------------------------------------------
