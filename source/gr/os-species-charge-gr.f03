#include "os-preprocess.fpp"
#include "os-config.h"

module m_species_charge_gr

#include "memory/memory.h"

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
  real(p_k_fld), dimension(-1:1)  :: St
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
  real(p_k_fld), dimension(-1:1)  :: St
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
  real(p_k_fld), dimension(-1:1)  :: St
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
  real(p_k_fld), dimension(-1:1)  :: St
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
  real(p_k_fld), dimension(-1:1)  :: St
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
  type( t_vdf ),     intent(inout) :: rho
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
      ERROR('Not implemented')
      call abort_program( p_err_notimplemented )

  end select
    

end subroutine deposit_density_species_gr
!-------------------------------------------------------------------------------