#include "os-preprocess.fpp"
#include "os-config.h"

module m_emf_mur_gr

#include "memory/memory.h"

use m_parameters
use m_vdf_define, only: t_vdf
use m_geometry_gr

implicit none

private

type :: t_mur_gr

  real( p_double ) :: dt
  class( t_geometry_gr ), pointer :: geometry
  type( t_vdf ) :: fsave

end type t_mur_gr

interface setup
  module procedure setup_mur
end interface

interface cleanup
  module procedure cleanup_mur
end interface

interface update_b_2d_mur
  module procedure update_b_2d_mur
end interface

interface update_e_2d_mur
  module procedure update_e_2d_mur
end interface

public :: t_mur_gr
public :: setup, cleanup
public :: update_b_2d_mur, update_e_2d_mur

contains

!-------------------------------------------------------------------------------
subroutine setup_mur( this, geometry, grid, dt, restart, restart_handle )

    use m_grid_define, only : t_grid
    use m_grid, only : my_nx
    use m_restart
  
    implicit none
  
    ! dummy variables
    class ( t_mur_gr ), intent( inout ) :: this
    class( t_geometry_gr ), intent( inout ), pointer :: geometry
    class( t_grid ), intent(in) :: grid
    real( p_double ), intent( in ) :: dt
    logical, intent(in) :: restart
    type( t_restart_handle ), intent(in) :: restart_handle
  
    ! local variables
    integer :: nt
    integer, dimension(p_x_dim) :: nx
    integer, dimension(2, 2) :: gc_num_f2
    real( p_double ), dimension(2) :: dx_vec
  
    if ( restart ) then
      ! do nothing
      ! fsave was already read from restart file
    else
      ! executable statements
      nx = my_nx(grid)
      nt = nx(p_tc_dim)

      ! guard cells for f2s
      ! should not be required, but are set for completeness
      gc_num_f2(:, 1) = 0
      gc_num_f2(p_lower, 2) = 1
      gc_num_f2(p_upper, 2) = 2

      ! dummy array to create vdfs
      dx_vec = (/0,0/)
      
      ! allocate vdf to store boundary fields
      ! x_dim = 2
      ! f_dim = 3 (1: br, 2: et, 3: ep)
      ! nx = 2 x nt cells
      call this%fsave%new( 2, 3, (/2, nt/), gc_num_f2, dx_vec, .true. )
    endif

    ! save dt
    this%dt = dt

    ! save pointer to geometry
    this%geometry => geometry

  end subroutine

!-------------------------------------------------------------------------------
subroutine cleanup_mur( this )

  implicit none

  ! dummy variables
  class ( t_mur_gr ), intent( inout )  ::  this

  ! local variables
  
  ! executable statements
  call this%fsave%cleanup()

end subroutine

!-------------------------------------------------------------------------------
subroutine update_b_2d_mur( this, b_vdf, step )

  use m_vdf_define, only : t_vdf

  implicit none

  ! dummy variables
  class ( t_mur_gr ), intent( inout )  ::  this
  class ( t_vdf ), intent( inout ) :: b_vdf
  integer, intent( in ) :: step

  integer :: i1, i2
  real(p_double) :: dr, dt
  real( p_k_fld ), pointer :: r(:,:)
  real( p_k_fld ), pointer :: b(:,:,:), fsave(:, :, :)
  real( p_k_fld ), pointer :: alpha(:,:)
  
  b     => b_vdf%f2
  r     => this%geometry%r%f1
  alpha => this%geometry%alpha%f1
  fsave => this%fsave%f2

  ! get index of last cell
  i1 = b_vdf%nx_(p_rc_dim)

  dt = this%dt
  dr = r(1, i1+1) - r(1, i1)
  
  do i2 = 0, b_vdf%nx_(p_tc_dim)
    
    if ( step == 1 ) then
      b(1, i1+1, i2) = &
        - alpha(1, i1+1) * b(1, i1, i2) / alpha(1, i1) - &
        fsave(1, 2, i2) * ( &
          - 1. + &
          alpha(1, i1+1)**2 * dt / dr + &
          alpha(1, i1+1) * dt / (2. * r(1, i1+1)) &
        ) - &
        fsave(1, 1, i2) * ( &
          - alpha( 1, i1+1 ) / alpha( 1, i1 ) - &
          alpha(1, i1+1) * alpha(1, i1) * dt / dr + &
          alpha(1, i1+1) * dt / (2. * r(1, i1)) &
        )

    else
      b(1, i1+1, i2) = &
        ( &
          fsave(1, 2, i2) / alpha(1, i1+1) + &
          fsave(1, 1, i2) / alpha(1, i1) - &
          b( 1, i1, i2) * ( &
            1. / alpha( 1, i1 ) - &
            dt * alpha(1, i1) / dr + &
            dt / (2. * r(1, i1)) &
          ) &
        ) / &
        ( &
          1. / alpha(1, i1+1 ) + &
          alpha( 1, i1+1 ) * dt / dr + &
          dt / (2. * r(1, i1+1)) &
        )

    end if

    ! inside the conductor:
    b(2, i1+1, i2) = 0.
    b(3, i1+1, i2) = 0.

    ! save values
    fsave(1, 1, i2) = b(1, i1  , i2)
    fsave(1, 2, i2) = b(1, i1+1, i2)

  enddo

  ! clear pointers
  nullify(b)
  nullify(r)
  nullify(alpha)
  nullify(fsave)


end subroutine

!-------------------------------------------------------------------------------
subroutine update_e_2d_mur( this, e_vdf )

  use m_vdf_define, only : t_vdf

  implicit none

  ! dummy variables
  class ( t_mur_gr ), intent( inout )  ::  this
  class ( t_vdf ), intent( inout ) :: e_vdf

  integer :: i1, i2
  real(p_double) :: dr, dt
  real( p_k_fld ), pointer :: r(:,:)
  real( p_k_fld ), pointer :: e(:,:,:), fsave(:, :, :)
  real( p_k_fld ), pointer :: alpha(:,:)
  
  e     => e_vdf%f2
  r     => this%geometry%r%f1
  alpha => this%geometry%alpha%f1
  fsave => this%fsave%f2

  ! get index of last cell
  i1 = e_vdf%nx_(p_rc_dim)

  dt = this%dt
  dr = r(1, i1+1) - r(1, i1)
  
  do i2 = 0, e_vdf%nx_(p_tc_dim)

  
    e(2, i1+1, i2) = ( &
        - ( &
          - 1. / (alpha(1, i1+1) * dt) + &
          alpha(1, i1+1) / dr + &
          1. / (2. * r(1, i1+1)) &
        ) * fsave(2, 2, i2) - &
        ( &
          1. / (alpha(1, i1) * dt) - &
          alpha(1, i1) / dr + &
          1. / (2. * r(1, i1)) &
        ) * e(2, i1, i2 ) - &
        ( &
          - 1. / (alpha(1, i1) * dt) - &
          alpha(1, i1) / dr + &
          1. / (2. * r(1, i1)) &
        ) * fsave(2, 1, i2) &
      ) / &
      ( &
        1. / (alpha(1, i1+1) * dt) + &
        alpha(1, i1+1) / dr + &
        1. / (2. * r(1, i1+1)) &
      )

    e(3, i1+1, i2) = ( &
        - ( &
          - 1. / (alpha(1, i1+1) * dt) + &
          alpha(1, i1+1) / dr + &
          1. / (2. * r(1, i1+1)) &
        ) * fsave(3, 2, i2) - &
        ( &
          1. / (alpha(1, i1) * dt) - &
          alpha(1, i1) / dr + &
          1. / (2. * r(1, i1)) &
        ) * e(3, i1, i2 ) - &
        ( &
          - 1. / (alpha(1, i1) * dt) - &
          alpha(1, i1) / dr + &
          1. / (2. * r(1, i1)) &
        ) * fsave(3, 1, i2) &
      ) / &
      ( &
        1. / (alpha(1, i1+1) * dt) + &
        alpha(1, i1+1) / dr + &
        1. / (2. * r(1, i1+1)) &
      )

    ! fill ghost cells outside the simulation domain
    e(1, i1+1, i2) = 0.

    ! save values
    fsave(2, 1, i2) = e(2, i1  , i2)
    fsave(2, 2, i2) = e(2, i1+1, i2)
    fsave(3, 1, i2) = e(3, i1  , i2)
    fsave(3, 2, i2) = e(3, i1+1, i2)

  enddo

  ! clear pointers
  nullify(e)
  nullify(r)
  nullify(alpha)
  nullify(fsave)


end subroutine

end module m_emf_mur_gr
