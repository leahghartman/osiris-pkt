#include "os-preprocess.fpp"
#include "os-config.h"

module m_emf_pec_gr

#include "memory/memory.h"

use m_parameters
use m_geometry_gr

implicit none

private

interface update_b_2d_pec
  module procedure update_b_2d_pec
end interface

interface update_e_2d_pec
  module procedure update_e_2d_pec
end interface

public :: update_b_2d_pec, update_e_2d_pec

contains

!-------------------------------------------------------------------------------
subroutine update_b_2d_pec( geometry, b_vdf, bnd )

  use m_vdf_define, only : t_vdf

  implicit none

  ! dummy variables
  class( t_geometry_gr ), intent(inout) :: geometry
  class ( t_vdf ), intent( inout ) :: b_vdf
  integer, intent( in ) :: bnd

  integer :: i1, i2
  real( p_k_fld ), pointer :: b(:,:,:)

  b => b_vdf%f2

  select case (bnd)

  case (p_lower)
    ! inner radial boundary
    do i2 = 0, b_vdf%nx_(p_tc_dim) + 1
      ! br (at the surface):
      b( 1, 1, i2 ) = 0.

      ! inside the conductor:
      b( 1, 0, i2 ) = 0.
      b( 2, 0, i2 ) = 0.
      b( 3, 0, i2 ) = 0.
    enddo
  case (p_upper)
    ! outer radial boundary
    i1 = b_vdf%nx_(p_rc_dim)

    do i2 = 0, b_vdf%nx_(p_tc_dim)
      ! br (at the surface):
      b( 1, i1+1, i2) = 0.

      ! inside the conductor:
      b( 2, i1+1, i2) = 0.
      b( 3, i1+1, i2) = 0.
    enddo
  end select

  ! clear pointers
  nullify(b)

end subroutine

!-------------------------------------------------------------------------------
subroutine update_e_2d_pec( geometry, e_vdf, bnd )

  use m_vdf_define, only : t_vdf

  implicit none

  ! dummy variables
  class( t_geometry_gr ), intent(inout) :: geometry
  class ( t_vdf ), intent( inout ) :: e_vdf
  integer, intent( in ) :: bnd

  integer :: i1, i2
  real( p_k_fld ), pointer :: e(:,:,:)

  e => e_vdf%f2

  select case (bnd)

  case (p_lower)
    ! inner radial boundary
    do i2 = 0, e_vdf%nx_(p_tc_dim) + 1
      ! et (at the surface)
      e( 2, 1, i2 ) = 0.

      ! ep
      e( 3, 1, i2 ) = 0.

      ! inside the conductor:
      e( 1, 0, i2 ) = 0.
      e( 2, 0, i2 ) = 0.
      e( 3, 0, i2 ) = 0.
    enddo
    
  case (p_upper)
    ! outer radial boundary
    i1 = e_vdf%nx_(p_rc_dim)

    do i2 = 0, e_vdf%nx_(p_tc_dim) + 1
      ! et (at the surface):
      e( 2, i1+1, i2 ) = 0.

      ! ep
      e( 3, i1+1, i2 ) = 0.

      ! inside the conductor:
      e( 1, i1+1, i2 ) = 0.
    enddo
  end select

  ! clear pointers
  nullify(e)

end subroutine

end module m_emf_pec_gr