#include "os-preprocess.fpp"
#include "os-config.h"

module m_emf_axi_gr

#include "memory/memory.h"

use m_parameters
use m_geometry_gr

implicit none

private

interface update_b_2d_axi
  module procedure update_b_2d_axi
end interface

interface update_e_2d_axi
  module procedure update_e_2d_axi
end interface

public :: update_b_2d_axi, update_e_2d_axi

contains

!-------------------------------------------------------------------------------
subroutine update_b_2d_axi( geometry, b_vdf, bnd )

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
    ! northern polar boundary
    do i1 = 0, b_vdf%nx_(p_rc_dim) + 1
      ! br (outside domain, mirror first component)
      b( 1, i1, 0 ) = b( 1, i1, 1 )

      ! bt (on boundary, set to zero)
      b( 2, i1, 1 ) = 0.

      ! bt (outside domain, mirror second component)
      b( 2, i1, 0 ) = b( 2, i1, 2 )

      ! bp (outside domain, mirror first component)
      b( 3, i1, 0 ) = -b( 3, i1, 1 )

    enddo

    ! handle inner radial corner
    b( 1, 0, 0 ) =  b( 1, 0, 1 )
    b( 2, 0, 0 ) = -b( 2, 0, 2 )
    b( 3, 0, 0 ) = -b( 3, 0, 1 )
    
  case (p_upper)
    ! southern polar boundary
    i2 = b_vdf%nx_(p_tc_dim)
    do i1 = 1, b_vdf%nx_(p_rc_dim) + 1
      ! br (outside domain, mirror last component)
      b( 1, i1, i2+1 ) = b( 1, i1, i2 )

      ! bt (on boundary, set to zero)
      b( 1, i1, i2+1 ) = 0.

      ! bp (outside domain, mirror last component)
      b( 3, i1, i2+1 ) = -b( 3, i1, i2 )

    enddo
  end select

  ! clear pointers
  nullify(b)

end subroutine

!-------------------------------------------------------------------------------
subroutine update_e_2d_axi( geometry, e_vdf, bnd )

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
    ! northern polar boundary
    do i1 = 1, e_vdf%nx_(p_rc_dim) + 1
      ! er (outside domain, mirror second component)
      e( 1, i1, 0 ) = e( 1, i1, 2 )

      ! et (outside domain, mirror first component)
      e( 2, i1, 0 ) = -e( 2, i1, 1 )

      ! ep (on boundary, set to zero)
      e( 3, i1, 1 ) = 0.

      ! ep (outside domain, mirror second component)
      e( 3, i1, 0 ) = -e( 3, i1, 2 )
    enddo
    
  case (p_upper)
    ! southern polar boundary
    i2 = e_vdf%nx_(p_tc_dim)
    do i1 = 1, e_vdf%nx_(p_rc_dim) + 1
      ! et (outside domain, mirror last component)
      e( 2, i1, i2+1 ) = -e( 2, i1, i2 )

      ! ep (on boundary, set to zero)
      e( 3, i1, i2+1 ) = 0.
    enddo
  end select

  ! clear pointers
  nullify(e)

end subroutine

end module m_emf_axi_gr