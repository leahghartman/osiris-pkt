module m_term_chars

implicit none

! terminal characters for grammar description of 2d space filling curves
public :: up, right, down, left
! terminal characters for grammar description of 3d space filling curves
public :: ux, dx, uy, dy, uz, dz

contains

!-------------------------------------------------------------------------------
! moves up one cell and fills in the correct value
!-------------------------------------------------------------------------------
subroutine up( hi2, x, t )

  implicit none

  integer, dimension(:,:), intent(inout) :: hi2
  integer, dimension(2), intent(inout) :: x
  integer, intent(inout) :: t

  x(2) = x(2) + 1
  t = t+1
  hi2( x(1), x(2) ) = t

end subroutine up
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! moves right one cell and fills in the correct value
!-------------------------------------------------------------------------------
subroutine right( hi2, x, t )

  implicit none

  integer, dimension(:,:), intent(inout) :: hi2
  integer, dimension(2), intent(inout) :: x
  integer, intent(inout) :: t

  x(1) = x(1) + 1
  t = t+1
  hi2( x(1), x(2) ) = t

end subroutine right
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! moves down one cell and fills in the correct value
!-------------------------------------------------------------------------------
subroutine down( hi2, x, t )

  implicit none

  integer, dimension(:,:), intent(inout) :: hi2
  integer, dimension(2), intent(inout) :: x
  integer, intent(inout) :: t

  x(2) = x(2) - 1
  t = t+1
  hi2( x(1), x(2) ) = t

end subroutine down
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! moves left one cell and fills in the correct value
!-------------------------------------------------------------------------------
subroutine left( hi2, x, t )

  implicit none

  integer, dimension(:,:), intent(inout) :: hi2
  integer, dimension(2), intent(inout) :: x
  integer, intent(inout) :: t

  x(1) = x(1) - 1
  t = t+1
  hi2( x(1), x(2) ) = t

end subroutine left
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! moves up one cell in x direction
!-------------------------------------------------------------------------------
subroutine ux( hi3, x, t )

  implicit none

  integer, dimension(:,:,:), intent(inout) :: hi3
  integer, dimension(3), intent(inout) :: x
  integer, intent(inout) :: t

  x(1) = x(1) + 1
  t = t+1
  hi3( x(1), x(2), x(3) ) = t

end subroutine ux
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! moves down one cell in x direction
!-------------------------------------------------------------------------------
subroutine dx( hi3, x, t )

  implicit none

  integer, dimension(:,:,:), intent(inout) :: hi3
  integer, dimension(3), intent(inout) :: x
  integer, intent(inout) :: t

  x(1) = x(1) - 1
  t = t+1
  hi3( x(1), x(2), x(3) ) = t

end subroutine dx
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! moves up one cell in y direction
!-------------------------------------------------------------------------------
subroutine uy( hi3, x, t )

  implicit none

  integer, dimension(:,:,:), intent(inout) :: hi3
  integer, dimension(3), intent(inout) :: x
  integer, intent(inout) :: t

  x(2) = x(2) + 1
  t = t+1
  hi3( x(1), x(2), x(3) ) = t

end subroutine uy
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! moves down one cell in y direction
!-------------------------------------------------------------------------------
subroutine dy( hi3, x, t )

  implicit none

  integer, dimension(:,:,:), intent(inout) :: hi3
  integer, dimension(3), intent(inout) :: x
  integer, intent(inout) :: t

  x(2) = x(2) - 1
  t = t+1
  hi3( x(1), x(2), x(3) ) = t

end subroutine dy
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! moves up one cell in z direction
!-------------------------------------------------------------------------------
subroutine uz( hi3, x, t )

  implicit none

  integer, dimension(:,:,:), intent(inout) :: hi3
  integer, dimension(3), intent(inout) :: x
  integer, intent(inout) :: t

  x(3) = x(3) + 1
  t = t+1
  hi3( x(1), x(2), x(3) ) = t

end subroutine uz
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! moves down one cell in z direction
!-------------------------------------------------------------------------------
subroutine dz( hi3, x, t )

  implicit none

  integer, dimension(:,:,:), intent(inout) :: hi3
  integer, dimension(3), intent(inout) :: x
  integer, intent(inout) :: t

  x(3) = x(3) - 1
  t = t+1
  hi3( x(1), x(2), x(3) ) = t

end subroutine dz
!-------------------------------------------------------------------------------

end module m_term_chars
