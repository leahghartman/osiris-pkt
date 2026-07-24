! #include "os-config.h"
#include "os-preprocess.fpp"

module m_hilbert

! #include "memory/memory.h"

use m_term_chars
use m_system

implicit none

public :: fill_hi2_hilbert
public :: fill_hi3_hilbert

private :: H, A, B, C
private :: xyz, xy_z_, yzx, yz_x_, zxy, zx_y_, x_yz_, x_y_z, y_zx_, y_z_x, z_xy_, z_x_y
private :: SH, SH_

contains

!-------------------------------------------------------------------------------
! non-terminal character for 2d hilbert grammar
!-------------------------------------------------------------------------------
recursive subroutine H( hi2, x, t, n )

  integer, dimension(:,:), intent(inout) :: hi2
  integer, dimension(2), intent(inout) :: x
  integer, intent(inout) :: t
  integer, intent(in) :: n

  if ( n==2 ) then
    call up(hi2,x,t)
    call right(hi2,x,t)
    call down(hi2,x,t)
  else
    call A(hi2,x,t,n/2)
    call up(hi2,x,t)
    call H(hi2,x,t,n/2)
    call right(hi2,x,t)
    call H(hi2,x,t,n/2)
    call down(hi2,x,t)
    call B(hi2,x,t,n/2)
  endif

end subroutine H
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! non-terminal character for 2d hilbert grammar
!-------------------------------------------------------------------------------
recursive subroutine A( hi2, x, t, n )

  integer, dimension(:,:), intent(inout) :: hi2
  integer, dimension(2), intent(inout) :: x
  integer, intent(inout) :: t
  integer, intent(in) :: n

  if ( n==2) then
    call right(hi2,x,t)
    call up(hi2,x,t)
    call left(hi2,x,t)
  else
    call H(hi2,x,t,n/2)
    call right(hi2,x,t)
    call A(hi2,x,t,n/2)
    call up(hi2,x,t)
    call A(hi2,x,t,n/2)
    call left(hi2,x,t)
    call C(hi2,x,t,n/2)
  endif

end subroutine A
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! non-terminal character for 2d hilbert grammar
!-------------------------------------------------------------------------------
recursive subroutine B( hi2, x, t, n )

  integer, dimension(:,:), intent(inout) :: hi2
  integer, dimension(2), intent(inout) :: x
  integer, intent(inout) :: t
  integer, intent(in) :: n

  if ( n==2 ) then
    call left(hi2,x,t)
    call down(hi2,x,t)
    call right(hi2,x,t)
  else
    call C(hi2,x,t,n/2)
    call left(hi2,x,t)
    call B(hi2,x,t,n/2)
    call down(hi2,x,t)
    call B(hi2,x,t,n/2)
    call right(hi2,x,t)
    call H(hi2,x,t,n/2)
  endif

end subroutine B
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! non-terminal character for 2d hilbert grammar
!-------------------------------------------------------------------------------
recursive subroutine C( hi2, x, t, n )

  integer, dimension(:,:), intent(inout) :: hi2
  integer, dimension(2), intent(inout) :: x
  integer, intent(inout) :: t
  integer, intent(in) :: n

  if ( n==2 ) then
    call down(hi2,x,t)
    call left(hi2,x,t)
    call up(hi2,x,t)
  else
    call B(hi2,x,t,n/2)
    call down(hi2,x,t)
    call C(hi2,x,t,n/2)
    call left(hi2,x,t)
    call C(hi2,x,t,n/2)
    call up(hi2,x,t)
    call A(hi2,x,t,n/2)
  endif

end subroutine C
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! fill hi2 with g_til_aid's in hilbert curve order
! Rectangular simulation box sizes are allowed but must satisfy:
!   - smallest dimension must have number of tiles equal to a power of two
!   - all other dimesions must be a multiple of the smallest dimension
!-------------------------------------------------------------------------------
subroutine fill_hi2_hilbert( hi2, n, t )

  implicit none

  integer, dimension(:,:), intent(inout) :: hi2
  integer, dimension(:), intent(in)      :: n
  integer, intent(inout)                 :: t

  integer, dimension(2) :: x   ! starting position
  integer :: i, lb, ub

  if ( IAND( minval(n), minval(n)-1 )/=0 ) then
    ERROR('Number of tiles in smallest dimension must be a power of 2')
    call abort_program()
  elseif ( minval(n)==1 .and. any(n /= 1) ) then
    ERROR('If number of tiles in smallest dimension is 1,')
    ERROR('all number of tiles must be 1. Try topology="snake"')
    call abort_program()
  endif

  if ( n(1)==n(2) ) then
    if ( n(1)==1 ) then
      hi2 = 1
    else
      x   = 1
      hi2 = t
      call H( hi2,x,t,n(1) )
    endif
  elseif ( (minval(n)==n(1)) .and. (mod(n(2),n(1))==0) ) then
    do i = 1, n(2)/n(1)
      lb  = (i-1)*n(1)+1
      ub  = i*n(1)
      x   = 1
      hi2(:,lb:ub) = n(1)*n(1)*(i-1) + 1
      call A( hi2(:,lb:ub),x,t,n(1) )
      t = t+1
    enddo
  elseif ( (minval(n)==n(2)) .and. (mod(n(1),n(2))==0) ) then
    do i = 1, n(1)/n(2)
      lb  = (i-1)*n(2)+1
      ub  = i*n(2)
      x   = 1
      hi2(lb:ub,:) = n(2)*n(2)*(i-1) + 1
      call H( hi2(lb:ub,:),x,t,n(2) )
      t = t+1
    enddo
  else
    ERROR('Number of tiles in largest dimension must be multiple of number')
    ERROR('of tiles in smallest dimension')
    call abort_program()
  endif

end subroutine fill_hi2_hilbert
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! non-terminal character for 3d hilbert grammar
!-------------------------------------------------------------------------------
recursive subroutine xyz( hi3, x, t, n )

  integer, dimension(:,:,:), intent(inout) :: hi3
  integer, dimension(3), intent(inout) :: x
  integer, intent(inout) :: t
  integer, intent(in) :: n

  if ( n==2 ) then
    call uy(hi3,x,t)
    call ux(hi3,x,t)
    call dy(hi3,x,t)
    call uz(hi3,x,t)
    call uy(hi3,x,t)
    call dx(hi3,x,t)
    call dy(hi3,x,t)
  else
		call yzx(hi3,x,t,n/2)
    call uy(hi3,x,t)
		call zxy(hi3,x,t,n/2)
    call ux(hi3,x,t)
		call zxy(hi3,x,t,n/2)
    call dy(hi3,x,t)
		call x_y_z(hi3,x,t,n/2)
    call uz(hi3,x,t)
		call x_y_z(hi3,x,t,n/2)
    call uy(hi3,x,t)
		call z_xy_(hi3,x,t,n/2)
    call dx(hi3,x,t)
		call z_xy_(hi3,x,t,n/2)
    call dy(hi3,x,t)
		call yz_x_(hi3,x,t,n/2)
  endif

end subroutine xyz
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! non-terminal character for 3d hilbert grammar
!-------------------------------------------------------------------------------
recursive subroutine xy_z_( hi3, x, t, n )

  integer, dimension(:,:,:), intent(inout) :: hi3
  integer, dimension(3), intent(inout) :: x
  integer, intent(inout) :: t
  integer, intent(in) :: n

  if ( n==2 ) then
    call dy(hi3,x,t)
    call ux(hi3,x,t)
    call uy(hi3,x,t)
    call dz(hi3,x,t)
    call dy(hi3,x,t)
    call dx(hi3,x,t)
    call uy(hi3,x,t)
  else
		call yz_x_(hi3,x,t,n/2)
    call dy(hi3,x,t)
		call zx_y_(hi3,x,t,n/2)
    call ux(hi3,x,t)
		call zx_y_(hi3,x,t,n/2)
    call uy(hi3,x,t)
		call x_yz_(hi3,x,t,n/2)
    call dz(hi3,x,t)
		call x_yz_(hi3,x,t,n/2)
    call dy(hi3,x,t)
		call z_x_y(hi3,x,t,n/2)
    call dx(hi3,x,t)
		call z_x_y(hi3,x,t,n/2)
    call uy(hi3,x,t)
		call yzx(hi3,x,t,n/2)
  endif

end subroutine xy_z_
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! non-terminal character for 3d hilbert grammar
!-------------------------------------------------------------------------------
recursive subroutine yzx( hi3, x, t, n )

  integer, dimension(:,:,:), intent(inout) :: hi3
  integer, dimension(3), intent(inout) :: x
  integer, intent(inout) :: t
  integer, intent(in) :: n

  if ( n==2 ) then
    call ux(hi3,x,t)
    call uz(hi3,x,t)
    call dx(hi3,x,t)
    call uy(hi3,x,t)
    call ux(hi3,x,t)
    call dz(hi3,x,t)
    call dx(hi3,x,t)
  else
		call zxy(hi3,x,t,n/2)
    call ux(hi3,x,t)
		call xyz(hi3,x,t,n/2)
    call uz(hi3,x,t)
		call xyz(hi3,x,t,n/2)
    call dx(hi3,x,t)
		call y_zx_(hi3,x,t,n/2)
    call uy(hi3,x,t)
		call y_zx_(hi3,x,t,n/2)
    call ux(hi3,x,t)
		call xy_z_(hi3,x,t,n/2)
    call dz(hi3,x,t)
		call xy_z_(hi3,x,t,n/2)
    call dx(hi3,x,t)
		call z_x_y(hi3,x,t,n/2)
  endif

end subroutine yzx
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! non-terminal character for 3d hilbert grammar
!-------------------------------------------------------------------------------
recursive subroutine yz_x_( hi3, x, t, n )

  integer, dimension(:,:,:), intent(inout) :: hi3
  integer, dimension(3), intent(inout) :: x
  integer, intent(inout) :: t
  integer, intent(in) :: n

  if ( n==2 ) then
    call ux(hi3,x,t)
    call dz(hi3,x,t)
    call dx(hi3,x,t)
    call dy(hi3,x,t)
    call ux(hi3,x,t)
    call uz(hi3,x,t)
    call dx(hi3,x,t)
  else
		call zx_y_(hi3,x,t,n/2)
    call ux(hi3,x,t)
		call xy_z_(hi3,x,t,n/2)
    call dz(hi3,x,t)
		call xy_z_(hi3,x,t,n/2)
    call dx(hi3,x,t)
		call y_z_x(hi3,x,t,n/2)
    call dy(hi3,x,t)
		call y_z_x(hi3,x,t,n/2)
    call ux(hi3,x,t)
		call xyz(hi3,x,t,n/2)
    call uz(hi3,x,t)
		call xyz(hi3,x,t,n/2)
    call dx(hi3,x,t)
		call z_xy_(hi3,x,t,n/2)
  endif

end subroutine yz_x_
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! non-terminal character for 3d hilbert grammar
!-------------------------------------------------------------------------------
recursive subroutine zxy( hi3, x, t, n )

  integer, dimension(:,:,:), intent(inout) :: hi3
  integer, dimension(3), intent(inout) :: x
  integer, intent(inout) :: t
  integer, intent(in) :: n

  if ( n==2 ) then
    call uz(hi3,x,t)
    call uy(hi3,x,t)
    call dz(hi3,x,t)
    call ux(hi3,x,t)
    call uz(hi3,x,t)
    call dy(hi3,x,t)
    call dz(hi3,x,t)
  else
		call xyz(hi3,x,t,n/2)
    call uz(hi3,x,t)
		call yzx(hi3,x,t,n/2)
    call uy(hi3,x,t)
		call yzx(hi3,x,t,n/2)
    call dz(hi3,x,t)
		call zx_y_(hi3,x,t,n/2)
    call ux(hi3,x,t)
		call zx_y_(hi3,x,t,n/2)
    call uz(hi3,x,t)
		call y_z_x(hi3,x,t,n/2)
    call dy(hi3,x,t)
		call y_z_x(hi3,x,t,n/2)
    call dz(hi3,x,t)
		call x_yz_(hi3,x,t,n/2)
  endif

end subroutine zxy
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! non-terminal character for 3d hilbert grammar
!-------------------------------------------------------------------------------
recursive subroutine zx_y_( hi3, x, t, n )

  integer, dimension(:,:,:), intent(inout) :: hi3
  integer, dimension(3), intent(inout) :: x
  integer, intent(inout) :: t
  integer, intent(in) :: n

  if ( n==2 ) then
    call dz(hi3,x,t)
    call dy(hi3,x,t)
    call uz(hi3,x,t)
    call ux(hi3,x,t)
    call dz(hi3,x,t)
    call uy(hi3,x,t)
    call uz(hi3,x,t)
  else
		call xy_z_(hi3,x,t,n/2)
    call dz(hi3,x,t)
		call yz_x_(hi3,x,t,n/2)
    call dy(hi3,x,t)
		call yz_x_(hi3,x,t,n/2)
    call uz(hi3,x,t)
		call zxy(hi3,x,t,n/2)
    call ux(hi3,x,t)
		call zxy(hi3,x,t,n/2)
    call dz(hi3,x,t)
		call y_zx_(hi3,x,t,n/2)
    call uy(hi3,x,t)
		call y_zx_(hi3,x,t,n/2)
    call uz(hi3,x,t)
		call x_y_z(hi3,x,t,n/2)
  endif

end subroutine zx_y_
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! non-terminal character for 3d hilbert grammar
!-------------------------------------------------------------------------------
recursive subroutine x_yz_( hi3, x, t, n )

  integer, dimension(:,:,:), intent(inout) :: hi3
  integer, dimension(3), intent(inout) :: x
  integer, intent(inout) :: t
  integer, intent(in) :: n

  if ( n==2 ) then
    call uy(hi3,x,t)
    call dx(hi3,x,t)
    call dy(hi3,x,t)
    call dz(hi3,x,t)
    call uy(hi3,x,t)
    call ux(hi3,x,t)
    call dy(hi3,x,t)
  else
		call y_zx_(hi3,x,t,n/2)
    call uy(hi3,x,t)
		call z_xy_(hi3,x,t,n/2)
    call dx(hi3,x,t)
		call z_xy_(hi3,x,t,n/2)
    call dy(hi3,x,t)
		call xy_z_(hi3,x,t,n/2)
    call dz(hi3,x,t)
		call xy_z_(hi3,x,t,n/2)
    call uy(hi3,x,t)
		call zxy(hi3,x,t,n/2)
    call ux(hi3,x,t)
		call zxy(hi3,x,t,n/2)
    call dy(hi3,x,t)
		call y_z_x(hi3,x,t,n/2)
  endif

end subroutine x_yz_
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! non-terminal character for 3d hilbert grammar
!-------------------------------------------------------------------------------
recursive subroutine x_y_z( hi3, x, t, n )

  integer, dimension(:,:,:), intent(inout) :: hi3
  integer, dimension(3), intent(inout) :: x
  integer, intent(inout) :: t
  integer, intent(in) :: n

  if ( n==2 ) then
    call dy(hi3,x,t)
    call dx(hi3,x,t)
    call uy(hi3,x,t)
    call uz(hi3,x,t)
    call dy(hi3,x,t)
    call ux(hi3,x,t)
    call uy(hi3,x,t)
  else
		call y_z_x(hi3,x,t,n/2)
    call dy(hi3,x,t)
		call z_x_y(hi3,x,t,n/2)
    call dx(hi3,x,t)
		call z_x_y(hi3,x,t,n/2)
    call uy(hi3,x,t)
		call xyz(hi3,x,t,n/2)
    call uz(hi3,x,t)
		call xyz(hi3,x,t,n/2)
    call dy(hi3,x,t)
		call zx_y_(hi3,x,t,n/2)
    call ux(hi3,x,t)
		call zx_y_(hi3,x,t,n/2)
    call uy(hi3,x,t)
		call y_zx_(hi3,x,t,n/2)
  endif

end subroutine x_y_z
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! non-terminal character for 3d hilbert grammar
!-------------------------------------------------------------------------------
recursive subroutine y_zx_( hi3, x, t, n )

  integer, dimension(:,:,:), intent(inout) :: hi3
  integer, dimension(3), intent(inout) :: x
  integer, intent(inout) :: t
  integer, intent(in) :: n

  if ( n==2 ) then
    call dx(hi3,x,t)
    call dz(hi3,x,t)
    call ux(hi3,x,t)
    call uy(hi3,x,t)
    call dx(hi3,x,t)
    call uz(hi3,x,t)
    call ux(hi3,x,t)
  else
		call z_xy_(hi3,x,t,n/2)
    call dx(hi3,x,t)
		call x_yz_(hi3,x,t,n/2)
    call dz(hi3,x,t)
		call x_yz_(hi3,x,t,n/2)
    call ux(hi3,x,t)
		call yzx(hi3,x,t,n/2)
    call uy(hi3,x,t)
		call yzx(hi3,x,t,n/2)
    call dx(hi3,x,t)
		call x_y_z(hi3,x,t,n/2)
    call uz(hi3,x,t)
		call x_y_z(hi3,x,t,n/2)
    call ux(hi3,x,t)
		call zx_y_(hi3,x,t,n/2)
  endif

end subroutine y_zx_
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! non-terminal character for 3d hilbert grammar
!-------------------------------------------------------------------------------
recursive subroutine y_z_x( hi3, x, t, n )

  integer, dimension(:,:,:), intent(inout) :: hi3
  integer, dimension(3), intent(inout) :: x
  integer, intent(inout) :: t
  integer, intent(in) :: n

  if ( n==2 ) then
    call dx(hi3,x,t)
    call uz(hi3,x,t)
    call ux(hi3,x,t)
    call dy(hi3,x,t)
    call dx(hi3,x,t)
    call dz(hi3,x,t)
    call ux(hi3,x,t)
  else
		call z_x_y(hi3,x,t,n/2)
    call dx(hi3,x,t)
		call x_y_z(hi3,x,t,n/2)
    call uz(hi3,x,t)
		call x_y_z(hi3,x,t,n/2)
    call ux(hi3,x,t)
		call yz_x_(hi3,x,t,n/2)
    call dy(hi3,x,t)
		call yz_x_(hi3,x,t,n/2)
    call dx(hi3,x,t)
		call x_yz_(hi3,x,t,n/2)
    call dz(hi3,x,t)
		call x_yz_(hi3,x,t,n/2)
    call ux(hi3,x,t)
		call zxy(hi3,x,t,n/2)
  endif

end subroutine y_z_x
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! non-terminal character for 3d hilbert grammar
!-------------------------------------------------------------------------------
recursive subroutine z_xy_( hi3, x, t, n )

  integer, dimension(:,:,:), intent(inout) :: hi3
  integer, dimension(3), intent(inout) :: x
  integer, intent(inout) :: t
  integer, intent(in) :: n

  if ( n==2 ) then
    call dz(hi3,x,t)
    call uy(hi3,x,t)
    call uz(hi3,x,t)
    call dx(hi3,x,t)
    call dz(hi3,x,t)
    call dy(hi3,x,t)
    call uz(hi3,x,t)
  else
		call x_yz_(hi3,x,t,n/2)
    call dz(hi3,x,t)
		call y_zx_(hi3,x,t,n/2)
    call uy(hi3,x,t)
		call y_zx_(hi3,x,t,n/2)
    call uz(hi3,x,t)
		call z_x_y(hi3,x,t,n/2)
    call dx(hi3,x,t)
		call z_x_y(hi3,x,t,n/2)
    call dz(hi3,x,t)
		call yz_x_(hi3,x,t,n/2)
    call dy(hi3,x,t)
		call yz_x_(hi3,x,t,n/2)
    call uz(hi3,x,t)
		call xyz(hi3,x,t,n/2)
  endif

end subroutine z_xy_
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! non-terminal character for 3d hilbert grammar
!-------------------------------------------------------------------------------
recursive subroutine z_x_y( hi3, x, t, n )

  integer, dimension(:,:,:), intent(inout) :: hi3
  integer, dimension(3), intent(inout) :: x
  integer, intent(inout) :: t
  integer, intent(in) :: n

  if ( n==2 ) then
    call uz(hi3,x,t)
    call dy(hi3,x,t)
    call dz(hi3,x,t)
    call dx(hi3,x,t)
    call uz(hi3,x,t)
    call uy(hi3,x,t)
    call dz(hi3,x,t)
  else
		call x_y_z(hi3,x,t,n/2)
    call uz(hi3,x,t)
		call y_z_x(hi3,x,t,n/2)
    call dy(hi3,x,t)
		call y_z_x(hi3,x,t,n/2)
    call dz(hi3,x,t)
		call z_xy_(hi3,x,t,n/2)
    call dx(hi3,x,t)
		call z_xy_(hi3,x,t,n/2)
    call uz(hi3,x,t)
		call yzx(hi3,x,t,n/2)
    call uy(hi3,x,t)
		call yzx(hi3,x,t,n/2)
    call dz(hi3,x,t)
		call xy_z_(hi3,x,t,n/2)
  endif

end subroutine z_x_y
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Non terminal character for hilbert-snake space filling curve
!  fills in rectangularly shaped hi3 by going through hi3 in snake pattern
!  filling in cubes in 3d hilbert curve order
!-------------------------------------------------------------------------------
recursive subroutine SH( hi3, t, lb, ub, n )

  implicit none

  integer, dimension(:,:,:), intent(inout) :: hi3
  integer, intent(inout) :: t      ! g_til_aid parameter
  integer, dimension(2), intent(inout) :: lb, ub ! keeps track of which part of hi3 were on
  integer, dimension(3), intent(in) :: n

  integer :: n_min
  integer, dimension(3) :: x

  select case (minloc(n,1))
    case (1)
      n_min = n(1)
      if ( ub(2)==n(3) ) then
        ! fill in window
        x = 1
        hi3(:,lb(1):ub(1),lb(2):ub(2)) = t
        call yzx( hi3(:,lb(1):ub(1),lb(2):ub(2)), x, t, n_min )
        if ( t==product(n) ) return
        t = t+1
        ! move window up y
        lb(1) = lb(1) + n_min
        ub(1) = ub(1) + n_min
        ! fill in window
        x = 1
        hi3(:,lb(1):ub(1),lb(2):ub(2)) = t
        call yzx( hi3(:,lb(1):ub(1),lb(2):ub(2)), x, t, n_min )
        if ( t==product(n) ) return
        t = t+1
        ! move window down z
        lb(2) = lb(2) - n_min
        ub(2) = ub(2) - n_min
        ! switch to other non-terminal character
        call SH_(hi3,t,lb,ub,n)
      else
        ! fill in window
        x = 1
        hi3(:,lb(1):ub(1),lb(2):ub(2)) = t
        call xyz( hi3(:,lb(1):ub(1),lb(2):ub(2)), x, t, n_min )
        if ( t==product(n) ) return
        t = t+1
        ! move window up z
        lb(2) = lb(2) + n_min
        ub(2) = ub(2) + n_min
        ! call this non-terminal character again
        call SH(hi3,t,lb,ub,n)
      endif
    case (2)
      n_min = n(2)
      if ( ub(2)==n(3) ) then
        ! fill in window
        x = 1
        hi3(lb(1):ub(1),:,lb(2):ub(2)) = t
        call zxy( hi3(lb(1):ub(1),:,lb(2):ub(2)), x, t, n_min )
        if ( t==product(n) ) return
        t = t+1
        ! move window up x
        lb(1) = lb(1) + n_min
        ub(1) = ub(1) + n_min
        ! fill in window
        x = 1
        hi3(lb(1):ub(1),:,lb(2):ub(2)) = t
        call zxy( hi3(lb(1):ub(1),:,lb(2):ub(2)), x, t, n_min )
        if ( t==product(n) ) return
        t = t+1
        ! move window down z
        lb(2) = lb(2) - n_min
        ub(2) = ub(2) - n_min
        ! switch to other non-terminal character
        call SH_(hi3,t,lb,ub,n)
      else
        ! fill in window
        x = 1
        hi3(lb(1):ub(1),:,lb(2):ub(2)) = t
        call xyz( hi3(lb(1):ub(1),:,lb(2):ub(2)), x, t, n_min )
        if ( t==product(n) ) return
        t = t+1
        ! move window up z
        lb(2) = lb(2) + n_min
        ub(2) = ub(2) + n_min
        ! call this non-terminal character again
        call SH(hi3,t,lb,ub,n)
      endif
    case (3)
      n_min = n(3)
      if ( ub(2)==n(2) ) then
        ! fill in window
        x = 1
        hi3(lb(1):ub(1),lb(2):ub(2),:) = t
        call zxy( hi3(lb(1):ub(1),lb(2):ub(2),:), x, t, n_min )
        if ( t==product(n) ) return
        t = t+1
        ! move window up x
        lb(1) = lb(1) + n_min
        ub(1) = ub(1) + n_min
        ! fill in window
        x = 1
        hi3(lb(1):ub(1),lb(2):ub(2),:) = t
        call zxy( hi3(lb(1):ub(1),lb(2):ub(2),:), x, t, n_min )
        if ( t==product(n) ) return
        t = t+1
        ! move window down z
        lb(2) = lb(2) - n_min
        ub(2) = ub(2) - n_min
        ! switch to other non-terminal character
        call SH_(hi3,t,lb,ub,n)
      else
        ! fill in window
        x = 1
        hi3(lb(1):ub(1),lb(2):ub(2),:) = t
        call yzx( hi3(lb(1):ub(1),lb(2):ub(2),:), x, t, n_min )
        if ( t==product(n) ) return
        t = t+1
        ! move window up z
        lb(2) = lb(2) + n_min
        ub(2) = ub(2) + n_min
        ! call this non-terminal character again
        call SH(hi3,t,lb,ub,n)
      endif
    case default
      ERROR('Error drawing hilbert curve')
      call abort_program()
  end select

end subroutine SH
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Non terminal character for hilbert-snake space filling curve
!  fills in rectangularly shaped hi3 by going through hi3 in snake pattern
!  filling in cubes in 3d hilbert curve order
!-------------------------------------------------------------------------------
recursive subroutine SH_( hi3, t, lb, ub, n )

  implicit none

  integer, dimension(:,:,:), intent(inout) :: hi3
  integer, intent(inout) :: t      ! g_til_aid parameter
  integer, dimension(2), intent(inout) :: lb, ub ! keeps track of which part of hi3 were on
  integer, dimension(3), intent(in) :: n

  integer :: n_min
  integer, dimension(3) :: x

  select case (minloc(n,1))
    case (1)
      n_min = n(1)
      if ( lb(2)==1 ) then
        ! fill in window
        x = (/ 1, n_min, n_min /)
        hi3(:,lb(1):ub(1),lb(2):ub(2)) = t
        call xy_z_( hi3(:,lb(1):ub(1),lb(2):ub(2)), x, t, n_min )
        if ( t==product(n) ) return
        t = t+1
        ! move window up y
        lb(1) = lb(1) + n_min
        ub(1) = ub(1) + n_min
        ! switch to other non-terminal character
        call SH(hi3,t,lb,ub,n)
      else
        ! fill in window
        x = (/ 1, n_min, n_min /)
        hi3(:,lb(1):ub(1),lb(2):ub(2)) = t
        call xy_z_( hi3(:,lb(1):ub(1),lb(2):ub(2)), x, t, n_min )
        if ( t==product(n) ) return
        t = t+1
        ! move window down z
        lb(2) = lb(2) - n_min
        ub(2) = ub(2) - n_min
        ! call this non-terminal character again
        call SH_(hi3,t,lb,ub,n)
      endif
    case (2)
      n_min = n(2)
      if ( lb(2)==1 ) then
        ! fill in window
        x = (/ n_min, 1, n_min /)
        hi3(lb(1):ub(1),:,lb(2):ub(2)) = t
        call x_yz_( hi3(lb(1):ub(1),:,lb(2):ub(2)), x, t, n_min )
        if ( t==product(n) ) return
        t = t+1
        ! move window up x
        lb(1) = lb(1) + n_min
        ub(1) = ub(1) + n_min
        ! switch to other non-terminal character
        call SH(hi3,t,lb,ub,n)
      else
        ! fill in window
        x = (/ n_min, 1, n_min /)
        hi3(lb(1):ub(1),:,lb(2):ub(2)) = t
        call x_yz_( hi3(lb(1):ub(1),:,lb(2):ub(2)), x, t, n_min )
        if ( t==product(n) ) return
        t = t+1
        ! move window down z
        lb(2) = lb(2) - n_min
        ub(2) = ub(2) - n_min
        ! call this non-terminal character again
        call SH_(hi3,t,lb,ub,n)
      endif
    case (3)
      n_min = n(3)
      if ( lb(2)==1 ) then
        ! fill in window
        x = (/ n_min, n_min, 1 /)
        hi3(lb(1):ub(1),lb(2):ub(2),:) = t
        call y_z_x( hi3(lb(1):ub(1),lb(2):ub(2),:), x, t, n_min )
        if ( t==product(n) ) return
        t = t+1
        ! move window up x
        lb(1) = lb(1) + n_min
        ub(1) = ub(1) + n_min
        ! switch to other non-terminal character
        call SH(hi3,t,lb,ub,n)
      else
        ! fill in window
        x = (/ n_min, n_min, 1 /)
        hi3(lb(1):ub(1),lb(2):ub(2),:) = t
        call y_z_x( hi3(lb(1):ub(1),lb(2):ub(2),:), x, t, n_min )
        if ( t==product(n) ) return
        t = t+1
        ! move window down y
        lb(2) = lb(2) - n_min
        ub(2) = ub(2) - n_min
        ! call this non-terminal character again
        call SH_(hi3,t,lb,ub,n)
      endif
    case default
      ERROR('Error drawing hilbert curve')
      call abort_program()
  end select

end subroutine SH_
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! sequentialize elements of hi3 by traversing array in hilbert curve order
! Rectangular simulation box sizes are allowed but must satisfy:
!   - smallest dimension must have number of tiles equal to a power of two
!   - all other dimesions must be a multiple of the smallest dimension
!-------------------------------------------------------------------------------
subroutine fill_hi3_hilbert( hi3, n )

  implicit none

  integer, dimension(:,:,:), intent(inout) :: hi3
  integer, dimension(:), intent(in)        :: n

  integer               :: t
  integer, dimension(2) :: lb, ub
  integer, dimension(3) :: x

  if ( IAND( minval(n), minval(n)-1 )/=0 ) then
    ERROR('Number of tiles in smallest dimension must be a power of 2')
    call abort_program()
  elseif ( minval(n)==1 .and. any(n /= 1) ) then
    ERROR('If number of tiles in smallest dimension is 1,')
    ERROR('all number of tiles must be 1. Try topology="snake"')
    call abort_program()
  endif

  if ( n(1)==n(2) .and. n(2)==n(3) ) then
    if ( n(1)==1 ) then
      hi3 = 1
    else
      t   = 1
      x   = 1
      hi3 = t
      call xyz(hi3,x,t,n(1))
    endif
  elseif ( (minloc(n,1)==1) .and. (mod(n(2),n(1))==0 .and. mod(n(3),n(1))==0) ) then
      t   = 1
      hi3 = 1
      lb  = 1
      ub  = n(1)
      call SH( hi3, t, lb, ub, n )
  elseif ( (minloc(n,1)==2) .and. (mod(n(1),n(2))==0 .and. mod(n(3),n(2))==0) ) then
      t   = 1
      hi3 = 1
      lb  = 1
      ub  = n(2)
      call SH( hi3, t, lb, ub, n )
  elseif ( (minloc(n,1)==3) .and. (mod(n(1),n(3))==0 .and. mod(n(2),n(3))==0) ) then
      t   = 1
      hi3 = 1
      lb  = 1
      ub  = n(3)
      call SH( hi3, t, lb, ub, n )
  else
    ERROR('Number of tiles in other dimensions must be multiple of number')
    ERROR('of tiles in smallest dimension')
    call abort_program()
  endif

end subroutine fill_hi3_hilbert
!-------------------------------------------------------------------------------

end module m_hilbert
