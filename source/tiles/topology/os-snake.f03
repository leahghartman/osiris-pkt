#include "os-preprocess.fpp"
! #include "os-config.h"
! #include "memory/memory.h"

module m_snake

use m_term_chars
use m_system

implicit none

public :: fill_hi2_snake
public :: fill_hi3_snake

private :: S, S_
private :: yx, y_x, yx_, y_x_

contains

!-------------------------------------------------------------------------------
! non-terminal character for 2d snake grammer
!-------------------------------------------------------------------------------
recursive subroutine S( hi2, x, t, t0, n )

  integer, dimension(:,:), intent(inout) :: hi2
  integer, dimension(2), intent(inout)   :: x
  integer, intent(inout)                 :: t
  integer, intent(in)                    :: t0
  integer, dimension(2), intent(in)      :: n

  if ( t-t0+1==product(n) ) then
    t = t+1
    return
  elseif ( x(2)==n(2) ) then
    call right(hi2,x,t)
    call S_(hi2,x,t,t0,n)
  else
    call up(hi2,x,t)
    call S(hi2,x,t,t0,n)
  endif

end subroutine S
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! non-terminal character for 2d snake grammar
!-------------------------------------------------------------------------------
recursive subroutine S_( hi2, x, t, t0, n )

  integer, dimension(:,:), intent(inout) :: hi2
  integer, dimension(2), intent(inout)   :: x
  integer, intent(inout)                 :: t
  integer, intent(in)                    :: t0
  integer, dimension(2), intent(in)      :: n

  if ( t-t0+1==product(n) ) then
    t = t+1
    return
  elseif ( x(2)==1 ) then
    call right(hi2,x,t)
    call S(hi2,x,t,t0,n)
  else
    call down(hi2,x,t)
    call S_(hi2,x,t,t0,n)
  endif

end subroutine S_
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! sequentialize elements of hi3 by traversing array in snake order
!-------------------------------------------------------------------------------
subroutine fill_hi2_snake( hi2, n, t )

  implicit none

  integer, dimension(:,:), intent(inout)  :: hi2
  integer, dimension(:), intent(in)       :: n
  integer, intent(inout)                  :: t

  integer :: t0                ! starting count
  integer, dimension(2) :: x   ! starting position

  x = 1
  t0 = t
  hi2 = t
  call S(hi2,x,t,t0,n)

end subroutine fill_hi2_snake
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! non-terminal character for 3d snake grammer
!-------------------------------------------------------------------------------
recursive subroutine yx( hi3, x, t, t0, n )

	integer, dimension(:,:,:), intent(inout) :: hi3
  integer, dimension(3), intent(inout)     :: x
  integer, intent(inout)                   :: t
  integer, intent(in)                      :: t0
  integer, dimension(3), intent(in)        :: n

  if ( t-t0+1==product(n) ) then
    t = t+1
    return
  elseif ( x(2)==n(2) ) then
		if ( x(1)==n(1) ) then
			! if you reach the corner, move up and turn around
			call uz(hi3,x,t)
			call y_x_(hi3,x,t,t0,n)
		else
			! if you reach an edge, turn around
			call ux(hi3,x,t)
			call y_x(hi3,x,t,t0,n)
		endif
  else
    call uy(hi3,x,t)
    call yx(hi3,x,t,t0,n)
  endif

end subroutine yx
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! non-terminal character for 3d snake grammer
!-------------------------------------------------------------------------------
recursive subroutine y_x( hi3, x, t, t0, n )

	integer, dimension(:,:,:), intent(inout) :: hi3
  integer, dimension(3), intent(inout)     :: x
  integer, intent(inout)                   :: t
  integer, intent(in)                      :: t0
  integer, dimension(3), intent(in)        :: n

  if ( t-t0+1==product(n) ) then
    t = t+1
    return
  elseif ( x(2)==1 ) then
		if ( x(1)==n(1) ) then
			! if you reach the corner, move up and turn around
			call uz(hi3,x,t)
			call yx_(hi3,x,t,t0,n)
		else
			! if you reach an edge, turn around
			call ux(hi3,x,t)
			call yx(hi3,x,t,t0,n)
		endif
  else
    call dy(hi3,x,t)
    call y_x(hi3,x,t,t0,n)
  endif

end subroutine y_x
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! non-terminal character for 3d snake grammer
!-------------------------------------------------------------------------------
recursive subroutine yx_( hi3, x, t, t0, n )

	integer, dimension(:,:,:), intent(inout) :: hi3
  integer, dimension(3), intent(inout)     :: x
  integer, intent(inout)                   :: t
  integer, intent(in)                      :: t0
  integer, dimension(3), intent(in)        :: n

  if ( t-t0+1==product(n) ) then
    t = t+1
    return
  elseif ( x(2)==n(2) ) then
		if ( x(1)==1 ) then
			! if you reach the corner, move up and turn around
			call uz(hi3,x,t)
			call y_x(hi3,x,t,t0,n)
		else
			! if you reach an edge, turn around
			call dx(hi3,x,t)
			call y_x_(hi3,x,t,t0,n)
		endif
  else
    call uy(hi3,x,t)
    call yx_(hi3,x,t,t0,n)
  endif

end subroutine yx_
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! non-terminal character for 3d snake grammer
!-------------------------------------------------------------------------------
recursive subroutine y_x_( hi3, x, t, t0, n )

	integer, dimension(:,:,:), intent(inout) :: hi3
  integer, dimension(3), intent(inout)     :: x
  integer, intent(inout)                   :: t
  integer, intent(in)                      :: t0
  integer, dimension(3), intent(in)        :: n

  if ( t-t0+1==product(n) ) then
    t = t+1
    return
  elseif ( x(2)==1 ) then
		if ( x(1)==1 ) then
			! if you reach the corner, move up and turn around
			call uz(hi3,x,t)
			call yx(hi3,x,t,t0,n)
		else
			! if you reach an edge, turn around
			call dx(hi3,x,t)
			call yx_(hi3,x,t,t0,n)
		endif
  else
    call dy(hi3,x,t)
    call y_x_(hi3,x,t,t0,n)
  endif

end subroutine y_x_
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! sequentialize elements of hi3 by traversing array in snake order
!-------------------------------------------------------------------------------
subroutine fill_hi3_snake( hi3, n, t )

  implicit none

  integer, dimension(:,:,:), intent(inout) :: hi3
  integer, dimension(:), intent(in)        :: n
  integer, intent(inout)                   :: t

  integer :: t0                ! starting count
  integer, dimension(3) :: x

  x = 1
  t0 = t
  hi3 = t
  call yx(hi3,x,t,t0,n)

end subroutine fill_hi3_snake
!-------------------------------------------------------------------------------

end module m_snake
