#include "os-config.h"
#include "os-preprocess.fpp"

module m_wall_comm_tiles

#include "memory/memory.h"

use m_parameters

use m_wall_define, only : t_wall, lbound, ubound
use m_vdf_comm, only : t_vdf_msg, irecv, isend
use m_wall_comm, only : pack_wall_msg, unpack_wall_msg
use m_node_conf, only : t_node_conf, neighbor
use m_node_conf_tiles, only : t_node_conf_tiles

private

public :: loc_recv_wall, loc_send_wall
public :: irecv_wall_tiles, isend_wall_tiles

contains

!-----------------------------------------------------------------------------------------
! Pack send buffer for local boundary communication
!-----------------------------------------------------------------------------------------
subroutine loc_send_wall( wall, update_type, dim, bnd, nx_move, send_msg )

  implicit none

  type( t_wall ), intent(in) :: wall
  integer, intent(in) :: update_type, dim, bnd

  integer, intent(in) :: nx_move
  type(t_vdf_msg), dimension(2), intent(inout) :: send_msg

  integer, dimension(2,p_max_dim) :: ns
  integer :: i, bsize, msg_size

  ! Get range of cells to send
  msg_size = wall%f_dim
  do i = 1, wall%x_dim
    if ( i == dim ) then
       if ( update_type == p_vdf_add ) then
         select case ( bnd )
           case ( p_upper )
             ns( p_lower, i ) = wall%nx(i) + 1 - wall%gc_num( p_lower, i )
             ns( p_upper, i ) = wall%nx(i) + wall%gc_num( p_upper, i )
           case ( p_lower )
             ns( p_lower, i ) = 1 - wall%gc_num( p_lower, i )
             ns( p_upper, i ) = wall%gc_num( p_upper, i )
         end select
       else
         select case ( bnd )
           case ( p_upper )
             ns( p_lower, i ) = wall%nx(i) + 1 - wall%gc_num( p_lower, i )
             ns( p_upper, i ) = wall%nx(i) - nx_move
           case ( p_lower )
             ns( p_lower, i ) = 1 - nx_move
             ns( p_upper, i ) = wall%gc_num( p_upper, i )
         end select
       endif
    else
      ! full cell range along other directions
      ns( p_lower, i ) = lbound( wall, i+1 )
      ns( p_upper, i ) = ubound( wall, i+1 )
    endif
    msg_size = msg_size * (ns( p_upper, i ) - ns( p_lower,i ) + 1 )
  enddo

  ! Only set member data needed for packing (none needed for sending)
  send_msg(bnd)%range( :, 1 : wall%x_dim ) = ns ( :, 1 : wall%x_dim )

  ! grow message buffer if necessary
  if ( associated( send_msg(bnd)%buffer ) ) then
    bsize = size( send_msg(bnd)%buffer )
  else
    bsize = 0
  endif

  if ( msg_size > bsize ) then
    call freemem( send_msg(bnd)%buffer )
    call alloc( send_msg(bnd)%buffer, (/msg_size/) )
  endif

  ! pack message data
  call pack_wall_msg( send_msg(bnd), wall )

end subroutine loc_send_wall
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! post send message for boundary communication
!-----------------------------------------------------------------------------------------
subroutine isend_wall_tiles( wall, update_type, dim, bnd, no_co, nx_move, send_msg )

  implicit none

  type( t_wall ), intent(in) :: wall
  integer, intent(in) :: update_type, dim, bnd
  class( t_node_conf ), intent(in) :: no_co

  integer, intent(in) :: nx_move
  type(t_vdf_msg), dimension(2), intent(inout) :: send_msg

  integer, dimension(2,p_max_dim) :: ns
  integer :: i, bsize, msg_size

  ! sanity check
  if ( send_msg(bnd)%handle /= MPI_REQUEST_NULL ) then
    ERROR( 'isend_wall_tiles called when message still waiting to complete' )
    call abort_program( p_err_invalid )
  endif

  ! Get range of cells to send
  msg_size = wall%f_dim
  do i = 1, wall%x_dim
    if ( i == dim ) then
       if ( update_type == p_vdf_add ) then
         select case ( bnd )
           case ( p_upper )
             ns( p_lower, i ) = wall%nx(i) + 1 - wall%gc_num( p_lower, i )
             ns( p_upper, i ) = wall%nx(i) + wall%gc_num( p_upper, i )
           case ( p_lower )
             ns( p_lower, i ) = 1 - wall%gc_num( p_lower, i )
             ns( p_upper, i ) = wall%gc_num( p_upper, i )
         end select
       else
         select case ( bnd )
           case ( p_upper )
             ns( p_lower, i ) = wall%nx(i) + 1 - wall%gc_num( p_lower, i )
             ns( p_upper, i ) = wall%nx(i) - nx_move
           case ( p_lower )
             ns( p_lower, i ) = 1 - nx_move
             ns( p_upper, i ) = wall%gc_num( p_upper, i )
         end select
       endif
    else
      ! full cell range along other directions
      ns( p_lower, i ) = lbound( wall, i+1 )
      ns( p_upper, i ) = ubound( wall, i+1 )
    endif
    msg_size = msg_size * (ns( p_upper, i ) - ns( p_lower,i ) + 1 )
  enddo

  ! Use vdf-comm send buffers
  send_msg(bnd)%type                       = update_type
  send_msg(bnd)%range( :, 1 : wall%x_dim ) = ns ( :, 1 : wall%x_dim )
  send_msg(bnd)%msg_size                   = msg_size
  send_msg(bnd)%node                       = neighbor( no_co, bnd, dim )
  select type(no_co); class is(t_node_conf_tiles)
    send_msg(bnd)%tag = (no_co%ti_co%neighbor_g_til( bnd, dim )-1)*2*wall%x_dim+(2-bnd)+2*(dim-1)
  end select

  ! grow message buffer if necessary
  if ( associated( send_msg(bnd)%buffer ) ) then
    bsize = size( send_msg(bnd)%buffer )
  else
    bsize = 0
  endif

  if ( msg_size > bsize ) then
    call freemem( send_msg(bnd)%buffer )
    call alloc( send_msg(bnd)%buffer, (/msg_size/) )
  endif

  ! pack message data
  call pack_wall_msg( send_msg(bnd), wall )

  ! post send using vdf-comm isend
  call isend( send_msg(bnd), no_co )

end subroutine isend_wall_tiles
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! post recv message for boundary communication
!-----------------------------------------------------------------------------------------
subroutine loc_recv_wall( wall, update_type, dim, bnd, nx_move, recv_msg )

  implicit none

  type(t_wall), intent(inout) :: wall
  integer, intent(in) :: update_type, dim, bnd
  integer, intent(in) :: nx_move
  type(t_vdf_msg), dimension(2), intent(inout) :: recv_msg

  integer, dimension(2,p_max_dim) :: nr
  integer :: i, msg_size

  ! Get range of cells to receive
  msg_size = wall%f_dim
  do i = 1, wall%x_dim
    if ( i == dim ) then
       if ( update_type == p_vdf_add ) then
         select case ( bnd )
           case ( p_upper )
             nr( p_lower, i ) = wall%nx(i) + 1 - wall%gc_num( p_lower, i )
             nr( p_upper, i ) = wall%nx(i) + wall%gc_num( p_upper, i )
           case ( p_lower )
             nr( p_lower, i ) = 1 - wall%gc_num( p_lower, i )
             nr( p_upper, i ) = wall%gc_num( p_upper, i )
         end select
       else
         select case ( bnd )
           case ( p_upper )
             nr( p_lower, i ) = wall%nx(i) + 1 - nx_move
             nr( p_upper, i ) = wall%nx(i) + wall%gc_num( p_upper, i )
           case ( p_lower )
             nr( p_lower, i ) = 1 - wall%gc_num( p_lower, i )
             nr( p_upper, i ) = 0 - nx_move
         end select
       endif
    else
      ! full cell range along other directions
      nr( p_lower, i ) = lbound( wall, i+1 )
      nr( p_upper, i ) = ubound( wall, i+1 )
    endif
    msg_size = msg_size * (nr( p_upper, i ) - nr( p_lower,i ) + 1 )
  enddo

  ! (*debug*)
  if ( msg_size == 0 ) then
    ERROR('Invalid message size, aborting.')
    call abort_program( p_err_invalid )
  endif

  ! Fill in only msg member data needed to unpack (none needed for receiving)
  recv_msg(bnd)%type                       = update_type
  recv_msg(bnd)%range( :, 1 : wall%x_dim ) = nr ( :, 1 : wall%x_dim )

  ! unpack message data
  call unpack_wall_msg( recv_msg(bnd), wall )

end subroutine loc_recv_wall
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! post recv message for boundary communication
!-----------------------------------------------------------------------------------------
subroutine irecv_wall_tiles( wall, update_type, dim, bnd, no_co, nx_move, recv_msg )

  implicit none

  type(t_wall), intent(inout) :: wall
  integer, intent(in) :: update_type, dim, bnd
  class( t_node_conf ), intent(in) :: no_co
  integer, intent(in) :: nx_move
  type(t_vdf_msg), dimension(2), intent(inout) :: recv_msg

  integer, dimension(2,p_max_dim) :: nr
  integer :: i, bsize, msg_size

  ! Check if message is still waiting to complete
  if ( recv_msg( bnd )%handle /= MPI_REQUEST_NULL ) then
    ERROR( 'irecv_wall_tiles called when message still waiting to complete' )
    call abort_program( p_err_invalid )
  endif

  ! Get range of cells to receive
  msg_size = wall%f_dim
  do i = 1, wall%x_dim
    if ( i == dim ) then
       if ( update_type == p_vdf_add ) then
         select case ( bnd )
           case ( p_upper )
             nr( p_lower, i ) = wall%nx(i) + 1 - wall%gc_num( p_lower, i )
             nr( p_upper, i ) = wall%nx(i) + wall%gc_num( p_upper, i )
           case ( p_lower )
             nr( p_lower, i ) = 1 - wall%gc_num( p_lower, i )
             nr( p_upper, i ) = wall%gc_num( p_upper, i )
         end select
       else
         select case ( bnd )
           case ( p_upper )
             nr( p_lower, i ) = wall%nx(i) + 1 - nx_move
             nr( p_upper, i ) = wall%nx(i) + wall%gc_num( p_upper, i )
           case ( p_lower )
             nr( p_lower, i ) = 1 - wall%gc_num( p_lower, i )
             nr( p_upper, i ) = 0 - nx_move
         end select
       endif
    else
      ! full cell range along other directions
      nr( p_lower, i ) = lbound( wall, i+1 )
      nr( p_upper, i ) = ubound( wall, i+1 )
    endif
    msg_size = msg_size * (nr( p_upper, i ) - nr( p_lower,i ) + 1 )
  enddo

  ! (*debug*)
  if ( msg_size == 0 ) then
    ERROR('Invalid message size, aborting.')
    call abort_program( p_err_invalid )
  endif

  ! Use vdf-comm recv buffers
  recv_msg(bnd)%type                       = update_type
  recv_msg(bnd)%range( :, 1 : wall%x_dim ) = nr ( :, 1 : wall%x_dim )
  recv_msg(bnd)%msg_size                   = msg_size
  recv_msg(bnd)%node                       = neighbor( no_co, bnd, dim )
  select type(no_co); class is(t_node_conf_tiles)
    recv_msg(bnd)%tag = (no_co%ti_co%g_til_aid-1)*2*wall%x_dim+(bnd-1)+2*(dim-1)
  end select

  ! grow message buffer if necessary
  if ( associated( recv_msg(bnd)%buffer ) ) then
    bsize = size( recv_msg(bnd)%buffer )
  else
    bsize = 0
  endif

  if ( msg_size > bsize ) then
    call freemem( recv_msg(bnd)%buffer )
    call alloc( recv_msg(bnd)%buffer, (/msg_size/) )
  endif

  ! post receive
  call irecv( recv_msg(bnd), no_co )

end subroutine irecv_wall_tiles
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Update boundary values of a wall on communication and parallel boundaries
! First pass when all communications were over mpi
!-----------------------------------------------------------------------------------------
! subroutine update_boundary_tiles_wall( wall, update_type, no_co, move_num, send_msg, recv_msg )

!   implicit none

!   type( t_wall ), intent(inout) :: wall
!   integer,                 intent(in) :: update_type
!   class( t_node_conf ),             intent(in) :: no_co
!   integer, dimension(:), intent(in) :: move_num
!   type(t_vdf_msg), dimension(2), intent(inout) :: send_msg, recv_msg

!   integer :: i, my_node, lneighbor, uneighbor

!   my_node = no_co % my_aid()

!   do i = 1, wall%x_dim

!     ! skip the direction perpendicular to the wall
!     if (i /= wall%idir) then

!       if ( my_node == neighbor(no_co,p_lower,i) ) then

!         !single node periodic
!         call update_periodic_wall( wall, i, update_type )

!       else

!         ! communication with other nodes
!         lneighbor = neighbor( no_co, p_lower, i )
!         uneighbor = neighbor( no_co, p_upper, i )

!         ! post receives
!         if ( lneighbor > 0 ) call irecv_wall_tiles( wall, update_type, i, p_lower, no_co, move_num(i), recv_msg )
!         if ( uneighbor > 0 ) call irecv_wall_tiles( wall, update_type, i, p_upper, no_co, move_num(i), recv_msg )

!         ! post sends
!         ! Note: sends MUST be posted in the opposite order of receives to account for a 2 node
!         ! periodic partition, where 1 node sends 2 messages to the same node
!         if ( uneighbor > 0 ) call isend_wall_tiles( wall, update_type, i, p_upper, no_co, move_num(i), send_msg )
!         if ( lneighbor > 0 ) call isend_wall_tiles( wall, update_type, i, p_lower, no_co, move_num(i), send_msg )

!         ! wait for receives and unpack data
!         if ( lneighbor > 0 ) call wait_irecv_wall( wall, p_lower, recv_msg )
!         if ( uneighbor > 0 ) call wait_irecv_wall( wall, p_upper, recv_msg )

!         ! wait for sends to complete
!         if ( uneighbor > 0 ) call wait_isend_wall( p_upper, send_msg )
!         if ( lneighbor > 0 ) call wait_isend_wall( p_lower, send_msg )

!       endif

!     endif

!   enddo

! end subroutine update_boundary_tiles_wall
!-----------------------------------------------------------------------------------------

end module m_wall_comm_tiles