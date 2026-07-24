!-----------------------------------------------------------------------------------------
! Species communication module for tiles
!-----------------------------------------------------------------------------------------

#include "os-config.h"
#include "os-preprocess.fpp"

module m_species_comm_tiles

#include "memory/memory.h"

use m_system
use m_parameters
use m_node_conf, only : t_node_conf, comm, neighbor
use m_node_conf_tiles, only : t_node_conf_tiles, p_spec_buf_block_t
use m_species_define, only : t_species, t_spec_msg, t_part_idx
use m_species_comm, only : int_pack_size, wait_spec_msg, unpack_particle_data_2
use m_species_comm, only : particle_pack_size, pack_particle_data, irecv_1_msg, isend_1_msg

#ifdef __HAS_TRACKS__
use m_species_tracks, only : missing_particles
#endif

implicit none

! extra boundary indexes for tile boundaries
integer, parameter :: p_lower_phys = 4
integer, parameter :: p_upper_phys = 8

! Modified buffer size for use in tile-to-tile comms
! Note, bnd_grp uses the p_max_buffer1_size parameter since it is presumably passing a
! large number of particles between entire nodes.
! Size in integer(p_byte)s, large enough for ~512 particles per tile
integer, parameter :: p_max_buffer1_size_t = 32768

! minimal block size to be used when growing particle buffers in a tile
! units are particle number
! this parameter is now imported from os-nconf-tiles so it can be set in the input deck
! integer, parameter :: p_spec_buf_block_t = 1024

interface set_spec_msg
  module procedure set_spec_msg
end interface

interface loc_send
  module procedure loc_send_spec
  module procedure loc_send_msg
end interface

interface isend_1_tiles
  module procedure isend_1_spec_tiles
end interface

interface irecv_1
  module procedure irecv_1_spec_group
end interface

interface irecv_1_tiles
  module procedure irecv_1_spec_tiles
end interface

interface irecv_2
  module procedure irecv_2_spec_group
end interface

interface loc_recv
  module procedure loc_recv_spec
  module procedure loc_recv_msg
end interface

interface unpack
  module procedure unpack_spec_msg_group
end interface

interface pack_group
  module procedure pack_spec_group
end interface

interface wait_send_tiles
  module procedure wait_send_tiles
end interface

public :: p_lower_phys, p_upper_phys, p_max_buffer1_size_t
public :: set_spec_msg, loc_send, isend_1_tiles, irecv_1, irecv_1_tiles, irecv_2
public :: loc_recv, unpack, pack_group, wait_send_tiles

contains

!-----------------------------------------------------------------------------------------
! initialize t_spec_msg data
!-----------------------------------------------------------------------------------------
subroutine set_spec_msg( send_spec, spec )

  implicit none

  type(t_spec_msg), intent(inout) :: send_spec
  class(t_species), intent(in) :: spec

  send_spec%size1 = 2*int_pack_size
  send_spec%size2 = 0
  send_spec%n_part2 = 0
  send_spec%n_tils1 = 0
  send_spec%particle_size = particle_pack_size( spec )

end subroutine set_spec_msg
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Prepare message buffer for intra-node tile-to-tile communication
!-----------------------------------------------------------------------------------------
subroutine loc_send_spec( spec, i, bnd, no_co, list, send_msg )

  implicit none

  class( t_species ), intent(inout) :: spec
  integer, intent(in) :: i, bnd
  class( t_node_conf ), intent(in) :: no_co

  type( t_part_idx ), intent(in) :: list
  type(t_spec_msg), dimension(2), intent(inout) :: send_msg

  integer, dimension(p_max_dim) :: shift

  send_msg( bnd )%comm = comm( no_co )
  send_msg( bnd )%node = neighbor( no_co, bnd, i )
  select type(no_co); class is(t_node_conf_tiles)
    send_msg( bnd )%tag = (no_co%ti_co%neighbor_g_til( bnd, i )-1)*2*p_x_dim+(2-bnd)+2*(i-1)
  end select
  send_msg( bnd )%particle_size = particle_pack_size( spec )

  ! if sending particles over the global simulation edge
  ! (parallel periodic) correct particle
  ! indexes
  shift = 0
  if ( no_co%on_edge( i, bnd ) ) then
    if ( bnd == p_lower ) then
      shift(i) = + spec%g_nx(i)
    else
      shift(i) = - spec%g_nx(i)
    endif
  endif

  ! Pack data and post send message
  call loc_send_msg( send_msg( bnd ), spec, list%idx, list%nidx, shift )

end subroutine loc_send_spec
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Prepare message buffer for intra-node tile-to-tile communication
!-----------------------------------------------------------------------------------------
subroutine loc_send_msg( msg, spec, par_idx, n_idx, shift )

  implicit none

  type( t_spec_msg ), intent(inout) :: msg
  class( t_species ), intent(inout) :: spec

  integer, dimension(:), pointer :: par_idx
  integer, intent(in) :: n_idx
  integer, dimension(:), intent(in) :: shift

  integer :: position, bsize, max_n_part1, cur_idx, ierr

  ! Get number of particles in each buffer
  max_n_part1 = ( msg%max_buffer1_size - int_pack_size ) / msg%particle_size

  ! sanity check
  if ( max_n_part1 < 1 ) then
    ERROR("msg%max_buffer1_size is too small")
    ERROR("please recompile with larger value")
    call abort_program()
  endif

  if ( n_idx > max_n_part1 ) then
    msg%n_part1 = max_n_part1
    msg%n_part2 = n_idx - msg%n_part1
  else
    msg%n_part1 = n_idx
    msg%n_part2 = 0
  endif

  ! Get actual message 1 size
  msg%size1 = int_pack_size + msg%n_part1 * msg%particle_size

  ! This is just a sanity check
  if ( msg%size1 > msg%max_buffer1_size ) then
    ERROR(" msg%size1 is too large, this should never happen" )
    call abort_program()
  endif

  ! Check if message buffer1 is big enough
  if ( associated( msg%buffer1 ) ) then
    bsize = size( msg%buffer1 )
  else
    bsize = 0
  endif

  if ( bsize < msg%size1 ) then
    call freemem( msg%buffer1 )
    call alloc( msg%buffer1, (/ msg%size1 /) )
  endif

  ! pack number of particles to send
  position = 0
  call mpi_pack( n_idx, 1, MPI_INTEGER, msg%buffer1, msg%size1, position, &
    mpi_comm_world, ierr )

  ! pack first set of particles if any
  if ( msg%n_part1 > 0 ) then
    cur_idx = 1
    call pack_particle_data( msg%buffer1, msg%size1, position, msg%n_part1, &
      spec, par_idx, cur_idx, shift )
  endif

  ! If necessary pack data for 2nd send message
  if ( msg%n_part2 > 0 ) then
    msg%size2 = msg%n_part2 * msg%particle_size

    ! Check if message buffer1 is big enough
    if ( associated( msg%buffer2 ) ) then
      bsize = size( msg%buffer2 )
    else
      bsize = 0
    endif
    if ( bsize < msg%size2 ) then
      call freemem( msg%buffer2 )
      call alloc( msg%buffer2, (/ msg%size2 /) )
    endif

    ! Pack particle data
    position = 0
    call pack_particle_data( msg%buffer2, msg%size2, position, msg%n_part2, &
      spec, par_idx, cur_idx, shift )

  else
    msg%size2 = 0
  endif

#ifdef __HAS_TRACKS__

! if tracking move tracks that were present and left to the missing list
  if ( ( spec%diag%ndump_fac_tracks > 0 ) .and. ( n_idx > 0 ) ) then
    call missing_particles( spec%diag%tracks, par_idx, n_idx, spec%tag )
  endif

#endif

end subroutine loc_send_msg
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Post send (step 1) for normal (boundary) species communication using local buffers
!-----------------------------------------------------------------------------------------
subroutine isend_1_spec_tiles( spec, dim, bnd, no_co, list, send_msg )

  implicit none

  class( t_species ), intent(inout) :: spec
  integer, intent(in) :: dim, bnd
  class( t_node_conf ), intent(in) :: no_co

  type( t_part_idx ), intent(in) :: list
  type(t_spec_msg), dimension(2), intent(inout) :: send_msg

  integer, dimension(p_max_dim) :: shift

  send_msg( bnd )%comm = comm( no_co )
  send_msg( bnd )%node = neighbor( no_co, bnd, dim )
  select type(no_co); class is(t_node_conf_tiles)
    send_msg( bnd )%tag  = (no_co%ti_co%neighbor_g_til( bnd, dim )-1)*&
                            2*p_x_dim+(2-bnd)+2*(dim-1)
  end select
  send_msg( bnd )%particle_size = particle_pack_size( spec )

  ! if sending particles over the global simulation edge (parallel periodic) correct
  ! particle indexes
  shift = 0
  if ( no_co%on_edge( dim, bnd ) ) then
    if ( bnd == p_lower ) then
      shift(dim) = + spec%g_nx(dim)
    else
      shift(dim) = - spec%g_nx(dim)
    endif
  endif

  ! Pack data and post send message
  call isend_1_msg( send_msg( bnd ), spec, list%idx, list%nidx, shift )

end subroutine isend_1_spec_tiles
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Post receive (step 1) for group (boundary) species communication using local buffers
!-----------------------------------------------------------------------------------------
subroutine irecv_1_spec_group( spec, node, no_co, recv_msg )

  implicit none

  class( t_species ), intent(in) :: spec
  integer, intent(in) :: node
  class( t_node_conf ), intent(in) :: no_co
  type(t_spec_msg), intent(inout) :: recv_msg

  recv_msg%comm = comm( no_co )
  recv_msg%node = node
  recv_msg%tag  = 0
  recv_msg%particle_size = particle_pack_size( spec )

  call irecv_1_msg( recv_msg )

end subroutine irecv_1_spec_group
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Post receive (step 1) for normal (boundary) species communication using local buffers
!-----------------------------------------------------------------------------------------
subroutine irecv_1_spec_tiles( spec, i, bnd, no_co, recv_msg )

  implicit none

  class( t_species ), intent(in) :: spec
  integer, intent(in) :: i, bnd
  class( t_node_conf ), intent(in) :: no_co
  type(t_spec_msg), dimension(2), intent(inout) :: recv_msg

  recv_msg( bnd )%comm = comm( no_co )
  recv_msg( bnd )%node = neighbor( no_co, bnd, i )
  select type(no_co); class is(t_node_conf_tiles)
    recv_msg( bnd )%tag  = (no_co%ti_co%g_til_aid-1)*2*p_x_dim+(bnd-1)+2*(i-1)
  end select
  recv_msg( bnd )%particle_size = particle_pack_size( spec )

  call irecv_1_msg( recv_msg( bnd ) )


end subroutine irecv_1_spec_tiles
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Start 2nd part of grouped receive message communication
!-----------------------------------------------------------------------------------------
subroutine irecv_2_spec_group( recv_spec, garbage )

  implicit none

  type( t_spec_msg ), intent(inout) :: recv_spec
  ! unused dummy argument for unique interface
  integer, intent(in) :: garbage

  integer :: bsize, position, ierr

  ! If 1st message block is still active wait for it to complete
  call wait_spec_msg( recv_spec )

  ! Get size of 2nd message
  position = 0
  call mpi_unpack( recv_spec%buffer1, recv_spec%size1, position, recv_spec%size2, 1, &
    MPI_INTEGER, mpi_comm_world, ierr )
  if (ierr/=0) then
   ERROR("MPI error")
   call abort_program( p_err_mpi )
  endif

  ! If 2nd message is required then
  if ( recv_spec%size2 > 0 ) then

    ! check if buffer is big enough
    if ( associated( recv_spec%buffer2 ) ) then
      bsize = size( recv_spec%buffer2 )
    else
      bsize = 0
    endif
    if ( recv_spec%size2 > bsize ) then
      call freemem( recv_spec%buffer2 )
      call alloc( recv_spec%buffer2, (/ recv_spec%size2 /) )
    endif

    ! post receive for second msg
    call mpi_irecv( recv_spec%buffer2, recv_spec%size2, MPI_PACKED, recv_spec%node - 1, &
      recv_spec%tag, recv_spec%comm, recv_spec%request, ierr )
    if (ierr/=0) then
       ERROR("MPI error")
       call abort_program( p_err_mpi )
    endif

  else
    recv_spec%request = MPI_REQUEST_NULL
  endif

end subroutine irecv_2_spec_group
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Grab message buffer from neighboring tile (same node) and unpack the data
!-----------------------------------------------------------------------------------------
subroutine loc_recv_spec( spec, i, bnd, no_co, recv_msg, bnd_cross )

  implicit none

  class( t_species ), intent(inout) :: spec
  integer, intent(in) :: i, bnd
  class( t_node_conf ), intent(in) :: no_co
  type(t_spec_msg), dimension(2), intent(inout) :: recv_msg
  type( t_part_idx ), dimension(:), intent(inout) :: bnd_cross

  recv_msg( bnd )%comm = comm( no_co )
  recv_msg( bnd )%node = neighbor( no_co, bnd, i )
  select type(no_co); class is(t_node_conf_tiles)
    recv_msg( bnd )%tag  = (no_co%ti_co%g_til_aid-1)*2*p_x_dim+(bnd-1)+2*(i-1)
  end select
  recv_msg( bnd )%particle_size = particle_pack_size( spec )

  call loc_recv_msg( recv_msg( bnd ) )

  select case (bnd)
    case ( p_lower )
       call loc_unpack_particles( recv_msg( p_lower ), spec, bnd_cross )
    case ( p_upper )
       call loc_unpack_particles( recv_msg( p_upper ), spec, bnd_cross )
  end select

end subroutine loc_recv_spec
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Begin to process msg buffer from a tile in the same node
!-----------------------------------------------------------------------------------------
subroutine loc_recv_msg( msg )

  implicit none

  type( t_spec_msg ), intent(inout) :: msg

  integer :: position, npart, ierr

  ! Get total number of particles to receive
  position = 0
  call mpi_unpack( msg%buffer1, msg%size1, position, npart, 1, MPI_INTEGER, mpi_comm_world, ierr )
  if (ierr/=0) then
   ERROR("MPI error")
   call abort_program( p_err_mpi )
  endif

  if ( int_pack_size + npart * msg%particle_size <= msg%max_buffer1_size ) then
    ! 1 message was enough
    msg%n_part1 = npart
    msg%n_part2 = 0
  else
    ! A 2nd message is required
    msg%n_part1 = ( msg%max_buffer1_size - int_pack_size ) / msg%particle_size
    msg%n_part2 = npart - msg%n_part1
  endif

  ! If 2nd message is required then
  if ( .not. ( msg%n_part2>0 ) ) then
    msg%request = MPI_REQUEST_NULL
    msg%size2   = 0
  endif

end subroutine loc_recv_msg
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! unpack single message from group
!-----------------------------------------------------------------------------------------
subroutine unpack_spec_msg_group( spec, recv_buff, recv_size, position, bnd_cross )

  implicit none

  class(t_species), intent(inout)                 :: spec
  integer(p_byte), dimension(:), intent(in)       :: recv_buff
  integer, intent(in)                             :: recv_size
  integer, intent(inout)                          :: position
  type( t_part_idx ), dimension(:), intent(inout) :: bnd_cross

  integer :: npart, new_particles, ierr

  ! unpack number of particles this tile will be receiving
  call mpi_unpack( recv_buff, recv_size, position, npart, 1, &
    MPI_INTEGER, mpi_comm_world, ierr )

  ! Check if it is necessary to grow particle buffer
  ! Number of new particles will be particles being received minus the number of
  ! holes left in the particle buffer

  new_particles = npart - ((bnd_cross(1)%nidx - bnd_cross(1)%start + 1 ) + &
    (bnd_cross(2)%nidx - bnd_cross(2)%start + 1 ))

  if ( spec%num_par + new_particles > spec%num_par_max ) then
    ! Decided against growing by 10%, now grow by specific amount
    call spec%grow_buffer( spec%num_par + new_particles + p_spec_buf_block_t )
    ! call spec%grow_buffer( ceiling( 1.1*( spec%num_par+new_particles ) ) )
  endif

  ! Unpack particle data
  call unpack_particle_data_2( recv_buff, recv_size, position, &
    npart, spec, bnd_cross )

end subroutine unpack_spec_msg_group
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Unpacks message buffer into species object for 2 hole lists
!-----------------------------------------------------------------------------------------
subroutine loc_unpack_particles( msg, species, bnd_cross )

  implicit none

  type( t_spec_msg ), intent(inout) :: msg
  class( t_species ), intent(inout) :: species
  type( t_part_idx ), dimension(:), intent(inout) :: bnd_cross


  integer :: npart, new_particles, position, ierr

  ! Number of new particles will be particles being received minus the number of
  ! holes left in the particle buffer
  new_particles = (msg%n_part1 + msg%n_part2) - &
                    ((bnd_cross(1)%nidx - bnd_cross(1)%start + 1 ) + &
                      (bnd_cross(2)%nidx - bnd_cross(2)%start + 1 ))

  if ( species%num_par + new_particles > species%num_par_max ) then
    ! Decided against growing by 10%, now grow by specific amount
    call species%grow_buffer( species%num_par + new_particles + p_spec_buf_block_t )
    ! call species%grow_buffer( ceiling( 1.1*( species%num_par+new_particles ) ) )
  endif

  ! unpack buffer1
  if ( msg%n_part1 > 0 ) then
    ! jump over total number of particles in buffer1
    position = 0
    call mpi_unpack( msg%buffer1, msg%size1, position, &
      npart, 1, MPI_INTEGER,  &
      mpi_comm_world, ierr )

    ! Unpack particle data
    call unpack_particle_data_2( msg%buffer1, msg%size1, position, msg%n_part1, &
      species, bnd_cross )
  endif

  ! unpack buffer2
  if ( msg%n_part2 > 0 ) then
    position = 0
    call unpack_particle_data_2( msg%buffer2, msg%size2, position, msg%n_part2, &
      species, bnd_cross )
  endif

end subroutine loc_unpack_particles
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Pack send buffer for grouped particle update boundary
!-----------------------------------------------------------------------------------------
subroutine pack_spec_group( spec, no_co, send_buff, send_size, bnd_cross, til_recv, &
                            npart, position, bnd, dim )

  implicit none

  class(t_species), intent(inout)              :: spec
  class(t_node_conf), intent(in)               :: no_co
  integer(p_byte), dimension(:), intent(inout) :: send_buff
  integer, intent(in)                          :: send_size
  type( t_part_idx ), dimension(2)             :: bnd_cross
  integer, intent(in)                          :: til_recv
  integer, intent(in)                          :: npart
  integer, intent(inout)                       :: position
  integer, intent(in)                          :: bnd
  integer, intent(in)                          :: dim

  integer, dimension(2) :: preface
  integer, dimension(p_max_dim) :: shift
  integer :: cur_idx, ierr

  ! Pack preface of til's msg
  preface = (/til_recv, npart/)
  call mpi_pack( preface, 2, MPI_INTEGER, send_buff, send_size, position, &
    mpi_comm_world, ierr )

  ! if sending particles over the global simulation edge (parallel periodic) correct
  ! particle indexes
  shift = 0
  if ( no_co%on_edge( dim, bnd ) ) then
    if ( bnd == p_lower ) then
      shift(dim) = + spec%g_nx(dim)
    else
      shift(dim) = - spec%g_nx(dim)
    endif
  endif

  ! Pack the particles
  cur_idx = 1
  call pack_particle_data( send_buff, send_size, position, npart, spec, &
    bnd_cross(bnd)%idx, cur_idx, shift )

end subroutine pack_spec_group
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!
!-----------------------------------------------------------------------------------------
subroutine wait_send_tiles( bnd, send_msg )

  implicit none

  integer, intent(in) :: bnd
  type(t_spec_msg), dimension(2), intent(inout) :: send_msg

  select case (bnd)

    case ( p_lower )
       call wait_spec_msg( send_msg( p_lower ) )

    case ( p_upper )
       call wait_spec_msg( send_msg( p_upper ) )

    case ( p_lower+p_upper )
       call wait_spec_msg( send_msg(p_lower), send_msg(p_upper) )

    case ( p_lower_phys+p_upper )
       call wait_spec_msg( send_msg( p_upper ) )

    case ( p_lower+p_upper_phys )
       call wait_spec_msg( send_msg( p_lower ) )

  end select


end subroutine wait_send_tiles
!-----------------------------------------------------------------------------------------

end module m_species_comm_tiles