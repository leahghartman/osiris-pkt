# 1 "tiles/os-vdf-comm-tiles.f90"
# 1 "<built-in>" 1
# 1 "<built-in>" 3
# 467 "<built-in>" 3
# 1 "<command line>" 1
# 1 "<built-in>" 2
# 1 "tiles/os-vdf-comm-tiles.f90" 2
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
# 2 "tiles/os-vdf-comm-tiles.f90" 2
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
# 3 "tiles/os-vdf-comm-tiles.f90" 2

module m_vdf_comm_tiles

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
# 7 "tiles/os-vdf-comm-tiles.f90" 2

use m_parameters
use m_vdf_define, only : t_vdf, t_vdf_report, t_vdf_report_item
use m_node_conf, only : t_node_conf, comm, neighbor
use m_node_conf_tiles, only : t_node_conf_tiles
use m_grid_define, only : t_grid
use m_space, only : t_space
use m_vdf_comm, only : t_vdf_msg, pack_vdf_msg, unpack_vdf_msg
use m_time_step, only : t_time_step

implicit none

private

! Parameter for passing all the messages
integer, parameter :: p_msg = 1
! Parameter for calling mpi wait
integer, parameter :: p_wait = 2

integer, parameter :: p_send_msg = 1
integer, parameter :: p_recv_msg = 2

type :: t_bnd_msg

  ! node to send/receive messages from
  integer :: node = -1

  ! number of tiles that communicate with this node
  integer :: n_tils = -1

  ! integer to add up message size
  ! dim 1 corresponds to send or recv
  ! dim 2 corresponds to tiles, with one extra for total
  ! see m_group_boundary for parameter definitions
  integer, dimension(:,:), pointer :: msg_size => null()

  ! keep track of current position in message.
  ! This may go in send_vdf/recv_vdf, or as a local variable in the routine
  integer :: vdf_pos = 0

  ! - Comments pertaining to t_bnd_patt%msg:
  ! This array contains information about tiles communicating with other node
  ! size(3,n_communications_with_node)
  ! dim 1 contains data ( my_til, upper/lower, other_til )
  ! - Comments pertaining to t_bnd_patt%msg_loc:
  ! This array contains info about tiles communicating within a node
  ! size(3,n_local_communications)
  ! dim 1 contains data ( smaller g_til_aid in comm, upper or lower from perspective
  ! of smaller g_til_aid tile, larger g_til_aid in comm )
  integer, dimension(:,:), pointer :: tils => null()

  ! cell range involved in communication
  ! First index is lower/upper, second is dim, third is tiles
  integer, dimension(:,:,:), pointer :: send_vdf_range => null(), recv_vdf_range => null()

  contains

  procedure :: map_tile_to_idx
  procedure :: total_size
  procedure :: cleanup => cleanup_bnd_msg

end type t_bnd_msg

interface isend_size
  module procedure isend_size_vdf
end interface

interface loc_send_vdf
  module procedure loc_send_vdf
end interface

interface irecv_size
  module procedure irecv_size_vdf
end interface

interface loc_recv_vdf
  module procedure loc_recv_vdf
end interface

interface fld_pack_size
  module procedure fld_pack_size
end interface

interface unpack_group
  module procedure unpack_vdf_msg_group
end interface

interface pack_group
  module procedure pack_vdf_msg_group
end interface

interface isend_tiles
  module procedure isend_vdf_tiles
  module procedure isend_vdf_msg_tiles
end interface

interface irecv_tiles
  module procedure irecv_vdf_group
  module procedure irecv_vdf_tiles
  module procedure irecv_vdf_msg_tiles
end interface

interface get_send_range
  module procedure get_send_range
end interface

interface get_recv_range
  module procedure get_recv_range
end interface

interface pack_fld_msg_dlb_group
  module procedure pack_fld_msg_dlb_group
end interface

interface unpack_fld_msg_dlb_group
  module procedure unpack_fld_msg_dlb_group
end interface

public :: p_msg, p_wait, p_send_msg, p_recv_msg
public :: t_bnd_msg
public :: isend_size, irecv_size, loc_send_vdf, loc_recv_vdf
public :: irecv_tiles, isend_tiles
public :: fld_pack_size, unpack_group, pack_group
public :: get_send_range, get_recv_range
public :: pack_fld_msg_dlb_group, unpack_fld_msg_dlb_group

contains

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
function map_tile_to_idx( this, tile, bnd )

  class(t_bnd_msg), intent(in) :: this
  integer, intent(in) :: tile, bnd

  integer :: map_tile_to_idx
  integer :: i

  map_tile_to_idx = 0

  do i = 1, this%n_tils
    if (this%tils(1,i) == tile .and. this%tils(2,i) == bnd) then
      map_tile_to_idx = i
      exit
    endif
  enddo

  ! Sanity check
  if (map_tile_to_idx == 0) then
    write(err_buf__,*) 'Tile not found in map_tile_to_idx';call err__("tiles/os-vdf-comm-tiles.f90",156)
  endif

end function map_tile_to_idx
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine total_size( this )

  class(t_bnd_msg), intent(in) :: this

  this%msg_size(p_send_msg,this%n_tils+1) = &
    sum(this%msg_size(p_send_msg,1:this%n_tils))

  this%msg_size(p_recv_msg,this%n_tils+1) = &
    sum(this%msg_size(p_recv_msg,1:this%n_tils))

end subroutine total_size
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine cleanup_bnd_msg( this )

  class(t_bnd_msg), intent(in) :: this

  if ( associated(this%msg_size) ) call freemem(this%msg_size,"tiles/os-vdf-comm-tiles.f90",183)
  if ( associated(this%tils) ) call freemem(this%tils,"tiles/os-vdf-comm-tiles.f90",184)
  if ( associated(this%send_vdf_range) ) call freemem(this%send_vdf_range,"tiles/os-vdf-comm-tiles.f90",185)
  if ( associated(this%recv_vdf_range) ) call freemem(this%recv_vdf_range,"tiles/os-vdf-comm-tiles.f90",186)

end subroutine cleanup_bnd_msg
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! calculate grouped send message size for boundary communication
!-----------------------------------------------------------------------------------------
subroutine isend_size_vdf( vdf, update_type, dim, bnd, nx_move, msg, til_idx, vdf_b )

  implicit none

  type( t_vdf ), intent(in) :: vdf
  type( t_vdf ), intent(in), optional :: vdf_b
  integer, intent(in) :: update_type, dim, bnd
  integer, intent(in) :: nx_move, til_idx
  class( t_bnd_msg ), intent(inout) :: msg

  integer, dimension(2,p_max_dim) :: ns
  integer :: msg_size

  ! Get range of cells to send
  call get_send_range( vdf, msg_size, ns, bnd, update_type, dim, nx_move )

  ! If second vdf is present double message size
  if ( present( vdf_b ) ) msg_size = msg_size * 2

  ! Add 1 for tile number and 1 for upper/lower
  msg_size = msg_size + 2

  ! Store ns for later use
  msg%send_vdf_range( :, 1 : vdf%x_dim_, til_idx ) = ns ( :, 1 : vdf%x_dim_ )

  ! Store message size
  msg%msg_size(p_send_msg,til_idx) = msg_size

end subroutine isend_size_vdf
!-----------------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! prepare send message for boundary communication between tiles on same node
!-------------------------------------------------------------------------------
subroutine loc_send_vdf( vdf, update_type, dim, bnd, no_co, nx_move, send_msg, &
                         vdf_b )

  implicit none

  type( t_vdf ), intent(in) :: vdf
  type( t_vdf ), intent(in), optional :: vdf_b

  integer, intent(in) :: update_type, dim, bnd
  class( t_node_conf ), intent(in) :: no_co

  integer, intent(in) :: nx_move
  type(t_vdf_msg), dimension(2), intent(inout) :: send_msg

  integer, dimension(2,p_max_dim) :: ns
  integer :: i, bsize, msg_size


  if ( neighbor( no_co, bnd, dim ) < 0 ) then
    write(0,*) 'irecv_vdf called without a valid neighbor, dim, bnd = ',dim, bnd, &
               trim("tiles/os-vdf-comm-tiles.f90"), ':', 248
    call abort_program( p_err_invalid )
  endif

  ! Get range of cells to send
  call get_send_range( vdf, msg_size, ns, bnd, update_type, dim, nx_move )

  ! If second vdf is present double message size
  if ( present( vdf_b ) ) msg_size = msg_size * 2

  send_msg(bnd)%type = update_type ! unnecessary
  send_msg(bnd)%range( :, 1:vdf%x_dim_ ) = ns ( :, 1 : vdf%x_dim_ )
  send_msg(bnd)%msg_size = msg_size
  send_msg(bnd)%node = neighbor( no_co, bnd, dim )
  select type(no_co); class is(t_node_conf_tiles)
    send_msg(bnd)%tag = (no_co%ti_co%neighbor_g_til( bnd, dim )-1)*2*vdf%x_dim_+(2-bnd)+2*(dim-1)
  end select

  ! grow message buffer if necessary
  if ( associated( send_msg(bnd)%buffer ) ) then
    bsize = size( send_msg(bnd)%buffer )
  else
    bsize = 0
  endif

  if ( msg_size > bsize ) then
    call freemem(send_msg(bnd)%buffer,"tiles/os-vdf-comm-tiles.f90",274)
    call alloc(send_msg(bnd)%buffer, (/msg_size/),"tiles/os-vdf-comm-tiles.f90",275)
  endif

  ! pack message data
  i = 1
  call pack_vdf_msg( send_msg(bnd), vdf, i )
  if ( present( vdf_b ) ) call pack_vdf_msg( send_msg(bnd), vdf_b, i )

  ! post send
  ! call isend_vdf_msg( send_msg(bnd), no_co )

end subroutine loc_send_vdf
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! post grouped recv message for boundary communication
!-----------------------------------------------------------------------------------------
subroutine irecv_vdf_group( update_type, node, no_co, recv_msg )

  implicit none

  integer, intent(in) :: update_type, node
  class( t_node_conf ), intent(in) :: no_co
  type(t_vdf_msg), intent(inout) :: recv_msg

  ! Set recv_msg parameters
  recv_msg%type = update_type
  ! msg_size set in check_buffer_size_vdf
  recv_msg%node = node
  recv_msg%tag = 0

  ! post receive
  call irecv_vdf_msg_tiles( recv_msg, no_co )

end subroutine irecv_vdf_group
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! calculate grouped recv message size for boundary communication
!-----------------------------------------------------------------------------------------
subroutine irecv_size_vdf( vdf, update_type, dim, bnd, nx_move, msg, til_idx, vdf_b )

  implicit none

  type( t_vdf ), intent(in) :: vdf
  type( t_vdf ), intent(in), optional :: vdf_b
  integer, intent(in) :: update_type, dim, bnd
  integer, intent(in) :: nx_move, til_idx
  class( t_bnd_msg ), intent(inout) :: msg

  integer, dimension(2,p_max_dim) :: nr
  integer :: msg_size

  ! Get range of cells to receive
  call get_recv_range( vdf, msg_size, nr, bnd, update_type, dim, nx_move )

  ! (* debug *)
  if ( msg_size == 0 ) then
    write(err_buf__,*) 'Invalid message size, aborting.';call err__("tiles/os-vdf-comm-tiles.f90",333)
    call abort_program( p_err_invalid )
  endif

  ! If second vdf is present double message size
  if ( present( vdf_b ) ) msg_size = msg_size * 2

  ! Add 1 for tile number and 1 for upper/lower
  msg_size = msg_size + 2

  ! Store nr for later use
  msg%recv_vdf_range( :, 1 : vdf%x_dim_, til_idx ) = nr ( :, 1 : vdf%x_dim_ )

  ! Store message size
  msg%msg_size(p_recv_msg,til_idx) = msg_size

end subroutine irecv_size_vdf
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Wait for a BC vdf recv message to complete
!-----------------------------------------------------------------------------------------
subroutine loc_recv_vdf( vdf, bnd, recv_msg, update_type, dim, &
                         nx_move, vdf_b )

  implicit none

  type( t_vdf ), intent(inout) :: vdf
  integer, intent(in) :: bnd
  type(t_vdf_msg), dimension(2), intent(inout) :: recv_msg
  integer, intent(in) :: update_type, dim
  integer, intent(in) :: nx_move
  type( t_vdf ), intent(inout), optional :: vdf_b

  integer, dimension(2,p_max_dim) :: nr
  integer :: i, msg_size

  ! Get range of cells to receive
  call get_recv_range( vdf, msg_size, nr, bnd, update_type, dim, nx_move )

  recv_msg(bnd)%type = update_type
  recv_msg(bnd)%range( :, 1 : vdf%x_dim_ ) = nr ( :, 1 : vdf%x_dim_ )

  ! unpack message data
  i = 1
  call unpack_vdf_msg( recv_msg( bnd ), vdf, i )
  if ( present( vdf_b ) ) call unpack_vdf_msg( recv_msg( bnd ), vdf_b, i )

end subroutine loc_recv_vdf
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Returns the size of the buffer for a single field point
!-----------------------------------------------------------------------------------------
function fld_pack_size( )

  implicit none

  integer :: fld_pack_size_, fld_pack_size

  integer :: fld_mpi_type, ierr

  fld_mpi_type = mpi_real_type( p_k_fld )

  call mpi_pack_size( 1, fld_mpi_type, mpi_comm_world, fld_pack_size_, ierr)

  if ( ierr /= 0 ) then
    write(err_buf__,*) "MPI error: fld_pack_size";call err__("tiles/os-vdf-comm-tiles.f90",400)
    call abort_program( p_err_mpi )
  endif
  fld_pack_size = fld_pack_size_

end function fld_pack_size
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Unapack grouped message buffer data into vdf
!-----------------------------------------------------------------------------------------
subroutine unpack_vdf_msg_group( msg, vdf, k, rng )

  implicit none

  type( t_vdf_msg ), intent(in) :: msg
  type( t_vdf ), intent(inout) :: vdf
  integer, intent(inout) :: k
  integer, dimension(2,p_max_dim), intent(in) :: rng

  integer :: i, i1, i2, i3
  real(p_k_fld), dimension(:,:), pointer :: f1
  real(p_k_fld), dimension(:,:,:), pointer :: f2
  real(p_k_fld), dimension(:,:,:,:), pointer :: f3
  real(p_k_fld), dimension(:), pointer :: buffer

  buffer => msg%buffer

  if (vdf%f_dim_ == 3) then

    select case ( msg%type )

    case ( p_vdf_add )

      ! Add message values to local values

      select case ( vdf%x_dim_ )
      case (1)
        f1 => vdf%f1
        do i1 = rng( p_lower, 1 ), rng( p_upper, 1 )
          f1( 1, i1 ) = f1( 1, i1 ) + buffer( k )
          f1( 2, i1 ) = f1( 2, i1 ) + buffer( k+1 )
          f1( 3, i1 ) = f1( 3, i1 ) + buffer( k+2 )
          k = k + 3
        enddo

      case (2)
        f2 => vdf%f2
        do i2 = rng( p_lower, 2 ), rng( p_upper, 2 )
          do i1 = rng( p_lower, 1 ), rng( p_upper, 1 )
            f2( 1, i1, i2 ) = f2( 1, i1, i2 ) + buffer( k )
            f2( 2, i1, i2 ) = f2( 2, i1, i2 ) + buffer( k+1 )
            f2( 3, i1, i2 ) = f2( 3, i1, i2 ) + buffer( k+2 )
            k = k + 3
          enddo
        enddo

      case (3)
        f3 => vdf%f3
        do i3 = rng( p_lower, 3 ), rng( p_upper, 3 )
          do i2 = rng( p_lower, 2 ), rng( p_upper, 2 )
            do i1 = rng( p_lower, 1 ), rng( p_upper, 1 )
              f3( 1, i1, i2, i3 ) = f3( 1, i1, i2, i3 ) + buffer( k )
              f3( 2, i1, i2, i3 ) = f3( 2, i1, i2, i3 ) + buffer( k+1 )
              f3( 3, i1, i2, i3 ) = f3( 3, i1, i2, i3 ) + buffer( k+2 )
              k = k + 3
            enddo
          enddo
        enddo
      end select

    case ( p_vdf_replace )

      ! replace local values with message values

      select case ( vdf%x_dim_ )
      case (1)
        f1 => vdf%f1
        do i1 = rng( p_lower, 1 ), rng( p_upper, 1 )
          f1( 1, i1 ) = buffer( k )
          f1( 2, i1 ) = buffer( k+1 )
          f1( 3, i1 ) = buffer( k+2 )
          k = k + 3
        enddo

      case (2)
        f2 => vdf%f2
        do i2 = rng( p_lower, 2 ), rng( p_upper, 2 )
          do i1 = rng( p_lower, 1 ), rng( p_upper, 1 )
            f2( 1, i1, i2 ) = buffer( k )
            f2( 2, i1, i2 ) = buffer( k+1 )
            f2( 3, i1, i2 ) = buffer( k+2 )
            k = k + 3
          enddo
        enddo

      case (3)
        f3 => vdf%f3
        do i3 = rng( p_lower, 3 ), rng( p_upper, 3 )
          do i2 = rng( p_lower, 2 ), rng( p_upper, 2 )
            do i1 = rng( p_lower, 1 ), rng( p_upper, 1 )
              f3( 1, i1, i2, i3 ) = buffer( k )
              f3( 2, i1, i2, i3 ) = buffer( k+1 )
              f3( 3, i1, i2, i3 ) = buffer( k+2 )
              k = k + 3
            enddo
          enddo
        enddo
      end select

    case default
      write(err_buf__,*) 'unpack_vdf_msg_group called with an invalid msg%type';call err__("tiles/os-vdf-comm-tiles.f90",511)
      call abort_program()

    end select

  else

    select case ( msg%type )

    case ( p_vdf_add )

      ! Add message values to local values

      select case ( vdf%x_dim_ )
      case (1)
        f1 => vdf%f1
        do i1 = rng( p_lower, 1 ), rng( p_upper, 1 )
          do i = 1, vdf%f_dim_
            f1( i, i1 ) = f1( i, i1 ) + buffer( k )
            k = k + 1
          enddo
        enddo

      case (2)
        f2 => vdf%f2
        do i2 = rng( p_lower, 2 ), rng( p_upper, 2 )
          do i1 = rng( p_lower, 1 ), rng( p_upper, 1 )
            do i = 1, vdf%f_dim_
              f2( i, i1, i2 ) = f2( i, i1, i2 ) + buffer( k )
              k = k + 1
            enddo
          enddo
        enddo

      case (3)
        f3 => vdf%f3
        do i3 = rng( p_lower, 3 ), rng( p_upper, 3 )
          do i2 = rng( p_lower, 2 ), rng( p_upper, 2 )
            do i1 = rng( p_lower, 1 ), rng( p_upper, 1 )
              do i = 1, vdf%f_dim_
                f3( i, i1, i2, i3 ) = f3( i, i1, i2, i3 ) + buffer( k )
                k = k + 1
              enddo
            enddo
          enddo
        enddo
      end select

    case ( p_vdf_replace )

      ! replace local values with message values

      select case ( vdf%x_dim_ )
      case (1)
        f1 => vdf%f1
        do i1 = rng( p_lower, 1 ), rng( p_upper, 1 )
          do i = 1, vdf%f_dim_
            f1( i, i1 ) = buffer( k )
            k = k + 1
          enddo
        enddo

      case (2)
        f2 => vdf%f2
        do i2 = rng( p_lower, 2 ), rng( p_upper, 2 )
          do i1 = rng( p_lower, 1 ), rng( p_upper, 1 )
            do i = 1, vdf%f_dim_
              f2( i, i1, i2 ) = buffer( k )
              k = k + 1
            enddo
          enddo
        enddo

      case (3)
        f3 => vdf%f3
        do i3 = rng( p_lower, 3 ), rng( p_upper, 3 )
          do i2 = rng( p_lower, 2 ), rng( p_upper, 2 )
            do i1 = rng( p_lower, 1 ), rng( p_upper, 1 )
              do i = 1, vdf%f_dim_
                f3( i, i1, i2, i3 ) = buffer( k )
                k = k + 1
              enddo
            enddo
          enddo
        enddo
      end select

    case default
      write(err_buf__,*) 'unpack_vdf_msg_group called with an invalid msg%type';call err__("tiles/os-vdf-comm-tiles.f90",599)
      call abort_program()

    end select

  endif

end subroutine unpack_vdf_msg_group
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Pack grouped vdf data into message buffer
!-----------------------------------------------------------------------------------------
subroutine pack_vdf_msg_group( msg, vdf, k, rng, tile, bnd )

  implicit none

  type( t_vdf_msg ), intent(inout) :: msg
  type( t_vdf ), intent(in) :: vdf
  integer, intent(inout) :: k
  integer, dimension(2,p_max_dim), intent(in) :: rng
  integer, intent(in), optional :: tile, bnd

  integer :: i, i1, i2, i3
  real(p_k_fld), dimension(:,:), pointer :: f1
  real(p_k_fld), dimension(:,:,:), pointer :: f2
  real(p_k_fld), dimension(:,:,:,:), pointer :: f3
  real(p_k_fld), dimension(:), pointer :: buffer

  buffer => msg%buffer

  ! Pack the preface, consisting of tile number and upper/lower
  if (present(tile) .and. present(bnd)) then
    buffer( k ) = real(tile,kind=p_k_fld)
    k = k + 1
    buffer( k ) = real(bnd,kind=p_k_fld)
    k = k + 1
  endif

  if (vdf%f_dim_ == 3) then

    select case ( vdf%x_dim_ )
    case (1)
      f1 => vdf%f1
      do i1 = rng( p_lower, 1 ), rng( p_upper, 1 )
        buffer( k ) = f1( 1, i1 )
        buffer( k+1 ) = f1( 2, i1 )
        buffer( k+2 ) = f1( 3, i1 )
        k = k + 3
      enddo

    case (2)
      f2 => vdf%f2
      do i2 = rng( p_lower, 2 ), rng( p_upper, 2 )
        do i1 = rng( p_lower, 1 ), rng( p_upper, 1 )
          buffer( k ) = f2( 1, i1, i2 )
          buffer( k+1 ) = f2( 2, i1, i2 )
          buffer( k+2 ) = f2( 3, i1, i2 )
          k = k + 3
        enddo
      enddo

    case (3)
      f3 => vdf%f3
      do i3 = rng( p_lower, 3 ), rng( p_upper, 3 )
        do i2 = rng( p_lower, 2 ), rng( p_upper, 2 )
          do i1 = rng( p_lower, 1 ), rng( p_upper, 1 )
            buffer( k ) = f3( 1, i1, i2, i3 )
            buffer( k+1 ) = f3( 2, i1, i2, i3 )
            buffer( k+2 ) = f3( 3, i1, i2, i3 )
            k = k + 3
          enddo
        enddo
      enddo
    end select

  else

    select case ( vdf%x_dim_ )
    case (1)
      f1 => vdf%f1
      do i1 = rng( p_lower, 1 ), rng( p_upper, 1 )
        do i = 1, vdf%f_dim_
          buffer( k ) = f1( i, i1 )
          k = k + 1
        enddo
      enddo

    case (2)
      f2 => vdf%f2
      do i2 = rng( p_lower, 2 ), rng( p_upper, 2 )
        do i1 = rng( p_lower, 1 ), rng( p_upper, 1 )
          do i = 1, vdf%f_dim_
            buffer( k ) = f2( i, i1, i2 )
            k = k + 1
          enddo
        enddo
      enddo

    case (3)
      f3 => vdf%f3
      do i3 = rng( p_lower, 3 ), rng( p_upper, 3 )
        do i2 = rng( p_lower, 2 ), rng( p_upper, 2 )
          do i1 = rng( p_lower, 1 ), rng( p_upper, 1 )
            do i = 1, vdf%f_dim_
              buffer( k ) = f3( i, i1, i2, i3 )
              k = k + 1
            enddo
          enddo
        enddo
      enddo
    end select

  endif


end subroutine pack_vdf_msg_group
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! post send message for boundary communication
!-----------------------------------------------------------------------------------------
subroutine isend_vdf_tiles( vdf, update_type, dim, bnd, no_co, nx_move, send_msg, vdf_b )

  implicit none

  type( t_vdf ), intent(in) :: vdf
  type( t_vdf ), intent(in), optional :: vdf_b

  integer, intent(in) :: update_type, dim, bnd
  class( t_node_conf ), intent(in) :: no_co

  integer, intent(in) :: nx_move
  type(t_vdf_msg), dimension(2), intent(inout) :: send_msg

  integer, dimension(2,p_max_dim) :: ns
  integer :: i, bsize, msg_size

  ! Sanity checks
  if ( send_msg(bnd)%handle /= MPI_REQUEST_NULL ) then
    write(err_buf__,*) 'isend_msg_vdf called when message still waiting to complete';call err__("tiles/os-vdf-comm-tiles.f90",739)
    call abort_program( p_err_invalid )
  endif

  if ( neighbor( no_co, bnd, dim ) < 0 ) then
    write(0,*) 'irecv_vdf called without a valid neighbor, dim, bnd = ',dim, bnd, &
               trim("tiles/os-vdf-comm-tiles.f90"), ':', 745
    call abort_program( p_err_invalid )
  endif

  ! Get range of cells to send
  call get_send_range( vdf, msg_size, ns, bnd, update_type, dim, nx_move )

  ! If second vdf is present double message size
  if ( present( vdf_b ) ) msg_size = msg_size * 2

  send_msg(bnd)%type = update_type ! unnecessary
  send_msg(bnd)%range( :, 1 : vdf%x_dim_ ) = ns ( :, 1 : vdf%x_dim_ )
  send_msg(bnd)%msg_size = msg_size
  send_msg(bnd)%node = neighbor( no_co, bnd, dim )
  select type(no_co); class is(t_node_conf_tiles)
    send_msg(bnd)%tag = (no_co%ti_co%neighbor_g_til(bnd,dim)-1)*2*vdf%x_dim_+(2-bnd)+2*(dim-1)
  end select

  ! grow message buffer if necessary
  if ( associated( send_msg(bnd)%buffer ) ) then
    bsize = size( send_msg(bnd)%buffer )
  else
    bsize = 0
  endif

  if ( msg_size > bsize ) then
    call freemem(send_msg(bnd)%buffer,"tiles/os-vdf-comm-tiles.f90",771)
    call alloc(send_msg(bnd)%buffer, (/msg_size/),"tiles/os-vdf-comm-tiles.f90",772)
  endif

  ! pack message data
  i = 1
  call pack_vdf_msg( send_msg(bnd), vdf, i )
  if ( present( vdf_b ) ) call pack_vdf_msg( send_msg(bnd), vdf_b, i )

  ! post send
  call isend_vdf_msg_tiles( send_msg(bnd), no_co )

end subroutine isend_vdf_tiles
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! post recv message for boundary communication
!-----------------------------------------------------------------------------------------
subroutine irecv_vdf_tiles( vdf, update_type, dim, bnd, no_co, nx_move, recv_msg, vdf_b )

  implicit none

  type( t_vdf ), intent(in) :: vdf
  integer, intent(in) :: update_type, dim, bnd
  class( t_node_conf ), intent(in) :: no_co
  integer, intent(in) :: nx_move
  type(t_vdf_msg), dimension(2), intent(inout) :: recv_msg
  type( t_vdf ), intent(in), optional :: vdf_b


  integer, dimension(2,p_max_dim) :: nr
  integer :: bsize, msg_size

  ! Check if message is still waiting to complete
  if ( recv_msg( bnd )%handle /= MPI_REQUEST_NULL ) then
    write(err_buf__,*) 'irecv_vdf called when message still waiting to complete';call err__("tiles/os-vdf-comm-tiles.f90",806)
    call abort_program( p_err_invalid )
  endif

  ! Sanity check
  if ( neighbor( no_co, bnd, dim ) < 0 ) then
    write(0,*) 'irecv_vdf called without a valid neighbor, dim, bnd = ',dim, bnd, &
                    trim("tiles/os-vdf-comm-tiles.f90"), ':', 813
    call abort_program( p_err_invalid )
  endif

  ! Get range of cells to receive
  call get_recv_range( vdf, msg_size, nr, bnd, update_type, dim, nx_move )

  ! (* debug *)
  if ( msg_size == 0 ) then
    write(err_buf__,*) 'Invalid message size, aborting.';call err__("tiles/os-vdf-comm-tiles.f90",822)
    call abort_program( p_err_invalid )
  endif

  ! If second vdf is present double message size
  if ( present( vdf_b ) ) msg_size = msg_size * 2

  ! Get send buffer and grow it if necessary
  recv_msg(bnd)%type = update_type
  recv_msg(bnd)%range( :, 1 : vdf%x_dim_ ) = nr ( :, 1 : vdf%x_dim_ )
  recv_msg(bnd)%msg_size = msg_size
  recv_msg(bnd)%node = neighbor( no_co, bnd, dim )
  select type(no_co); class is(t_node_conf_tiles)
    recv_msg(bnd)%tag = (no_co%ti_co%g_til_aid-1)*2*vdf%x_dim_+(bnd-1)+2*(dim-1)
  end select

  ! grow message buffer if necessary
  if ( associated( recv_msg(bnd)%buffer ) ) then
    bsize = size( recv_msg(bnd)%buffer )
  else
    bsize = 0
  endif

  if ( msg_size > bsize ) then
    if ( bsize > 0 ) call freemem(recv_msg(bnd)%buffer,"tiles/os-vdf-comm-tiles.f90",846)
    call alloc(recv_msg(bnd)%buffer, (/msg_size/),"tiles/os-vdf-comm-tiles.f90",847)
  endif

  ! post receive
  call irecv_vdf_msg_tiles( recv_msg(bnd), no_co )

end subroutine irecv_vdf_tiles
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
! Post receive for message, ignoring if recv node is same as my node
!-----------------------------------------------------------------------------------------
subroutine irecv_vdf_msg_tiles( recv, no_co )

  implicit none

  type (t_vdf_msg), intent(inout) :: recv
  class (t_node_conf ), intent(in) :: no_co

  integer :: ierr
  integer :: mpi_type

  mpi_type = mpi_real_type( p_k_fld )

  call mpi_irecv( recv%buffer, recv%msg_size, mpi_type, recv%node - 1, recv%tag, comm( no_co ), &
    recv%handle, ierr )

  if (ierr/=0) then
    write(err_buf__,*) "MPI error";call err__("tiles/os-vdf-comm-tiles.f90",876)
    call abort_program( p_err_mpi )
  endif
  ! print *, '[', my_node, '] <- ', recv%node , ' ', recv%msg_size, ' reals, handle = ', recv%handle, ' tag = ', recv%tag

end subroutine irecv_vdf_msg_tiles
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Post send for message, ignoring if send node is same as my node
!-----------------------------------------------------------------------------------------
subroutine isend_vdf_msg_tiles( send, no_co )

  implicit none

  type (t_vdf_msg), intent(inout) :: send
  class (t_node_conf ), intent(in) :: no_co

  integer :: ierr
  integer :: mpi_type

  mpi_type = mpi_real_type( p_k_fld )

  call mpi_isend( send%buffer, send%msg_size, mpi_type, send%node-1, send%tag, comm(no_co), &
    send%handle, ierr )
  if (ierr/=0) then
    write(err_buf__,*) "MPI error";call err__("tiles/os-vdf-comm-tiles.f90",902)
    call abort_program( p_err_mpi )
  endif
  ! print *, '[', my_node, '] -> ', send%node , ' ', send%msg_size, ' reals, handle = ', send%handle, ' tag = ', send%tag

end subroutine isend_vdf_msg_tiles
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Get range of cells to send
!-----------------------------------------------------------------------------------------
subroutine get_send_range( vdf, msg_size, ns, bnd, update_type, dim, nx_move )

  implicit none

  type( t_vdf ), intent(in) :: vdf
  integer, intent(inout) :: msg_size
  integer, dimension(2,p_max_dim), intent(inout) :: ns
  integer, intent(in) :: bnd, update_type, dim, nx_move

  integer :: i

  msg_size = vdf%f_dim_
  do i = 1, vdf%x_dim_
    if ( i == dim ) then
      if ( update_type == p_vdf_add ) then
        select case ( bnd )
        case ( p_upper )
          ns( p_lower, i ) = vdf%nx_(i) + 1 - vdf%gc_num_( p_lower, i )
          ns( p_upper, i ) = vdf%nx_(i) + vdf%gc_num_( p_upper, i )
        case ( p_lower )
          ns( p_lower, i ) = 1 - vdf%gc_num_( p_lower, i )
          ns( p_upper, i ) = vdf%gc_num_( p_upper, i )
        end select
      else
        select case ( bnd )
        case ( p_upper )
          ns( p_lower, i ) = vdf%nx_(i) + 1 - vdf%gc_num_( p_lower, i )
          ns( p_upper, i ) = vdf%nx_(i) - nx_move
        case ( p_lower )
          ns( p_lower, i ) = 1 - nx_move
          ns( p_upper, i ) = vdf%gc_num_( p_upper, i )
        end select
      endif
    else
      ns( p_lower, i ) = vdf%lbound( i+1 )
      ns( p_upper, i ) = vdf%ubound( i+1 )
    endif
    msg_size = msg_size * (ns( p_upper, i ) - ns( p_lower,i ) + 1 )
  enddo

end subroutine get_send_range
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Get range of cells to receive
!-----------------------------------------------------------------------------------------
subroutine get_recv_range( vdf, msg_size, nr, bnd, update_type, dim, nx_move )

  implicit none

  type( t_vdf ), intent(in) :: vdf
  integer, intent(inout) :: msg_size
  integer, dimension(2,p_max_dim), intent(inout) :: nr
  integer, intent(in) :: bnd, update_type, dim, nx_move

  integer :: i

  msg_size = vdf%f_dim_
  do i = 1, vdf%x_dim_
    if ( i == dim ) then
      if ( update_type == p_vdf_add ) then
        select case ( bnd )
        case ( p_upper )
          nr( p_lower, i ) = vdf%nx_(i) + 1 - vdf%gc_num_( p_lower, i )
          nr( p_upper, i ) = vdf%nx_(i) + vdf%gc_num_( p_upper, i )
        case ( p_lower )
          nr( p_lower, i ) = 1 - vdf%gc_num_( p_lower, i )
          nr( p_upper, i ) = vdf%gc_num_( p_upper, i )
        end select
      else
        select case ( bnd )
        case ( p_upper )
          nr( p_lower, i ) = vdf%nx_(i) + 1 - nx_move
          nr( p_upper, i ) = vdf%nx_(i) + vdf%gc_num_( p_upper, i )
        case ( p_lower )
          nr( p_lower, i ) = 1 - vdf%gc_num_( p_lower, i )
          nr( p_upper, i ) = 0 - nx_move
        end select
      endif
    else
      nr( p_lower, i ) = vdf%lbound( i+1 )
      nr( p_upper, i ) = vdf%ubound( i+1 )
    endif
    msg_size = msg_size * (nr( p_upper, i ) - nr( p_lower,i ) + 1 )
  enddo

end subroutine get_recv_range
!-----------------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
subroutine pack_fld_msg_dlb_group( f1, f2, f3, fld_size, lbound_fld, ubound_fld, buffer, &
                                  bsize, position)

  implicit none

  real(p_k_fld), dimension(:,:), pointer, intent(in) :: f1
  real(p_k_fld), dimension(:,:,:), pointer, intent(in) :: f2
  real(p_k_fld), dimension(:,:,:,:), pointer, intent(in) :: f3
  integer, intent(in) :: fld_size
  integer, dimension(:), intent(in) :: lbound_fld, ubound_fld
  integer(p_byte), dimension(:), pointer, intent(in) :: buffer
  integer, intent(in) :: bsize
  integer, intent(inout) :: position


  real(p_k_fld), dimension(:), pointer :: temp_buffer => null()
  integer :: fld_mpi_type, ierr
  integer :: i, i1, i2, i3, k

  fld_mpi_type = mpi_real_type( p_k_fld )

  k = 1

  call alloc(temp_buffer, (/fld_size/),"tiles/os-vdf-comm-tiles.f90",1027)

  select case ( p_x_dim )
  case (1)
    do i1 = lbound_fld(2), ubound_fld(2)
      do i = lbound_fld(1), ubound_fld(1)
        temp_buffer( k ) = f1( i, i1 )
        k = k + 1
      enddo
    enddo
  case (2)
    do i2 = lbound_fld(3), ubound_fld(3)
      do i1 = lbound_fld(2), ubound_fld(2)
        do i = lbound_fld(1), ubound_fld(1)
          temp_buffer( k ) = f2( i, i1, i2 )
          k = k + 1
        enddo
      enddo
    enddo
  case (3)
    do i3 = lbound_fld(4), ubound_fld(4)
      do i2 = lbound_fld(3), ubound_fld(3)
        do i1 = lbound_fld(2), ubound_fld(2)
          do i = lbound_fld(1), ubound_fld(1)
            temp_buffer( k ) = f3( i, i1, i2, i3 )
            k = k + 1
          enddo
        enddo
      enddo
    enddo
  end select

  call mpi_pack( temp_buffer, fld_size, fld_mpi_type, buffer, &
    bsize, position, mpi_comm_world, ierr )

  call freemem(temp_buffer,"tiles/os-vdf-comm-tiles.f90",1062)

end subroutine pack_fld_msg_dlb_group
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
subroutine unpack_fld_msg_dlb_group(f1, f2, f3, fld_size, lbound_fld, ubound_fld, buffer, &
                                    bsize, position)

  implicit none

  real(p_k_fld), dimension(:,:), pointer, intent(in) :: f1
  real(p_k_fld), dimension(:,:,:), pointer, intent(in) :: f2
  real(p_k_fld), dimension(:,:,:,:), pointer, intent(in) :: f3
  integer, intent(in) :: fld_size
  integer, dimension(:), intent(in) :: lbound_fld, ubound_fld
  integer(p_byte), dimension(:), pointer, intent(in) :: buffer
  integer, intent(in) :: bsize
  integer, intent(inout) :: position


  real(p_k_fld), dimension(:), pointer :: temp_buffer => null()
  integer :: fld_mpi_type, ierr
  integer :: i, i1, i2, i3, k

  fld_mpi_type = mpi_real_type( p_k_fld )

  k = 1

  call alloc(temp_buffer, (/fld_size/),"tiles/os-vdf-comm-tiles.f90",1093)

  call mpi_unpack( buffer, bsize, position, temp_buffer, fld_size, &
    fld_mpi_type, mpi_comm_world, ierr )

  select case ( p_x_dim )
  case (1)
    do i1 = lbound_fld(2), ubound_fld(2)
      do i = lbound_fld(1), ubound_fld(1)
        f1( i, i1 ) = temp_buffer( k )
        k = k + 1
      enddo
    enddo
  case (2)
    do i2 = lbound_fld(3), ubound_fld(3)
      do i1 = lbound_fld(2), ubound_fld(2)
        do i = lbound_fld(1), ubound_fld(1)
          f2( i, i1, i2 ) = temp_buffer( k )
          k = k + 1
        enddo
      enddo
    enddo
  case (3)
    do i3 = lbound_fld(4), ubound_fld(4)
      do i2 = lbound_fld(3), ubound_fld(3)
        do i1 = lbound_fld(2), ubound_fld(2)
          do i = lbound_fld(1), ubound_fld(1)
            f3( i, i1, i2, i3 ) = temp_buffer( k )
            k = k + 1
          enddo
        enddo
      enddo
    enddo
  end select

  call freemem(temp_buffer,"tiles/os-vdf-comm-tiles.f90",1128)

end subroutine unpack_fld_msg_dlb_group
!-------------------------------------------------------------------------------


end module m_vdf_comm_tiles
