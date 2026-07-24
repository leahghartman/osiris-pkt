module m_group_comm

#include "os-config.h"
#include "os-preprocess.fpp"
#include "memory/memory.h"

use m_system
use m_parameters

use m_grid_define,         only: t_grid, t_msg_patt, t_msg
use m_grid,                only: coordinates
use m_grid_tiles,          only: t_grid_group
use m_space,               only: t_space, extend
use m_simulation_group,    only: t_simulation_group, t_til_container
use m_simulation_tiles,    only: t_simulation_tiles
use m_simulation,          only: t_simulation
use m_node_conf,           only: t_node_conf
use m_node_conf_tiles,     only: t_node_conf_tiles, p_spec_buf_block_t
use m_tile_conf,           only: p_topology_viktor, g_tgp
use m_vdf_define,          only: t_vdf, t_vdf_arr
use m_vdf_comm,            only: t_vdf_msg, wait, irecv, isend, update_periodic_vdf, irecv_wait, isend_wait, wait
use m_vdf_math,            only: add
use m_vdf_comm_tiles,      only: t_bnd_msg, fld_pack_size, pack_group, unpack_group
use m_vdf_comm_tiles,      only: irecv_tiles, isend_tiles, loc_send_vdf, get_recv_range
use m_vdf_comm_tiles,      only: loc_recv_vdf, get_send_range
use m_vdf_comm_tiles,      only: p_send_msg, p_recv_msg
use m_vdf_comm_tiles,      only: pack_fld_msg_dlb_group, unpack_fld_msg_dlb_group
use m_species_comm,        only: comm, int_pack_size, particle_pack_size, neighbor, wait_spec_msg, wait_spec_msg_1, pack_single_particle, unpack_single_particle
use m_species_comm_tiles,  only: pack_group, unpack
use m_species_define,      only: t_species, t_spec_msg, t_spec_arr
use m_wall_define,         only: t_wall, ubound, lbound
use m_emf_define,          only: t_emf, t_vpml
use m_species_charge,      only: norm_charge_cyl, deposit_rho
use m_particles_define,    only: t_particles
use m_species_diagnostics, only: p_norm, p_charge


implicit none

private

interface irecv_til_size
  module procedure irecv_til_size
end interface

interface isend_til_size
  module procedure isend_til_size
end interface

interface irecv_til_msg
  module procedure irecv_til_msg_dlb
  module procedure irecv_til_msg_diag
end interface

interface isend_til_msg
  module procedure isend_til_msg_dlb
  module procedure isend_til_msg_diag
end interface

interface irecv_wait
  module procedure wait_irecv_vdf_group
end interface

interface isend
  module procedure isend_vdf_group
end interface

interface isend_size_spec
  module procedure isend_size_spec_group
end interface

interface isend_1
  module procedure isend_1_spec_group
end interface

interface wait_send
  module procedure wait_send_msg
end interface

interface wait_free_send
  module procedure wait_free_send_diag
end interface

interface wait_free_recv
  module procedure wait_free_recv_diag
end interface

interface unpack
  module procedure unpack_tiles_dlb
  ! module procedure unpack_tiles_diag
  module procedure unpack_spec_group
end interface

interface get_size
  module procedure get_size_preface
  module procedure get_size_emf
  module procedure get_size_part
end interface

interface deposit_quant_1
  module procedure deposit_quant_1_tiles
end interface

interface deposit_quant_2
  module procedure deposit_quant_2_tiles
end interface

interface deposit_charge_1
  module procedure deposit_charge_part_1
end interface

interface deposit_charge_2
  module procedure deposit_charge_part_2
end interface

interface loc_fill
  module procedure loc_fill_vdf_tiles
end interface

public :: irecv_til_size, isend_til_size, irecv_til_msg, isend_til_msg, wait_send
public :: wait_free_send, wait_free_recv, unpack, irecv_wait, isend
public :: deposit_quant_1, deposit_quant_2
public :: deposit_charge_1, deposit_charge_2, isend_size_spec
public :: isend_1, loc_fill

contains

!-------------------------------------------------------------------------------
! receives integers containing the total size in memory of msgs it will
! receive from each node.  Puts this into recv_bsize
!-------------------------------------------------------------------------------
subroutine irecv_til_size( sim, patt, bsize )

  implicit none

  class( t_simulation_group ), intent(in)       :: sim
  type( t_msg_patt ), intent(inout)             :: patt
  integer, dimension(:), pointer, intent(inout) :: bsize

  integer :: ierr
  integer :: i

  call alloc( bsize, (/patt%n_recv/) )

  ierr  = 0
  bsize = 0

  do i = 1, patt%n_recv

    patt%recv(i)%request = MPI_REQUEST_NULL

    ! receive size of ith message, store in buffer bsize
    call mpi_irecv( bsize(i), 1, MPI_INTEGER, patt%recv(i)%node-1, 0, &
                    comm(sim%no_co), patt%recv(i)%request, ierr )

    if (ierr/=0) then
      ERROR("MPI error")
      call abort_program( p_err_mpi )
    endif
  enddo

end subroutine irecv_til_size
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! sends int to each node containing msg size needed to send tiles
!-------------------------------------------------------------------------------
subroutine isend_til_size( sim, patt, bsize)

  implicit none

  class( t_simulation_group ), intent(in)       :: sim
  type( t_msg_patt ), intent(inout)             :: patt
  integer, dimension(:), pointer, intent(inout) :: bsize

  class(t_simulation), pointer :: til
  integer, dimension(p_x_dim) :: tgp
  integer :: ierr
  integer :: i, j, k
  integer :: g_til_aid_

  call alloc( bsize, (/patt%n_send/) )

  ierr  = 0
  bsize = 0

  select type(no_co => sim%no_co); class is(t_node_conf_tiles)

  do i = 1, patt%n_send
    ! compute size of ith message
    do j = 1, size( patt%send(i)%tils )
      select case ( no_co%ti_co%topology_type )
        case (p_topology_viktor)
          do k = 1, p_x_dim
            tgp(k) = g_tgp( no_co%ti_co, patt%send(i)%tils(j), k )
          enddo
          if ( p_x_dim==2 ) g_til_aid_ = no_co%ti_co%hi2_old( tgp(1), tgp(2) )
          if ( p_x_dim==3 ) g_til_aid_ = no_co%ti_co%hi3_old( tgp(1), tgp(2), tgp(3) )
        case default
          g_til_aid_ = patt%send(i)%tils(j)
      end select

      til => sim % tiles(g_til_aid_) % t

      call get_size( til%part%num_species, bsize(i) )
      call get_size( til%emf, bsize(i) )
      call get_size( til%part, bsize(i) )
    enddo

    if (bsize(i) < 0) then
      SCR_MPINODE('Invalid value for tile buffer size, possible integer overflow.')
      SCR_MPINODE('Message to node ',patt%send(i)%node-1)
      call abort_program( p_err_invalid )
    endif

    patt%send(i)%request = MPI_REQUEST_NULL

    ! send ith message to proper process
    call mpi_isend( bsize(i), 1, MPI_INTEGER, patt%send(i)%node-1, &
                    0, comm(no_co), patt%send(i)%request, ierr )

    if (ierr/=0) then
      ERROR("MPI error")
      call abort_program( p_err_mpi )
    endif

  enddo

  end select

end subroutine isend_til_size
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
subroutine irecv_til_msg_dlb( patt, bsize, tiles )

  implicit none

  type( t_msg_patt ), intent(inout)                           :: patt
  integer, dimension(:), pointer, intent(inout)               :: bsize
  class( t_til_container ), dimension(:), pointer, intent(in) :: tiles

  integer :: i, ierr

  ! Wait for message size to be received
  call wait_msg_array( patt%recv, patt%n_recv )

  do i = 1, patt%n_recv

    ! Allocate message buffers
    call alloc( patt%recv(i)%buffer, (/ bsize(i) /) )

    patt%recv(i)%request = MPI_REQUEST_NULL

    ! Post receive
    call mpi_irecv( patt%recv(i)%buffer, bsize(i), MPI_PACKED, patt%recv(i)%node-1, &
                    0, comm(tiles(lbound(tiles,1))%t%no_co), patt%recv(i)%request, ierr )

    if (ierr/=0) then
      ERROR("MPI error")
      call abort_program( p_err_mpi )
    endif

  enddo

end subroutine irecv_til_msg_dlb
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
subroutine isend_til_msg_dlb( no_co, patt, bsize, tiles )

  implicit none

  class( t_node_conf ), intent(in)                            :: no_co
  type( t_msg_patt ), intent(inout)                           :: patt
  integer, dimension(:), pointer, intent(inout)               :: bsize
  class( t_til_container ), dimension(:), pointer, intent(in) :: tiles

  integer :: position, g_til_aid_, num_spec, fld_mpi_type, part_mpi_type, fld_size
  integer :: i, j, k, ierr, np
  integer, dimension(p_x_dim) :: tgp
  integer, dimension(:), pointer :: preface
  integer, dimension(p_max_dim) :: lshift
  class(t_species), pointer :: spec
  type(t_vdf), pointer :: vdf
  type(t_vpml), pointer :: vpml
  type(t_wall), pointer :: wall

  call wait_msg_array( patt%send, patt%n_send )

  num_spec = tiles(lbound(tiles,1)) % t % part % num_species
  ! Array containing g_til_aid and num_par for each species
  call alloc( preface, (/ num_spec+1 /) )

  fld_mpi_type = mpi_real_type( p_k_fld )
  part_mpi_type = mpi_real_type( p_k_part )

  do i = 1, patt%n_send

    ! Allocate send buffer
    call alloc( patt%send(i)%buffer, (/ bsize(i) /) )

    ! Pack tiles
    position = 0

    ! Pack all tiles going to the same node serially, one after the other
    do j = 1, size( patt%send(i)%tils )

      ! We want to use the post dlb g_til_aid index (only important in viktor topology)
      ! in the preface since this gets sent in the mpi message and is what the 
      ! receiving node uses to know which tile it is receiving
      g_til_aid_ = patt%send(i)%tils(j)
      preface(1) = g_til_aid_

      ! from this point on we want to use the pre dlb g_til_aid index (only important
      ! in viktor topology), since we use it to grab from the old tile array which is 
      ! indexed in the old way
      select type(no_co); class is(t_node_conf_tiles)
      select case ( no_co%ti_co%topology_type )
        case (p_topology_viktor)
          ! for viktor topology, this makes a difference
          do k = 1, p_x_dim
            tgp(k) = g_tgp( no_co%ti_co, patt%send(i)%tils(j), k )
          enddo
          if ( p_x_dim==2 ) g_til_aid_ = no_co%ti_co%hi2_old( tgp(1), tgp(2) )
          if ( p_x_dim==3 ) g_til_aid_ = no_co%ti_co%hi3_old( tgp(1), tgp(2), tgp(3) )
        case default
          ! for all other topologies, pre and post dlb g_til_aid's are the same always
      end select
      end select

      ! Add number of particles for each species to preface
      spec => tiles(g_til_aid_) % t % part % species
      k = 1
      do
        if (.not. associated(spec)) exit
        k = k + 1
        preface(k) = spec % num_par
        spec => spec % next
      enddo

      ! First pack preface, containing g_til_aid and num_par for each species
      call mpi_pack( preface, num_spec+1, MPI_INTEGER, patt%send(i)%buffer, &
        bsize(i), position, mpi_comm_world, ierr )

      ! Pack e
      vdf => tiles(g_til_aid_)%t%emf%e
      fld_size = vdf%f_dim() * product(vdf%size())

      call pack_fld_msg_dlb_group( vdf%f1, vdf%f2, vdf%f3, fld_size, vdf%lbound_vdf(), &
                                   vdf%ubound_vdf(), patt%send(i)%buffer, bsize(i), position )

      ! Pack b
      vdf => tiles(g_til_aid_)%t%emf%b
      fld_size = vdf%f_dim() * product(vdf%size())

      call pack_fld_msg_dlb_group( vdf%f1, vdf%f2, vdf%f3, fld_size, vdf%lbound_vdf(), &
                                   vdf%ubound_vdf(), patt%send(i)%buffer, bsize(i), position )

#ifdef __HAS_PML__
      ! Pack vpml boundary fields if necessary
      vpml => tiles(g_til_aid_)%t%emf%bnd_con%vpml_all

      do k = 1, vpml%n_vpml
        ! Pack wall_e
        wall => vpml%wall_array_e(k)
        fld_size = product(ubound(wall)-lbound(wall)+1)

        call pack_fld_msg_dlb_group( wall%f1, wall%f2, wall%f3, fld_size, lbound(wall), &
                                     ubound(wall), patt%send(i)%buffer, bsize(i), position )

        ! Pack wall_b
        wall => vpml%wall_array_b(k)
        fld_size = product(ubound(wall)-lbound(wall)+1)

        call pack_fld_msg_dlb_group( wall%f1, wall%f2, wall%f3, fld_size, lbound(wall), &
                                     ubound(wall), patt%send(i)%buffer, bsize(i), position )
      enddo
#endif

      ! Pack each species
      spec => tiles(g_til_aid_) % t % part % species

      do k = 1, p_x_dim
        lshift(k) = spec%my_nx_p(p_lower, k)
      enddo

      do
        if (.not. associated(spec)) exit
        np = spec % num_par

        do k = 1, np
          call pack_single_particle( patt%send(i)%buffer, bsize(i), position, spec, k, lshift )
        enddo

        spec => spec % next
      enddo

    enddo

    ! Send message to node
    call mpi_isend( patt%send(i)%buffer, bsize(i), MPI_PACKED, patt%send(i)%node-1, &
                    0, comm(tiles(g_til_aid_)%t%no_co), patt%send(i)%request, ierr )

  enddo

  call freemem( bsize )
  call freemem( preface )

end subroutine isend_til_msg_dlb
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Post recv for incoming tiles for diagnostics
!-------------------------------------------------------------------------------
subroutine irecv_til_msg_diag( group_vdf, no_co, grid, recv, n_recv, n_sub, fc_msg, tot_tils )

  implicit none

  type( t_vdf ), intent(inout)               :: group_vdf
  class( t_node_conf ), intent(in)           :: no_co
  class( t_grid ), intent(in)                :: grid
  type( t_msg ), dimension(:), intent(inout) :: recv
  integer, intent(in)                        :: n_recv
  integer, intent(in)                        :: n_sub
  integer, dimension(2,n_sub), intent(in)    :: fc_msg
  integer, intent(in)                        :: tot_tils

  ! (til_id, 1:p_x_dim): first coord is field comp, last p_x_dim coords are spatial coords
  integer, dimension(1+p_x_dim) :: til_sizes, til_starts

  ! first coord is field comp, last p_x_dim coords are spatial coords
  integer, dimension(1+p_x_dim) :: group_size

  ! allocate arrays as if we were sending every tile, probably won't need all this space
  integer, dimension(tot_tils*n_sub)                        :: block_lens
  integer(kind=MPI_ADDRESS_KIND), dimension(tot_tils*n_sub) :: block_disps

  integer, dimension(2,p_max_dim) :: gc_num

  integer :: n_tils, idx, idx_til, til, ierr, fld_mpi_type, struct_type, fc_off
  integer :: i, j, k, l

  gc_num = group_vdf%gc_num_

  select type(no_co); class is (t_node_conf_tiles)
  select type(grid);  class is (t_grid_group)

    fld_mpi_type = mpi_real_type( p_k_fld )

    ! will always have 1 element in each block (each tile has own mpi datatype)
    do i = 1, tot_tils*n_sub
      block_lens(i)  = 1
      block_disps(i) = 0
    enddo

    ! Get group size (depends on which fc's we're reporting)
    group_size = (/ group_vdf%f_dim_, grid%my_nx(3,1:p_x_dim) /)
    do l = 1, p_x_dim
      group_size(l+1) = group_size(l+1) + gc_num(1,l) + gc_num(2,l)
    enddo

    ! Allow for case where group_vdf has only one fc and til_vdf has many, but is sending
    ! only one fc to group_vdf
    if ( group_vdf%f_dim_ == 1 .and. fc_msg(1,1) > 1 ) then
      ! sanity check
      if ( n_sub > 1 .or. fc_msg(1,1) /= fc_msg(2,1) ) then
        ERROR('Problem passing single field component into group_vdf')
      endif
      fc_off = 1 - fc_msg(1,1)
    else
      fc_off = 0
    endif

    ! Loop over messages to be received
    do i = 1, n_recv

      n_tils = size( recv(i)%tils )

      ! Allocate receive type array to hold committed MPI subarray types
      call alloc( recv(i)%sub_type, (/ n_sub * n_tils /) )

      ! Loop over tiles to be received from a node
      idx = 0
      do j = 1, n_tils
        til = recv(i)%tils(j)

        ! Loop over subarray types
        do k = 1, n_sub
          idx = idx+1

          ! Get tile's position in group buffer and total size (field components)
          til_starts(1) = fc_msg(1,k) - 1 + fc_off
          til_sizes(1)  = fc_msg(2,k) - fc_msg(1,k) + 1

          ! Get tile's position in group buffer and total size (spatial coords)
          idx_til = 0
          do l = 1, p_x_dim
            til_starts(l+1) = grid%nx_til( 1, idx_til + g_tgp(no_co%ti_co, til, l) ) - &
                              grid%my_nx(1,l)

            til_sizes(l+1)  = grid%nx_til( 3, idx_til + g_tgp(no_co%ti_co, til, l) ) + &
                              gc_num(1,l) + gc_num(2,l)
            
            idx_til = idx_til + no_co%ti_co%g_til(l)
          enddo

          ! Create MPI subarray type for this tile (stored in recv(i)%sub_type(idx)),
          call mpi_type_create_subarray( 1+p_x_dim, group_size, til_sizes, &
            til_starts, MPI_ORDER_FORTRAN, fld_mpi_type, recv(i)%sub_type(idx), ierr )

          call mpi_type_commit( recv(i)%sub_type(idx), ierr )
        enddo
      enddo

      ! Create MPI custom type for the array of types we just made
      call mpi_type_create_struct( n_tils*n_sub, block_lens(1:n_tils*n_sub), &
        block_disps(1:n_tils*n_sub), recv(i)%sub_type, struct_type, ierr ) 

      call mpi_type_commit( struct_type, ierr )

      ! Post receive w custom data type for all tiles coming in
      ! unpack directly into the field buffer like an absolute boss
      ! TODO change these to blocking recv's
      select case(p_x_dim)
      case (1)
        call mpi_irecv( group_vdf % f1, 1, struct_type, recv(i)%node-1, 1, &
          comm(no_co), recv(i)%request, ierr )
      case (2)
        call mpi_irecv( group_vdf % f2, 1, struct_type, recv(i)%node-1, 1, &
          comm(no_co), recv(i)%request, ierr )
      case (3)
        call mpi_irecv( group_vdf % f3, 1, struct_type, recv(i)%node-1, 1, &
          comm(no_co), recv(i)%request, ierr )
      end select

    enddo

  end select
  end select

end subroutine irecv_til_msg_diag
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Pack and post sends for outgoing tiles for diags
!-------------------------------------------------------------------------------
subroutine isend_til_msg_diag( vdf_arr, til1_id, no_co, send, n_send, n_sub, fc_msg, &
                               n_fc, gc_num )

  implicit none

  type( t_vdf_arr ), dimension(:), intent(in) :: vdf_arr
  integer, intent(in)                         :: til1_id
  class( t_node_conf ), intent(in)            :: no_co
  type( t_msg ), dimension(:), intent(inout)  :: send
  integer, intent(in)                         :: n_send
  integer, intent(in)                         :: n_sub
  integer, dimension(:,:), intent(in)         :: fc_msg
  integer, intent(in)                         :: n_fc
  integer, dimension(:,:), intent(in)         :: gc_num

  select case (p_x_dim)
  case(1)
    call isend_til_msg_diag_1d( vdf_arr, til1_id, no_co, send, n_send, n_sub, fc_msg, &
      n_fc, gc_num )
  case(2)
    call isend_til_msg_diag_2d( vdf_arr, til1_id, no_co, send, n_send, n_sub, fc_msg, &
      n_fc, gc_num )
  case(3)
    call isend_til_msg_diag_3d( vdf_arr, til1_id, no_co, send, n_send, n_sub, fc_msg, &
      n_fc, gc_num )
  case default
    ERROR('Invalid value for vdf % x_dim_')
    call abort_program( p_err_invalid )
  end select

end subroutine isend_til_msg_diag
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
subroutine isend_til_msg_diag_1d( vdf_arr, til1_id, no_co, send, n_send, n_sub, fc_msg, &
                                  n_fc, gc_num )

  implicit none

  type( t_vdf_arr ), dimension(:), intent(in) :: vdf_arr
  integer, intent(in)                         :: til1_id
  class( t_node_conf ), intent(in)            :: no_co
  type( t_msg ), dimension(:), intent(inout)  :: send
  integer, intent(in)                         :: n_send
  integer, intent(in)                         :: n_sub
  integer, dimension(:,:), intent(in)         :: fc_msg
  integer, intent(in)                         :: n_fc
  integer, dimension(:,:), intent(in)         :: gc_num

  ! allocate arrays as if we were sending every tile, probably won't need all this space
  integer, dimension(size(vdf_arr))                        :: block_lens
  integer(kind=MPI_ADDRESS_KIND), dimension(size(vdf_arr)) :: block_disps

  real(p_k_fld), dimension(:), pointer     :: buffer
  real(p_k_fld), dimension(:,:), pointer   :: f1

  integer :: til_sizes, til_starts, ierr, struct_type
  integer :: til, n_tils, idx, bsize, fld_mpi_type
  integer :: i, j, k, fc, bidx
  integer :: ii, ii0, ii1

  fld_mpi_type = mpi_real_type( p_k_fld )

  ! will always have 1 element in each block (each tile has own mpi datatype)
  do i = 1, size(vdf_arr)
    block_lens(i)  = 1
    block_disps(i) = 0
  enddo

  ! Loop over messages to be sent
  do i = 1, n_send

    buffer => null()

    ! reset buffer index
    bidx = 1

    n_tils = size( send(i)%tils )

    ! Get size of comm buffer
    bsize = 0
    do j = 1, n_tils
      til = send(i)%tils(j) - til1_id + 1
      bsize = bsize + (vdf_arr(til) % v % nx_(1) + gc_num(1,1) + gc_num(2,1))
    enddo
    bsize = bsize * n_fc

    ! Allocate comm_buffer
    call alloc( send(i)%vdf_diag_buffer, (/bsize/) )
    buffer => send(i)%vdf_diag_buffer

    ! Allocate receive type array to hold committed MPI subarray types
    call alloc( send(i)%sub_type, (/ n_tils /) )

    ! Loop over tiles to be sent to a node and pack them into buffer
    idx = 0
    do j = 1, n_tils
      idx = idx+1
      til = send(i)%tils(j) - til1_id + 1

      til_starts = bidx - 1
      til_sizes  = (vdf_arr(til) % v % nx_(1) + gc_num(1,1) + gc_num(2,1)) * &
                   n_fc

      f1 => vdf_arr(til) % v % f1

      do k = 1, n_sub

        ii0 = 1 - gc_num(1,1)
        ii1 = vdf_arr(til) % v % nx_(1) + gc_num(2,1)
        do ii = ii0, ii1
          do fc = fc_msg(1,k), fc_msg(2,k)
            buffer(bidx) = f1(fc, ii)
            bidx = bidx+1
          enddo
        enddo

      enddo

      ! Create MPI subarray type for this tile (stored in send(i)%sub_type(j)),
      call mpi_type_create_subarray( 1, (/bsize/), (/til_sizes/), (/til_starts/), &
        MPI_ORDER_FORTRAN, fld_mpi_type, send(i)%sub_type(idx), ierr )

      call mpi_type_commit( send(i)%sub_type(idx), ierr )

    enddo

    ! Create MPI custom type for the array of types we just made
    call mpi_type_create_struct( n_tils, block_lens(1:n_tils), &
      block_disps(1:n_tils), send(i)%sub_type, struct_type, ierr ) 

    call mpi_type_commit( struct_type, ierr )
    
    ! Post send on comm buffer with custom mpi type
    call mpi_isend( buffer, 1, struct_type, send(i)%node-1, 1, &
      comm(no_co), send(i)%request, ierr )

  enddo

end subroutine isend_til_msg_diag_1d
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
subroutine isend_til_msg_diag_2d( vdf_arr, til1_id, no_co, send, n_send, n_sub, fc_msg, &
                                  n_fc, gc_num )

  implicit none

  type( t_vdf_arr ), dimension(:), intent(in) :: vdf_arr
  integer, intent(in)                         :: til1_id
  class( t_node_conf ), intent(in)            :: no_co
  type( t_msg ), dimension(:), intent(inout)  :: send
  integer, intent(in)                         :: n_send
  integer, intent(in)                         :: n_sub
  integer, dimension(:,:), intent(in)         :: fc_msg
  integer, intent(in)                         :: n_fc
  integer, dimension(:,:), intent(in)         :: gc_num

  ! allocate arrays as if we were sending every tile, probably won't need all this space
  integer, dimension(size(vdf_arr))                        :: block_lens
  integer(kind=MPI_ADDRESS_KIND), dimension(size(vdf_arr)) :: block_disps

  real(p_k_fld), dimension(:), pointer       :: buffer
  real(p_k_fld), dimension(:,:,:), pointer   :: f2

  integer :: til_sizes, til_starts, ierr, struct_type
  integer :: til, n_tils, idx, bsize, fld_mpi_type
  integer :: i, j, k, fc, bidx
  integer :: ii, ii0, ii1, jj, jj0, jj1

  fld_mpi_type = mpi_real_type( p_k_fld )

  ! will always have 1 element in each block (each tile has own mpi datatype)
  do i = 1, size(vdf_arr)
    block_lens(i)  = 1
    block_disps(i) = 0
  enddo

  ! Loop over messages to be sent
  do i = 1, n_send

    buffer => null()

    ! reset buffer index
    bidx = 1

    n_tils = size( send(i)%tils )

    ! Get size of comm buffer
    bsize = 0
    do j = 1, n_tils
      til = send(i)%tils(j) - til1_id + 1
      bsize = bsize + (vdf_arr(til) % v % nx_(1) + gc_num(1,1) + gc_num(2,1)) * &
                      (vdf_arr(til) % v % nx_(2) + gc_num(1,2) + gc_num(2,2))
    enddo
    bsize = bsize * n_fc

    ! Allocate comm_buffer
    call alloc( send(i)%vdf_diag_buffer, (/bsize/) )
    buffer => send(i)%vdf_diag_buffer

    ! Allocate receive type array to hold committed MPI subarray types
    call alloc( send(i)%sub_type, (/ n_tils /) )

    ! Loop over tiles to be sent to a node and pack them into buffer
    idx = 0
    do j = 1, n_tils
      idx = idx+1
      til = send(i)%tils(j) - til1_id + 1

      til_starts = bidx - 1
      til_sizes  = (vdf_arr(til) % v % nx_(1) + gc_num(1,1) + gc_num(2,1)) * &
                   (vdf_arr(til) % v % nx_(2) + gc_num(1,2) + gc_num(2,2)) * n_fc

      f2 => vdf_arr(til) % v % f2

      do k = 1, n_sub

        jj0 = 1 - gc_num(1,2)
        ii0 = 1 - gc_num(1,1)
        jj1 = vdf_arr(til) % v % nx_(2) + gc_num(2,2)
        ii1 = vdf_arr(til) % v % nx_(1) + gc_num(2,1)
        do jj = jj0, jj1
          do ii = ii0, ii1
            do fc = fc_msg(1,k), fc_msg(2,k)
              buffer(bidx) = f2(fc, ii, jj)
              bidx = bidx+1
            enddo
          enddo
        enddo

      enddo

      ! Create MPI subarray type for this tile (stored in send(i)%sub_type(j)),
      call mpi_type_create_subarray( 1, (/bsize/), (/til_sizes/), (/til_starts/), &
        MPI_ORDER_FORTRAN, fld_mpi_type, send(i)%sub_type(idx), ierr )

      call mpi_type_commit( send(i)%sub_type(idx), ierr )

    enddo

    ! Create MPI custom type for the array of types we just made
    call mpi_type_create_struct( n_tils, block_lens(1:n_tils), &
      block_disps(1:n_tils), send(i)%sub_type, struct_type, ierr ) 

    call mpi_type_commit( struct_type, ierr )
    
    ! Post send on comm buffer with custom mpi type
    call mpi_isend( buffer, 1, struct_type, send(i)%node-1, 1, &
      comm(no_co), send(i)%request, ierr )

  enddo

end subroutine isend_til_msg_diag_2d
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
subroutine isend_til_msg_diag_3d( vdf_arr, til1_id, no_co, send, n_send, n_sub, fc_msg, &
                                  n_fc, gc_num )

  implicit none

  type( t_vdf_arr ), dimension(:), intent(in) :: vdf_arr
  integer, intent(in)                         :: til1_id
  class( t_node_conf ), intent(in)            :: no_co
  type( t_msg ), dimension(:), intent(inout)  :: send
  integer, intent(in)                         :: n_send
  integer, intent(in)                         :: n_sub
  integer, dimension(:,:), intent(in)         :: fc_msg
  integer, intent(in)                         :: n_fc
  integer, dimension(:,:), intent(in)         :: gc_num

  ! allocate arrays as if we were sending every tile, probably won't need all this space
  integer, dimension(size(vdf_arr))                        :: block_lens
  integer(kind=MPI_ADDRESS_KIND), dimension(size(vdf_arr)) :: block_disps

  real(p_k_fld), dimension(:), pointer       :: buffer
  real(p_k_fld), dimension(:,:,:,:), pointer :: f3

  integer :: til_sizes, til_starts, ierr, struct_type
  integer :: til, n_tils, idx, bsize, fld_mpi_type
  integer :: i, j, k, fc, bidx
  integer :: ii, ii0, ii1, jj, jj0, jj1, kk, kk0, kk1

  fld_mpi_type = mpi_real_type( p_k_fld )

  ! will always have 1 element in each block (each tile has own mpi datatype)
  do i = 1, size(vdf_arr)
    block_lens(i)  = 1
    block_disps(i) = 0
  enddo

  ! Loop over messages to be sent
  do i = 1, n_send

    buffer => null()

    ! reset buffer index
    bidx = 1

    n_tils = size( send(i)%tils )

    ! Get size of comm buffer
    bsize = 0
    do j = 1, n_tils
      til = send(i)%tils(j) - til1_id + 1
      bsize = bsize + (vdf_arr(til) % v % nx_(1) + gc_num(1,1) + gc_num(2,1)) * &
                      (vdf_arr(til) % v % nx_(2) + gc_num(1,2) + gc_num(2,2)) * &
                      (vdf_arr(til) % v % nx_(3) + gc_num(1,3) + gc_num(2,3))
    enddo
    bsize = bsize * n_fc

    ! Allocate comm_buffer
    call alloc( send(i)%vdf_diag_buffer, (/bsize/) )
    buffer => send(i)%vdf_diag_buffer

    ! Allocate receive type array to hold committed MPI subarray types
    call alloc( send(i)%sub_type, (/ n_tils /) )

    ! Loop over tiles to be sent to a node and pack them into buffer
    idx = 0
    do j = 1, n_tils
      idx = idx+1
      til = send(i)%tils(j) - til1_id + 1

      til_starts = bidx - 1
      til_sizes  = (vdf_arr(til) % v % nx_(1) + gc_num(1,1) + gc_num(2,1)) * &
                   (vdf_arr(til) % v % nx_(2) + gc_num(1,2) + gc_num(2,2)) * &
                   (vdf_arr(til) % v % nx_(3) + gc_num(1,3) + gc_num(2,3)) * n_fc

      f3 => vdf_arr(til) % v % f3

      do k = 1, n_sub

        kk0 = 1 - gc_num(1,3)
        jj0 = 1 - gc_num(1,2)
        ii0 = 1 - gc_num(1,1)
        kk1 = vdf_arr(til) % v % nx_(3) + gc_num(2,3)
        jj1 = vdf_arr(til) % v % nx_(2) + gc_num(2,2)
        ii1 = vdf_arr(til) % v % nx_(1) + gc_num(2,1)
        do kk = kk0, kk1
          do jj = jj0, jj1
            do ii = ii0, ii1
              do fc = fc_msg(1,k), fc_msg(2,k)
                buffer(bidx) = f3(fc, ii, jj, kk)
                bidx = bidx+1
              enddo
            enddo
          enddo
        enddo

      enddo

      ! Create MPI subarray type for this tile (stored in send(i)%sub_type(j)),
      call mpi_type_create_subarray( 1, (/bsize/), (/til_sizes/), (/til_starts/), &
        MPI_ORDER_FORTRAN, fld_mpi_type, send(i)%sub_type(idx), ierr )

      call mpi_type_commit( send(i)%sub_type(idx), ierr )

    enddo

    ! Create MPI custom type for the array of types we just made
    call mpi_type_create_struct( n_tils, block_lens(1:n_tils), &
      block_disps(1:n_tils), send(i)%sub_type, struct_type, ierr ) 

    call mpi_type_commit( struct_type, ierr )
    
    ! Post send on comm buffer with custom mpi type
    call mpi_isend( buffer, 1, struct_type, send(i)%node-1, 1, &
      comm(no_co), send(i)%request, ierr )

  enddo

end subroutine isend_til_msg_diag_3d
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
subroutine wait_free_recv_diag( patt )

  implicit none

  type( t_msg_patt ), intent(inout) :: patt

  integer :: ierr
  integer :: i, j

  if ( patt%n_recv > 0 ) then
    ! Wait for message to be completed
    call wait_msg_array( patt%recv, patt%n_recv )
  endif

  ! Loop over messages to be received
  do i = 1, patt%n_recv

    ! Free MPI types in patt%recv(i)%recv_type
    do j = 1, size(patt%recv(i)%sub_type)
      call mpi_type_free( patt%recv(i)%sub_type(j), ierr )
    enddo

    ! Free memory at the end
    call freemem( patt%recv(i)%sub_type )

  enddo

end subroutine wait_free_recv_diag
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
subroutine wait_free_send_diag( patt )

  implicit none

  type( t_msg_patt ), intent(inout) :: patt

  integer :: ierr
  integer :: i, j

  if ( patt%n_send > 0 ) then
    ! Wait for message to be completed
    call wait_msg_array( patt%send, patt%n_send )
  endif

  ! Loop over messages to be sent
  do i = 1, patt%n_send

    ! Free MPI types in patt%send(i)%send_type
    do j = 1, size(patt%send(i)%sub_type)
      call mpi_type_free( patt%send(i)%sub_type(j), ierr )
    enddo

    ! Free memory at the end
    call freemem( patt%send(i)%sub_type )
    call freemem( patt%send(i)%vdf_diag_buffer )

  enddo

end subroutine wait_free_send_diag
!-------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Wait for a grouped vdf recv message to complete and unpack it
!-----------------------------------------------------------------------------------------
subroutine wait_irecv_vdf_group( vdf_arr_a, til1_id, msg, recv_msg, vdf_arr_b )

  implicit none

  type( t_vdf_arr ), dimension(:), intent(inout) :: vdf_arr_a
  integer, intent(in) :: til1_id
  class( t_bnd_msg ), intent(inout) :: msg
  type(t_vdf_msg), intent(inout) :: recv_msg
  type( t_vdf_arr ), dimension(:), intent(inout), optional :: vdf_arr_b

  integer :: i, tile, bnd, t_idx_msg, t_idx_vdf, vdf_pos
  integer, dimension(:), pointer :: msg_sum

  ! wait for message to complete
  call wait( recv_msg )

  ! create cumulative sum of message size array, accounting for a random order of tiles
  call alloc( msg_sum, (/msg%n_tils/) )
  msg_sum(1) = 1
  do i = 1, msg%n_tils-1
    tile = int(recv_msg%buffer(msg_sum(i)))
    bnd = 3 - int(recv_msg%buffer(msg_sum(i)+1))
    t_idx_msg = msg%map_tile_to_idx(tile,bnd)
    msg_sum(i+1) = msg_sum(i) + msg%msg_size(p_recv_msg, t_idx_msg)
  enddo

  ! unpack message data
  !$omp parallel do schedule(static) private(vdf_pos,tile,bnd,t_idx_msg,t_idx_vdf)
  do i = 1, msg%n_tils

    vdf_pos = msg_sum(i)

    ! Get tile number (convert from real)
    tile = int(recv_msg%buffer(vdf_pos))
    vdf_pos = vdf_pos + 1

    ! Get upper/lower boundary (convert from real)
    bnd = 3 - int(recv_msg%buffer(vdf_pos))
    vdf_pos = vdf_pos + 1

    t_idx_msg = msg%map_tile_to_idx(tile,bnd)
    t_idx_vdf = tile - til1_id + 1

    ! unpack data
    call unpack_group( recv_msg, vdf_arr_a(t_idx_vdf)%v, vdf_pos, &
                        msg%recv_vdf_range(:,:,t_idx_msg))
    if (present(vdf_arr_b)) then
      call unpack_group( recv_msg, vdf_arr_b(t_idx_vdf)%v, vdf_pos, &
                          msg%recv_vdf_range(:,:,t_idx_msg) )
    endif
  enddo
  !$omp end parallel do

  call freemem(msg_sum)

end subroutine wait_irecv_vdf_group
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! pack and post grouped send message for vdf boundary communication
!-----------------------------------------------------------------------------------------
subroutine isend_vdf_group( vdf_arr_a, til1_id, no_co, msg, send_vdf, vdf_arr_b )

  implicit none

  type( t_vdf_arr ), dimension(:), intent(in) :: vdf_arr_a
  integer, intent(in) :: til1_id
  class( t_node_conf ), intent(in) :: no_co
  class( t_bnd_msg ), intent(inout) :: msg
  type(t_vdf_msg), intent(inout) :: send_vdf
  type( t_vdf_arr ), dimension(:), intent(in), optional :: vdf_arr_b

  integer :: t_idx, tile_out, bnd, i, vdf_pos
  integer, dimension(:), pointer :: msg_sum

  ! update type only matters for receiving message
  ! msg_size set in check_buffer_size_vdf
  send_vdf%node = msg%node
  send_vdf%tag  = 0

  ! create cumulative sum of message size array
  call alloc(msg_sum, (/msg%n_tils/) )
  msg_sum(1) = 1
  do i = 1, msg%n_tils-1
    msg_sum(i+1) = msg_sum(i) + msg%msg_size(p_send_msg, i)
  enddo

  ! pack message data
  !$omp parallel do schedule(static) private(t_idx,bnd,tile_out,vdf_pos)
  do i = 1, msg%n_tils
    t_idx = msg%tils(1,i) - til1_id + 1
    bnd = msg%tils(2,i)
    tile_out = msg%tils(3,i)
    vdf_pos = msg_sum(i)

    call pack_group( send_vdf, vdf_arr_a(t_idx)%v, vdf_pos, msg%send_vdf_range(:,:,i), &
                      tile_out, bnd )

    if (present(vdf_arr_b)) then
      call pack_group( send_vdf, vdf_arr_b(t_idx)%v, vdf_pos, msg%send_vdf_range(:,:,i) )
    endif
  enddo
  !$omp end parallel do

  ! post send
  call isend( send_vdf, no_co )

  call freemem(msg_sum)

end subroutine isend_vdf_group
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Gather the size of a tile's msg
!-----------------------------------------------------------------------------------------
subroutine isend_size_spec_group( sim, msg )

  implicit none

  class(t_simulation_group), intent(in) :: sim
  class(t_bnd_msg), intent(inout) :: msg

  integer :: j, til_idx, bnd

  ! count size of each msg (only need to count send size, recv size fixed)
  do j = 1, msg%n_tils
    til_idx = msg%tils(1,j)
    bnd = msg%tils(2,j)

    ! set msg size - it is just the number of particles
    msg%msg_size(p_send_msg,j) = sim%tiles(til_idx)%t%bnd%bnd_cross(bnd)%nidx
  enddo

end subroutine isend_size_spec_group
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine isend_1_spec_group( sim, spec, send_spec, msg, dim )

  implicit none

  class( t_simulation_group ), intent(inout) :: sim
  type(t_spec_arr), dimension(sim%til_id(1):sim%til_id(2)), intent(inout) :: spec
  type(t_spec_msg), intent(inout)            :: send_spec
  class( t_bnd_msg ), intent(in)             :: msg
  integer, intent(in)                        :: dim

  integer, dimension(2) :: preface
  integer :: position, til_recv, til_send, npart, bnd, size_max
  integer :: j, ierr
  class(t_simulation), pointer :: til
  integer, dimension(:), pointer :: pos_scan

  ! Pack the preface
  position = 0
  preface = (/ send_spec%size2, send_spec%n_tils1 /)
  call mpi_pack( preface, 2, MPI_INTEGER, send_spec%buffer1, send_spec%size1, position, &
    mpi_comm_world, ierr )

  ! Prepare to pack particles in parallel
  size_max = send_spec%size1
  if ( send_spec%size2 > size_max ) size_max = send_spec%size2
  if ( size_max > 0 ) call alloc( pos_scan, (/ size_max /) )

  ! Pack 1st msg buffer
  if ( send_spec%size1 > 0 ) then

    ! Get start position of each tile
    pos_scan(1) = position
    do j = 1, send_spec%n_tils1-1
      pos_scan(j+1) = pos_scan(j) + 2*int_pack_size + &
                    msg%msg_size(p_send_msg,j) * send_spec%particle_size
    enddo

    ! Pack buffer
    !$omp parallel do schedule(dynamic) private(til_send, bnd, til_recv, npart, til, position)
    do j = 1, send_spec%n_tils1
      til_send = msg%tils(1,j)
      bnd      = msg%tils(2,j)
      til_recv = msg%tils(3,j)
      npart    = msg%msg_size(p_send_msg,j)
      til     => sim%tiles(til_send)%t
      position = pos_scan(j)
      call pack_group( spec( til_send )%s, til%no_co, send_spec%buffer1, &
        send_spec%size1, til%bnd%bnd_cross, til_recv, npart, position, bnd, dim )
    enddo
    !$omp end parallel do
  endif

  ! Initialize rest of communication parameters
  send_spec%comm = comm( sim%no_co )
  send_spec%node = msg%node
  send_spec%tag  = 0

  ! Post send for 1st message block
  call mpi_isend( send_spec%buffer1, send_spec%size1, MPI_PACKED, send_spec%node - 1, &
    send_spec%tag, send_spec%comm, send_spec%request, ierr )

  if (ierr/=0) then
    ERROR("MPI error")
    call abort_program( p_err_mpi )
  endif

  ! If necessary pack 2nd msg buffer
  if ( send_spec%size2 > 0 ) then

    ! Get start position of each tile
    pos_scan(1) = 0
    do j = send_spec%n_tils1+1, msg%n_tils-1
      pos_scan(j-send_spec%n_tils1+1) = pos_scan(j-send_spec%n_tils1) + 2*int_pack_size +&
                    msg%msg_size(p_send_msg,j) * send_spec%particle_size
    enddo

    ! Pack buffer
    !$omp parallel do schedule(dynamic) private(til_send, bnd, til_recv, npart, til, position)
    do j = send_spec%n_tils1+1, msg%n_tils
      til_send = msg%tils(1,j)
      bnd      = msg%tils(2,j)
      til_recv = msg%tils(3,j)
      npart    = msg%msg_size(p_send_msg,j)
      til     => sim%tiles(til_send)%t
      position = pos_scan(j-send_spec%n_tils1)
      call pack_group( spec( til_send )%s, til%no_co, send_spec%buffer2, &
        send_spec%size2, til%bnd%bnd_cross, til_recv, npart, position, bnd, dim )
    enddo
    !$omp end parallel do
  endif

  if ( size_max > 0 ) call freemem( pos_scan )

#ifdef __HAS_TRACKS__
  ! ERROR("TRACKS not implemented")
  ! call abort_program( p_err_invalid )
  ! TODO
  ! if tracking move tracks that were present and left to the missing list
  ! if ( ( spec%diag%ndump_fac_tracks > 0 ) .and. ( n_idx > 0 ) ) then
  !   call missing_particles( spec%diag%tracks, par_idx, n_idx, spec%tag )
  ! endif
#endif

end subroutine isend_1_spec_group
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! unpack incoming group of msgs for update boundary
!-----------------------------------------------------------------------------------------
subroutine unpack_spec_group( sim, spec, recv_spec, n_tils )

  implicit none

  class(t_simulation_group), intent(inout) :: sim
  type(t_spec_arr), dimension(sim%til_id(1):sim%til_id(2)), intent(inout) :: spec
  type( t_spec_msg ), intent(inout)        :: recv_spec
  integer, intent(in)                      :: n_tils

  integer :: position, ierr, i, j, til_recv, size_max, npart, my_til
  integer, dimension(:), pointer :: pos_scan, tile_nums
  logical, dimension(:), pointer :: tile_duplicate

  ! If message hasn't arrived wait for it
  call wait_spec_msg_1( recv_spec )

  ! Unpack preface
  ! skip the first int (size of second buffer) and unpack number of tiles in first msg
  position = int_pack_size
  call mpi_unpack( recv_spec%buffer1, recv_spec%size1, position, &
    recv_spec%n_tils1, 1, MPI_INTEGER, mpi_comm_world, ierr )

  ! Prepare for omp loop
  size_max = recv_spec%n_tils1
  if ( n_tils - size_max > size_max ) size_max = n_tils - size_max
  if ( size_max > 0 ) then
    call alloc( pos_scan, (/ size_max /) )
    call alloc( tile_nums, (/ size_max /) )
    allocate( tile_duplicate(size_max) )
  endif

  do j = 1, recv_spec%n_tils1
    pos_scan(j) = position
    ! Unpack tile id for determining if we need to skip duplicates
    call mpi_unpack( recv_spec%buffer1, recv_spec%size1, position, tile_nums(j), 1, &
      MPI_INTEGER, mpi_comm_world, ierr )
    ! Get number of particles for this tile
    call mpi_unpack( recv_spec%buffer1, recv_spec%size1, position, npart, 1, &
      MPI_INTEGER, mpi_comm_world, ierr )
    position = position + npart * recv_spec%particle_size
  enddo

  ! Calculate which tiles have duplicate entries, change those values to true
  tile_duplicate = .false.
  do j = 1, recv_spec%n_tils1
    if (.not. tile_duplicate(j)) then
      my_til = tile_nums(j)
      do i = j+1, recv_spec%n_tils1
        if (my_til == tile_nums(i)) then
          tile_duplicate(j) = .true.
          tile_duplicate(i) = .true.
        endif
      enddo
    endif
  enddo

  ! Unpack msgs from non-duplicate tiles in first buffer
  !$omp parallel do schedule(dynamic) private(position, til_recv, ierr)
  do j = 1, recv_spec%n_tils1

    if ( .not. tile_duplicate(j) ) then
      position = pos_scan(j)

      ! unpack g_til_aid of recipient
      call mpi_unpack( recv_spec%buffer1, recv_spec%size1, position, til_recv, 1, &
        MPI_INTEGER, mpi_comm_world, ierr )

      call unpack( spec(til_recv)%s, recv_spec%buffer1, recv_spec%size1, position, &
        sim%tiles(til_recv)%t%bnd%bnd_cross )
    endif

  enddo
  !$omp end parallel do

  ! Unpack msgs from duplicate tiles in first buffer
  ! This must be done outside an omp loop since bnd_cross is modified in unpack
  do j = 1, recv_spec%n_tils1

    if ( tile_duplicate(j) ) then
      position = pos_scan(j)

      ! unpack g_til_aid of recipient
      call mpi_unpack( recv_spec%buffer1, recv_spec%size1, position, til_recv, 1, &
        MPI_INTEGER, mpi_comm_world, ierr )

      call unpack( spec(til_recv)%s, recv_spec%buffer1, recv_spec%size1, position, &
        sim%tiles(til_recv)%t%bnd%bnd_cross )
    endif

  enddo

  if ( recv_spec%n_tils1 < n_tils ) then

    ! reset position
    position = 0

    do j = 1, n_tils-recv_spec%n_tils1
      pos_scan(j) = position
      ! Unpack tile id for determining if we need to skip duplicates
      call mpi_unpack( recv_spec%buffer2, recv_spec%size2, position, tile_nums(j), 1, &
        MPI_INTEGER, mpi_comm_world, ierr )
      ! Get number of particles for this tile
      call mpi_unpack( recv_spec%buffer2, recv_spec%size2, position, npart, 1, &
        MPI_INTEGER, mpi_comm_world, ierr )
      position = position + npart * recv_spec%particle_size
    enddo

    ! Calculate which tiles have duplicate entries, change those values to true
    tile_duplicate = .false.
    do j = 1, n_tils-recv_spec%n_tils1
      if (.not. tile_duplicate(j)) then
        my_til = tile_nums(j)
        do i = j+1, n_tils-recv_spec%n_tils1
          if (my_til == tile_nums(i)) then
            tile_duplicate(j) = .true.
            tile_duplicate(i) = .true.
          endif
        enddo
      endif
    enddo

    ! Unpack msgs from non-duplicate tiles in second buffer
    !$omp parallel do schedule(dynamic) private(position, til_recv, ierr)
    do j = 1, n_tils-recv_spec%n_tils1

      if ( .not. tile_duplicate(j) ) then
        position = pos_scan(j)

        ! unpack g_til_aid of recipient
        call mpi_unpack( recv_spec%buffer2, recv_spec%size2, position, til_recv, 1, &
          MPI_INTEGER, mpi_comm_world, ierr )

        call unpack( spec(til_recv)%s, recv_spec%buffer2, recv_spec%size2, position, &
          sim%tiles(til_recv)%t%bnd%bnd_cross )
      endif

    enddo
    !$omp end parallel do

    ! Unpack msgs from duplicate tiles in second buffer
    do j = 1, n_tils-recv_spec%n_tils1

      if ( tile_duplicate(j) ) then
        position = pos_scan(j)

        ! unpack g_til_aid of recipient
        call mpi_unpack( recv_spec%buffer2, recv_spec%size2, position, til_recv, 1, &
          MPI_INTEGER, mpi_comm_world, ierr )

        call unpack( spec(til_recv)%s, recv_spec%buffer2, recv_spec%size2, position, &
          sim%tiles(til_recv)%t%bnd%bnd_cross )
      endif

    enddo

  endif

  if ( size_max > 0 ) then
    call freemem( pos_scan )
    call freemem( tile_nums )
    deallocate( tile_duplicate )
  endif

end subroutine
!-----------------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
subroutine wait_send_msg( patt )

  implicit none

  type( t_msg_patt ), intent(inout) :: patt

  call wait_msg_array( patt%send, patt%n_send )

end subroutine wait_send_msg
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
subroutine wait_til_msg( msg )

  implicit none

  type(t_msg), intent(inout) :: msg

  integer, dimension( MPI_STATUS_SIZE ) :: status
  integer :: ierr

  if ( msg%request /= MPI_REQUEST_NULL ) then
    call MPI_WAIT( msg%request, status, ierr )
    if ( ierr /= 0 ) then
      ERROR("MPI error")
      call abort_program( p_err_mpi )
    endif
  endif

end subroutine wait_til_msg
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
subroutine wait_msg_array( msg, n_msg )

  implicit none

  type(t_msg), dimension(:), pointer :: msg
  integer, intent(in)                :: n_msg

  integer :: ierr
  integer, dimension(:,:), pointer :: status
  integer, dimension(n_msg) :: request

  integer :: i,j

  j = 0
  do i = 1, n_msg
    if (msg(i)%request /= MPI_REQUEST_NULL) then
      j = j+1
      request(j) = msg(i)%request
    endif
  enddo

  if (j > 0) then
    call alloc ( status, (/MPI_STATUS_SIZE, j/) )

    call MPI_WAITALL(j, request, status, ierr)
    if (ierr/=0) then
      ERROR("MPI error")
      call abort_program( p_err_mpi )
    endif

    call freemem( status )

  endif

end subroutine wait_msg_array
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! unpack tile messages
!-------------------------------------------------------------------------------
subroutine unpack_tiles_dlb( patt, bsize, tiles )

  implicit none

  type( t_msg_patt ), intent(inout)                           :: patt
  integer, dimension(:), pointer, intent(inout)               :: bsize
  class( t_til_container ), dimension(:), pointer, intent(in) :: tiles

  integer :: position, til, num_spec, fld_mpi_type, part_mpi_type, vdf_size
  integer :: i, j, k, l, ierr, np, bseg_size
  integer, dimension(:), pointer :: preface
  class(t_species), pointer :: spec
  type(t_vdf), pointer :: vdf
  type(t_vpml), pointer :: vpml
  type(t_wall), pointer :: wall

  num_spec = tiles(lbound(tiles,1)) % t % part % num_species
  call alloc( preface, (/ num_spec+1 /) )

  fld_mpi_type = mpi_real_type( p_k_fld )
  part_mpi_type = mpi_real_type( p_k_part )

  do i = 1, patt%n_recv

    ! Wait for message to be received
    call wait_til_msg( patt%recv(i) )

    ! Unpack tiles
    position = 0

    ! Pack all tiles going to the same node serially, one after the other
    do j = 1, size( patt%recv(i)%tils )

      ! First unpack preface, containing g_til_aid and num_par for each species
      call mpi_unpack( patt%recv(i)%buffer, bsize(i), position, preface, num_spec+1,  &
        MPI_INTEGER, mpi_comm_world, ierr )

      til = preface(1)

      ! Unpack e
      vdf => tiles(til) % t % emf % e
      vdf_size = vdf%f_dim() * product(vdf%size())

      call unpack_fld_msg_dlb_group( vdf%f1, vdf%f2, vdf%f3, vdf_size, vdf%lbound_vdf(), &
                                     vdf%ubound_vdf(), patt%recv(i)%buffer, bsize(i), position )

      ! Unpack b
      vdf => tiles(til) % t % emf % b
      vdf_size = vdf%f_dim() * product(vdf%size())

      call unpack_fld_msg_dlb_group( vdf%f1, vdf%f2, vdf%f3, vdf_size, vdf%lbound_vdf(), &
                                     vdf%ubound_vdf(), patt%recv(i)%buffer, bsize(i), position )

#ifdef __HAS_PML__
      ! Unpack vpml boundary fields if necessary
      vpml => tiles(til)%t%emf%bnd_con%vpml_all

      do k = 1, vpml%n_vpml
        ! Unpack wall_e
        wall => vpml%wall_array_e(k)
        vdf_size = product(ubound(wall)-lbound(wall)+1)

        call unpack_fld_msg_dlb_group( wall%f1, wall%f2, wall%f3, vdf_size, lbound(wall), &
                                       ubound(wall), patt%recv(i)%buffer, bsize(i), position )

        ! Unpack wall_b
        wall => vpml%wall_array_b(k)
        vdf_size = product(ubound(wall)-lbound(wall)+1)

        call unpack_fld_msg_dlb_group( wall%f1, wall%f2, wall%f3, vdf_size, lbound(wall), &
                                       ubound(wall), patt%recv(i)%buffer, bsize(i), position )
      enddo
#endif

      ! Unpack each species
      spec => tiles(til) % t % part % species
      k = 1
      do
        if (.not. associated(spec)) exit
        k = k+1
        np = preface(k)

        ! allocate needed buffers
        ! sanity check
        if (associated(spec%x)) then
          ERROR(" --in unpack_tiles, species particle arrays already allocated",til)
          call abort_program( p_err_invalid )
        endif

        ! initialize particle buffers with some extra space if needed
        if ( np >= spec%num_par_max ) then
          bseg_size = np + p_spec_buf_block_t
        else
          bseg_size = spec%num_par_max
        endif
        call spec%init_buffer( bseg_size )

        do l = 1, np
          call unpack_single_particle( patt%recv(i)%buffer, bsize(i), position, spec, l )
        enddo

        spec % num_par = np

        spec => spec % next
      enddo

    enddo

  enddo

  call freemem(preface)
  call freemem(bsize)

end subroutine unpack_tiles_dlb
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! msg for each tile will be prefaceded with a number of integers:
! g_til_aid and number of particles in each species
!-------------------------------------------------------------------------------
subroutine get_size_preface( num_species, bsize )

  implicit none

  integer, intent(in)    :: num_species
  integer, intent(inout) :: bsize

  integer :: preface_size

  ! int_pack_size comes from os-spec-comm
  preface_size = int_pack_size + num_species*int_pack_size

  call check_bsize_overflow( bsize, preface_size )

  bsize = bsize + preface_size

end subroutine get_size_preface
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! get emf contribution to size of message
!-------------------------------------------------------------------------------
subroutine get_size_emf( emf, bsize )

  implicit none

  class( t_emf ), intent(in) :: emf
  integer, intent(inout)     :: bsize

  type( t_vpml ) :: vpml_all
  integer :: i
  integer :: emf_size

  vpml_all = emf%bnd_con%vpml_all

  ! e and b contributions
  emf_size = fld_pack_size()*product(emf%e%ubound()-emf%e%lbound()+1)
  call check_bsize_overflow( bsize, emf_size )
  bsize = bsize + emf_size

  emf_size = fld_pack_size()*product(emf%b%ubound()-emf%b%lbound()+1)
  call check_bsize_overflow( bsize, emf_size )
  bsize = bsize + emf_size

  ! vpml contribution to size of message
  ! only add e and b field data in wall arrays, rest of member data gets initialized
  do i = 1, vpml_all%n_vpml
    emf_size = fld_pack_size()*product(ubound(vpml_all%wall_array_e(i))-&
                                            lbound(vpml_all%wall_array_e(i))+1)
    call check_bsize_overflow( bsize, emf_size )
    bsize = bsize + emf_size

    emf_size = fld_pack_size()*product(ubound(vpml_all%wall_array_b(i))-&
                                            lbound(vpml_all%wall_array_b(i))+1)
    call check_bsize_overflow( bsize, emf_size )
    bsize = bsize + emf_size
  enddo

end subroutine get_size_emf
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! get particle contribution to size of message
!-------------------------------------------------------------------------------
subroutine get_size_part( part, bsize )

  implicit none

  class(t_particles), pointer :: part
  integer, intent(inout)      :: bsize

  integer :: spec_size
  class(t_species), pointer :: species

  species => part % species

  do
    if (.not. associated(species)) exit
    spec_size = species%num_par * particle_pack_size(species)
    call check_bsize_overflow( bsize, spec_size )
    bsize   =  bsize + spec_size
    species => species % next
  enddo

end subroutine get_size_part
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
subroutine check_bsize_overflow( bsize, pack_size )

  implicit none

  integer, intent(in) :: bsize
  integer, intent(in) :: pack_size

  integer :: bsize_max

  bsize_max = huge(bsize)

  if ( bsize > bsize_max - pack_size ) then
    SCR_MPINODE('Invalid value for tile buffer size. Integer overflow.')
    call abort_program( p_err_invalid )
  endif

end subroutine check_bsize_overflow
!-------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! First half of deposit species charge on a grid for tiled run
!-----------------------------------------------------------------------------------------
subroutine deposit_quant_1_tiles( spec, grid, charge, quant )

  implicit none

  class( t_species ), intent(in) :: spec
  class( t_grid ), intent(in) :: grid
  type( t_vdf ), intent(inout) :: charge
  integer, intent(in) :: quant

  ! local variables
  integer, parameter :: p_part_block = 4096
  real(p_k_part), dimension(p_part_block) :: q
  integer :: i1, i2

  integer, dimension (2, p_x_dim) :: gc_num
  integer :: lquant

  ! create vdf

  ! The number of guard cells required depends on the interpolation level
  gc_num(p_lower,:) = spec%interpolation
  gc_num(p_upper,:) = spec%interpolation + 1

  ! note that this automatically sets the vdf value to 0.0
  call charge% new( p_x_dim, 1, grid%my_nx(3,:), gc_num , spec%dx, .true.)

  ! deposit given quantity
  if ( quant == p_norm ) then
    lquant = p_charge
  else
    lquant = quant
  endif

  do i1 = 1, spec%num_par, p_part_block
    i2 = i1 + p_part_block - 1
    if ( i2 > spec%num_par ) i2 = spec%num_par
    call spec % get_quant( i1, i2, lquant, q )
    call spec % deposit_density( charge, i1, i2, q  )
  enddo

end subroutine deposit_quant_1_tiles
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Second half of deposit species charge on a grid for tiled run
!-----------------------------------------------------------------------------------------
subroutine deposit_quant_2_tiles( spec, grid, charge, quant )

  implicit none

  class( t_species ), intent(in) :: spec
  class( t_grid ), intent(in) :: grid
  type( t_vdf ), intent(inout) :: charge
  integer, intent(in) :: quant

  ! local variables
  integer, parameter :: p_part_block = 4096
  integer :: i1, i2, i3


  ! normalize charge for cylindrical coordinates
  if ( coordinates( grid ) == p_cylindrical_b ) then
     call norm_charge_cyl( charge, grid%my_nx( 1, p_r_dim ), &
                           real( spec%dx( p_r_dim ), p_k_fld ) )
  endif

  ! get p_norm quantity requires special treatment
  if ( quant == p_norm ) then
    select case( p_x_dim )
    case (1)
      do i1 = lbound( charge%f1, 2 ), ubound( charge%f1, 2 )
        if ( charge%f1( 1, i1 ) /= 0.0 ) then
          charge%f1( 1, i1 ) = abs( 1.0 / charge%f1( 1, i1 ) )
        else
          charge%f1( 1, i1 ) = 1.0
        endif
      enddo

    case (2)
      do i2 = lbound( charge%f2, 3 ), ubound( charge%f2, 3 )
        do i1 = lbound( charge%f2, 2 ), ubound( charge%f2, 2 )
          if ( charge%f2( 1, i1, i2 ) /= 0.0 ) then
            charge%f2( 1, i1, i2 ) = abs( 1.0 / charge%f2( 1, i1, i2 ) )
          else
            charge%f2( 1, i1, i2 ) = 1.0
          endif
        enddo
      enddo

    case (3)
      do i3 = lbound( charge%f3, 4 ), ubound( charge%f3, 4 )
        do i2 = lbound( charge%f3, 3 ), ubound( charge%f3, 3 )
          do i1 = lbound( charge%f3, 2 ), ubound( charge%f3, 2 )
            if ( charge%f3( 1, i1, i2, i3 ) /= 0.0 ) then
              charge%f3( 1, i1, i2, i3 ) = abs( 1.0 / charge%f3( 1, i1, i2, i3 ) )
            else
              charge%f3( 1, i1, i2, i3 ) = 1.0
            endif
          enddo
        enddo
      enddo
    end select
  endif

end subroutine deposit_quant_2_tiles
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Deposit charge density of particles. Used for charge conservation and diagnostics
!-----------------------------------------------------------------------------------------
subroutine deposit_charge_part_1( this, g_space, grid, charge )

  implicit none

  class( t_particles ), intent(in) :: this
  class( t_grid ), intent(in) :: grid
  type( t_space ), intent(in) :: g_space
  type( t_vdf ), intent(inout) :: charge

  integer, dimension (2, p_max_dim) :: lgc_num
  real(p_double), dimension(p_x_dim) :: ldx
  type( t_vdf ) :: tmp_charge
  class(t_species), pointer :: species

  ! Create charge vdf if required
  if ( charge%x_dim_ < 1 ) then
    ldx = extend(g_space)/ grid % g_nx(1 : p_x_dim)
    lgc_num(1,:) = this%interpolation
    lgc_num(2,:) = this%interpolation+1
    call charge % new( grid%x_dim, 1, grid%my_nx(3,:), lgc_num , ldx, .false.)
  endif

  ! set array for charge density to zero before depositing charge for each species
  call charge % zero()

  if ( this%num_species > 0 ) then

    species => this % species
    call deposit_rho( species, charge )

    if ( this%num_species > 1 ) then

      if ( this%charge%low_roundoff ) then

        ! This minimizes roundoff in charge deposition by depositing 1 species at a time
        call tmp_charge % new( charge )

        do
          species => species % next
          if (.not. associated(species)) exit
          call tmp_charge % zero()
          call deposit_rho( species, tmp_charge )
          call add( charge, tmp_charge )
        enddo

        call tmp_charge % cleanup()

      else

        do
          species => species % next
          if (.not. associated(species)) exit
          call deposit_rho( species, charge )
        enddo

      endif

    endif
  endif

end subroutine deposit_charge_part_1
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Deposit charge density of particles. Used for charge conservation and diagnostics
!-----------------------------------------------------------------------------------------
subroutine deposit_charge_part_2( this, grid, charge )

  implicit none

  class( t_particles ), intent(in) :: this
  class( t_grid ),intent(in) :: grid
  type( t_vdf ), intent(inout) :: charge

  if ( this%num_species > 0 ) then
    ! normalize charge for cylindrical coordinates
    if ( grid%coordinates == p_cylindrical_b ) then
      call norm_charge_cyl( charge, grid%my_nx( 1, p_r_dim ), real( charge%dx(p_r_dim), p_k_fld ) )
    endif
  endif

end subroutine deposit_charge_part_2
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! send/recv e and b ghost cells between two tiles
!-----------------------------------------------------------------------------------------
subroutine loc_fill_vdf_tiles( vdfs, update_type, dim, bnd, nx_move )

  implicit none

  type( t_vdf_arr ), dimension(:,:), intent(inout) :: vdfs
  integer, intent(in) :: update_type, dim, bnd, nx_move

  integer, dimension(2,p_max_dim) :: l_nr, l_ns, u_nr, u_ns
  integer :: i, j, i1, i2, i3, msg_size
  integer, dimension(p_max_dim) :: off
  real(p_k_fld), dimension(:,:), pointer :: lf1, uf1
  real(p_k_fld), dimension(:,:,:), pointer :: lf2, uf2
  real(p_k_fld), dimension(:,:,:,:), pointer :: lf3, uf3

  ! get ranges for sending/receiving data
  call get_recv_range( vdfs(p_lower,1)%v, msg_size, l_nr, bnd, update_type, dim, nx_move )
  call get_send_range( vdfs(p_lower,1)%v, msg_size, l_ns, bnd, update_type, dim, nx_move )
  call get_recv_range( vdfs(p_upper,1)%v, msg_size, u_nr, 3-bnd, update_type, dim, nx_move )
  call get_send_range( vdfs(p_upper,1)%v, msg_size, u_ns, 3-bnd, update_type, dim, nx_move )

  if (vdfs(p_lower,1)%v%f_dim_ == 3) then

    select case( update_type )
    case( p_vdf_add )

      do j = 1, size(vdfs,2)

        ! add data from upper tile to lower
        ! compute offset once
        off = u_ns(p_lower,:) - l_nr(p_lower,:)
        select case ( vdfs(p_lower,j)%v%x_dim_ )
        case (1)
          lf1 => vdfs(p_lower,j)%v%f1; uf1 => vdfs(p_upper,j)%v%f1
          do i1 = l_nr( p_lower, 1 ), l_nr( p_upper, 1 )
            lf1( 1, i1 ) = lf1( 1, i1 ) + &
                                uf1( 1, i1+off(1) )
            lf1( 2, i1 ) = lf1( 2, i1 ) + &
                                uf1( 2, i1+off(1) )
            lf1( 3, i1 ) = lf1( 3, i1 ) + &
                                uf1( 3, i1+off(1) )
          enddo

        case (2)
          lf2 => vdfs(p_lower,j)%v%f2; uf2 => vdfs(p_upper,j)%v%f2
          do i2 = l_nr( p_lower, 2 ), l_nr( p_upper, 2 )
            do i1 = l_nr( p_lower, 1 ), l_nr( p_upper, 1 )
              lf2( 1, i1, i2 ) = lf2( 1, i1, i2 ) + &
                                      uf2( 1, i1+off(1), i2+off(2) )
              lf2( 2, i1, i2 ) = lf2( 2, i1, i2 ) + &
                                      uf2( 2, i1+off(1), i2+off(2) )
              lf2( 3, i1, i2 ) = lf2( 3, i1, i2 ) + &
                                      uf2( 3, i1+off(1), i2+off(2) )
            enddo
          enddo

        case (3)
          lf3 => vdfs(p_lower,j)%v%f3; uf3 => vdfs(p_upper,j)%v%f3
          do i3 = l_nr( p_lower, 3 ), l_nr( p_upper, 3 )
            do i2 = l_nr( p_lower, 2 ), l_nr( p_upper, 2 )
              do i1 = l_nr( p_lower, 1 ), l_nr( p_upper, 1 )
                lf3( 1, i1, i2, i3 ) = lf3( 1, i1, i2, i3 ) + &
                                            uf3( 1, i1+off(1), i2+off(2), i3+off(3) )
                lf3( 2, i1, i2, i3 ) = lf3( 2, i1, i2, i3 ) + &
                                            uf3( 2, i1+off(1), i2+off(2), i3+off(3) )
                lf3( 3, i1, i2, i3 ) = lf3( 3, i1, i2, i3 ) + &
                                            uf3( 3, i1+off(1), i2+off(2), i3+off(3) )
              enddo
            enddo
          enddo
        end select

        ! copy data from lower tile to upper
        ! compute offset once
        off = l_ns(p_lower,:) - u_nr(p_lower,:)
        select case ( vdfs(p_upper,j)%v%x_dim_ )
        case (1)
          do i1 = u_nr( p_lower, 1 ), u_nr( p_upper, 1 )
            uf1( 1, i1 ) = lf1( 1, i1+off(1) )
            uf1( 2, i1 ) = lf1( 2, i1+off(1) )
            uf1( 3, i1 ) = lf1( 3, i1+off(1) )
          enddo

        case (2)
          do i2 = u_nr( p_lower, 2 ), u_nr( p_upper, 2 )
            do i1 = u_nr( p_lower, 1 ), u_nr( p_upper, 1 )
              uf2( 1, i1, i2 ) = lf2( 1, i1+off(1), i2+off(2) )
              uf2( 2, i1, i2 ) = lf2( 2, i1+off(1), i2+off(2) )
              uf2( 3, i1, i2 ) = lf2( 3, i1+off(1), i2+off(2) )
            enddo
          enddo

        case (3)
          do i3 = u_nr( p_lower, 3 ), u_nr( p_upper, 3 )
            do i2 = u_nr( p_lower, 2 ), u_nr( p_upper, 2 )
              do i1 = u_nr( p_lower, 1 ), u_nr( p_upper, 1 )
                uf3( 1, i1, i2, i3 ) = lf3( 1, i1+off(1), i2+off(2), i3+off(3) )
                uf3( 2, i1, i2, i3 ) = lf3( 2, i1+off(1), i2+off(2), i3+off(3) )
                uf3( 3, i1, i2, i3 ) = lf3( 3, i1+off(1), i2+off(2), i3+off(3) )
              enddo
            enddo
          enddo
        end select

      enddo

    case( p_vdf_replace )

      do j = 1, size(vdfs,2)

        ! copy data from upper tile to lower
        ! compute offset once
        off = u_ns(p_lower,:) - l_nr(p_lower,:)
        select case ( vdfs(p_lower,j)%v%x_dim_ )
        case (1)
          lf1 => vdfs(p_lower,j)%v%f1; uf1 => vdfs(p_upper,j)%v%f1
          do i1 = l_nr( p_lower, 1 ), l_nr( p_upper, 1 )
            lf1( 1, i1 ) = uf1( 1, i1+off(1) )
            lf1( 2, i1 ) = uf1( 2, i1+off(1) )
            lf1( 3, i1 ) = uf1( 3, i1+off(1) )
          enddo

        case (2)
          lf2 => vdfs(p_lower,j)%v%f2; uf2 => vdfs(p_upper,j)%v%f2
          do i2 = l_nr( p_lower, 2 ), l_nr( p_upper, 2 )
            do i1 = l_nr( p_lower, 1 ), l_nr( p_upper, 1 )
              lf2( 1, i1, i2 ) = uf2( 1, i1+off(1), i2+off(2) )
              lf2( 2, i1, i2 ) = uf2( 2, i1+off(1), i2+off(2) )
              lf2( 3, i1, i2 ) = uf2( 3, i1+off(1), i2+off(2) )
            enddo
          enddo

        case (3)
          lf3 => vdfs(p_lower,j)%v%f3; uf3 => vdfs(p_upper,j)%v%f3
          do i3 = l_nr( p_lower, 3 ), l_nr( p_upper, 3 )
            do i2 = l_nr( p_lower, 2 ), l_nr( p_upper, 2 )
              do i1 = l_nr( p_lower, 1 ), l_nr( p_upper, 1 )
                lf3( 1, i1, i2, i3 ) = uf3( 1, i1+off(1), i2+off(2), i3+off(3) )
                lf3( 2, i1, i2, i3 ) = uf3( 2, i1+off(1), i2+off(2), i3+off(3) )
                lf3( 3, i1, i2, i3 ) = uf3( 3, i1+off(1), i2+off(2), i3+off(3) )
              enddo
            enddo
          enddo
        end select

        ! copy data from lower tile to upper
        ! compute offset once
        off = l_ns(p_lower,:) - u_nr(p_lower,:)
        select case ( vdfs(p_upper,j)%v%x_dim_ )
        case (1)
          do i1 = u_nr( p_lower, 1 ), u_nr( p_upper, 1 )
            uf1( 1, i1 ) = lf1( 1, i1+off(1) )
            uf1( 2, i1 ) = lf1( 2, i1+off(1) )
            uf1( 3, i1 ) = lf1( 3, i1+off(1) )
          enddo

        case (2)
          do i2 = u_nr( p_lower, 2 ), u_nr( p_upper, 2 )
            do i1 = u_nr( p_lower, 1 ), u_nr( p_upper, 1 )
              uf2( 1, i1, i2 ) = lf2( 1, i1+off(1), i2+off(2) )
              uf2( 2, i1, i2 ) = lf2( 2, i1+off(1), i2+off(2) )
              uf2( 3, i1, i2 ) = lf2( 3, i1+off(1), i2+off(2) )
            enddo
          enddo

        case (3)
          do i3 = u_nr( p_lower, 3 ), u_nr( p_upper, 3 )
            do i2 = u_nr( p_lower, 2 ), u_nr( p_upper, 2 )
              do i1 = u_nr( p_lower, 1 ), u_nr( p_upper, 1 )
                uf3( 1, i1, i2, i3 ) = lf3( 1, i1+off(1), i2+off(2), i3+off(3) )
                uf3( 2, i1, i2, i3 ) = lf3( 2, i1+off(1), i2+off(2), i3+off(3) )
                uf3( 3, i1, i2, i3 ) = lf3( 3, i1+off(1), i2+off(2), i3+off(3) )
              enddo
            enddo
          enddo
        end select

      enddo

    end select

  else

    select case( update_type )
    case( p_vdf_add )

      do j = 1, size(vdfs,2)

        ! add data from upper tile to lower
        ! compute offset once
        off = u_ns(p_lower,:) - l_nr(p_lower,:)
        select case ( vdfs(p_lower,j)%v%x_dim_ )
        case (1)
          lf1 => vdfs(p_lower,j)%v%f1; uf1 => vdfs(p_upper,j)%v%f1
          do i1 = l_nr( p_lower, 1 ), l_nr( p_upper, 1 )
            do i = 1, vdfs(p_lower,j)%v%f_dim_
              lf1( i, i1 ) = lf1( i, i1 ) + &
                                  uf1( i, i1+off(1) )
            enddo
          enddo

        case (2)
          lf2 => vdfs(p_lower,j)%v%f2; uf2 => vdfs(p_upper,j)%v%f2
          do i2 = l_nr( p_lower, 2 ), l_nr( p_upper, 2 )
            do i1 = l_nr( p_lower, 1 ), l_nr( p_upper, 1 )
              do i = 1, vdfs(p_lower,j)%v%f_dim_
                lf2( i, i1, i2 ) = lf2( i, i1, i2 ) + &
                                        uf2( i, i1+off(1), i2+off(2) )
              enddo
            enddo
          enddo

        case (3)
          lf3 => vdfs(p_lower,j)%v%f3; uf3 => vdfs(p_upper,j)%v%f3
          do i3 = l_nr( p_lower, 3 ), l_nr( p_upper, 3 )
            do i2 = l_nr( p_lower, 2 ), l_nr( p_upper, 2 )
              do i1 = l_nr( p_lower, 1 ), l_nr( p_upper, 1 )
                do i = 1, vdfs(p_lower,j)%v%f_dim_
                  lf3( i, i1, i2, i3 ) = lf3( i, i1, i2, i3 ) + &
                                              uf3( i, i1+off(1), i2+off(2), i3+off(3) )
                enddo
              enddo
            enddo
          enddo
        end select

        ! copy data from lower tile to upper
        ! compute offset once
        off = l_ns(p_lower,:) - u_nr(p_lower,:)
        select case ( vdfs(p_upper,j)%v%x_dim_ )
        case (1)
          do i1 = u_nr( p_lower, 1 ), u_nr( p_upper, 1 )
            do i = 1, vdfs(p_upper,j)%v%f_dim_
              uf1( i, i1 ) = lf1( i, i1+off(1) )
            enddo
          enddo

        case (2)
          do i2 = u_nr( p_lower, 2 ), u_nr( p_upper, 2 )
            do i1 = u_nr( p_lower, 1 ), u_nr( p_upper, 1 )
              do i = 1, vdfs(p_upper,j)%v%f_dim_
                uf2( i, i1, i2 ) = lf2( i, i1+off(1), i2+off(2) )
              enddo
            enddo
          enddo

        case (3)
          do i3 = u_nr( p_lower, 3 ), u_nr( p_upper, 3 )
            do i2 = u_nr( p_lower, 2 ), u_nr( p_upper, 2 )
              do i1 = u_nr( p_lower, 1 ), u_nr( p_upper, 1 )
                do i = 1, vdfs(p_upper,j)%v%f_dim_
                  uf3( i, i1, i2, i3 ) = lf3( i, i1+off(1), i2+off(2), i3+off(3) )
                enddo
              enddo
            enddo
          enddo
        end select

      enddo

    case( p_vdf_replace )

      do j = 1, size(vdfs,2)

        ! copy data from upper tile to lower
        ! compute offset once
        off = u_ns(p_lower,:) - l_nr(p_lower,:)
        select case ( vdfs(p_lower,j)%v%x_dim_ )
        case (1)
          lf1 => vdfs(p_lower,j)%v%f1; uf1 => vdfs(p_upper,j)%v%f1
          do i1 = l_nr( p_lower, 1 ), l_nr( p_upper, 1 )
            do i = 1, vdfs(p_lower,j)%v%f_dim_
              lf1( i, i1 ) = uf1( i, i1+off(1) )
            enddo
          enddo

        case (2)
          lf2 => vdfs(p_lower,j)%v%f2; uf2 => vdfs(p_upper,j)%v%f2
          do i2 = l_nr( p_lower, 2 ), l_nr( p_upper, 2 )
            do i1 = l_nr( p_lower, 1 ), l_nr( p_upper, 1 )
              do i = 1, vdfs(p_lower,j)%v%f_dim_
                lf2( i, i1, i2 ) = uf2( i, i1+off(1), i2+off(2) )
              enddo
            enddo
          enddo

        case (3)
          lf3 => vdfs(p_lower,j)%v%f3; uf3 => vdfs(p_upper,j)%v%f3
          do i3 = l_nr( p_lower, 3 ), l_nr( p_upper, 3 )
            do i2 = l_nr( p_lower, 2 ), l_nr( p_upper, 2 )
              do i1 = l_nr( p_lower, 1 ), l_nr( p_upper, 1 )
                do i = 1, vdfs(p_lower,j)%v%f_dim_
                  lf3( i, i1, i2, i3 ) = uf3( i, i1+off(1), i2+off(2), i3+off(3) )
                enddo
              enddo
            enddo
          enddo
        end select

        ! copy data from lower tile to upper
        ! compute offset once
        off = l_ns(p_lower,:) - u_nr(p_lower,:)
        select case ( vdfs(p_upper,j)%v%x_dim_ )
        case (1)
          do i1 = u_nr( p_lower, 1 ), u_nr( p_upper, 1 )
            do i = 1, vdfs(p_upper,j)%v%f_dim_
              uf1( i, i1 ) = lf1( i, i1+off(1) )
            enddo
          enddo

        case (2)
          do i2 = u_nr( p_lower, 2 ), u_nr( p_upper, 2 )
            do i1 = u_nr( p_lower, 1 ), u_nr( p_upper, 1 )
              do i = 1, vdfs(p_upper,j)%v%f_dim_
                uf2( i, i1, i2 ) = lf2( i, i1+off(1), i2+off(2) )
              enddo
            enddo
          enddo

        case (3)
          do i3 = u_nr( p_lower, 3 ), u_nr( p_upper, 3 )
            do i2 = u_nr( p_lower, 2 ), u_nr( p_upper, 2 )
              do i1 = u_nr( p_lower, 1 ), u_nr( p_upper, 1 )
                do i = 1, vdfs(p_upper,j)%v%f_dim_
                  uf3( i, i1, i2, i3 ) = lf3( i, i1+off(1), i2+off(2), i3+off(3) )
                enddo
              enddo
            enddo
          enddo
        end select

      enddo

    end select

  endif

end subroutine loc_fill_vdf_tiles
!-----------------------------------------------------------------------------------------

end module m_group_comm

!-----------------------------------------------------------------------------------------
! move tiles to proper process in dynamic load balance
!-----------------------------------------------------------------------------------------
subroutine send_recv_tils_dlb( this, patt, tiles )

  use m_simulation_group, only : t_simulation_group, t_til_container, cleanup_tils_out
  use m_grid_define, only : t_msg_patt
  use m_group_comm, only : irecv_til_size, isend_til_size, irecv_til_msg, isend_til_msg
  use m_group_comm, only : unpack, wait_send
  use m_species_current, only: tmp_xold
  use m_system

  implicit none

  class( t_simulation_group ), intent(inout)      :: this
  type( t_msg_patt ), intent(inout)               :: patt
  class( t_til_container ), dimension(:), pointer :: tiles

  integer, dimension(:), pointer :: recv_bsize
  integer, dimension(:), pointer :: send_bsize
  logical :: alloc_buf

  recv_bsize => null()
  send_bsize => null()
  alloc_buf = associated( tmp_xold )

  ! post receives for size of message
  if (patt%n_recv>0) call irecv_til_size( this, patt, recv_bsize )

  ! calculate message sizes and post sends for size of message
  if (patt%n_send>0) call isend_til_size( this, patt, send_bsize )

  ! wait for recv_bsize, allocate recv buffer and post recv for message
  if (patt%n_recv>0) call irecv_til_msg( patt, recv_bsize, tiles )

  ! wait for send_bsize, allocate send buffer, pack, and post sends
  if (patt%n_send>0) call isend_til_msg( this%no_co, patt, send_bsize, this%tiles )

  ! cleanup outgoing tiles since they have been packed
  ! and we don't want to run out of memory
  call cleanup_tils_out( this, patt )
  ! deallocate tile container
  deallocate( this%tiles )

  ! wait to receive actual message data, unpack
  if (patt%n_recv>0) call unpack( patt, recv_bsize, tiles )

  ! wait for all sends to complete
  if (patt%n_send>0) call wait_send( patt )

  ! cleanup on tiles deallocates module variables that we need, reallocate them
  ! TODO: come up with a better solution for this
  if ( alloc_buf .and. (.not. associated(tmp_xold)) ) &
    call tiles(this%til_id(1))%t%part%species%bnd_con%init_tmp_buf_current()

  ! temp tile array becomes the new tile array
  this%tiles => tiles

end subroutine send_recv_tils_dlb
!-----------------------------------------------------------------------------------------
