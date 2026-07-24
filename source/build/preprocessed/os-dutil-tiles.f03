# 1 "tiles/os-dutil-tiles.f03"
# 1 "<built-in>" 1
# 1 "<built-in>" 3
# 467 "<built-in>" 3
# 1 "<command line>" 1
# 1 "<built-in>" 2
# 1 "tiles/os-dutil-tiles.f03" 2
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
# 2 "tiles/os-dutil-tiles.f03" 2
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
# 3 "tiles/os-dutil-tiles.f03" 2

module m_diag_utilities_tiles

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
# 7 "tiles/os-dutil-tiles.f03" 2

use m_system
use m_parameters
use m_grid_define, only : t_grid, t_msg, t_msg_patt, alloc
use m_grid_tiles, only : t_grid_tiles, t_grid_group, unique_count_1d
use m_node_conf, only : t_node_conf
use m_node_conf_tiles, only : t_node_conf_tiles, node_of_til
use m_tile_conf, only : g_tgp
use m_vdf_define, only : t_vdf, t_vdf_arr, t_vdf_report
use m_vdf_memory
use m_emf_define, only : t_emf, p_extfld_none
use m_emf_diag
use m_time_step, only : t_time_step
use m_vdf_report, only : report_vdf
use m_current_define, only : t_current
use m_current_diag
use m_space, only : t_space
use m_group_comm, only : wait_free_recv, wait_free_send, irecv_til_msg, isend_til_msg

implicit none

private

interface send_recv
  module procedure send_recv_vdf_diag
end interface

interface loc_pack
  module procedure loc_pack_vdf_diag
end interface

interface init_dutil
  module procedure init_diag_utils_emf_group
  module procedure init_diag_utils_jay_group
end interface

interface report_vdf
  module procedure report_vdf_tiles
end interface

interface write_vdf
  module procedure write_vdf_tiles
end interface

public :: send_recv, loc_pack
public :: new_diag_patt_recvs, new_diag_patt_sends
public :: init_dutil, report_vdf, write_vdf

contains

!-----------------------------------------------------------------------------------------
! Send/recv tiles to proper process, pack tile data into group buffers
!-----------------------------------------------------------------------------------------
subroutine send_recv_vdf_diag( group_vdf, til_vdfs, til1_id, no_co, grid, patt, fc_rep )

  implicit none

  type( t_vdf ), intent(inout) :: group_vdf
  type( t_vdf_arr ), dimension(:), intent(in) :: til_vdfs
  integer, intent(in) :: til1_id
  class( t_node_conf ), intent(in) :: no_co
  class( t_grid ), intent(in) :: grid
  type( t_msg_patt ), intent(inout) :: patt
  logical, dimension(:), intent(in) :: fc_rep

  integer, dimension(2,2) :: fc_msg ! ( lower:upper, n_msg ). Maximum of two messages
  integer :: i, n_sub, n_fc, tot_tils
  logical :: in_msg

  ! Get num field components we'll need to report
  n_fc = 0
  do i = 1, size(fc_rep)
    if ( fc_rep(i) ) n_fc = n_fc+1
  enddo

  ! Get the total number of tiles that could be coming in
  select type(no_co); class is( t_node_conf_tiles )
    tot_tils = product( no_co % ti_co % my_nt_io(3,1:p_x_dim) )
  end select

  ! Count number of subarray types necessary for each tile
  ! Loop over the fc_rep array, process message boundary if fc_rep(i) is different from in_msg
  ! Yields fc_msg array, which contains the upper/lower field components for either the
  ! one or two subarray MPI types necessary for each tile.
  ! If field size is not three, then just include all field components in one message from
  ! the minimum needed to the maximum needed
  if ( size(fc_rep) == 3 ) then

    in_msg = .false.
    n_sub = 0
    do i = 1, size(fc_rep)
      if ( fc_rep(i) .neqv. in_msg ) then
        if ( in_msg ) then
          fc_msg( p_upper, n_sub ) = i - 1
          in_msg = .false.
        else
          n_sub = n_sub + 1
          fc_msg( p_lower, n_sub ) = i
          in_msg = .true.
        endif
      endif
    enddo
    if ( in_msg ) fc_msg( p_upper, n_sub ) = size(fc_rep)

  else

    n_sub = 1
    do i = 1, size(fc_rep)
      if (fc_rep(i)) then
        fc_msg( p_lower, n_sub ) = i
        exit
      endif
    enddo

    do i = size(fc_rep), 1, -1
      if (fc_rep(i)) then
        fc_msg( p_upper, n_sub ) = i
        exit
      endif
    enddo

  endif

  ! commit user-defined MPI types and post recvs
  if ( patt%n_recv > 0 ) then
    call irecv_til_msg( group_vdf, no_co, grid, patt%recv, patt%n_recv, &
      n_sub, fc_msg(:,1:n_sub), tot_tils )
  endif

  ! commit user-defined MPI types and post sends
  if ( patt%n_send > 0 ) then
    call isend_til_msg( til_vdfs, til1_id, no_co, patt%send, patt%n_send, &
      n_sub, fc_msg(:,1:n_sub), n_fc, group_vdf%gc_num_ )
  endif

end subroutine send_recv_vdf_diag
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Pack local tile data into group buffers
!-----------------------------------------------------------------------------------------
subroutine loc_pack_vdf_diag( group_vdf, til_vdfs, til1_id, no_co, grid, fc_rep )

  implicit none

  type( t_vdf ), intent(inout) :: group_vdf
  type( t_vdf_arr ), dimension(:), intent(in) :: til_vdfs
  integer, intent(in) :: til1_id
  class( t_node_conf ), intent(in), target :: no_co
  class( t_grid ), intent(in) :: grid
  logical, dimension(:), intent(in) :: fc_rep

  integer, dimension(:,:), pointer :: my_nt
  integer, dimension(p_x_dim) :: tgp
  integer :: t

  select type(no_co); class is (t_node_conf_tiles)
  select type(grid); class is (t_grid_group)

    ! get boundaries for the tiles node is responsible for
    my_nt => no_co%ti_co%my_nt_io(:,1:p_x_dim)

    do t = 1, size(til_vdfs)
      ! check if tile is ours for diags
      tgp = no_co%ti_co%h(:,t+til1_id-1)
      if ( .not.( any(tgp<my_nt(1,:)) .or. any(tgp>my_nt(2,:)) ) ) then

        ! Pack tile data for reports into group buffers
        call pack_fields( group_vdf, til_vdfs(t)%v, no_co, grid, fc_rep, t+til1_id-1 )

      endif
    enddo

  end select
  end select

end subroutine loc_pack_vdf_diag
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! pack tiles' (proper use of "'"?) field data into group's field array
!-----------------------------------------------------------------------------------------
subroutine pack_fields( group_vdf, til_vdf, no_co, grid, fc_rep, t )

  implicit none

  type( t_vdf ), intent(inout) :: group_vdf
  type( t_vdf ), intent(inout) :: til_vdf
  class( t_node_conf_tiles ), intent(in) :: no_co
  class( t_grid_group ), intent(in) :: grid
  logical, dimension(:), intent(in) :: fc_rep
  integer, intent(in) :: t

  ! nx(1:2,dim) = (lower cell bnd, upper cell bnd)
  integer, dimension(2,p_x_dim) :: nx

  integer :: strt, stpp, step, idx_til, fc, fc_off
  integer :: i, j, k

  idx_til = 0
  do i = 1, p_x_dim
    ! include lower/upper guard cells of group_vdf
    nx(1,i) = grid%nx_til( 1, idx_til + g_tgp(no_co%ti_co, t, i) ) - grid%my_nx(1,i) + 1 &
              - group_vdf%gc_num(1,i)
    nx(2,i) = grid%nx_til( 2, idx_til + g_tgp(no_co%ti_co, t, i) ) - grid%my_nx(1,i) + 1 &
              + group_vdf%gc_num(2,i)

    ! increment idx_til (for indexing nx_til)
    idx_til = idx_til + no_co%ti_co%g_til(i)
  enddo

  ! figure out do loop over field components based on fc array
  if ( size(fc_rep) == 3 ) then

    strt = 0; stpp = -1; step = 1
    do i = 1, size(fc_rep)
      if ( fc_rep(i) .and. strt < 1 ) then
        strt = i; stpp = i
      elseif ( fc_rep(i) ) then
        if ( strt == stpp ) step = i - strt
        stpp = i
      endif
    enddo

  else

    step = 1
    do i = 1, size(fc_rep)
      if (fc_rep(i)) then
        strt = i
        exit
      endif
    enddo

    do i = size(fc_rep), 1, -1
      if (fc_rep(i)) then
        stpp = i
        exit
      endif
    enddo

  endif

  ! if group_vdf has only one fc and til_vdf has many, offset group_vdf fc
  if ( group_vdf%f_dim_ < size(fc_rep) ) then
    ! sanity check
    if ( group_vdf%f_dim_ /= 1 ) then
      write(err_buf__,*) 'Field dimensions differ, but group_vdf field dimension is not 1';call err__("tiles/os-dutil-tiles.f03",254)
    endif
    fc_off = 1 - strt
  else
    fc_off = 0
  endif

  ! copy tile field data into group field array
  select case (p_x_dim)
  case(1)

    ! group_vdf indexing starts at nx(1,1), til_vdf indexing starts at 1
    do i = nx(1,1), nx(2,1)
      do fc = strt, stpp, step
        group_vdf % f1( fc+fc_off, i ) &
          = til_vdf % f1( fc, i-nx(1,1)+1 )
      enddo
    enddo

  case(2)

    ! group_vdf indexing starts at nx(1,1), til_vdf indexing starts at 1
    do j = nx(1,2), nx(2,2)
      do i = nx(1,1), nx(2,1)
        do fc = strt, stpp, step
          group_vdf % f2( fc+fc_off, i, j ) &
            = til_vdf % f2( fc, i-nx(1,1)+1, j-nx(1,2)+1 )
        enddo
      enddo
    enddo

  case(3)

    ! group_vdf indexing starts at nx(1,1), til_vdf indexing starts at 1
    do k = nx(1,3), nx(2,3)
      do j = nx(1,2), nx(2,2)
        do i = nx(1,1), nx(2,1)
          do fc = strt, stpp, step
            group_vdf % f3( fc+fc_off, i, j, k ) &
              = til_vdf % f3( fc, i-nx(1,1)+1, j-nx(1,2)+1, k-nx(1,3)+1 )
          enddo
        enddo
      enddo
    enddo

  end select

end subroutine pack_fields
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Find which nodes we need to recv tiles from
! Need to make sure send/recv node pair order tiles in the same way
!-----------------------------------------------------------------------------------------
subroutine new_diag_patt_recvs( no_co, patt )

  implicit none

  class( t_node_conf_tiles ), intent(in), target :: no_co
  type( t_msg_patt ), intent(inout) :: patt ! who's patt?

  integer, dimension(:,:), pointer :: my_nt

  ! msg_nodes(1:2,i)=(node of sending til, sending til)
  integer, dimension(:,:), pointer :: msg_nodes

  integer :: til, node
  integer :: i, j, k, c

  my_nt => no_co%ti_co%my_nt_io(:,1:p_x_dim)

  ! compute max number of cells we could have to recv (j)
  j = 1
  do i = 1, p_x_dim
    j = j * ( my_nt(2,i)-my_nt(1,i)+1)
  enddo

  ! worst case scenario we pass all local tiles
  call alloc(msg_nodes, (/2,j/),"tiles/os-dutil-tiles.f03",332)

  ! Identify tiles not stored locally that we need for diags
  select case ( p_x_dim )
  case (1)

    c = 0
    ! loop over tiles we need for diags
    do i = my_nt(1,1), my_nt(2,1)
      til = no_co%ti_co%hi1(i)
      node = node_of_til( no_co, til )

      ! check if tile currently stored on different node
      if ( node /= no_co%my_aid() ) then
        c = c+1
        ! store til and node to send it to in msg_nodes
        msg_nodes(1,c) = node
        msg_nodes(2,c) = til
      endif
    enddo

  case (2)

    c = 0
    ! loop over tiles we need for diags
    do j = my_nt(1,2), my_nt(2,2)
      do i = my_nt(1,1), my_nt(2,1)
        til = no_co%ti_co%hi2(i,j)
        node = node_of_til( no_co, til )

        ! check if tile currently stored on different node
        if ( node /= no_co%my_aid() ) then
          c = c+1
          ! store til and node to send it to in msg_nodes
          msg_nodes(1,c) = node
          msg_nodes(2,c) = til
        endif
      enddo
    enddo

  case (3)

    c = 0
    ! loop over tiles we need for diags
    do k = my_nt(1,3), my_nt(2,3)
      do j = my_nt(1,2), my_nt(2,2)
        do i = my_nt(1,1), my_nt(2,1)
          til = no_co%ti_co%hi3(i,j,k)
          node = node_of_til( no_co, til )

          ! check if tile currently stored on different node
          if ( node /= no_co%my_aid() ) then
            c = c+1
            ! store til and node to send it to in msg_nodes
            msg_nodes(1,c) = node
            msg_nodes(2,c) = til
          endif
        enddo
      enddo
    enddo

  end select

  ! Find unique nodes to recv tils from, which tils are coming from those nodes
  call new_diag_msg( patt%recv, patt%n_recv, msg_nodes, c )

  ! (*debug*)
  ! do i = 1, patt%n_recv
  !
  ! enddo

  ! Cleanup
  call freemem(msg_nodes,"tiles/os-dutil-tiles.f03",404)

end subroutine new_diag_patt_recvs
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Find which nodes we need to send tiles to for diags and which tiles to send
! Need to make sure send/recv node pair order tiles in the same way
!-----------------------------------------------------------------------------------------
subroutine new_diag_patt_sends( no_co, til_id, patt )

  implicit none

  class( t_node_conf_tiles ), intent(in), target :: no_co
  integer, dimension(3), intent(in) :: til_id
  type( t_msg_patt ), intent(inout) :: patt ! who's patt?

  integer, dimension(:,:), pointer :: my_nt

  ! msg_nodes(1:2,i)=(node of sending til, sending til)
  integer, dimension(:,:), pointer :: msg_nodes

  ! window(1:2,dim) = (lower bnd, upper bnd)
  integer, dimension(2,p_x_dim) :: window

  integer, dimension(p_x_dim) :: ngp,tgp

  integer :: q, r, til, node
  integer :: i, j, k, l, c

  my_nt => no_co%ti_co%my_nt_io(:,1:p_x_dim)

  ! worst case scenario we pass all local tiles
  call alloc(msg_nodes, (/2,til_id(3)/),"tiles/os-dutil-tiles.f03",437)

  ! find min and max coords in hid in each dimension of our tiles
  ! this way we don't have to loop through the entire hid array
  window(1,:) = no_co%ti_co%h(:,til_id(1))
  window(2,:) = no_co%ti_co%h(:,til_id(1))
  do til = til_id(1)+1, til_id(2)
    tgp = no_co%ti_co%h(:,til)
    do i = 1, p_x_dim
      if ( tgp(i)<window(1,i) ) then
        ! set a new min of window
        window(1,i) = tgp(i)
      elseif ( tgp(i)>window(2,i) ) then
        ! set a new max of window
        window(2,i) = tgp(i)
      endif
    enddo
  enddo

  ! Identify tiles stored locally that we need to send to a different node for diags
  select case ( p_x_dim )
  case (1)

    c = 0
    do i = window(1,1), window(2,1)
      til = no_co%ti_co%hi1(i)
      tgp = no_co%ti_co%h(:,til)

      ! check if til on this node AND need to go to different node for diags
      if ( ( til>=til_id(1) .and. til<=til_id(2) ) .and. &
           ( any( tgp<my_nt(1,:) ) .or. any( tgp>my_nt(2,:) ) ) ) then

        ! find node til needs to be on
        do l = 1, p_x_dim
          q = no_co%ti_co%g_til(l)/no_co%nx(l)
          r = mod(no_co%ti_co%g_til(l),no_co%nx(l))

          if ( tgp(l) <= r*q + r ) then
            ngp(l) = (tgp(l)-1)/(q+1) + 1
          else
            ngp(l) = (tgp(l)-(r*q+r)-1)/q + r + 1
          endif
        enddo
        node = no_co % coords_rank( ngp-1 ) + 1

        c = c+1

        ! store til and node to send it to in msg_nodes
        msg_nodes(1,c) = node
        msg_nodes(2,c) = til

      endif
    enddo

  case (2)

    c = 0
    do j = window(1,2), window(2,2)
      do i = window(1,1), window(2,1)
        til = no_co%ti_co%hi2(i,j)
        tgp = no_co%ti_co%h(:,til)

        ! check if til on this node AND need to go to different node for diags
        if ( ( til>=til_id(1) .and. til<=til_id(2) ) .and. &
             ( any( tgp<my_nt(1,:) ) .or. any( tgp>my_nt(2,:) ) ) ) then

          ! find node til needs to be on
          do l = 1, p_x_dim
            q = no_co%ti_co%g_til(l)/no_co%nx(l)
            r = mod(no_co%ti_co%g_til(l),no_co%nx(l))

            if ( tgp(l) <= r*q + r ) then
              ngp(l) = (tgp(l)-1)/(q+1) + 1
            else
              ngp(l) = (tgp(l)-(r*q+r)-1)/q + r + 1
            endif
          enddo
          node = no_co % coords_rank( ngp-1 ) + 1

          c = c+1

          ! store til and node to send it to in msg_nodes
          msg_nodes(1,c) = node
          msg_nodes(2,c) = til

        endif
      enddo
    enddo

  case (3)

    c = 0
    do k = window(1,3), window(2,3)
      do j = window(1,2), window(2,2)
        do i = window(1,1), window(2,1)
          til = no_co%ti_co%hi3(i,j,k)
          tgp = no_co%ti_co%h(:,til)

          ! check if til on this node AND need to go to different node for diags
          if ( ( til>=til_id(1) .and. til<=til_id(2) ) .and. &
               ( any( tgp<my_nt(1,:) ) .or. any( tgp>my_nt(2,:) ) ) ) then

            ! find node til needs to be on
            do l = 1, p_x_dim
              q = no_co%ti_co%g_til(l)/no_co%nx(l)
              r = mod(no_co%ti_co%g_til(l),no_co%nx(l))

              if ( tgp(l) <= r*q + r ) then
                ngp(l) = (tgp(l)-1)/(q+1) + 1
              else
                ngp(l) = (tgp(l)-(r*q+r)-1)/q + r + 1
              endif
            enddo
            node = no_co % coords_rank( ngp-1 ) + 1

            c = c+1

            ! store til and node to send it to in msg_nodes
            msg_nodes(1,c) = node
            msg_nodes(2,c) = til

          endif
        enddo
      enddo
    enddo

  end select

  ! Find unique nodes to send tils to, which tils are going to those nodes
  call new_diag_msg( patt%send, patt%n_send, msg_nodes, c )

  ! (*debug*)
  ! do i = 1, patt%n_send
  !
  ! enddo

  ! Cleanup
  call freemem(msg_nodes,"tiles/os-dutil-tiles.f03",574)

end subroutine new_diag_patt_sends
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine new_diag_msg( msg, n_msg, msg_nodes, c )

  implicit none

  type(t_msg), dimension(:), pointer :: msg
  integer, intent(inout) :: n_msg
  integer, dimension(:,:), pointer :: msg_nodes
  integer, intent(in) :: c

  integer, dimension(:,:), pointer :: unique

  integer :: cc, node
  integer :: i, j

  unique => null()

  if ( c==0 ) then
    n_msg = 0
  else
    call unique_count_1d( msg_nodes(1,1:c), unique )

    n_msg = size(unique,1)
    call alloc(msg, (/ n_msg /),"tiles/os-dutil-tiles.f03",603)

    ! loop over unique nodes, figure out which tiles they need
    do i = 1, n_msg
      node = unique(i,1)
      msg(i)%node = node

      ! number of tiles going to that node
      call alloc(msg(i)%tils, (/unique(i,2)/),"tiles/os-dutil-tiles.f03",611)

      ! loop through all out of bounds tiles, pack ones for this node into msg patt
      cc = 0
      do j = 1, c
        if ( msg_nodes(1,j)==node ) then
          cc = cc+1
          msg(i)%tils(cc) = msg_nodes(2,j)
        endif
      enddo

    enddo
  endif

  if ( associated(unique) ) call freemem(unique,"tiles/os-dutil-tiles.f03",625)

end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Allocate emf objects for diagnostics, setting up fields to hold emf data for reports
!-----------------------------------------------------------------------------------------
subroutine init_diag_utils_emf_group( this, grid, dx, part_fld_alloc, interpolation )

  implicit none

  class( t_emf ), intent(inout), target :: this
  class( t_grid ), intent(in) :: grid
  real(p_double), dimension(:), intent(in) :: dx
  logical, intent(in) :: part_fld_alloc
  integer, intent(in) :: interpolation

  type( t_vdf_report ), pointer :: rep
  integer, dimension( 2, p_x_dim ) :: min_gc

  ! dumping divergence, psi, and poynting diagnostics requires one guard cell
  min_gc = 1

  this%part_fld_alloc = .false.

  rep => this%diag%reports
  do
    if ( .not. associated( rep ) ) exit

    select case ( rep%quant )
    case ( p_e1, p_e2, p_e3, p_ene_e1, p_ene_e2, p_ene_e3, p_ene_e, p_div_e, p_psi )
      if ( this%e%f_dim_ < 1 ) then
        call this%e%new( p_x_dim, p_f_dim, grid%my_nx(3,:), min_gc, dx, .false. )
      endif

    case ( p_b1, p_b2, p_b3, p_ene_b1, p_ene_b2, p_ene_b3, p_ene_b, p_div_b )
      if ( this%b%f_dim_ < 1 ) then
        call this%b%new( p_x_dim, p_f_dim, grid%my_nx(3,:), min_gc, dx, .false. )
      endif

    case ( p_ext_e1, p_ext_e2, p_ext_e3 )
      if ( this%ext_e%f_dim_ < 1 .and. this%ext_fld /= p_extfld_none ) then
        call this%ext_e%new( p_x_dim, p_f_dim, grid%my_nx(3,:), min_gc, dx, .false. )
      endif

    case ( p_ext_b1, p_ext_b2, p_ext_b3 )
      if ( this%ext_b%f_dim_ < 1 .and. this%ext_fld /= p_extfld_none ) then
        call this%ext_b%new( p_x_dim, p_f_dim, grid%my_nx(3,:), min_gc, dx, .false. )
      endif

    case ( p_part_e1, p_part_e2, p_part_e3, p_part_b1, p_part_b2, p_part_b3 )
      if ( part_fld_alloc ) then
        this%part_fld_alloc = .true.
        if ( .not. associated(this%e_part) ) then
          call alloc(this%e_part,"tiles/os-dutil-tiles.f03",680)
          call alloc(this%b_part,"tiles/os-dutil-tiles.f03",681)
          call this%e_part%new( p_x_dim, p_f_dim, grid%my_nx(3,:), min_gc, dx, .false. )
          call this%b_part%new( p_x_dim, p_f_dim, grid%my_nx(3,:), min_gc, dx, .false. )
        endif
      else
        this%part_fld_alloc = .false.
        if ( this%e%f_dim_ < 1 ) then
          call this%e%new( p_x_dim, p_f_dim, grid%my_nx(3,:), min_gc, dx, .false. )
        endif
        if ( this%b%f_dim_ < 1 ) then
          call this%b%new( p_x_dim, p_f_dim, grid%my_nx(3,:), min_gc, dx, .false. )
        endif
        this%e_part => this%e
        this%b_part => this%b
      endif

    case ( p_ene_emf, p_s1, p_s2, p_s3 )
      if ( this%e%f_dim_ < 1 ) then
        call this%e%new( p_x_dim, p_f_dim, grid%my_nx(3,:), min_gc, dx, .false. )
      endif
      if ( this%b%f_dim_ < 1 ) then
        call this%b%new( p_x_dim, p_f_dim, grid%my_nx(3,:), min_gc, dx, .false. )
      endif

    case ( p_charge_cons )
      write(err_buf__,*) "Charge conservation diagnostic temporarily disabled";call err__("tiles/os-dutil-tiles.f03",706)
      call abort_program()

    case default
      ! unknown quantity, must belong to a subclass
      continue

    end select

    rep => rep%next
  enddo

  call this % diag % init( this%ext_fld == p_extfld_none, this%part_fld_alloc, &
                           interpolation )

end subroutine init_diag_utils_emf_group
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Allocate emf objects for diagnostics, setting up fields to hold emf data for reports
!-----------------------------------------------------------------------------------------
subroutine init_diag_utils_jay_group( this, grid, dx )

  implicit none

  class( t_current ), intent(inout) :: this
  class( t_grid ), intent(in) :: grid
  real(p_double), dimension(:), intent(in) :: dx

  type( t_vdf_report ), pointer :: rep
  integer, dimension( 2, p_x_dim ) :: min_gc

  ! don't need to dump guard cells so don't need to allocate space for them
  ! except that divergence diagnostics require one guard cell
  min_gc = 1

  ! Always allocate this array since it is called in report_current
  call alloc(this%pf, (/1/),"tiles/os-dutil-tiles.f03",743)

  rep => this%diag%reports
  do
    if ( .not. associated( rep ) ) exit

    select case ( rep%quant )
    case ( p_j1, p_j2, p_j3, p_div_j )
      if ( this%pf(1)%f_dim_ < 1 ) then
        call this%pf(1)%new( p_x_dim, p_f_dim, grid%my_nx(3,:), min_gc, dx, .false. )
      endif

    case default
      ! unknown quantity, must belong to a subclass
      continue

    end select

    rep => rep%next
  enddo

  call this % diag % init( this%interpolation )

end subroutine init_diag_utils_jay_group
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Process all reports for the quantity specified in report
!-----------------------------------------------------------------------------------------
subroutine report_vdf_tiles( report, vdf_arr, til1_id, fc, g_space, grid, no_co, patt, tstep, t )

  implicit none

  type(t_vdf_report), pointer :: report
  type(t_vdf_arr), dimension(:), intent(in) :: vdf_arr
  integer, intent(in) :: til1_id, fc
  type(t_space), intent(in) :: g_space
  class(t_grid), intent(in) :: grid
  class( t_node_conf ), intent(in) :: no_co
  type(t_msg_patt), intent(inout) :: patt
  type(t_time_step), intent(in) :: tstep
  real(p_double), intent(in) :: t

  integer, dimension( 2, p_x_dim ) :: min_gc
  type(t_vdf) :: vdf_data
  logical, dimension( vdf_arr(1)%v%f_dim_ ) :: fc_rep

  ! Sanity check
  if ( vdf_arr(1)%v%x_dim_ /= p_x_dim ) then
    write(err_buf__,*) 'report_vdf_tiles not built for vdfs with different';call err__("tiles/os-dutil-tiles.f03",792)
    write(err_buf__,*) 'dimension than p_x_dim.  Aborting...';call err__("tiles/os-dutil-tiles.f03",793)
    call abort_program(p_err_notimplemented)
  endif

  ! Create vdf on each node corresponding to regularized OSIRIS grid
  ! No guard cells necessary
  min_gc = 0
  call vdf_data % new( p_x_dim, 1, grid%my_nx(3,:), min_gc, vdf_arr(1)%v%dx_(1:p_x_dim), &
                        .false. )

  fc_rep = .false.
  fc_rep(fc) = .true.

  ! Post non-blocking send/recv messages to proper process using defined MPI types
  call send_recv( vdf_data, vdf_arr, til1_id, no_co, grid, patt, fc_rep )

  ! Pack local tile data into group vdf
  call loc_pack( vdf_data, vdf_arr, til1_id, no_co, grid, fc_rep )

  ! wait for receive and free MPI types
  call wait_free_recv( patt )

  ! wait for send and free MPI types
  call wait_free_send( patt )

  call report_vdf( report, vdf_data, fc, g_space, grid, no_co, tstep, t )

  call vdf_data % cleanup()

end subroutine report_vdf_tiles
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Save vdf grid to file
!-----------------------------------------------------------------------------------------
subroutine write_vdf_tiles( vdf_arr, report, til1_id, fc, g_space, grid, no_co, patt )

  implicit none

  type(t_vdf_arr), dimension(:), intent(in) :: vdf_arr
  type(t_vdf_report), intent(in) :: report
  integer, intent(in) :: til1_id, fc
  type(t_space), intent(in) :: g_space
  class(t_grid), intent(in) :: grid
  class( t_node_conf ), intent(in) :: no_co
  type(t_msg_patt), intent(inout) :: patt

  integer, dimension( 2, p_x_dim ) :: min_gc
  type(t_vdf) :: vdf_data
  logical, dimension( vdf_arr(1)%v%f_dim_ ) :: fc_rep

  ! Sanity check
  if ( vdf_arr(1)%v%x_dim_ /= p_x_dim ) then
    write(err_buf__,*) 'write_vdf_tiles not built for vdfs with different';call err__("tiles/os-dutil-tiles.f03",846)
    write(err_buf__,*) 'dimension than p_x_dim.  Aborting...';call err__("tiles/os-dutil-tiles.f03",847)
    call abort_program(p_err_notimplemented)
  endif

  ! Create vdf on each node corresponding to regularized OSIRIS grid
  ! No guard cells necessary
  min_gc = 0
  call vdf_data % new( p_x_dim, 1, grid%my_nx(3,:), min_gc, vdf_arr(1)%v%dx_(1:p_x_dim), &
                        .false. )

  fc_rep = .false.
  fc_rep(fc) = .true.

  ! Post non-blocking send/recv messages to proper process using defined MPI types
  call send_recv( vdf_data, vdf_arr, til1_id, no_co, grid, patt, fc_rep )

  ! Pack local tile data into group vdf
  call loc_pack( vdf_data, vdf_arr, til1_id, no_co, grid, fc_rep )

  ! wait for receive and free MPI types
  call wait_free_recv( patt )

  ! wait for send and free MPI types
  call wait_free_send( patt )

  call vdf_data % write( report, fc, g_space, grid, no_co )

  call vdf_data % cleanup()

end subroutine write_vdf_tiles
!-----------------------------------------------------------------------------------------


end module m_diag_utilities_tiles


!-----------------------------------------------------------------------------------------
! report on electro-magnetic field diagnostics for tiles code
!-----------------------------------------------------------------------------------------
subroutine report_emf_group( sim )

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
# 889 "tiles/os-dutil-tiles.f03" 2

  use m_system
  use m_simulation_group, only : t_simulation_group
  use m_emf_diag
  use m_time, only : t
  use m_logprof, only : begin_event, end_event
  use m_diag_utilities_tiles, only : send_recv, loc_pack
  use m_vdf_define, only : t_vdf_arr, t_vdf_report, t_vdf
  use m_vdf_report, only : if_report
  use m_group_comm, only : wait_free_recv, wait_free_send

  implicit none

  class( t_simulation_group ), intent(inout) :: sim

  logical, dimension(p_max_dim), target :: fc_e, fc_b, fc_ext_e, fc_ext_b
  logical, dimension(p_max_dim), target :: fc_part_e, fc_part_b
  logical, dimension(:), pointer :: fc_rep
  logical :: val
  type( t_vdf_report ), pointer :: rep
  integer :: id, i
  type( t_vdf_arr ), dimension(sim%til_id(1):sim%til_id(2)) :: til_vdfs
  type( t_vdf ), pointer :: group_vdf

  call begin_event( diag_emf_ev )

  ! Get needed field components for e and b by looping over report types
  fc_e = .false.; fc_b = .false.
  fc_ext_e = .false.; fc_ext_b = .false.
  fc_part_e = .false.; fc_part_b = .false.
  rep => sim%emf%diag%reports
  do
    if ( .not. associated( rep ) ) exit
    if ( if_report( rep, sim%tstep ) ) then

      select case ( rep%quant )
      case ( p_e1, p_e2, p_e3 )
        fc_e( rep%quant - p_e1 + 1 ) = .true.

      case ( p_b1, p_b2, p_b3 )
        fc_b( rep%quant - p_b1 + 1 ) = .true.

      case ( p_ext_e1, p_ext_e2, p_ext_e3 )
        fc_ext_e( rep%quant - p_ext_e1 + 1 ) = .true.

      case ( p_ext_b1, p_ext_b2, p_ext_b3 )
        fc_ext_b( rep%quant - p_ext_b1 + 1 ) = .true.

      case ( p_part_e1, p_part_e2, p_part_e3 )
        if ( sim%emf%part_fld_alloc ) then
          fc_part_e( rep%quant - p_part_e1 + 1 ) = .true.
        else
          fc_e( rep%quant - p_part_e1 + 1 ) = .true.
        endif

      case ( p_part_b1, p_part_b2, p_part_b3 )
        if ( sim%emf%part_fld_alloc ) then
          fc_part_b( rep%quant - p_part_b1 + 1 ) = .true.
        else
          fc_b( rep%quant - p_part_b1 + 1 ) = .true.
        endif

      case ( p_ene_e1, p_ene_e2, p_ene_e3 )
        fc_e( rep%quant - p_ene_e1 + 1 ) = .true.

      case ( p_ene_b1, p_ene_b2, p_ene_b3 )
        fc_b( rep%quant - p_ene_b1 + 1 ) = .true.

      case ( p_ene_e )
        fc_e = .true.

      case ( p_ene_b )
        fc_b = .true.

      case ( p_ene_emf )
        fc_e = .true.
        fc_b = .true.

      case ( p_div_e )
        fc_e = .true.

      case ( p_div_b )
        fc_b = .true.

      case ( p_charge_cons )
        ! print *, 'Reporting charge conservation '
        ! call report_vdf( rep, this%f, 1, g_space, grid, no_co, tstep, t )
        write(err_buf__,*) "Charge conservation diagnostic temporarily disabled";call err__("tiles/os-dutil-tiles.f03",976)
        call abort_program()

      case ( p_psi )
        fc_e = .true.

      case ( p_s1, p_s2, p_s3 )
        val = fc_e( rep%quant - p_s1 + 1 )
        fc_e = .true.
        fc_e( rep%quant - p_s1 + 1 ) = val
        val = fc_b( rep%quant - p_s1 + 1 )
        fc_b = .true.
        fc_b( rep%quant - p_s1 + 1 ) = val

      case default
        ! unknown quantity, must belong to a subclass
        ! print *, 'In report_diag_emf, unknown quantity : ', trim(rep%name)
        continue

      end select

    endif
    rep => rep%next
  enddo

  ! For each vdf necessary to put into regular partition, first post receives/sends for
  ! external communications, then pack in local tiles, then unpack/pack external communications
  do id = p_e1, p_part_b1, 3

    select case( id )
    case( p_e1 )

      do i = sim%til_id(1), sim%til_id(2)
        til_vdfs(i) % v => sim % tiles(i) % t % emf % e
      enddo
      group_vdf => sim % emf % e
      fc_rep => fc_e

    case( p_b1 )

      do i = sim%til_id(1), sim%til_id(2)
        til_vdfs(i) % v => sim % tiles(i) % t % emf % b
      enddo
      group_vdf => sim % emf % b
      fc_rep => fc_b

    case( p_ext_e1 )

      do i = sim%til_id(1), sim%til_id(2)
        til_vdfs(i) % v => sim % tiles(i) % t % emf % ext_e
      enddo
      group_vdf => sim % emf % ext_e
      fc_rep => fc_ext_e

    case( p_ext_b1 )

      do i = sim%til_id(1), sim%til_id(2)
        til_vdfs(i) % v => sim % tiles(i) % t % emf % ext_b
      enddo
      group_vdf => sim % emf % ext_b
      fc_rep => fc_ext_b

    case( p_part_e1 )

      do i = sim%til_id(1), sim%til_id(2)
        til_vdfs(i) % v => sim % tiles(i) % t % emf % e_part
      enddo
      group_vdf => sim % emf % e_part
      fc_rep => fc_part_e

    case( p_part_b1 )

      do i = sim%til_id(1), sim%til_id(2)
        til_vdfs(i) % v => sim % tiles(i) % t % emf % b_part
      enddo
      group_vdf => sim % emf % b_part
      fc_rep => fc_part_b

    end select

    if ( .not. any(fc_rep) ) cycle

    ! Post non-blocking send/recv messages to proper process using defined MPI types
    call send_recv( group_vdf, til_vdfs, sim%til_id(1), sim%no_co, &
                                  sim%grid, sim%diag_patt, fc_rep )

    ! Pack local tile data into group vdf
    call loc_pack( group_vdf, til_vdfs, sim%til_id(1), sim%no_co, sim%grid, fc_rep )

    ! wait for receive and free MPI types
    call wait_free_recv( sim%diag_patt )

    ! wait for send and free MPI types
    call wait_free_send( sim%diag_patt )

  enddo

!
  ! Call the regular report emf function
  call sim % emf % report( sim%g_space, sim%grid, sim%no_co, sim%tstep, t(sim%time), &
    sim%bnd%send_vdf, sim%bnd%recv_vdf )

  call end_event( diag_emf_ev )

end subroutine report_emf_group
!-----------------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! report particles, had to put everything at the top level
!-------------------------------------------------------------------------------
subroutine report_part_group( sim )

  use m_system
  use m_parameters
  use m_species_define, only : t_spec_arr
  use m_particles_define, only : p_charge, p_charge_htc, p_dcharge_dt, diag_part_ev
  use m_particles_charge, only : get_charge_htc, get_dcharge_dt
  use m_simulation_group, only : t_simulation_group
  use m_logprof, only : begin_event, end_event
  use m_time, only : t
  use m_vdf_define, only : t_vdf_arr, t_vdf_report
  use m_vdf_report, only : if_report
  use m_species_phasespace_tiles, only : report
  use m_time_step, only : test_if_report, n, ndump, dt
  use m_group_boundary, only : deposit_quant
  use m_species_rawdiag_tiles, only : write_raw
  use m_diag_utilities_tiles, only : report_vdf
  use m_species_diagnostics, only : p_norm
  use m_vdf_math, only : mult
  use m_neutral, only : t_neutral
  use m_diag_neutral, only : p_ion_charge, p_neut_den

  implicit none

  class(t_simulation_group), intent(inout) :: sim

  integer :: i, j
  type( t_vdf_report ), pointer :: rep
  type( t_vdf_arr ), dimension(:), pointer :: charge
  type( t_spec_arr ), dimension(:), pointer :: spec_arr
  logical :: needs_update_charge

  spec_arr => null()
  charge => null()

  call begin_event(diag_part_ev)

  ! Species diagnostics (phasespaces, energy, temperature, heat flux, RAW, tracks, etc.)
  allocate( spec_arr( sim%til_id(3) ) )
  do i = 1, sim%til_id(3)
    spec_arr(i)%s => sim%tiles(i+sim%til_id(1)-1)%t%part%species
  enddo
  do
    if (.not. associated(spec_arr(1)%s)) exit
    call report_spec_group( spec_arr, sim )
    do i = 1, sim%til_id(3)
      spec_arr(i)%s => spec_arr(i)%s%next
    enddo
  enddo
  deallocate(spec_arr)


  ! Allocate array of vdfs, one per tile
  if ( .not. associated(charge) ) then
    allocate( charge(sim%til_id(1):sim%til_id(2)) )
  endif

  ! Report each neutral, one at a time
  do i=1, sim%tiles(sim%til_id(1))%t%part%num_neutral

    do j = sim%til_id(1), sim%til_id(2)
      charge(j)%v => sim%tiles(j)%t%part%neutral(i)%multi_ion
    enddo

    call report_neutral_group( sim%tiles(sim%til_id(1))%t%part%neutral(i), charge, sim )

    do j = sim%til_id(1), sim%til_id(2)
      charge(j)%v => null()
    enddo

  enddo



  ! Regular particle diagnostics
  ! Check if charge deposit is required
  needs_update_charge = .false.
  rep => sim%tiles(sim%til_id(1))%t%part%reports
  do
    if ( .not. associated( rep ) .or. needs_update_charge ) exit

    needs_update_charge = if_report( rep, sim%tstep ) .or. &
                 ( ( rep%quant == p_charge_htc .or. rep%quant == p_dcharge_dt ) .and. &
                 if_report( rep, sim%tstep, iter = +1 ) )

    rep => rep%next
  enddo

  ! Deposit charge if required
  if ( needs_update_charge ) then
    call sim % update_charge()
  endif

  ! Do selected reports
  rep => sim%tiles(sim%til_id(1))%t%part%reports
  do
    if (.not. associated(rep)) exit

    if ( if_report( rep, sim%tstep ) ) then

      ! Allocate array of vdfs, one per tile
      if ( .not. associated(charge) ) then
        allocate( charge(sim%til_id(1):sim%til_id(2)) )
      endif

      ! deposit present charge

      select case ( rep%quant )
      case ( p_charge )

        do i = sim%til_id(1), sim%til_id(2)
          charge(i)%v => sim%tiles(i)%t%part%charge%current
        enddo
        call report_vdf( rep, charge, sim%til_id(1), 1, sim%g_space, sim%grid, &
                          sim%no_co, sim%diag_patt, sim%tstep, t(sim%time) )
        do i = sim%til_id(1), sim%til_id(2)
          charge(i)%v => null()
        enddo

      case ( p_charge_htc )

        do i = sim%til_id(1), sim%til_id(2)
          allocate( charge(i)%v )
          call get_charge_htc( sim%tiles(i)%t%part, charge(i)%v )
        enddo
        call report_vdf( rep, charge, sim%til_id(1), 1, sim%g_space, sim%grid, &
                          sim%no_co, sim%diag_patt, sim%tstep, t(sim%time) )
        do i = sim%til_id(1), sim%til_id(2)
          call charge(i) % v % cleanup()
          deallocate( charge(i)%v )
        enddo

      case ( p_dcharge_dt )

        do i = sim%til_id(1), sim%til_id(2)
          allocate( charge(i)%v )
          call get_dcharge_dt( sim%tiles(i)%t%part, dt( sim%tstep ), charge(i)%v )
        enddo
        call report_vdf( rep, charge, sim%til_id(1), 1, sim%g_space, sim%grid, &
                          sim%no_co, sim%diag_patt, sim%tstep, t(sim%time) )
        do i = sim%til_id(1), sim%til_id(2)
          call charge(i) % v % cleanup()
          deallocate( charge(i)%v )
        enddo

      end select

    endif

    rep => rep % next

  enddo

  if ( associated(charge) ) deallocate(charge)

  call end_event(diag_part_ev)

  contains

  subroutine report_spec_group( spec_arr, sim )

    implicit none

    type( t_spec_arr ), dimension(:), pointer :: spec_arr
    class( t_simulation_group ), intent(inout) :: sim

    type( t_vdf_arr ), dimension(:), pointer :: charge, norm
    type( t_vdf_report ), pointer :: rep
    integer :: i

    charge => null(); norm => null()

    ! Density diagnostics
    rep => spec_arr(1)%s%diag%reports
    do
      if (.not. associated(rep)) exit

      if ( if_report( rep, sim%tstep ) ) then

        ! Allocate array of vdfs, one per tile
        if ( .not. associated(charge) ) allocate( charge(1:sim%til_id(3)) )

        ! deposit quantity (this also creates the required charge vdf)
        call deposit_quant( spec_arr, charge, rep%quant, sim )

        ! do all reports for current quantity
        call report_vdf( rep, charge, sim%til_id(1), 1, sim%g_space, sim%grid, &
                            sim%no_co, sim%diag_patt, sim%tstep, t(sim%time) )

        ! cleanup temp. grid
        do i = 1, sim%til_id(3)
          call charge(i) % v % cleanup()
          deallocate( charge(i)%v )
        enddo

      endif

      rep => rep % next

    enddo

    ! Momentum distribution diagnostics are not currently implemented

    ! Cell average diagnostics
    rep => spec_arr(1)%s%diag%rep_cell_avg
    do
      if ( .not. associated( rep ) ) exit

      if ( if_report( rep, sim%tstep ) ) then

        ! deposit normalization value
        if ( .not. associated(norm) ) then
          allocate( norm(1:sim%til_id(3)) )
          call deposit_quant( spec_arr, norm, p_norm, sim )
        endif

        ! Allocate array of vdfs, one per tile
        if ( .not. associated(charge) ) allocate( charge(1:sim%til_id(3)) )

        ! deposit quantity and normalize it
        call deposit_quant( spec_arr, charge, rep%quant, sim )
        do i = 1, sim%til_id(3)
          call mult( charge(i)%v, norm(i)%v )
        enddo

        ! do all reports for current quantity
        call report_vdf( rep, charge, sim%til_id(1), 1, sim%g_space, sim%grid, &
                            sim%no_co, sim%diag_patt, sim%tstep, t(sim%time) )

        ! cleanup temp. grid
        do i = 1, sim%til_id(3)
          call charge(i) % v % cleanup()
          deallocate( charge(i)%v )
        enddo
      endif

      rep => rep%next
    enddo
    if (associated(norm)) then
      do i = 1, sim%til_id(3)
        call norm(i) % v % cleanup()
        deallocate( norm(i)%v )
      enddo
      deallocate(norm)
    endif

    ! Phasespaces
    call report( spec_arr(1)%s%diag%phasespaces, spec_arr, sim%no_co, sim%g_space, &
                  sim%tstep, t(sim%time) )

    ! Raw particle data diagnostics
    if (test_if_report( sim%tstep, spec_arr(1)%s%diag%ndump_fac_raw )) then

      call write_raw( spec_arr, sim%no_co, n(sim%tstep), t(sim%time), &
                      n(sim%tstep) / ndump(sim%tstep) )

    endif

    ! Temperature and heat flux diagnostics not currently implemented

    ! Particle tracking diagnostics not currently implemented

    if ( associated(charge) ) deallocate( charge )

  end subroutine report_spec_group


  subroutine report_neutral_group( neut, multi_ion, sim )

    implicit none

    class( t_neutral ), intent(inout) :: neut
    ! multi_ion is indexed from sim%til_id(1):sim%til_id(2) here
    type( t_vdf_arr ), dimension(:), pointer :: multi_ion
    class( t_simulation_group ), intent(inout) :: sim

    type( t_vdf_arr ), dimension(:), pointer :: neut_den, neut_default
    type( t_vdf_report ), pointer :: rep
    integer :: i

    neut_den => null(); neut_default => null()

    rep => neut%diag%reports
    do
      if ( .not. associated( rep ) ) exit

      if ( if_report( rep, sim%tstep ) ) then
        select case ( rep%quant )
        case ( p_ion_charge )
          if ( mpi_node() == 0 ) then
            write(0,*) "Warning, ion_charge diagnostic may be incorrect."
          endif
          call report_vdf( rep, multi_ion, sim%til_id(1), neut%ion_idx, sim%g_space, &
                            sim%grid, sim%no_co, sim%diag_patt, sim%tstep, t(sim%time) )

        case ( p_neut_den )

          if ( .not. associated(neut_den) ) then
            allocate( neut_den(sim%til_id(1):sim%til_id(2)) )
            allocate( neut_default(sim%til_id(1):sim%til_id(2)) )
          endif

          ! get real neutral density
          do i = sim%til_id(1), sim%til_id(2)
            allocate( neut_den(i)%v, neut_default(i)%v )
            call neut_den(i)%v % new( multi_ion(i)%v, copy = .true., f_dim = 1, fc = 1 )
            call neut_default(i)%v % new( multi_ion(i)%v, copy = .true., f_dim = 1, &
                                          fc = neut%neut_idx )
            call mult( neut_den(i)%v, neut_default(i)%v )
          enddo

          ! save result
          call report_vdf( rep, neut_den, sim%til_id(1), 1, sim%g_space, sim%grid, &
                            sim%no_co, sim%diag_patt, sim%tstep, t(sim%time) )

          ! free used memory
          do i = sim%til_id(1), sim%til_id(2)
            call neut_den(i)%v % cleanup()
            call neut_default(i)%v % cleanup()
            deallocate( neut_den(i)%v, neut_default(i)%v )
          enddo

        end select
      endif

      rep => rep%next
    enddo

    if ( associated(neut_den) ) deallocate( neut_den, neut_default )

  end subroutine report_neutral_group


end subroutine report_part_group
!-------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! report on current - diagnostic
!-----------------------------------------------------------------------------------------
subroutine report_current_group( sim )

  use m_system
  use m_parameters
  use m_time, only : t
  use m_current_diag
  use m_logprof, only : begin_event, end_event
  use m_vdf_define, only : t_vdf, t_vdf_arr, t_vdf_report
  use m_vdf_report, only : if_report
  use m_diagnostic_utilities, only : t_diag_file
  use m_diag_utilities_tiles, only : send_recv, loc_pack
  use m_simulation_group, only : t_simulation_group
  use m_group_comm, only : wait_free_recv, wait_free_send

  implicit none

  class( t_simulation_group ), intent(inout), target :: sim

  logical, dimension(p_max_dim), target :: fc_jay
  type( t_vdf_report ), pointer :: rep
  integer :: i
  type( t_vdf_arr ), dimension(sim%til_id(1):sim%til_id(2)) :: til_vdfs
  type( t_vdf ), pointer :: group_vdf

  call begin_event( diag_current_ev )

  ! Get needed field components for jay by looping over report types
  fc_jay = .false.
  rep => sim%jay%diag%reports
  do
    if ( .not. associated( rep ) ) exit
    if ( if_report( rep, sim%tstep ) ) then

      select case ( rep%quant )
      case ( p_j1, p_j2, p_j3 )
        fc_jay( rep%quant - p_j1 + 1 ) = .true.

      case ( p_div_j )
        fc_jay = .true.

      case default
        ! unknown quantity, must belong to a subclass
        continue

      end select

    endif
    rep => rep%next
  enddo

  ! If required to put jay into regular partition, first post receives/sends for external
  ! communications, then pack in local tiles, then unpack/pack external communications
  if ( any(fc_jay) ) then

    do i = sim%til_id(1), sim%til_id(2)
      til_vdfs(i) % v => sim % tiles(i) % t % jay % pf(1)
    enddo
    group_vdf => sim % jay % pf(1)

    ! Post non-blocking send/recv messages to proper process using defined MPI types
    call send_recv( group_vdf, til_vdfs, sim%til_id(1), sim%no_co, &
                                  sim%grid, sim%diag_patt, fc_jay )

    ! Pack local tile data into group vdf
    call loc_pack( group_vdf, til_vdfs, sim%til_id(1), sim%no_co, sim%grid, fc_jay )

    ! wait for receive and free MPI types
    call wait_free_recv( sim%diag_patt )

    ! wait for send and free MPI types
    call wait_free_send( sim%diag_patt )

  endif

  ! Call the regular report jay function
  call sim % jay % report( sim%g_space, sim%grid, sim%no_co, sim%tstep, t(sim%time) )

  call end_event( diag_current_ev )

end subroutine report_current_group
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine report_energy_emf_group( sim )

  use m_simulation_group, only : t_simulation_group
  use m_time, only : t, tmin
  use m_time_step
  use m_emf_diag
  use m_logprof
  use m_parameters
  use m_node_conf
  use m_vdf_math

  implicit none

  class( t_simulation_group ), intent(inout) :: sim

  real(p_double), dimension(2*p_f_dim) :: temp_int, temp_int_tot

  integer :: ierr, i

  ! file name and path
  character(len=256) :: full_name, path

  call begin_event( diag_emf_ev )

  ! reports on integrated field energy
  if (test_if_report( sim%tstep, sim%emf%diag%ndump_fac_ene_int ) ) then

    call sim%emf%fill_data( sim%tstep )

    ! get local e and b field integrals
    do i = 1, 2*p_f_dim
      temp_int_tot(i) = 0
    enddo

    do i = sim%til_id(1), sim%til_id(2)
      call total( sim%tiles(i)%t%emf%b, temp_int, pow = 2 )
      call total( sim%tiles(i)%t%emf%e, temp_int(p_f_dim+1:), pow = 2 )
      temp_int_tot = temp_int_tot + temp_int
    enddo

    ! sum up results from all nodes
    call reduce_array( sim%no_co, temp_int_tot, operation = p_sum)

    ! save the data
    if ( root(sim%no_co) ) then
      ! normalize energies
      temp_int_tot = temp_int_tot * sim % tiles(sim%til_id(1)) % t % emf % e % dvol() &
                       * 0.5_p_double

      ! setup path and file names
      path = trim(path_hist)
      full_name = trim(path) // 'fld_ene'

      ! Open file and position at the last record
      if ( t(sim%time) == tmin(sim%time) ) then

        call mkdir( path, ierr )

        open (unit=file_id_fldene, file=full_name, status = 'REPLACE' , &
              form='formatted')

        ! Write Header
        write(file_id_fldene, '(A)' ) '! EM field energy per field component'
        write(file_id_fldene,'( A6, 1X,A15,6(1X,A23) )') &
                                'Iter','Time    ','B1','B2','B3',&
                                'E1','E2','E3'
      else

        open (unit=file_id_fldene, file=full_name, position = 'append', &
              form='formatted')
      endif

      write(file_id_fldene, '( I6, 1X,g15.8, 6(1X,es23.16) )') n(sim%tstep), t(sim%time), temp_int_tot
      close(file_id_fldene)

    endif

  endif

  call end_event( diag_emf_ev )

end subroutine report_energy_emf_group
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine report_energy_part_group( sim )

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
# 1598 "tiles/os-dutil-tiles.f03" 2

  use m_simulation_group, only : t_simulation_group
  use m_particles_define
  use m_time_step
  use m_logprof
  use m_parameters, only: file_id_parene, p_err_diagfile, path_hist
  use m_node_conf
  use m_species_define, only : t_species, t_spec_arr
  use m_system
  use m_time, only : t, tmin

  implicit none

  class( t_simulation_group ), intent(inout) :: sim

  integer :: i, j, ierr
  integer, parameter :: p_header_length = 23
  real(p_double), dimension(:), pointer :: ene_all
  character(len = p_header_length), dimension(:), allocatable :: header
  character(len = 64) :: fmt
  character(len = 80) :: path, full_name

  character(len=*), parameter :: p_err_msg = "Unable to create directory for particle energy diagnostic"

  class( t_species ), pointer :: species

  type( t_spec_arr ), dimension(:), pointer :: spec_arr

  call begin_event( diag_part_ev )

  ! Report global particle energy
  if ( test_if_report( sim%tstep, sim%part%ndump_fac_ene ) ) then

    call alloc(ene_all, (/ sim%part%num_species + 1 /),"tiles/os-dutil-tiles.f03",1631)

    do i = 1, sim%part%num_species + 1
      ene_all(i) = 0
    enddo

    ! Get normalized kinetic energy for all species and calculate total part. energy
    do j = sim%til_id(1), sim%til_id(2)
      species => sim % tiles(j) % t % part % species
      i = 1
      do
        if (.not. associated(species)) exit
        ene_all(i+1) = ene_all(i+1) + species % get_energy()
        ene_all(1) = ene_all(1) + species % get_energy()
        species => species % next
        i = i + 1
      enddo
    enddo

    ! Gather results on root node
    call reduce_array( sim%no_co, ene_all )

    ! Write output to file
    if ( root(sim%no_co) ) then

      ! prepare path and file names

      path = trim(path_hist)

      full_name = trim(path) // 'par_ene'

      ! Open file and position at the last record
      if ( t(sim%time) == tmin(sim%time) ) then

        call mkdir( path, ierr )
        call check_error(ierr,p_err_msg,p_err_diagfile,"tiles/os-dutil-tiles.f03",1666)

        open (unit=file_id_parene, file=full_name, status = 'REPLACE' , form='formatted')

        ! Write Header

        write(file_id_parene, '(A)' ) '! Total and per species particle kinetic energy'

        allocate( header( sim%part%num_species ) )

        write( fmt, '(A,I0,A,I0,A)' ) '( A6, (1X,A15), ', sim%part%num_species + 1, &
                                        '(1X, A', p_header_length, ') )'

        i = 1
        species => sim % tiles(sim%til_id(1)) % t % part % species
        do
          if (.not. associated(species)) exit
          header(i) = trim(species%name)
          header(i) = adjustr(header(i))
          species => species % next
          i = i+1
        enddo

        write(file_id_parene, fmt ) 'Iter','Time    ', 'Total', header
        deallocate( header )
      else

        open (unit=file_id_parene, file=full_name, position = 'append', &
        form='formatted')

      endif

      ! Writes the particle energies to file
      write( fmt, '(A,I0,A)' ) '( I6, (1X,g15.8), ', sim%part%num_species + 1, '(1X, es23.16) )'
      write(file_id_parene, fmt) n(sim%tstep), t(sim%time), ene_all

      ! Close file
      close(file_id_parene)

    endif

    call freemem(ene_all,"tiles/os-dutil-tiles.f03",1707)

  endif

  ! Report individual species energy
  allocate( spec_arr( sim%til_id(3) ) )
  do i = 1, sim%til_id(3)
    spec_arr(i)%s => sim%tiles(i+sim%til_id(1)-1)%t%part%species
  enddo
  do
    if (.not. associated(spec_arr(1)%s)) exit
    call report_energy_spec_group( spec_arr, sim )
    do i = 1, sim%til_id(3)
      spec_arr(i)%s => spec_arr(i)%s%next
    enddo
  enddo
  deallocate(spec_arr)

  call end_event( diag_part_ev )


  contains


  !----------
  ! Energy diagnostic is performed in double precision
  ! - Time centered kinetic energy has been previously calculated during the push
  !----------
  subroutine report_energy_spec_group( spec_arr, sim )

    implicit none

    type( t_spec_arr ), dimension(:), pointer :: spec_arr
    class( t_simulation_group ), intent(inout) :: sim

    character(80) :: path, full_name
    integer :: ierr, i
    integer(p_int64) :: total_par
    real(p_double) :: energy

    character(len=*), parameter :: p_err_msg_spec = "Unable to create directory for species energy diagnostic"

    if (test_if_report( sim%tstep, spec_arr(1)%s%diag%ndump_fac_ene )) then

      ! Get energy on local node (sum over all local tiles)
      energy = 0
      do i = 1, sim%til_id(3)
        energy = energy + spec_arr(i) % s % get_energy()
      enddo

      ! sums the results from all nodes
      call reduce(sim%no_co, energy, operation = p_sum)

      ! get the total number of particles in the simulation box (sum over all local tiles)
      total_par = 0
      do i = 1, sim%til_id(3)
        total_par = total_par + spec_arr(i) % s % num_par
      enddo

      ! sums the results from all nodes
      call reduce(sim%no_co, total_par, operation = p_sum)

      if ( root(sim%no_co) ) then

        ! File path
        path = trim(path_hist)

        ! file name including path
        write( full_name, '(A,A,I0.2,A)' ) trim(path), 'par', spec_arr(1)%s%sp_id, '_ene'


        ! Open file and position at the last record
        if ( t(sim%time) == tmin(sim%time) ) then
          call mkdir( path, ierr )
          call check_error(ierr,p_err_msg_spec,p_err_diagfile,"tiles/os-dutil-tiles.f03",1781)

          open (unit=file_id_parene, file=full_name, status = 'REPLACE' , &
            form='formatted')

          ! Write Header
          write(file_id_parene, '(A,A,A)' ) '! Kinetic energy for species "', &
                trim( spec_arr(1)%s%name ), '"'
          write(file_id_parene,'( A6, 2(1X,A15), (1X,A23) )') &
                'Iter','Time    ','Total Par.','Kin. Energy'

        else
          open (unit=file_id_parene, file=full_name, position = 'append', &
            form='formatted')
        endif

        ! Writes the particle energies to file
        write(file_id_parene, '( I6,1X,g15.8,1X,I14,1(1X,es23.16) )') &
        n(sim%tstep), t(sim%time), total_par, energy
        close(file_id_parene)

      endif

    endif


  end subroutine report_energy_spec_group
  !----------


end subroutine report_energy_part_group
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Allocate objects for diagnostics
!-----------------------------------------------------------------------------------------
subroutine init_diag_utils_group( sim, dx )

  use m_simulation_group, only : t_simulation_group
  use m_parameters, only : p_double
  use m_diagnostic_utilities, only : p_max_diag_buffer_size, init_dutil
  use m_diag_utilities_tiles, only : init_dutil

  implicit none

  class( t_simulation_group ), intent(inout) :: sim
  real(p_double), dimension(:), intent(in) :: dx

  ! Set up group objects for reporting emf
  call init_dutil( sim%emf, sim%grid, dx, sim%tiles(sim%til_id(1))%t%emf%part_fld_alloc, &
                   sim%part%interpolation )

  ! Set up group objects for reporting jay
  call init_dutil( sim%jay, sim%grid, dx )

  ! No initialization necessary for particles

  ! Call regular osiris init_dutil
  call init_dutil( sim%options, p_max_diag_buffer_size )

end subroutine init_diag_utils_group
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! find outgoing/incoming tiles and their destinations
!-----------------------------------------------------------------------------------------
subroutine new_diag_patt( sim )

  use m_node_conf_tiles, only : t_node_conf_tiles
  use m_simulation_group, only : t_simulation_group
  use m_diag_utilities_tiles, only : new_diag_patt_recvs, new_diag_patt_sends

  implicit none

  class( t_simulation_group ), intent(inout) :: sim

  call sim % diag_patt % cleanup()

  select type(no_co => sim%no_co); class is (t_node_conf_tiles)

    ! Find which nodes we need to recv tiles from for diags
    call new_diag_patt_recvs( no_co, sim%diag_patt )

    ! Find which nodes we need to send tiles to for diags and which tiles to send
    call new_diag_patt_sends( no_co, sim%til_id, sim%diag_patt )

  end select

end subroutine new_diag_patt
!-----------------------------------------------------------------------------------------
