# 1 "tiles/os-group-boundary.f03"
# 1 "<built-in>" 1
# 1 "<built-in>" 3
# 467 "<built-in>" 3
# 1 "<command line>" 1
# 1 "<built-in>" 2
# 1 "tiles/os-group-boundary.f03" 2
module m_group_boundary

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
# 4 "tiles/os-group-boundary.f03" 2
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
# 5 "tiles/os-group-boundary.f03" 2
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
# 6 "tiles/os-group-boundary.f03" 2

use m_system
use m_parameters

use m_simulation_group, only : t_simulation_group
use m_simulation, only : t_simulation
use m_node_conf, only : t_node_conf, neighbor, n_threads
use m_node_conf_tiles, only : t_node_conf_tiles
use m_vdf_define, only : t_vdf_arr
use m_vdf_comm, only : t_vdf_msg, wait, update_periodic_vdf
use m_vdf_comm_tiles, only : t_bnd_msg, p_send_msg, p_recv_msg, irecv_size, isend_size
use m_vdf_comm_tiles, only : p_msg, p_wait, irecv_tiles
use m_emf_define, only : t_lindman, t_vpml
use m_wall_define, only : t_wall
use m_wall_comm, only : update_periodic_wall, wait_irecv_wall, wait_isend_wall
use m_wall_comm_tiles, only : loc_recv_wall, loc_send_wall, irecv_wall_tiles, isend_wall_tiles
use m_species_define, only : t_spec_msg, t_spec_arr, t_species
use m_species_comm, only : remove_particles, wait_spec_msg, isend_2, int_pack_size
use m_species_comm_tiles, only : p_lower_phys, p_upper_phys, irecv_2, loc_recv, irecv_1, loc_send
use m_species_comm_tiles, only : set_spec_msg
use m_species_boundary, only : check_boundary, check_node_cross, process_phys_boundary
use m_bnd_tiles, only : t_bnd_group, t_bnd_patt
use m_group_comm, only : unpack, isend_size_spec, isend_1, irecv_wait, isend, loc_fill
use m_group_comm, only : deposit_quant_1, deposit_quant_2


implicit none

private

! For vpml updates with tiles
integer, parameter :: p_wall_e = 1
integer, parameter :: p_wall_b = 2

interface send_recv
  module procedure send_recv_vdf_ub
end interface

interface wait_all
  module procedure wait_all_vdf_ub
end interface

interface loc_update_boundary
  module procedure loc_update_boundary_vdf
end interface

interface update_boundary
  module procedure update_boundary_lindman
  module procedure update_boundary_wall
  module procedure update_boundary_vpml
  module procedure update_boundary_vdf_arr_tiles
end interface

interface deposit_quant
  module procedure deposit_quant_tiles
end interface

interface copy_spec_msg
  module procedure copy_spec_msg
end interface

interface check_buffer_size
  module procedure check_buffer_size_vdf
  module procedure check_buffer_size_part
end interface

interface check_ncross_grp
  module procedure check_node_cross_group
end interface

interface check_bnd_grp
  module procedure check_boundary_group
end interface

interface loc_ub_spec_1
  module procedure loc_update_bnd_spec_1_group
end interface

interface loc_ub_spec_2
  module procedure loc_update_bnd_spec_2_group
end interface

interface ext_ub_spec_1
  module procedure ext_update_bnd_spec_1_group
end interface

interface ext_ub_spec_2
  module procedure ext_update_bnd_spec_2_group
end interface

interface ext_ub_spec_3
  module procedure ext_update_bnd_spec_3_group
end interface

interface rmv_part_grp
  module procedure remove_particles_group
end interface

public :: send_recv, wait_all, loc_update_boundary
public :: update_boundary, deposit_quant, copy_spec_msg, check_buffer_size
public :: check_bnd_grp, check_ncross_grp
public :: loc_ub_spec_1, loc_ub_spec_2
public :: ext_ub_spec_1, ext_ub_spec_2, ext_ub_spec_3
public :: rmv_part_grp
public :: p_wall_e, p_wall_b

contains

!-----------------------------------------------------------------------------------------
! Post sends and receives for vdf boundary values for tiles on communication boundaries
!-----------------------------------------------------------------------------------------
subroutine send_recv_vdf_ub( vdf_arr_a, update_type, dim, nx_move, til1_id, &
                                no_co, patt, send_msg, recv_msg, vdf_arr_b )

  implicit none

  type( t_vdf_arr ), dimension(:), intent(in) :: vdf_arr_a
  integer, intent(in) :: update_type, dim, nx_move, til1_id
  class( t_node_conf ), intent(in) :: no_co
  class( t_bnd_patt ), intent(inout) :: patt
  type( t_vdf_msg ), intent(inout), target :: send_msg, recv_msg
  type( t_vdf_arr ), dimension(:), intent(in), optional :: vdf_arr_b

  type(t_vdf_msg), pointer :: snd_vdf, rcv_vdf
  integer :: i, j, t_idx, bnd
  class(t_bnd_msg), pointer :: msg

  snd_vdf => send_msg
  rcv_vdf => recv_msg

  ! Post sends and receives for all external communications
  do i = 1, patt%n_msg

    msg => patt%msg(i)

    ! Zero out message size
    msg%msg_size(:,:) = 0

    ! Count size of each message
    do j = 1, msg%n_tils
      t_idx = msg%tils(1,j) - til1_id + 1
      bnd = msg%tils(2,j)

      ! Count recv size
      if (present(vdf_arr_b)) then
        call irecv_size( vdf_arr_a(t_idx)%v, update_type, dim, bnd, &
          nx_move, msg, j, vdf_arr_b(t_idx)%v )
      else
        call irecv_size( vdf_arr_a(t_idx)%v, update_type, dim, bnd, &
          nx_move, msg, j )
      endif

      ! Count send size
      if (present(vdf_arr_b)) then
        call isend_size( vdf_arr_a(t_idx)%v, update_type, dim, bnd, &
          nx_move, msg, j, vdf_arr_b(t_idx)%v )
      else
        call isend_size( vdf_arr_a(t_idx)%v, update_type, dim, bnd, &
          nx_move, msg, j )
      endif
    enddo

    ! Add up all message sizes
    call msg % total_size()

    ! Reallocate buffer if necessary, both for send and receive
    call check_buffer_size( msg, rcv_vdf, snd_vdf )

    ! Post receive
    call irecv_tiles( update_type, msg%node, no_co, rcv_vdf )

    ! Pack data and post send
    call isend( vdf_arr_a, til1_id, no_co, msg, snd_vdf, vdf_arr_b )

    ! point to next element in linked list
    rcv_vdf => rcv_vdf%next
    snd_vdf => snd_vdf%next

    ! TODO: process lindman boundaries if needed
  enddo

end subroutine send_recv_vdf_ub
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Wait for and unpack vdf boundary values for tiles on communication boundaries
!-----------------------------------------------------------------------------------------
subroutine wait_all_vdf_ub( vdf_arr_a, til1_id, patt, send_msg, recv_msg, vdf_arr_b )

  implicit none

  type( t_vdf_arr ), dimension(:), intent(inout) :: vdf_arr_a
  integer, intent(in) :: til1_id
  class( t_bnd_patt ), intent(inout) :: patt
  type( t_vdf_msg ), intent(inout), target :: send_msg, recv_msg
  type( t_vdf_arr ), dimension(:), intent(inout), optional :: vdf_arr_b

  type(t_vdf_msg), pointer :: snd_vdf, rcv_vdf
  integer :: i
  class(t_bnd_msg), pointer :: msg

  ! Wait for and unpack all external recvs
  rcv_vdf => recv_msg
  do i = 1, patt%n_msg

    msg => patt%msg(i)

    if (present(vdf_arr_b)) then
      call irecv_wait( vdf_arr_a, til1_id, msg, rcv_vdf, vdf_arr_b )
    else
      call irecv_wait( vdf_arr_a, til1_id, msg, rcv_vdf )
    endif
    rcv_vdf => rcv_vdf%next

  enddo

  ! Wait for all external sends to complete
  snd_vdf => send_msg
  do i = 1, patt%n_msg

    call wait( snd_vdf )
    snd_vdf => snd_vdf%next

  enddo

end subroutine wait_all_vdf_ub
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Update vdf boundary values on local tile boundaries
!-----------------------------------------------------------------------------------------
subroutine loc_update_boundary_vdf( vdf_arr_a, update_type, dim, nx_move, til1_id, &
                                      msg, vdf_arr_b )

  implicit none

  type( t_vdf_arr ), dimension(:), intent(in) :: vdf_arr_a
  integer, intent(in) :: update_type, dim, nx_move, til1_id
  class( t_bnd_msg ), intent(inout) :: msg
  type( t_vdf_arr ), dimension(:), intent(in), optional :: vdf_arr_b

  type( t_vdf_arr ), dimension(2,2) :: vdfs ! ( p_lower/p_upper, a/b )
  integer :: j, bnd, d2

  ! Get size of second dimension of vdfs to pass to loc_fill
  if (present(vdf_arr_b)) then
    d2 = 2
  else
    d2 = 1
  endif

  !$omp parallel do schedule(static) private(vdfs,bnd)
  do j = 1, msg%n_tils

    ! l_til has the lower g_til_aid of the two, not necessary the p_lower bnd
    vdfs(p_lower,1)%v => vdf_arr_a( msg%tils(1,j) - til1_id + 1 )%v
    vdfs(p_upper,1)%v => vdf_arr_a( msg%tils(3,j) - til1_id + 1 )%v
    if (present(vdf_arr_b)) then
      vdfs(p_lower,2)%v => vdf_arr_b( msg%tils(1,j) - til1_id + 1 )%v
      vdfs(p_upper,2)%v => vdf_arr_b( msg%tils(3,j) - til1_id + 1 )%v
    endif

    bnd = msg%tils(2,j)

    call loc_fill( vdfs(:,1:d2), update_type, dim, bnd, nx_move )

  enddo
  !$omp end parallel do

end subroutine loc_update_boundary_vdf
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Calls all functions necessary for updating boundaries of a vdf array
!-----------------------------------------------------------------------------------------
subroutine update_boundary_vdf_arr_tiles( vdf_arr_a, update_type, til1_id, no_co, &
                                no_co_til, patt, send_msg, recv_msg, vdf_arr_b, nx_move_ )

  implicit none

  type( t_vdf_arr ), dimension(:), intent(inout) :: vdf_arr_a
  integer, intent(in) :: update_type, til1_id
  class( t_node_conf ), intent(in) :: no_co, no_co_til
  class( t_bnd_patt ), dimension(:), intent(inout) :: patt
  type( t_vdf_msg ), intent(inout), target :: send_msg, recv_msg
  type( t_vdf_arr ), dimension(:), intent(inout), optional :: vdf_arr_b
  integer, dimension(:), intent(in), optional :: nx_move_

  integer, dimension(p_max_dim) :: nx_move
  integer :: my_node, dim, i
  integer, dimension(p_x_dim) :: g_til

  if (present(nx_move_)) then
    nx_move( 1:vdf_arr_a(1)%v%x_dim_ ) = nx_move_( 1:vdf_arr_a(1)%v%x_dim_ )
  else
    nx_move = 0
  endif

  my_node = no_co%my_aid()
  select type(no_co_til); class is(t_node_conf_tiles)
    g_til = no_co_til%ti_co%g_til
  end select

  do dim = 1, p_x_dim

    ! Check if single node periodic
    if ( g_til(dim) == 1 .and. my_node == neighbor(no_co_til,p_lower,dim) ) then

      do i = 1, size(vdf_arr_a)
        call update_periodic_vdf( vdf_arr_a(i)%v, dim, update_type )
        if (present(vdf_arr_b)) then
          call update_periodic_vdf( vdf_arr_b(i)%v, dim, update_type )
        endif
      enddo

    else

      if (present(vdf_arr_b)) then

        call send_recv( vdf_arr_a, update_type, dim, nx_move(dim), til1_id, &
                        no_co, patt(dim), send_msg, recv_msg, vdf_arr_b )

        call loc_update_boundary( vdf_arr_a, update_type, dim, nx_move(dim), til1_id, &
                                  patt(dim)%msg_loc, vdf_arr_b )

        call wait_all( vdf_arr_a, til1_id, patt(dim), send_msg, recv_msg, vdf_arr_b )

      else

        call send_recv( vdf_arr_a, update_type, dim, nx_move(dim), til1_id, &
                        no_co, patt(dim), send_msg, recv_msg )

        call loc_update_boundary( vdf_arr_a, update_type, dim, nx_move(dim), til1_id, &
                                  patt(dim)%msg_loc )

        call wait_all( vdf_arr_a, til1_id, patt(dim), send_msg, recv_msg )

      endif

    endif

  enddo

end subroutine update_boundary_vdf_arr_tiles
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Deposit species charge on a grid for all tiles
!-----------------------------------------------------------------------------------------
subroutine deposit_quant_tiles( spec_arr, charge, quant, sim )

  implicit none

  type( t_spec_arr ), dimension(:), intent(in) :: spec_arr
  type( t_vdf_arr ), dimension(:), intent(inout) :: charge
  integer, intent(in) :: quant
  class( t_simulation_group ), intent(inout) :: sim

  integer :: i, t

  do i = 1, sim%til_id(3)
    t = i + sim%til_id(1) - 1
    if (.not. associated(charge(i)%v)) allocate( charge(i)%v )
    call deposit_quant_1( spec_arr(i)%s, sim%tiles(t)%t%grid, charge(i)%v, quant )
  enddo

  ! Call update boundary on all the tiles
  call update_boundary( charge, p_vdf_add, sim%til_id(1), sim%no_co, &
    sim%tiles(sim%til_id(1))%t%no_co, sim%bnd_patt, sim%bnd_grp%send_vdf, &
    sim%bnd_grp%recv_vdf )

  ! second half of deposit quantity
  do i = 1, sim%til_id(3)
    t = i + sim%til_id(1) - 1
    call deposit_quant_2( spec_arr(i)%s, sim%tiles(t)%t%grid, charge(i)%v, quant )
  enddo

end subroutine deposit_quant_tiles
!-----------------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine update_boundary_lindman(this, no_co, nx_move, send_msg, recv_msg, step)
!-------------------------------------------------------------------------------

  implicit none

  type (t_lindman), intent(inout) :: this
  class(t_node_conf),intent(in) :: no_co
  integer, dimension(:), intent(in):: nx_move
  type(t_vdf_msg), dimension(2), intent(inout) :: send_msg, recv_msg
  integer, intent(in) :: step


  ! only update boundaries if boundary is active
  if (this%idir > 0) then
    write(err_buf__,*) "Lindman not yet allowed with tiling";call err__("tiles/os-group-boundary.f03",401)
    ! call update_boundary(this%b_buffer, p_vdf_replace, no_co, nx_move, send_msg, recv_msg, step)
    ! call update_boundary(this%e_buffer, p_vdf_replace, no_co, nx_move, send_msg, recv_msg, step)
  end if

end subroutine update_boundary_lindman
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Update boundary values of a wall on communication and parallel boundaries
!-------------------------------------------------------------------------------
subroutine update_boundary_wall( wall, update_type, no_co, move_num, send_msg, recv_msg, &
                                 step, dim, sim )

  implicit none

  type( t_wall ), intent(inout) :: wall
  integer, intent(in) :: update_type
  class( t_node_conf ), intent(in) :: no_co
  integer, dimension(:), intent(in) :: move_num
  type(t_vdf_msg), dimension(2), intent(inout) :: send_msg, recv_msg
  integer, intent(in) :: step, dim
  class(t_simulation_group), intent(inout) :: sim

  integer :: my_node, l_node, u_node
  type(t_vdf_msg), dimension(2) :: loc_msg

  my_node = no_co%my_aid()

  ! skip the direction perpendicular to the wall
  if (dim /= wall%idir) then

    select type(no_co); class is(t_node_conf_tiles)

    if ( no_co%ti_co%g_til(dim) == 1 .and. my_node == neighbor(no_co,p_lower,dim) ) then

       !single node periodic
       call update_periodic_wall( wall, dim, update_type )

    else

      ! communication with other nodes
      l_node = neighbor( no_co, p_lower, dim )
      u_node = neighbor( no_co, p_upper, dim )

      select case(step)
      case(p_msg)

        ! post receives
        if ( l_node > 0 .and. l_node/=my_node ) then
          call irecv_wall_tiles( wall, update_type, dim, p_lower, no_co, move_num(dim), recv_msg )
        endif
        if ( u_node > 0 .and. u_node/=my_node ) then
          call irecv_wall_tiles( wall, update_type, dim, p_upper, no_co, move_num(dim), recv_msg )
        endif

        ! post sends
        ! Note: sends MUST be posted in the opposite order of receives to account for
        ! a 2 node periodic partition, where 1 node sends 2 messages to the same node
        if ( u_node > 0 .and. u_node/=my_node ) then
          call isend_wall_tiles( wall, update_type, dim, p_upper, no_co, move_num(dim), send_msg )
        elseif ( u_node > 0 ) then
          call loc_send_wall( wall, update_type, dim, p_upper, move_num(dim), send_msg )
        endif

        if ( l_node > 0 .and. l_node/=my_node ) then
          call isend_wall_tiles( wall, update_type, dim, p_lower, no_co, move_num(dim), send_msg )
        elseif( l_node > 0 ) then
          call loc_send_wall( wall, update_type, dim, p_lower, move_num(dim), send_msg )
        endif

      case(p_wait)

        ! wait for receives and unpack data
        if ( l_node > 0 .and. l_node/=my_node ) then
          call wait_irecv_wall( wall, p_lower, recv_msg )
        elseif ( l_node > 0 ) then
          loc_msg(p_lower)%buffer => sim%tiles( no_co%ti_co%neighbor_g_til(p_lower,dim) ) &
            %t%bnd%send_vdf(p_upper)%buffer
          call loc_recv_wall( wall, update_type, dim, p_lower, move_num(dim), loc_msg )
        endif
        if ( u_node > 0 .and. u_node/=my_node ) then
          call wait_irecv_wall( wall, p_upper, recv_msg )
        elseif ( u_node > 0 ) then
          loc_msg(p_upper)%buffer => sim%tiles( no_co%ti_co%neighbor_g_til(p_upper,dim) ) &
            %t%bnd%send_vdf(p_lower)%buffer
          call loc_recv_wall( wall, update_type, dim, p_upper, move_num(dim), loc_msg )
        endif

        ! wait for sends to complete
        if ( u_node > 0 .and. u_node/=my_node ) then
          call wait_isend_wall( p_upper, send_msg )
        endif
        if ( l_node > 0 .and. l_node/=my_node ) then
          call wait_isend_wall( p_lower, send_msg )
        endif

      end select

    endif

    end select

  endif


end subroutine update_boundary_wall
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine update_boundary_vpml( this, no_co, if_move, nx_move, send_msg, recv_msg, &
                                step, dim, eb, sim )
!-------------------------------------------------------------------------------

  implicit none

  type (t_vpml), intent(inout) :: this
  class(t_node_conf),intent(in) :: no_co
  logical, dimension(:), intent(in):: if_move
  integer, dimension(:), intent(in):: nx_move
  type(t_vdf_msg), dimension(2), intent(inout) :: send_msg, recv_msg
  integer, intent(in) :: step
  integer, intent(in) :: dim ! dimension
  integer, intent(in) :: eb ! e wall or b wall
  class(t_simulation_group), intent(inout) :: sim

  ! local variables
  integer :: i_vpml, j, move_gc
  integer, dimension( 2, p_x_dim ) :: lgc_num

  ! The guard cell values on vpml walls (which are created by the new_vpml_wall routine) don't have
  ! exactly the same meaning as in normal walls, which breaks the communication routines.

  ! This needs to change, but in the mean time just temporarily overwrite the gc_num member

  ! variable with the standard gc_num value of gc_num(p_lower,:) = 1 and gc_num(p_upper,:) = 2
  ! (if doing a moving window then gc_num(p_lower,1) = 2)

  do while (this%c <= this%n_vpml)

    i_vpml = this%c

    if (dim == this%wall_array_e(i_vpml)%idir) then

      this%c = i_vpml + 1

    else

      ! Store original gc_num
      lgc_num = this%wall_array_e(i_vpml) % gc_num ( :, 1:p_x_dim )

      ! Set it to standard EMF gc_num
      do j = 1, this%wall_array_e(i_vpml)%x_dim

        if ( if_move( j ) ) then
          move_gc = 1
        else
          move_gc = 0
        endif

        this%wall_array_e(i_vpml) % gc_num( p_lower, j ) = 1 + move_gc
        this%wall_array_e(i_vpml) % gc_num( p_upper, j ) = 2

        this%wall_array_b(i_vpml) % gc_num( p_lower, j ) = 1 + move_gc
        this%wall_array_b(i_vpml) % gc_num( p_upper, j ) = 2
      enddo

      ! Do a normal wall update boundary
      if (eb==p_wall_e) then
        call update_boundary(this%wall_array_e(i_vpml), p_vdf_replace, no_co, nx_move, send_msg, recv_msg, step, dim, sim )
      elseif (eb==p_wall_b) then
        call update_boundary(this%wall_array_b(i_vpml), p_vdf_replace, no_co, nx_move, send_msg, recv_msg, step, dim, sim )
      endif

      ! Set it back to original value
      this%wall_array_e(i_vpml) % gc_num ( :, 1:p_x_dim ) = lgc_num
      this%wall_array_b(i_vpml) % gc_num ( :, 1:p_x_dim ) = lgc_num

      if (step == p_wait) this%c = i_vpml + 1
      exit

    endif

  enddo

end subroutine update_boundary_vpml
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Copy function replaces an mpi communication when species communication
! is local to node
!-------------------------------------------------------------------------------
subroutine copy_spec_msg ( source_msg, dest_msg )

  implicit none

  type(t_spec_msg), intent(in) :: source_msg
  type(t_spec_msg), intent(inout) :: dest_msg

  dest_msg % n_part = source_msg % n_part
  dest_msg % n_part1 = source_msg % n_part1
  dest_msg % n_part2 = source_msg % n_part2
  dest_msg % size1 = source_msg % size1
  dest_msg % size2 = source_msg % size2
  dest_msg % max_buffer1_size = source_msg % max_buffer1_size

  dest_msg % buffer1 => source_msg % buffer1
  dest_msg % buffer2 => source_msg % buffer2

end subroutine copy_spec_msg
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
subroutine check_buffer_size_vdf( msg, recv_vdf, send_vdf )

  implicit none

  class(t_bnd_msg), intent(in) :: msg
  type(t_vdf_msg), intent(inout) :: recv_vdf, send_vdf

  integer :: bsize, msg_size

  ! Check if recv buffer is large enough
  msg_size = msg%msg_size(p_recv_msg,msg%n_tils+1)

  ! (* debug *)
  if ( msg_size == 0 ) then
    write(err_buf__,*) 'Invalid message size, aborting.';call err__("tiles/os-group-boundary.f03",629)
    call abort_program( p_err_invalid )
  endif

  ! grow message buffer if necessary
  if ( associated( recv_vdf%buffer ) ) then
    bsize = size( recv_vdf%buffer )
  else
    bsize = 0
  endif

  if ( msg_size > bsize ) then
    if ( bsize > 0 ) then
      call freemem(recv_vdf%buffer,"tiles/os-group-boundary.f03",642)
    endif
    call alloc(recv_vdf%buffer, (/msg_size/),"tiles/os-group-boundary.f03",644)
  endif

  ! Set message size
  recv_vdf%msg_size = msg_size

  ! Check if send buffer is large enough
  msg_size = msg%msg_size(p_send_msg,msg%n_tils+1)

  ! (* debug *)
  if ( msg_size == 0 ) then
    write(err_buf__,*) 'Invalid message size, aborting.';call err__("tiles/os-group-boundary.f03",655)
    call abort_program( p_err_invalid )
  endif

  ! grow message buffer if necessary
  if ( associated( send_vdf%buffer ) ) then
    bsize = size( send_vdf%buffer )
  else
    bsize = 0
  endif
  if ( msg_size > bsize ) then
    if ( bsize > 0 ) then
      call freemem(send_vdf%buffer,"tiles/os-group-boundary.f03",667)
    endif
    call alloc(send_vdf%buffer, (/msg_size/),"tiles/os-group-boundary.f03",669)
  endif

  ! Set message size
  send_vdf%msg_size = msg_size

end subroutine check_buffer_size_vdf
!-------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Determine how to pack msgs into the two t_spec_msg buffers
!-----------------------------------------------------------------------------------------
subroutine check_buffer_size_part( msg, send_spec, spec )

  implicit none

  class(t_bnd_msg), intent(in) :: msg
  type( t_spec_msg ), intent(inout), pointer :: send_spec
  class(t_species), intent(inout) :: spec

  integer :: j, msg_size_j, bsize

  ! initialize send_spec member data
  call set_spec_msg( send_spec, spec )

  ! get size for first message
  do j = 1, msg%n_tils
    ! each message contains 2 ints:
    ! one for the g_til_aid of tile it's going to, one for number of particles
    msg_size_j = 2*int_pack_size + msg%msg_size(p_send_msg,j)*send_spec%particle_size

    ! move on to second buffer if it will overflow
    if ( send_spec%size1+msg_size_j > send_spec%max_buffer1_size ) exit

    send_spec%size1 = send_spec%size1 + msg_size_j
    ! incremement number of tils in first msg
    send_spec%n_tils1 = send_spec%n_tils1 + 1
  enddo

  ! get size for second message
  do j = send_spec%n_tils1+1, msg%n_tils
    ! each message contains 2 ints:
    ! one for the g_til_aid of tile it's going to, one for number of particles
    msg_size_j = 2*int_pack_size + msg%msg_size(p_send_msg,j)*send_spec%particle_size

    ! put it into the second buffer
    send_spec%size2 = send_spec%size2 + msg_size_j

    ! set n_part2 > 0 so that we can use isend_2_msg in os-spec-comm
    send_spec%n_part2 = 1
  enddo

  ! Check if message buffer1 is big enough, resize if necessary
  if ( associated( send_spec%buffer1 ) ) then
    bsize = size( send_spec%buffer1 )
  else
    bsize = 0
  endif

  if ( bsize < send_spec%size1 ) then
    write(err_buf__,*) 'buffer1 should never be larger than max_buffer1_size';call err__("tiles/os-group-boundary.f03",729)
    call freemem(send_spec%buffer1,"tiles/os-group-boundary.f03",730)
    call alloc(send_spec%buffer1, (/ send_spec%size1 /),"tiles/os-group-boundary.f03",731)
  endif

  ! If 2nd msg necessary
  if ( send_spec%size2 > 0 ) then
    ! Check if message buffer2 is big enough, resize if necessary
    if ( associated( send_spec%buffer2 ) ) then
      bsize = size( send_spec%buffer2 )
    else
      bsize = 0
    endif

    if ( bsize < send_spec%size2 ) then
      call freemem(send_spec%buffer2,"tiles/os-group-boundary.f03",744)
      call alloc(send_spec%buffer2, (/ send_spec%size2 /),"tiles/os-group-boundary.f03",745)
    endif
  endif

end subroutine check_buffer_size_part
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Check node cross for this species for all tiles
! We call this node cross because we reuse the old check_node_cross routines, but
! with tiles, this really just amounts to checking which particles have left a given tile
! Determination of which crossings are local to a node, and which are external
! is done later in update_boundary_part_group
!-----------------------------------------------------------------------------------------
subroutine check_node_cross_group( sim, spec )

  implicit none

  class(t_simulation_group), intent(inout) :: sim
  type(t_spec_arr), dimension(sim%til_id(1):sim%til_id(2)), intent(inout) :: spec

  class(t_simulation), pointer :: til
  class(t_node_conf_tiles), pointer :: no_co_tiles
  integer :: j, t

  select type(no_co => sim%no_co); class is(t_node_conf_tiles)
    no_co_tiles => no_co
    ! do omp over light tiles
    if (no_co_tiles%omp_patt%n_light > 0) then
      !$omp parallel do schedule(dynamic) private(t, til)
      do j = 1, no_co_tiles%omp_patt%n_light
        t = no_co_tiles%omp_patt%light_tils(j)
        til => sim%tiles(t)%t
        call check_node_cross( spec(t)%s, n_threads(til%no_co), til%bnd%node_cross )
        spec(t)%s%num_par_save = spec(t)%s%num_par
      enddo
      !$omp end parallel do
    endif

    ! do omp inside heavy tiles
    do j = 1, no_co_tiles%omp_patt%n_heavy
      t = no_co_tiles%omp_patt%heavy_tils(j)
      til => sim%tiles(t)%t
      call check_node_cross( spec(t)%s, n_threads(til%no_co), til%bnd%node_cross )
      spec(t)%s%num_par_save = spec(t)%s%num_par
    enddo
  end select

end subroutine check_node_cross_group
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine check_boundary_group( sim, spec, dim )

  implicit none

  class(t_simulation_group), intent(inout) :: sim
  integer, intent(in) :: dim
  type(t_spec_arr), dimension(sim%til_id(1):sim%til_id(2)), intent(inout) :: spec

  class(t_simulation), pointer :: til
  class(t_node_conf_tiles), pointer :: no_co_tiles
  integer :: min_npar, j, t

  select type(no_co => sim%no_co); class is(t_node_conf_tiles)
    no_co_tiles => no_co

    ! do omp over all tiles, heavy tiles first to get them started
    !$omp parallel do schedule(dynamic) private(t, til, min_npar)
    do j = 1, sim%til_id(3)
      if ( j <= no_co_tiles%omp_patt%n_heavy ) then
        ! process heavy tiles first
        t = no_co_tiles%omp_patt%heavy_tils(j)
      else
        ! process light tiles last
        t = no_co_tiles%omp_patt%light_tils( j - no_co_tiles%omp_patt%n_heavy )
      endif
      til => sim%tiles(t)%t

      ! check which particles crossed node bndry and store their indexes
      min_npar = spec(t)%s%num_par_save

      if ( spec(t)%s%num_par < min_npar ) then
        spec(t)%s%num_par_save = spec(t)%s%num_par
        min_npar = spec(t)%s%num_par_save
      endif

      call check_boundary( spec(t)%s, min_npar, dim, til%bnd%bnd_cross, til%bnd%node_cross )
    enddo
    !$omp end parallel do

  end select

end subroutine check_boundary_group
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Local update boundary routine (for comms btwn tiles on same node)
! Pack local communication buffers
!-----------------------------------------------------------------------------------------
subroutine loc_update_bnd_spec_1_group( sim, spec, comms_int, dim )

  implicit none

  class(t_simulation_group), intent(inout) :: sim
  integer, intent(in) :: dim
  type(t_spec_arr), dimension(sim%til_id(1):sim%til_id(2)), intent(inout) :: spec
  integer, dimension(sim%til_id(1):sim%til_id(2)), intent(inout) :: comms_int

  class(t_simulation), pointer :: til
  class(t_node_conf_tiles), pointer :: no_co_tiles
  integer :: g_til_dim, l_node, u_node, my_node, t, j

  my_node = sim%no_co%my_aid()

  ! Can't use "select type"/"class is" inside omp loop
  select type(no_co => sim%no_co); class is(t_node_conf_tiles)
    no_co_tiles => no_co
    g_til_dim = no_co_tiles%ti_co%g_til(dim)

    ! Do heavy tiles first to get them started
    !$omp parallel do schedule(dynamic) private(t, til, l_node, u_node)
    do j = 1, sim%til_id(3)

      if ( j <= no_co_tiles%omp_patt%n_heavy ) then
        ! process heavy tiles first
        t = no_co_tiles%omp_patt%heavy_tils(j)
      else
        ! process light tiles last
        t = no_co_tiles%omp_patt%light_tils( j - no_co_tiles%omp_patt%n_heavy )
      endif

      til => sim%tiles(t)%t

      l_node = neighbor( til%no_co, p_lower, dim )
      u_node = neighbor( til%no_co, p_upper, dim )

      ! result is setting comms_int
      comms_int(t) = 0
      if ( l_node==my_node ) then
        if (g_til_dim==1) then
          comms_int(t) = comms_int(t) + p_lower_phys
        else
          comms_int(t) = comms_int(t) + p_lower
        endif
      elseif ( l_node < 0 ) then
        comms_int(t) = comms_int(t) + p_lower_phys
      endif

      if ( u_node==my_node ) then
        if (g_til_dim==1) then
          comms_int(t) = comms_int(t) + p_upper_phys
        else
          comms_int(t) = comms_int(t) + p_upper
        endif
      elseif ( u_node < 0 ) then
        comms_int(t) = comms_int(t) + p_upper_phys
      endif

      ! critical statements are necessary because the current buffers are module variables
      select case ( comms_int(t) )
      case (p_lower)
        ! prepare message buffer for lower tile in same node
        call loc_send( spec(t)%s, dim, p_lower, til%no_co, til%bnd%bnd_cross(p_lower), til%bnd%send_spec )
      case (p_upper)
        ! prepare message buffer for upper tile in same node
        call loc_send( spec(t)%s, dim, p_upper, til%no_co, til%bnd%bnd_cross( p_upper ), til%bnd%send_spec )
      case (p_lower+p_upper)
        ! prepare message buffer for upper and lower tiles in same node
        call loc_send( spec(t)%s, dim, p_upper, til%no_co, til%bnd%bnd_cross( p_upper ), til%bnd%send_spec )
        call loc_send( spec(t)%s, dim, p_lower, til%no_co, til%bnd%bnd_cross( p_lower ), til%bnd%send_spec )
      case (p_lower_phys)
        ! process lower physical boundary, no mpi communication with upper node
        !$omp critical
        call process_phys_boundary( spec(t)%s, til%jay, spec(t)%s%dt, dim, p_lower, til%bnd%bnd_cross(p_lower) )
        !$omp end critical
      case (p_lower_phys+p_upper)
        ! prepare message buffer for upper tile in same node
        call loc_send( spec(t)%s, dim, p_upper, til%no_co, til%bnd%bnd_cross( p_upper ), til%bnd%send_spec )
        ! process lower physical boundary
        !$omp critical
        call process_phys_boundary( spec(t)%s, til%jay, spec(t)%s%dt, dim, p_lower, til%bnd%bnd_cross(p_lower) )
        !$omp end critical
      case (p_upper_phys)
        ! process upper physical boundary, no mpi communication with lower node
        !$omp critical
        call process_phys_boundary( spec(t)%s, til%jay, spec(t)%s%dt, dim, p_upper, til%bnd%bnd_cross(p_upper) )
        !$omp end critical
      case (p_lower+p_upper_phys)
        ! prepare message buffer for lower tile in same node
        call loc_send( spec(t)%s, dim, p_lower, til%no_co, til%bnd%bnd_cross(p_lower), til%bnd%send_spec )
        ! process upper physical boundary
        !$omp critical
        call process_phys_boundary( spec(t)%s, til%jay, spec(t)%s%dt, dim, p_upper, til%bnd%bnd_cross(p_upper) )
        !$omp end critical
      case (p_lower_phys+p_upper_phys)
        ! process both lower and upper physical boundaries
        !$omp critical
        call process_phys_boundary( spec(t)%s, til%jay, spec(t)%s%dt, dim, p_lower, til%bnd%bnd_cross(p_lower) )
        call process_phys_boundary( spec(t)%s, til%jay, spec(t)%s%dt, dim, p_upper, til%bnd%bnd_cross(p_upper) )
        !$omp end critical
      end select

    enddo
    !$omp end parallel do

  end select

end subroutine loc_update_bnd_spec_1_group
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Local update boundary routine (for comms btwn tiles on same node)
! Gather and unpack local communication buffers
!-----------------------------------------------------------------------------------------
subroutine loc_update_bnd_spec_2_group( sim, spec, comms_int, dim )

  implicit none

  class(t_simulation_group), intent(inout) :: sim
  integer, intent(in) :: dim
  type(t_spec_arr), dimension(sim%til_id(1):sim%til_id(2)), intent(inout) :: spec
  integer, dimension(sim%til_id(1):sim%til_id(2)), intent(inout) :: comms_int

  type(t_spec_msg), dimension(2) :: loc_msg_buff
  class(t_simulation), pointer :: til, l_til, u_til
  class(t_node_conf), pointer :: no_co

  integer, dimension(2,sim%til_id(1):sim%til_id(2)) :: nbr_til
  integer :: t

  ! need nbr_til as an array because can't use "select type"/"class is" inside omp loop
  do t = sim%til_id(1), sim%til_id(2)
    no_co => sim%tiles(t)%t%no_co
    select type(no_co); class is(t_node_conf_tiles)
      nbr_til(:,t) = no_co%ti_co%neighbor_g_til(:,dim)
    end select
  enddo

  !$omp parallel do schedule(dynamic) private(loc_msg_buff, til, l_til, u_til)
  do t = sim%til_id(1), sim%til_id(2)

    til => sim%tiles(t)%t

    ! copy message data to local message buffer
    select case (comms_int(t))
    case (p_lower,p_lower+p_upper_phys)
      ! nbr_til = no_co%ti_co%neighbor_g_til(p_lower,dim)
      l_til => sim%tiles(nbr_til(p_lower,t))%t
      call copy_spec_msg( l_til%bnd%send_spec(p_upper), loc_msg_buff(p_lower) )
    case (p_upper,p_lower_phys+p_upper)
      ! nbr_til = no_co%ti_co%neighbor_g_til(p_upper,dim)
      u_til => sim%tiles(nbr_til(p_upper,t))%t
      call copy_spec_msg( u_til%bnd%send_spec(p_lower), loc_msg_buff(p_upper) )
    case (p_lower+p_upper)
      ! nbr_til = no_co%ti_co%neighbor_g_til(p_lower,dim)
      l_til => sim%tiles(nbr_til(p_lower,t))%t
      call copy_spec_msg( l_til%bnd%send_spec(p_upper), loc_msg_buff(p_lower) )
      ! nbr_til = no_co%ti_co%neighbor_g_til(p_upper,dim)
      u_til => sim%tiles(nbr_til(p_upper,t))%t
      call copy_spec_msg( u_til%bnd%send_spec(p_lower), loc_msg_buff(p_upper) )
    end select

    ! receive and unpack local message buffer
    select case ( comms_int(t) )
    case (p_lower,p_lower+p_upper_phys)
      call loc_recv( spec(t)%s, dim, p_lower, til%no_co, loc_msg_buff, til%bnd%bnd_cross )
    case (p_upper,p_lower_phys+p_upper)
      call loc_recv( spec(t)%s, dim, p_upper, til%no_co, loc_msg_buff, til%bnd%bnd_cross )
    case (p_lower + p_upper)
      call loc_recv( spec(t)%s, dim, p_lower, til%no_co, loc_msg_buff, til%bnd%bnd_cross )
      call loc_recv( spec(t)%s, dim, p_upper, til%no_co, loc_msg_buff, til%bnd%bnd_cross )
    end select

  enddo
  !$omp end parallel do

end subroutine loc_update_bnd_spec_2_group
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! External update boundary routine for species/particles (for comms btwn nodes)
! Pack buffers for external comms, send buff1, post recv for buff1
!-----------------------------------------------------------------------------------------
subroutine ext_update_bnd_spec_1_group( sim, spec, dim )

  implicit none

  class(t_simulation_group), intent(inout) :: sim
  integer, intent(in) :: dim
  type(t_spec_arr), dimension(sim%til_id(1):sim%til_id(2)), intent(inout) :: spec

  type(t_spec_msg), pointer :: snd_spc, rcv_spc
  class(t_bnd_patt), pointer :: patt
  class(t_bnd_msg), pointer :: msg
  integer :: i

  patt => sim%bnd_patt(dim)
  snd_spc => sim%bnd_grp%send_spec
  rcv_spc => sim%bnd_grp%recv_spec

  ! loop over msgs (each msg corresponds to a single node)
  do i = 1, patt%n_msg

    msg => patt%msg(i)

    ! zero out msg size
    msg%msg_size(:,:) = 0

    ! Post receive
    call irecv_1( spec(sim%til_id(1))%s, msg%node, sim%no_co, rcv_spc )

    ! count size of each msg (only need to count send size, recv size fixed)
    call isend_size_spec( sim, msg )

    ! check buff size and reallocate if necessary
    call check_buffer_size( msg, snd_spc, spec(sim%til_id(1))%s )

    ! pack data and post send
    call isend_1( sim, spec, snd_spc, msg, dim )

    ! point to next element in linked list
    snd_spc => snd_spc%next
    rcv_spc => rcv_spc%next

  enddo

end subroutine ext_update_bnd_spec_1_group
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! External update boundary routine for species/particles (for comms btwn nodes)
! Send/post recv for buff2 if necessary
!-----------------------------------------------------------------------------------------
subroutine ext_update_bnd_spec_2_group( sim, dim )

  implicit none

  class(t_simulation_group), intent(inout) :: sim
  integer, intent(in) :: dim

  type(t_spec_msg), pointer :: snd_spc, rcv_spc
  class(t_bnd_patt), pointer :: patt
  class(t_bnd_msg), pointer :: msg
  integer :: i

  patt => sim%bnd_patt(dim)
  snd_spc => sim%bnd_grp%send_spec
  rcv_spc => sim%bnd_grp%recv_spec

  ! loop over msgs (each msg corresponds to a single node)
  do i = 1, patt%n_msg

    msg => patt%msg(i)

    ! wait for 1st message, then post recv for 2nd msg (if necessary)
    ! dummy integer for interface purposes
    call irecv_2( rcv_spc, 0 )

    ! post send for second msg (if necessary)
    call isend_2( snd_spc )

    ! point to next element of linked list
    snd_spc => snd_spc%next
    rcv_spc => rcv_spc%next

  enddo

end subroutine ext_update_bnd_spec_2_group
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! External update boundary routine for species/particles (for comms btwn nodes)
! Unpack msgs and wait for all msgs to send
!-----------------------------------------------------------------------------------------
subroutine ext_update_bnd_spec_3_group( sim, spec, dim )

  implicit none

  class(t_simulation_group), intent(inout) :: sim
  integer, intent(in) :: dim
  type(t_spec_arr), dimension(sim%til_id(1):sim%til_id(2)), intent(inout) :: spec

  type(t_spec_msg), pointer :: snd_spc, rcv_spc
  class(t_bnd_patt), pointer :: patt
  class(t_bnd_msg), pointer :: msg
  integer :: i

  patt => sim%bnd_patt(dim)
  rcv_spc => sim%bnd_grp%recv_spec

  ! Unpack messages
  do i = 1, patt%n_msg
    msg => patt%msg(i)
    call unpack( sim, spec, rcv_spc, msg%n_tils )
    rcv_spc => rcv_spc%next
  enddo

  ! Wait for messages to send
  snd_spc => sim%bnd_grp%send_spec
  do i = 1, patt%n_msg
    call wait_spec_msg( snd_spc )
    snd_spc => snd_spc%next
  enddo

end subroutine ext_update_bnd_spec_3_group
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! remove particles whose index is still in par_idx_l or par_idx_u
! and pack the spec(t)%s buffer
!-----------------------------------------------------------------------------------------
subroutine remove_particles_group( sim, spec )

  implicit none

  class(t_simulation_group), intent(inout) :: sim
  type(t_spec_arr), dimension(sim%til_id(1):sim%til_id(2)), intent(inout) :: spec

  integer :: j, t
  class(t_node_conf_tiles), pointer :: no_co_tiles

  select type(no_co => sim%no_co); class is(t_node_conf_tiles)
    no_co_tiles => no_co

    ! do omp over all tiles, heavy tiles first to get them started
    !$omp parallel do schedule(dynamic) private(t)
    do j = 1, sim%til_id(3)
      if ( j <= no_co_tiles%omp_patt%n_heavy ) then
        ! process heavy tiles first
        t = no_co_tiles%omp_patt%heavy_tils(j)
      else
        ! process light tiles last
        t = no_co_tiles%omp_patt%light_tils( j - no_co_tiles%omp_patt%n_heavy )
      endif

      call remove_particles( spec(t)%s, sim%tiles(t)%t%bnd%bnd_cross )
    enddo
    !$omp end parallel do

  end select

end subroutine remove_particles_group
!-----------------------------------------------------------------------------------------

end module m_group_boundary


!-------------------------------------------------------------------------------
! updates emf at boundaries
!-------------------------------------------------------------------------------
subroutine update_boundary_emf_group( sim )
!-------------------------------------------------------------------------------

  use m_parameters
  use m_simulation_group, only : t_simulation_group, emfboundev_glob, emfboundev_loc
  use m_group_boundary, only : p_wall_e, p_wall_b, update_boundary, check_buffer_size
  use m_group_boundary, only : send_recv, wait_all, loc_update_boundary
  use m_logprof, only : begin_event, end_event
  use m_space, only : x_dim, nx_move, if_move
  use m_node_conf, only : t_node_conf, neighbor
  use m_node_conf_tiles, only : t_node_conf_tiles
  use m_vdf_define, only : t_vdf_arr
  use m_vdf_comm, only : update_periodic_vdf
  use m_vdf_comm_tiles, only : p_msg, p_wait
  use m_emf_bound, only : emfboundev


  use m_vpml, only : is_active


  implicit none

  class(t_simulation_group), intent(inout) :: sim

  integer :: i, j, my_node
  integer, dimension(p_max_dim) :: lnx_move
  logical, dimension(p_max_dim) :: lif_move
  integer :: n_vpml_max, i_vpml, dim, eb
  class(t_node_conf), pointer :: no_co
  integer, dimension(p_x_dim) :: g_til
  type(t_vdf_arr), dimension(:), pointer :: e_arr, b_arr

  call begin_event(emfboundev)

  lnx_move( 1 : x_dim( sim%g_space ) ) = nx_move( sim%g_space )
  lif_move( 1 : x_dim( sim%g_space ) ) = if_move( sim%g_space )

  ! TODO: shared-memoryize
  ! update to the lindman boundaries must occur before updating em fields to
  ! get the corner value right
  do i = p_msg, p_wait
    do j = sim%til_id(1), sim%til_id(2)
      call update_boundary(sim%tiles(j)%t%emf%bnd_con%lindman(p_lower), &
        sim%tiles(j)%t%no_co, lnx_move, sim%tiles(j)%t%bnd%send_vdf, &
        sim%tiles(j)%t%bnd%recv_vdf, i )
      call update_boundary(sim%tiles(j)%t%emf%bnd_con%lindman(p_upper), &
        sim%tiles(j)%t%no_co, lnx_move, sim%tiles(j)%t%bnd%send_vdf, &
        sim%tiles(j)%t%bnd%recv_vdf, i )
    enddo
  enddo



  ! Update VPML walls if necessary
  n_vpml_max = 0
  do j = sim%til_id(1), sim%til_id(2)
    n_vpml_max = max(n_vpml_max,sim%tiles(j)%t%emf%bnd_con%vpml_all%n_vpml)
  enddo

  do dim = 1, p_x_dim

    ! Loop over e walls first, then b walls
    do eb = p_wall_e, p_wall_b

      ! Set current vpml to 1 (will be ignored for those with n_vpml==0)
      do j = sim%til_id(1), sim%til_id(2)
        sim%tiles(j)%t%emf%bnd_con%vpml_all%c = 1
      enddo

      ! Make sure to loop through at least the maximum number of vpml
      do i_vpml = 1, n_vpml_max
        do i = p_msg, p_wait
          do j = sim%til_id(1), sim%til_id(2)
            if ( is_active(sim%tiles(j)%t%emf%bnd_con%vpml_all) ) then
              call update_boundary(sim%tiles(j)%t%emf%bnd_con%vpml_all, &
                sim%tiles(j)%t%no_co, lif_move, lnx_move, sim%tiles(j)%t%bnd%send_vdf, &
                sim%tiles(j)%t%bnd%recv_vdf, i, dim, eb, sim )
            endif
          enddo
        enddo
      enddo

    enddo

  enddo



  my_node = sim%no_co%my_aid()
  no_co => sim%tiles(sim%til_id(1))%t%no_co;
  select type(no_co); class is(t_node_conf_tiles)
    g_til = no_co%ti_co%g_til
  end select

  allocate( e_arr(sim%til_id(1):sim%til_id(2)), b_arr(sim%til_id(1):sim%til_id(2)) )
  do i = sim%til_id(1), sim%til_id(2)
    e_arr(i)%v => sim%tiles(i)%t%emf%e
    b_arr(i)%v => sim%tiles(i)%t%emf%b
  enddo

  ! Update other boundary types
  do dim = 1, p_x_dim

    ! Check if single-node and single-tile periodic
    if ( g_til(dim) == 1 .and. &
      my_node == neighbor(sim%tiles(sim%til_id(1))%t%no_co,p_lower,dim) ) then

      do j = sim%til_id(1), sim%til_id(2)
        call update_periodic_vdf( e_arr(j)%v, dim, p_vdf_replace )
        call update_periodic_vdf( b_arr(j)%v, dim, p_vdf_replace )
      enddo

    else

      call begin_event(emfboundev_glob)

      if ( sim%bnd_patt(dim)%n_msg>0 ) then
        call send_recv( e_arr, p_vdf_replace, dim, lnx_move(dim), sim%til_id(1), &
                        sim%no_co, sim%bnd_patt(dim), sim%bnd_grp%send_vdf, &
                        sim%bnd_grp%recv_vdf, b_arr )
      endif

      call end_event(emfboundev_glob)
      call begin_event(emfboundev_loc)

      call loc_update_boundary( e_arr, p_vdf_replace, dim, lnx_move(dim), sim%til_id(1), &
                                sim%bnd_patt(dim)%msg_loc, b_arr )

      call end_event(emfboundev_loc)
      call begin_event(emfboundev_glob)

      if ( sim%bnd_patt(dim)%n_msg>0 ) then
        call wait_all( e_arr, sim%til_id(1), sim%bnd_patt(dim), sim%bnd_grp%send_vdf, &
                        sim%bnd_grp%recv_vdf, b_arr )
      endif

      call end_event(emfboundev_glob)

    endif

  enddo

  deallocate( e_arr, b_arr )

  call end_event(emfboundev)

end subroutine update_boundary_emf_group
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
! updates the particle data at the boundaries
!-------------------------------------------------------------------------------
subroutine update_boundary_part_group( sim )

  use m_parameters
  use m_particles_define, only : partboundev
  use m_simulation_group, only : t_simulation_group, partboundev_glob, partboundev_loc
  use m_group_boundary, only : check_ncross_grp, check_bnd_grp, ext_ub_spec_1, ext_ub_spec_2
  use m_group_boundary, only : ext_ub_spec_3, loc_ub_spec_1, loc_ub_spec_2, rmv_part_grp
  use m_group_boundary, only : update_boundary
  use m_logprof, only : begin_event, end_event
  use m_species_define, only : t_spec_arr
  use m_vdf_define, only : t_vdf_arr
  use m_space, only : nx_move

  implicit none

  class(t_simulation_group), intent(inout) :: sim

  type(t_spec_arr), dimension(sim%til_id(1):sim%til_id(2)) :: spec
  integer, dimension(sim%til_id(1):sim%til_id(2)) :: comms_int
  type(t_vdf_arr), dimension(:), pointer :: multi_ion_arr
  integer :: t, dim, i

  call begin_event(partboundev)

  do t = sim%til_id(1), sim%til_id(2)
    spec(t)%s => sim%tiles(t)%t%part%species
  enddo

  do
    if (.not. associated(spec(sim%til_id(1))%s)) exit

    ! Check node cross for group of tiles
    call check_ncross_grp( sim, spec )

    do dim = 1, p_x_dim

      ! Check bnd cross for group of tiles
      call check_bnd_grp( sim, spec, dim )

      call begin_event(partboundev_glob)
      ! Pack buffers for external comms, send buff1, post recv for buff1
      call ext_ub_spec_1( sim, spec, dim )
      call end_event(partboundev_glob)

      call begin_event(partboundev_loc)
      ! Pack local communication buffers
      call loc_ub_spec_1( sim, spec, comms_int, dim )

      ! Gather and unpack local communication buffers
      call loc_ub_spec_2( sim, spec, comms_int, dim )
      call end_event(partboundev_loc)

      call begin_event(partboundev_glob)
      ! Send/post recv for buff2 if necessary
      call ext_ub_spec_2( sim, dim )

      ! Unpack msgs and wait for all msgs to send
      call ext_ub_spec_3( sim, spec, dim )
      call end_event(partboundev_glob)

      ! Remove particles no longer on this node
      call rmv_part_grp( sim, spec )

    enddo

    ! advance species pointers
    do t = sim%til_id(1),sim%til_id(2)
      spec(t)%s => spec(t)%s%next
    enddo

  enddo



  allocate( multi_ion_arr(sim%til_id(1):sim%til_id(2)) )

  do i = 1, sim%tiles(sim%til_id(1))%t%part%num_neutral
    do t = sim%til_id(1), sim%til_id(2)
      multi_ion_arr(t)%v => sim%tiles(t)%t%part%neutral(i)%multi_ion
    enddo

    call update_boundary( multi_ion_arr, p_vdf_replace, sim%til_id(1), sim%no_co, &
                          sim%tiles(sim%til_id(1))%t%no_co, sim%bnd_patt, &
                          sim%bnd_grp%send_vdf, sim%bnd_grp%recv_vdf, &
                          nx_move_ = nx_move(sim%g_space) )
  enddo

  deallocate( multi_ion_arr )



  call end_event(partboundev)

end subroutine update_boundary_part_group
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
! Update boundary of electrical current
!-------------------------------------------------------------------------------
subroutine update_boundary_current_group( sim )

  use m_parameters
  use m_node_conf, only : t_node_conf, neighbor
  use m_node_conf_tiles, only : t_node_conf_tiles
  use m_logprof, only : begin_event, end_event
  use m_simulation_group, only : t_simulation_group, jayboundev_glob, jayboundev_loc
  use m_simulation, only : t_simulation
  use m_group_boundary, only : check_buffer_size, send_recv, wait_all, loc_update_boundary
  use m_vdf_define, only : t_vdf_arr
  use m_vdf_comm, only : update_periodic_vdf, t_vdf_msg
  use m_current_boundary, only : spec_bc_1d, spec_bc_2d, spec_bc_3d
  use m_current_define, only : jayboundev

  implicit none

  class(t_simulation_group), intent(inout) :: sim

  class(t_simulation), pointer :: til
  integer :: i, j, dim
  integer :: my_node
  class(t_node_conf), pointer :: no_co
  integer, dimension(p_x_dim) :: g_til
  type(t_vdf_arr), dimension(:), pointer :: jay_arr

  call begin_event( jayboundev )

  ! it is not necessary to use move_num(space) since the values for the current
  ! are recalculated at each time-step so no data needs to be shifted

  my_node = sim%no_co%my_aid()
  no_co => sim%tiles(sim%til_id(1))%t%no_co;
  select type(no_co); class is(t_node_conf_tiles)
    g_til = no_co%ti_co%g_til
  end select

  allocate( jay_arr(sim%til_id(1):sim%til_id(2)) )
  do i = sim%til_id(1), sim%til_id(2)
    jay_arr(i)%v => sim%tiles(i)%t%jay%pf(1)
  enddo

  ! Update other boundary types
  do dim = 1, p_x_dim

    ! Check if single-node and single-tile periodic
    if ( g_til(dim) == 1 .and. &
      my_node == neighbor(sim%tiles(sim%til_id(1))%t%no_co,p_lower,dim) ) then

      do j = sim%til_id(1), sim%til_id(2)
        call update_periodic_vdf( jay_arr(j)%v, dim, p_vdf_add )
      enddo

    else

      call begin_event(jayboundev_glob)

      if ( sim%bnd_patt(dim)%n_msg>0 ) then
        call send_recv( jay_arr, p_vdf_add, dim, 0, sim%til_id(1), sim%no_co, &
                        sim%bnd_patt(dim), sim%bnd_grp%send_vdf, sim%bnd_grp%recv_vdf )
      endif

      call end_event(jayboundev_glob)
      call begin_event(jayboundev_loc)

      call loc_update_boundary( jay_arr, p_vdf_add, dim, 0, sim%til_id(1), &
                                sim%bnd_patt(dim)%msg_loc )

      call end_event(jayboundev_loc)
      call begin_event(jayboundev_glob)

      if ( sim%bnd_patt(dim)%n_msg>0 ) then
        call wait_all( jay_arr, sim%til_id(1), sim%bnd_patt(dim), sim%bnd_grp%send_vdf, &
                        sim%bnd_grp%recv_vdf )
      endif

      call end_event(jayboundev_glob)

    endif

  enddo

  deallocate( jay_arr )


  ! Process specular boundaries
  do j = sim%til_id(1), sim%til_id(2)
    til => sim%tiles(j)%t

    select case ( p_x_dim )
      case (1)
        if ( til%jay%bc_type( p_lower, 1 ) == p_bc_specular ) then
          call spec_bc_1d( til%jay%pf(1), p_lower, til%jay%interpolation )
        endif
        if ( til%jay%bc_type( p_upper, 1 ) == p_bc_specular ) then
          call spec_bc_1d( til%jay%pf(1), p_upper, til%jay%interpolation )
        endif
      case (2)
        do i = 1, 2
          if ( til%jay%bc_type( p_lower, i ) == p_bc_specular ) then
            call spec_bc_2d( til%jay%pf(1), i, p_lower, til%jay%interpolation )
          endif
          if ( til%jay%bc_type( p_upper, i ) == p_bc_specular ) then
            call spec_bc_2d( til%jay%pf(1), i, p_upper, til%jay%interpolation )
          endif
        enddo
      case (3)
        do i = 1, 3
          if ( til%jay%bc_type( p_lower, i ) == p_bc_specular ) then
            call spec_bc_3d( til%jay%pf(1), i, p_lower, til%jay%interpolation )
          endif
          if ( til%jay%bc_type( p_upper, i ) == p_bc_specular ) then
            call spec_bc_3d( til%jay%pf(1), i, p_upper, til%jay%interpolation )
          endif
        enddo
    end select

  enddo


  call end_event( jayboundev )


end subroutine update_boundary_current_group
!-------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Update the global charge
!-----------------------------------------------------------------------------------------
subroutine update_charge_part_group( sim )

  use m_parameters
  use m_simulation_group, only : t_simulation_group
  use m_particles_define, only : t_particles
  use m_simulation, only : t_simulation
  use m_vdf_define, only : t_vdf_arr
  use m_group_comm, only : deposit_charge_1, deposit_charge_2
  use m_group_boundary, only : update_boundary

  implicit none

  class(t_simulation_group), intent(inout) :: sim

  class( t_particles ), pointer :: particles
  class(t_simulation), pointer :: til
  integer :: i
  type(t_vdf_arr), dimension(:), pointer :: vdf_arr

  allocate( vdf_arr(sim%til_id(1):sim%til_id(2)) )

  ! first half of update charge
  do i = sim%til_id(1), sim%til_id(2)
    til => sim%tiles(i)%t
    particles => til%part

    ! If charge is already up to date return silently
    if ( particles%n_current == particles%charge%n_last_update ) return

    ! If keeping the previous charge rotate the pointers
    if ( particles%charge%keep_previous ) then

      if ( associated( particles%charge%current, particles%charge%charge0 ) ) then
        particles%charge%current => particles%charge%charge1
        particles%charge%previous => particles%charge%charge0
      else
        particles%charge%current => particles%charge%charge0
        particles%charge%previous => particles%charge%charge1
      endif

      particles%charge%n_prev_update = particles%charge%n_last_update

    else

      particles%charge%current => particles%charge%charge0

    endif

    ! Deposit current
    call deposit_charge_1( particles, sim%g_space, til%grid, particles%charge%current )

    vdf_arr(i)%v => particles%charge%current
  enddo

  call update_boundary( vdf_arr, p_vdf_add, sim%til_id(1), sim%no_co, &
    sim%tiles(sim%til_id(1))%t%no_co, sim%bnd_patt, sim%bnd_grp%send_vdf, &
    sim%bnd_grp%recv_vdf )
  deallocate( vdf_arr )

  ! second half of update charge
  do i = sim%til_id(1), sim%til_id(2)
    til => sim%tiles(i)%t
    particles => til%part

    ! Deposit current
    call deposit_charge_2( particles, til%grid, particles%charge%current )

    ! Update the last deposit counter
    particles%charge%n_last_update = particles%n_current
  enddo

end subroutine update_charge_part_group
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Gathers information for update boundary (fills t_bnd_patt)
!-----------------------------------------------------------------------------------------
subroutine new_bnd_patt( this )

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
# 1658 "tiles/os-group-boundary.f03" 2
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
# 1659 "tiles/os-group-boundary.f03" 2
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
# 1660 "tiles/os-group-boundary.f03" 2

  use m_parameters
  use m_simulation_group, only : t_simulation_group
  use m_node_conf, only : t_node_conf, neighbor
  use m_node_conf_tiles, only : t_node_conf_tiles
  use m_grid_tiles, only : unique_count_1d

  implicit none

  class( t_simulation_group ), intent(inout) :: this

  ! nbrs(1:4,i)=(node, sending til, upper/lower, receiving til)
  integer, dimension(:,:), pointer :: nbrs
  integer, dimension(:,:), pointer :: unique
  integer :: my_node, node
  integer :: til, nbr_til, n_tils
  integer :: dim, j, k, c, cc
  class(t_node_conf), pointer :: no_co

  my_node = this%no_co%my_aid()

  call alloc(nbrs, (/4, 2*this%til_id(3)/),"tiles/os-group-boundary.f03",1681)

  do dim = 1, p_x_dim

    ! cleanup old bnd_patt info
    call this%bnd_patt(dim)%cleanup()

    ! ---- Fill this%bnd_patt%msg
    ! find external neighbors in ith dimension of all tiles on my node
    ! (dim.e.not including local neighbors)
    cc = 0
    do til = this%til_id(1), this%til_id(2)
      do k = p_lower, p_upper
        node = neighbor( this%tiles(til)%t%no_co, k, dim )
        if ( node/=my_node .and. node>0 ) then
          cc = cc+1
          no_co => this%tiles(til)%t%no_co;
          select type(no_co); class is(t_node_conf_tiles)
            nbr_til = no_co%ti_co%neighbor_g_til( k, dim )
          end select
          nbrs(:,cc) = (/node, til, k, nbr_til /)
        endif
      enddo
    enddo

    if ( cc==0 ) then
      this%bnd_patt(dim)%n_msg = 0
    else
      ! find unique nodes and how many msgs to each node
      unique => null()
      call unique_count_1d( nbrs(1,1:cc), unique )
      this%bnd_patt(dim)%n_msg = size(unique,1)
      allocate( this%bnd_patt(dim)%msg( this%bnd_patt(dim)%n_msg ) )

      ! fill in msg%node, msg%tils, allocate other member data
      do k = 1, this%bnd_patt(dim)%n_msg

        ! node id for current message
        node = unique(k,1)
        this%bnd_patt(dim)%msg(k)%node = node

        ! number of messages going to that node
        n_tils = unique(k,2)
        this%bnd_patt(dim)%msg(k)%n_tils = n_tils
        call alloc(this%bnd_patt(dim)%msg(k)%tils, (/3,n_tils/),"tiles/os-group-boundary.f03",1725)

        ! allocate arrays for later use
        call alloc(this%bnd_patt(dim)%msg(k)%msg_size, (/2,n_tils+1/),"tiles/os-group-boundary.f03",1728)
        call alloc(this%bnd_patt(dim)%msg(k)%send_vdf_range, (/2,p_max_dim,n_tils/),"tiles/os-group-boundary.f03",1729)
        call alloc(this%bnd_patt(dim)%msg(k)%recv_vdf_range, (/2,p_max_dim,n_tils/),"tiles/os-group-boundary.f03",1730)

        ! find g_til_aid's and pack them into the msg
        ! tiles can be repeated up to two times for p_lower and p_upper
        c = 0
        do j = 1, cc
          if ( nbrs(1,j)==node ) then
            c = c+1
            this%bnd_patt(dim)%msg(k) % tils(1,c) = nbrs(2,j)
            this%bnd_patt(dim)%msg(k) % tils(2,c) = nbrs(3,j)
            this%bnd_patt(dim)%msg(k) % tils(3,c) = nbrs(4,j)
          endif
        enddo
      enddo

      call freemem(unique,"tiles/os-group-boundary.f03",1745)
    endif
    ! ----

    ! ---- fill this%bnd_patt%msg_loc
    ! find total number of local comms
    n_tils = 0
    do til = this%til_id(1), this%til_id(2)
      do j = p_lower, p_upper

        ! check for lower local comm
        node = neighbor( this%tiles(til)%t%no_co, j, dim )
        no_co => this%tiles(til)%t%no_co;
        select type(no_co); class is(t_node_conf_tiles)
          nbr_til = no_co%ti_co%neighbor_g_til( j, dim )
        end select

        ! only process if you have smaller g_til_aid to avoid double counting
        if ( node==my_node .and. til<nbr_til ) then
          n_tils = n_tils+1
        endif

      enddo
    enddo

    ! fill msg_loc
    if ( n_tils==0 ) then
      allocate( this%bnd_patt(dim)%msg_loc )
      this%bnd_patt(dim)%msg_loc%n_tils = 0
    else
      allocate( this%bnd_patt(dim)%msg_loc )
      call alloc(this%bnd_patt(dim)%msg_loc%tils, (/3,n_tils/),"tiles/os-group-boundary.f03",1776)

      this%bnd_patt(dim)%msg_loc%n_tils = n_tils

      n_tils = 0
      do til = this%til_id(1), this%til_id(2)
        do j = p_lower, p_upper

          ! check for lower local comm
          node = neighbor( this%tiles(til)%t%no_co, j, dim )
          no_co => this%tiles(til)%t%no_co;
          select type(no_co); class is(t_node_conf_tiles)
            nbr_til = no_co%ti_co%neighbor_g_til( j, dim )
          end select

          ! only process if you have smaller g_til_aid to avoid double counting
          if ( node==my_node .and. til<nbr_til ) then
            n_tils = n_tils+1
            this%bnd_patt(dim)%msg_loc % tils(1,n_tils) = til
            this%bnd_patt(dim)%msg_loc % tils(2,n_tils) = j
            this%bnd_patt(dim)%msg_loc % tils(3,n_tils) = nbr_til
          endif

        enddo
      enddo

    endif
    ! ----

  enddo

  call freemem(nbrs,"tiles/os-group-boundary.f03",1807)

end subroutine new_bnd_patt
!-----------------------------------------------------------------------------------------
