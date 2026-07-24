#include "os-config.h"
#include "os-preprocess.fpp"

module m_sim_group_cuda

#include "memory/memory.h"

! #define __DEBUG__ 1

use m_parameters
use m_simulation_group
use m_sim_tiles_cuda
use m_node_conf, only: comm
use m_node_conf_tiles
use m_node_conf_cuda
use m_time
use m_logprof
use m_particles_define
use m_species_define
use m_vdf_define
use m_time_step, only : dt, test_if_report, advance
use m_vdf_report, only : if_report
use m_space, only : t_space, nx_move
use m_emf_define, only : t_emf
use m_grid_define, only: t_msg_patt, t_grid_arr
use m_restart, only: if_restart_write
use m_time_step, only: n, ndump

implicit none

private

type, extends( t_simulation_group ) :: t_sim_group_cuda

contains

  procedure :: allocate_objs      => allocate_objs_group_cuda
  procedure :: allocate_tiles     => allocate_tiles_cuda
  procedure :: init               => init_sim_group_cuda
  procedure :: advance_deposit    => advance_deposit_group_cuda
  procedure :: move_part_window   => move_window_part_group_cuda
  procedure :: cleanup            => cleanup_sim_group_cuda

  procedure :: update_boundary_part     => update_boundary_part_group_cuda
  procedure :: update_boundary_part_cpu => update_boundary_part_cpu_group_cuda
  procedure :: report_dry_part          => report_dry_part_group_cuda
  procedure :: report_energy_dry_part   => report_energy_dry_part_group_cuda

  ! Load balance procedures
  procedure :: allocate_tiles_dlb => allocate_tiles_dlb_cuda
  procedure :: send_recv          => send_recv_tils_dlb_cuda
  procedure :: init_dlb           => init_dlb_cuda

end type t_sim_group_cuda

! Add interface so we can call this from subclass without type-casting back to superclass
interface
  subroutine init_sim_group( sim, grid_gc_min )
  import t_simulation_group, p_x_dim
  class( t_simulation_group ), intent(inout) ::  sim
  integer, dimension( 2, p_x_dim ), intent(in) :: grid_gc_min
  end subroutine
end interface

! Add interface so we can call this from subclass without type-casting back to superclass
interface
  subroutine init_dlb_group( sim, patt, tiles_temp, dif, diff, difff )
    import t_simulation_group, t_msg_patt, t_til_container
    class( t_simulation_group ), intent(inout)      :: sim
    type( t_msg_patt ), intent(inout)               :: patt
    class( t_til_container ), dimension(:), pointer :: tiles_temp
    integer, dimension(:), optional                 :: dif
    integer, dimension(:,:), optional               :: diff
    integer, dimension(:,:,:), optional             :: difff
  end subroutine init_dlb_group
end interface

! Interfaces to functions defined in C
interface
  subroutine init_node_conf( n_tils, my_aid, no_num, neighbor_til_id, shift ) bind(C)
    use iso_c_binding
    integer(c_int), value :: n_tils
    integer(c_int), value :: my_aid
    integer(c_int), value :: no_num
    type(c_ptr), value    :: neighbor_til_id
    type(c_ptr), value    :: shift
  end subroutine
end interface

interface
  subroutine init_chunk_pool( chunk_pool_size, chunk_size, n_tils ) bind(C)
    use iso_c_binding
    integer(c_int), value :: chunk_pool_size, chunk_size, n_tils
  end subroutine
end interface

interface
  subroutine init_mem_transfer_manager( com_buf_size ) bind(C)
    use iso_c_binding
    integer(c_int), value :: com_buf_size
  end subroutine
end interface

interface
  subroutine flush_mem_transfer_manager() bind(C)
    use iso_c_binding
  end subroutine
end interface

interface
  subroutine init_flds_group() bind(C)
    use iso_c_binding
  end subroutine
end interface

interface
  subroutine init_flds_tile( til_id, jay_nx, jay_gc_num, jay, emf_nx, emf_gc_num, e, &
                             b, if_last_til ) bind(C)
    use iso_c_binding
    integer(c_int), value :: til_id
    type(c_ptr), value    :: jay_nx, jay_gc_num, jay, emf_nx, emf_gc_num, e, b
    logical(c_bool), value:: if_last_til
  end subroutine
end interface

interface
  subroutine init_spec_group( spec_id, rqm, free_stream, dx, interpolation, grid_center, &
                              dt, ihole_size, ihole_size_frac, nchunks_pad, nchunks_pad_frac, &
                              chunk_growth_factor, nblocks_requested, tiles_per_node ) bind(C)
    use iso_c_binding
    use m_parameters
    integer(c_int), value  :: spec_id
    real(c_k_part), value  :: rqm
    logical(c_bool), value :: free_stream
    type(c_ptr), value     :: dx
    integer(c_int), value  :: interpolation
    logical(c_bool), value :: grid_center
    real(c_double), value  :: dt
    integer(c_int), value  :: ihole_size
    real(c_double), value  :: ihole_size_frac
    integer(c_int), value  :: nchunks_pad
    real(c_double), value  :: nchunks_pad_frac
    real(c_double), value  :: chunk_growth_factor
    integer(c_int), value  :: nblocks_requested
    integer(c_int), value  :: tiles_per_node
  end subroutine
end interface

interface
  subroutine init_spec_tile( til_id, spec_id, num_par_max, num_par, energy, my_nx_p, x, &
                             p, q, ix, nchunks_max_init, if_last_til ) bind(C)
    use iso_c_binding
    use m_parameters
    integer(c_int), value :: til_id, spec_id, num_par_max, num_par
    type(c_ptr), value    :: energy, my_nx_p, x, p, q, ix
    integer(c_int), value :: nchunks_max_init
    logical(c_bool), value:: if_last_til
  end subroutine
end interface

interface
  subroutine copy_particles_to_device() bind(C)
    use iso_c_binding
  end subroutine
end interface

interface
  subroutine pre_advance_deposit_cuda() bind(C)
    use iso_c_binding
  end subroutine
end interface

interface
  subroutine advance_deposit_cuda( report_energy ) bind(C)
    use iso_c_binding
    type(c_ptr), value :: report_energy
  end subroutine
end interface

interface
  subroutine post_advance_deposit_cuda() bind(C)
    use iso_c_binding
  end subroutine
end interface

interface
  subroutine move_window_part_cuda( nmove, num_par ) bind(C)
    use iso_c_binding
    type(c_ptr), value :: nmove, num_par
  end subroutine
end interface

interface
  subroutine update_boundary_part_cuda_1( num_par ) bind(C)
    use iso_c_binding
    type(c_ptr), value :: num_par
  end subroutine
end interface

interface
  subroutine copy_particles_ext_from_device( num_par) bind(C)
    use iso_c_binding
    type(c_ptr), value :: num_par
  end subroutine
end interface

interface
  subroutine update_boundary_part_cuda_2( num_par ) bind(C)
    use iso_c_binding
    type(c_ptr), value :: num_par
  end subroutine
end interface

interface
  subroutine copy_report_from_device( copy_spec, report_energy ) bind(C)
    use iso_c_binding
    type(c_ptr), value :: copy_spec, report_energy
  end subroutine
end interface

interface
  subroutine init_dlb_cuda_( n_tils_new, tils_stnry, tils_stnry_old, n_stnry, &
       tils_send_old, n_send, ihole_size, ihole_size_frac, nchunks_pad, nchunks_pad_frac ) bind(C)
    use iso_c_binding
    integer(c_int), value :: n_tils_new
    type(c_ptr), value    :: tils_stnry
    type(c_ptr), value    :: tils_stnry_old
    integer(c_int), value :: n_stnry
    type(c_ptr), value    :: tils_send_old
    integer(c_int), value :: n_send
    integer(c_int), value :: ihole_size
    real(c_double), value  :: ihole_size_frac
    integer(c_int), value :: nchunks_pad
    real(c_double), value  :: nchunks_pad_frac
  end subroutine
end interface

interface
  subroutine init_tile_last_dlb() bind(C)
    use iso_c_binding
  end subroutine
end interface

interface
  subroutine copy_particles_to_device_dlb( n_tils_recv, tils_recv ) bind(C)
    use iso_c_binding
    integer(c_int), value :: n_tils_recv
    type(c_ptr), value    :: tils_recv
  end subroutine
end interface

interface
  subroutine copy_particles_from_device_dlb( n_tils_send, tils_send ) bind(C)
    use iso_c_binding
    integer(c_int), value :: n_tils_send
    type(c_ptr), value    :: tils_send
  end subroutine
end interface

interface
  subroutine init_node_conf_dlb( n_tils, my_aid, no_num, neighbor_til_id, shift ) bind(C)
    use iso_c_binding
    integer(c_int), value :: n_tils
    integer(c_int), value :: my_aid
    integer(c_int), value :: no_num
    type(c_ptr), value    :: neighbor_til_id
    type(c_ptr), value    :: shift
  end subroutine
end interface

interface
  subroutine cleanup_cuda() bind(C)
    use iso_c_binding
  end subroutine
end interface

public :: t_sim_group_cuda

contains

!-----------------------------------------------------------------------------------------
! Initialize simulation objects
!-----------------------------------------------------------------------------------------
subroutine allocate_objs_group_cuda( sim )

  implicit none

  class( t_sim_group_cuda ), intent(inout) :: sim

  if ( .not. associated( sim%no_co ) ) then
    allocate( t_node_conf_cuda :: sim%no_co )
  endif

  call sim % t_simulation_group % allocate_objs()

end subroutine allocate_objs_group_cuda

!-----------------------------------------------------------------------------------------
! Initialize simulation objects
! Necessary to override this routine so we can allocate tiles as t_sim_tiles_cuda
!-----------------------------------------------------------------------------------------
subroutine allocate_tiles_cuda( sim )

  implicit none

  class( t_sim_group_cuda ), intent(inout) :: sim

  integer :: i

  ! Allocate tile containter
  if (.not. associated(sim%tiles)) then
    allocate( sim%tiles( sim%til_id(1):sim%til_id(2) ) )

    ! Allocate tiles
    do i = sim%til_id(1),sim%til_id(2)
      allocate( t_sim_tiles_cuda :: sim%tiles(i)%t )
    enddo
  endif

  call sim % t_simulation_group % allocate_tiles()

end subroutine allocate_tiles_cuda
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Initialize new tiles array for dynamic load balance
! Necessary to override this routine so we can allocate tiles as t_sim_tiles_cuda
!-----------------------------------------------------------------------------------------
subroutine allocate_tiles_dlb_cuda( sim, patt, tiles_temp )

  implicit none

  class( t_sim_group_cuda ), intent(inout)        :: sim
  type( t_msg_patt ), intent(inout)               :: patt
  class( t_til_container ), dimension(:), pointer :: tiles_temp

  integer :: i, j, til

  ! Allocate temp tile container for new tile array
  allocate( tiles_temp( sim%til_id(1):sim%til_id(2) ) )

  ! Allocate, allocate objects, read input for incoming tiles
  if ( patt%n_recv>0 ) then

    do i = 1, patt%n_recv
      do j = 1, size( patt%recv(i)%tils )
        til = patt%recv(i)%tils(j)
        allocate( t_sim_tiles_cuda :: tiles_temp(til)%t )
      enddo
    enddo

    !$omp parallel do schedule(dynamic)
    do i = 1, patt%n_recv
      !$omp critical
      do j = 1, size( patt%recv(i)%tils )
        til = patt%recv(i)%tils(j)
        call tiles_temp(til)%t % allocate_objs()
      enddo
      !$omp end critical
    enddo
    !$omp end parallel do

    ! Read input for incoming tiles
    ! set input file position to same as it was when reading for tiles initially
    do i = 1, patt%n_recv
      do j = 1, size( patt%recv(i)%tils )
        til = patt%recv(i)%tils(j)
        sim%input_file%pos = sim%pos
        call tiles_temp(til)%t % read_copy_sim( sim, sim%input_file, til )
      enddo
    enddo

  endif

end subroutine allocate_tiles_dlb_cuda
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Initialize simulation objects
!-----------------------------------------------------------------------------------------
subroutine init_sim_group_cuda( sim, grid_gc_min )

  use m_tile_conf, only: t_tile_conf

  implicit none

  class ( t_sim_group_cuda ), intent(inout) :: sim
  integer, dimension( 2, p_x_dim ), intent(in) :: grid_gc_min

  integer :: n_tils, sp_id, til, til_id, i, tiles_per_node
  class(t_species), pointer :: spec
  class(t_vdf), pointer :: jay, e, b
  integer, dimension(:,:), pointer :: neighbor_til_id
  integer, dimension(:,:,:), pointer :: shift

  type( t_grid_arr ), dimension(:), pointer :: grid_arr
  type( t_tile_conf ), pointer :: ti_co

  integer :: chunk_pool_size, chunk_size, com_buf_size, ihole_size, nchunks_pad
  integer :: nchunks_max_init, nblocks_requested
  real(p_double) :: ihole_size_frac, nchunks_pad_frac, chunk_growth_factor

  ! call super class
  call init_sim_group(sim, grid_gc_min)

  SCR_ROOT('')
  SCR_ROOT('Initializing CUDA simulation objects:')

  n_tils = sim%til_id(3)

  select type( no_co => sim%no_co )
  class is( t_node_conf_cuda )
    ti_co              => no_co%ti_co
    chunk_pool_size     = no_co%chunk_pool_size
    com_buf_size        = no_co%com_buf_size
    chunk_size          = no_co%chunk_size
    ihole_size          = no_co%ihole_size
    ihole_size_frac     = no_co%ihole_size_frac
    nchunks_pad         = no_co%nchunks_pad
    nchunks_pad_frac    = no_co%nchunks_pad_frac
    nchunks_max_init    = no_co%nchunks_max_init
    chunk_growth_factor = no_co%chunk_growth_factor
    nblocks_requested   = no_co%nblocks_requested
  end select

  tiles_per_node = ti_co%g_til_num / sim%no_co%no_num ! no_co%numnah

  ! initialize cuda objects, groups first
  call alloc( neighbor_til_id, (/ 3**p_x_dim, n_tils /) )
  call alloc( shift, (/ p_x_dim, 3**p_x_dim, n_tils /) )

  allocate( grid_arr( n_tils ) )
  do i = 1, n_tils
    grid_arr(i)%g => sim%tiles(i+sim%til_id(1)-1)%t%grid
  enddo

  select type( no_co => sim%no_co )
  class is( t_node_conf_tiles )
    call populate_neighbors( no_co, grid_arr, sim%til_id, neighbor_til_id, shift )
  end select
  deallocate(grid_arr)

  call init_node_conf( n_tils, sim%no_co%my_aid_, sim%no_co%no_num, &
                       c_loc(neighbor_til_id), c_loc(shift) )

  call freemem( neighbor_til_id )
  call freemem( shift )

  SCR_ROOT('  Setting up chunk pool')
  call init_chunk_pool( chunk_pool_size, chunk_size, n_tils )

  SCR_ROOT('  Setting up memory transfer manager')
  call init_mem_transfer_manager( com_buf_size )

  SCR_ROOT('  Setting up emf and particles')
  ! initialize cuda fields
  call init_flds_group()

  sp_id = 1
  spec => sim%tiles(sim%til_id(1))%t%part%species
  do
    if (.not. associated(spec)) exit
    call init_spec_group( sp_id, spec%rqm, logical(spec%free_stream,kind=c_bool), &
                          c_loc(spec%dx), spec%interpolation, &
                          logical(spec%grid_center,kind=c_bool), dt(sim%tstep), ihole_size, &
                          ihole_size_frac, nchunks_pad, nchunks_pad_frac, chunk_growth_factor, &
                          nblocks_requested, tiles_per_node )
    spec => spec % next
    sp_id = sp_id + 1
  enddo

  ! initialize cuda objects, tile by tile
  do til = 1, n_tils
    til_id = til + sim%til_id(1) - 1
    jay => sim%tiles(til_id)%t%jay%pf(1)
    e   => sim%tiles(til_id)%t%emf%e_part
    b   => sim%tiles(til_id)%t%emf%b_part
    call init_flds_tile( til, c_loc(jay%nx_), c_loc(jay%gc_num_), c_loc(jay%buffer), &
                         c_loc(e%nx_), c_loc(e%gc_num_), c_loc(e%buffer), c_loc(b%buffer), &
                         logical(.false.,kind=c_bool) )
  enddo

  call flush_mem_transfer_manager();

  do til = 1, n_tils
    til_id = til + sim%til_id(1) - 1

    sp_id = 1
    spec => sim%tiles(til_id)%t%part%species
    do
      if (.not. associated(spec)) exit
      call init_spec_tile( til, sp_id, spec%num_par_max, spec%num_par,c_loc(spec%energy),&
                           c_loc(spec%my_nx_p), c_loc(spec%x), c_loc(spec%p), &
                           c_loc(spec%q), c_loc(spec%ix), nchunks_max_init, logical(.false.,kind=c_bool) )
      spec => spec % next
      sp_id = sp_id + 1
    enddo
  enddo

  SCR_ROOT('  Copying particles to device')
  call copy_particles_to_device()


  SCR_ROOT('')
  SCR_ROOT('CUDA initialization complete!')

end subroutine init_sim_group_cuda
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Fill in tile ids for neighbors of each tile in all directions (including corners).
! These ids are referenced to the list of tiles on this MPI rank, with 0 being the first.
!-----------------------------------------------------------------------------------------
subroutine populate_neighbors( no_co, grid_arr, til_id, neighbor_til_id, shift )

  implicit none

  class( t_node_conf_tiles ), intent(in) :: no_co
  type(t_grid_arr), dimension(:), pointer :: grid_arr
  integer, dimension(3), intent(in) :: til_id
  integer, dimension(:,:), pointer :: neighbor_til_id
  integer, dimension(:,:,:), pointer :: shift

  integer :: i, j, k, x, b, n
  integer, dimension(p_max_dim) :: ind

  neighbor_til_id = 0
  shift = 0

  do k = 1, til_id(3)

    ! Loop over neighbors of tile k
    do j = 1, 3**p_x_dim

      ! Get coordinates of tile k in the global grid of tiles
      ind(1:no_co%x_dim) = no_co%ti_co%h(:,k+til_id(1)-1)

      ! Get coordinates of neighbor j
      x = j-1
      b = 3 ! converting j to base 3
      i = 1 ! start with index 1
      do while (x>0)
        n=int(x-b*int(x/b))
        ! Do nothing if n==0
        if (n==1) then
          ind(i) = ind(i) - 1
        else if (n==2) then
          ind(i) = ind(i) + 1
        endif
        x=int(x/b)
        i = i + 1
      enddo

      ! Now, set neighbor_til_id and shift for neighbor j of tile k
      ! neighbor_til_id to -1 for physical boundary or external neighbor, else just til_id
      ! set shift same as in CPU update boundary for local neighbors, 0 for external neighbors
      ! or phys boundaries (since these will be handled on the CPCPU)

      do i = 1, p_x_dim

        ! Check if the lower bound has gone outside domain
        if (ind(i)<1) then
          if (no_co%ifpr_(i)) then
            ind(i) = no_co%ti_co%g_til(i)
            shift(i,j,k) = +grid_arr(1)%g%g_nx(i)
          else
            neighbor_til_id(j,k) = -1
          endif
        endif

        ! Check if the upper bound has gone outside domain
        if (ind(i)>no_co%ti_co%g_til(i)) then
          if (no_co%ifpr_(i)) then
            ind(i) = 1
            shift(i,j,k) = -grid_arr(1)%g%g_nx(i)
          else
            neighbor_til_id(j,k) = -1
          endif
        endif

      enddo

      if (neighbor_til_id(j,k) /= -1) then

        ! Fill in til_id, zero indexed to this rank's tiles
        select case( no_co%x_dim )
        case(1)
          neighbor_til_id(j,k) = no_co%ti_co%hi1( ind(1) ) - til_id(1)
        case(2)
          neighbor_til_id(j,k) = no_co%ti_co%hi2( ind(1), ind(2) ) - til_id(1)
        case(3)
          neighbor_til_id(j,k) = no_co%ti_co%hi3( ind(1), ind(2), ind(3) ) - til_id(1)
        end select

        if ( neighbor_til_id(j,k) < 0 .or. neighbor_til_id(j,k) > til_id(3)-1 ) then
          ! Set til_id to -1 if neighbor j belongs to different rank from tile k
          neighbor_til_id(j,k) = -1
        else
          ! Only set shift if neighbor j belongs to same rank as tile k
          shift(1:p_x_dim,j,k) = shift(1:p_x_dim,j,k) &
                                 + grid_arr(k)%g%my_nx(p_lower, 1:p_x_dim) &
                                 - grid_arr(neighbor_til_id(j,k)+1)%g%my_nx(p_lower, 1:p_x_dim)
        endif

      endif

      ! shift should be zero if the neighbor is external
      if (neighbor_til_id(j,k) == -1) then
        shift(1:p_x_dim,j,k) = 0
      endif

    enddo

  enddo

end subroutine populate_neighbors
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine advance_deposit_group_cuda( sim )

  implicit none

  class( t_sim_group_cuda ), intent(inout) :: sim

  integer :: til, til_id, sp_id
  integer, target :: num_par(sim%til_id(3),sim%part%num_species)
  class(t_species), pointer :: spec
  logical(kind=c_bool), dimension(sim%part%num_species), target :: report_energy

  call begin_event( push_outside_ev )

  do til = 1, sim%til_id(3)
    til_id = til + sim%til_id(1) - 1
    ! call sim%tiles(til_id)%t%jay%pf(1)%zero()

    sp_id = 1
    spec => sim%tiles(til_id)%t%part%species
    do
      if (.not. associated(spec)) exit
      num_par(til,sp_id) = spec%num_par
      spec => spec % next
      sp_id = sp_id + 1
    enddo
  enddo

  report_energy(:) = sim % report_energy_dry_part()

  call pre_advance_deposit_cuda()

  call begin_event( pushev )
  call advance_deposit_cuda( c_loc(report_energy) )
  call end_event( pushev )

  call post_advance_deposit_cuda()

  do til = 1, sim%til_id(3)
    til_id = til + sim%til_id(1) - 1
    sim%tiles(til_id)%t%part%n_current = sim%tiles(til_id)%t%part%n_current + 1
  enddo

  call end_event( push_outside_ev )

end subroutine advance_deposit_group_cuda
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Note, no support yet for cathode or neutral species
!-----------------------------------------------------------------------------------------
subroutine move_window_part_group_cuda( sim )

  implicit none

  class( t_sim_group_cuda ), intent(inout) :: sim

  integer :: til, til_id, sp_id, nmove_tot
  integer, target :: nmove(p_max_dim), num_par(sim%til_id(3),sim%part%num_species)
  class(t_particles), pointer :: part
  class(t_species), pointer :: spec

  ! Get nmove for each dimension, return if not moving at all
  nmove = 0
  nmove_tot = 0
  do til = 1, p_x_dim
    nmove(til) = nx_move( sim%g_space, til )
    nmove_tot = nmove_tot + abs(nmove(til))
  enddo
  if (nmove_tot==0) return

  ! Start the move_window process
  do til = sim%til_id(1), sim%til_id(2)

    ! Set num_par equal to zero for each species so no particles are edited
    part => sim%tiles(til)%t%part
    spec => part%species
    sp_id = 1
    do
      if (.not. associated(spec)) exit
      spec%num_par = 0
      spec => spec%next
      sp_id = sp_id + 1
    enddo

    ! Perform CPU moving window (mainly just injection)
    ! Note: after this, spec%num_par is equal to the number of particles injected
    call part % move_window( sim%g_space, sim%tiles(til)%t%grid, sim%tiles(til)%t%jay, &
                             sim%tiles(til)%t%no_co, sim%tiles(til)%t%bnd )

    ! store number of injected particles (to be passed to C/CUDA)
    spec => part%species
    sp_id = 1
    do
      if (.not. associated(spec)) exit
      num_par(til-sim%til_id(1)+1,sp_id) = spec%num_par
      spec => spec%next
      sp_id = sp_id + 1
    enddo

  enddo

  ! Edit GPU particle indices, then pass injected particles from CPU to GPU
  ! Note: after this is called, num_par is equal to the total number of particles on each tile
  call move_window_part_cuda( c_loc(nmove), c_loc(num_par) )

  ! update num_par for each tile for each species to reflect actual number of particles
  do til = sim%til_id(1), sim%til_id(2)
    part => sim%tiles(til)%t%part
    spec => part%species
    sp_id = 1
    do
      if (.not. associated(spec)) exit
      spec%num_par = num_par(til-sim%til_id(1)+1, sp_id)
      spec => spec%next
      sp_id = sp_id + 1
    enddo
  enddo

  ! Copy particles from GPU to CPU for reporting if necessary this time step
  call prep_report_move_window( sim )

end subroutine move_window_part_group_cuda
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! updates the particle data at the boundaries
!-----------------------------------------------------------------------------------------
subroutine update_boundary_part_group_cuda( this )

  use m_particles, only : validate

  implicit none

  class(t_sim_group_cuda), intent(inout) :: this

  integer, target :: num_par(this%til_id(3),this%part%num_species)
  class(t_species), pointer :: spec
  integer :: til, til_id, sp_id

  call begin_event(partboundev)

  ! Sort particles and gather information about holes in the particle buffer in CUDA
  call update_boundary_part_cuda_1( c_loc(num_par) )

  ! Set num_par to reflect only relevant particles. And make sure CPU buffers are big enough
  do til = 1, this%til_id(3)
    til_id = til + this%til_id(1) - 1
    sp_id = 1
    spec => this%tiles(til_id)%t%part%species
    do
      if (.not. associated(spec)) exit

      spec%num_par = num_par(til,sp_id)

      if ( spec%num_par > spec%num_par_max ) then
        call spec%grow_buffer( spec%num_par + p_spec_buf_block_t )
      endif

      spec => spec % next
      sp_id = sp_id + 1
    enddo
  enddo

  ! Now that there's space, copy relevant particles from device to host
  call copy_particles_ext_from_device( c_loc(num_par) )


#ifdef __DEBUG__
  do til = this%til_id(1), this%til_id(2)
    call validate( this%tiles(til)%t%part, 'after update_boundary_part_cuda_1', over=.true. )
  enddo
  call this%no_co%barrier()
#endif


  ! call CPU update boundary routines (mainly this passes particles via MPI to other nodes)
  call this % update_boundary_part_cpu()

  ! update num_par to reflect number of particles which need to be processed by GPU
  do til = 1, this%til_id(3)
    til_id = til + this%til_id(1) - 1
    sp_id = 1
    spec => this%tiles(til_id)%t%part%species
    do
      if (.not. associated(spec)) exit
      num_par(til,sp_id) = spec%num_par
      spec => spec % next
      sp_id = sp_id + 1
    enddo
  enddo

  ! Fill holes in each tile's particle buffer in CUDA
  call update_boundary_part_cuda_2( c_loc(num_par) )

  ! update num_par for each tile for each species to reflect actual number of particles
  do til = 1, this%til_id(3)
    til_id = til + this%til_id(1) - 1
    sp_id = 1
    spec => this%tiles(til_id)%t%part%species
    do
      if (.not. associated(spec)) exit
      spec%num_par = num_par(til,sp_id)
      spec => spec % next
      sp_id = sp_id + 1
    enddo
  enddo

  call prep_report_update_bound( this )

  call end_event(partboundev)

end subroutine update_boundary_part_group_cuda
!-----------------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! updates the particle data at the boundaries on the CPU side
! Notes:
!  * Naively one would think local comm routines are not needed here, since they should be 
!    handled in CUDA. But actually, they are still needed for the case where a particle
!    crosses both an external boundary and at least one local boundary. 
!  * The only difference btwn this routine and update_boundary_part_group() is that here
!    we don't call begin_event(partboundev), since this is already done inside 
!    update_boundary_part_group_cuda. A huge chunk of duplicated code is a big price to pay
!    for that but, short of implementing some sort of wrapper, worth it
!-------------------------------------------------------------------------------
subroutine update_boundary_part_cpu_group_cuda( sim )

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

  class(t_sim_group_cuda), intent(inout) :: sim

  type(t_spec_arr), dimension(sim%til_id(1):sim%til_id(2)) :: spec
  integer, dimension(sim%til_id(1):sim%til_id(2)) :: comms_int
  type(t_vdf_arr), dimension(:), pointer :: multi_ion_arr
  integer :: t, dim, i

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

#ifdef __HAS_IONIZATION__
  if (sim%tiles(sim%til_id(1))%t%part%num_neutral>0) then
    ERROR("Neutrals not yet implemented for CUDA (although it wouldn't be too much work...maybe you, user, should be the one to do it)")
    call abort_program( p_err_notimplemented )
  endif
#endif

end subroutine update_boundary_part_cpu_group_cuda
!-------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
! Check if particles/energy are needed for reporting on next time step
! If so, copy them from device to host
!-----------------------------------------------------------------------------------------
subroutine prep_report_update_bound( this )

  implicit none

  class(t_sim_group_cuda), intent(inout) :: this

  logical(kind=c_bool), dimension(this%part%num_species), target :: copy_spec, report_energy
  class(t_emf), pointer :: emf
  logical :: if_move_next
  class(t_species), pointer :: spec

  integer :: til, til_id, sp_id

  ! Should particles be communicated to CPU for reporting next time step?
  call advance( this%tstep, 1 )
  emf => this%tiles(this%til_id(1))%t%emf
  if_move_next = move_boundary_dry( this%g_space, emf%dx(), dt(this%tstep) )
  if (if_move_next) then
    copy_spec = .false.
  else
    copy_spec(:) = this % report_dry_part()
  endif
  report_energy(:) = this % report_energy_dry_part()
  call advance( this%tstep, -1 )

  ! Are we doing a restart at the end of this timestep? If so particles should be 
  ! comm'd to the CPU
  if ( if_restart_write( this%restart, n(this%tstep), ndump(this%tstep), &
                                  comm(this%no_co), this%no_co%ngp_id() ) ) then
    copy_spec(:) = .true.
  end if

  ! Of species that need to be communicated to CPU for reporting, check if any buffers
  ! need to be resized
  do til = 1, this%til_id(3)
    til_id = til + this%til_id(1) - 1
    sp_id = 1
    spec => this%tiles(til_id)%t%part%species
    do
      if (.not. associated(spec)) exit
      if ( copy_spec(sp_id) ) then
        if ( spec%num_par > spec%num_par_max ) then
          call spec%grow_buffer( spec%num_par + p_spec_buf_block_t )
        endif
      endif
      spec => spec % next
      sp_id = sp_id + 1
    enddo
  enddo

  ! Copy energy and particles to CPU for reporting
  call copy_report_from_device( c_loc(copy_spec), c_loc(report_energy) )

end subroutine prep_report_update_bound
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
! Check if particles are needed for reporting on next time step, if so, copy them from
! device to host
!-----------------------------------------------------------------------------------------
subroutine prep_report_move_window( this )

  implicit none

  class(t_sim_group_cuda), intent(inout) :: this

  logical(kind=c_bool), dimension(this%part%num_species), target :: copy_spec, report_energy
  integer, dimension(this%part%num_species), target :: num_par_max
  class(t_emf), pointer :: emf
  logical :: if_move_next
  class(t_species), pointer :: spec

  integer :: til, til_id, sp_id

  ! Should particles be communicated to CPU for reporting next time step?
  copy_spec(:) = this % report_dry_part()

  ! energy never needs to be communicated to CPU at this stage of the loop
  report_energy(:) = .false.

  ! Of species that need to be communicated to CPU for reporting, check if any buffers
  ! need to be resized
  do til = 1, this%til_id(3)
    til_id = til + this%til_id(1) - 1
    sp_id = 1
    spec => this%tiles(til_id)%t%part%species
    do
      if (.not. associated(spec)) exit
      if ( copy_spec(sp_id) ) then
        if ( spec%num_par > spec%num_par_max ) then
          call spec%grow_buffer( spec%num_par + p_spec_buf_block_t )
        endif
      endif
      spec => spec % next
      sp_id = sp_id + 1
    enddo
  enddo

  ! Copy energy and particles to CPU for reporting
  call copy_report_from_device( c_loc(copy_spec), c_loc(report_energy) )

end subroutine prep_report_move_window
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
! dry run of move_boundary, to see if particles will need to be copied back after
! move_window instead
!-----------------------------------------------------------------------------------------
function move_boundary_dry( space, dx, dt )

  implicit none

  type(t_space), intent(in) :: space
  real(p_double), dimension(:), intent(in) :: dx
  real(p_double),               intent(in) :: dt
  logical :: move_boundary_dry

  integer :: i_dim
  real(p_double), dimension(2,p_max_dim) :: x_bnd_old, xmoved

  move_boundary_dry = .false.

  x_bnd_old = space%x_bnd

  do i_dim=1, space%x_dim
    if ( space%if_move(i_dim) ) then

      if ( space%move_vel > 0.0_p_double ) then
        ! move boundary at specified velocity
        xmoved(:,i_dim) = space%xmoved(:,i_dim) + space%move_vel*dt
      else
        ! move boundary at c
        xmoved(:,i_dim) = space%xmoved(:,i_dim) + dt
      endif

      if ( xmoved(1,i_dim) >= dx(i_dim) ) then
        move_boundary_dry = .true.
      endif
    endif
  enddo

end function move_boundary_dry
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! dry run of report particles, check to see if we need to send species to CPU
!-----------------------------------------------------------------------------------------
function report_dry_part_group_cuda( sim )

  implicit none

  class(t_sim_group_cuda), intent(inout) :: sim

  logical, dimension(sim%part%num_species) :: report_dry_part_group_cuda

  integer :: j
  type( t_vdf_report ), pointer :: rep
  class( t_species ), pointer :: spec
  logical :: needs_update_charge

  report_dry_part_group_cuda(:) = .false.

  ! Species diagnostics (phasespaces, energy, temperature, heat flux, RAW, tracks, etc.)
  spec => sim%tiles(sim%til_id(1))%t%part%species
  j = 1
  do
    if (.not. associated(spec)) exit
    if (.not. report_dry_part_group_cuda(j)) then
      report_dry_part_group_cuda(j) = report_dry_spec_group_cuda( spec, sim )
    endif
    spec => spec%next
    j = j + 1
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
    report_dry_part_group_cuda(:) = .true.
  endif

  contains

  function report_dry_spec_group_cuda( spec, sim )

    implicit none

    class( t_species ), intent(in) :: spec
    class( t_simulation_group ), intent(inout) :: sim

    logical :: report_dry_spec_group_cuda

    type( t_vdf_report ), pointer :: rep

    report_dry_spec_group_cuda = .false.

    ! Density diagnostics
    rep => spec%diag%reports
    do
      if (.not. associated(rep)) exit

      if ( if_report( rep, sim%tstep ) ) then

        report_dry_spec_group_cuda = .true.

      endif

      rep => rep % next

    enddo

    ! Momentum distribution diagnostics are not currently implemented

    ! Cell average diagnostics
    rep => spec%diag%rep_cell_avg
    do
      if ( .not. associated( rep ) ) exit

      if ( if_report( rep, sim%tstep ) ) then

        report_dry_spec_group_cuda = .true.

      endif

      rep => rep%next
    enddo

    ! Phasespaces
    if (test_if_report( sim%tstep, spec%diag%phasespaces%ndump_fac )) then

      report_dry_spec_group_cuda = .true.

    endif

    ! Raw particle data diagnostics
    if (test_if_report( sim%tstep, spec%diag%ndump_fac_raw )) then

      report_dry_spec_group_cuda = .true.

    endif

    ! Temperature and heat flux diagnostics not currently implemented

    ! Particle tracking diagnostics not currently implemented

  end function report_dry_spec_group_cuda

end function report_dry_part_group_cuda
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! dry run of report particles, check to see if we need to send species to CPU
!-----------------------------------------------------------------------------------------
function report_energy_dry_part_group_cuda( sim )

  implicit none

  class(t_sim_group_cuda), intent(inout) :: sim

  logical, dimension(sim%part%num_species) :: report_energy_dry_part_group_cuda

  integer :: j
  class( t_species ), pointer :: spec

  report_energy_dry_part_group_cuda(:) = .false.

  ! Species diagnostics (phasespaces, energy, temperature, heat flux, RAW, tracks, etc.)
  spec => sim%tiles(sim%til_id(1))%t%part%species
  j = 1
  do
    if (.not. associated(spec)) exit
    if (.not. report_energy_dry_part_group_cuda(j)) then
      report_energy_dry_part_group_cuda(j) = report_energy_dry_spec_group_cuda( spec, sim )
    endif
    spec => spec%next
    j = j + 1
  enddo

  ! Particle energy diagnostics
  if ( test_if_report( sim%tstep, sim%part%ndump_fac_ene ) ) then
    report_energy_dry_part_group_cuda(:) = .true.
  endif

  contains

  function report_energy_dry_spec_group_cuda( spec, sim )

    implicit none

    class( t_species ), intent(in) :: spec
    class( t_simulation_group ), intent(inout) :: sim

    logical :: report_energy_dry_spec_group_cuda

    report_energy_dry_spec_group_cuda = .false.

    ! Energy diagnostics
    if (test_if_report( sim%tstep, spec%diag%ndump_fac_ene )) then

      report_energy_dry_spec_group_cuda = .true.

    endif

  end function report_energy_dry_spec_group_cuda

end function report_energy_dry_part_group_cuda
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine cleanup_sim_group_cuda( sim )

  implicit none

  class( t_sim_group_cuda ), intent(inout) :: sim

  ! call superclass
  call sim % t_simulation_group % cleanup()

  call cleanup_cuda()

end subroutine cleanup_sim_group_cuda
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Initialize new tiles array for dynamic load balance
!-----------------------------------------------------------------------------------------
subroutine init_dlb_cuda( sim, patt, tiles_temp, dif, diff, difff )

  use m_node_conf
  use m_species_cuda

  implicit none

  class( t_sim_group_cuda ), intent(inout)        :: sim
  type( t_msg_patt ), intent(inout)               :: patt
  class( t_til_container ), dimension(:), pointer :: tiles_temp
  integer, dimension(:), optional                 :: dif
  integer, dimension(:,:), optional               :: diff
  integer, dimension(:,:,:), optional             :: difff

  class( t_species ), pointer :: spec
  class(t_node_conf), pointer :: no_co

  integer :: n_tils_send, n_tils_stnry, n_tils_new
  integer, dimension(:), pointer :: tils_send_old, tils_stnry, tils_stnry_old
  integer, dimension(3) :: sim_til_id_old ! to store the pre-dlb value of sim%til_id
  integer :: ihole_size, nchunks_pad
  real(p_double) :: ihole_size_frac, nchunks_pad_frac
  
  integer :: i, j
  integer :: til_id


  ! copy tiles to be sent via MPI from device to host
  ! have to do this here, not in send_recv(), because we delete these tiles in the next step
  no_co => sim%no_co
  select type(no_co); class is(t_node_conf_tiles)
    sim_til_id_old = no_co%ti_co % g_til_aid_min_max_old(:,no_co%my_aid())
    call get_tils_send( no_co, n_tils_send, tils_send_old, patt, sim_til_id_old(1) )
  end select

  ! Make room for incoming data: resize particle buffers for this tile if necessary
  if ( n_tils_send > 0 ) then
    do i = 1, patt%n_send
      do j = 1, size( patt%send(i)%tils )
        til_id = patt%send(i)%tils(j)
        spec => sim%tiles(til_id)%t%part%species
        do
          if (.not. associated(spec)) exit
          if ( spec%num_par > spec%num_par_max ) then
            select type( spec ); class is( t_species_cuda )
              spec%til_id_min = sim_til_id_old(1)
            end select
            call spec%grow_buffer( spec%num_par + p_spec_buf_block_t )
          endif
          spec => spec % next
        enddo
      enddo
    enddo

    ! only need to copy particles, fields are already on CPU
    call copy_particles_from_device_dlb( n_tils_send, c_loc(tils_send_old) )

  endif


  ! call superclass init function
  call init_dlb_group( sim, patt, tiles_temp, dif, diff, difff )

  ! Update species member data for all tiles about the new partition
  do i = sim%til_id(1), sim%til_id(2)
    spec => tiles_temp(i)%t%part%species
    no_co => tiles_temp(i)%t%no_co
    select type(no_co); class is(t_node_conf_tiles)
    do
      if (.not. associated(spec)) exit
        select type( spec ); class is( t_species_cuda )
          spec % til_id     = no_co % ti_co % g_til_aid 
          spec % til_id_min = no_co % ti_co % g_til_aid_min_max( 1, no_co%local_rank()+1 )
        end select
      spec => spec % next
    enddo
    end select
  enddo


  ! Reinitialize simulation objects in C/CUDA
  n_tils_new = sim%til_id(3)

  ! List of tiles staying on this node, using the post- and pre-dlb tile id, respectively.
  ! At most, all tiles were stationary.
  call alloc( tils_stnry, (/ n_tils_new /) )
  call alloc( tils_stnry_old, (/ n_tils_new /) )

  ! obtain a list of all stationary tiles
  no_co => sim%no_co
  select type(no_co); class is(t_node_conf_tiles)
    call get_tils_stnry( no_co, sim%til_id, sim_til_id_old, tils_stnry, tils_stnry_old, &
                         n_tils_stnry, dif, diff, difff )
  end select

  ! unpack nchunks_pad and nchunks_pad_frac
  select type( no_co => sim%no_co )
  class is( t_node_conf_cuda )
    ihole_size          = no_co%ihole_size
    ihole_size_frac     = no_co%ihole_size_frac
    nchunks_pad         = no_co%nchunks_pad
    nchunks_pad_frac    = no_co%nchunks_pad_frac
  end select
  
  ! reallocate tiled arrays to prepare for incoming tiles, make new pointers to stationary tiles
  call init_dlb_cuda_( n_tils_new, c_loc(tils_stnry), c_loc(tils_stnry_old), n_tils_stnry, &
                       c_loc(tils_send_old), n_tils_send, ihole_size, ihole_size_frac, & 
                       nchunks_pad, nchunks_pad_frac )

  call freemem( tils_stnry )
  call freemem( tils_stnry_old )
  if (n_tils_send > 0) call freemem( tils_send_old )


end subroutine init_dlb_cuda
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine get_tils_send( no_co, n_tils_send, tils_send_old, patt, til_id_old_1 )

  use m_tile_conf

  implicit none

  class(t_node_conf_tiles), intent(in)            :: no_co
  integer, intent(out)                            :: n_tils_send
  integer, dimension(:), pointer, intent(inout)   :: tils_send_old
  type( t_msg_patt ), intent(inout)               :: patt
  integer, intent(in)                             :: til_id_old_1

  integer, dimension(p_x_dim) :: tgp
  integer :: i, j, k, til_id

  ! aggregate list of all tiles to be sent, regardless of to which node, using the old
  ! (pre-dlb) tile id
  n_tils_send = 0
  do i = 1, patt%n_send
    n_tils_send = n_tils_send + size( patt%send(i)%tils )
  enddo

  if ( n_tils_send > 0 ) call alloc( tils_send_old, (/ n_tils_send /) )

  n_tils_send = 0
  do i = 1, patt%n_send
    do j = 1, size( patt%send(i)%tils )
      n_tils_send = n_tils_send + 1

      ! post-dlb tile index
      til_id = patt%send(i)%tils(j)

      ! we want the pre-dlb tile index so we can grab the tiles from the old array
      ! (only important in viktor topology, otherwise these are the same)
      select case ( no_co%ti_co%topology_type )
        case (p_topology_viktor)
          ! for viktor topology, this makes a difference
          do k = 1, p_x_dim
            tgp(k) = g_tgp( no_co%ti_co, til_id, k )
          enddo
          if ( p_x_dim==2 ) til_id = no_co%ti_co%hi2_old( tgp(1), tgp(2) )
          if ( p_x_dim==3 ) til_id = no_co%ti_co%hi3_old( tgp(1), tgp(2), tgp(3) )
        case default
          ! for all other topologies, pre and post dlb g_til_aid's are the same always
      end select

      tils_send_old( n_tils_send ) = til_id - til_id_old_1
    enddo
  enddo

end subroutine get_tils_send
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
! get a list of all stationary tiles (those that remained on this node over dynamic load balance)
!-----------------------------------------------------------------------------------------
subroutine get_tils_stnry( no_co, til_id, til_id_old, tils_stnry, tils_stnry_old, n_tils_stnry, &
                           dif, diff, difff )

  use m_tile_conf
  use m_node_conf

  implicit none

  class(t_node_conf_tiles), intent(in)            :: no_co
  integer, dimension(:), intent(in)               :: til_id
  integer, dimension(:), intent(in)               :: til_id_old
  integer, dimension(:), intent(inout)            :: tils_stnry
  integer, dimension(:), intent(inout)            :: tils_stnry_old
  integer, intent(inout)                          :: n_tils_stnry
  integer, dimension(:), optional                 :: dif
  integer, dimension(:,:), optional               :: diff
  integer, dimension(:,:,:), optional             :: difff

  integer, dimension(p_x_dim) :: lb, ub
  integer :: til, til_old
  integer :: i, j, k

  ! no_co => sim%no_co
  ! select type(no_co); class is(t_node_conf_tiles)

  n_tils_stnry = 0

  select case ( no_co%ti_co%topology_type )
    case (p_topology_viktor)
      select case( p_x_dim )
      case (2)
        ! dimension 1 lb and ub
        lb(1) = no_co%ti_co%cart_til_mm(1)%ii( 1, my_ngp(no_co,1), 1, 1 )
        ub(1) = no_co%ti_co%cart_til_mm(1)%ii( 2, my_ngp(no_co,1), 1, 1 )
        ! dimension 2 lb and ub
        lb(2) = no_co%ti_co%cart_til_mm(2)%ii( 1, my_ngp(no_co,1), my_ngp(no_co,2), 1 )
        ub(2) = no_co%ti_co%cart_til_mm(2)%ii( 2, my_ngp(no_co,1), my_ngp(no_co,2), 1 )
        ! get old tile data
        do i = lb(1), ub(1)
          do j = lb(2), ub(2)
            if ( diff(i,j)==0 ) then
              ! post load balance g_til_aid, used to index the new tile array
              til = no_co%ti_co%hi2( i, j )
              ! pre load balance g_til_aid, used to index the old tile array
              til_old = no_co%ti_co%hi2_old(i,j)

              n_tils_stnry = n_tils_stnry + 1

              tils_stnry( n_tils_stnry ) = til - til_id(1)
              tils_stnry_old( n_tils_stnry ) = til_old - til_id_old(1)
            endif
          enddo
        enddo
      case (3)
        ! dimension 1 lb and ub
        lb(1) = no_co%ti_co%cart_til_mm(1)%ii( 1, my_ngp(no_co,1), 1, 1 )
        ub(1) = no_co%ti_co%cart_til_mm(1)%ii( 2, my_ngp(no_co,1), 1, 1 )
        ! dimension 2 lb and ub
        lb(2) = no_co%ti_co%cart_til_mm(2)%ii( 1, my_ngp(no_co,1), my_ngp(no_co,2), 1 )
        ub(2) = no_co%ti_co%cart_til_mm(2)%ii( 2, my_ngp(no_co,1), my_ngp(no_co,2), 1 )
        ! dimension 3 lb and ub
        lb(3) = no_co%ti_co%cart_til_mm(3)%ii( 1, my_ngp(no_co,1), my_ngp(no_co,2), &
                                             my_ngp(no_co,3) )
        ub(3) = no_co%ti_co%cart_til_mm(3)%ii( 2, my_ngp(no_co,1), my_ngp(no_co,2), &
                                             my_ngp(no_co,3) )
        ! get old tile data
        do i = lb(1), ub(1)
          do j = lb(2), ub(2)
            do k = lb(3), ub(3)
              if ( difff(i,j,k)==0 ) then
                ! post load balance g_til_aid, used to index the new tile array
                til = no_co%ti_co%hi3( i, j, k )
                ! pre load balance g_til_aid, used to index the old tile array
                til_old = no_co%ti_co%hi3_old(i,j,k)

                n_tils_stnry = n_tils_stnry + 1

                tils_stnry( n_tils_stnry ) = til - til_id(1)
                tils_stnry_old( n_tils_stnry ) = til_old - til_id_old(1)
              endif
            enddo
          enddo
        enddo
      end select
    case default
      ! for all other topologies it's much more straight forward
      do i = til_id(1), til_id(2)
        if (dif(i)==0) then
          til = i
          n_tils_stnry = n_tils_stnry + 1
          tils_stnry( n_tils_stnry ) = til - til_id(1)
          tils_stnry_old( n_tils_stnry ) = til - til_id_old(1)
        endif
      enddo
  end select

end subroutine get_tils_stnry
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! move tiles to proper process in ye olde dynamic load-e balance
! and copy received tiles to device
!-----------------------------------------------------------------------------------------
subroutine send_recv_tils_dlb_cuda( this, patt, tiles )

  use m_simulation_group, only : t_simulation_group, t_til_container, cleanup_tils_out
  use m_grid_define, only : t_msg_patt

  implicit none

  class( t_sim_group_cuda ), intent(inout)        :: this
  type( t_msg_patt ), intent(inout)               :: patt ! who's patt?
  class( t_til_container ), dimension(:), pointer :: tiles

  integer :: n_tils, sp_id, til_id, i, j, nchunks_max_init, n_tils_recv
  integer, dimension(:), pointer     :: tils_recv
  integer, dimension(:,:), pointer   :: neighbor_til_id
  integer, dimension(:,:,:), pointer :: shift
  class(t_species), pointer :: spec
  class(t_vdf), pointer :: jay, e, b
  logical :: if_last_til

  type( t_grid_arr ), dimension(:), pointer :: grid_arr

  ! call superclass, send recv tiles via MPI as normal
  call this % t_simulation_group % send_recv( patt, tiles )

  ! copy newly recvd tiles from host to device (as at initializatipon)

  n_tils = this%til_id(3)

  select type( no_co => this%no_co )
  class is( t_node_conf_cuda )
    nchunks_max_init = no_co%nchunks_max_init
  end select

  call alloc( neighbor_til_id, (/ 3**p_x_dim, n_tils /) )
  call alloc( shift, (/ p_x_dim, 3**p_x_dim, n_tils /) )

  allocate( grid_arr( n_tils ) )
  do i = 1, n_tils
    grid_arr(i)%g => this%tiles(i+this%til_id(1)-1)%t%grid
  enddo

  select type( no_co => this%no_co )
  class is( t_node_conf_tiles )
    call populate_neighbors( no_co, grid_arr, this%til_id, neighbor_til_id, shift )
  end select
  deallocate(grid_arr)

  call init_node_conf_dlb( n_tils, this%no_co%my_aid_, this%no_co%no_num, &
                           c_loc(neighbor_til_id), c_loc(shift) )

  call freemem( neighbor_til_id )
  call freemem( shift )

  ! note: chunk_pool, mem_transfer_manager, and group fields and spec are already g2g 
  ! (good to go)

  ! initialize cuda objects, tile by tile, for newly received tiles only
  n_tils_recv = 0
  do i = 1, patt%n_recv
    n_tils_recv = n_tils_recv + size( patt%recv(i)%tils )
  enddo

  if ( n_tils_recv > 0 ) then

    call alloc( tils_recv, (/ n_tils_recv /) )

    n_tils_recv = 0
    if_last_til = .false.
    do i = 1, patt%n_recv
      do j = 1, size( patt%recv(i)%tils )
        ! some special initialization things need to happen once we've passed all the tiles to cuda
        if ( i==patt%n_recv .and. j==size(patt%recv(i)%tils) ) if_last_til = .true.

        til_id = patt%recv(i)%tils(j)

        ! pack til_id for recv'd tiles, 0 indexed for use in C
        n_tils_recv = n_tils_recv + 1
        tils_recv(n_tils_recv) = til_id - this%til_id(1)

        jay => this%tiles(til_id)%t%jay%pf(1)
        e   => this%tiles(til_id)%t%emf%e_part
        b   => this%tiles(til_id)%t%emf%b_part
        call init_flds_tile( til_id-this%til_id(1)+1, c_loc(jay%nx_), c_loc(jay%gc_num_), &
                             c_loc(jay%buffer), c_loc(e%nx_), c_loc(e%gc_num_), &
                             c_loc(e%buffer), c_loc(b%buffer), logical(if_last_til,kind=c_bool) )
      enddo
    enddo

    call flush_mem_transfer_manager();

    if_last_til = .false.
    do i = 1, patt%n_recv
      do j = 1, size( patt%recv(i)%tils )
        ! some special initialization things need to happen once we've passed all the tiles to cuda
        if ( i==patt%n_recv .and. j==size(patt%recv(i)%tils) ) if_last_til = .true.

        til_id = patt%recv(i)%tils(j)

        sp_id = 1
        spec => this%tiles(til_id)%t%part%species
        do
          if (.not. associated(spec)) exit
          call init_spec_tile( til_id-this%til_id(1)+1, sp_id, spec%num_par_max, spec%num_par, &
                               c_loc(spec%energy), c_loc(spec%my_nx_p), c_loc(spec%x), &
                               c_loc(spec%p), c_loc(spec%q), c_loc(spec%ix), nchunks_max_init, &
                               logical(if_last_til,kind=c_bool) )
          spec => spec % next
          sp_id = sp_id + 1
        enddo
      enddo
    enddo

    call copy_particles_to_device_dlb( n_tils_recv, c_loc(tils_recv) )

    call freemem(tils_recv)

  else
    ! Complete the initialization for tiles
    call init_tile_last_dlb()
  endif

end subroutine send_recv_tils_dlb_cuda
!-----------------------------------------------------------------------------------------

end module m_sim_group_cuda
