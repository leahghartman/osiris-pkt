! #define DEBUG_FILE

module m_simulation_group

#include "os-config.h"
#include "os-preprocess.fpp"
#include "memory/memory.h"

use m_system
use m_parameters
use m_simulation, only : t_simulation, loopev, dynlbev, dynlb_set_int_load_ev, &
                         dynlb_new_lb_ev, dynlb_reshape_ev, dynlb_reshape_part_ev, &
                         dynlb_reshape_grid_ev, restart_write_ev, report_load_ev
use m_simulation_tiles, only : t_simulation_tiles
use m_input_file
use m_node_conf, only : t_node_conf, my_ngp, comm, root, n_threads, no_num
use m_node_conf_tiles
use m_tile_conf, only : p_topology_viktor, g_tgp, g_tgp_old, set_ti_co_old
use m_grid_define, only : t_grid, t_msg_patt
use m_grid_parallel, only : if_dynamic_lb, clear_int_load, use_part_load
use m_grid_tiles, only : t_grid_tiles, init_tile_load, add_load, conv_load_2_int_load, &
                         new_dlb_patt, t_grid_group
use m_restart
use m_random, only : init_genrand, list_algorithm_rnd
use m_space, only : xmin_initial, xmax_initial, x_bnd, if_move, xmax, xmin, move_boundary
use m_space, only : restart_write, setup
use m_time_step, only : n, dt, ndump, setup, advance, restart_write
use m_time, only : t, setup, tmin, update, restart_write, tmax
use m_emf_define, only : t_emf, p_extfld_none
use m_emf_diag, only : p_e1, p_e2, p_e3
use m_particles_define, only : t_particles
use m_particles, only : get_min_gc, get_grid_center, interpolation, &
                        ionize_neutral, sort_collide
use m_particles_tiles, only : add_load
use m_antenna_array, only : t_antenna_array, antenna_on
use m_bnd_tiles, only : t_bnd_group, t_bnd_patt
use m_diagnostic_utilities, only : cleanup_dutil
use m_species_comm_tiles, only : p_max_buffer1_size_t
use m_species_define, only : t_species
use m_logprof
use m_vdf_define, only : t_vdf_report
use m_dynamic_loadbalance, only : sim_imbalance
use m_random, only : rng, write_checkpoint_rnd


implicit none

private

! this extra type necessary to make an array of pointers in fortran
type :: t_til_container

  class( t_simulation_tiles ), pointer :: t => null()

end type t_til_container

type, extends( t_simulation ) :: t_simulation_group

  ! array of pointers to objects of type t_simulation
  class( t_til_container ), dimension(:), pointer :: tiles => null()

  ! Lower/upper/number global tile id's for this node
  integer, dimension(3) :: til_id

  ! input file pos
  integer :: pos

  ! Boundary information
  class( t_bnd_group ), pointer :: bnd_grp => null()
  ! bnd_patt will be allocated with dimension p_x_dim
  class( t_bnd_patt ), dimension(:), pointer :: bnd_patt => null()

  ! diagnostic information
  type( t_msg_patt ) :: diag_patt


contains

  procedure :: read_input       => read_input_sim_group
  procedure :: allocate_objs    => allocate_objs_group
  procedure :: allocate_tiles   => allocate_tiles_
  procedure :: init             => init_sim_group
  procedure :: list_algorithm   => list_algorithm_sim_group
  procedure :: iter             => iter_sim_group
  procedure :: cleanup          => cleanup_sim_group
  procedure :: write_checkpoint => write_checkpoint_group

  procedure :: advance_deposit  => advance_deposit_group
  procedure :: move_part_window => move_window_part_group

  ! Update boundary procedures
  procedure :: update_boundary_emf     => update_boundary_emf_group
  procedure :: update_boundary_part    => update_boundary_part_group
  procedure :: update_boundary_current => update_boundary_current_group
  procedure :: update_charge           => update_charge_part_group

  ! Diagnostic procedures
  procedure :: report_part        => report_part_group
  procedure :: report_emf         => report_emf_group
  procedure :: report_current     => report_current_group
  procedure :: init_dutil         => init_diag_utils_group
  procedure :: new_diag_patt
  procedure :: report_energy_emf  => report_energy_emf_group
  procedure :: report_energy_part => report_energy_part_group

  ! Load balance procedures
  procedure :: allocate_tiles_dlb => allocate_tiles_dlb_
  procedure :: init_dlb           => init_dlb_group
  procedure :: set_int_load       => set_int_load_group
  procedure :: dlb                => dlb_group
  procedure :: reshape            => reshape_group
  procedure :: send_recv          => send_recv_tils_dlb
  procedure :: report_load        => report_load_group
  procedure :: new_omp_patt
  procedure :: new_bnd_patt

end type t_simulation_group

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

interface
  subroutine update_boundary_emf_group( this )
    import t_simulation_group
    class(t_simulation_group), intent(inout) :: this
  end subroutine
end interface

interface
  subroutine update_boundary_part_group( this )
    import t_simulation_group
    class(t_simulation_group), intent(inout) :: this
  end subroutine
end interface

interface
  subroutine update_boundary_current_group( this )
    import t_simulation_group
    class(t_simulation_group), intent(inout) :: this
  end subroutine
end interface

interface
  subroutine update_charge_part_group( this )
    import t_simulation_group
    class(t_simulation_group), intent(inout) :: this
  end subroutine
end interface

interface
  subroutine report_part_group( this )
    import t_simulation_group
    class(t_simulation_group), intent(inout) :: this
  end subroutine
end interface

interface
  subroutine report_emf_group( sim )
    import t_simulation_group
    class( t_simulation_group ), intent(inout) :: sim
  end subroutine
end interface

interface
  subroutine init_diag_utils_group( sim, dx )
    import t_simulation_group, p_double
    class( t_simulation_group ), intent(inout) :: sim
    real(p_double), dimension(:), intent(in) :: dx
  end subroutine
end interface

interface
  subroutine new_diag_patt( sim )
    import t_simulation_group
    class( t_simulation_group ), intent(inout) :: sim
  end subroutine
end interface

interface
  subroutine report_current_group( sim )
    import t_simulation_group
    class( t_simulation_group ), intent(inout), target :: sim
  end subroutine
end interface

interface
  subroutine report_energy_emf_group( sim )
    import t_simulation_group
    class( t_simulation_group ), intent(inout) :: sim
  end subroutine
end interface

interface
  subroutine report_energy_part_group( sim )
    import t_simulation_group
    class ( t_simulation_group ), intent(inout) :: sim
  end subroutine
end interface

interface
  subroutine send_recv_tils_dlb( this, patt, tiles )
    import t_simulation_group, t_msg_patt, t_til_container
    class( t_simulation_group ), intent(inout)      :: this
    type( t_msg_patt ), intent(inout)               :: patt
    class( t_til_container ), dimension(:), pointer :: tiles
  end subroutine
end interface

interface
  subroutine report_load_group( this )
    import t_simulation_group
    class(t_simulation_group), intent(inout) :: this
  end subroutine
end interface

interface
  subroutine new_bnd_patt( this )
    import t_simulation_group
    class(t_simulation_group), intent(inout) :: this
  end subroutine
end interface

interface cleanup_tils_out
  module procedure cleanup_tils_out_group
end interface

integer :: emfboundev_glob, emfboundev_loc
integer :: partboundev_glob, partboundev_loc
integer :: jayboundev_glob, jayboundev_loc
integer :: omp_patt_ev
integer :: fsolver_outside_ev, push_outside_ev, smooth_outside_ev, sort_outside_ev

public :: t_simulation_group, t_til_container
public :: cleanup_tils_out
public :: emfboundev_glob, emfboundev_loc, partboundev_glob
public :: partboundev_loc, jayboundev_glob, jayboundev_loc
public :: omp_patt_ev
public :: fsolver_outside_ev, push_outside_ev, smooth_outside_ev, sort_outside_ev
public :: set_init_part

contains

!-----------------------------------------------------------------------------------------
! Read simulation parameters from input file
!-----------------------------------------------------------------------------------------
subroutine read_input_sim_group( sim, input_file )

  implicit none

  class ( t_simulation_group ), intent(inout) :: sim
  class( t_input_file ), intent(inout) :: input_file

  integer :: ierr, pos_nconf

  pos_nconf = input_file%pos
  ! Advance input file position to after node_conf for reading successive tiles
  call get_namelist( input_file, "nl_node_conf", ierr )

  ! Store input file position for later reference in allocate_tiles
  sim % pos = input_file%pos
  input_file%pos = pos_nconf

  call sim % t_simulation % read_input( input_file )

  ! Turn off displaying output while reading so that info is not printed for every tile
  call no_disp_out(input_file)

end subroutine read_input_sim_group

!-----------------------------------------------------------------------------------------
! Initialize simulation objects
!-----------------------------------------------------------------------------------------
subroutine allocate_objs_group( sim )

  implicit none

  class( t_simulation_group ), intent(inout) :: sim

  if ( .not. associated( sim%no_co ) ) then
    allocate( t_node_conf_tiles :: sim%no_co )
  endif

  if ( .not. associated( sim%grid ) ) then
    allocate( t_grid_group :: sim%grid )
  endif

  if ( .not. associated( sim%bnd_patt ) ) then
    allocate( sim%bnd_patt(p_x_dim) )
  endif

  if ( .not. associated( sim%bnd_grp ) ) then
    allocate( sim%bnd_grp )
  endif

  call sim % t_simulation % allocate_objs()

end subroutine allocate_objs_group

!-----------------------------------------------------------------------------------------
! Initialize simulation objects
!-----------------------------------------------------------------------------------------
subroutine allocate_tiles_( sim )

  implicit none

  class( t_simulation_group ), intent(inout) :: sim

  integer :: i

  ! Allocate tile containter
  if (.not. associated(sim%tiles)) then
    allocate( sim%tiles( sim%til_id(1):sim%til_id(2) ) )

    ! Allocate tiles
    do i = sim%til_id(1),sim%til_id(2)
      allocate( t_simulation_tiles :: sim%tiles(i)%t )
    enddo
  endif

  ! Allocate objects for tiles
  !$omp parallel do schedule(dynamic)
  do i=sim%til_id(1),sim%til_id(2)
    !$omp critical
    call sim%tiles(i)%t % allocate_objs()
    !$omp end critical
  enddo
  !$omp end parallel do

  do i=sim%til_id(1),sim%til_id(2)
    sim%input_file%pos = sim%pos
    call sim % tiles(i) % t % read_copy_sim( sim, sim%input_file, i )
  enddo

end subroutine allocate_tiles_
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Initialize new tiles array for dynamic load balance
!-----------------------------------------------------------------------------------------
subroutine allocate_tiles_dlb_( sim, patt, tiles_temp )

  implicit none

  class( t_simulation_group ), intent(inout)      :: sim
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
        allocate( t_simulation_tiles :: tiles_temp(til)%t )
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

end subroutine allocate_tiles_dlb_
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! List algorithm from the first tile only
!-----------------------------------------------------------------------------------------
subroutine list_algorithm_sim_group( sim )

  implicit none

  class( t_simulation_group ), intent(in) :: sim

  call list_algorithm( sim%tiles(sim%til_id(1))%t%options )
  call list_algorithm_rnd()
  call sim%tiles(sim%til_id(1))%t % emf % list_algorithm( )
  call sim%tiles(sim%til_id(1))%t % part % list_algorithm( )
  call sim%tiles(sim%til_id(1))%t % jay % list_algorithm( )

end subroutine list_algorithm_sim_group

!-----------------------------------------------------------------------------------------
! Simulation iteration
!-----------------------------------------------------------------------------------------
subroutine iter_sim_group( sim )

  implicit none

  class ( t_simulation_group ), intent(inout) :: sim

  class(t_emf), pointer :: emf
  class(t_node_conf_tiles), pointer :: no_co_tiles
  integer :: i, j

  ! set the initial timer for the loop
  call begin_event(loopev)

  ! advance time step counter and time
  DEBUG('Before advance time')
  call advance( sim%tstep, 1 )
  call update( sim%time, dt(sim%tstep), n(sim%tstep) )
  do i = sim%til_id(1), sim%til_id(2)
    call advance( sim%tiles(i)%t%tstep, 1 )
    call update( sim%tiles(i)%t%time, dt(sim%tiles(i)%t%tstep), n(sim%tiles(i)%t%tstep) )
  enddo

  if ( root(sim%no_co) ) then
     print '(A,I10,A,F0.3)', 'n = ', n(sim%tstep), ', t = ', t(sim%time)
  endif

  ! Electric current diagnostics  must occur before dynamic load
  ! balancing, since it destroys the current values
  DEBUG('Before report current')
  call sim % report_current()

  DEBUG('Before zpulse launch')
  do i = sim%til_id(1), sim%til_id(2)
    call sim%tiles(i)%t%zpulse_list%launch( sim%tiles(i)%t%emf, sim%g_space, &
      sim%tiles(i)%t%grid%my_nx(p_lower, :), sim%tiles(i)%t%grid%g_nx, &
      t(sim%tiles(i)%t%time), dt(sim%tiles(i)%t%tstep), sim%tiles(i)%t%no_co )
  enddo

  ! do this here because we may be setting fields -- want to inject them before the field solve
  DEBUG('Before delayed particle injection')
  do i = sim%til_id(1), sim%til_id(2)
    call sim%tiles(i)%t%part%delayed_injection( sim%jay, sim%emf, sim%grid, sim%no_co, sim%bnd, t(sim%time) )
  enddo

  ! ********************* move window ******************************

  ! test if%tiles(i) simulation window needs to be moved - use
  ! global space to avoid numerical inconsistencies
  ! between move on different nodes and tiles
  DEBUG('Before move window sim%g_space')
  emf => sim%tiles(sim%til_id(1))%t%emf
  call move_boundary( sim%g_space, emf%dx(), dt(sim%tiles(sim%til_id(1))%t%tstep) )

  ! move boundaries of other objects as required
  DEBUG('Before move window sim%emf')
  do i = sim%til_id(1), sim%til_id(2)
    call sim%tiles(i)%t%emf % move_window( sim%g_space, sim%tiles(i)%t%grid%my_nx( p_lower, : ), sim%tiles(i)%t%no_co%on_edge( 1, p_upper ) )
  enddo

  DEBUG('Before move window sim%part')
  call sim % move_part_window()

  ! jay does not need to be moved since it is recalculated
  ! from scratch during each interation (in advance_deposit)
  ! use antenna if necessary
  DEBUG('Before antenna')
  do i = sim%til_id(1), sim%til_id(2)
    if(antenna_on(sim%tiles(i)%t%antenna_array)) then
       call sim%tiles(i)%t%antenna_array % launch( sim%tiles(i)%t%emf, &
          dt(sim%tiles(i)%t%tstep), t(sim%tiles(i)%t%time), &
          sim%g_space, sim%tiles(i)%t%grid%my_nx( p_lower, : ))
    end if
  enddo
  ! update boundaries of em-fields and calculate
  ! time averaged fields if necessary
  DEBUG('Before update boundary sim%emf')
  call sim%update_boundary_emf()

  ! Update the fields used for particle interpolation. This includes spatial smoothing,
  ! external fields and subcycling
  DEBUG('Before update_particle_fld sim%emf')
  !$omp parallel do schedule(dynamic)
  do i = sim%til_id(1), sim%til_id(2)
    call sim%tiles(i)%t%emf%update_particle_fld( n(sim%tiles(i)%t%tstep), &
      t(sim%tiles(i)%t%time), dt(sim%tiles(i)%t%tstep), sim%g_space, &
      sim%tiles(i)%t%grid%my_nx( p_lower, : ))
  enddo
  !$omp end parallel do

  ! diagnostic of em-field and particles
  ! currents can only be diagnosed after advance_deposit and are taken care of at the
  ! beggining of the loop
  DEBUG('Before report emf')
  if ( sim%tiles(sim%til_id(1))%t%emf%if_charge_cons( sim%tiles(sim%til_id(1))%t%tstep ) ) then
    call sim % update_charge()
    ERROR( "Charge conservation diagnostic is temporarily disabled" )
    call abort_program()
    ! call update_charge_cons( sim%tiles(i)%t%emf, sim%tiles(i)%t%part%charge%current, sim%tiles(i)%t%bnd%send_vdf, sim%tiles(i)%t%bnd%recv_vdf )
  endif

  call sim%report_emf()

  DEBUG('Before report part')
  call sim % report_part()

  DEBUG('Before calculate omp_patt')
  call begin_event(omp_patt_ev)
  call sim % new_omp_patt()
  call end_event(omp_patt_ev)

  ! move particles a timestep further and get new current
  DEBUG('Before advance_deposit')
  call sim % advance_deposit()

  ! take care of boundaries for particles
  DEBUG('Before update boundary sim%part')
  call sim % update_boundary_part()

  ! report time centered energy
  call sim % report_energy_emf()
  call sim % report_energy_part()

  ! --------------------------------------------------
  ! begin injection of particles (cathode, ionization)
  ! --------------------------------------------------

  ! injection of particles must be done outside the
  ! move_boundary(part) <...> update_boundary(part) section
  ! of the code or else some paralell partitions will fail on
  ! moving window runs

  ! note that particles injected in this section will only
  ! be pushed by the EMF in the next iteration

  ! also note that particle injection here should contribute
  ! (if necessary) to the current (e.g. cathode)

  DEBUG('Before cathode')
  do i = sim%til_id(1), sim%til_id(2)
    call sim%tiles(i)%t%part%inject_cathode( sim%tiles(i)%t%jay, &
      t(sim%tiles(i)%t%time), dt(sim%tiles(i)%t%tstep), sim%tiles(i)%t%no_co, &
      sim%tiles(i)%t%grid%coordinates )
  enddo

  DEBUG('Before ionization')
  ! --------------------------------------------------
  ! end injection of particles
  do i = sim%til_id(1), sim%til_id(2)
    call ionize_neutral(sim%tiles(i)%t%part, sim%tiles(i)%t%emf, &
      dt(sim%tiles(i)%t%tstep), sim%tiles(i)%t%grid%coordinates)
  enddo
  ! --------------------------------------------------

  ! --------------------------------------------------
  ! begin particle sorting / collisions
  ! --------------------------------------------------

  ! sorting and collisions are done here, outside the
  ! move_boundary() <...> update_boundary() section
  ! of the loop to avoid some issues related to the moving
  ! window algorithm.
  ! Outside this section all the data (grids/particles)
  ! correctly placed in the local space
  ! sorting and collisions are done together since the
  ! collision algorithm requires sorting the particles

  DEBUG('Before sort')
  call begin_event(sort_outside_ev)
  select type(no_co => sim%no_co); class is(t_node_conf_tiles)
    no_co_tiles => no_co
    ! do omp over light tiles
    if (no_co_tiles%omp_patt%n_light > 0) then
      !$omp parallel do schedule(dynamic) private(i)
      do j = 1, no_co_tiles%omp_patt%n_light
        i = no_co_tiles%omp_patt%light_tils(j)
        call sort_collide( sim%tiles(i)%t%part, n(sim%tiles(i)%t%tstep), &
          t(sim%tiles(i)%t%time), n_threads(sim%tiles(i)%t%no_co) )
      enddo
      !$omp end parallel do
    endif

    ! do omp inside of heavy tiles
    do j = 1, no_co%omp_patt%n_heavy
      i = no_co%omp_patt%heavy_tils(j)
      call sort_collide( sim%tiles(i)%t%part, n(sim%tiles(i)%t%tstep), &
        t(sim%tiles(i)%t%time), n_threads(sim%tiles(i)%t%no_co) )
    enddo
    call end_event(sort_outside_ev)
  end select
  ! --------------------------------------------------
  ! end particle sorting / collisions
  ! --------------------------------------------------

  ! normalize the current using the volume elements for cylindrical coordinates
  ! this also takes care of axial boundary
  if ( sim%tiles(sim%til_id(1))%t%grid%coordinates == p_cylindrical_b ) then
    do i = sim%til_id(1), sim%til_id(2)
      DEBUG('Before normalize sim%jay')
      call sim%tiles(i)%t%jay%normalize_cyl( )
    enddo
  endif

  ! take care of boundaries for currents
  DEBUG('Before update boundary sim%jay')
  call sim%update_boundary_current()

  ! smooth currents
  ! - smoothing has to be after update boundary
  DEBUG('Before smooth sim%jay')
  call begin_event(smooth_outside_ev)
  !$omp parallel do schedule(dynamic)
  do i = sim%til_id(1), sim%til_id(2)
    call sim%tiles(i)%t % jay % smooth()
  enddo
  !$omp end parallel do
  call end_event(smooth_outside_ev)
#if 0
  ! This is currently disabled

  launch current source based EM waves
  do i = sim%til_id(1), sim%til_id(2)
    call sim%tiles(i)%t%zpulse%launch( sim%tiles(i)%t%jay, sim%g_space, &
      sim%tiles(i)%t%grid%my_nx(p_lower, :), sim%tiles(i)%t%grid%g_nx, &
      t(sim%tiles(i)%t%time), dt(sim%tiles(i)%t%tstep), sim%tiles(i)%t%no_co )
  enddo
#endif

  ! calculate new electro-magnetic fields
  DEBUG('Before advance sim%emf')
  call begin_event(fsolver_outside_ev)
  !$omp parallel do schedule(dynamic)
  do i = sim%til_id(1), sim%til_id(2)
    call sim%tiles(i)%t % emf % advance( sim%tiles(i)%t%jay, &
      dt(sim%tiles(i)%t%tstep), sim%tiles(i)%t%bnd%send_vdf, sim%tiles(i)%t%bnd%recv_vdf )
  enddo
  !$omp end parallel do
  call end_event( fsolver_outside_ev )

  ! set the final time for the loop
  call end_event(loopev)

end subroutine iter_sim_group
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Cleanup simulation objects
!-----------------------------------------------------------------------------------------
subroutine cleanup_sim_group( sim )

  implicit none

  class( t_simulation_group ), intent(inout) :: sim

  integer :: i

  ! cleanup regular simulation objects
  call sim % grid % cleanup()

  ! cleanup particles
  call sim % part % cleanup( )

  ! cleanup current object
  call sim % jay % cleanup()

  ! cleanup emf object
  call sim % emf % cleanup()

  ! cleanup laser pulses
  call sim % zpulse_list % cleanup( )

  ! cleanup antennas
  call sim % antenna_array % cleanup()

  ! cleanup vdf and spec comm buffers
  call sim % bnd % cleanup()

  ! cleanup diagnostic utilities
  call cleanup_dutil()
  call sim % diag_patt % cleanup()

  ! cleanup tiles
  do i = sim%til_id(1), sim%til_id(2)
    call sim%tiles(i)%t % cleanup()
  enddo

  call sim % bnd_grp % cleanup()

  do i = 1, p_x_dim
    call sim % bnd_patt(i) % cleanup()
  enddo

  call sim % no_co % cleanup()

  ! deallocate polymorphic objects
  deallocate( sim%input_file )
  deallocate( sim%no_co )
  deallocate( sim%emf )
  deallocate( sim%part )
  deallocate( sim%grid )
  deallocate( sim%jay )
  deallocate( sim%zpulse_list )
  deallocate( sim%bnd )
  do i=sim%til_id(1),sim%til_id(2)
    deallocate( sim%tiles(i)%t )
  enddo
  deallocate( sim%tiles ) ! tile container
  deallocate( sim%bnd_grp )
  deallocate( sim%bnd_patt )

end subroutine cleanup_sim_group
!---------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------
! Write checkpoint data
!---------------------------------------------------------------------------------------
subroutine write_checkpoint_group( sim, restart_handle )

  implicit none

  class( t_simulation_group ), intent(in) :: sim
  type( t_restart_handle ), intent(inout) :: restart_handle

  integer :: i

  call sim % no_co % write_checkpoint( restart_handle )
  call sim % grid % write_checkpoint( restart_handle )
  call write_checkpoint_rnd( rng  , restart_handle )
  call restart_write( sim%tstep   , restart_handle )
  call restart_write( sim%restart , restart_handle )
  call restart_write( sim%g_space , restart_handle )
  call restart_write( sim%time    , restart_handle )

  do i = sim%til_id(1), sim%til_id(2)
    call sim % tiles(i) % t % grid % write_checkpoint( restart_handle )
  enddo

  do i = sim%til_id(1), sim%til_id(2)
    call sim % tiles(i) % t % emf % write_checkpoint( restart_handle )
  enddo

  do i = sim%til_id(1), sim%til_id(2)
    call sim % tiles(i) % t % jay % write_checkpoint( restart_handle )
  enddo

  do i = sim%til_id(1), sim%til_id(2)
    call sim % tiles(i) % t % part % write_checkpoint( restart_handle )
  enddo

  do i = sim%til_id(1), sim%til_id(2)
    call sim % tiles(i) % t % antenna_array % write_checkpoint( restart_handle )
  enddo

end subroutine write_checkpoint_group
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
subroutine advance_deposit_group( sim )

  implicit none

  class( t_simulation_group ), intent(inout) :: sim

  class(t_node_conf_tiles), pointer :: no_co_tiles
  integer :: i, j

  call begin_event( push_outside_ev )

  select type(no_co => sim%no_co); class is(t_node_conf_tiles)
    no_co_tiles => no_co
    ! do omp over light tiles
    if (no_co_tiles%omp_patt%n_light > 0) then
      !$omp parallel do schedule(dynamic) private(i)
      do j = 1, no_co_tiles%omp_patt%n_light
        i = no_co_tiles%omp_patt%light_tils(j)
        call sim%tiles(i)%t % part % advance_deposit( sim%tiles(i)%t%emf, &
          sim%tiles(i)%t%jay, sim%tiles(i)%t%tstep, t(sim%tiles(i)%t%time), &
          sim%tiles(i)%t%no_co, sim%tiles(i)%t%options )
      enddo
      !$omp end parallel do
    endif

    ! do omp inside of heavy tiles
    do j = 1, no_co%omp_patt%n_heavy
      i = no_co%omp_patt%heavy_tils(j)
      call sim%tiles(i)%t % part % advance_deposit( sim%tiles(i)%t%emf, &
        sim%tiles(i)%t%jay, sim%tiles(i)%t%tstep, t(sim%tiles(i)%t%time), &
        sim%tiles(i)%t%no_co, sim%tiles(i)%t%options )
    enddo
  end select

  call end_event( push_outside_ev )

end subroutine advance_deposit_group
!-------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine move_window_part_group( sim )

  implicit none

  class( t_simulation_group ), intent(inout) :: sim

  integer :: i

  do i = sim%til_id(1), sim%til_id(2)
    call sim%tiles(i)%t%part % move_window( sim%g_space, sim%tiles(i)%t%grid, &
      sim%tiles(i)%t%jay, sim%tiles(i)%t%no_co, sim%tiles(i)%t%bnd )
  enddo

end subroutine move_window_part_group
!-----------------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Sets the load on each tile
! - The result is saved on a grid object, which may not be the main grid object
!-------------------------------------------------------------------------------
subroutine set_int_load_group( sim, grid, n_ )

  implicit none

  class( t_simulation_group ), intent(inout) :: sim
  class( t_grid ), intent(inout) :: grid
  integer, intent(in) :: n_

  integer :: i, topology_type

  call begin_event(dynlb_set_int_load_ev)

  ! initialize int_load array
  call init_tile_load( grid, sim%no_co )

  select type(no_co => sim%no_co); class is(t_node_conf_tiles)
    topology_type = no_co%ti_co%topology_type
  end select

  select type(grid); class is(t_grid_group)

  ! add particle load if required
  ! density data will be used for iteration 0
  if ( use_part_load( grid ) ) then
    if ( n_ == 0 ) then
      ! We do the density calculation once, not for every tile
      call add_load( sim % part, grid, sim%no_co, n_ )
    else
      ! Every tile needs to report its number of particles
      do i=sim%til_id(1),sim%til_id(2)
        call add_load( sim%tiles(i)%t%part, grid, sim%tiles(i)%t%no_co, n_ )
      enddo
    endif
  endif

  ! add grid / expression load if required
  if ( n_ == 0 ) then
    call add_load( grid, sim%no_co, n_ )
  else
    ! Every tile needs to report its number of cells
    do i=sim%til_id(1),sim%til_id(2)
      select case( topology_type )
      case(p_topology_viktor)
        select case(p_x_dim)
        case(2)
          call add_load( sim%tiles(i)%t%grid, sim%tiles(i)%t%no_co, n_, &
            til_load_2=grid%til_load_2 )
        case(3)
          call add_load( sim%tiles(i)%t%grid, sim%tiles(i)%t%no_co, n_, &
            til_load_3=grid%til_load_3 )
        end select
      case default
        call add_load( sim%tiles(i)%t%grid, sim%tiles(i)%t%no_co, n_, &
          int_load=grid%int_load )
      end select
    enddo
  endif

  end select

  ! gather results from all nodes
  call grid % gather_load( sim%no_co )

  ! convert to integral
  if (topology_type/=p_topology_viktor) then
    call conv_load_2_int_load( grid%int_load )
  endif

  call end_event(dynlb_set_int_load_ev)

end subroutine set_int_load_group
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!  dynamically load balance the simulation (for tiles algorithm)
!-------------------------------------------------------------------------------
subroutine dlb_group( this )

  implicit none

  class( t_simulation_group ), intent(inout) ::  this

  class( t_grid ), pointer :: new_lb
  logical :: do_dlb

  if (if_dynamic_lb(this%grid, n(this%tstep), this%no_co)) then

    DEBUG('Attempting dynamic load balance')

    call begin_event(dynlbev)

    ! Check if the simulation imbalance is over the set threshold. (The imbalance is defined as
    ! the ratio between the maximum load and the average load)
    if ( this%grid%max_imbalance >= 1 ) then
      do_dlb = sim_imbalance( this%grid, this%part, this%no_co ) > this%grid%max_imbalance
    else
      do_dlb = .true.
    endif

    if ( do_dlb ) then

      if ( root( this%no_co ) ) then
        print *, ''
        print *, ' Dynamically load balancing the simulation.'
      endif

      allocate( t_grid_group :: new_lb )

      ! copy data to new_lb
      call new_lb % copy( this%grid )

      ! determine the current integral particle load
      call this % set_int_load( new_lb, n(this%tstep) )

      ! Save old tile configuration data
      select type(no_co => this%no_co); class is(t_node_conf_tiles)
        call set_ti_co_old( no_co%ti_co )
      end select

      ! determine optimal distribution of tiles on nodes
      call begin_event( dynlb_new_lb_ev )
      select type(new_lb); class is(t_grid_group)
        call new_lb % parallel_partition_group( this%no_co, n(this%tstep) )
      end select
      call end_event( dynlb_new_lb_ev )

      ! clear the memory used by the integral particle load
      call clear_int_load( new_lb )

! #if 0
!       ! (* debug *)

!       call wait_for_all( this%no_co  )
!       SCR_ROOT("")
!       SCR_ROOT(" Verifying dynamic load balance (A)...")

!       ! (* debug *)
!       ! Check if reshaping a vdf with the same shape as the e field works ok
!       call test_reshape_vdf( this%emf%e, this%grid, new_lb, this%no_co )

!       ! check particles are ok for debug purposes
!       call validate( this%part, "Veryfying dynamic load balance (A)" )

! #endif

      call begin_event( dynlb_reshape_ev )

      ! Load balance the simulation: exchange tiles btwn nodes
      call this % reshape()

      call end_event( dynlb_reshape_ev )

      ! Cleanup temp. objects
      call new_lb % cleanup()
      deallocate( new_lb )

      ! synchronize all nodes
      ! Doesn't seem to be necessary, remove wait_for_all
      ! call wait_for_all( this%no_co  )

! #if 0
!       ! (* debug *)

!       SCR_ROOT("")
!       SCR_ROOT(" Verifying dynamic load balance (B)...")

!       ! check particles are ok for debug purposes
!       call validate( this%part, "Veryfying dynamic load balance (B)" )

!       ! check all nodes have the same grid object
!       call check_consistency( this%grid, this%no_co )

!       ! synchronize all nodes for debug purposes
!       call wait_for_all( this%no_co  )
!       SCR_ROOT(" No errors found!")
!       SCR_ROOT("")

!       ! (* end debug *)
! #endif

      if ( root( this%no_co ) ) then
        print *, ''
      endif

    endif

    call end_event(dynlbev)

  endif

end subroutine dlb_group
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! move tiles to proper process
!-------------------------------------------------------------------------------
subroutine reshape_group( this )

  implicit none

  class( t_simulation_group ), intent(inout) :: this

  class( t_til_container ), dimension(:), pointer :: tiles_temp
  type( t_msg_patt )                              :: msg_patt

  integer, dimension(:), pointer :: diff1
  integer, dimension(:,:), pointer :: diff2
  integer, dimension(:,:,:), pointer :: diff3

  select type(no_co => this%no_co); class is(t_node_conf_tiles)

  ! check if load balance needs to be done
  ! TODO: We should check if there is a significant global load balance
  ! advantage and only change it in that case
  if ( all(no_co%ti_co%g_til_aid_min_max_old==no_co%ti_co%g_til_aid_min_max) ) then
    SCR_ROOT( ' Skipping dynamic load balance' )
    return
  endif

  ! ! (*debug*)
  ! if (root(no_co)) then
  !   print *, '--g_til_aid_min_max_old'
  !   print '(3I6)', no_co%ti_co%g_til_aid_min_max_old
  !   print *, '--g_til_aid_min_max'
  !   print '(3I6)', no_co%ti_co%g_til_aid_min_max
  ! endif

  select case (no_co%ti_co%topology_type)
  case (p_topology_viktor)
    ! ! (*debug*)
    ! if ( root(no_co) ) then
    !   do i = 1, p_x_dim
    !     print *, 'dim',i
    !     print *, 'cart_til_mm_old'
    !     print '(3I6)', no_co%ti_co%cart_til_mm_old(i)%ii
    !     print *, 'cart_til_mm'
    !     print '(3I6)', no_co%ti_co%cart_til_mm(i)%ii
    !   enddo
    ! endif
    select case ( p_x_dim )
    case (2)
      call alloc(diff2, (/ no_co%ti_co%g_til(1), no_co%ti_co%g_til(2) /) )
      ! compute diff2
      call prep_dlb_v( this, diff=diff2 )
      ! fill message pattern
      call new_dlb_patt( this%no_co, msg_patt, diff=diff2 )
      ! initialize new tile array
      call this % init_dlb( msg_patt, tiles_temp, diff=diff2 )
      call freemem(diff2)
    case (3)
      call alloc(diff3, (/ no_co%ti_co%g_til(1), no_co%ti_co%g_til(2), no_co%ti_co%g_til(3) /) )
      ! compute diff3
      call prep_dlb_v( this, difff=diff3 )
      ! fill message pattern
      call new_dlb_patt( this%no_co, msg_patt, difff=diff3 )
      ! initialize new tile array
      call this % init_dlb( msg_patt, tiles_temp, difff=diff3 )
      call freemem(diff3)
    end select
  case default
    call alloc(diff1, (/ product(no_co%ti_co%g_til) /) )
    ! prepare local variables
    call prep_dlb( this, diff1 )
    ! get new message pattern: outgoing/incoming tiles and their destinations
    call new_dlb_patt( msg_patt, this%no_co, no_co%ti_co%g_til_aid_min_max_old, &
      no_co%ti_co%g_til_aid_min_max, diff1 )
    ! initialize new tile array
    call this % init_dlb( msg_patt, tiles_temp, dif=diff1 )
    call freemem(diff1)
  end select

  end select

  ! move tiles to proper process, cleanup outgoing tiles, deallocate this%tiles, and
  ! make temp tile array the new tile array
  call this % send_recv( msg_patt, tiles_temp )

  ! Set up tile boundary information
  call this % new_bnd_patt()
  call this % bnd_grp % init( this%bnd_patt )

  ! have lunch

  ! cleanup
  call msg_patt % cleanup()

end subroutine reshape_group
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! prepare variables to be used in dynamic load balance, viktor's algorithm
!-------------------------------------------------------------------------------
subroutine prep_dlb_v( this, diff, difff )

  implicit none

  class( t_simulation_group ), intent(in)            :: this
  integer, dimension(:,:), intent(inout), optional   :: diff
  integer, dimension(:,:,:), intent(inout), optional :: difff

  integer :: g_til_num, num_nodes
  integer :: i, node
  integer, dimension(3) :: lb, ub

  select type(no_co => this%no_co); class is(t_node_conf_tiles)

  g_til_num = product( no_co%ti_co%g_til )
  num_nodes = no_num( no_co )

  lb = 1
  do i = 1, p_x_dim
    ub(i) = no_co%ti_co%g_til( i )
  enddo

  select case ( p_x_dim )
    case (2)
      ! call alloc( diff, lb, ub )

      ! calculate diff
      do i = 1, product( no_co%ti_co%g_til )
        ! find node that g_til_aid=i belongs to
        node = node_of_til( no_co, i )
        ! set diff of tiles g_tgp
        diff( g_tgp(no_co%ti_co,i,1), g_tgp(no_co%ti_co,i,2) ) = node
      enddo
      do i = 1, product( no_co%ti_co%g_til )
        ! find node that g_til_aid=i used to belong to
        node = node_of_til_old( no_co, i )
        ! set diff of tile's old g_tgp to difference btwn new node id, old tile id
        diff( g_tgp_old(no_co%ti_co,i,1), g_tgp_old(no_co%ti_co,i,2) ) &
          = diff( g_tgp_old(no_co%ti_co,i,1), g_tgp_old(no_co%ti_co,i,2) ) - node
      enddo

      ! if ( root(no_co) ) then
      !   ! print tile topology
      !   call vis_hid(no_co%ti_co)
      !   print *, '--'
      !   ! print node tile topology
      !   print *, 'diff:'
      !   do i = g_til(no_co,2),1,-1
      !     print *, diff(:,i)
      !   enddo
      ! endif
    case (3)
      ! call alloc( difff, lb, ub )

      ! calculate diff
      do i = 1, product(no_co%ti_co%g_til)
        ! find node that g_til_aid=i belongs to
        node = node_of_til( no_co, i )
        ! fill tile's element of inv. space filling op. with node of tile
        difff( g_tgp(no_co%ti_co,i,1), g_tgp(no_co%ti_co,i,2), g_tgp(no_co%ti_co,i,3) ) = node
      enddo
      do i = 1, product(no_co%ti_co%g_til)
        ! find node that g_til_aid=i belongs to
        node = node_of_til_old( no_co, i )
        ! fill tile's element of inv. space filling op. with node of tile
        difff( g_tgp_old(no_co%ti_co,i,1), g_tgp_old(no_co%ti_co,i,2), &
          g_tgp_old(no_co%ti_co,i,3) ) = difff( g_tgp_old(no_co%ti_co,i,1), &
          g_tgp_old(no_co%ti_co,i,2), g_tgp_old(no_co%ti_co,i,3) ) - node
      enddo

      ! if ( root(no_co) ) then
      !   ! print tile topology
      !   call vis_hid(no_co%ti_co)
      !   print *, '--'
      !   ! print node tile topology
      !   print *, 'diff:'
      !   do j = g_til(no_co,3),1,-1
      !     print *, 'k=',j
      !     do i = g_til(no_co,2),1,-1
      !       print *, difff(:,i,j)
      !     enddo
      !   enddo
      ! endif
  end select

  end select

end subroutine prep_dlb_v
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! prepare variables to be used in dynamic load balance
!-------------------------------------------------------------------------------
subroutine prep_dlb( this, dif )

  implicit none

  class( t_simulation_group ), intent(in) :: this
  integer, dimension(:), intent(inout)    :: dif

  integer, dimension(:), pointer :: tils_node_old, tils_node_new

  integer :: i, lb, ub

  select type(no_co => this%no_co); class is(t_node_conf_tiles)

  allocate(tils_node_old(product(no_co%ti_co%g_til)))
  allocate(tils_node_new(product(no_co%ti_co%g_til)))

  do i = 1, no_num( no_co )
    lb = no_co%ti_co%g_til_aid_min_max_old(1,i)
    ub = no_co%ti_co%g_til_aid_min_max_old(2,i)
    tils_node_old(lb:ub) = i

    lb = no_co%ti_co%g_til_aid_min_max(1,i)
    ub = no_co%ti_co%g_til_aid_min_max(2,i)
    tils_node_new(lb:ub) = i
  enddo

  dif = tils_node_new - tils_node_old

  deallocate(tils_node_old)
  deallocate(tils_node_new)

  end select

end subroutine prep_dlb
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! cleanup tiles leaving node in dynamic load balance
!-------------------------------------------------------------------------------
subroutine cleanup_tils_out_group( this, msg_patt )

  implicit none

  class( t_simulation_group ), intent(inout) :: this
  type( t_msg_patt ), intent(in)             :: msg_patt

  integer :: i, j, k
  integer :: til
  integer, dimension(p_x_dim) :: tgp

  select type(no_co => this%no_co); class is(t_node_conf_tiles)

    select case (no_co%ti_co%topology_type)
    case (p_topology_viktor)
      if ( msg_patt%n_send>0 ) then
        do i = 1, msg_patt%n_send
          do j = 1, size(msg_patt%send(i)%tils)
            ! g_tgp of outoing tile
            do k = 1, p_x_dim
              tgp(k) = g_tgp( no_co%ti_co, msg_patt%send(i)%tils(j), k )
            enddo
            ! this is the pre load balance g_til_aid of this outgoing tile
            if ( p_x_dim==2 ) til = no_co%ti_co%hi2_old( tgp(1), tgp(2) )
            if ( p_x_dim==3 ) til = no_co%ti_co%hi3_old( tgp(1), tgp(2), tgp(3) )
            ! cleanup the outgoing tile
            call this%tiles(til)%t % cleanup()
            deallocate( this%tiles(til)%t )
          enddo
        enddo
      endif
    case default
      if ( msg_patt%n_send>0 ) then
        do i = 1, msg_patt%n_send
          do j = 1, size(msg_patt%send(i)%tils)
            til = msg_patt%send(i)%tils(j)
            call this%tiles(til)%t % cleanup()
            deallocate( this%tiles(til)%t )
          enddo
        enddo
      endif
    end select

  end select

end subroutine cleanup_tils_out_group
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! When initializing tiles in dlb, don't initialize particles
!-------------------------------------------------------------------------------
subroutine set_init_part( part, bool )

  implicit none

  class( t_particles ), intent(inout) :: part
  logical, intent(in) :: bool

  class( t_species ), pointer :: spec

  spec => part%species
  do
    if (.not. associated(spec)) exit
    spec % init_part = bool
    spec => spec % next
  enddo

end subroutine set_init_part
!-------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Determine tiles with load greater than/less than the total load per thread
! Heavy and light tiles will then be processed differently with openmp parallelism
!-----------------------------------------------------------------------------------------
subroutine new_omp_patt( this )

  implicit none

  class(t_simulation_group), intent(inout) :: this

  class(t_species), pointer :: species

  integer, dimension(this%til_id(1):this%til_id(2)) :: til_load
  integer :: n_threads, n_heavy, n_light
  integer :: tot_load
  real(p_double) :: threshold
  integer :: i, t

  ! use extended version of t_node_conf, t_node_conf_tiles
  select type(no_co => this%no_co); class is(t_node_conf_tiles)

    ! allocate space to hold g_til_aid's of heavy and light tiles
    if ( associated(no_co%omp_patt%heavy_tils) ) then
      ! free and reallocate if not large enough
      if ( size(no_co%omp_patt%heavy_tils) < this%til_id(3) ) then
        call freemem( no_co%omp_patt%heavy_tils )
        call freemem( no_co%omp_patt%light_tils )
        call alloc( no_co%omp_patt%heavy_tils, (/this%til_id(3)/) )
        call alloc( no_co%omp_patt%light_tils, (/this%til_id(3)/) )
      endif
    else
      ! allocate if not allocated
      call alloc( no_co%omp_patt%heavy_tils, (/this%til_id(3)/) )
      call alloc( no_co%omp_patt%light_tils, (/this%til_id(3)/) )
    endif

    ! number of threads available to the node
    n_threads = no_co%n_threads

    select case( no_co%omp_patt%force_omp )
    case( p_heavy )

      n_heavy = this%til_id(3)
      n_light = 0

    case( p_light )

      n_heavy = 0
      n_light = this%til_id(3)

    case( p_default )

      ! get total number of particles on the node
      tot_load = 0
      til_load = 0
      do t = this%til_id(1), this%til_id(2)
        species => this%tiles(t)%t%part%species
        do
          if (.not. associated(species)) exit
          tot_load = tot_load + species%num_par
          til_load(t) = til_load(t) + species%num_par
          species => species % next
        enddo
      enddo

      ! calculate num tiles with load greater than/less than the total load per thread
      n_heavy = 0
      n_light = 0
      threshold = real(tot_load,p_double)/real(n_threads,p_double)

      do t = this%til_id(1), this%til_id(2)
        if ( til_load(t) > threshold ) then
          n_heavy = n_heavy + 1
        else
          n_light = n_light + 1
        endif
      enddo

      ! if have fewer light tils than threads process all tils as heavy tils
      if ( n_light < n_threads ) then
        n_heavy = this%til_id(3)
        n_light = 0
      endif

    end select

    ! fill heavy and light til g_til_aid's in omp_patt
    if ( n_heavy==0 ) then
      n_light = 0
      do t = this%til_id(1), this%til_id(2)
        n_light = n_light + 1
        no_co%omp_patt%light_tils(n_light) = t
      enddo
    elseif ( n_light==0 ) then
      n_heavy = 0
      do t = this%til_id(1), this%til_id(2)
        n_heavy = n_heavy + 1
        no_co%omp_patt%heavy_tils(n_heavy) = t
      enddo
    else
      n_heavy = 0
      n_light = 0
      do t = this%til_id(1), this%til_id(2)
        if ( til_load(t) > threshold ) then
          n_heavy = n_heavy + 1
          no_co%omp_patt%heavy_tils(n_heavy) = t
        else
          n_light = n_light + 1
          no_co%omp_patt%light_tils(n_light) = t
        endif
      enddo
    endif

    ! set n_threads to tot number of threads available for heavy tils
    do i = 1, n_heavy
      t = no_co%omp_patt%heavy_tils(i)
      this%tiles(t)%t%no_co%n_threads = n_threads
    enddo
    ! set n_threads to 1 for light tils
    do i = 1, n_light
      t = no_co%omp_patt%light_tils(i)
      this%tiles(t)%t%no_co%n_threads = 1
    enddo

    ! store number of heavy/light tiles for future reference
    no_co%omp_patt%n_heavy = n_heavy
    no_co%omp_patt%n_light = n_light

  end select

end subroutine new_omp_patt
!-----------------------------------------------------------------------------------------


end module m_simulation_group


!-----------------------------------------------------------------------------------------
! Initialize simulation objects
!-----------------------------------------------------------------------------------------
subroutine init_sim_group( sim, grid_gc_min )

  use m_simulation_group
  use m_parameters
  use m_simulation, only : t_simulation, loopev, dynlbev, dynlb_set_int_load_ev, &
                           dynlb_new_lb_ev, dynlb_reshape_ev, dynlb_reshape_part_ev, &
                           dynlb_reshape_grid_ev, restart_write_ev, report_load_ev
  use m_node_conf, only : t_node_conf, comm
  use m_node_conf_tiles
  use m_grid_define, only : t_grid, t_msg_patt
  use m_grid_parallel, only : if_dynamic_lb, clear_int_load, use_part_load
  use m_grid_tiles, only : t_grid_tiles, init_tile_load, add_load, conv_load_2_int_load, &
                           new_dlb_patt, t_grid_group
  use m_restart
  use m_random, only : init_genrand, list_algorithm_rnd
  use m_space, only : xmin_initial, xmax_initial, x_bnd, if_move, xmax, xmin, move_boundary
  use m_space, only : restart_write, setup
  use m_time_step, only : n, dt, ndump, setup, advance, restart_write
  use m_time, only : t, setup, tmin, update, restart_write, tmax
  use m_emf_define, only : t_emf, p_extfld_none
  use m_particles, only : get_min_gc, get_grid_center, interpolation, &
                          ionize_neutral, sort_collide
  use m_species_comm_tiles, only : p_max_buffer1_size_t
  use m_species_define, only : t_species
  use m_logprof

  implicit none

  class( t_simulation_group ), intent(inout) ::  sim
  integer, dimension( 2, p_x_dim ), intent(in) :: grid_gc_min
  
  real( p_double ), dimension( p_x_dim ) :: ldx
  integer :: i_dim, i
  integer, dimension( 2, p_x_dim ) :: min_gc, min_gc_tmp
  logical :: restart
  type( t_restart_handle )    ::  restart_handle
  class( t_species ), pointer :: species

  ! these pointers are added to fix a compiler error in GCC 7 wherein GCC could not
  !   correctly resolve long deference chains that ended in functions that return arrays.
  class( t_emf ), pointer :: emf
  class(t_grid), pointer :: grid
  class(t_node_conf), pointer :: no_co
  class(t_simulation), pointer :: simulation


  ! initialize timing structure
  loopev                = create_event('loop')
  dynlbev               = create_event('dynamic load balance (total)')
  dynlb_set_int_load_ev = create_event('dynlb get int load')
  dynlb_new_lb_ev       = create_event('dynlb get new lb')
  dynlb_reshape_ev      = create_event('dynlb reshape simulation')
  dynlb_reshape_part_ev = create_event('dynlb reshape simulation (particles)')
  dynlb_reshape_grid_ev = create_event('dynlb reshape simulation (grids)')
  restart_write_ev      = create_event('writing restart information')
  report_load_ev        = create_event('load diagnostics')

  emfboundev_glob = create_event('emf bound glob')
  emfboundev_loc = create_event('emf bound loc')
  partboundev_glob = create_event('part bound glob')
  partboundev_loc = create_event('part bound loc')
  jayboundev_glob = create_event('jay bound glob')
  jayboundev_loc = create_event('jay bound loc')
  omp_patt_ev = create_event('heavy/light tile calculation (omp_patt)')

  fsolver_outside_ev = create_event('field solver outside omp loop')
  push_outside_ev = create_event('advance deposit outside omp loop')
  smooth_outside_ev = create_event('current smooth outside omp loop')
  sort_outside_ev = create_event('particle sort (total) outside omp loop')

  ! use restart files if restart is set in the input file
  restart = if_restart_read(sim%restart)

  SCR_ROOT('')
  if ( restart ) then
    SCR_ROOT('**********************************************')
    SCR_ROOT('* Restarting simulation from checkpoint data *')
    SCR_ROOT('**********************************************')
  else
    SCR_ROOT('Initializing simulation:')
  endif

  !---------------------------------------------------------------------------------------
  ! Set up group objects
  !---------------------------------------------------------------------------------------

  SCR_ROOT(" Setting up group no_co")
  call sim % no_co % init()

  if ( restart ) then
    ! open restart file (this needs to happen after setup no_co)
    call restart_read_open( sim%restart, comm(sim%no_co), sim%no_co%ngp_id(), file_id_rst, restart_handle )

    ! read restart information for no_co
    ! This used to only check for inconsistencies in the input file vs. restart files,
    ! but now sets extra parameters in no_co
    call sim % no_co % restart_read( restart_handle )
  endif

  ! setup the sim group grid - want superclass init
  SCR_ROOT(" Setting up group grid")
  call sim % grid % init( sim%no_co, grid_gc_min, x_bnd( sim%g_space ), restart, &
    restart_handle )

  ! if not restarting generate initial load balance partition
  if ( .not. restart ) then

    if ( sim % grid % needs_int_load( 0 ) ) then
      SCR_ROOT(' - generating initial integral tile load...')
      ! create even partition for calculating particle (and grid) loads
      call sim % grid % parallel_partition( sim%no_co, -2 )
      call sim % set_int_load( sim%grid, 0 )
    endif

    SCR_ROOT(' - setting group grid boundaries...')
    select type(grid => sim%grid); class is(t_grid_group)
      call grid % parallel_partition_group( sim%no_co, 0 )
    end select
    call clear_int_load( sim%grid )

  endif

  SCR_ROOT(" - setting up group grid boundaries for diagnostics...")
  ! Initialize member data of group grid needed for diagnostics
  select type(grid => sim%grid); class is (t_grid_group)
    call grid % parallel_partition_io( sim%no_co )
  end select

  ! Initialize group I/O object
  call sim % grid % io % init( sim%grid, sim%no_co )

  ! initialization of the random number seed
  ! only do it once because it's a module variable
  ! TODO: double check all of this
  SCR_ROOT(" Setting up random number seed")
  call init_genrand( sim%options%random_class, sim%options%enforce_rng_constancy, &
                     sim%options%random_seed, sim%no_co%ngp_id() + sim%options%random_seed, &
                     sim%grid, restart, restart_handle  )

  ! determine cell size
  ldx = (xmax( sim%g_space ) - xmin( sim%g_space ))/ sim%grid%g_nx( 1:p_x_dim )

  SCR_ROOT(" Setting up group tstep")
  call setup( sim%tstep, -1, restart, restart_handle )

  SCR_ROOT(" Setting up group restart")
  call setup( sim%restart, restart, restart_handle )

  SCR_ROOT(" Setting up group g_space")
  call setup( sim%g_space, ldx, sim%grid%coordinates, restart, restart_handle )

  SCR_ROOT(" Setting up group time")
  call setup( sim%time, tmin(sim%time) - dt(sim%tstep), &
              dt(sim%tstep) , restart, restart_handle )

  SCR_ROOT(' Allocating simulation_tiles objects.')
  select type(nc_til => sim%no_co); class is(t_node_conf_tiles)
    sim % til_id = nc_til % ti_co % g_til_aid_min_max(:,nc_til%my_aid())
  end select
  ! Allocate tiles and read input into all tiles
  call sim % allocate_tiles()

  !---------------------------------------------------------------------------------------
  ! Set up tiles objects
  !---------------------------------------------------------------------------------------

  SCR_ROOT(" Setting up tiles no_co")
  do i = sim%til_id(1), sim%til_id(2)
    no_co => sim%tiles(i)%t%no_co; select type(no_co); class is(t_node_conf_tiles)
      call no_co % setup_tiles( sim%no_co )
    end select
  enddo

  ! setup the sim tiles grid
  SCR_ROOT(" Setting up tiles grid")
  do i = sim%til_id(1), sim%til_id(2)
    ! This is to make the min_gc_sim call work for some gcc compilers
    simulation => sim%tiles(i)%t
    call simulation % grid % init( simulation%no_co, &
      simulation%min_gc_sim(), x_bnd( sim%g_space ), restart, restart_handle )
  enddo

  ! Get partition information for all tiles
  SCR_ROOT(' - setting tile grid boundaries...')
  select type(group_grid => sim%grid); class is(t_grid_group)

    do i = sim%til_id(1), sim%til_id(2)
      grid => sim%tiles(i)%t%grid
      select type(grid); class is(t_grid_tiles)
        call grid % parallel_partition_tiles( sim%tiles(i)%t%no_co, group_grid%nx_til )
      end select
    enddo

  end select

  ! SCR_ROOT(" Setting up tstep")
  do i = sim%til_id(1), sim%til_id(2)
    sim%tiles(i)%t%tstep = sim%tstep
  enddo

  ! SCR_ROOT(" Setting up restart")
  do i = sim%til_id(1), sim%til_id(2)
    sim%tiles(i)%t%restart = sim%restart
  enddo

  ! We only have one g_space to avoid roundoff error in moving window

  ! SCR_ROOT(" Setting up time")
  do i = sim%til_id(1), sim%til_id(2)
    sim%tiles(i)%t%time = sim%time
  enddo

  SCR_ROOT(" Setting up boundary buffers")
  do i = sim%til_id(1), sim%til_id(2)
    call sim % tiles(i) % t % bnd % init( p_max_buffer1_size_t )
  enddo

  SCR_ROOT(" Setting up emf")

  min_gc = get_min_gc( sim%tiles(sim%til_id(1))%t%part )

  ! KYLE AND ROMAN MUFFIN LOBSTER
  ! #banter
  do i = sim%til_id(1), sim%til_id(2)
    call sim%tiles(i)%t%emf%init( get_grid_center( sim%tiles(i)%t%part ), &
                                  interpolation( sim%tiles(i)%t%part ), &
                                  sim%g_space, sim%tiles(i)%t%grid, &
                                  min_gc, ldx, sim%tiles(i)%t%tstep, &
                                  tmin(sim%tiles(i)%t%time), tmax(sim%tiles(i)%t%time), &
                                  sim%tiles(i)%t%no_co, &
                                  sim%tiles(i)%t%bnd%send_vdf, &
                                  sim%tiles(i)%t%bnd%recv_vdf, &
                                  restart, &
                                  restart_handle, sim%tiles(i)%t%options )
  enddo

  SCR_ROOT(" Setting up current")

  ! The field solver may require than particles ( but the cells for emf smoothing don't
  ! need to be taken into account )
  emf => sim%tiles(sim%til_id(1))%t%emf
  min_gc_tmp = emf % min_gc( )


  do i_dim = 1, p_x_dim
    if ( min_gc_tmp( 1, i_dim ) > min_gc( 1, i_dim ) ) min_gc( 1, i_dim ) = min_gc_tmp( 1, i_dim )
    if ( min_gc_tmp( 2, i_dim ) > min_gc( 2, i_dim ) ) min_gc( 2, i_dim ) = min_gc_tmp( 2, i_dim )
  enddo

  do i = sim%til_id(1), sim%til_id(2)
    emf => sim%tiles(i)%t%emf
    call sim%tiles(i)%t%jay%init( sim%tiles(i)%t%grid, min_gc, ldx, &
                sim%tiles(i)%t%no_co, if_move(sim%g_space), &
                emf%bc_type(), interpolation( sim%tiles(i)%t%part ), &
                restart, restart_handle, sim%tiles(i)%t%options )
  enddo
  ! Set interpolation for group diagnostic initialization
  sim%jay%interpolation = sim%tiles(sim%til_id(1))%t%jay%interpolation

  SCR_ROOT(" Setting up zpulses")
  do i = sim%til_id(1), sim%til_id(2)
    call sim%tiles(i)%t % zpulse_list % init( restart, t(sim%tiles(i)%t%time), &
              dt(sim%tiles(i)%t%tstep), sim%tiles(i)%t%emf, sim%g_space, &
              sim%tiles(i)%t%grid%my_nx(p_lower,:), sim%tiles(i)%t%no_co, &
              interpolation( sim%tiles(i)%t%part ) )
  enddo

  SCR_ROOT(" Setting up particles")
  do i = sim%til_id(1), sim%til_id(2)
    call sim%tiles(i)%t%part%init( sim%g_space, sim%tiles(i)%t%jay, sim%tiles(i)%t%emf, &
                sim%tiles(i)%t%grid, sim%tiles(i)%t%no_co, sim%tiles(i)%t%bnd, &
                sim%tiles(i)%t%zpulse_list, ndump(sim%tiles(i)%t%tstep), restart, &
                restart_handle, t(sim%tiles(i)%t%time), sim%tiles(i)%t%tstep, &
                tmin(sim%tiles(i)%t%time), tmax(sim%tiles(i)%t%time), &
                sim%tiles(i)%t%options  )
  enddo

  SCR_ROOT(" Setting up antenna array")
  do i = sim%til_id(1), sim%til_id(2)
    call sim%tiles(i)% t % antenna_array % init( restart, restart_handle )
  enddo


  SCR_ROOT(" Setting up group diagnostics utilities")
  ! Setup objects for diagnostics
  call sim % init_dutil( ldx )
  ! Get new diag msg pattern - sending/receiving tiles
  call sim % new_diag_patt()

  if ( restart ) then
    call restart_read_close( sim%restart, restart_handle )
  endif

  SCR_ROOT(" Setting up tile boundary information")
  call sim % new_bnd_patt()
  call sim % bnd_grp % init( sim%bnd_patt )

  ! Set n_current_* params here since the corresponding sim objects don't participate
  ! in restart
  sim % jay % n_current_jay = n(sim%tstep) + 1
  sim % emf % n_current_emf = n(sim%tstep) + 1
  species => sim%part%species
  do
    if (.not. associated(species)) exit
    species % n_current_spec = n(sim%tstep) + 1
    species => species % next
  enddo

  SCR_ROOT('')
  SCR_ROOT('Initialization complete!')

end subroutine init_sim_group
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Initialize new tiles array for dynamic load balance
!-----------------------------------------------------------------------------------------
subroutine init_dlb_group( sim, patt, tiles_temp, dif, diff, difff )

  use m_simulation_group
  use m_system
  use m_parameters
  use m_simulation, only : t_simulation, loopev, dynlbev, dynlb_set_int_load_ev, &
                           dynlb_new_lb_ev, dynlb_reshape_ev, dynlb_reshape_part_ev, &
                           dynlb_reshape_grid_ev, restart_write_ev, report_load_ev
  use m_simulation_tiles, only : t_simulation_tiles
  use m_input_file
  use m_node_conf, only : t_node_conf, my_ngp, comm, root, n_threads, no_num
  use m_node_conf_tiles
  use m_tile_conf, only : p_topology_viktor, g_tgp, g_tgp_old, set_ti_co_old
  use m_grid_define, only : t_grid, t_msg_patt
  use m_grid_parallel, only : if_dynamic_lb, clear_int_load, use_part_load
  use m_grid_tiles, only : t_grid_tiles, init_tile_load, add_load, conv_load_2_int_load, &
                           new_dlb_patt, t_grid_group
  use m_restart
  use m_random, only : init_genrand, list_algorithm_rnd
  use m_space, only : xmin_initial, xmax_initial, x_bnd, if_move, xmax, xmin, move_boundary
  use m_space, only : restart_write, setup
  use m_time_step, only : n, dt, ndump, setup, advance, restart_write
  use m_time, only : t, setup, tmin, update, restart_write, tmax
  use m_emf_define, only : t_emf, p_extfld_none
  use m_emf_diag, only : p_e1, p_e2, p_e3
  use m_particles_define, only : t_particles
  use m_particles, only : get_min_gc, get_grid_center, interpolation, &
                          ionize_neutral, sort_collide, ionize_neutral, sort_collide
  use m_particles_tiles, only : add_load
  use m_antenna_array, only : t_antenna_array, antenna_on
  use m_bnd_tiles, only : t_bnd_group, t_bnd_patt
  use m_diagnostic_utilities, only : cleanup_dutil
  use m_species_comm_tiles, only : p_max_buffer1_size_t
  use m_species_define, only : t_species
  use m_logprof
  use m_vdf_define, only : t_vdf_report
  use m_dynamic_loadbalance, only : sim_imbalance
  use m_random, only : rng, write_checkpoint_rnd

  implicit none

  class( t_simulation_group ), intent(inout)      :: sim
  type( t_msg_patt ), intent(inout)               :: patt
  class( t_til_container ), dimension(:), pointer :: tiles_temp
  integer, dimension(:), optional                 :: dif
  integer, dimension(:,:), optional               :: diff
  integer, dimension(:,:,:), optional             :: difff

  real( p_double ), dimension( p_x_dim ) :: ldx
  type( t_restart_handle )               :: restart_handle
  logical                                :: restart
  integer, dimension( 2, p_x_dim )       :: min_gc, min_gc_tmp
  integer                                :: i_dim
  integer                                :: til, til_old
  integer                                :: til_id1
  integer, dimension(p_x_dim)            :: lb, ub
  integer                                :: i, j, k
  real(p_double), dimension(2,p_x_dim)   :: x_bnd_initial
  class(t_node_conf), pointer            :: no_co, no_co_temp
  class(t_grid), pointer                 :: grid

  ! these pointers are added to fix a compiler error in GCC 7 wherein GCC could not
  !   correctly resolve long deference chains that ended in functions that return arrays.
  class(t_simulation), pointer           :: simulation
  class(t_emf), pointer                  :: emf

  ! Store sim%til_id(1) for later use initializing time/tstep properly
  til_id1 = sim % til_id(1)

  no_co => sim%no_co
  select type(no_co); class is(t_node_conf_tiles)

  ! Set new, post-dlb g_til_aid range
  sim % til_id = no_co%ti_co % g_til_aid_min_max(:,sim%no_co%my_aid())

  ! Allocate temporary tile array
  call sim % allocate_tiles_dlb( patt, tiles_temp )

  ! Copy data for stationary tiles (tiles that will remain on node after dlb)
  select case ( no_co%ti_co%topology_type )
    case (p_topology_viktor)
      select case( p_x_dim )
      case (2)
        ! dimension 1 lb and ub
        lb(1) = no_co%ti_co%cart_til_mm(1)%ii( 1, my_ngp(sim%no_co,1), 1, 1 )
        ub(1) = no_co%ti_co%cart_til_mm(1)%ii( 2, my_ngp(sim%no_co,1), 1, 1 )
        ! dimension 2 lb and ub
        lb(2) = no_co%ti_co%cart_til_mm(2)%ii( 1, my_ngp(sim%no_co,1), my_ngp(sim%no_co,2), 1 )
        ub(2) = no_co%ti_co%cart_til_mm(2)%ii( 2, my_ngp(sim%no_co,1), my_ngp(sim%no_co,2), 1 )
        ! get old tile data
        do i = lb(1), ub(1)
          do j = lb(2), ub(2)
            if ( diff(i,j)==0 ) then
              ! post load balance g_til_aid, used to index the new tile array
              til = no_co%ti_co%hi2( i, j )
              ! pre load balance g_til_aid, used to index the old tile array
              til_old = no_co%ti_co%hi2_old(i,j)
              ! put stationary tile in proper place in new tile array
              tiles_temp(til)%t => sim%tiles(til_old)%t
              ! g_til_aid's may change upon dlb with viktor's alg
              no_co_temp => tiles_temp(til)%t%no_co;
              select type(no_co_temp); class is(t_node_conf_tiles)
                no_co_temp%ti_co%g_til_aid = til
              end select
            endif
          enddo
        enddo
      case (3)
        ! dimension 1 lb and ub
        lb(1) = no_co%ti_co%cart_til_mm(1)%ii( 1, my_ngp(sim%no_co,1), 1, 1 )
        ub(1) = no_co%ti_co%cart_til_mm(1)%ii( 2, my_ngp(sim%no_co,1), 1, 1 )
        ! dimension 2 lb and ub
        lb(2) = no_co%ti_co%cart_til_mm(2)%ii( 1, my_ngp(sim%no_co,1), my_ngp(sim%no_co,2), 1 )
        ub(2) = no_co%ti_co%cart_til_mm(2)%ii( 2, my_ngp(sim%no_co,1), my_ngp(sim%no_co,2), 1 )
        ! dimension 3 lb and ub
        lb(3) = no_co%ti_co%cart_til_mm(3)%ii( 1, my_ngp(sim%no_co,1), my_ngp(sim%no_co,2), &
                                             my_ngp(sim%no_co,3) )
        ub(3) = no_co%ti_co%cart_til_mm(3)%ii( 2, my_ngp(sim%no_co,1), my_ngp(sim%no_co,2), &
                                             my_ngp(sim%no_co,3) )
        ! get old tile data
        do i = lb(1), ub(1)
          do j = lb(2), ub(2)
            do k = lb(3), ub(3)
              if ( difff(i,j,k)==0 ) then
                ! post load balance g_til_aid, used to index the new tile array
                til = no_co%ti_co%hi3( i, j, k )
                ! pre load balance g_til_aid, used to index the old tile array
                til_old = no_co%ti_co%hi3_old(i,j,k)
                ! put stationary tile in proper place in new tile array
                tiles_temp(til)%t => sim%tiles(til_old)%t
                ! g_til_aid's may change upon dlb with viktor's alg
                no_co_temp => tiles_temp(til)%t%no_co;
                select type(no_co_temp); class is(t_node_conf_tiles)
                  no_co_temp%ti_co%g_til_aid = til
                end select
              endif
            enddo
          enddo
        enddo
      end select
    case default
      ! for all other topologies it's much more straight forward
      ! any given tile's id does not change over the course of a load balance
      do i = sim%til_id(1), sim%til_id(2)
        if (dif(i)==0) then
          tiles_temp(i)%t => sim%tiles(i)%t
        endif
      enddo
  end select

  end select

  ! Setup incoming tiles as at time-step 0, see init_sim_group()

  ! never use restart for dlb initialization
  restart = .false.

  ! intialization of noco, need to do this for all tiles, not just incoming
  ! tiles, so that neighbors get updated for post-dlb partition
  do i = sim%til_id(1), sim%til_id(2)
    no_co => tiles_temp(i)%t%no_co; select type(no_co); class is(t_node_conf_tiles)
      call no_co % setup_tiles( sim%no_co )
    end select
  enddo

  do i = 1, patt%n_recv
    do j = 1, size( patt%recv(i)%tils )
      til = patt%recv(i)%tils(j)
      x_bnd_initial(p_lower,:) = xmin_initial( sim%g_space )
      x_bnd_initial(p_upper,:) = xmax_initial( sim%g_space )
      ! This is to make the min_gc_sim call work for some gcc compilers
      simulation => tiles_temp(til)%t
      call simulation % grid % init( simulation%no_co, &
                  simulation%min_gc_sim(), x_bnd_initial, restart, restart_handle )
    enddo
  enddo

  ! Get partition information for all tiles
  select type(group_grid => sim%grid); class is(t_grid_group)

    do i = 1, patt%n_recv
      do j = 1, size( patt%recv(i)%tils )

        til = patt%recv(i)%tils(j)
        grid => tiles_temp(til)%t%grid
        select type(grid); class is(t_grid_tiles)
          call grid % parallel_partition_tiles( tiles_temp(til)%t%no_co, &
                                                group_grid%nx_til )
        end select

      enddo
    enddo

  end select
  
  ! determine cell size
  ldx = ( xmax_initial( sim%g_space ) - xmin_initial( sim%g_space ) ) / &
        tiles_temp(sim%til_id(1))%t%grid%g_nx( 1:p_x_dim )

  do i = 1, patt%n_recv
    do j = 1, size( patt%recv(i)%tils )
      til = patt%recv(i)%tils(j)
      tiles_temp(til)%t%tstep = sim%tstep
    enddo
  enddo

  do i = 1, patt%n_recv
    do j = 1, size( patt%recv(i)%tils )
      til = patt%recv(i)%tils(j)
      tiles_temp(til)%t%restart = sim%restart
    enddo
  enddo

  ! We only have one g_space to avoid roundoff error in moving window

  do i = 1, patt%n_recv
    do j = 1, size( patt%recv(i)%tils )
      til = patt%recv(i)%tils(j)
      tiles_temp(til)%t%time = sim%time
    enddo
  enddo

  do i = 1, patt%n_recv
    do j = 1, size( patt%recv(i)%tils )
      til = patt%recv(i)%tils(j)
      call tiles_temp(til) % t % bnd % init( p_max_buffer1_size_t )
    enddo
  enddo

  min_gc = get_min_gc( tiles_temp(sim%til_id(1))%t%part )

  do i = 1, patt%n_recv
    do j = 1, size( patt%recv(i)%tils )
      til = patt%recv(i)%tils(j)
      ! TODO: We may have problems with set_fld_values and using sim%g_space

      ! KYLE AND ROMAN MUFFIN LOBSTER
      call tiles_temp(til)%t%emf%init( get_grid_center( tiles_temp(til)%t%part ), &
                                     interpolation( tiles_temp(til)%t%part ), &
                                     sim%g_space, &
                                     tiles_temp(til)%t%grid, min_gc, ldx, &
                                     tiles_temp(til)%t%tstep, &
                                     tmin(tiles_temp(til)%t%time), &
                                     tmax(tiles_temp(til)%t%time), &
                                     tiles_temp(til)%t%no_co, &
                                     tiles_temp(til)%t%bnd%send_vdf, &
                                     tiles_temp(til)%t%bnd%recv_vdf, &
                                     restart, &
                                     restart_handle, tiles_temp(til)%t%options )
    enddo
  enddo

  ! The field solver may require than particles ( but the cells for emf smoothing don't
  ! need to be taken into account )
  emf => tiles_temp(sim%til_id(1))%t%emf
  min_gc_tmp = emf % min_gc()

  do i_dim = 1, p_x_dim
    if ( min_gc_tmp(1, i_dim) > min_gc(1, i_dim) ) min_gc(1, i_dim) &
      = min_gc_tmp(1, i_dim)
    if ( min_gc_tmp(2, i_dim) > min_gc(2, i_dim) ) min_gc(2, i_dim) &
      = min_gc_tmp(2, i_dim)
  enddo

  do i = 1, patt%n_recv
    do j = 1, size( patt%recv(i)%tils )
      til = patt%recv(i)%tils(j)
      emf => tiles_temp(til)%t%emf
      call tiles_temp(til)%t%jay%init( tiles_temp(til)%t%grid, min_gc, ldx, &
                  tiles_temp(til)%t%no_co, if_move(sim%g_space), &
                  emf%bc_type(), interpolation( tiles_temp(til)%t%part ), &
                  restart, restart_handle, tiles_temp(til)%t%options )
    enddo
  enddo

  do i = 1, patt%n_recv
    do j = 1, size( patt%recv(i)%tils )
      til = patt%recv(i)%tils(j)
      ! don't launch if already past the time
      call tiles_temp(til)%t % zpulse_list % init( .true., t(tiles_temp(til)%t%time), &
              dt(tiles_temp(til)%t%tstep), tiles_temp(til)%t%emf, sim%g_space, &
              tiles_temp(til)%t%grid%my_nx(p_lower,:), tiles_temp(til)%t%no_co, &
              interpolation( tiles_temp(til)%t%part ) )
    enddo
  enddo

  do i = 1, patt%n_recv
    do j = 1, size( patt%recv(i)%tils )
      til = patt%recv(i)%tils(j)
      ! Don't initialize particles because they will be recieved over MPI
      call set_init_part( tiles_temp(til)%t%part, .false. )
      call tiles_temp(til)%t%part%init( sim%g_space, &
                                      tiles_temp(til)%t%jay, tiles_temp(til)%t%emf, &
                                      tiles_temp(til)%t%grid, tiles_temp(til)%t%no_co, &
                                      tiles_temp(til)%t%bnd, &
                                      tiles_temp(til)%t%zpulse_list, &
                                      ndump(tiles_temp(til)%t%tstep), restart, &
                                      restart_handle, t(tiles_temp(til)%t%time), &
                                      tiles_temp(til)%t%tstep, &
                                      tmin(tiles_temp(til)%t%time), &
                                      tmax(tiles_temp(til)%t%time), &
                                      tiles_temp(til)%t%options )
      ! Set this to true for proper moving window behavior later on
      call set_init_part( tiles_temp(til)%t%part, .true. )
    enddo
  enddo

  do i = 1, patt%n_recv
    do j = 1, size( patt%recv(i)%tils )
      til = patt%recv(i)%tils(j)
      call tiles_temp(til)%t%antenna_array % init( restart, restart_handle )
    enddo
  enddo

  ! Get new diag msg pattern - sending/receiving tiles
  ! Need to do this every dlb since tiles move, but all other diag setup remains const
  call sim % new_diag_patt()
  ! Do we need init_dutil? I don't think so.
  ! call sim % init_dutil( ldx )

  ! Set up tile boundary information after tiles_temp is now sim%tiles

end subroutine init_dlb_group
!-----------------------------------------------------------------------------------------

