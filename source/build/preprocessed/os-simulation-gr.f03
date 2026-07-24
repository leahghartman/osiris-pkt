# 1 "gr/os-simulation-gr.f03"
# 1 "<built-in>" 1
# 1 "<built-in>" 3
# 467 "<built-in>" 3
# 1 "<command line>" 1
# 1 "<built-in>" 2
# 1 "gr/os-simulation-gr.f03" 2
module m_simulation_gr

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
# 4 "gr/os-simulation-gr.f03" 2
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
# 5 "gr/os-simulation-gr.f03" 2

use m_simulation
use m_system
use m_parameters

use m_input_file, only : t_input_file
use m_restart, only : t_restart_handle
use m_geometry_gr, only : t_geometry_gr

! this is just meant to prevent re-exporting m_simulation symbols, everything is public
private

type, extends( t_simulation ) :: t_simulation_gr

  type( t_geometry_gr ), pointer :: geometry => null()

contains

  procedure :: iter => iter_sim_gr
  procedure :: allocate_objs => allocate_objs_gr
  procedure :: init => init_gr
  procedure :: read_input => read_input_gr
  procedure :: cleanup => cleanup_gr

  procedure :: write_checkpoint => write_checkpoint_gr

end type t_simulation_gr

public :: t_simulation_gr

contains

!-----------------------------------------------------------------------------------------
! Simulation iteration
!-----------------------------------------------------------------------------------------
subroutine iter_sim_gr( sim )

  use m_time_step, only : advance, dt, n
  use m_time, only : update, t, tmin
  use m_node_conf, only : root, n_threads
  use m_emf
  use m_particles
  use m_logprof
  use m_particles_charge
  use m_space
  use m_vdf_define
  use m_particles_gr, only : sort

  implicit none

  class ( t_simulation_gr ), intent(inout) :: sim

  ! set the initial timer for the loop
  call begin_event(loopev)

  ! advance time step counter and time

  call advance( sim%tstep, 1 )
  call update( sim%time, dt(sim%tstep), n(sim%tstep) )

  call sim % geometry % advance(t(sim%time))

  if ( root(sim%no_co) ) then
     print '(A,I10,A,F0.3)', 'n = ', n(sim%tstep), ', t = ', t(sim%time)
  endif

  ! Electric current diagnostics must occur before dynamic load
  ! balancing, since it destroys the current values

  call sim%jay%report(sim%g_space, sim%grid, sim%no_co, sim%tstep, t(sim%time) )

  ! update boundaries of em-fields and calculate
  ! time averaged fields if necessary

  call sim%emf%update_boundary( sim%g_space, sim%no_co, sim%bnd%send_vdf, sim%bnd%recv_vdf )

  ! Update the fields used for particle interpolation. This includes spatial smoothing,
  ! external fields and subcycling

  call sim%emf%update_particle_fld( n(sim%tstep), t(sim%time), dt(sim%tstep), sim%g_space, &
                            sim%grid%my_nx( p_lower, : ))

  ! diagnostic of em-field and particles
  ! currents can only be diagnosed after advance_deposit and are taken care of at the
  ! beggining of the loop

  if ( sim % emf % if_charge_cons( sim%tstep ) ) then
    call update_charge( sim%part, sim%g_space, sim%grid, sim%no_co, sim%bnd%send_vdf, sim%bnd%recv_vdf )
    write(err_buf__,*) "Charge conservation diagnostic is temporarily disabled";call err__("gr/os-simulation-gr.f03",93)
    call abort_program()
    ! call update_charge_cons( sim%emf, sim%part%charge%current, sim%bnd%send_vdf, sim%bnd%recv_vdf )
  endif

  call sim % emf % report( sim%g_space, sim%grid, sim%no_co, sim%tstep, t(sim%time), sim%bnd%send_vdf, sim%bnd%recv_vdf )


  call sim % part % report( sim%emf, sim%g_space, sim%grid, sim%no_co, &
                    sim%tstep, t(sim%time), tmin(sim%time), sim%bnd%send_vdf, sim%bnd%recv_vdf )

  ! move particles a timestep further and get new current

  call sim % part % advance_deposit( sim%emf, sim%jay, sim%tstep, t(sim%time), &
                                     sim%no_co, sim%options )

  ! take care of boundaries for particles

  call sim % part % update_boundary( sim%jay, sim%g_space, sim%no_co, dt(sim%tstep), sim%bnd )

  ! report time centered energy
  call sim % emf % report_energy( sim%no_co, sim%tstep, t(sim%time), tmin(sim%time))
  call sim % part % report_energy( sim%no_co, sim%tstep, t(sim%time), tmin(sim%time))

  ! particle sorting
  call sort( sim%part, n(sim%tstep), t(sim%time) )

  ! take care of boundaries for currents

  call sim%jay%update_boundary( sim%no_co, sim%bnd%send_vdf, sim%bnd%recv_vdf )

  ! smooth currents
  ! - smoothing has to be after update boundary

  call sim % jay % smooth()

  ! calculate new electro-magnetic fields

  call sim % emf % advance( sim%jay, dt(sim%tstep), sim%bnd%send_vdf, sim%bnd%recv_vdf )

  ! set the final time for the loop
  call end_event(loopev)

end subroutine iter_sim_gr

!-----------------------------------------------------------------------------------------
! Allocate algorithm specific simulation objects
!-----------------------------------------------------------------------------------------
subroutine allocate_objs_gr( sim )

  use m_geometry_gr, only : t_geometry_gr
  use m_emf_define_gr, only : t_emf_gr
  use m_particles_define_gr, only : t_particles_gr

  implicit none

  class( t_simulation_gr ), intent(inout) :: sim

  if (mpi_node()==0) print *,'Allocating gr geometry, emf and particle objects.'

  allocate( t_geometry_gr :: sim%geometry )

  allocate( t_emf_gr :: sim%emf )

  allocate( t_particles_gr :: sim%part )

  ! allocate any remaining objects in superclass
  call sim % t_simulation % allocate_objs()

end subroutine allocate_objs_gr

!-----------------------------------------------------------------------------------------
! Initialize simulation objects
!-----------------------------------------------------------------------------------------
subroutine init_gr( sim, grid_gc_min )

  use m_restart, only : t_restart_handle, if_restart_read, &
                        restart_read_open, restart_read_close
  use m_node_conf, only : comm
  use m_space, only : x_bnd, xmin, xmax, if_move
  use m_time, only : setup, tmin, tmax, t
  use m_time_step, only : dt, setup, ndump
  use m_emf_define, only : t_emf
  use m_emf_define_gr, only: t_emf_gr
  use m_particles_define_gr, only : t_particles_gr

  use m_grid_parallel
  use m_random
  use m_logprof
  use m_particles
  use m_emf
  use m_diagnostic_utilities, only : init_dutil

  implicit none

  ! dummy variables
  class ( t_simulation_gr ), intent(inout) :: sim

  real( p_double ), dimension( p_x_dim ) :: ldx

  integer :: i_dim
  integer, dimension( 2, p_x_dim ) :: min_gc, min_gc_tmp
  integer, dimension( 2, p_x_dim ), intent(in) :: grid_gc_min


  logical :: restart
  type( t_restart_handle ) :: restart_handle
  integer :: diag_buffer_size

  ! initialize timing structure
  loopev = create_event('loop')
  dynlbev = create_event('dynamic load balance (total)')
  dynlb_set_int_load_ev = create_event('dynlb get int load')
  dynlb_new_lb_ev = create_event('dynlb get new lb')
  dynlb_reshape_ev = create_event('dynlb reshape simulation')
  dynlb_reshape_part_ev = create_event('dynlb reshape simulation (particles)')
  dynlb_reshape_grid_ev = create_event('dynlb reshape simulation (grids)')
  restart_write_ev = create_event('writing restart information')
  report_load_ev = create_event('load diagnostics')

  ! use restart files if restart is set in the input file
  restart = if_restart_read(sim%restart)

  if (mpi_node()==0) print *,''
  if ( restart ) then
    if (mpi_node()==0) print *,'**********************************************'
    if (mpi_node()==0) print *,'* Restarting simulation from checkpoint data *'
    if (mpi_node()==0) print *,'**********************************************'
  else
    if (mpi_node()==0) print *,'Initializing simulation:'
  endif

  if (mpi_node()==0) print *," Setting up no_co"
  call sim % no_co % init()

  if ( restart ) then
    ! open restart file (this needs to happen after setup no_co)
    call restart_read_open( sim%restart, comm(sim%no_co), sim%no_co%ngp_id(), file_id_rst, restart_handle )

    ! read restart information for no_co
    ! This only checks for inconsistencies in the input file vs. restart files, all setup
    ! is completed above
    call sim % no_co % restart_read( restart_handle )
  endif

  ! setup the simulation grid
  if (mpi_node()==0) print *," Setting up grid"

  call sim%grid%init( sim%no_co, grid_gc_min, x_bnd( sim%g_space ), &
              restart, restart_handle )

  ! if not restarting generate initial load balance partition
  if ( .not. restart ) then

     if ( sim % grid % needs_int_load( 0 ) ) then
       if (mpi_node()==0) print *,' - generating initial integral load...'
       ! create even partition for calculating particle (and grid) loads
       call sim % grid % parallel_partition( sim%no_co, -2 )
       call sim % set_int_load( sim%grid, 0 )
     endif

     if (mpi_node()==0) print *,' - setting parallel grid boundaries...'
     call sim % grid % parallel_partition( sim%no_co, 0 )
     call clear_int_load( sim%grid )

  endif

  ! Initialize I/O object
  call sim % grid % io % init( sim%grid, sim%no_co )

  ! initialization of the random number seed
  if (mpi_node()==0) print *," Setting up random number seed"
  call init_genrand( sim%options%random_class, sim%options%enforce_rng_constancy, &
                     sim%options%random_seed, sim%no_co%ngp_id() + sim%options%random_seed, &
                     sim%grid, restart, restart_handle )

  ! determine cell size
  ldx = (xmax( sim%g_space ) - xmin( sim%g_space ))/ sim%grid%g_nx( 1:p_x_dim )

  if (mpi_node()==0) print *," Setting up tstep"
  call setup( sim%tstep, -1 , restart, restart_handle )

  if (mpi_node()==0) print *," Setting up restart"
  call setup( sim%restart, restart, restart_handle )

  if (mpi_node()==0) print *," Setting up g_space"
  call setup( sim%g_space, ldx, sim%grid%coordinates , restart, restart_handle )

  if (mpi_node()==0) print *," Setting up time"
  call setup( sim%time, tmin(sim%time) - dt(sim%tstep), &
              dt(sim%tstep) , restart, restart_handle )

  if (mpi_node()==0) print *," Setting up geometry"
  call sim % geometry % init(sim%g_space, sim%grid)

  if (mpi_node()==0) print *," Setting up boundary buffers"
  call sim % bnd % init()

  if (mpi_node()==0) print *," Setting up emf"
  min_gc = get_min_gc( sim%part )

  select type( emf => sim%emf )
    class is ( t_emf_gr )
      emf%geometry => sim%geometry
    class default
      write(err_buf__,*) 'pointer to geometry must be set with a t_emf_gr object';call err__("gr/os-simulation-gr.f03",298)
      call abort_program( p_err_invalid )
  end select

  call sim%emf%init( get_grid_center( sim%part ), interpolation( sim%part ), &
                sim%g_space, sim%grid, min_gc, ldx, sim%tstep, tmin(sim%time), &
                tmax(sim%time), sim%no_co, sim%bnd%send_vdf, sim%bnd%recv_vdf, &
                restart, restart_handle, sim%options )

  if (mpi_node()==0) print *," Setting up current"
  ! The field solver may require than particles ( but the cells for emf smoothing don't
  ! need to be taken into account )
  min_gc_tmp = sim%emf%min_gc()

  do i_dim = 1, p_x_dim
    if ( min_gc_tmp( 1, i_dim ) > min_gc( 1, i_dim ) ) min_gc( 1, i_dim ) = min_gc_tmp( 1, i_dim )
    if ( min_gc_tmp( 2, i_dim ) > min_gc( 2, i_dim ) ) min_gc( 2, i_dim ) = min_gc_tmp( 2, i_dim )
  enddo

  call sim%jay%init( sim%grid, min_gc, ldx, &
              sim%no_co, if_move(sim%g_space), &
              sim % emf % bc_type(), interpolation( sim%part ), &
              restart, restart_handle, sim%options )

  if (mpi_node()==0) print *," Setting up particles"

  select type( part => sim%part )
    class is ( t_particles_gr )
    part%geometry => sim%geometry
    class default
      write(err_buf__,*) 'pointer to geometry must be set with a t_particles_gr object';call err__("gr/os-simulation-gr.f03",328)
      call abort_program( p_err_invalid )
  end select

  call sim%part%init( sim%g_space, sim%jay, sim%emf, &
              sim%grid, sim%no_co, sim%bnd, sim%zpulse_list, ndump(sim%tstep), &
              restart, restart_handle, t(sim%time), sim%tstep, tmin(sim%time), &
              tmax(sim%time), sim%options)

  if (mpi_node()==0) print *," Setting up diagnostics utilities"
  diag_buffer_size = 0
  call sim % part % get_diag_buffer_size( sim%grid%g_nx, diag_buffer_size )
  call sim % emf % get_diag_buffer_size( sim%grid%g_nx, diag_buffer_size )
  call sim % jay % get_diag_buffer_size( sim%grid%g_nx, diag_buffer_size )
  call init_dutil( sim%options, diag_buffer_size )

  if ( restart ) then
    call restart_read_close( sim%restart, restart_handle )
  endif

  if (mpi_node()==0) print *,''
  if (mpi_node()==0) print *,'Initialization complete!'

end subroutine init_gr

!-----------------------------------------------------------------------------------------
! Initialize simulation objects
!-----------------------------------------------------------------------------------------
subroutine read_input_gr( sim, input_file )

  use m_node_conf, only: periodic
  use m_time_step, only: read_nml, dt
  use m_space, only: read_nml, if_move, xmin, xmax
  use m_time, only: read_nml

  implicit none

  class( t_simulation_gr ), intent(inout) :: sim
  class( t_input_file ), intent(inout) :: input_file

  real( p_double ), dimension( p_x_dim ) :: ldx

  if (mpi_node()==0) print *,'Reading parallel node configuration... '
  call sim % no_co % read_input( input_file, p_x_dim )

  if (mpi_node()==0) print *,'Reading grid configuration...'
  call sim % grid % read_input( input_file, sim % no_co % nx )

  if (mpi_node()==0) print *,'Reading tstep configuration... '
  call read_nml( sim%tstep, input_file )

  if (mpi_node()==0) print *,'Reading restart configuration... '
  call read_nml( sim%restart, input_file, sim%options%restart )

  if (mpi_node()==0) print *,'Reading g_space configuration... '
  call read_nml( sim%g_space, input_file, periodic(sim%no_co), sim%grid%coordinates )

  if (mpi_node()==0) print *,'Reading time configuration... '
  call read_nml( sim%time, input_file )

  if (mpi_node()==0) print *,'Reading geometry configuration...'
  call sim % geometry % read_input( input_file )

  if (mpi_node()==0) print *,'Reading emf configuration... '
  ! determine cell size
  ldx = (xmax( sim%g_space ) - xmin( sim%g_space ))/ sim%grid%g_nx( 1:p_x_dim )
  call sim % emf % read_input( input_file, periodic(sim%no_co), if_move(sim%g_space), &
                                sim%grid, ldx, dt(sim%tstep), sim%options%gamma )

  if (mpi_node()==0) print *,'Reading part configuration... '
  call sim%part%read_input( input_file, periodic(sim%no_co), &
                  if_move(sim%g_space), sim%grid, dt(sim%tstep), &
                  sim%options )

  if (mpi_node()==0) print *,'Reading current configuration...'
  call sim % jay % read_input( input_file )

end subroutine read_input_gr

!-----------------------------------------------------------------------------------------
! Cleanup simulation objects
!-----------------------------------------------------------------------------------------
subroutine cleanup_gr( sim )

  use m_vdf_comm
  use m_wall_comm
  use m_particles
  use m_diagnostic_utilities

  implicit none

  class( t_simulation_gr ), intent(inout) :: sim

  call sim % grid % cleanup()

  ! cleanup particles
  call sim % part % cleanup( )

  ! cleanup current object
  call sim%jay%cleanup()

  ! cleanup emf object
  call sim%emf%cleanup()

  ! cleanup vdf and spec comm buffers
  call sim % bnd % cleanup()

  ! cleanup diagnostic utilities
  call cleanup_dutil()

  ! Cleanup MPI system
  call sim % no_co % cleanup()

  ! cleanup geometry
  call sim % geometry % cleanup()

  ! deallocate polymorphic objects
  deallocate( sim%input_file )
  deallocate( sim%no_co )
  deallocate( sim%grid )
  deallocate( sim%emf )
  deallocate( sim%part )
  deallocate( sim%jay )
  deallocate( sim%bnd )
  deallocate( sim%geometry )

end subroutine cleanup_gr

!---------------------------------------------------------------------------------------
! Write checkpoint data
!---------------------------------------------------------------------------------------
subroutine write_checkpoint_gr( sim, restart_handle )

  implicit none

  class ( t_simulation_gr ), intent(in) :: sim
  type ( t_restart_handle ), intent(inout) :: restart_handle

  ! write superclass checkpoint data first
  call sim % t_simulation % write_checkpoint( restart_handle )

  ! for now its placed here
  call sim % geometry % write_checkpoint( restart_handle )

  !do nothing for now

end subroutine write_checkpoint_gr

end module m_simulation_gr
