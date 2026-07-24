# 1 "cyl_modes/os-neutral-cyl-modes.f03"
# 1 "<built-in>" 1
# 1 "<built-in>" 3
# 467 "<built-in>" 3
# 1 "<command line>" 1
# 1 "<built-in>" 2
# 1 "cyl_modes/os-neutral-cyl-modes.f03" 2
! - BSI and BSI random ionization models
! - Vacuum ionization
! - Impact ionization
!
!-----------------------------------------------------------------------------------------


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
# 9 "cyl_modes/os-neutral-cyl-modes.f03" 2
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
# 10 "cyl_modes/os-neutral-cyl-modes.f03" 2



module m_neutral_cyl_modes

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
# 16 "cyl_modes/os-neutral-cyl-modes.f03" 2

use m_parameters

use m_restart, only: t_restart_handle, restart_io_write, restart_io_read
use m_vdf, only: move_window
use m_vdf_define, only: t_vdf, t_vdf_report
use m_vdf_comm, only: update_boundary, t_vdf_msg, reshape_copy, reshape_nocopy
use m_vdf_report, only: new
use m_vdf_memory
use m_node_conf, only: t_node_conf
use m_grid_define, only: t_grid
use m_space, only: t_space, nx_move
use m_species_cyl_modes_define, only: t_species_cyl_modes
use m_emf_cyl_modes, only: t_emf_cyl_modes
use m_species_define, only: t_species, p_cell_low, p_cell_near
use m_emf_define, only: t_emf
use m_neutral, only: t_neutral
use m_diag_neutral, only: setup, cleanup, report
use m_random, only: rng
use m_cross, only: setup, cleanup, restart_read, restart_write
use m_time_step, only: t_time_step
use m_math, only: pi

implicit none

private

! string to id restart data
character(len=*), parameter :: p_neutral_cyl_modes_rst_id = "neutral_cyl_modes rst data - 0x0002"

type, extends( t_neutral ) :: t_neutral_cyl_modes

  type( t_vdf ), dimension(:), pointer :: w_cyl => null()

  ! Ionization level previous
  type( t_vdf ), dimension(:), pointer :: ion_level_old_cyl => null()

  ! Multi-level densities
  type( t_vdf ), dimension(:), pointer :: multi_ion_cyl => null()

  ! small cosine tables for the calculating ionization rates (mode,x1,x2,theta)
  real(p_k_fld), dimension(:,:,:,:), pointer :: cosine => null(), sine => null()

  ! stores random theta offset for each cell
  type( t_vdf ) :: rand_theta

  !number particles in circle
  integer :: num_par_theta

contains

  procedure :: adk_field_rates => adk_field_rates_neutral_cyl_modes
  procedure :: init => init_neutral_cyl_modes
  procedure :: cleanup => cleanup_neutral_cyl_modes
  procedure :: advance => advance_density_cyl_modes
  procedure :: inject => inject_particles_cyl_modes
  procedure :: ionize => ionize_neutral_cyl_modes
  procedure :: reshape_obj => reshape_neutral_cyl_modes
  procedure :: add_particle_load => add_particle_load_neutral_cyl_modes
  procedure :: move_window => move_window_neutral_cyl_modes
  procedure :: update_boundary => update_boundary_neutral_cyl_modes
  procedure :: report => report_neutral_cyl_modes
  procedure :: restart_write => restart_write_neutral_cyl_modes
  procedure :: restart_read => restart_read_neutral_cyl_modes

end type t_neutral_cyl_modes

public :: t_neutral_cyl_modes

!-----------------------------------------------------------------------------------------
contains

!-----------------------------------------------------------------------------------------
subroutine init_neutral_cyl_modes( this, neutral_id, emf, nx_p_min, g_space, &
                          restart, restart_handle, sim_options )
!-----------------------------------------------------------------------------------------

  implicit none

  class( t_neutral_cyl_modes ), intent(inout) :: this
  integer, intent(in) :: neutral_id
  class( t_emf),intent(in) :: emf
  integer, dimension(:), intent(in) :: nx_p_min
  type( t_space ), intent(in) :: g_space
  logical, intent(in) :: restart
  type( t_restart_handle ), intent(in) :: restart_handle
  type( t_options ), intent(in) :: sim_options


  integer, dimension(p_x_dim) :: nlbound, nubound
  integer, dimension(2,p_x_dim) :: lgc_num
  integer :: num_par_theta, num_modes, mode, theta, i, j
  real(p_double) :: actual_theta
  logical :: rand_theta

  if ( restart ) then

    call this % restart_read( restart_handle )

  else

    ! set the internal neutral id (index in neutral array of particle object)
    this%neutral_id = neutral_id

    !set up cross section profile
    if (this%if_impact) then
      call setup( this%cross_section)
    endif

    ! check if associated species are ok
    if (.not. associated(this%species1)) then
      print *, "(*error*) species1 is not associated in setup_neutral, in os-neutral-cyl-modes.f03"
      print *, "(*error*) unable to setup neutral"
      print *, "(*error*) bailing out..."
      call abort_program()
    endif

    if (.not. associated(this%species2) .and. this%if_mov_ions) then
      print *, "(*error*) species2 is not associated in setup_neutral, in os-neutral-cyl-modes.f03"
      print *, "(*error*) unable to setup neutral"
      print *, "(*error*) bailing out..."
      call abort_program()
    endif

    select type( species1 => this%species1 )
    class is ( t_species_cyl_modes )
      num_par_theta = species1%num_par_theta
      rand_theta = species1%rand_theta
    end select

    ! Get random offset in theta for each cell
    ! We set guard cells to one so that the dlb communication doesn't error out
    lgc_num = 1
    call this % rand_theta % new( p_x_dim, 1, emf%e%nx_(1:p_x_dim), lgc_num, emf%e%dx(), .true. )
    if (rand_theta) then
      do j = 1, emf%e%nx_(2)
        do i = 1, emf%e%nx_(1)
          call rng % harvest_real2( this%rand_theta%f2(1,i,j), (/ i, j /) )
        enddo
      enddo
    endif

    call alloc(this%ion_level_old_cyl, (/ num_par_theta /),"cyl_modes/os-neutral-cyl-modes.f03",158)
    call alloc(this%w_cyl, (/ num_par_theta /),"cyl_modes/os-neutral-cyl-modes.f03",159)
    call alloc(this%multi_ion_cyl, (/ num_par_theta /),"cyl_modes/os-neutral-cyl-modes.f03",160)

    this%num_par_theta = num_par_theta

    do theta=1,num_par_theta

      ! Create ionization level and density vdf structures, based in vdf e
      ! copy = .false. causes vdf to be initialized with 0.
      call this % ion_level_old_cyl(theta) % new(emf%e, copy = .false., f_dim = 1)

      ! vdf with ion densities, neutral density (1st component) and
      ! total ion density (last component)
      ! copy = .false. causes vdf to be initialized with 0.
      call this % multi_ion_cyl(theta) % new( emf%e, copy = .false., f_dim = this%ion_idx)

      call this % w_cyl(theta) % new( emf%e , copy = .false. , f_dim = this%ion_idx)

      ! Fill neutral density (1st component of multi_ion) according to given neutral profile
      ! Since value at position uses guard cells we need to initialize them as well
      lgc_num = this%multi_ion_cyl(theta) % gc_num()
      nlbound = 1 - lgc_num( p_lower, 1:p_x_dim )
      nubound = emf%e % nx_(1:p_x_dim) + lgc_num( p_upper, 1:p_x_dim )
      call this % t_neutral % set_neutden_values( this%multi_ion_cyl(theta), this%den_neutral, this%neut_idx, this%den_min, &
        nx_p_min, g_space, nlbound, nubound)

    enddo

  endif

  ! set values not present in restart info
  this%omega_p = sim_options%omega_p0

  select type( emf_cyl => emf )
  class is ( t_emf_cyl_modes )
    num_modes = emf_cyl%n_cyl_modes
  end select

  ! We need at least the regular sine and cosine values
  if ( num_modes == 0 ) then
    num_modes = 1
  endif

  call alloc(this%cosine, (/ num_modes, emf%e%nx_(1), emf%e%nx_(2), this%num_par_theta /),"cyl_modes/os-neutral-cyl-modes.f03",202)
  call alloc(this%sine, (/ num_modes, emf%e%nx_(1), emf%e%nx_(2), this%num_par_theta /),"cyl_modes/os-neutral-cyl-modes.f03",203)

  do theta = 1, this%num_par_theta
    do j = 1, emf%e%nx_(2)
      do i = 1, emf%e%nx_(1)
        actual_theta = ( real(theta-1,p_k_fld)/real(this%num_par_theta,p_k_fld) &
                         + this%rand_theta%f2(1,i,j) ) * 2.0_p_k_fld * pi
        do mode = 1, num_modes
          this%cosine(mode,i,j,theta) = cos( mode*actual_theta )
          this% sine(mode,i,j,theta) = sin( mode*actual_theta )
        enddo
      enddo
    enddo
  enddo

  call setup( this%diag, this%name, this%species1%interpolation )

  ! warn about charge normalization error near the axis if not injecting in the middle of
  ! a cell in r
  if ( .not. this%inject_line ) then
    if ( mpi_node() == 0 ) then
      write(0,"(A)") " (*WARNING*) inject_line has been set to false for neutral."
      write(0,"(A)") " (*WARNING*) This will cause a small error in the charge density"
      write(0,"(A)") " (*WARNING*) near the axis. Set to true to inject in middle of cell"
      write(0,"(A)") " (*WARNING*) and eliminate this error."
    endif
  endif

end subroutine init_neutral_cyl_modes
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
subroutine cleanup_neutral_cyl_modes( this )
!-----------------------------------------------------------------------------------------

  implicit none

  class( t_neutral_cyl_modes ), intent(inout) :: this
  integer :: i

  call this % t_neutral % cleanup()

  do i=1,this%num_par_theta
    call this%ion_level_old_cyl(i) % cleanup()
    call this%multi_ion_cyl(i) % cleanup()
    call this%w_cyl(i) % cleanup()
  enddo

  call freemem(this%ion_level_old_cyl,"cyl_modes/os-neutral-cyl-modes.f03",252)
  call freemem(this%multi_ion_cyl,"cyl_modes/os-neutral-cyl-modes.f03",253)
  call freemem(this%w_cyl,"cyl_modes/os-neutral-cyl-modes.f03",254)

  call freemem(this%cosine,"cyl_modes/os-neutral-cyl-modes.f03",256)
  call freemem(this%sine,"cyl_modes/os-neutral-cyl-modes.f03",257)

  call this % rand_theta % cleanup()

end subroutine cleanup_neutral_cyl_modes
!-----------------------------------------------------------------------------------------



!-----------------------------------------------------------------------------------------
! report on ionization level diagnostic
!-----------------------------------------------------------------------------------------
subroutine report_neutral_cyl_modes( this, g_space, grid, no_co, tstep, t )
!-----------------------------------------------------------------------------------------

  implicit none

  class( t_neutral_cyl_modes ), intent(inout) :: this
  type( t_space ), intent(in) :: g_space
  class( t_grid ), intent(in) :: grid
  class( t_node_conf ), intent(in) :: no_co
  type( t_time_step ), intent(in) :: tstep
  real(p_double), intent(in) :: t

  call report(this%diag, this%multi_ion_cyl(1), this%ion_idx, this%neut_idx, g_space, &
              grid, no_co, tstep, t)
end subroutine report_neutral_cyl_modes
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
! Update boundaries of neutral object (communication only)
!-----------------------------------------------------------------------------------------
subroutine update_boundary_neutral_cyl_modes( this, nx_move, no_co, send_msg, recv_msg )

  implicit none

  class( t_neutral_cyl_modes ), intent( inout ) :: this
  integer :: theta
  integer, dimension(:), intent(in) :: nx_move
  class( t_node_conf ), intent(in) :: no_co
  type(t_vdf_msg), dimension(2), intent(inout) :: send_msg, recv_msg

  ! update boundaries ion_level vdf
  do theta =1, this%num_par_theta
    call update_boundary(this%multi_ion_cyl(theta), p_vdf_replace, no_co, send_msg, recv_msg, nx_move)
  enddo

end subroutine update_boundary_neutral_cyl_modes
!-----------------------------------------------------------------------------------------



!-----------------------------------------------------------------------------------------
! Move simulation window for the neutral object
!-----------------------------------------------------------------------------------------
subroutine move_window_neutral_cyl_modes( this, nx_p_min, g_space , need_den_val )

  implicit none

  class( t_neutral_cyl_modes ), intent(inout) :: this
  integer, dimension(:), intent(in) :: nx_p_min
  type( t_space ), intent(in) :: g_space
  logical , intent(in) :: need_den_val

  integer, dimension( 2, p_x_dim ) :: init_bnd
  integer :: i, num_par_theta, theta

  ! initalize ionization values where required
  ! note that we do not use values from another node, we simply
  ! recalculate them (which ends up being faster that having the
  ! extra communication)
  num_par_theta = this%num_par_theta
  do theta =1,num_par_theta

    ! move ion_level vdf
    call move_window( this%multi_ion_cyl(theta), g_space )
    call move_window( this%ion_level_old_cyl(theta), g_space )
    do i = 1, p_x_dim

      if ( nx_move( g_space, i ) > 0) then


        if(need_den_val) then

          init_bnd(p_lower,1:p_x_dim) = &
              1 - this%multi_ion_cyl(theta)%gc_num_( p_lower, 1:p_x_dim )


          init_bnd(p_upper,1:p_x_dim) = this%multi_ion_cyl(theta)%nx_(1:p_x_dim ) + &
          this%multi_ion_cyl(theta)%gc_num_( p_upper, 1:p_x_dim )

          init_bnd(p_lower,i) = this%multi_ion_cyl(theta)%nx_(i) - nx_move( g_space, i ) + this%multi_ion_cyl(theta)%gc_num_( p_lower, i )

          call this % t_neutral % set_neutden_values( this%multi_ion_cyl(theta), &
                                this%den_neutral, this%neut_idx, this%den_min, nx_p_min, &
                                g_space, init_bnd(p_lower,:), init_bnd(p_upper,:) )
        endif
      endif
    enddo
  enddo

end subroutine move_window_neutral_cyl_modes
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Reshape neutral object for dynamic load balance
!-----------------------------------------------------------------------------------------
subroutine reshape_neutral_cyl_modes( this, old_lb, new_lb, no_co, send_msg, recv_msg )
!-----------------------------------------------------------------------------------------

  implicit none

  class(t_neutral_cyl_modes), intent(inout) :: this
  class(t_grid), intent(in) :: old_lb, new_lb
  class( t_node_conf ), intent(in) :: no_co
  type(t_vdf_msg), dimension(2), intent(inout) :: send_msg, recv_msg

  integer :: theta, num_par_theta, num_modes, mode, i, j
  real(p_double) :: actual_theta
  type(t_vdf_report), pointer :: report

  num_par_theta = this%num_par_theta
  do theta = 1, num_par_theta

    call reshape_copy( this%ion_level_old_cyl(theta), old_lb, new_lb, no_co, send_msg, recv_msg )
    call reshape_copy( this%multi_ion_cyl(theta), old_lb, new_lb, no_co, send_msg, recv_msg )
    call reshape_nocopy( this%w_cyl(theta), new_lb )

  enddo

  ! Reshape other cyl_modes-specific data
  call reshape_copy( this%rand_theta, old_lb, new_lb, no_co, send_msg, recv_msg )

  num_modes = size( this%cosine, 1 )
  call freemem(this%cosine,"cyl_modes/os-neutral-cyl-modes.f03",392)
  call freemem(this%sine,"cyl_modes/os-neutral-cyl-modes.f03",393)
  call alloc(this%cosine, (/ num_modes, this%w_cyl(1)%nx_(1), this%w_cyl(1)%nx_(2), this%num_par_theta /),"cyl_modes/os-neutral-cyl-modes.f03",394)
  call alloc(this%sine, (/ num_modes, this%w_cyl(1)%nx_(1), this%w_cyl(1)%nx_(2), this%num_par_theta /),"cyl_modes/os-neutral-cyl-modes.f03",395)

  do theta = 1, num_par_theta
    do j = 1, this%w_cyl(1)%nx_(2)
      do i = 1, this%w_cyl(1)%nx_(1)
        actual_theta = ( real(theta-1,p_k_fld)/real(this%num_par_theta,p_k_fld) &
                         + this%rand_theta%f2(1,i,j) ) * 2.0_p_k_fld * pi
        do mode = 1, num_modes
          this%cosine(mode,i,j,theta) = cos( mode*actual_theta )
          this% sine(mode,i,j,theta) = sin( mode*actual_theta )
        enddo
      enddo
    enddo
  enddo

  ! Reshape tavg_data
  report => this%diag%reports
  do
    if ( .not. associated(report) ) exit
    if ( report%tavg_data%x_dim_ > 0 ) then
      call reshape_copy( report%tavg_data, old_lb, new_lb, no_co, send_msg, recv_msg )
    endif
    report => report%next
  enddo

end subroutine reshape_neutral_cyl_modes
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Add ionization calculation load
!-----------------------------------------------------------------------------------------
subroutine add_particle_load_neutral_cyl_modes( this, grid )
!-----------------------------------------------------------------------------------------

  implicit none

  ! dummy vars
  class(t_neutral_cyl_modes), intent(in) :: this
  class(t_grid), intent(inout) :: grid

  integer, dimension( p_max_dim ) :: dim_off
  real( p_double ), dimension(:), pointer :: int_load
  integer :: i, j, i1, i2

  ! base_load is the basic load for calculating adk rates, doing injection, etc.
  ! multi_fac is the additional load for each ionization level (mostly in advance_density)
  ! note that the real load will depend on the nearby field strength, but that is ignored
  ! these parameters were calculated based on a few sample simulations by K. Miller
  real(p_double), parameter :: base_load = 0.35, multi_fac = 0.15
  real(p_double) :: load

  ! Get the offset of the int_load array for each direction
  dim_off(1) = 0
  do i = 2, grid%x_dim
    dim_off(i) = dim_off(i-1) + grid%g_nx(i-1)
  enddo

  ! Add the correction for global cell positions to dim_off
  do i = 1, grid%x_dim
    dim_off(i) = dim_off(i) + ( this%species1%my_nx_p( p_lower, i ) - 1 )
  enddo

  ! just for clarity
  int_load => grid%int_load

  ! Approximate load for each cell
  load = ( base_load + this%multi_max * multi_fac ) * this%num_par_theta

  do j = 1, this%multi_ion_cyl(1)%nx_(2)
    do i = 1, this%multi_ion_cyl(1)%nx_(1)
      if ( this%multi_ion_cyl(1)%f2(this%neut_idx,i,j) > this%den_min ) then
        i1 = dim_off(1) + i
        i2 = dim_off(2) + j
        int_load(i1) = int_load(i1) + load
        int_load(i2) = int_load(i2) + load
      endif
    enddo
  enddo

end subroutine add_particle_load_neutral_cyl_modes
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
! Advance ionization level density
!-----------------------------------------------------------------------------------------
subroutine advance_density_cyl_modes(this , dt)
!-----------------------------------------------------------------------------------------

  class(t_neutral_cyl_modes), intent(inout) :: this ! neutral array
  real( p_k_fld ), intent(in) :: dt

  if(.not. this%if_analytic) then
    call advance_density_cyl_modes_std(this, dt)
  else
    call advance_density_cyl_modes_analytic(this, dt)
  endif

end subroutine advance_density_cyl_modes
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
subroutine advance_density_cyl_modes_std(this , dt)
!-----------------------------------------------------------------------------------------

  implicit none

  class(t_neutral_cyl_modes), intent(inout) :: this ! neutral array
  real( p_k_fld ), intent(in) :: dt
  integer :: i,j, m, ion_idx,theta
  real(p_k_fld) :: dens_temp, cons, inj
  logical :: shoot
  real(p_k_fld) :: lab_dt

  shoot = .false.
  ion_idx = this%ion_idx
  lab_dt = dt


  do theta = 1, this%num_par_theta

    !$omp parallel do private(i,cons,shoot,m,inj,dens_temp)
    do j=1, this%multi_ion_cyl(theta) % nx_( 2 )
      do i=1, this%multi_ion_cyl(theta)% nx_( 1 )

        if ( this%multi_ion_cyl(theta)%f2(this%neut_idx,i,j) > this%den_min ) then

          ! 1st order method
          ! cons = this(n)%multi_ion%f2(1,i,j) * w_ion(1) * dt

          ! 2nd order Runge-Kutta
          cons = this%multi_ion_cyl(theta)%f2(1,i,j) * this%w_cyl(theta)%f2(1,i,j) * lab_dt / &
                ( 1.0_p_k_fld + 0.5_p_k_fld * this%w_cyl(theta)%f2(1,i,j) * lab_dt )

          ! overshoot
          if (cons > this%multi_ion_cyl(theta)%f2(1,i,j)) then
            shoot = .true.
            cons = this%multi_ion_cyl(theta)%f2(1,i,j)




          endif

          ! subtract from 0.
          this%multi_ion_cyl(theta)%f2(1,i,j) = this%multi_ion_cyl(theta)%f2(1,i,j) - cons

          do m=2, this%multi_max

            ! remember how much charge will be added to m-1. level
            inj = cons

            ! overshoot in previous level (we go to time centered densities)
            if (shoot) then
              dens_temp = inj * 0.5_p_k_fld
              shoot = .false.
            else
              dens_temp = 0.0_p_k_fld
            endif

            ! ionizing m-1. level (index m) adding to m. level

            ! 1st order method
            !cons = (this(n)%multi_ion%f2(m,i,j) + dens_temp)* w_ion(m) * dt

            ! 2nd order Runge-Kutta
            cons = (this%multi_ion_cyl(theta)%f2(m,i,j) + dens_temp)* this%w_cyl(theta)%f2(m,i,j) * lab_dt / &
                  ( 1.0_p_k_fld + 0.5_p_k_fld * this%w_cyl(theta)%f2(m,i,j) * lab_dt )

            ! overshoot in current level (temp might come from previous overshoot)
            if (cons > this%multi_ion_cyl(theta)%f2(m,i,j) + dens_temp) then
              shoot = .true.
              cons = this%multi_ion_cyl(theta)%f2(m,i,j) + dens_temp




            endif

            ! update density for m-1. level (index m)
            this%multi_ion_cyl(theta)%f2(m,i,j) = min(this%multi_ion_cyl(theta)%f2(m,i,j) - &
                                                  cons + inj, 1.0_p_k_fld)
            ! min only due to round off

          enddo ! all levels

          ! update the last level
          this%multi_ion_cyl(theta)%f2(this%multi_max+1,i,j) = &
                this%multi_ion_cyl(theta)%f2(this%multi_max+1,i,j) + cons

          ! calculate total charge
          this%multi_ion_cyl(theta)%f2(ion_idx,i,j) = 0.0_p_k_fld
          do m=2, this%multi_max+1
            this%multi_ion_cyl(theta)%f2(ion_idx,i,j) = &
                  this%multi_ion_cyl(theta)%f2(ion_idx,i,j) + &
                  (m-1)*(this%multi_ion_cyl(theta)%f2(m,i,j))
          enddo

          ! again: just round off
          if (this%multi_ion_cyl(theta)%f2(ion_idx,i,j) > this%multi_max) then
            this%multi_ion_cyl(theta)%f2(ion_idx,i,j) = this%multi_max
          endif

        endif

      enddo
    enddo
    !$omp end parallel do

  enddo

end subroutine advance_density_cyl_modes_std
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine advance_density_cyl_modes_analytic(this , dt)
!-----------------------------------------------------------------------------------------

  implicit none

  class(t_neutral_cyl_modes), intent(inout) :: this ! neutral array
  real( p_k_fld ), intent(in) :: dt
  integer :: i, j, k, m, l, ion_idx, theta
  real(p_k_fld) :: csum
  real(p_k_fld) :: lab_dt
  real(p_k_fld), allocatable, dimension(:) :: w_exp, coefs, multi_ion_next

  ion_idx = this%ion_idx
  lab_dt = dt


  ! allocate temporary buffers
  allocate(w_exp(this%multi_max + 1))
  allocate(coefs(this%multi_max + 1))
  allocate(multi_ion_next(this%multi_max + 1))

  do theta = 1, this%num_par_theta

    !$omp parallel do private(i,m, l, w_exp, coefs, multi_ion_next)
    do j=1, this%multi_ion_cyl(theta) % nx_( 2 )
      do i=1, this%multi_ion_cyl(theta)% nx_( 1 )
        if ( this%multi_ion_cyl(theta)%f2(this%neut_idx,i,j) > this%den_min) then

          ! reset temp array level m
          multi_ion_next(1:this%multi_max+1) = 0.0_p_k_part
          l = 1

          do m= 1, this%multi_max + 1
            if(m > 1 .and. this%w_cyl(theta)%f2(m,i,j) == 0.0_p_k_part .and. this%w_cyl(theta)%f2(m-1,i,j) == 0.0_p_k_part ) then
              l = m + 1
              CYCLE
            endif

            ! calculate exp(-rate[level m] * dt) and store (only done once)
            w_exp(m) = exp(- this%w_cyl(theta)%f2(m,i,j) * lab_dt )

            ! coef sum from 1 to m-1
            csum = 0.0

            do k= l, m-1
              ! update coef of exp(-rate[k] * dt)
              coefs(k) = coefs(k) * this%w_cyl(theta)%f2(m-1,i,j)/ &
              (this%w_cyl(theta)%f2(m,i,j)- this%w_cyl(theta)%f2(k,i,j))

              ! multiple coef by exponential exp(-rate[k] * dt) and sum to current level m
              multi_ion_next(m) = multi_ion_next(m) + w_exp(k) * coefs(k)
              csum = csum + coefs(k)
            enddo

            coefs(m) = this%multi_ion_cyl(theta)%f2(m,i,j) - csum
            multi_ion_next(m) = multi_ion_next(m) + w_exp(m) * coefs(m)
            this%multi_ion_cyl(theta)%f2(m,i,j) = multi_ion_next(m)

          enddo ! all levels

          ! calculate total charge
          this%multi_ion_cyl(theta)%f2(ion_idx,i,j) = 0.0_p_k_fld
          do m=2, this%multi_max+1
            this%multi_ion_cyl(theta)%f2(ion_idx,i,j) = &
                  this%multi_ion_cyl(theta)%f2(ion_idx,i,j) + &
                  (m-1)*(this%multi_ion_cyl(theta)%f2(m,i,j))
          enddo

          ! again: just round off
          if (this%multi_ion_cyl(theta)%f2(ion_idx,i,j) > this%multi_max) then
            this%multi_ion_cyl(theta)%f2(ion_idx,i,j) = this%multi_max
          endif


        endif

      enddo
    enddo
    !$omp end parallel do


  enddo

  ! free temp buffers
  deallocate( w_exp )
  deallocate( coefs )
  deallocate( multi_ion_next )

end subroutine advance_density_cyl_modes_analytic
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
subroutine inject_particles_cyl_modes(this, coordinates)
!-----------------------------------------------------------------------------------------

  implicit none

  integer, parameter :: rank = 2

  class(t_neutral_cyl_modes), intent(inout) :: this
  integer, intent(in) :: coordinates
  integer, dimension(rank) :: num_par !=num_par_x
  integer :: n_p_c, n_p_add, i, j, ix, ion_idx

  real(p_k_part), dimension(4) :: xnewpart
  integer, dimension(rank) :: ixnewpart

  real(p_k_part), dimension(p_p_dim) :: pnewpart
  real(p_k_part) :: norm1, norm2, qnewpart1, qnewpart2
  real(p_k_part) :: den_center
  real(p_k_part) :: r, dr
  integer :: shift_ix2, theta
  class(t_species_cyl_modes), pointer :: spec1, spec2

  ! charge normalization coefficient for injecting particles
  real(p_k_part) :: alpha

  select type( spec1_temp => this%species1 )
  class is ( t_species_cyl_modes )
    spec1 => spec1_temp
  class default
    ! This must never happen, inject_particles_cyl_modes must use a species1 that is a
    ! t_species_cyl_modes object
    call abort_program( p_err_invalid )
  end select

  select type( spec2_temp => this%species2 )
  class is ( t_species_cyl_modes )
    ! if castable
    spec2 => spec2_temp
  class default
    ! if not castable then only abort if_mov_ions
    if(this%if_mov_ions) then
      call abort_program( p_err_invalid )
    endif
  end select

  pnewpart = 0.0_p_k_part
  ion_idx = this%ion_idx
  num_par(1) = spec1%num_par_x(1)
  num_par(2) = spec1%num_par_x(2)
  n_p_c = num_par(1) * num_par(2)
  norm1 = -1.0_p_k_part / (n_p_c*this%num_par_theta)
  dr = real( spec1%dx(2), p_k_part )
  shift_ix2 = spec1%my_nx_p(p_lower, 2) - 2

  ! Note that the number of particles per cell of the ion species is
  ! being ignored (which is correct )
  if (this%if_mov_ions) norm2 = +1.0_p_k_part / (n_p_c*this%num_par_theta)

  do theta = 1, this%num_par_theta

    do j=1, this%multi_ion_cyl(theta) % nx_( 2 )

      ! If injected into the middle of the box, we can correct charge normalization
      ! If injected randomly, there will be a small error that we cannot correct
      if ( this%inject_line ) then
        r = spec1%g_box( p_lower , p_r_dim ) + ( j + shift_ix2 ) * dr
        alpha = spec1 % norm_charge_axis( j, r, ppc2_in=1 )
      else
        alpha = 1.0_p_k_part
      endif

      do i=1, this%multi_ion_cyl(theta) % nx_( 1 )

        if ( this%multi_ion_cyl(theta)%f2(ion_idx,i,j) > this%multi_min ) then

          ! determine neutral density
          den_center = this%multi_ion_cyl(theta)%f2(this%neut_idx,i,j)

          if (den_center > this%den_min) then

            n_p_add = int(this%multi_ion_cyl(theta)%f2(ion_idx,i,j)*n_p_c/this%multi_max + 0.5_p_k_fld) - &
                      int(this%ion_level_old_cyl(theta)%f2(1,i,j)*n_p_c/this%multi_max + 0.5_p_k_fld)

            ! determine charge of the new particle
            qnewpart1 = real( this%multi_max*den_center*norm1*alpha, p_k_part )

            if (this%if_mov_ions) qnewpart2 = real( den_center*norm2*alpha, p_k_part )

            ! place particles in cell
            do ix = 0, n_p_add-1

              ixnewpart(1) = i
              ixnewpart(2) = j

              if ( this%inject_line ) then
                xnewpart(1) = (ix + 0.5)/n_p_add - 0.5_p_k_part
                xnewpart(2) = 0
              else
                ! using genrand_real3() makes sure that the particle is
                ! never injected in the cell boundary
                call rng % harvest_real3( xnewpart(1), ixnewpart )
                call rng % harvest_real3( xnewpart(2), ixnewpart )

                xnewpart(1) = xnewpart(1) - 0.5_p_k_part
                xnewpart(2) = xnewpart(2) - 0.5_p_k_part
              endif

              ! get radial position normalized to cell size
              ! since g_box is adjusted based on interpolation order, this will always
              ! give the correct particle position

              ! get radial position
              r = spec1%g_box( p_lower , p_r_dim ) + &
                  ( ( ixnewpart(2) + shift_ix2 ) + xnewpart(2) ) * dr

              ! Only inject inside the box. This could be optimized for near cell positions
              ! by only injecting starting from cell 2
              if ( r > 0 ) then

                xnewpart(3) = r*this%cosine(1,i,j,theta)
                xnewpart(4) = r*this% sine(1,i,j,theta)

                ! add particles to the corresponding buffers
                call spec1%create_particle_single_cell_p_cyl( ixnewpart, xnewpart, &
                                                              pnewpart, qnewpart1*r )

                if (this%if_mov_ions) then
                  ! add particles to the ion buffer
                  call spec2%create_particle_single_cell_p_cyl( ixnewpart, xnewpart, &
                                                                pnewpart, qnewpart2*r )

                endif
              endif

            enddo ! n_p_add
          endif ! den_center > den_min
        endif
      enddo ! nx1
    enddo ! nx2
  enddo !theta

end subroutine inject_particles_cyl_modes
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
! Ionizes a background source
!-----------------------------------------------------------------------------------------
subroutine ionize_neutral_cyl_modes(this, species, emf, gdt, coordinates)
!-----------------------------------------------------------------------------------------

  implicit none

  class(t_neutral_cyl_modes) :: this ! neutral obj
  class(t_species), pointer :: species ! species array (impact ionization only)
  class( t_emf ), intent(in) :: emf ! EMF data
  real( p_double ), intent(in) :: gdt ! time step

  integer, intent(in) :: coordinates

  real( p_k_fld ) :: dt
  logical :: if_impact
  integer :: theta , num_par_theta

  dt = real( gdt, p_k_fld )
  if_impact = .false.

  num_par_theta = this%num_par_theta
  ! Clear ionization rates
  !$omp parallel do
  do theta =1,num_par_theta
    call this%w_cyl(theta)%zero()

    ! total ion density is the last component of multi_ion vdf
    call this%ion_level_old_cyl(theta)%copy(this%multi_ion_cyl(theta), this%ion_idx)

  enddo
  !$omp end parallel do

  ! field ionization rates
  if ( this%if_tunnel ) then
      call this% adk_field_rates( emf )
  endif

  ! impact ionization rates
  ! if ( this(n)%if_impact ) call impact_rates( this(n), species, dt )
  ! Advance densities
  !uncomment later
  call this % advance( dt )

  ! Inject particles
  call this % inject(coordinates)


end subroutine ionize_neutral_cyl_modes
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
! Get ADK ionization rates on each cell in 2D
!-----------------------------------------------------------------------------------------
subroutine adk_field_rates_neutral_cyl_modes( this, emf )

  implicit none
  class( t_neutral_cyl_modes ), intent(inout) :: this ! neutral
  class( t_emf ), intent(in) :: emf ! electrical field

  type(t_vdf), pointer :: e_re, e_im
  integer :: i, j, l
  real(p_k_fld) :: den_center, eij, e1, e2, e3
  integer :: num_par_theta, theta , mode, num_modes

  select type ( emf )
  class is ( t_emf_cyl_modes )

    ! load values
    num_par_theta = this%num_par_theta
    num_modes = emf%n_cyl_modes

    do theta = 1,num_par_theta

      if ( this%species1%pos_type == p_cell_low ) then

        !$omp parallel do private(i,den_center,e_re,e1,e2,e3,mode,e_im,eij,l)
        do j = 1, emf%e_cyl_m%pf_re(0)%nx_( 2 )
          do i = 1, emf%e_cyl_m%pf_re(0)%nx_( 1 )

            den_center = this%multi_ion_cyl(theta)%f2(this%neut_idx,i,j)

            if (den_center > this%den_min) then
              ! Interpolate at center of the cell
              e_re => emf%e_part_cyl_m%pf_re(0)
              e1 = 0.5_p_k_fld*( e_re%f2(1,i,j) + e_re%f2(1,i,j+1) )
              e2 = 0.5_p_k_fld*( e_re%f2(2,i,j) + e_re%f2(2,i+1,j) )
              e3 = 0.25_p_k_fld*( e_re%f2(3,i,j) + e_re%f2(3,i+1,j) + &
                                  e_re%f2(3,i,j+1) + e_re%f2(3,i+1,j+1))

              do mode = 1, num_modes
                e_re => emf%e_part_cyl_m%pf_re(mode)
                e_im => emf%e_part_cyl_m%pf_im(mode)
                e1 = e1 + 0.5_p_k_fld*( this%cosine(mode,i,j,theta)*( e_re%f2(1,i,j) + e_re%f2(1,i,j+1) ) &
                                       + this%sine(mode,i,j,theta)* ( e_im%f2(1,i,j) + e_im%f2(1,i,j+1) ) )
                e2 = e2 + 0.5_p_k_fld*( this%cosine(mode,i,j,theta)*( e_re%f2(2,i,j) + e_re%f2(2,i+1,j) ) &
                                       + this%sine(mode,i,j,theta)* ( e_im%f2(2,i,j) + e_im%f2(2,i+1,j) ) )
                e3 = e3 + 0.25_p_k_fld*( this%cosine(mode,i,j,theta)*( e_re%f2(3,i,j) + e_re%f2(3,i+1,j) + &
                                                                       e_re%f2(3,i,j+1) + e_re%f2(3,i+1,j+1)) &
                                       + this%sine(mode,i,j,theta)* ( e_im%f2(3,i,j) + e_im%f2(3,i+1,j) + &
                                                                       e_im%f2(3,i,j+1) + e_im%f2(3,i+1,j+1)) )
              enddo

              eij = sqrt(e1**2 + e2**2 + e3**2)*this%omega_p*1.704e-12
              if (eij > this%e_min) then
                ! w = r1 * eij^-r3 * EXP(-r2/eij) * 1/wp
                do l = 1, this%multi_max
                  this%w_cyl(theta)%f2(l,i,j) = &
                      this%rate_param(1,l) * &
                      eij**(-this%rate_param(3,l)) * exp(-this%rate_param(2,l)/(eij)) / &
                      this%omega_p
                enddo
              endif
            else
              this%w_cyl(theta)%f2(1:this%multi_max,i,j) = 0.0
            endif

          enddo
        enddo
        !$omp end parallel do

      else

        !$omp parallel do private(i,den_center,e_re,e1,e2,e3,mode,e_im,eij,l)
        do j = 1, emf%e_cyl_m%pf_re(0)%nx_( 2 )
          do i = 1, emf%e_cyl_m%pf_re(0)%nx_( 1 )

            den_center = this%multi_ion_cyl(theta)%f2(this%neut_idx,i,j)

            if (den_center > this%den_min) then
              ! Interpolate at cell corner
              e_re => emf%e_part_cyl_m%pf_re(0)
              e1 = 0.5_p_k_fld*( e_re%f2(1,i-1,j) + e_re%f2(1,i,j) )
              e2 = 0.5_p_k_fld*( e_re%f2(2,i,j-1) + e_re%f2(2,i,j) )
              e3 = e_re%f2(3,i,j)

              do mode = 1, num_modes
                e_re => emf%e_part_cyl_m%pf_re(mode)
                e_im => emf%e_part_cyl_m%pf_im(mode)
                e1 = e1 + 0.5_p_k_fld*( this%cosine(mode,i,j,theta)*( e_re%f2(1,i-1,j) + e_re%f2(1,i,j) ) &
                                      + this%sine(mode,i,j,theta)* ( e_im%f2(1,i-1,j) + e_im%f2(1,i,j) ) )
                e2 = e2 + 0.5_p_k_fld*( this%cosine(mode,i,j,theta)*( e_re%f2(2,i,j-1) + e_re%f2(2,i,j) ) &
                                      + this%sine(mode,i,j,theta)* ( e_im%f2(2,i,j-1) + e_im%f2(2,i,j) ) )
                e3 = e3 + this%cosine(mode,i,j,theta)*e_re%f2(3,i,j) + this%sine(mode,i,j,theta)*e_im%f2(3,i,j)
              enddo

              eij = sqrt(e1**2 + e2**2 + e3**2)*this%omega_p*1.704e-12
              if (eij > this%e_min) then
                ! w = r1 * eij^-r3 * EXP(-r2/eij) * 1/wp
                do l = 1, this%multi_max
                  this%w_cyl(theta)%f2(l,i,j) = &
                      this%rate_param(1,l) * &
                      eij**(-this%rate_param(3,l)) * exp(-this%rate_param(2,l)/(eij)) / &
                      this%omega_p
                enddo
              endif
            else
              this%w_cyl(theta)%f2(1:this%multi_max,i,j) = 0.0
            endif

          enddo
        enddo
        !$omp end parallel do

      endif

    enddo

  end select

end subroutine adk_field_rates_neutral_cyl_modes
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
! Read checkpoint information from file
!-----------------------------------------------------------------------------------------
subroutine restart_read_neutral_cyl_modes( this, restart_handle )
!-----------------------------------------------------------------------------------------


  implicit none

  class( t_neutral_cyl_modes ), intent(inout) :: this
  type( t_restart_handle ), intent(in) :: restart_handle

  character(len=*), parameter :: err_msg = 'error reading restart data for neutral object.'
  character(len=len(p_neutral_cyl_modes_rst_id)) :: rst_id
  integer :: ierr, n


  call restart_io_read(rst_id, restart_handle, ierr)
  call check_error(ierr,err_msg,p_err_rstrd,"cyl_modes/os-neutral-cyl-modes.f03",1043)

  ! check if restart file is compatible
  if ( rst_id /= p_neutral_cyl_modes_rst_id) then
    write(err_buf__,*) 'Corrupted restart file, or restart file, from incompatible binary (neutral)';call err__("cyl_modes/os-neutral-cyl-modes.f03",1047)
    call abort_program(p_err_rstrd)
  endif

  ! diagnostic this%diag is not read, so that the user
  ! can change the dump factors before restarting

  call restart_io_read(this%neutral_id, restart_handle, ierr)
  call check_error(ierr,err_msg,p_err_rstrd,"cyl_modes/os-neutral-cyl-modes.f03",1055)

  call restart_io_read(this%name, restart_handle, ierr)
  call check_error(ierr,err_msg,p_err_rstrd,"cyl_modes/os-neutral-cyl-modes.f03",1058)

  call restart_io_read(this%if_mov_ions, restart_handle, ierr)
  call check_error(ierr,err_msg,p_err_rstrd,"cyl_modes/os-neutral-cyl-modes.f03",1061)

  call restart_io_read(this%omega_p, restart_handle, ierr)
  call check_error(ierr,err_msg,p_err_rstrd,"cyl_modes/os-neutral-cyl-modes.f03",1064)

  call restart_io_read(this%den_min, restart_handle, ierr)
  call check_error(ierr,err_msg,p_err_rstrd,"cyl_modes/os-neutral-cyl-modes.f03",1067)

  call restart_io_read(this%if_tunnel, restart_handle, ierr)
  call check_error(ierr,err_msg,p_err_rstrd,"cyl_modes/os-neutral-cyl-modes.f03",1070)

  call restart_io_read(this%if_impact, restart_handle, ierr)
  call check_error(ierr,err_msg,p_err_rstrd,"cyl_modes/os-neutral-cyl-modes.f03",1073)

  call restart_io_read(this%multi_max, restart_handle, ierr)
  call check_error(ierr,err_msg,p_err_rstrd,"cyl_modes/os-neutral-cyl-modes.f03",1076)

  call restart_io_read(this%num_par_theta, restart_handle, ierr)
  call check_error(ierr,err_msg,p_err_rstrd,"cyl_modes/os-neutral-cyl-modes.f03",1079)

  call alloc(this%multi_ion_cyl, (/ this%num_par_theta /),"cyl_modes/os-neutral-cyl-modes.f03",1081)
  call alloc(this%ion_level_old_cyl, (/ this%num_par_theta /),"cyl_modes/os-neutral-cyl-modes.f03",1082)
  call alloc(this%w_cyl, (/ this%num_par_theta /),"cyl_modes/os-neutral-cyl-modes.f03",1083)
  do n=1, this%num_par_theta
    call this % multi_ion_cyl(n) % read_checkpoint( restart_handle )
    call this % ion_level_old_cyl(n) % read_checkpoint( restart_handle )
    call this % w_cyl(n) % read_checkpoint( restart_handle )
  enddo

  if ( this%if_impact ) then
    call restart_read( this%cross_section, restart_handle )
  endif

  call this % rand_theta % read_checkpoint( restart_handle )

end subroutine restart_read_neutral_cyl_modes
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! write object information into a restart file
!-----------------------------------------------------------------------------------------
subroutine restart_write_neutral_cyl_modes( this, restart_handle )
!-----------------------------------------------------------------------------------------

  implicit none

  class( t_neutral_cyl_modes ), intent(in) :: this
  type( t_restart_handle ), intent(inout) :: restart_handle

  character(len=*), parameter :: err_msg = 'error writing restart data for neutral object.'
  integer :: ierr, n

  call restart_io_write("p_neutral_cyl_modes_rst_id", p_neutral_cyl_modes_rst_id, restart_handle, ierr)
  call check_error(ierr,err_msg,p_err_rstwrt,"cyl_modes/os-neutral-cyl-modes.f03",1114)

  ! diagnostic this%diag is not saved, so that the user
  ! can change the dump factors before restarting
  call restart_io_write("this%neutral_id", this%neutral_id, restart_handle, ierr)
  call check_error(ierr,err_msg,p_err_rstwrt,"cyl_modes/os-neutral-cyl-modes.f03",1119)

  call restart_io_write("this%name", this%name, restart_handle, ierr)
  call check_error(ierr,err_msg,p_err_rstwrt,"cyl_modes/os-neutral-cyl-modes.f03",1122)

  call restart_io_write("this%if_mov_ions", this%if_mov_ions, restart_handle, ierr)
  call check_error(ierr,err_msg,p_err_rstwrt,"cyl_modes/os-neutral-cyl-modes.f03",1125)

  call restart_io_write("this%omega_p", this%omega_p, restart_handle, ierr)
  call check_error(ierr,err_msg,p_err_rstwrt,"cyl_modes/os-neutral-cyl-modes.f03",1128)

  call restart_io_write("this%den_min", this%den_min, restart_handle, ierr)
  call check_error(ierr,err_msg,p_err_rstwrt,"cyl_modes/os-neutral-cyl-modes.f03",1131)

  call restart_io_write("this%if_tunnel", this%if_tunnel, restart_handle, ierr)
  call check_error(ierr,err_msg,p_err_rstwrt,"cyl_modes/os-neutral-cyl-modes.f03",1134)

  call restart_io_write("this%if_impact", this%if_impact, restart_handle, ierr)
  call check_error(ierr,err_msg,p_err_rstwrt,"cyl_modes/os-neutral-cyl-modes.f03",1137)

  call restart_io_write("this%multi_max", this%multi_max, restart_handle, ierr)
  call check_error(ierr,err_msg,p_err_rstwrt,"cyl_modes/os-neutral-cyl-modes.f03",1140)

  call restart_io_write("this%num_par_theta", this%num_par_theta, restart_handle, ierr)
  call check_error(ierr,err_msg,p_err_rstwrt,"cyl_modes/os-neutral-cyl-modes.f03",1143)

  do n=1, this%num_par_theta
    call this % multi_ion_cyl(n) % write_checkpoint( restart_handle )
    call this % ion_level_old_cyl(n) % write_checkpoint( restart_handle )
    call this % w_cyl(n) % write_checkpoint( restart_handle )
  enddo

  if ( this%if_impact ) then
    call restart_write( this%cross_section, restart_handle)
  endif

  call this % rand_theta % write_checkpoint( restart_handle )

end subroutine restart_write_neutral_cyl_modes
!-----------------------------------------------------------------------------------------

end module m_neutral_cyl_modes
