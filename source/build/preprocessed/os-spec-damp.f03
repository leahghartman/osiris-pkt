# 1 "overdense/os-spec-damp.f03"
# 1 "<built-in>" 1
# 1 "<built-in>" 3
# 467 "<built-in>" 3
# 1 "<command line>" 1
# 1 "<built-in>" 2
# 1 "overdense/os-spec-damp.f03" 2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Module for particle damper
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
# 6 "overdense/os-spec-damp.f03" 2
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
# 7 "overdense/os-spec-damp.f03" 2

module m_spec_damp

use m_system
use m_parameters
use m_species_define
use m_random
use m_input_file
use m_grid_define
use m_spec_fluid_moments
use m_node_conf
use m_space, only : t_space, xmin, xmax
use m_emf_define, only : t_emf
use m_current_define, only : t_current
use m_vdf_define, only : t_vdf
use m_vdf_comm, only : t_vdf_msg
use m_restart, only : t_restart_handle
use m_time_step
use m_fparser, only : t_fparser, p_max_expr_len, setup

implicit none

private

integer, parameter :: p_all = 0
integer, parameter :: p_forward = 1
integer, parameter :: p_backward = 2

!-----------------------------------------------------------------------------------------
! Particle damper 'class'
!-----------------------------------------------------------------------------------------
type t_damper

  integer :: damper_direction
  integer :: orientation ! forward or backward

  real(p_k_part) :: mfp ! mean free path, used for the standard input
  real(p_double) :: residual ! fraction of particles that make it through the absorber
  logical :: custom_input ! if yes, the user is manually setting x_start, x_full,
                          ! stopping_dist_start and stopping_dist_full

  real(p_k_part) :: x_start ! the x value at which absorption begins
  real(p_k_part) :: x_full ! the x value at whiche absorption reaches it maximum value
  real(p_k_part) :: x_range ! the difference of the two above

  ! the gamma cutoff, below which particles are left untouched
  real(p_k_part) :: gamma_start ! the cutoff gamma at the start of the damper
  real(p_k_part) :: gamma_full ! the cutoff gamma in the maximum damper region
                                 ! note this is usually less than gamma_start, confusingly
  real(p_k_part) :: gamma_range ! the difference of the two above

  logical :: local_vth_cutoff ! define the damping cutoff as a multiple of the local temperature
  real(p_k_part) :: vth_mult_start
  real(p_k_part) :: vth_mult_full
  real(p_k_part) :: vth_mult_range

  logical :: reemit_with_local_temperature

  logical :: fourth_root_temperature ! use the fourth root to calculate the temperature
                                     ! affects both local_vth_cutoff and reemit_with_local_temperature

  real(p_k_part), dimension(p_p_dim) :: reemit_temp_start ! the initial reemision temperature
  real(p_k_part), dimension(p_p_dim) :: reemit_temp_full ! the reemision temp in the maximum absorption region
                                                        ! again less than reemite_temp_full
  real(p_k_part), dimension(p_p_dim) :: reemit_temp_range ! the difference of the two above

  real(p_k_part) :: stopping_dist_start ! the (target) stopping distance at the damper start
  real(p_k_part) :: stopping_dist_full ! the stopping distance in the maximum absortion region
                                      ! again, a bad name. you think I'd fix that.
  real(p_k_part) :: stopping_dist_range ! the difference of the two above

  integer :: n_damp ! the frequency of absorption events, as a multiple of actual timesteps
                      ! this is only to decrease computational cost, the stopping distance is
                      ! adjusted to correct for this

  logical :: if_cold_perp ! flag for whether to re-emit thermally in all three dimensions (.false.)
                          ! or only in the physical directions of the simulation (.true.)

  logical :: if_damp ! if this species has an damper anywhere
                       ! both set automatically and can be overridden by the input deck
  logical :: if_damp_local ! if the damper is on this node, again just to save on computation time
                             ! (Though if you don't load balance correctly, the nodes w/o dampers
                             ! (will just be waiting idlely for the nodes with them,
                             ! (but it's the thought that counts )
  logical :: linear ! whether to use damping based on velocity magnitude (.true.) or based
                    ! on the hazard function (.false.). It is recommended linear=.true.
  integer :: stop_dir, reemit_dir ! stop/reemit particles traveling in a preferred direction
  real(p_double) :: damp_start_time ! time to start damping

  logical :: use_local_vth_func ! whether to use (.true.) or not use (.false.) the below
                                ! function specifying where to/not to calculate the local
                                ! thermal velocity
  type(t_fparser) :: local_vth_func ! function (of space) defining where to use (>0) or
                                    ! not use (<=0) the local_vth_cutoff and
                                    ! reemit_with_local_temperature parameters

  type(t_damper), pointer :: next => null()

end type t_damper

interface init
  module procedure init_damper
end interface

interface read_nml
  module procedure read_nml_damper
end interface

public :: t_damper, p_all, p_forward, p_backward
public :: if_damp, init, read_nml

contains

!-----------------------------------------------------------------------------------------
function if_damp( damper, n, t )

  type( t_damper ), intent(in) :: damper
  integer, intent(in) :: n
  real(p_double), intent(in) :: t
  logical :: if_damp

  if_damp = .false.

  if( t < damper%damp_start_time ) return

  if( .not. damper%if_damp ) return

  if( .not. damper%if_damp_local ) return

  if( damper%n_damp > 0 ) then
    if( mod(n, damper%n_damp) .ne. 0 ) return
  endif

  if_damp = .true.
  return

end function if_damp
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Read information from input file
!-----------------------------------------------------------------------------------------
subroutine read_nml_damper( this, input_file )

  type( t_damper ), intent(inout) :: this
  class( t_input_file ), intent(inout) :: input_file

  logical :: if_damp
  integer :: n_damp
  integer :: damper_direction
  character(len=16) :: orientation
  logical :: custom_input
  real(p_k_part) :: mfp
  real(p_double) :: residual
  real(p_k_part) :: x_start
  real(p_k_part) :: x_full
  real(p_k_part) :: gamma_start
  real(p_k_part) :: gamma_full
  logical :: local_vth_cutoff
  real(p_k_part) :: vth_mult_start
  real(p_k_part) :: vth_mult_full
  logical :: reemit_with_local_temperature, fourth_root_temperature
  real(p_k_part), dimension(p_p_dim) :: reemit_temp_start
  real(p_k_part), dimension(p_p_dim) :: reemit_temp_full
  real(p_k_part) :: stopping_dist_start
  real(p_k_part) :: stopping_dist_full
  logical :: if_cold_perp
  logical :: linear
  character(len=16) :: stop_dir, reemit_dir
  real(p_double) :: damp_start_time
  logical :: use_local_vth_func
  character(len=p_max_expr_len) :: local_vth_func_expr

  integer :: ierr, i

  namelist /nl_damper/ if_damp, n_damp, damper_direction, x_start, x_full, gamma_start,&
                        gamma_full, local_vth_cutoff, vth_mult_start, vth_mult_full, &
                        reemit_with_local_temperature, fourth_root_temperature, &
                        reemit_temp_start, reemit_temp_full, stopping_dist_start, &
                        stopping_dist_full, if_cold_perp, linear, stop_dir, reemit_dir, &
                        damp_start_time, use_local_vth_func, local_vth_func_expr, &
                        orientation, custom_input, mfp, residual

  if_damp = .true.
  n_damp = 0
  damper_direction = 1
  orientation = "forward"
  mfp = -1.0
  residual = 1.0d-7
  custom_input = .false.
  x_start = 0.0
  x_full = 0.0
  gamma_start = -1.0
  gamma_full = -1.0
  local_vth_cutoff = .true.
  vth_mult_start = 6.0
  vth_mult_full = 6.0
  reemit_with_local_temperature = .true.
  fourth_root_temperature = .false.
  reemit_temp_start = -1.0
  reemit_temp_full = -1.0
  stopping_dist_start = 0.0
  stopping_dist_full = 0.0
  if_cold_perp = .false.
  linear = .true.
  stop_dir = "all"
  reemit_dir = "all"
  damp_start_time = -1.0_p_double
  use_local_vth_func = .false.
  local_vth_func_expr = "NO_FUNCTION_SUPPLIED!"

  call get_namelist(input_file, "nl_damper", ierr)

  if ( ierr == 0 ) then
    read (input_file%nml_text, nml = nl_damper, iostat = ierr)
    if (ierr /= 0) then
      if ( mpi_node() == 0 ) then
        write(0,*) ""
        write(0,*) "   Error reading damper parameters"
        write(0,*) "   aborting..."
      endif
      stop
    endif
  else
    this%if_damp = .false.
    return
  endif

  if( .not. if_damp ) then
    this%if_damp = .false.
    return
  endif

  this%if_damp = if_damp
  this%damper_direction = damper_direction

  select case ( trim(orientation) )
  case ("forward")
    this%orientation = p_forward
  case ("backward")
    this%orientation = p_backward
  case default
    if ( mpi_node() == 0 ) then
      write(0,*) ''
      write(0,*) '   Error in damper parameters'
      write(0,*) '   orientation must be either "forward" or "backward"'
      write(0,*) '   aborting...'
    endif
    stop
  end select

  select case ( trim(stop_dir) )
  case ("all")
    this%stop_dir = p_all
  case ("toward boundary")
    this%stop_dir = p_forward
  case ("away from boundary")
    this%stop_dir = p_backward
  case default
    if ( mpi_node() == 0 ) then
      write(0,*) ''
      write(0,*) '   Error in damper parameters'
      write(0,*) '   stop_dir must be either "all", "toward boundary" or "away from boundary"'
      write(0,*) '   aborting...'
    endif
    stop
  end select

  select case ( trim(reemit_dir) )
  case ("all")
    this%reemit_dir = p_all
  case ("toward boundary")
    this%reemit_dir = p_forward
  case ("away from boundary")
    this%reemit_dir = p_backward
  case default
    if ( mpi_node() == 0 ) then
      write(0,*) ''
      write(0,*) '   Error in damper parameters'
      write(0,*) '   reemit_dir must be either "all", "toward boundary" or "away from boundary"'
      write(0,*) '   aborting...'
    endif
    stop
  end select

  this%custom_input = custom_input
  if ( custom_input ) then

    this%x_start = x_start
    this%x_full = x_full

    this%stopping_dist_start = stopping_dist_start
    this%stopping_dist_full = stopping_dist_full

    if( (this%x_full - this%x_start) .le. 0.0 .and. this%orientation == p_forward ) then
      write(err_buf__,*) "Damper x_range must be > 0 for forward orientation";call err__("overdense/os-spec-damp.f03",302)
      call abort_program( p_err_invalid )
    elseif( (this%x_full - this%x_start) .ge. 0.0 .and. this%orientation == p_backward ) then
      write(err_buf__,*) "Damper x_range must be < 0 for backward orientation";call err__("overdense/os-spec-damp.f03",305)
      call abort_program( p_err_invalid )
    endif

    if( this%stopping_dist_start .le. 0.0 ) then
      write(err_buf__,*) "Damper stopping distances must be positive";call err__("overdense/os-spec-damp.f03",310)
      call abort_program( p_err_invalid )
    endif

    if( this%stopping_dist_full .le. 0.0 ) then
      write(err_buf__,*) "Damper stopping distances must be positive";call err__("overdense/os-spec-damp.f03",315)
      call abort_program( p_err_invalid )
    endif

  else

    ! This is a custom absorber profile picked by Kyle that should work well in most cases
    ! Here the mean free path is set to be approximately 35% of the entire absorber length
    ! Users simply specify their desired mean free path, and what fraction of particles
    ! they are okay with making it entirely through the boundary. Kyle recommends this to
    ! be less than 1e-6.
    this%mfp = mfp
    this%residual = residual

    if( this%mfp .le. 0.0 ) then
      write(err_buf__,*) "Damper mean free path must be positive";call err__("overdense/os-spec-damp.f03",330)
      call abort_program( p_err_invalid )
    endif

    if( this%residual .le. 0.0 .or. this%residual .ge. 0.001 ) then
      write(err_buf__,*) "Damper residual must be 0 < residual < 0.001";call err__("overdense/os-spec-damp.f03",335)
      call abort_program( p_err_invalid )
    endif

  endif

  this%gamma_start = gamma_start
  this%gamma_full = gamma_full

  this%local_vth_cutoff = local_vth_cutoff
  this%vth_mult_start = vth_mult_start
  this%vth_mult_full = vth_mult_full

  this%reemit_with_local_temperature = reemit_with_local_temperature
  this%fourth_root_temperature = fourth_root_temperature
  this%reemit_temp_start = reemit_temp_start
  this%reemit_temp_full = reemit_temp_full

  this%if_cold_perp = if_cold_perp
  this%linear = linear
  this%damp_start_time = damp_start_time
  this%n_damp = n_damp

  this%use_local_vth_func = use_local_vth_func
  ! Setup and check math function
  if ( use_local_vth_func ) then
    select case (p_x_dim)
    case (1)
      call setup(this%local_vth_func, trim(local_vth_func_expr), &
                (/'x1'/), ierr)
    case (2)
      call setup(this%local_vth_func, trim(local_vth_func_expr), &
                (/'x1','x2'/), ierr)
    case (3)
      call setup(this%local_vth_func, trim(local_vth_func_expr), &
                (/'x1','x2','x3'/), ierr)
    end select

    if (ierr /= 0) then
      if ( mpi_node() == 0 ) then
        write(0,*) "(*error*) Invalid function supplied : '", &
                    trim(local_vth_func_expr), "'"
      endif
      stop
    endif
  endif

  if( (this%damper_direction .lt. 1) .or. (this%damper_direction .gt. 3) ) then
    write(err_buf__,*) "Damper direction not 1, 2, or 3";call err__("overdense/os-spec-damp.f03",383)
    call abort_program( p_err_invalid )
  endif

  if( this%damper_direction .gt. p_x_dim ) then
    write(err_buf__,*) "Damper direction higher than simulation dimensionality";call err__("overdense/os-spec-damp.f03",388)
    call abort_program( p_err_invalid )
  endif

  if(this%local_vth_cutoff) then
    if( (this%vth_mult_start .le. 0) .or. (this%vth_mult_full .le. 0) ) then
      write(err_buf__,*) "damper vth multipliers must be positive";call err__("overdense/os-spec-damp.f03",394)
      call abort_program( p_err_invalid)
    endif
  endif

  if( (.not. this%local_vth_cutoff) .or. this%use_local_vth_func ) then
    if( (this%gamma_full .lt. 1.0) .or. (this%gamma_start .lt. 1.0) ) then
      write(err_buf__,*) "Damper gammas must be larger than 1.0 when not using";call err__("overdense/os-spec-damp.f03",401)
      write(err_buf__,*) "a local_vth_cutoff or when using a local_vth_func.";call err__("overdense/os-spec-damp.f03",402)
      call abort_program( p_err_invalid )
    endif
  endif

  if( (.not. this%reemit_with_local_temperature) .or. this%use_local_vth_func ) then
    do i = 1, p_p_dim
      if ( this%reemit_temp_start(i) .lt. 0.0 ) then
        write(err_buf__,*) "Damper reemision temperatures must be nonnegative";call err__("overdense/os-spec-damp.f03",410)
        write(err_buf__,*) "when not reemitting with local temperature or";call err__("overdense/os-spec-damp.f03",411)
        write(err_buf__,*) "when using a local_vth_func.";call err__("overdense/os-spec-damp.f03",412)
        call abort_program( p_err_invalid )
      endif
      if ( this%reemit_temp_full(i) .lt. 0.0 ) then
        write(err_buf__,*) "Damper reemision temperatures must be nonnegative";call err__("overdense/os-spec-damp.f03",416)
        write(err_buf__,*) "when not reemitting with local temperature or";call err__("overdense/os-spec-damp.f03",417)
        write(err_buf__,*) "when using a local_vth_func.";call err__("overdense/os-spec-damp.f03",418)
        call abort_program( p_err_invalid )
      endif
    enddo
  endif

  if( this%n_damp .le. 0 ) then
    write(err_buf__,*) "n_damp must be positive or if_damp set to false";call err__("overdense/os-spec-damp.f03",425)
    call abort_program( p_err_invalid )
  endif

end subroutine read_nml_damper
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine init_damper( this, fluid_moments, grid, no_co, dx, g_space )

  type( t_damper ), intent(inout) :: this
  type( t_fluid_moments ), intent(inout) :: fluid_moments
  class( t_grid ), intent(in) :: grid
  class( t_node_conf ), intent(in) :: no_co
  real(p_double), dimension(:), intent(in) :: dx
  type( t_space ), intent(in) :: g_space

  integer :: dir, i
  real(p_double) :: l_xmin, l_xmax

  if( .not. this%if_damp) then
    this%if_damp_local = .false.
    return
  endif

  dir = this%damper_direction

  if( .not. this%custom_input ) then
    ! Set up x and stopping_dist
    ! Damping region is set to be three times the mean free path
    ! x_full is one mean free path from the boundary
    ! Initial stopping distance is set to twice the mean free path
    ! Final stopping distance is determined by the residual
    if ( this%orientation == p_forward ) then
      this%x_start = real( xmax(g_space,dir), p_k_part ) - 3*this%mfp
      this%x_full = real( xmax(g_space,dir), p_k_part ) - this%mfp
    else
      this%x_start = real( xmin(g_space,dir), p_k_part ) + 3*this%mfp
      this%x_full = real( xmin(g_space,dir), p_k_part ) + this%mfp
    endif
    this%stopping_dist_start = 2*this%mfp
    this%stopping_dist_full = calc_stop_full( this%residual ) * 3*this%mfp
  endif

  if( this%reemit_with_local_temperature .or. this%local_vth_cutoff) then
    if ( this%fourth_root_temperature ) then
      call setup(fluid_moments, 'fourth_root_abs_proper_velocity', grid, no_co, dx)
    else
      call setup(fluid_moments, 'abs_proper_velocity', grid, no_co, dx)
    endif
  endif

  if ( this%x_start < xmin(g_space,dir) ) then
    write(err_buf__,*) "Damper x_start must be inside simulation region";call err__("overdense/os-spec-damp.f03",478)
    call abort_program( p_err_invalid )
  endif

  if ( this%x_full > xmax(g_space,dir) ) then
    write(err_buf__,*) "Damper x_full must be inside simulation region";call err__("overdense/os-spec-damp.f03",483)
    call abort_program( p_err_invalid )
  endif

  l_xmin = xmin(g_space, dir) + dx(dir)*grid%my_nx(p_lower, dir)
  l_xmax = xmin(g_space, dir) + dx(dir)*grid%my_nx(p_upper, dir)
  if( (l_xmax .le. this%x_start .and. this%orientation==p_forward) .or. &
      (l_xmin .ge. this%x_start .and. this%orientation==p_backward) ) then
    this%if_damp_local = .false.
  else
    this%if_damp_local = .true.
  endif

  this%x_range = this%x_full - this%x_start

  this%gamma_range = (-1*this%gamma_full) + this%gamma_start

  this%vth_mult_range = (-1*this%vth_mult_full) + this%vth_mult_start

  this%stopping_dist_range = (-1*this%stopping_dist_full) + this%stopping_dist_start

  do i = 1, p_p_dim
    this%reemit_temp_range(i) = (-1*this%reemit_temp_full(i)) + this%reemit_temp_start(i)
  enddo

end subroutine init_damper
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Function designed by Kyle Miller to calculate the stopping distance at x_full based on
! the input of how many particles make it through the absorbing region (residual). Check
! Kyle's dissertation for more information.
!-----------------------------------------------------------------------------------------
function calc_stop_full( residual )

  real(p_double), intent(in) :: residual
  real(p_k_part) :: calc_stop_full

  real(p_double) :: log2res, loglog2res, prob_full

  log2res = log(2*residual)
  loglog2res = log(-log2res)

  prob_full = -(residual/(3*log2res**3))*(6*(1 + log2res + log2res**2)*loglog2res &
                                          -3*(3+log2res)*loglog2res**2 + 2*loglog2res**3 &
                                          + 6*log2res**3*log(-2*residual*log2res))

  calc_stop_full = real( -1/(2*log(prob_full)), p_k_part )

end function calc_stop_full
!-----------------------------------------------------------------------------------------

end module m_spec_damp
