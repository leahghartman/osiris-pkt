#include "os-config.h"
#include "os-preprocess.fpp"

module m_zpulse_flyfoc_cyl_modes

#include "memory/memory.h"

use m_zpulse_flyfoc, only : p_astrl, p_gaussian_new, p_laguerre_new, p_astrl_chromatic, p_astrl_laguerre
use m_zpulse_flyfoc, only : p_astrl_analytic, p_astrl_discrete, p_astrl_discrete_num_pulses_max
use m_zpulse_flyfoc, only : interp_3d
use m_zpulse_mov_wall_cyl_modes
use m_zpulse_std, only : p_forward, p_backward

use m_parameters
use m_math
use m_fparser 
use m_node_conf,     only: t_node_conf
use m_zpulse_mov_wall, only : t_zpulse_mov_wall 
use m_input_file,    only: t_input_file, get_namelist, p_max_nml_section_name
use m_emf_define,    only: t_emf_bound, t_emf
use m_emf_cyl_modes, only: t_emf_cyl_modes
use m_space,         only: t_space, x_bnd, get_x_bnd, xmin, xmax 
use m_vdf_define,    only: t_vdf
use m_grid_define,   only: t_grid
use m_grid_cyl_modes,only: t_grid_cyl_modes
use m_fft_new 

use m_zpulse_std,    only: p_zpulse_bnormal
use m_zpulse_std,    only: p_max_chirp_order, p_polynomial, p_sin2, p_gaussian, p_func
use m_zpulse_std,    only: p_const
use m_zpulse_std,    only: p_hermite_gaussian, p_hermite_gaussian_astigmatic
use m_zpulse_std,    only: p_laguerre_gaussian, p_plane, p_gaussian_asym, t_zpulse, p_forward, p_plane

use m_math


!-----------------------------------------------------------------------------------------
type, extends( t_zpulse_mov_wall_cyl_modes ) :: t_zpulse_flyfoc_cyl_modes
  
    ! whether we should precompute the envelope and use this for injection 
    logical :: if_interp_env

    ! shared flying focus options 
    integer        :: env_type
    real(p_double) :: interp_dx2, interp_dt 
    real(p_double) :: vff 

    ! ! chromatic flying focus params
    ! real(p_double) :: chirp_strength
    ! logical        :: if_use_diffractive_lens

    ! astrl params
    real(p_double) :: z1, z2, z_lead, zr, initial_delay, final_delay, vg, deltak_max
    integer :: integral_size
    integer :: focal_weight_type
    integer :: focal_delay_type
    integer :: focal_phase_type
    real(p_double) :: env_norm
    real(p_double) :: dz0

    type(t_fparser) :: focal_weight_math_func
    type(t_fparser) :: focal_delay_math_func
    type(t_fparser) :: focal_phase_math_func

    ! math funcs for astrl analytic
    type(t_fparser) :: s0_math_func
    type(t_fparser) :: w0_math_func
    type(t_fparser) :: xperp0_math_func
    type(t_fparser) :: a0_math_func
    type(t_fparser) :: pol_math_func
    type(t_fparser) :: phase_math_func
    type(t_fparser) :: phase_derivative_math_func
    
    ! discrete astrl params
    integer :: astrl_discrete_num_pulses
    real(p_double), dimension(:), pointer :: astrl_discrete_focal_points => null()  
    real(p_double), dimension(:), pointer :: astrl_discrete_focal_delays => null()  
    real(p_double), dimension(:), pointer :: astrl_discrete_focal_weights => null()  
    real(p_double), dimension(:), pointer :: astrl_discrete_focal_phases => null()  
    integer :: astrl_discrete_norm_idx 

    ! ! fft params
    ! real(p_double) :: fft_omega_max, fft_omega_delta, fft_k2_max, fft_k2_delta

    complex(p_double), dimension(:,:,:), pointer :: e_env_buffer => null()
  
    ! temporary computation buffer for astrl integral 
    complex(p_double), dimension(:), pointer :: buf_1d => null()

    ! variables used for the prepcomputed envelope for any pulse type
    real(p_double), dimension(4,2) :: e_env_grid_bounds 
    integer, dimension(4) :: e_env_grid_shape
    real(p_double), dimension(4) :: e_env_grid_dx

contains

  procedure :: read_input => read_input_zpulse_flyfoc_cyl_modes
  procedure :: init_env_norm => init_env_norm_zpulse_flyfoc_cyl_modes
  procedure :: init_extra_vars => init_extra_vars_zpulse_flyfoc_cyl_modes
  procedure :: cleanup    => cleanup_zpulse_flyfoc_cyl_modes 

  procedure :: t_duration => t_duration_zpulse_flyfoc_cyl_modes
  procedure :: t_envelope => t_env_zpulse_flyfoc_cyl_modes 

  procedure :: get_env => get_env_zpulse_flyfoc_cyl_modes

  procedure :: get_env_gaussian => get_env_gaussian_zpulse_flyfoc_cyl_modes
  procedure :: get_env_laguerre => get_env_laguerre_zpulse_flyfoc_cyl_modes
  procedure :: get_env_astrl => get_env_astrl_zpulse_flyfoc_cyl_modes
  procedure :: get_env_astrl_laguerre => get_env_astrl_laguerre_zpulse_flyfoc_cyl_modes
  procedure :: get_env_astrl_chromatic => get_env_astrl_chromatic_zpulse_flyfoc_cyl_modes
  procedure :: get_env_astrl_analytic => get_env_astrl_analytic_zpulse_flyfoc_cyl_modes
  procedure :: get_env_astrl_discrete => get_env_astrl_discrete_zpulse_flyfoc_cyl_modes

  procedure :: precompute_envelope => precompute_envelope_zpulse_flyfoc_cyl_modes    

  procedure :: set_fields => set_fields_zpulse_flyfoc_cyl_modes
  procedure :: set_fields_t0 => set_fields_t0_zpulse_flyfoc_cyl_modes

end type t_zpulse_flyfoc_cyl_modes
!-----------------------------------------------------------------------------------------

contains 


!-----------------------------------------------------------------------------------------
subroutine read_input_zpulse_flyfoc_cyl_modes( this, input_file, g_space, bnd_con, periodic, grid, sim_options )

  use m_input_file

  implicit none

  class( t_zpulse_flyfoc_cyl_modes ), intent(inout) :: this
  class( t_input_file ), intent(inout) :: input_file
  type( t_space ), intent(in) :: g_space
  class (t_emf_bound), intent(in) :: bnd_con
  logical, dimension(:), intent(in) :: periodic
  class( t_grid ), intent(in) :: grid
  type( t_options ), intent(in) :: sim_options

  ! local variables
  real(p_double) :: a0
  real(p_double) :: omega0
  real(p_double) :: k0 
  real(p_double) :: phase
  integer        :: pol_type
  real(p_double) :: pol
  character(len=16) :: propagation
  integer         :: direction

  ! Temporal profile
  character(len=16) :: tenv_type
  real(p_double) :: tenv_fwhm
  real(p_double) :: tenv_rise, tenv_flat, tenv_fall
  real(p_double) :: tenv_duration, tenv_range
  character(len = p_max_expr_len) :: tenv_math_func
  real(p_double)  :: xi0

  ! shared flying focus params
  ! real(p_double) :: focal_length 
  character(len=32) :: env_type 
  logical :: if_interp_env
  logical :: if_launch_t0
  real(p_double) :: interp_dx2, interp_dt 
  real(p_double) :: vff 

  ! astrl params
  real(p_double) :: z1, z2, initial_delay, final_delay, vg
  integer :: integral_size
  character(len=16) :: focal_weight_type  
  character(len=16) :: focal_delay_type  
  character(len=16) :: focal_phase_type  
  character(len = p_max_expr_len) :: focal_weight_math_func
  character(len = p_max_expr_len) :: focal_delay_math_func
  character(len = p_max_expr_len) :: focal_phase_math_func

  ! astrl analytic math funcs
  character(len = p_max_expr_len) :: s0_math_func
  character(len = p_max_expr_len) :: w0_math_func
  character(len = p_max_expr_len) :: xperp0_math_func
  character(len = p_max_expr_len) :: a0_math_func
  character(len = p_max_expr_len) :: pol_math_func
  character(len = p_max_expr_len) :: phase_math_func
  character(len = p_max_expr_len) :: phase_derivative_math_func

  ! chromatic flying focus params
  real(p_double) :: deltak_max 
  
  ! discrete astrl params
  integer :: astrl_discrete_num_pulses
  real(p_double), dimension(p_astrl_discrete_num_pulses_max) :: astrl_discrete_focal_points
  real(p_double), dimension(p_astrl_discrete_num_pulses_max) :: astrl_discrete_focal_delays
  real(p_double), dimension(p_astrl_discrete_num_pulses_max) :: astrl_discrete_focal_weights
  real(p_double), dimension(p_astrl_discrete_num_pulses_max) :: astrl_discrete_focal_phases
  logical :: astrl_discrete_use_linspace
  integer :: astrl_discrete_norm_idx

  ! launch parameters
  real(p_double) :: launch_time
  logical        :: if_launch

  ! Beam profile
  ! character(len=16) :: per_type
  real(p_double), dimension(2) :: per_center
  real(p_double), dimension(2) :: per_w0
  integer, dimension(2)        :: per_tem_mode

  ! wall params
  real(p_double)  :: wall_pos  ! initial wall position
  real(p_double)  :: wall_vel  
  integer         :: ncells 
  logical         :: if_launch_from_wall
  logical         :: if_clean_fields
  logical         :: if_apply_divcorr 
  logical         :: if_add_fields

  ! tmp variables
  integer :: N, i 
  real(p_k_fparse), dimension(1) :: fparser_arr 

  namelist /nl_zpulse_flyfoc/ &
        a0, omega0, k0, phase, pol_type, pol, propagation, direction, &
        tenv_type, tenv_fwhm, tenv_rise, tenv_flat, tenv_fall, &
        tenv_duration, tenv_range, & ! tenv_launch_duration, 
        tenv_math_func, xi0, &
        per_center, per_w0, per_tem_mode, &
        env_type, &
        if_interp_env, if_launch_t0, interp_dx2, interp_dt, & 
        vff, vg, z1, z2, initial_delay, final_delay, integral_size, &
        focal_weight_type, focal_weight_math_func, &
        focal_delay_type,  focal_delay_math_func, &
        focal_phase_type,  focal_phase_math_func, &
        s0_math_func, w0_math_func, xperp0_math_func, a0_math_func, pol_math_func, &
        phase_math_func, phase_derivative_math_func, & 
        astrl_discrete_num_pulses, &
        astrl_discrete_focal_delays, astrl_discrete_focal_points, &
        astrl_discrete_focal_weights, astrl_discrete_focal_phases, & 
        astrl_discrete_use_linspace, astrl_discrete_norm_idx, & 
        deltak_max, &
        wall_pos, wall_vel, ncells, if_launch_from_wall, &
        if_clean_fields, if_apply_divcorr, if_add_fields, &
        launch_time, if_launch

  integer :: ierr

  if (disp_out(input_file)) then
    SCR_ROOT(" - Reading zpulse_flyfoc configuration...")
  endif

  if_launch = .true.

  a0             = 1.0_p_double
  omega0         = 10.0_p_double
  k0             = 0.0_p_double
  phase          = 0.0_p_double
  pol_type       = 0
  pol            = 90.0_p_double

  propagation    = "forward"
  direction      = 1

  tenv_type       = "gaussian"
  tenv_fwhm       = -1.0_p_double
  tenv_rise       = 0.0_p_double
  tenv_flat       = 0.0_p_double
  tenv_fall       = 0.0_p_double
  xi0             = 0.0_p_double

  tenv_duration   = 0.0_p_double
  tenv_range      = -huge(1.0_p_double)
  ! tenv_launch_duration      = -huge(1.0_p_double)
  tenv_math_func  = "NO_FUNCTION_SUPPLIED!"

  ! per_type       = "gaussian"
  per_center     = 0.0_p_double
  per_tem_mode   = 0.0_p_double

  ! shared ff options 
  env_type             = "" 
  if_interp_env            = .true.
  if_launch_t0           = .false.
  interp_dx2             = -1.0_p_double
  interp_dt              = -1.0_p_double 
  
  ! continuous astrl
  z1 = 0.0_p_double
  z2 = 0.0_p_double
  initial_delay = 0.0_p_double
  final_delay = 0.0_p_double
  vg = 1.0_p_double  ! group velocity 
  integral_size = 1
  
  ! chromatic astrl
  deltak_max = 0.0_p_double


  focal_weight_type  = "const"
  focal_weight_math_func  = "NO_FUNCTION_SUPPLIED!"
  focal_delay_type  = "const"
  focal_delay_math_func  = "NO_FUNCTION_SUPPLIED!"
  focal_phase_type  = "const"
  focal_phase_math_func  = "NO_FUNCTION_SUPPLIED!"

  ! discrete astrl 
  astrl_discrete_num_pulses = 1 
  astrl_discrete_focal_delays(:) = 0.0_p_double 
  astrl_discrete_focal_points(:) = 0.0_p_double 
  astrl_discrete_focal_weights(:) = 1.0_p_double 
  astrl_discrete_focal_phases(:) = 0.0_p_double 
  astrl_discrete_use_linspace = .true.
  astrl_discrete_norm_idx = 1 

  ! analytic astrl
  s0_math_func = "0"
  w0_math_func = "1"
  xperp0_math_func = "0"
  a0_math_func = "1"
  pol_math_func = "0"
  phase_math_func = "0"
  phase_derivative_math_func = "0"
  
  ! general options 
  launch_time           = 0.0_p_double
  if_launch             = .true.

  ! moving wall params
  wall_vel = 0.0_p_double
  wall_pos = -huge(1.0d0)
  ncells = 4
  if_launch_from_wall = .false. 
  if_clean_fields = .false.
  if_apply_divcorr = .true.
  if_add_fields = .false.

  if( p_x_dim .eq. 1 ) then 
    SCR_ROOT(   "ERROR: zpulse_flyfoc only works for spatial dimensions > 1")
    stop 
  endif

  ! read data from file
  call get_namelist( input_file, "nl_zpulse_flyfoc", ierr )

  if ( ierr /= 0 ) then
    if ( mpi_node() == 0 ) then
      if (ierr < 0) then
        SCR_ROOT( "Error reading zpulse_flyfoc parameters" )
      else
        SCR_ROOT( "Error: zpulse_flyfoc parameters missing" )
      endif
      SCR_ROOT( "aborting...")
    endif
    stop
  endif

  read (input_file%nml_text, nml = nl_zpulse_flyfoc, iostat = ierr)
  if (ierr /= 0) then
    if ( mpi_node() == 0 ) then
      write(0,*)  ""
      write(0,*)  "   Error reading zpulse_flyfoc parameters"
      write(0,*)  "   aborting..."
    endif
    stop
  endif

  this%a0         = a0
  this%omega0     = omega0
  if (k0 /= 0.0) then
     this%k0 = k0
  else
     this%k0 = omega0
  endif
  ! this%wavelength = 2 * pi / omega0 
  this%phase0     = real( phase * pi_180, p_double )

  ! Polarization parameters
  if ((pol_type <-1) .or. (pol_type >1)) then
    if ( mpi_node() == 0 ) then
      write(0,*)  ""
      write(0,*)  "   Error in zpulse_flyfoc parameters"
      write(0,*)  "   pol must be in the range [-1,0,+1]"
      write(0,*)  "   aborting..."
    endif
    stop
  endif
  this%pol_type  = pol_type
  this%pol = real( pol * pi_180, p_double )
  
  ! Propagation direction
  select case ( trim(propagation) )
  case ("forward")
     this%propagation = p_forward
     this%prop_sign = 1.0_p_double
  case ("backward")
     this%propagation = p_backward
     this%prop_sign = -1.0_p_double 
  case default
    if ( mpi_node() == 0 ) then
      write(0,*)  ''
      write(0,*)  '   Error in zpulse_flyfoc parameters'
      write(0,*)  '   propagation must be either "forward" or "backward"'
      write(0,*)  '   aborting...'
    endif
    stop
  end select

  ! Launch direction
  if ((direction < 1) .or. (direction > p_x_dim)) then
    if ( mpi_node() == 0 ) then
      write(0,*)  ""
      write(0,*)  "   Error in zpulse_flyfoc parameters"
      write(0,*)  "   direction must be in the range [1 .. x_dim]"
      write(0,*)  "   aborting..."
    endif
    stop
  endif

  this%direction = direction

    ! this%tenv_launch_duration      = tenv_launch_duration

  this % xi0 = xi0

  this%per_w0 = per_w0
  this%per_tem_mode = per_tem_mode

  this%per_center     = per_center

  this%launch_time    = launch_time
  this%if_launch = if_launch

  this%wall_pos = wall_pos
  this%wall_vel = wall_vel
  this%ncells   = ncells
  this%if_launch_from_wall = if_launch_from_wall
  this%if_clean_fields = if_clean_fields
  this%if_apply_divcorr = if_apply_divcorr
  this%if_add_fields = if_add_fields

  ! interp params
  this%if_interp_env = if_interp_env
  this%if_launch_t0 = if_launch_t0
  this%interp_dx2    = interp_dx2 
  this%interp_dt     = interp_dt 
  
  ! shared ff params
  this%z1 = z1
  this%z2 = z2 
  this%initial_delay = initial_delay
  this%final_delay = final_delay
  this%vff = vff 
  this%vg = vg 

  ! chromatic ff params
  this%deltak_max = deltak_max

  this%zr = 0.5_p_double * omega0 * this%per_w0(1)**2 

! get type of envelope to use 
  select case( env_type ) 

  case( 'gaussian' )
    this%env_type = p_gaussian_new

  case( 'laguerre' )
    this%env_type = p_laguerre_new
    this%per_type = p_laguerre_gaussian

  case( 'astrl' )
    this%env_type = p_astrl

  case( 'astrl_chromatic' )
    this%env_type = p_astrl_chromatic

  case( 'astrl_analytic' )
     this%env_type = p_astrl_analytic 

  case( 'astrl_laguerre' )
    this%env_type = p_astrl_laguerre
    this%per_type = p_laguerre_gaussian

  case( 'astrl_discrete' )
    this%env_type = p_astrl_discrete

  case default  
    SCR_ROOT( 'ERROR in os-zpulse-flyfoc: invalid pulse type requested.')
    SCR_ROOT( 'Options are: gaussian, astrl, astrl_discrete ')
    stop 

  end select 

  if( this%env_type .ne. p_astrl_analytic ) then 
     
     select case ( trim(tenv_type) )
        
     case ("polynomial")
        this%tenv_type = p_polynomial

        if ( tenv_fwhm > 0.0_p_double ) then
           this%tenv_rise       = tenv_fwhm/2.0_p_double
           this%tenv_flat       = 0.0_p_double
           this%tenv_fall       = tenv_fwhm/2.0_p_double
        else
           this%tenv_rise       = tenv_rise
           this%tenv_flat       = tenv_flat
           this%tenv_fall       = tenv_fall
        endif

     case ("sin2")
        this%tenv_type = p_sin2

        if ( tenv_fwhm > 0.0_p_double ) then
           this%tenv_rise       = tenv_fwhm/2.0_p_double
           this%tenv_flat       = 0.0_p_double
           this%tenv_fall       = tenv_fwhm/2.0_p_double
        else
           this%tenv_rise       = tenv_rise
           this%tenv_flat       = tenv_flat
           this%tenv_fall       = tenv_fall
        endif

     case ("gaussian")
        this%tenv_type = p_gaussian

        if ( tenv_fwhm > 0.0_p_double  ) then
           this%tenv_duration   = tenv_fwhm/sqrt( 2.0_p_double * log(2.0_p_double))
        else
           this%tenv_duration   = tenv_duration
        endif
        this%tenv_range      = tenv_range

        if ( tenv_range == -huge(1.0_p_double) ) then
           if ( mpi_node() == 0 ) then
              write(0,*)  ''
              write(0,*)  '   Error in zpulse_flyfoc parameters'
              write(0,*)  '   When using "gaussian" temporal envelope, '
              write(0,*)  '   tenv_range must also be defined'
              write(0,*)  '   aborting...'
           endif
           stop
        endif

     case ("math")
        this%tenv_type = p_func
        this%tenv_range      = tenv_range

        if ( tenv_range == -huge(1.0_p_double) ) then
           if ( mpi_node() == 0 ) then
              write(0,*)  ''
              write(0,*)  '   Error in zpulse_flyfoc parameters'
              write(0,*)  '   When using "math" temporal envelopes, '
              write(0,*)  '   tenv_range must also be defined'
              write(0,*)  '   aborting...'
           endif
           stop
        endif

        call setup(this%tenv_math_func, trim(tenv_math_func), (/'t'/), ierr)
        if (ierr /= 0) then
           if ( mpi_node() == 0 ) then
              write(0,*)  ''
              write(0,*)  '   Error in zpulse_flyfoc parameters'
              write(0,*)  '   Supplied tenv_math_func failed to compile'
              write(0,*)  '   aborting...'
           endif
           stop
        endif

     case default
        if ( mpi_node() == 0 ) then
           write(0,*)  ''
           write(0,*)  '   Error in zpulse_flyfoc parameters'
           write(0,*)  '   tenv_type must be either "polynomial", "gaussian",'
           write(0,*)  '   "sin2", "math" or "hermite"'
           write(0,*)  '   aborting...'
        endif
        stop
     end select

  end if
  
  ! todo: restructure code so that axis fields are corrected when launching at t=0.
  ! currently, the code will crash because axis correction is only implemented for moving wall mode.
  if( this%if_launch_t0 ) then
    this%if_correct_axis_fields = .false.
  endif

  if( this%if_interp_env .and. &
        ( ( this%interp_dx2 .le. 0 ) .or. ( this%interp_dt .le. 0 ) ) ) then 
          SCR_ROOT( '   ERROR: when using if_interp_env = .true., interp_dx2 and interp_dt must be set')
          stop
  endif 

  
  if( this%env_type .eq. p_gaussian_new .or. this%env_type .eq. p_laguerre_new ) then 

    this%tenv_launch_duration = this%initial_delay + this%final_delay

  else if( this%env_type .eq. p_astrl .or. this%env_type .eq. p_astrl_laguerre &
                                      .or. this%env_type .eq. p_astrl_chromatic  ) then 

    this%integral_size = integral_size

    if( z2 <= z1 ) then
      SCR_ROOT( "ERROR in os-zpulse-flyfoc: z2 > z1 is required when using tunable_env_gaussian")
      stop
    endif

    if( integral_size .le. 0 ) then 
      SCR_ROOT( "ERROR in os-zpulse-flyfoc: integral size > 0 is required")
      stop
    endif

    this%tenv_launch_duration = abs( ( this%vg/this%vff - 1 ) * (z2 - z1) ) &
      + this%initial_delay + this%final_delay 

    select case( focal_weight_type ) 

    case( "const" ) 
      this%focal_weight_type = p_const

    case( "math" )
      this%focal_weight_type = p_func 

      call setup(this%focal_weight_math_func, &
                  trim(focal_weight_math_func), (/'z0'/), ierr)
      
      if (ierr /= 0) then
        if ( mpi_node() == 0 ) then
          write(0,*)  ''
          write(0,*)  '   Error in zpulse_flyfoc parameters'
          write(0,*)  '   Supplied focal_weight_math_func failed to compile'
          write(0,*)  '   aborting...'
        endif
        stop
      endif

    case default 
      SCR_ROOT( '   ERROR: invalid focal_weight_type: ', focal_weight_type )
      SCR_ROOT( '   Options are: const, math ')
      stop 
    
    end select 


    select case( focal_delay_type ) 

    case( "const" ) 
      this%focal_delay_type = p_const

    case( "math" )
      this%focal_delay_type = p_func 

      call setup(this%focal_delay_math_func, &
                  trim(focal_delay_math_func), (/'z0'/), ierr)
      
      if (ierr /= 0) then
        if ( mpi_node() == 0 ) then
          write(0,*)  ''
          write(0,*)  '   Error in zpulse_flyfoc parameters'
          write(0,*)  '   Supplied focal_delay_math_func failed to compile'
          write(0,*)  '   aborting...'
        endif
        stop
      endif

    case default 
      SCR_ROOT( '   ERROR: invalid focal_delay_type: ', focal_delay_type )
      SCR_ROOT( '   Options are: const, math ')
      stop 
    
    end select 


    select case( focal_phase_type ) 

    case( "const" ) 
      this%focal_phase_type = p_const

    case( "math" )
      this%focal_phase_type = p_func 

      call setup(this%focal_phase_math_func, &
                  trim(focal_phase_math_func), (/'z0'/), ierr)
      
      if (ierr /= 0) then
        if ( mpi_node() == 0 ) then
          write(0,*)  ''
          write(0,*)  '   Error in zpulse_flyfoc parameters'
          write(0,*)  '   Supplied focal_phase_math_func failed to compile'
          write(0,*)  '   aborting...'
        endif
        stop
      endif

    case default 
      SCR_ROOT( '   ERROR: invalid focal_phase_type: ', focal_phase_type )
      SCR_ROOT( '   Options are: const, math ')
      stop 
    
    end select 
  
  else if( this%env_type .eq. p_astrl_discrete ) then 

    N = astrl_discrete_num_pulses

    if( N > p_astrl_discrete_num_pulses_max ) then
      SCR_ROOT( '   ERROR: requested astrl_discrete_num_pulses = ', astrl_discrete_num_pulses )
      SCR_ROOT( '   Max num pulses = ', p_astrl_discrete_num_pulses_max )
      stop
    endif

    this%astrl_discrete_num_pulses = N
    this%astrl_discrete_norm_idx = astrl_discrete_norm_idx

    call alloc( this%astrl_discrete_focal_points, (/ N /))
    call alloc( this%astrl_discrete_focal_delays, (/ N /))
    call alloc( this%astrl_discrete_focal_weights, (/ N /))
    call alloc( this%astrl_discrete_focal_phases, (/ N /))

    ! astrl discrete -- focal points 
    if( astrl_discrete_use_linspace ) then 
      do i=1, N 
        this%astrl_discrete_focal_points(i) = this%z1 + (this%z2 - this%z1 ) / (N-1)* (i-1)  
      enddo
    else
      this%astrl_discrete_focal_points(1:N) = astrl_discrete_focal_points(1:N)
    endif

    ! astrl discrete --  focal delay 
    if( focal_delay_type .eq. 'math' ) then  

      call setup(this%focal_delay_math_func, &
                  trim(focal_delay_math_func), (/'z0'/), ierr)
      
      if (ierr /= 0) then
        if ( mpi_node() == 0 ) then
          write(0,*)  ''
          write(0,*)  '   Error in zpulse_flyfoc parameters'
          write(0,*)  '   Supplied focal_delay_math_func failed to compile'
          write(0,*)  '   aborting...'
        endif
        stop
      endif

      do i=1, N 
        fparser_arr(1) = this%astrl_discrete_focal_points(i)
        this%astrl_discrete_focal_delays(i) = eval( this%focal_delay_math_func, fparser_arr )
      enddo

    else
      this%astrl_discrete_focal_delays(1:N) = astrl_discrete_focal_delays(1:N) 
    endif

    ! astrl discrete -- focal weight 
    if( focal_weight_type .eq. 'math' ) then  

      call setup(this%focal_weight_math_func, &
                  trim(focal_weight_math_func), (/'z0'/), ierr)
      
      if (ierr /= 0) then
        if ( mpi_node() == 0 ) then
          write(0,*)  ''
          write(0,*)  '   Error in zpulse_flyfoc parameters'
          write(0,*)  '   Supplied focal_weight_math_func failed to compile'
          write(0,*)  '   aborting...'
        endif
        stop
      endif

      do i=1, N 
        fparser_arr(1) = this%astrl_discrete_focal_points(i)
        this%astrl_discrete_focal_weights(i) = eval( this%focal_weight_math_func, fparser_arr )
      enddo

    else
      this%astrl_discrete_focal_weights(1:N) = astrl_discrete_focal_weights(1:N) 
    endif

    ! astrl discrete --  focal phase 
    if( focal_phase_type .eq. 'math' ) then  

      call setup(this%focal_phase_math_func, &
                  trim(focal_phase_math_func), (/'z0'/), ierr)
      
      if (ierr /= 0) then
        if ( mpi_node() == 0 ) then
          write(0,*)  ''
          write(0,*)  '   Error in zpulse_flyfoc parameters'
          write(0,*)  '   Supplied focal_phase_math_func failed to compile'
          write(0,*)  '   aborting...'
        endif
        stop
      endif

      do i=1, N 
        fparser_arr(1) = this%astrl_discrete_focal_points(i)
        this%astrl_discrete_focal_phases(i) = eval( this%focal_phase_math_func, fparser_arr )
      enddo

    else
      this%astrl_discrete_focal_phases(1:N) = astrl_discrete_focal_phases(1:N) 
    endif

    this%tenv_launch_duration = (maxval(this%astrl_discrete_focal_delays) &
                              - minval(this%astrl_discrete_focal_delays)) &
                                + this%initial_delay + this%final_delay 

    ! astrl analytic 
  else if( this%env_type .eq. p_astrl_analytic ) then

     call setup(this%s0_math_func, trim(s0_math_func), (/'xi'/), ierr)

     if (ierr /= 0) then
        if ( mpi_node() == 0 ) then
           write(0,*)  ''
           write(0,*)  '   Error in zpulse_flyfoc parameters'
           write(0,*)  '   Supplied s0_math_func failed to compile'
           write(0,*)  '   aborting...'
        endif
        stop
     endif

     call setup(this%w0_math_func, trim(w0_math_func), (/'xi'/), ierr)

     if (ierr /= 0) then
        if ( mpi_node() == 0 ) then
           write(0,*)  ''
           write(0,*)  '   Error in zpulse_flyfoc parameters'
           write(0,*)  '   Supplied w0_math_func failed to compile'
           write(0,*)  '   aborting...'
        endif
        stop
     endif

     call setup(this%xperp0_math_func, trim(xperp0_math_func), (/'xi'/), ierr)

     if (ierr /= 0) then
        if ( mpi_node() == 0 ) then
           write(0,*)  ''
           write(0,*)  '   Error in zpulse_flyfoc parameters'
           write(0,*)  '   Supplied xperp0_math_func failed to compile'
           write(0,*)  '   aborting...'
        endif
        stop
     endif

     call setup(this%a0_math_func, trim(a0_math_func), (/'xi'/), ierr)

     if (ierr /= 0) then
        if ( mpi_node() == 0 ) then
           write(0,*)  ''
           write(0,*)  '   Error in zpulse_flyfoc parameters'
           write(0,*)  '   Supplied a0_math_func failed to compile'
           write(0,*)  '   aborting...'
        endif
        stop
     endif

     call setup(this%pol_math_func, trim(pol_math_func), (/'xi'/), ierr)

     if (ierr /= 0) then
        if ( mpi_node() == 0 ) then
           write(0,*)  ''
           write(0,*)  '   Error in zpulse_flyfoc parameters'
           write(0,*)  '   Supplied pol_math_func failed to compile'
           write(0,*)  '   aborting...'
        endif
        stop
     endif
    
     call setup(this%phase_math_func, trim(phase_math_func), (/'xi'/), ierr)

     if (ierr /= 0) then
        if ( mpi_node() == 0 ) then
           write(0,*)  ''
           write(0,*)  '   Error in zpulse_flyfoc parameters'
           write(0,*)  '   Supplied phase_math_func failed to compile'
           write(0,*)  '   aborting...'
        endif
        stop
     endif

     call setup(this%phase_derivative_math_func, trim(phase_derivative_math_func), (/'xi'/), ierr)

     if (ierr /= 0) then
        if ( mpi_node() == 0 ) then
           write(0,*)  ''
           write(0,*)  '   Error in zpulse_flyfoc parameters'
           write(0,*)  '   Supplied phase_derivative_math_func failed to compile'
           write(0,*)  '   aborting...'
        endif
        stop
     endif
    
  endif 


  ! Boost pulse
  this%if_boost = .false.
  
  if (sim_options%gamma > 1.0_p_double) then

    this % gamma          = sim_options%gamma
    this % gamma_times_beta = sqrt(this%gamma**2-1.0_p_double)

    if( this % propagation .eq. p_forward ) then
      this % boost_per_fac = 1.0_p_double / ( this%gamma + this%gamma_times_beta )
    else
      this % boost_per_fac = this%gamma + this%gamma_times_beta
    endif 
    SCR_ROOT('   boosting zpulse with gamma = ', this%gamma )

    this%if_boost = .true.
  else
    this % boost_per_fac = 1.0_p_double
  endif

  call this%init_extra_vars() 

  call this%init_env_norm( grid, g_space ) 

end subroutine read_input_zpulse_flyfoc_cyl_modes
!-----------------------------------------------------------------------------------------



!-----------------------------------------------------------------------------------------
subroutine init_extra_vars_zpulse_flyfoc_cyl_modes( this )

  implicit none 

  class( t_zpulse_flyfoc_cyl_modes ), intent(inout) :: this

  select case( this%env_type )

  case( p_astrl, p_astrl_laguerre, p_astrl_chromatic ) 

    call alloc( this%buf_1d, (/ this%integral_size /) )
    this%buf_1d(:) = 0

    if( this%focal_delay_type .eq. p_const ) then 
      if( this%vff > 0 .and. this%vff < this%vg ) then 
        this%z_lead = this%z1
      else
        this%z_lead = this%z2 
      endif       
    endif 

    this%dz0 = ( this%z2 - this%z1 ) / this%integral_size

  end select 

end subroutine init_extra_vars_zpulse_flyfoc_cyl_modes
!-----------------------------------------------------------------------------------------



!-----------------------------------------------------------------------------------------
subroutine init_env_norm_zpulse_flyfoc_cyl_modes( this, grid, g_space )

  implicit none 

  class( t_zpulse_flyfoc_cyl_modes ), intent(inout) :: this
  class( t_grid ), intent(in) :: grid
  type( t_space ),     intent(in) :: g_space

  real(p_double) :: default_env_norm 
  complex(p_double) :: tmp_env_norm 
  real(p_double ) :: t, z0, x2, delay
  real(p_k_fparse), dimension(1) :: fparser_arr 
  real(p_double) :: gxmin2, gxmax2, dx2
  integer :: nx2, i2
  real(p_double), dimension(:), pointer :: norm_buf => null()
  logical :: if_interp_env

  default_env_norm = this%boost_per_fac * this%a0

  x2 = 0.0_p_double

  gxmin2 = g_space%x_bnd_initial(p_lower,2)
  gxmax2 = g_space%x_bnd_initial(p_upper,2)
  nx2 = grid%g_nx(2)
  dx2 = ( gxmax2 - gxmin2 ) / nx2

  select case( this%env_type )

  case( p_astrl, p_astrl_laguerre, p_astrl_chromatic ) 
    
    ! set to 1 for calculation 
    this%env_norm = 1.0_p_double

    ! normalize envelope to the central value of z 
    z0 = ( this%z1 + this%z2 ) / 2 

    ! compute the time at which we will hit this value of z 
    select case( this%focal_delay_type )

    case( p_const ) 

      delay = ( z0 - this%z_lead ) * ( this%vg/this%vff - 1 ) 

    case( p_func )

      fparser_arr(1) = z0
      delay = eval( this%focal_delay_math_func, fparser_arr )

    case default

      SCR_ROOT( "ERROR: Should not reach here in init_env_norm_zpulse_flyfoc_cyl_modes")
      stop
      delay = 0

    end select

    ! must add lon_center so that the z0-focused pulse component is peaked at (z0,t(z0)) 
    t = ( z0 + delay ) + this%lon_center() 

    call alloc( norm_buf, (/nx2/))

    ! disable envelope interpolation for env norm initialization.
    if_interp_env = this%if_interp_env 
    this % if_interp_env = .false.

    do i2 = 1, nx2 
      x2 = gxmin2 + (i2-1) * dx2 
      call this%get_env( t, z0, x2, tmp_env_norm )
      norm_buf(i2) = abs(tmp_env_norm)
    enddo

    tmp_env_norm = maxval( norm_buf )

    ! restore if_interp_env to original value 
    this % if_interp_env = if_interp_env 
          
    call freemem( norm_buf )

    this%env_norm = default_env_norm / abs( tmp_env_norm )


  case( p_astrl_discrete )

    this%env_norm = 1.0_p_double

    z0 = this%astrl_discrete_focal_points( this%astrl_discrete_norm_idx )
    delay = this%astrl_discrete_focal_delays( this%astrl_discrete_norm_idx )
    
    t = z0 + delay + this%lon_center() 

    call this%get_env_astrl_discrete( t, z0, x2, tmp_env_norm )

    this%env_norm = default_env_norm / abs( tmp_env_norm )

  ! everything else 
  case default 
    
    this%env_norm = default_env_norm

  end select 

end subroutine init_env_norm_zpulse_flyfoc_cyl_modes
!-----------------------------------------------------------------------------------------



!-----------------------------------------------------------------------------------------
! Total duration of the laser pulse
!-----------------------------------------------------------------------------------------
function t_duration_zpulse_flyfoc_cyl_modes( this ) result( t_duration )

  implicit none

  class( t_zpulse_flyfoc_cyl_modes ), intent(in) :: this
  real(p_double) :: t_duration

  t_duration = this%tenv_launch_duration

  if( this % if_boost ) then
    t_duration = this % gamma * t_duration!  - this % gamma_times_beta * this % wall_pos
  endif

end function t_duration_zpulse_flyfoc_cyl_modes
!----------------------------------------------------------------------------------------



!-----------------------------------------------------------------------------------------
function t_env_zpulse_flyfoc_cyl_modes( this, t_sim, z_sim ) result(t_env)

  implicit none

  class( t_zpulse_flyfoc_cyl_modes ), intent(inout) :: this ! intent must be inout because of
                                                ! the eval( f_parser ) function
  real(p_double), intent(in) :: t_sim
  real(p_double), intent(in) :: z_sim 

  ! lab coordinates 
  real(p_double) :: z, t

  real(p_double) :: t_env

  ! local variables

  real(p_k_fparse), dimension(1) ::t_dbl

  t = t_sim
  z = z_sim 

  ! ignore pulse broadening 
  if( this%propagation .eq. p_forward ) then 
    t = t - z + this % xi0
  else
    t = t + z - this % xi0
  endif

  select case ( this%tenv_type )

  case( p_const )
    t_env = 1.0_p_double 

  case (p_polynomial)
    t_env = fenv_poly( t, this%tenv_rise, this%tenv_flat, this%tenv_fall )

  case (p_sin2)
    t_env = fenv_sin2( t, this%tenv_rise, this%tenv_flat, this%tenv_fall )

  case (p_gaussian)
    t_env = exp( - 2*((t - 0.5*this % tenv_range)/this%tenv_duration)**2 )

  case (p_func)

    if ( t > this % tenv_range ) then
      t_env = 0
    else
      t_dbl(1) = t
      t_env = eval( this%lon_math_func, t_dbl )
    endif

  case default
    t_env = 0

  end select

end function t_env_zpulse_flyfoc_cyl_modes
!-----------------------------------------------------------------------------------------



!-----------------------------------------------------------------------------------------
subroutine cleanup_zpulse_flyfoc_cyl_modes( this )

  implicit none 

  class( t_zpulse_flyfoc_cyl_modes ), intent(inout) :: this

  if ( associated(this%e_env_buffer) ) call freemem( this%e_env_buffer )

  if ( associated(this%buf_1d)) call freemem( this%buf_1d )

  if( associated( this%astrl_discrete_focal_delays )) &
      call freemem( this%astrl_discrete_focal_delays)
  if( associated( this%astrl_discrete_focal_points )) &
      call freemem( this%astrl_discrete_focal_points)
  if( associated( this%astrl_discrete_focal_weights )) &
      call freemem( this%astrl_discrete_focal_weights)
  if( associated( this%astrl_discrete_focal_phases )) &
      call freemem( this%astrl_discrete_focal_phases)

end subroutine cleanup_zpulse_flyfoc_cyl_modes
!-----------------------------------------------------------------------------------------



!-----------------------------------------------------------------------------------------
function gaussian_beam_q3d( z, r, omega, zr, z0 ) result( ret )

  implicit none 

  real(p_double), intent(in) :: z, r, omega, zr, z0
  complex(p_double) :: ret
  complex(p_double) :: q 

  q = cmplx( z - z0, zr, p_double )
  ret = exp( cmplx(0,-omega, p_double)  * ( r**2 / (2 * q ) ) ) * (zr/q) 

end function gaussian_beam_q3d
!-----------------------------------------------------------------------------------------



!-----------------------------------------------------------------------------------------
function laguerre_beam_q3d( z, r, omega, zr, w0, z0, ell, p ) result( ret )

  implicit none 

  real(p_double), intent(in) :: z, r, omega, zr, w0, z0
  integer, intent(in) :: ell, p
  real(p_double) :: zcent, w_of_z, Rinverse_of_z, guoy_shift
  complex(p_double) :: ret

  zcent  = z - z0 
  w_of_z = w0 * sqrt( 1 + (zcent / zr )**2 )
  Rinverse_of_z = zcent / ( zcent**2 + zr**2 )
  guoy_shift = ( 2*p + abs(ell) + 1 ) * atan2( zcent, z0 )

  ret = (w0 / w_of_z ) &
            * ( r * sqrt(2.0_p_double) / w_of_z ) ** abs(ell) & 
            * exp( - ( r / w_of_z )**2 ) &
            * laguerre( p, ell, 2 * ( r / w_of_z ) **2 ) & 
            * exp( cmplx( 0, 1) * ( - 0.5_p_double * omega * r**2 * Rinverse_of_z + guoy_shift ) ) 

end function laguerre_beam_q3d
!-----------------------------------------------------------------------------------------



!-----------------------------------------------------------------------------------------
subroutine precompute_envelope_zpulse_flyfoc_cyl_modes( this, e, g_space, nx_p_min, dt ) 

  implicit none 

  integer, parameter :: rank = 2

  class( t_zpulse_flyfoc_cyl_modes ), intent(inout) :: this
  type( t_vdf ), intent(in) :: e
  type( t_space ), intent(in) :: g_space
  integer, dimension(:), intent(in) :: nx_p_min
  real(p_double), intent(in) :: dt

  real(p_double), dimension(2,rank) :: g_x_range
  real(p_double), dimension(rank) :: ldx

  integer :: nt, nz, nx2
  integer :: i0, i1, i2

  ! bounds of the spacetime region that must be computed 
  real(p_double) :: tmin, tmax, zmin, zmax, x2min, x2max, delta_t, zwall 
  real(p_double) :: x2, z, t
  real(p_double) :: interp_dx1

  complex(p_double) :: env

  ! progress messages
  integer :: progress_counter, num_progress_intervals

  ! type(t_fft_manager) :: fft_manager 

  SCR_ROOT( "Computing pulse envelope..." )

  ! calculate spatial data 
  call get_x_bnd( g_space, g_x_range )

  ldx(1:2) = real( e%dx_(1:2), p_double )

  if( .not. this%if_launch_from_wall ) then 

    zmin = this%wall_pos - 2 * ldx(1) 
    zmax = this%wall_pos + (this%ncells + 2) * ldx(1)

  else

    if( this%propagation .eq. p_forward ) then 

      zwall = real( g_x_range(p_lower,1) )
      zmin = zwall - 2 * ldx(1)
      zmax = zwall + ( this%ncells + 2 ) * ldx(1) 

    else

      zwall = real( g_x_range(p_upper,1) )
      zmin = zwall - ( this%ncells + 2 ) * ldx(1)
      zmax = zwall + 2 * ldx(1) 

    endif

  endif

  x2min = real( g_x_range(p_lower,2) ) + ( lbound( e%f2, 3 ) + nx_p_min(2) - 3 ) * ldx(2) 
  x2max = real( g_x_range(p_lower,2) ) + ( ubound( e%f2, 3 ) + nx_p_min(2) + 3 ) * ldx(2) 

  select case( this%interpolation )
  case( p_linear, p_cubic )
    ! Nothing to change here
  case( p_quadratic, p_quartic )
    zmin = zmin + ldx(1)*0.5_p_double
    zmax = zmax + ldx(1)*0.5_p_double
  end select

  if( .not. this%if_boost ) then 
    tmin = 0 - this%interp_dt  
    tmax = this%tenv_launch_duration + this%interp_dt 
    interp_dx1 = ldx(1) 
  else 
    tmin = this%gamma_times_beta * this%wall_pos
    tmax = this%tenv_launch_duration + ldx(1) * this%gamma_times_beta
    zmin = this%gamma * zmin
    zmax = this%gamma * zmax 
    interp_dx1 = this%interp_dt
  endif

  ! temporarily disable envelope field interpolation so that we can get 
  ! the exact field values within this function call -- unset at end of this routine
  this%if_interp_env = .false.



  delta_t = this%interp_dt

  ! dimensions of the spatial parts of the buffer 
  nt = ceiling( (tmax - tmin )/ this%interp_dt ) + 1 
  nz = ceiling( (zmax - zmin) / interp_dx1 ) + 1 
  nx2 = ceiling( (x2max - x2min) / this%interp_dx2 ) + 1 

  ! print *, mpi_node(), 'x2min, x2max, nx2 ', zmin, zmax, interp_dx1, nz

  this%e_env_grid_shape(1:3) = (/ nt, nz, nx2 /)

  this%e_env_grid_bounds(1,1:2) = (/ tmin, tmax /)
  this%e_env_grid_bounds(2,1:2) = (/ zmin, zmax /)
  this%e_env_grid_bounds(3,1:2) = (/ x2min, x2max /)

  this%e_env_grid_dx(1) = this%interp_dt
  this%e_env_grid_dx(2) = interp_dx1
  this%e_env_grid_dx(3) = this%interp_dx2

  call alloc( this%e_env_buffer, (/ nt, nz, nx2 /) )
  this%e_env_buffer(:,:,:) = 0 

  progress_counter = 0
  num_progress_intervals = 10

  do i0 = 1, nt 

    if( i0 * 1.0_p_double / nt >= ( progress_counter * 1.0_p_double / num_progress_intervals ) ) then 
      SCR_ROOT( i0, " / ", nt, "..." )
      progress_counter = progress_counter + 1
    endif

    t = tmin + ( i0 - 1 ) * this%e_env_grid_dx(1) 

    do i1 = 1, nz
 
      z = zmin + ( i1 - 1 ) * this%e_env_grid_dx(2) 

      do i2 = 1, nx2  
        x2 = x2min + (i2 - 1 ) * this%e_env_grid_dx(3)

        call this%get_env( t, z, x2, env )

        this%e_env_buffer(i0, i1, i2 ) =  env 

      enddo 

    enddo

  enddo

  SCR_ROOT( "Successfully computed envelope." )

  ! revert if_interp_env -- now, the fields can be interpolated. 
  this%if_interp_env = .true.

end subroutine precompute_envelope_zpulse_flyfoc_cyl_modes
!-----------------------------------------------------------------------------------------



!-----------------------------------------------------------------------------------------
subroutine get_env_zpulse_flyfoc_cyl_modes( this, t, z, x2, env )

  implicit none 

  class( t_zpulse_flyfoc_cyl_modes ), intent(inout) :: this
  real(p_double), intent(in) :: t, z, x2 
  complex(p_double), intent(out) :: env 

  real(p_double), dimension(3) :: x, x1, dx
  integer, dimension(3) :: ix1

  if( this%if_interp_env ) then

    dx(1:3) = this%e_env_grid_dx(1:3)
    
    ix1(1) = floor( (t - this%e_env_grid_bounds(1,1)) / this%e_env_grid_dx(1) )
    x1(1) = dx(1) * ix1(1) + this%e_env_grid_bounds(1,1)
    x(1) = t
    
    ix1(2) = floor( (z - this%e_env_grid_bounds(2,1)) / this%e_env_grid_dx(2) )
    x1(2) = dx(2) * ix1(2) + this%e_env_grid_bounds(2,1)
    x(2) = z

    ix1(3) = floor( (x2 - this%e_env_grid_bounds(3,1)) / this%e_env_grid_dx(3) )
    x1(3) = dx(3) * ix1(3) + this%e_env_grid_bounds(3,1)
    x(3) = x2

    env = interp_3d( x, ix1, x1, dx, this%e_env_buffer )

  else

    select case( this%env_type )
  
    case( p_gaussian_new ) 
      call this%get_env_gaussian( t, z, x2, env )

    case( p_laguerre_new ) 
      call this%get_env_laguerre( t, z, x2, env )

    case( p_astrl )
      call this%get_env_astrl( t, z, x2, env ) 

    case( p_astrl_laguerre )
      call this%get_env_astrl_laguerre( t, z, x2, env ) 

    case( p_astrl_chromatic )
      call this%get_env_astrl_chromatic( t, z, x2, env ) 

   case( p_astrl_analytic )
      call this%get_env_astrl_analytic( t, z, x2, env )
      
    case( p_astrl_discrete ) 
      call this%get_env_astrl_discrete( t, z, x2, env ) 

    end select

  endif

end subroutine get_env_zpulse_flyfoc_cyl_modes
!-----------------------------------------------------------------------------------------



!-----------------------------------------------------------------------------------------
subroutine get_env_gaussian_zpulse_flyfoc_cyl_modes( this, t, z, x2, env )

  implicit none 

  class( t_zpulse_flyfoc_cyl_modes ), intent(inout) :: this
  real(p_double), intent(in) :: t, z, x2 
  complex(p_double), intent(out) :: env 

  env = gaussian_beam_q3d( z, x2, this%omega0, this%zr, this%z1 ) * this%t_envelope( t, z ) 

  env = env * this%env_norm 

end subroutine get_env_gaussian_zpulse_flyfoc_cyl_modes
!-----------------------------------------------------------------------------------------



!-----------------------------------------------------------------------------------------
subroutine get_env_laguerre_zpulse_flyfoc_cyl_modes( this, t, z, x2, env )

  implicit none 

  class( t_zpulse_flyfoc_cyl_modes ), intent(inout) :: this
  real(p_double), intent(in) :: t, z, x2 
  complex(p_double), intent(out) :: env 

  env = laguerre_beam_q3d( z, x2, this%omega0, this%zr, this%per_w0(1), this%z1, &
                           this%per_tem_mode(2), this%per_tem_mode(1) ) * this%t_envelope( t, z ) 

  env = env * this%env_norm 

end subroutine get_env_laguerre_zpulse_flyfoc_cyl_modes
!-----------------------------------------------------------------------------------------



!-----------------------------------------------------------------------------------------
subroutine get_env_astrl_zpulse_flyfoc_cyl_modes( this, t, z, x2, env )

  implicit none 

  class( t_zpulse_flyfoc_cyl_modes ), intent(inout) :: this
  real(p_double), intent(in) :: t, z, x2 
  complex(p_double), intent(out) :: env 

  integer :: i3, nz0 
  real(p_double) :: z0, delay 
  real(p_k_fparse), dimension(1) :: fparser_arr 

  nz0 = this%integral_size

  do i3 = 1, nz0
    
    z0 = this%z1 + (i3 - 1) * this%dz0

    delay = ( z0 - this%z_lead ) * ( this%vg/this%vff - 1 ) 

    this%buf_1d(i3) = gaussian_beam_q3d( z, x2, this%omega0, this%zr, z0 ) &
                      * this%t_envelope( t - delay, z ) 

    if( this%focal_weight_type .eq. p_func ) then 
      fparser_arr(1) = z0 
      this%buf_1d(i3) = this%buf_1d(i3) * eval( this%focal_weight_math_func, fparser_arr )
    endif 

    if( this%focal_phase_type .eq. p_func ) then 
      fparser_arr(1) = z0 
      this%buf_1d(i3) = this%buf_1d(i3) * exp( cmplx( 0, eval( this%focal_phase_math_func, fparser_arr ), p_double ) )
    endif 

  enddo

  env = this%dz0 * sum( this%buf_1d ) 

  env = env * this%env_norm 

end subroutine get_env_astrl_zpulse_flyfoc_cyl_modes
!-----------------------------------------------------------------------------------------



!-----------------------------------------------------------------------------------------
subroutine get_env_astrl_laguerre_zpulse_flyfoc_cyl_modes( this, t, z, x2, env )

  implicit none 

  class( t_zpulse_flyfoc_cyl_modes ), intent(inout) :: this
  real(p_double), intent(in) :: t, z, x2 
  complex(p_double), intent(out) :: env 

  integer :: i3, nz0 
  real(p_double) :: z0, delay 
  real(p_k_fparse), dimension(1) :: fparser_arr 

  nz0 = this%integral_size

  do i3 = 1, nz0
    
    z0 = this%z1 + (i3 - 1) * this%dz0

    delay = ( z0 - this%z_lead ) * ( this%vg/this%vff - 1 ) 

    this%buf_1d(i3) = laguerre_beam_q3d( z, x2, this%omega0, this%zr, this%per_w0(1), z0, &
                           this%per_tem_mode(2), this%per_tem_mode(1) ) &
                           * this%t_envelope( t - delay, z ) 

    if( this%focal_weight_type .eq. p_func ) then 
      fparser_arr(1) = z0 
      this%buf_1d(i3) = this%buf_1d(i3) * eval( this%focal_weight_math_func, fparser_arr )
    endif 

    if( this%focal_phase_type .eq. p_func ) then 
      fparser_arr(1) = z0 
      this%buf_1d(i3) = this%buf_1d(i3) * exp( cmplx( 0, eval( this%focal_phase_math_func, fparser_arr ), p_double ) )
    endif 

  enddo

  env = this%dz0 * sum( this%buf_1d ) 

  env = env * this%env_norm 

end subroutine get_env_astrl_laguerre_zpulse_flyfoc_cyl_modes
!-----------------------------------------------------------------------------------------



!-----------------------------------------------------------------------------------------
subroutine get_env_astrl_chromatic_zpulse_flyfoc_cyl_modes( this, t, z, x2, env )

  class( t_zpulse_flyfoc_cyl_modes ), intent(inout) :: this
  real(p_double), intent(in) :: t, z, x2 
  complex(p_double), intent(out) :: env 

  integer :: i3, nz0 
  real(p_double) :: z0, deltak, k, zr, delay
  real(p_k_fparse), dimension(1) :: fparser_arr 

  nz0 = this%integral_size

  do i3 = 1, nz0
    
    z0 = this%z1 + (i3 - 1) * this%dz0

    ! frequency shift due to chirp -- goes from -deltak_max to deltak_max
    deltak = this%deltak_max * 2.0_p_double * ( i3 - nz0/2.0_p_double) / nz0 

    ! effective k for this pulse 
    k = this%omega0 + deltak 

    ! use normal ASTRL flying focus delay profile 
    delay = ( z0 - this%z_lead ) * ( this%vg/this%vff - 1 ) 

    ! rayleigh range is changed because of frequency variation 
    zr = ( 1 + deltak/this%omega0 )**2 * this%zr

    ! last factor is chirp -- blue light leads 
    this%buf_1d(i3) = gaussian_beam_q3d( z, x2, k, zr, z0 ) &
                   * this%t_envelope( t - delay, z ) * exp( cmplx(0, deltak * (t-z), p_double ) )

    if( this%focal_weight_type .eq. p_func ) then 
      fparser_arr(1) = z0 
      this%buf_1d(i3) = this%buf_1d(i3) * eval( this%focal_weight_math_func, fparser_arr )
    endif 

    if( this%focal_phase_type .eq. p_func ) then 
      fparser_arr(1) = z0 
      this%buf_1d(i3) = this%buf_1d(i3) * exp( cmplx( 0, eval( this%focal_phase_math_func, fparser_arr ), p_double ) )
    endif 

  enddo

  env = this%dz0 * sum( this%buf_1d ) 

  env = env * this%env_norm 

end subroutine get_env_astrl_chromatic_zpulse_flyfoc_cyl_modes
!-----------------------------------------------------------------------------------------



!-----------------------------------------------------------------------------------------
subroutine get_env_astrl_analytic_zpulse_flyfoc_cyl_modes( this, t, z, x2, env )

  class( t_zpulse_flyfoc_cyl_modes ), intent(inout) :: this
  real(p_double), intent(in) :: t, z, x2 
  complex(p_double), intent(out) :: env 

  real(p_double) :: xi, s0, w0, xperp0, a0, zr, omega, phase
  real(p_k_fparse), dimension(1) :: fparser_arr 
  
  ! xi accounts for direction of pulse propagation. 
  xi = z - this%prop_sign * t 

  ! compute all paraxial beam properties for this particular xi slice 
  fparser_arr(1) = xi

  s0 = eval( this%s0_math_func, fparser_arr ) 
  w0 = eval( this%w0_math_func, fparser_arr ) 
  xperp0 = eval( this%xperp0_math_func, fparser_arr ) 
  a0 = eval( this%a0_math_func, fparser_arr )
  phase = eval( this%phase_math_func, fparser_arr ) 
  omega = this%omega0 + eval( this%phase_derivative_math_func, fparser_arr )

  zr = 0.5_p_double * this%omega0 * w0**2 

  env = a0 * gaussian_beam_q3d( z, x2, omega, zr, s0 ) * exp( cmplx( 0, -1 ) * phase ) 
  
end subroutine get_env_astrl_analytic_zpulse_flyfoc_cyl_modes
!-----------------------------------------------------------------------------------------



!-----------------------------------------------------------------------------------------
subroutine get_env_astrl_discrete_zpulse_flyfoc_cyl_modes( this, t, z, x2, env )

  implicit none 

  class( t_zpulse_flyfoc_cyl_modes ), intent(inout) :: this
  real(p_double), intent(in) :: t, z, x2 
  complex(p_double), intent(out) :: env 

  real(p_double) :: z0, delay, weight, phase
  integer :: i 

  env = 0.0_p_double

  do i = 1, this%astrl_discrete_num_pulses

    z0 = this%astrl_discrete_focal_points(i)
    delay = this%astrl_discrete_focal_delays(i)
    weight = this%astrl_discrete_focal_weights(i)
    phase = this%astrl_discrete_focal_phases(i)

    env = env + weight * gaussian_beam_q3d( z, x2, this%omega0, this%zr, z0 ) &
                        * this%t_envelope( t - delay, z ) &
                        * exp( cmplx( 0, phase, p_double ) )

  enddo

  env = env * this%env_norm 

end subroutine get_env_astrl_discrete_zpulse_flyfoc_cyl_modes
!-----------------------------------------------------------------------------------------




!-----------------------------------------------------------------------------------------
! B = zhat \times E ---> B_phi = prop_sign * E_r, B_r = - prop_sign * E_phi 
!-----------------------------------------------------------------------------------------
subroutine set_fields_zpulse_flyfoc_cyl_modes( this, e_cyl_m, b_cyl_m, g_space, &
                                               nx_p_min, g_nx, t_sim, dt, no_co )

  implicit none 

  class( t_zpulse_flyfoc_cyl_modes ), intent(inout) :: this
  type( t_cyl_modes ), intent(inout) :: e_cyl_m, b_cyl_m
  type( t_space ), intent(in) :: g_space
  integer, intent(in), dimension(:) :: nx_p_min, g_nx
  real(p_double), intent(in) :: t_sim
  real(p_double), intent(in) :: dt
  class( t_node_conf ), intent(in) :: no_co

  real(p_double) :: zmin, zmax, rmin, z, z_2, r, r_2, dr_2, dz_2, t, t_2
  real(p_double) :: z_sim, z_sim_2
  real(p_double) :: gpos, prop_sign, amp

  complex(p_double) :: env, field_0_0, field_0_2, field_2_0, field_2_2

  integer :: ipos, gipos 
  
  integer :: lbnd, rbnd
  logical :: lbnd_in_domain, rbnd_in_domain, domain_enclosed, lbnd_right_of_domain

  integer :: i1, i2 
  integer :: i1start, i1finish, clean_end 

  integer, parameter :: rank = 2
  real(p_double), dimension(2,rank) :: g_x_range
  
  real(p_double) :: cos_pol, sin_pol
  real(p_double), dimension(rank) :: ldx

  integer :: n, ell

  type(t_vdf), pointer :: e_re0

  if( this%if_launch .and. this%if_launch_t0 ) then 
    call this % set_fields_t0(e_cyl_m, b_cyl_m, g_space, nx_p_min, g_nx, t_sim, dt, no_co)
    this%if_launch = .false.
    return
  endif

  ! Turn off antenna if pulse has ended

  if ( t_sim > this % t_duration() ) then
    SCR_ROOT( "duration exceeded. " )
    SCR_ROOT( "this%t_duration(): ", this%t_duration() )
    SCR_ROOT( "t: ", t )
    this%if_launch = .false.
    return
  endif

  e_re0 => e_cyl_m%pf_re(0)

  if( this%if_interp_env .and. .not. associated(this%e_env_buffer ) ) then 
    call this % precompute_envelope( e_re0, g_space, nx_p_min, dt )
  endif

  amp = this%omega0

  ldx(1) = e_re0%dx(1)
  ldx(2) = e_re0%dx(2)

  ! get zmin, zmax, rmin in absolute coordinates 
  call get_x_bnd( g_space, g_x_range )
  zmin = real( g_x_range(p_lower,1), p_double ) 
  zmax = real( g_x_range(p_upper,1), p_double ) 
  rmin = real( g_x_range(p_lower,2), p_double ) 

  cos_pol = cos( this%pol )
  sin_pol = sin( this%pol )

  dz_2 = ldx(1)/2.0_p_double
  dr_2 = ldx(2)/2.0_p_double

  select case( this%interpolation )
  case( p_linear, p_cubic )
    ! Nothing to change here
  case( p_quadratic, p_quartic )
    zmin = zmin + dz_2
    zmax = zmax + dz_2
    ! rmin = rmin + dr_2
  end select

  if ( this%propagation == p_forward ) then
    prop_sign = 1.0_p_double
  else
    prop_sign = -1.0_p_double
  endif


  ! start near domain boundary 
  if( this%if_launch_from_wall ) then 
    
    if( this%propagation == p_forward ) then
      gipos = 1  
    else
      gipos = floor( ( zmax - zmin ) / ldx(1) ) - 1 
    endif
  
  ! else compute the wall position in global coordinates, stored as gpos 
  else 
    gpos = this%wall_pos + this%wall_vel * t_sim
    gipos = int((gpos - zmin)/ldx(1)) + 1

    if( (gpos - ldx(1) < zmin) .or. (gpos + (this%ncells + 1) * ldx(1) > zmax ) ) then 
      SCR_ROOT( "WARNING: zpulse_mov_wall hit the boundary before pulse injection" )
      SCR_ROOT( "was complete. Stopping injection.")
      this%if_launch = .false.
      return
    endif

  endif

  ipos = gipos - nx_p_min(1) + 1

  ! switch to local coordinates 
  zmin = zmin + real(nx_p_min(1), p_double)*ldx(1)
  rmin = rmin + real(nx_p_min(2), p_double)*ldx(2)

  ! left and right index bounds of the calculation in local coordinates 
  if( this%propagation == p_forward ) then 
    lbnd = ipos 
    rbnd = ipos + this%ncells
  else
    lbnd = ipos - this%ncells 
    rbnd = ipos 
  endif

  lbnd_in_domain  = (lbnd >= 1) .and. (lbnd  <= e_re0%nx_(1))
  lbnd_right_of_domain  = (lbnd  >= e_re0%nx_(1))
  rbnd_in_domain = (rbnd >= 1) .and. (rbnd <= e_re0%nx_(1))
  domain_enclosed  = (lbnd < 1)  .and. (rbnd >  e_re0%nx_(1))

  lbnd = max( lbnd, 1 )
  rbnd = min( rbnd, e_re0%nx_(1) )

  i1start = lbnd
  i1finish = rbnd

  ell = this % per_tem_mode(2)

  if( lbnd_in_domain .or. rbnd_in_domain .or. domain_enclosed ) then 

    if( .not. this%if_add_fields ) then 

      e_cyl_m%pf_re(0)%f2(:,i1start:i1finish,:) = 0
      b_cyl_m%pf_re(0)%f2(:,i1start:i1finish,:) = 0

      do n = 1, ubound( e_cyl_m%pf_re, 1 )

        e_cyl_m%pf_re(n)%f2(:,i1start:i1finish,:) = 0
        e_cyl_m%pf_im(n)%f2(:,i1start:i1finish,:) = 0
        b_cyl_m%pf_re(n)%f2(:,i1start:i1finish,:) = 0
        b_cyl_m%pf_im(n)%f2(:,i1start:i1finish,:) = 0

      enddo

    endif


    do i1 = i1start, i1finish 

      z_sim = zmin + real( i1 - 1, p_double ) * ldx(1) 
      z_sim_2 = z_sim + 0.5_p_double * ldx(1)

      ! convert to lab frame coords 
      ! note: z_sim_2 is z_sim shifted by a half cell
      ! in the lab frame, t_2 = t and z_2 = z_sim_2 
      ! in the boosted frame, (z_2,t_2) are the boosted coordinates of (t_sim, z_sim_2) 
      ! so t_2 should always we used when using z_2. 
      ! t_2 and z_2 are not the best variable names since they are different from the 
      ! usual convention (_2 usually means + dx/2), but i couldn't think of anything better.

      if( this % if_boost ) then 
        call this % apply_lorentz_transformation( z_sim, t_sim, z, t )
        call this % apply_lorentz_transformation( z_sim_2, t_sim, z_2, t_2 )
      else 
        z = z_sim
        z_2 = z_sim_2
        t = t_sim 
        t_2 = t_sim
      endif 

      do i2 = e_re0%lbound(3), e_re0%ubound(3)
      ! do i2 = 2, e_re%nx_(2)

        r = rmin + real(i2 - 1, p_double ) * ldx(2)
        r_2 = r + 0.5_p_double * ldx(2)

        ! z, r_2
        call this%get_env( t, z, r_2, env )
        field_0_2 = amp * env * this%fast_phase( t, z )

        ! z, r 
        call this%get_env( t, z, r, env ) 
        field_0_0 = amp * env * this%fast_phase( t, z )

        ! z_2, r 
        call this%get_env( t_2, z_2, r, env )
        field_2_0 = amp * env * this%fast_phase( t_2, z_2 )

        ! z_2, r_2
        call this%get_env( t_2, z_2, r_2, env )
        field_2_2 = amp * env * this%fast_phase( t_2, z_2 ) 

        if( ( this%env_type .ne. p_laguerre_new &
              .and. this%env_type .ne. p_astrl_laguerre )  &
              .or. ell .eq. 0 ) then 
    
          select case( this%pol_type ) 

          ! lin pol
          case( 0 )

            ! e_r: z, r_2
            e_cyl_m%pf_re(1)%f2(2,i1,i2) = + real( field_0_2 ) * cos( this%pol )

            e_cyl_m%pf_im(1)%f2(2,i1,i2) = + real( field_0_2 ) * sin( this%pol )

            ! e_phi: z, r
            e_cyl_m%pf_re(1)%f2(3,i1,i2) = + real( field_0_0 ) * sin( this%pol )

            e_cyl_m%pf_im(1)%f2(3,i1,i2) = - real( field_0_0 ) * cos( this%pol )

            ! b_r: z_2, r
            b_cyl_m%pf_re(1)%f2(2,i1,i2) = - real( field_2_0 ) * sin( this%pol ) * prop_sign

            b_cyl_m%pf_im(1)%f2(2,i1,i2) = + real( field_2_0 ) * cos( this%pol ) * prop_sign

            ! b_phi: z_2, r_2
            b_cyl_m%pf_re(1)%f2(3,i1,i2) = + real( field_2_2 ) * cos( this%pol ) * prop_sign

            b_cyl_m%pf_im(1)%f2(3,i1,i2) = + real( field_2_2 ) * sin( this%pol ) * prop_sign

          ! handle circ pol 
          ! same as general ell profile below 
          case( 1 )
            
            ! ! e_r 
            ! call this % per_envelope( z, r_2, per_env_amp, per_env_phase )
            ! total_phase = per_env_phase + fast_phase

            ! emf%e_cyl_m%pf_re(1)%f2(2,i1,i2) = + lenv * per_env_amp * cos( total_phase )
            ! emf%e_cyl_m%pf_im(1)%f2(2,i1,i2) = - lenv * per_env_amp * sin( total_phase )

            ! ! e_phi 
            ! call this % per_envelope( z, r, per_env_amp, per_env_phase )
            ! total_phase = per_env_phase + fast_phase
            
            ! emf%e_cyl_m%pf_re(1)%f2(3,i1,i2) = - lenv * per_env_amp * sin( total_phase )
            ! emf%e_cyl_m%pf_im(1)%f2(3,i1,i2) = - lenv * per_env_amp * cos( total_phase )

            ! ! b_r
            ! call this % per_envelope( z_2, r, per_env_amp, per_env_phase )
            ! total_phase = per_env_phase + fast_phase_2

            ! emf%b_cyl_m%pf_re(1)%f2(2,i1,i2) = + lenv_2 * per_env_amp * sin( total_phase ) * prop_sign
            ! emf%b_cyl_m%pf_im(1)%f2(2,i1,i2) = + lenv_2 * per_env_amp * cos( total_phase ) * prop_sign

            ! ! b_phi 
            ! call this % per_envelope( z_2, r_2, per_env_amp, per_env_phase )
            ! total_phase = per_env_phase + fast_phase_2
            
            ! emf%b_cyl_m%pf_re(1)%f2(3,i1,i2) = + lenv_2 * per_env_amp * cos( total_phase ) * prop_sign
            ! emf%b_cyl_m%pf_im(1)%f2(3,i1,i2) = - lenv_2 * per_env_amp * sin( total_phase ) * prop_sign

          ! same as general ell except with sign switches for imaginary parts
          case( -1 ) 

            ! ! e_r 
            ! call this % per_envelope( z, r_2, per_env_amp, per_env_phase )
            ! total_phase = per_env_phase + fast_phase

            ! emf%e_cyl_m%pf_re(1)%f2(2,i1,i2) = + lenv * per_env_amp * cos( total_phase )
            ! emf%e_cyl_m%pf_im(1)%f2(2,i1,i2) = + lenv * per_env_amp * sin( total_phase )

            ! ! e_phi 
            ! call this % per_envelope( z, r, per_env_amp, per_env_phase )
            ! total_phase = per_env_phase + fast_phase
            
            ! emf%e_cyl_m%pf_re(1)%f2(3,i1,i2) = + lenv * per_env_amp * sin( total_phase )      
            ! emf%e_cyl_m%pf_im(1)%f2(3,i1,i2) = - lenv * per_env_amp * cos( total_phase )      

            ! ! b_r
            ! call this % per_envelope( z_2, r, per_env_amp, per_env_phase )
            ! total_phase = per_env_phase + fast_phase_2

            ! emf%b_cyl_m%pf_re(1)%f2(2,i1,i2) = - lenv_2 * per_env_amp * sin( total_phase ) * prop_sign           
            ! emf%b_cyl_m%pf_im(1)%f2(2,i1,i2) = + lenv_2 * per_env_amp * cos( total_phase ) * prop_sign           

            ! ! b_phi 
            ! call this % per_envelope( z_2, r_2, per_env_amp, per_env_phase )
            ! total_phase = per_env_phase + fast_phase_2
            
            ! emf%b_cyl_m%pf_re(1)%f2(3,i1,i2) = + lenv_2 * per_env_amp * cos( total_phase ) * prop_sign
            ! emf%b_cyl_m%pf_im(1)%f2(3,i1,i2) = + lenv_2 * per_env_amp * sin( total_phase ) * prop_sign
            
          end select 
        
        ! else: is a laguerre pulse with ell /= 0 
        else

          select case( this%pol_type ) 

          ! lin pol 
          case( 0 )

            ! e_r: z, r_2
            e_cyl_m%pf_re(ell-1)%f2(2,i1,i2) = + 0.5 * real( field_0_2 * exp( cmplx( 0, this%pol, p_k_fld ) ) )

            if( ell-1 > 0 ) then
              e_cyl_m%pf_im(ell-1)%f2(2,i1,i2) = - 0.5 * aimag( field_0_2 * exp( cmplx( 0, this%pol, p_k_fld ) ) )
            endif

            e_cyl_m%pf_re(ell+1)%f2(2,i1,i2) = + 0.5 * real( field_0_2 * exp( cmplx( 0, -this%pol, p_k_fld ) ) )

            e_cyl_m%pf_im(ell+1)%f2(2,i1,i2) = - 0.5 * aimag( field_0_2 * exp( cmplx( 0, -this%pol, p_k_fld ) ) )


            ! e_phi: z, r
            e_cyl_m%pf_re(ell-1)%f2(3,i1,i2) = + 0.5 * aimag( field_0_0 * exp( cmplx( 0, this%pol, p_k_fld ) ) )

            if( ell-1 > 0 ) then
              e_cyl_m%pf_im(ell-1)%f2(3,i1,i2) = + 0.5 * real( field_0_0 * exp( cmplx( 0, this%pol, p_k_fld ) ) )
            endif

            e_cyl_m%pf_re(ell+1)%f2(3,i1,i2) = - 0.5 * aimag( field_0_0 * exp( cmplx( 0, -this%pol, p_k_fld ) ) )

            e_cyl_m%pf_im(ell+1)%f2(3,i1,i2) = - 0.5 * real( field_0_0 * exp( cmplx( 0, -this%pol, p_k_fld ) ) )


            ! b_r: z_2, r
            b_cyl_m%pf_re(ell-1)%f2(2,i1,i2) = - 0.5 * aimag( field_2_0 * exp( cmplx( 0, this%pol, p_k_fld ) ) ) * prop_sign

            if( ell-1 > 0 ) then
              b_cyl_m%pf_im(ell-1)%f2(2,i1,i2) = - 0.5 * real( field_2_0 * exp( cmplx( 0, this%pol, p_k_fld ) ) ) * prop_sign
            endif

            b_cyl_m%pf_re(ell+1)%f2(2,i1,i2) = + 0.5 * aimag( field_2_0 * exp( cmplx( 0, -this%pol, p_k_fld ) ) ) * prop_sign

            b_cyl_m%pf_im(ell+1)%f2(2,i1,i2) = + 0.5 * real( field_2_0 * exp( cmplx( 0, -this%pol, p_k_fld ) ) ) * prop_sign


            ! b_phi: z_2, r_2
            b_cyl_m%pf_re(ell-1)%f2(3,i1,i2) = + 0.5 * real( field_2_2 * exp( cmplx( 0, this%pol, p_k_fld ) ) ) * prop_sign

            if( ell-1 > 0 ) then
              b_cyl_m%pf_im(ell-1)%f2(3,i1,i2) = - 0.5 * aimag( field_2_2 * exp( cmplx( 0, this%pol, p_k_fld ) ) ) * prop_sign
            endif

            b_cyl_m%pf_re(ell+1)%f2(3,i1,i2) = + 0.5 * real( field_2_2 * exp( cmplx( 0, -this%pol, p_k_fld ) ) ) * prop_sign

            b_cyl_m%pf_im(ell+1)%f2(3,i1,i2) = - 0.5 * aimag( field_2_2 * exp( cmplx( 0, -this%pol, p_k_fld ) ) ) * prop_sign


          ! handle circ pol 
          case( 1 )
            
            ! ! e_r 
            ! call this % per_envelope_laguerre( z, r_2, per_env_amp, per_env_phase )
            ! total_phase = per_env_phase + fast_phase

            ! emf%e_cyl_m%pf_re(ell+1)%f2(2,i1,i2) = lenv * per_env_amp * cos( total_phase )

            ! emf%e_cyl_m%pf_im(ell+1)%f2(2,i1,i2) = - lenv * per_env_amp * sin( total_phase )

            ! ! e_phi 
            ! call this % per_envelope_laguerre( z, r, per_env_amp, per_env_phase )
            ! total_phase = per_env_phase + fast_phase
            
            ! emf%e_cyl_m%pf_re(ell+1)%f2(3,i1,i2) = - lenv * per_env_amp * sin( total_phase )

            ! emf%e_cyl_m%pf_im(ell+1)%f2(3,i1,i2) = - lenv * per_env_amp * cos( total_phase )

            ! ! b_r
            ! call this % per_envelope_laguerre( z_2, r, per_env_amp, per_env_phase )
            ! total_phase = per_env_phase + fast_phase_2

            ! emf%b_cyl_m%pf_re(ell+1)%f2(2,i1,i2) = + lenv_2 * per_env_amp * sin( total_phase ) * prop_sign

            ! emf%b_cyl_m%pf_im(ell+1)%f2(2,i1,i2) = + lenv_2 * per_env_amp * cos( total_phase ) * prop_sign

            ! ! b_phi 
            ! call this % per_envelope_laguerre( z_2, r_2, per_env_amp, per_env_phase )
            ! total_phase = per_env_phase + fast_phase_2
            
            ! emf%b_cyl_m%pf_re(ell+1)%f2(3,i1,i2) = + lenv_2 * per_env_amp * cos( total_phase ) * prop_sign

            ! emf%b_cyl_m%pf_im(ell+1)%f2(3,i1,i2) = - lenv_2 * per_env_amp * sin( total_phase ) * prop_sign


          case( -1 ) 

            ! ! e_r 
            ! call this % per_envelope_laguerre( z, r_2, per_env_amp, per_env_phase )
            ! total_phase = per_env_phase + fast_phase

            ! emf%e_cyl_m%pf_re(ell-1)%f2(2,i1,i2) = + lenv * per_env_amp * cos( total_phase + this%pol )
            
            ! if( ell-1 > 0 ) then 
            !   emf%e_cyl_m%pf_im(ell-1)%f2(2,i1,i2) = - lenv * per_env_amp * sin( total_phase + this%pol )
            ! endif          

            ! ! e_phi 
            ! call this % per_envelope_laguerre( z, r, per_env_amp, per_env_phase )
            ! total_phase = per_env_phase + fast_phase
            
            ! emf%e_cyl_m%pf_re(ell-1)%f2(3,i1,i2) = + lenv * per_env_amp * sin( total_phase + this%pol )
            
            ! if( ell-1 > 0 ) then 
            !   emf%e_cyl_m%pf_im(ell-1)%f2(3,i1,i2) = + lenv * per_env_amp * cos( total_phase + this%pol )
            ! endif          

            ! ! b_r
            ! call this % per_envelope_laguerre( z_2, r, per_env_amp, per_env_phase )
            ! total_phase = per_env_phase + fast_phase_2

            ! emf%b_cyl_m%pf_re(ell-1)%f2(2,i1,i2) = - lenv_2 * per_env_amp * sin( total_phase + this%pol ) * prop_sign          

            ! if( ell-1 > 0 ) then 
            !   emf%b_cyl_m%pf_im(ell-1)%f2(2,i1,i2) = - lenv_2 * per_env_amp * cos( total_phase + this%pol ) * prop_sign
            ! endif        

            ! ! b_phi 
            ! call this % per_envelope_laguerre( z_2, r_2, per_env_amp, per_env_phase )
            ! total_phase = per_env_phase + fast_phase_2
            
            ! emf%b_cyl_m%pf_re(ell-1)%f2(3,i1,i2) = + lenv_2 * per_env_amp * cos( total_phase + this%pol ) * prop_sign
            
            ! if( ell-1 > 0 ) then 
            !   emf%b_cyl_m%pf_im(ell-1)%f2(3,i1,i2) = - lenv_2 * per_env_amp * sin( total_phase + this%pol ) * prop_sign
            ! endif    

          end select 

        endif
      
      enddo ! i2 
    
    enddo ! i1 

  endif

  if( this%if_clean_fields ) then 
  
    if( lbnd_in_domain ) then 
      clean_end = i1start - 1
    else if( lbnd_right_of_domain ) then 
      return
    else
      return
    endif

    if( ell .eq. 0 ) then

      e_cyl_m%pf_re(1)%f2(:,1:clean_end,:) = 0
      e_cyl_m%pf_im(1)%f2(:,1:clean_end,:) = 0
      b_cyl_m%pf_re(1)%f2(:,1:clean_end,:) = 0
      b_cyl_m%pf_im(1)%f2(:,1:clean_end,:) = 0

    ! no imaginary parts for mode 0
    else if( ell .eq. 1 ) then

      e_cyl_m%pf_re(0)%f2(:,:clean_end,:) = 0
      b_cyl_m%pf_re(0)%f2(:,:clean_end,:) = 0

      e_cyl_m%pf_re(2)%f2(:,:clean_end,:) = 0
      e_cyl_m%pf_im(2)%f2(:,:clean_end,:) = 0
      b_cyl_m%pf_re(2)%f2(:,:clean_end,:) = 0
      b_cyl_m%pf_im(2)%f2(:,:clean_end,:) = 0

    else

      e_cyl_m%pf_re(ell-1)%f2(:,:clean_end,:) = 0
      e_cyl_m%pf_im(ell-1)%f2(:,:clean_end,:) = 0
      b_cyl_m%pf_re(ell-1)%f2(:,:clean_end,:) = 0
      b_cyl_m%pf_im(ell-1)%f2(:,:clean_end,:) = 0

      e_cyl_m%pf_re(ell+1)%f2(:,:clean_end,:) = 0
      e_cyl_m%pf_im(ell+1)%f2(:,:clean_end,:) = 0
      b_cyl_m%pf_re(ell+1)%f2(:,:clean_end,:) = 0
      b_cyl_m%pf_im(ell+1)%f2(:,:clean_end,:) = 0

    endif

  endif

end subroutine set_fields_zpulse_flyfoc_cyl_modes
!-----------------------------------------------------------------------------------------




!-----------------------------------------------------------------------------------------
! set fields at t=0 instead of at a fixed wall position 
! B = zhat \times E ---> B_phi = prop_sign * E_r, B_r = - prop_sign * E_phi 
!-----------------------------------------------------------------------------------------
subroutine set_fields_t0_zpulse_flyfoc_cyl_modes( this, e_cyl_m, b_cyl_m, g_space, &
                                                  nx_p_min, g_nx, t_sim, dt, no_co )

  implicit none 

  class( t_zpulse_flyfoc_cyl_modes ), intent(inout) :: this
  type( t_cyl_modes ), intent(inout) :: e_cyl_m, b_cyl_m
  type( t_space ), intent(in) :: g_space
  integer, intent(in), dimension(:) :: nx_p_min, g_nx
  real(p_double), intent(in) :: t_sim
  real(p_double), intent(in) :: dt
  class( t_node_conf ), intent(in) :: no_co

  real(p_double) :: zmin, zmax, rmin, z, z_2, r, r_2, dr_2, dz_2, t, t_2
  real(p_double) :: z_sim, z_sim_2
  real(p_double) :: prop_sign, amp 
  
  complex(p_double) :: env, field_0_0, field_0_2, field_2_0, field_2_2
  
  integer :: i1, i2 

  integer, parameter :: rank = 2
  real(p_double), dimension(2,rank) :: g_x_range
  
  real(p_double) :: cos_pol, sin_pol
  real(p_double), dimension(rank) :: ldx

  integer :: ell 

  type(t_vdf), pointer :: e_re0

  SCR_ROOT( 'Injecting flying focus fields...' )

  e_re0 => e_cyl_m%pf_re(0)
  ldx(1) = e_re0%dx(1)
  ldx(2) = e_re0%dx(2)

  ! get zmin, zmax, rmin in absolute coordinates 
  call get_x_bnd( g_space, g_x_range )
  zmin = real( g_x_range(p_lower,1), p_double ) 
  zmax = real( g_x_range(p_upper,1), p_double ) 
  rmin = real( g_x_range(p_lower,2), p_double ) 

  cos_pol = cos( this%pol )
  sin_pol = sin( this%pol )

  amp = this%omega0
  
  dz_2 = ldx(1)/2.0_p_double
  dr_2 = ldx(2)/2.0_p_double

  select case( this%interpolation )
  case( p_linear, p_cubic )
    ! Nothing to change here
  case( p_quadratic, p_quartic )
    zmin = zmin + dz_2
    zmax = zmax + dz_2
    ! rmin = rmin + dr_2
  end select

  if ( this%propagation == p_forward ) then
    prop_sign = 1.0_p_double
  else
    prop_sign = -1.0_p_double
  endif

  ! switch to local coordinates 
  zmin = zmin + real(nx_p_min(1), p_double)*ldx(1)
  rmin = rmin + real(nx_p_min(2), p_double)*ldx(2)

  ell = this % per_tem_mode(2)

  ! inject pulse 
  do i1 = e_re0%lbound(2), e_re0%ubound(2)

    z_sim = zmin + real( i1 - 1, p_double ) * ldx(1) 
    z_sim_2 = z_sim + 0.5_p_double * ldx(1)

    ! convert to lab frame coords 
    ! note: z_sim_2 is z_sim shifted by a half cell
    ! in the lab frame, t_2 = t and z_2 = z_sim_2 
    ! in the boosted frame, (z_2,t_2) are the boosted coordinates of (t_sim, z_sim_2) 
    ! so t_2 should always we used when using z_2. 
    ! t_2 and z_2 are not the best variable names since they are different from the 
    ! usual convention (_2 usually means + dx/2), but i couldn't think of anything better.

    if( this % if_boost ) then 
      call this % apply_lorentz_transformation( z_sim, t_sim, z, t )
      call this % apply_lorentz_transformation( z_sim_2, t_sim, z_2, t_2 )
    else 
      z = z_sim
      z_2 = z_sim_2
      t = t_sim 
      t_2 = t_sim
   endif

   do i2 = e_re0%lbound(3), e_re0%ubound(3)

      r = rmin + real(i2 - 1, p_double ) * ldx(2)
      r_2 = r + 0.5_p_double * ldx(2)

      ! z, r_2
      call this%get_env( t, z, r_2, env )
      field_0_2 = amp * env * this%fast_phase( t, z ) 

      ! z, r 
      call this%get_env( t, z, r, env ) 
      field_0_0 = amp * env * this%fast_phase( t, z ) 

      ! z_2, r 
      call this%get_env( t_2, z_2, r, env )
      field_2_0 = amp * env * this%fast_phase( t_2, z_2 )  

      ! z_2, r_2
      call this%get_env( t_2, z_2, r_2, env )
      field_2_2 = amp * env * this%fast_phase( t_2, z_2 ) 

      if( ( this%env_type .ne. p_laguerre_new &
            .and. this%env_type .ne. p_astrl_laguerre )  &
            .or. ell .eq. 0 ) then 
  
        select case( this%pol_type ) 

        ! lin pol
        case( 0 )

          ! e_r: z, r_2
          e_cyl_m%pf_re(1)%f2(2,i1,i2) = + real( field_0_2 ) * cos( this%pol )

          e_cyl_m%pf_im(1)%f2(2,i1,i2) = + real( field_0_2 ) * sin( this%pol )

          ! e_phi: z, r
          e_cyl_m%pf_re(1)%f2(3,i1,i2) = + real( field_0_0 ) * sin( this%pol )

          e_cyl_m%pf_im(1)%f2(3,i1,i2) = - real( field_0_0 ) * cos( this%pol )

          ! b_r: z_2, r
          b_cyl_m%pf_re(1)%f2(2,i1,i2) = - real( field_2_0 ) * sin( this%pol ) * prop_sign

          b_cyl_m%pf_im(1)%f2(2,i1,i2) = + real( field_2_0 ) * cos( this%pol ) * prop_sign

          ! b_phi: z_2, r_2
          b_cyl_m%pf_re(1)%f2(3,i1,i2) = + real( field_2_2 ) * cos( this%pol ) * prop_sign

          b_cyl_m%pf_im(1)%f2(3,i1,i2) = + real( field_2_2 ) * sin( this%pol ) * prop_sign

        ! handle circ pol 
        ! same as general ell profile below 
        case( 1 )
          
          ! ! e_r 
          ! call this % per_envelope( z, r_2, per_env_amp, per_env_phase )
          ! total_phase = per_env_phase + fast_phase

          ! emf%e_cyl_m%pf_re(1)%f2(2,i1,i2) = + lenv * per_env_amp * cos( total_phase )
          ! emf%e_cyl_m%pf_im(1)%f2(2,i1,i2) = - lenv * per_env_amp * sin( total_phase )

          ! ! e_phi 
          ! call this % per_envelope( z, r, per_env_amp, per_env_phase )
          ! total_phase = per_env_phase + fast_phase
          
          ! emf%e_cyl_m%pf_re(1)%f2(3,i1,i2) = - lenv * per_env_amp * sin( total_phase )
          ! emf%e_cyl_m%pf_im(1)%f2(3,i1,i2) = - lenv * per_env_amp * cos( total_phase )

          ! ! b_r
          ! call this % per_envelope( z_2, r, per_env_amp, per_env_phase )
          ! total_phase = per_env_phase + fast_phase_2

          ! emf%b_cyl_m%pf_re(1)%f2(2,i1,i2) = + lenv_2 * per_env_amp * sin( total_phase ) * prop_sign
          ! emf%b_cyl_m%pf_im(1)%f2(2,i1,i2) = + lenv_2 * per_env_amp * cos( total_phase ) * prop_sign

          ! ! b_phi 
          ! call this % per_envelope( z_2, r_2, per_env_amp, per_env_phase )
          ! total_phase = per_env_phase + fast_phase_2
          
          ! emf%b_cyl_m%pf_re(1)%f2(3,i1,i2) = + lenv_2 * per_env_amp * cos( total_phase ) * prop_sign
          ! emf%b_cyl_m%pf_im(1)%f2(3,i1,i2) = - lenv_2 * per_env_amp * sin( total_phase ) * prop_sign

        ! same as general ell except with sign switches for imaginary parts
        case( -1 ) 

          ! ! e_r 
          ! call this % per_envelope( z, r_2, per_env_amp, per_env_phase )
          ! total_phase = per_env_phase + fast_phase

          ! emf%e_cyl_m%pf_re(1)%f2(2,i1,i2) = + lenv * per_env_amp * cos( total_phase )
          ! emf%e_cyl_m%pf_im(1)%f2(2,i1,i2) = + lenv * per_env_amp * sin( total_phase )

          ! ! e_phi 
          ! call this % per_envelope( z, r, per_env_amp, per_env_phase )
          ! total_phase = per_env_phase + fast_phase
          
          ! emf%e_cyl_m%pf_re(1)%f2(3,i1,i2) = + lenv * per_env_amp * sin( total_phase )      
          ! emf%e_cyl_m%pf_im(1)%f2(3,i1,i2) = - lenv * per_env_amp * cos( total_phase )      

          ! ! b_r
          ! call this % per_envelope( z_2, r, per_env_amp, per_env_phase )
          ! total_phase = per_env_phase + fast_phase_2

          ! emf%b_cyl_m%pf_re(1)%f2(2,i1,i2) = - lenv_2 * per_env_amp * sin( total_phase ) * prop_sign           
          ! emf%b_cyl_m%pf_im(1)%f2(2,i1,i2) = + lenv_2 * per_env_amp * cos( total_phase ) * prop_sign           

          ! ! b_phi 
          ! call this % per_envelope( z_2, r_2, per_env_amp, per_env_phase )
          ! total_phase = per_env_phase + fast_phase_2
          
          ! emf%b_cyl_m%pf_re(1)%f2(3,i1,i2) = + lenv_2 * per_env_amp * cos( total_phase ) * prop_sign
          ! emf%b_cyl_m%pf_im(1)%f2(3,i1,i2) = + lenv_2 * per_env_amp * sin( total_phase ) * prop_sign
          
        end select 
      
      ! else: is a laguerre pulse with ell /= 0 
      else

        select case( this%pol_type ) 

        ! lin pol 
        case( 0 )

          ! e_r: z, r_2
          e_cyl_m%pf_re(ell-1)%f2(2,i1,i2) = + 0.5 * real( field_0_2 * exp( cmplx( 0, this%pol, p_k_fld ) ) )

          if( ell-1 > 0 ) then
            e_cyl_m%pf_im(ell-1)%f2(2,i1,i2) = - 0.5 * aimag( field_0_2 * exp( cmplx( 0, this%pol, p_k_fld ) ) )
          endif

          e_cyl_m%pf_re(ell+1)%f2(2,i1,i2) = + 0.5 * real( field_0_2 * exp( cmplx( 0, -this%pol, p_k_fld ) ) )

          e_cyl_m%pf_im(ell+1)%f2(2,i1,i2) = - 0.5 * aimag( field_0_2 * exp( cmplx( 0, -this%pol, p_k_fld ) ) )


          ! e_phi: z, r
          e_cyl_m%pf_re(ell-1)%f2(3,i1,i2) = + 0.5 * aimag( field_0_0 * exp( cmplx( 0, this%pol, p_k_fld ) ) )

          if( ell-1 > 0 ) then
            e_cyl_m%pf_im(ell-1)%f2(3,i1,i2) = + 0.5 * real( field_0_0 * exp( cmplx( 0, this%pol, p_k_fld ) ) )
          endif

          e_cyl_m%pf_re(ell+1)%f2(3,i1,i2) = - 0.5 * aimag( field_0_0 * exp( cmplx( 0, -this%pol, p_k_fld ) ) )

          e_cyl_m%pf_im(ell+1)%f2(3,i1,i2) = - 0.5 * real( field_0_0 * exp( cmplx( 0, -this%pol, p_k_fld ) ) )


          ! b_r: z_2, r
          b_cyl_m%pf_re(ell-1)%f2(2,i1,i2) = - 0.5 * aimag( field_2_0 * exp( cmplx( 0, this%pol, p_k_fld ) ) ) * prop_sign

          if( ell-1 > 0 ) then
            b_cyl_m%pf_im(ell-1)%f2(2,i1,i2) = - 0.5 * real( field_2_0 * exp( cmplx( 0, this%pol, p_k_fld ) ) ) * prop_sign
          endif

          b_cyl_m%pf_re(ell+1)%f2(2,i1,i2) = + 0.5 * aimag( field_2_0 * exp( cmplx( 0, -this%pol, p_k_fld ) ) ) * prop_sign

          b_cyl_m%pf_im(ell+1)%f2(2,i1,i2) = + 0.5 * real( field_2_0 * exp( cmplx( 0, -this%pol, p_k_fld ) ) ) * prop_sign


          ! b_phi: z_2, r_2
          b_cyl_m%pf_re(ell-1)%f2(3,i1,i2) = + 0.5 * real( field_2_2 * exp( cmplx( 0, this%pol, p_k_fld ) ) ) * prop_sign

          if( ell-1 > 0 ) then
            b_cyl_m%pf_im(ell-1)%f2(3,i1,i2) = - 0.5 * aimag( field_2_2 * exp( cmplx( 0, this%pol, p_k_fld ) ) ) * prop_sign
          endif

          b_cyl_m%pf_re(ell+1)%f2(3,i1,i2) = + 0.5 * real( field_2_2 * exp( cmplx( 0, -this%pol, p_k_fld ) ) ) * prop_sign

          b_cyl_m%pf_im(ell+1)%f2(3,i1,i2) = - 0.5 * aimag( field_2_2 * exp( cmplx( 0, -this%pol, p_k_fld ) ) ) * prop_sign


        ! handle circ pol 
        case( 1 )
          
          ! ! e_r 
          ! call this % per_envelope_laguerre( z, r_2, per_env_amp, per_env_phase )
          ! total_phase = per_env_phase + fast_phase

          ! emf%e_cyl_m%pf_re(ell+1)%f2(2,i1,i2) = lenv * per_env_amp * cos( total_phase )

          ! emf%e_cyl_m%pf_im(ell+1)%f2(2,i1,i2) = - lenv * per_env_amp * sin( total_phase )

          ! ! e_phi 
          ! call this % per_envelope_laguerre( z, r, per_env_amp, per_env_phase )
          ! total_phase = per_env_phase + fast_phase
          
          ! emf%e_cyl_m%pf_re(ell+1)%f2(3,i1,i2) = - lenv * per_env_amp * sin( total_phase )

          ! emf%e_cyl_m%pf_im(ell+1)%f2(3,i1,i2) = - lenv * per_env_amp * cos( total_phase )

          ! ! b_r
          ! call this % per_envelope_laguerre( z_2, r, per_env_amp, per_env_phase )
          ! total_phase = per_env_phase + fast_phase_2

          ! emf%b_cyl_m%pf_re(ell+1)%f2(2,i1,i2) = + lenv_2 * per_env_amp * sin( total_phase ) * prop_sign

          ! emf%b_cyl_m%pf_im(ell+1)%f2(2,i1,i2) = + lenv_2 * per_env_amp * cos( total_phase ) * prop_sign

          ! ! b_phi 
          ! call this % per_envelope_laguerre( z_2, r_2, per_env_amp, per_env_phase )
          ! total_phase = per_env_phase + fast_phase_2
          
          ! emf%b_cyl_m%pf_re(ell+1)%f2(3,i1,i2) = + lenv_2 * per_env_amp * cos( total_phase ) * prop_sign

          ! emf%b_cyl_m%pf_im(ell+1)%f2(3,i1,i2) = - lenv_2 * per_env_amp * sin( total_phase ) * prop_sign


        case( -1 ) 

          ! ! e_r 
          ! call this % per_envelope_laguerre( z, r_2, per_env_amp, per_env_phase )
          ! total_phase = per_env_phase + fast_phase

          ! emf%e_cyl_m%pf_re(ell-1)%f2(2,i1,i2) = + lenv * per_env_amp * cos( total_phase + this%pol )
          
          ! if( ell-1 > 0 ) then 
          !   emf%e_cyl_m%pf_im(ell-1)%f2(2,i1,i2) = - lenv * per_env_amp * sin( total_phase + this%pol )
          ! endif          

          ! ! e_phi 
          ! call this % per_envelope_laguerre( z, r, per_env_amp, per_env_phase )
          ! total_phase = per_env_phase + fast_phase
          
          ! emf%e_cyl_m%pf_re(ell-1)%f2(3,i1,i2) = + lenv * per_env_amp * sin( total_phase + this%pol )
          
          ! if( ell-1 > 0 ) then 
          !   emf%e_cyl_m%pf_im(ell-1)%f2(3,i1,i2) = + lenv * per_env_amp * cos( total_phase + this%pol )
          ! endif          

          ! ! b_r
          ! call this % per_envelope_laguerre( z_2, r, per_env_amp, per_env_phase )
          ! total_phase = per_env_phase + fast_phase_2

          ! emf%b_cyl_m%pf_re(ell-1)%f2(2,i1,i2) = - lenv_2 * per_env_amp * sin( total_phase + this%pol ) * prop_sign          

          ! if( ell-1 > 0 ) then 
          !   emf%b_cyl_m%pf_im(ell-1)%f2(2,i1,i2) = - lenv_2 * per_env_amp * cos( total_phase + this%pol ) * prop_sign
          ! endif        

          ! ! b_phi 
          ! call this % per_envelope_laguerre( z_2, r_2, per_env_amp, per_env_phase )
          ! total_phase = per_env_phase + fast_phase_2
          
          ! emf%b_cyl_m%pf_re(ell-1)%f2(3,i1,i2) = + lenv_2 * per_env_amp * cos( total_phase + this%pol ) * prop_sign
          
          ! if( ell-1 > 0 ) then 
          !   emf%b_cyl_m%pf_im(ell-1)%f2(3,i1,i2) = - lenv_2 * per_env_amp * sin( total_phase + this%pol ) * prop_sign
          ! endif    

        end select 

      endif
    
    enddo ! i2 
  
  enddo ! i1 



end subroutine set_fields_t0_zpulse_flyfoc_cyl_modes
!-----------------------------------------------------------------------------------------



end module m_zpulse_flyfoc_cyl_modes
