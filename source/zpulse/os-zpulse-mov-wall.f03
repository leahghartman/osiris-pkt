#include "os-config.h"
#include "os-preprocess.fpp"

module m_zpulse_mov_wall

use m_system
use m_parameters
use m_math

use m_zpulse_std
use m_zpulse_wall

use m_node_conf
use m_space
use m_vdf_define
use m_grid_define

use m_fparser

use m_emf_define, only: t_emf_bound, t_emf
use m_emf_bound, only: is_open

implicit none

private


!-----------------------------------------------------------------------------------------
type, extends(t_zpulse) :: t_zpulse_mov_wall

  ! whether to initialize fields fully at t=0 
  logical :: if_launch_t0

  integer :: tenv_type

  real(p_double)  :: inj_time

  real(p_double)  :: tenv_rise, tenv_flat, tenv_fall
  real(p_double)  :: tenv_duration, tenv_range, tenv_launch_duration
  type(t_fparser) :: tenv_math_func
  real(p_double)  :: xi0
  integer         :: ncells 

  logical         :: if_tenv_tilt
  real(p_double), dimension(2) :: tenv_tilt

  real(p_double)  :: wall_pos  ! initial wall position
  real(p_double)  :: wall_vel  ! velocity (must be >= 0 and <= 1)
  logical         :: if_launch_from_wall 

  logical         :: if_clean_fields 
  logical         :: if_clean_x2
  real(p_double)  :: min_x2_clean, max_x2_clean 

  real(p_double)  :: clean_fields_delta
  real(p_double)  :: clean_fields_time 

  logical         :: if_apply_divcorr
  logical         :: if_add_fields

contains

  procedure :: init       => init_zpulse_mov_wall
  procedure :: read_input => read_input_zpulse_mov_wall
  procedure :: launch => launch_zpulse_mov_wall
  procedure :: t_duration => t_duration_zpulse_mov_wall
  procedure :: t_envelope => t_env_zpulse_mov_wall
  ! procedure :: clean_fields => clean_fields_zpulse_mov_wall 
  procedure :: lon_center => lon_center_zpulse_mov_wall

  procedure :: launch_wall_1d => launch_wall_1d_zpulse_mov_wall
  procedure :: launch_wall_2d => launch_wall_2d_zpulse_mov_wall
  procedure :: launch_wall_3d => launch_wall_3d_zpulse_mov_wall

  procedure :: apply_divcorr_2d => apply_divcorr_2d_zpulse_mov_wall
  procedure :: apply_divcorr_3d => apply_divcorr_3d_zpulse_mov_wall
  procedure :: get_env_1d => get_env_1d_zpulse_mov_wall
  procedure :: get_env_2d => get_env_2d_zpulse_mov_wall
  procedure :: get_env_3d => get_env_3d_zpulse_mov_wall
  procedure :: fast_phase => fast_phase_zpulse

  ! procedure :: clean_fields_2d => clean_fields_2d_zpulse_mov_wall
  ! procedure :: clean_fields_3d => clean_fields_3d_zpulse_mov_wall


end type t_zpulse_mov_wall
!-----------------------------------------------------------------------------------------



public :: t_zpulse_mov_wall

contains

!-----------------------------------------------------------------------------------------
subroutine init_zpulse_mov_wall( this, restart, t, dt, b, g_space, nx_p_min, no_co, interpolation )

  implicit none

  class( t_zpulse_mov_wall ), intent(inout) :: this
  logical, intent(in) :: restart
  real( p_double ), intent(in) :: t, dt
  type( t_vdf ) ,  intent(in)  :: b
  type( t_space ), intent(in) :: g_space
  integer, intent(in), dimension(:) :: nx_p_min
  class( t_node_conf ), intent(in) :: no_co
  integer, intent(in) :: interpolation

  logical :: if_launch

  ! Need to preserve if_launch value for this zpulse
  if_launch = this%if_launch

  call this % t_zpulse % init( restart, t, dt, b, g_space, nx_p_min, no_co, interpolation )

  this%if_launch = if_launch

  if ( this%if_launch .and. ( t > this%launch_time + this % t_duration() ) ) then
    this%if_launch = .false.
  endif

end subroutine init_zpulse_mov_wall
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine read_input_zpulse_mov_wall( this, input_file, g_space, bnd_con, periodic, grid, sim_options )

  use m_input_file

  implicit none

  class( t_zpulse_mov_wall ), intent(inout) :: this
  class( t_input_file ), intent(inout) :: input_file
  type( t_space ), intent(in) :: g_space
  class (t_emf_bound), intent(in) :: bnd_con
  logical, dimension(:), intent(in) :: periodic
  class( t_grid ), intent(in) :: grid
  type( t_options ), intent(in) :: sim_options

  real(p_double) :: a0
  real(p_double) :: omega0
  real(p_double) :: k0
  real(p_double) :: phase
  integer        :: pol_type
  real(p_double) :: pol
  character(len=16) :: propagation
  integer :: direction 
  real(p_double) :: prop_sign

  integer            :: chirp_order
  real(p_double), dimension(p_max_chirp_order) :: chirp_coefs

  logical :: if_launch_t0

  ! Temporal profile
  character(len=16) :: tenv_type

  real(p_double) :: tenv_fwhm

  real(p_double) :: tenv_rise, tenv_flat, tenv_fall
  real(p_double) :: tenv_duration, tenv_range, tenv_launch_duration
  real(p_double) :: xi0 

  character(len = p_max_expr_len) :: tenv_math_func
  real(p_double), dimension(2) :: tenv_tilt

  ! Beam profile
  character(len=16) :: per_type
  real(p_double), dimension(2) :: per_center

  integer            :: per_chirp_order
  real(p_double), dimension(p_max_chirp_order) :: per_chirp_coefs

  real(p_double), dimension(2) :: per_w0, per_fwhm, per_focus
  integer, dimension(2) :: per_tem_mode

  real(p_double), dimension(2,2) :: per_w0_asym, per_fwhm_asym, per_focus_asym
  real(p_double), dimension(2) :: per_asym_trans

  integer         :: per_n, per_0clip
  real(p_double) :: per_kt

  real(p_double) :: launch_time
  logical        :: if_launch

  real(p_double) :: wall_pos, wall_vel, clean_fields_delta
  integer        :: ncells 
  logical        :: if_clean_fields, if_launch_from_wall, if_apply_divcorr, if_add_fields

  namelist /nl_zpulse_mov_wall/ &
    a0, omega0, k0, phase, pol_type, pol, propagation, direction, &
    chirp_order, chirp_coefs, xi0, tenv_fwhm, &
    if_launch_t0, tenv_type, tenv_rise, tenv_flat, tenv_fall, &
    tenv_duration, tenv_range, tenv_launch_duration, tenv_math_func, tenv_tilt, &
    per_type, per_center, per_w0, per_fwhm, per_focus, &
    per_chirp_order, per_chirp_coefs, &
    per_w0_asym, per_fwhm_asym, per_focus_asym, per_asym_trans, &
    per_n, per_0clip, per_kt, per_tem_mode, &
    launch_time, if_launch, wall_pos, wall_vel, if_launch_from_wall, ncells, &
    if_clean_fields, clean_fields_delta, if_apply_divcorr, if_add_fields

  integer :: i, ierr

  if (disp_out(input_file)) then
    SCR_ROOT(" - Reading zpulse_mov_wall configuration...")
  endif

  if_launch = .true.

  a0             = 1.0_p_double
  omega0         = 10.0_p_double
  k0             = 0.0_p_double
  phase          = 0.0_p_double
  pol_type       = 0
  pol            = 90.0_p_double

  chirp_order    = 0
  chirp_coefs    = 0.0_p_double

  propagation    = "forward"
  direction      = 1 

  if_launch_t0   = .false.
  tenv_type       = "polynomial"

  tenv_fwhm       = -1.0_p_double

  tenv_rise       = 0.0_p_double
  tenv_flat       = 0.0_p_double
  tenv_fall       = 0.0_p_double

  xi0             = 0.0_p_double


  tenv_duration   = 0.0_p_double
  tenv_range      = -huge(1.0_p_double)
  tenv_launch_duration = -huge(1.0_p_double)
  tenv_math_func  = "NO_FUNCTION_SUPPLIED!"

  tenv_tilt       = 0.0_p_double

  per_type       = "plane"
  per_center     = -huge(1.0_p_double)

  per_chirp_order    = 0
  per_chirp_coefs    = 0.0_p_double

  per_w0         = 0.0_p_double
  per_fwhm       = 0.0_p_double
  per_focus      = -huge(1.0_p_double)
  per_tem_mode = 0 ! tem mode, integer, default = 0,0

  per_w0_asym     = 0.0_p_double
  per_fwhm_asym   = 0.0_p_double
  per_focus_asym  = -huge(1.0_p_double)
  per_asym_trans  = 0.0_p_double

  per_n          = 0
  per_0clip      = 1
  per_kt         = 0.0_p_double

  wall_vel = 0.0_p_double
  wall_pos = -huge(1.0d0)
  ncells = 4
  if_clean_fields = .false. 
  clean_fields_delta = 0.0_p_double
  if_launch_from_wall = .false.

  if_apply_divcorr = .true. 
  if_add_fields = .false. 


  launch_time           = 0.0_p_double


  ! read data from file
  call get_namelist( input_file, "nl_zpulse_mov_wall", ierr )

  if ( ierr /= 0 ) then
    if ( mpi_node() == 0 ) then
      if (ierr < 0) then
        print *, "Error reading zpulse_mov_wall parameters"
      else
        print *, "Error: zpulse_mov_wall parameters missing"
      endif
      print *, "aborting..."
    endif
    stop
  endif

  read (input_file%nml_text, nml = nl_zpulse_mov_wall, iostat = ierr)
  if (ierr /= 0) then
    if ( mpi_node() == 0 ) then
      write(0,*)  ""
      write(0,*)  "   Error reading zpulse_mov_wall parameters"
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
  this%phase0     = real( phase * pi_180, p_double )

  ! Polarization parameters
  if ((pol_type <-1) .or. (pol_type >1)) then
    if ( mpi_node() == 0 ) then
      write(0,*)  ""
      write(0,*)  "   Error in zpulse_mov_wall parameters"
      write(0,*)  "   pol must be in the range [-1,0,+1]"
      write(0,*)  "   aborting..."
    endif
    stop
  endif
  this%pol_type  = pol_type
  this%pol = real( pol * pi_180, p_double )

  this % xi0 = xi0 


  ! Frequency chirp parameters
  if ((chirp_order < 0) .or. (chirp_order > p_max_chirp_order)) then
    if ( mpi_node() == 0 ) then
      write(0,*)  ""
      write(0,*)  "   Error in zpulse_mov_wall parameters"
      write(0,*)  "   chirp_order must be in the range [0,",p_max_chirp_order,"]"
      write(0,*)  "   aborting..."
    endif
    stop
  endif
  this%chirp_order    = chirp_order
  this%chirp_coefs    = chirp_coefs

  if ((per_chirp_order < 0) .or. (per_chirp_order > p_max_chirp_order)) then
    if ( mpi_node() == 0 ) then
      write(0,*)  ""
      write(0,*)  "   Error in zpulse_mov_wall parameters"
      write(0,*)  "   per_chirp_order must be in the range [0,",p_max_chirp_order,"]"
      write(0,*)  "   aborting..."
    endif
    stop
  endif
  this%per_chirp_order    = per_chirp_order
  this%per_chirp_coefs    = per_chirp_coefs

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
      write(0,*)  '   Error in zpulse_mov_wall parameters'
      write(0,*)  '   propagation must be either "forward" or "backward"'
      write(0,*)  '   aborting...'
    endif
    stop
  end select

  if ((direction < 1) .or. (direction > p_x_dim)) then
     if ( mpi_node() == 0 ) then
       write(0,*)  ""
       write(0,*)  "   Error in zpulse_mov_wall parameters"
       write(0,*)  "   direction must be in the range [1 .. x_dim]"
       write(0,*)  "   aborting..."
     endif
     stop
   endif

   this%direction = direction

  ! Temporal envelope tilt
  this%if_tenv_tilt = .false.
  if ( p_x_dim > 1 ) then
    do i = 1, p_x_dim-1
      if ( tenv_tilt(i) /= 0.0_p_double ) then
        if (( tenv_tilt(i) < -45.0_p_double ) .or. (tenv_tilt(i) > 45.0_p_double )) then
          if ( mpi_node() == 0 ) then
            write(0,*)  ""
            write(0,*)  "   Error in zpulse_mov_wall parameters"
            write(0,*)  "   tenv_tilt must be in the range [ -45 , 45 ] deg"
            write(0,*)  "   aborting..."
          endif
          stop
        endif
        this%tenv_tilt(i) = tan( real( tenv_tilt(i) * pi_180, p_double ))
        if ( this%propagation == p_forward ) this%tenv_tilt(i) = -this%tenv_tilt(i)
        this%if_tenv_tilt = .true.
      else
        this%tenv_tilt(i) = 0.0_p_double
      endif
    enddo
  endif


  ! Temporal envelope parameters
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
    this%tenv_launch_duration  = tenv_launch_duration

    if ( tenv_launch_duration == -huge(1.0_p_double) ) then
      if ( mpi_node() == 0 ) then
        write(0,*)  ''
        write(0,*)  '   Error in zpulse_mov_wall parameters'
        write(0,*)  '   When using "gaussian" temporal envelope, '
        write(0,*)  '   tenv_launch_duration must also be defined'
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
        write(0,*)  '   Error in zpulse_mov_wall parameters'
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
        write(0,*)  '   Error in zpulse_mov_wall parameters'
        write(0,*)  '   Supplied tenv_math_func failed to compile'
        write(0,*)  '   aborting...'
      endif
      stop
    endif

  case default
    if ( mpi_node() == 0 ) then
      write(0,*)  ''
      write(0,*)  '   Error in zpulse_mov_wall parameters'
      write(0,*)  '   lon_type must be either "polynomial", "gaussian",'
      write(0,*)  '   "sin2", "math" or "hermite"'
      write(0,*)  '   aborting...'
    endif
    stop
  end select



  ! perpendiular profile
  select case ( trim(per_type) )
  case ("plane")
    this%per_type       = p_plane

  case ("gaussian","hermite")

    this%per_type       = p_hermite_gaussian
    if ( p_x_dim > 1 ) then
      this%per_tem_mode = 0
      do i = 1, p_x_dim - 1
        this%per_tem_mode(i) = per_tem_mode(i)
        if ( per_tem_mode(i) > 7 ) then
          if (mpi_node() == 0) then
            write(0,*)  ''
            write(0,*)  '   Error in zpulse_mov_wall parameters'
            write(0,*)  '   Only TEM modes up to 7 have been implemented. Please contact the'
            write(0,*)  '   development team if you need higher order modes.'
            write(0,*)  '   aborting...'
          endif
          stop
        endif
      enddo
    endif

    ! get main spot size and focal plane
    this%per_w0 = 0.0_p_double
    if (per_fwhm(1) > 0.0_p_double) then
      this%per_w0 = per_fwhm(1)/sqrt(4.0_p_double * log(2.0_p_double))
    else
      this%per_w0 = per_w0(1)
    endif

    if ( this%per_w0(1) <= 0.0_p_double ) then
      if (mpi_node() == 0) then
        write(0,*)  ''
        write(0,*)  '   Error in zpulse_mov_wall parameters'
        write(0,*)  '   for gaussian pulses per_w0 or per_fwhm must be > 0'
        write(0,*)  '   aborting...'
      endif
      stop
    endif

    if ( per_focus(1) == -huge(1.0_p_double) ) then
      if (mpi_node() == 0) then
        write(0,*)  ''
        write(0,*)  '   Error in zpulse_mov_wall parameters'
        write(0,*)  '   for gaussian pulses per_focus must be defined'
        write(0,*)  '   aborting...'
      endif
      stop
    endif
    this%per_focus      = per_focus(1)

    ! check for astigmatic beam
    if (p_x_dim == 3) then
      if (per_fwhm(2) > 0.0_p_double) then
        this%per_w0(2) = per_fwhm(2)/sqrt(4.0_p_double * log(2.0_p_double))
      else if ( per_w0(2) > 0.0_p_double ) then
        this%per_w0(2) = per_w0(2)
      endif

      if ( per_focus(2) /= -huge(1.0_p_double) ) then
        this%per_focus(2)      = per_focus(2)
      else
        this%per_focus(2)      = this%per_focus(1)
      endif

      if ((this%per_w0(2) /= this%per_w0(1)) .or. (this%per_focus(2) /= this%per_focus(1))) then
        this%per_type = p_hermite_gaussian_astigmatic
      endif

    endif

  case ("laguerre")

    this%per_type       = p_laguerre_gaussian
    if ( p_x_dim > 1 ) then
      this%per_tem_mode = 0
      do i = 1, 2
        this%per_tem_mode(i) = per_tem_mode(i)
        if ( abs(per_tem_mode(i)) > 7 ) then
          if (mpi_node() == 0) then
            print *, ''
            print *, '   Error in zpulse_mov_wall parameters'
            print *, '   Only Laguerre modes up to 7 have been implemented. Please contact the'
            print *, '   development team if you need higher order modes.'
            print *, '   aborting...'
          endif
          stop
        endif
      enddo
    endif


    ! get main spot size and focal plane
    this%per_w0 = 0.0_p_double
    if (per_fwhm(1) > 0.0_p_double) then
      this%per_w0 = per_fwhm(1)/sqrt(4.0_p_double * log(2.0_p_double))
    else
      this%per_w0 = per_w0(1)
    endif

    if ( this%per_w0(1) <= 0.0_p_double ) then
      if (mpi_node() == 0) then
        print *, ''
        print *, '   Error in zpulse_mov_wall parameters'
        print *, '   for Laguerre-Gaussian pulses per_w0 or per_fwhm must be > 0'
        print *, '   aborting...'
      endif
      stop
    endif

    if ( per_focus(1) == -huge(1.0_p_double) ) then
      if (mpi_node() == 0) then
        print *, ''
        print *, '   Error in zpulse_mov_wall parameters'
        print *, '   for Laguerre-Gaussian pulses per_focus must be defined'
        print *, '   aborting...'
      endif
      stop
    endif
    this%per_focus      = per_focus(1)

    ! default parameters
    if (per_fwhm(1) > 0.0_p_double) then
      this%per_w0(1) = per_fwhm(1)/sqrt(4.0_p_double * log(2.0_p_double))
    else
      this%per_w0(1) = per_w0(1)
    endif
    this%per_focus(1)      = per_focus(1)

  case default
    if (mpi_node() == 0) then
      write(0,*)  ''
      write(0,*)  '   Error in zpulse_mov_wall parameters'
      write(0,*)  '   per_type must be one of the following:'
      write(0,*)  '   "plane", "gaussian"/"hermite"'
      write(0,*)  '   aborting...'
    endif
    stop
  end select

  ! Set default per_center values
  if ( p_x_dim > 1 ) then
    if ( per_center(1) == -huge(1.0_p_double) ) then
        per_center(1) = 0.5_p_double * ( xmin( g_space, 2 ) + xmax( g_space, 2 ) )
    endif

    if ( p_x_dim == 3 ) then
      if ( per_center(2) == -huge(1.0_p_double) ) then
          per_center(2) = 0.5_p_double * ( xmin( g_space, 3 ) + xmax( g_space, 3 ) )
      endif
    endif
  endif

  this%per_center     = per_center

  if ( wall_pos == -huge(1.0d0) ) then
    if ( mpi_node() == 0 ) then
      write(0,*)  ""
      write(0,*)  "   Error in zpulse_mov_wall parameters"
      write(0,*)  "   When using a moving wall, the initial position (wall_pos) must be set"
      write(0,*)  "   aborting..."
    endif
    stop
  endif

  this%if_launch_t0 = if_launch_t0

  this%wall_pos = wall_pos
  this%wall_vel = wall_vel
  this%ncells   = ncells
  this%if_clean_fields = if_clean_fields 
  this%clean_fields_delta = clean_fields_delta
  this%if_launch_from_wall = if_launch_from_wall
  this%if_apply_divcorr = if_apply_divcorr
  this%if_add_fields = if_add_fields

  this%launch_time    = launch_time
  this%if_launch = if_launch

  ! Boost pulse
  this%if_boost = .false.
  
  if (sim_options%gamma > 1.0_p_double) then

    this % gamma          = sim_options%gamma
    this % gamma_times_beta = sqrt(this%gamma**2-1.0_p_double)

    if( this % propagation .eq. p_forward ) then
      this % boost_per_fac = this % gamma - this % gamma_times_beta
    else
      this % boost_per_fac = this % gamma + this % gamma_times_beta
    endif 
    SCR_ROOT('   boosting zpulse with gamma = ', this%gamma )

    this%if_boost = .true.
  else
    this % boost_per_fac = 1.0_p_double
  endif

  this % if_clean_fields = if_clean_fields 

  this % clean_fields_time = this % t_duration() + this % clean_fields_delta 

end subroutine
!-----------------------------------------------------------------------------------------



!-----------------------------------------------------------------------------------------
function t_env_zpulse_mov_wall( this, t_sim, z_sim ) result(t_env)

  implicit none

  class( t_zpulse_mov_wall ), intent(inout) :: this ! intent must be inout because of
                                                ! the eval( f_parser ) function
  real(p_double), intent(in) :: t_sim
  real(p_double), intent(in) :: z_sim 

  ! lab coordinates 
  real(p_double) :: z, t

  real(p_double) :: t_env

  real(p_k_fparse), dimension(1) ::t_dbl

  ! convert to lab frame coords 
  if( this % if_boost ) then 
    call this % apply_lorentz_transformation( z_sim, t_sim, z, t )
  else 
    z = z_sim 
    t = t_sim 
  endif 

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

end function t_env_zpulse_mov_wall
!-----------------------------------------------------------------------------------------



!-----------------------------------------------------------------------------------------
! Total duration of the laser pulse
!-----------------------------------------------------------------------------------------
function t_duration_zpulse_mov_wall( this ) result( t_duration )

  implicit none

  class( t_zpulse_mov_wall ), intent(in) :: this
  real(p_double) :: t_duration

  select case ( this%tenv_type )
  case (p_polynomial, p_sin2)
    t_duration = this%tenv_rise + this%tenv_flat + this%tenv_fall

  case default
    t_duration = this%tenv_launch_duration

  end select

  if( this % if_boost ) then
    t_duration = this % gamma * t_duration!  - this % gamma_times_beta * this % wall_pos
  endif

end function t_duration_zpulse_mov_wall
!----------------------------------------------------------------------------------------



!----------------------------------------------------------------------------------------
function lon_center_zpulse_mov_wall( this ) result( lon_center )

  implicit none

  class( t_zpulse_mov_wall ), intent(in)  :: this
  real(p_double) :: lon_center

  select case ( this%tenv_type )

  case (p_polynomial, p_sin2)
    ! the center of the pulse is define as being
    ! half way in the flat region
    lon_center = this % tenv_rise + this % tenv_flat / 2

  case (p_gaussian, p_func)
    lon_center = 0.5*this % tenv_range

  case default
    lon_center = 0

  end select

  lon_center = lon_center

end function lon_center_zpulse_mov_wall
!----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!       returns the value of the perpendicular
!       envelope at the requested position
!-----------------------------------------------------------------------------------------
subroutine get_env_2d_zpulse_mov_wall( this, tin_sim, x1_sim, x2, env )

  implicit none

  class( t_zpulse_mov_wall ), intent(inout) :: this

  integer, parameter :: rank = 2

  ! x1 is the longitudinal coordinate, x2 the perpendicular one
  real(p_double), intent(in) :: x1_sim
  real(p_double), intent(in) :: x2
  real(p_double), intent(in) :: tin_sim

  ! lab frame coordinates 
  real(p_double) :: x1, t

  ! simulation time 
  real(p_double) :: t_sim 

  ! local variables
  real(p_double) :: z_center, r_center
  real(p_double) :: k, z, rWl2, Wl2, z0, rho, rho2
  real(p_double) :: gouy_shift, rWu, curv

  ! t_envelope
  real(p_double) :: tenv

  ! asymetric laser pulses variables
  real(p_double) :: w0_asym, focus_asym, shp

  complex(p_double), intent(out) :: env

  integer :: i

  ! ! Boost variables
  ! real(p_double) :: uvel, beta, rg1pb

  !if( present( tin_sim ) ) then 
  t_sim = tin_sim
  !else 
  !  t_sim = 0.0_p_double
  !endif 
  
  tenv = this % a0 * this % t_envelope(t_sim, x1_sim)

  if( this % if_boost ) then 
    call this % apply_lorentz_transformation( x1_sim, t_sim, x1, t )
  else 
    x1 = x1_sim 
    t = t_sim 
  endif 

  ! executable statements

  ! get pulse center for chirp and phase calculations
  z_center = this%lon_center()
  r_center = this%per_center(1)

  ! get wavenumber
  k = this%omega0
  ! add chirp
  if ( this%propagation == p_forward ) then
    do i = 1, this%chirp_order
      k = k + this%chirp_coefs(i) * ((x1 - z_center)**(i))
    enddo
    do i = 1, this%per_chirp_order
      k = k + this%per_chirp_coefs(i) * ((x2 - r_center)**(i))
    enddo
  else
    do i = 1, this%chirp_order
      k = k + this%chirp_coefs(i) * ((z_center - x1)**(i))
    enddo
    do i = 1, this%per_chirp_order
      k = k + this%per_chirp_coefs(i) * ((r_center - x2)**(i))
    enddo
  endif

  ! get envelope
  select case ( this%per_type )

  case (p_plane)
    env = tenv

  case (p_hermite_gaussian)

    ! get rayleigh range
    z0 = k * this%per_w0(1)**2 * 0.5_p_double

    ! calculate radius
    rho = x2 - r_center
    rho2 = rho**2

    ! calculate gaussian beam parameters
    z  = x1 - this%per_focus(1)

    ! In some situations (i.e. chirps) z0 may be 0 so z = 0 must be treated
    ! as a special case
    if ( z /= 0 ) then
      rWl2 = z0**2 / (z0**2 + z**2)
      curv = 0.5_p_double*rho2*z/(z**2 + z0**2)
      gouy_shift = ( this%per_tem_mode(1) + 1 ) * atan2( z, z0 )
    else
      rWl2 = 1
      curv = 0
      gouy_shift = 0
    endif
    rWu = sqrt(2 * rWl2) / this%per_w0(1)

    ! note that this is different in 2D and 3D
    env = tenv * sqrt( sqrt(rWl2) ) * &
    hermite( this%per_tem_mode(1), rho * rWu ) * &
    exp( - rho2 * rWl2/this%per_w0(1)**2 ) * &
    exp( cmplx(0,1.0_p_double)*( k * curv - gouy_shift) )

  case (p_gaussian_asym)

    ! THIS DOES NOT WORK WITH CHIRPED PULSES
    ! (there is a possibility of a division by 0 in that case)

    ! calculate radius
    rho = x2 - r_center
    rho2 = rho**2

    ! get w0 and focsus as a smooth transition between the 2 regimes
    shp = smooth_heaviside( rho/this%per_asym_trans(1) )
    w0_asym = this%per_w0_asym(1, 1) + shp*( this%per_w0_asym(2, 1) - this%per_w0_asym(1, 1))
    focus_asym = this%per_focus_asym(1,1) + &
    shp*(this%per_focus_asym(2,1) - this%per_focus_asym(1,1))

    ! get rayleigh range
    z0 = k * w0_asym**2 * 0.5_p_double

    ! calculate gaussian beam parameters
    z  = x1 - focus_asym

    Wl2 = 1.0_p_double + (z/z0)**2

    gouy_shift = atan2( z, z0 )

    ! note that this is different in 2D and 3D
    env = tenv * sqrt( 1.0_p_double /  sqrt(Wl2) ) * &
    exp( - rho2/(w0_asym**2 * Wl2) ) * &
    exp( cmplx(0,1.0_p_double) * ( 0.5_p_double*k*rho2*z/(z**2 + z0**2) - gouy_shift) )

  case default

    ! This should never happen
    env = 0

  end select

end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!       returns the value of the perpendicular
!       envelope at the requested position
!-----------------------------------------------------------------------------------------
subroutine get_env_3d_zpulse_mov_wall( this, tin_sim, x1_sim, x2, x3, env )

  implicit none

  integer, parameter :: rank = 3

  real(p_double) :: per_envelope_3d

  class( t_zpulse_mov_wall ), intent(inout) :: this

  ! x1, x2, x3 are organized as lon, per1, per2
  real(p_double), intent(in) :: x1_sim
  real(p_double), intent(in) :: x2, x3
  real(p_double), intent(in) :: tin_sim

  ! lab-frame x1 and t
  real(p_double) :: x1, t

  ! simulation time 
  real(p_double) :: t_sim 

  ! local variables
  real(p_double) :: z_center, rx_center, ry_center
  real(p_double) :: kbessel, k, z, rWl2, z0, rho2, phi, ox, oy
  real(p_double) :: gouy_shift, rWu, curv, theta
  real(p_single) :: ktr
  real(p_double) :: tenv

  ! variables for astigmatic gaussian beams
  real(p_double) :: r_wx, r_wy, wx, wy, zx, zy, z0x, z0y, r_Rx, r_Ry, ox2, oy2
  real(p_double) :: gouy_angx, gouy_angy

  ! variables for asymmetric gaussian beams
  real(p_double), dimension(2) :: w0_asym, focus_asym, shp

  complex(p_double), intent(out) :: env

  integer :: i

  !if( present( tin_sim ) ) then 
  t_sim = tin_sim
  !else 
  !  t_sim = 0.0_p_double
  !endif 

  tenv = this % a0 * this % t_envelope(t_sim, x1_sim)

  if( this % if_boost ) then 
    call this % apply_lorentz_transformation( x1_sim, t_sim, x1, t )
  else 
    x1 = x1_sim 
    t = t_sim 
  endif 

  ! get pulse center for chirp and phase calculations
  z_center  = this % lon_center()
  rx_center = this%per_center(1)
  ry_center = this%per_center(2)

  !calculate axial distances
  ox = x2 - rx_center
  oy = x3 - ry_center

  ! get wavenumber
  k = this%omega0
  ! add chirp
  if ( this%propagation == p_forward ) then
    do i = 1, this%chirp_order
      k = k + this%chirp_coefs(i) * ((x1 - z_center)**(i))
    enddo
    do i = 1, this%per_chirp_order
      k = k + this%per_chirp_coefs(i) * (ox**(i))
    enddo
  else
    do i = 1, this%chirp_order
      k = k + this%chirp_coefs(i) * ((z_center - x1)**(i))
    enddo
    do i = 1, this%per_chirp_order
      k = k + this%per_chirp_coefs(i) * ((-ox)**(i))
    enddo
  endif

  
  select case ( this%per_type )

  case (p_plane)
    env = tenv

  case (p_hermite_gaussian)

    !calculate radius
    rho2 = ox**2 + oy**2

    z  = x1 - this%per_focus(1)

    z0 = k * this%per_w0(1)**2 /2

    ! In some situations (i.e. chirps) z0 may be 0 so z = 0 must be treated
    ! as a special case
    if ( z /= 0 ) then
      rWl2 = z0**2 / (z0**2 + z**2)
      curv = 0.5_p_double*rho2*z/(z**2 + z0**2)
      gouy_shift = (this%per_tem_mode(1) + this%per_tem_mode(2) + 1) * atan2( z, z0 )
    else
      rWl2 = 1
      curv = 0
      gouy_shift = 0
    endif

    rWu = sqrt(2 * rWl2) / this%per_w0(1)

    ! for mode(1) = 0 and mode(2) = 0 this reproduces the gaussian beam
    env = tenv * sqrt(rWl2) * &
    hermite( this%per_tem_mode(1), ox * rWu ) * &
    hermite( this%per_tem_mode(2), oy * rWu ) * &
    exp( - rho2 * rWl2 /this%per_w0(1)**2 ) * &
    exp( cmplx(0,1.0_p_double) * (k * curv - gouy_shift) )

  case (p_laguerre_gaussian)

    !calculate radius
    rho2 = ox**2 + oy**2

    !calculate phi
    theta = atan2(oy,ox)

    z  = x1 - this%per_focus(1)

    z0 = k * this%per_w0(1)**2 /2

    ! In some situations (i.e. chirps) z0 may be 0 so z = 0 must be treated
    ! as a special case
    if ( z /= 0 ) then
      rWl2 = z0**2 / (z0**2 + z**2)
      curv = 0.5_p_double*rho2*z/(z**2 + z0**2)
      gouy_shift = ( 2*this%per_tem_mode(1) + abs(this%per_tem_mode(2)) + 1) * atan2( z, z0 )
    else
      rWl2 = 1
      curv = 0
      gouy_shift = 0
    endif

    rWu = sqrt(2 * rWl2) / this%per_w0(1)

    ! for mode(1) = 0 and mode(2) = 0 this reproduces the gaussian beam
    ! if (.not. this%if_boost) then
    env = tenv * sqrt(rWl2) * &
    ( sqrt(rho2)*sqrt(rWl2)/this%per_w0(1) )**abs(this%per_tem_mode(2)) * &
    laguerre( this%per_tem_mode(1), abs(this%per_tem_mode(2)), &
      2*rho2*rWl2/this%per_w0(1)**2  ) * &
    exp( - rho2 * rWl2/this%per_w0(1)**2  ) * &
    exp( cmplx(0,1.0_p_double) * (k * curv - gouy_shift + this%per_tem_mode(2) * theta) )


  case (p_hermite_gaussian_astigmatic )

    ! get rayleigh ranges
    z0x = k * this%per_w0(1)**2 *0.5_p_double ! Rayleigh range for x
    z0y = k * this%per_w0(2)**2 *0.5_p_double ! Rayleigh range for y

    ! calculate transverse coordinates squared
    ox = x2 - rx_center
    oy = x3 - ry_center
    ox2 = ox**2
    oy2 = oy**2

    ! calculate distance to focal planes
    zx  = x1 - this%per_focus(1) ! distance to x focal plane
    zy  = x1 - this%per_focus(2) ! distance to y focal plane

    ! In some situations (i.e. chirps) z0x or z0y may be 0 so zx = 0 and zy = 0 must be treated
    ! as a special case

    if ( zx /= 0 ) then
      r_wx = z0x / sqrt( zx**2 + z0x**2 ) / this%per_w0(1) ! 1/local spot, x direction
      r_Rx = zx/(zx**2 + z0x**2) ! 1 / ( x curvature radius )
      gouy_angx = atan2( zx, z0x )
    else
      r_wx = 1/this%per_w0(1)
      r_Rx = 0
      gouy_angx = 0
    endif

    if ( zy /= 0 ) then
      r_wy = z0y / sqrt( zy**2 + z0y**2 ) / this%per_w0(2) ! 1/local spot, y direction
      r_Ry = zy/(zy**2 + z0y**2) ! 1 / ( y curvature radius )
      gouy_angy = atan2( zy, z0y )
    else
      r_wy = 1/this%per_w0(2)
      r_Ry = 0
      gouy_angy = 0
    endif

    ! gouy phase shift
    ! Half of the Guoy phase comes from each transverse coord. See
    ! Siegman, Lasers, 16.4 Higher Order Gaussian Modes, Astigmatic
    ! mode functions
    gouy_shift = 0.5_p_double*( gouy_angx + gouy_angy ) * &
    (this%per_tem_mode(1) + this%per_tem_mode(2) + 1)

    ! asymmetric gaussian envelope
    env = tenv * sqrt(this%per_w0(1) * this%per_w0(2) * r_wx * r_wy) * &
    hermite( this%per_tem_mode(1), real(sqrt_2,p_double) * ox * r_wx ) * &
    hermite( this%per_tem_mode(2), real(sqrt_2,p_double) * oy * r_wy ) * &
    exp( - ox2 * r_wx**2 - oy2 * r_wy**2 ) * &
    exp( cmplx(0,1.0_p_double) * (0.5_p_double*k*(ox2*r_Rx + oy2*r_Ry) - gouy_shift) )


  ! case (p_bessel)
  !   ! calculate axial distances
  !   ox = x2 - rx_center
  !   oy = x3 - ry_center

  !   ! calculate Kt * rho
  !   ! currently bessel functions only work in single precision
  !   ktr = real(sqrt(ox**2+oy**2) * this%per_kt, p_single)

  !   ! get kbessel
  !   kbessel = sqrt(k**2 + this%per_kt**2) - abs(k)

  !   ! clip to the specified 0
  !   if (ktr >= this%per_clip_pos) then
  !     env = 0.0_p_double
  !   else
  !     !kbessel = k * sqrt(1.0 - (kt/k)**2), calculated at setup
  !     phi = this%per_n * atan2(oy,ox) + kbessel * (x1 - z_center)
  !     env = tenv * besseljn(this%per_n, ktr) * exp( cmplx(0,1.0_p_double) * phi )
  !   endif

  case (p_gaussian_asym)

    ! THIS DOES NOT WORK WITH CHIRPED PULSES
    ! (there is a possibility of a division by 0 in that case)

    ! calculate transverse coordinates squared
    ox = x2 - rx_center
    oy = x3 - ry_center
    ox2 = ox**2
    oy2 = oy**2

    ! detect which quadrant we are in:
    ! get w0 and focsus as a smooth transition between the 2 regimes
    shp(1) = smooth_heaviside( ox/this%per_asym_trans(1) )
    shp(2) = smooth_heaviside( oy/this%per_asym_trans(2) )

    w0_asym(1) = this%per_w0_asym(1, 1) + shp(1)*( this%per_w0_asym(2, 1) - this%per_w0_asym(1, 1))
    focus_asym(1) = this%per_focus_asym(1,1) + &
    shp(1)*(this%per_focus_asym(2,1) - this%per_focus_asym(1,1))
    w0_asym(2) = this%per_w0_asym(1, 2) + shp(2)*( this%per_w0_asym(2, 2) - this%per_w0_asym(1, 2))
    focus_asym(2) = this%per_focus_asym(1,2) + &
    shp(1)*(this%per_focus_asym(2,2) - this%per_focus_asym(1,2))

    ! get rayleigh ranges
    z0x = k * w0_asym(1)**2 *0.5_p_double ! Rayleigh range for x
    z0y = k * w0_asym(2)**2 *0.5_p_double ! Rayleigh range for y

    ! calculate distance to focal planes
    zx  = x1 - focus_asym(1) ! distance to x focal plane
    zy  = x1 - focus_asym(2) ! distance to y focal plane

    wx = w0_asym(1) *sqrt( 1.0_p_double + (zx/z0x)**2 ) ! local spot, x direction
    wy = w0_asym(2) *sqrt( 1.0_p_double + (zy/z0y)**2 ) ! local spot, y direction

    r_Rx = zx/(zx**2 + z0x**2) ! 1 / ( x curvature radius )
    r_Ry = zy/(zy**2 + z0y**2) ! 1 / ( y curvature radius )

    ! gouy phase shift
    ! Half of the Guoy phase comes from each transverse coord. See
    ! Siegman, Lasers, 16.4 Higher Order Gaussian Modes, Astigmatic
    ! mode functions
    gouy_shift = 0.5_p_double*(atan2( zx, z0x ) + atan2( zy, z0y ))

    ! asymmetric gaussian envelope
    env = tenv * sqrt( w0_asym(1) * w0_asym(2) / (wx * wy)) * &
    exp( - ox2/(wx**2) - oy2/(wy**2) ) * &
    exp( cmplx(0,1.0_p_double) * (0.5_p_double*k*(ox2*r_Rx + oy2*r_Ry) - gouy_shift) )

  case default

    ! This should never happen
    env = 0

  end select


end subroutine
!-----------------------------------------------------------------------------------------



!-----------------------------------------------------------------------------------------
function fast_phase_zpulse(this, tin_sim, x1_sim) result (phasor)
  implicit none
  class( t_zpulse_mov_wall ), intent(inout) :: this
  real(p_double), intent(in) :: tin_sim, x1_sim
  real(p_double) :: x1, t 
  complex(p_double) :: phasor
 
  if( this % if_boost ) then 
    call this % apply_lorentz_transformation( x1_sim, tin_sim, x1, t )
  else 
    x1 = x1_sim 
    t = tin_sim 
 endif
 
 phasor = exp(cmplx(0,1.0_p_double)*( this%prop_sign*this%omega0*t - this%k0*x1 + this%phase0))

end function fast_phase_zpulse
!-----------------------------------------------------------------------------------------



!-----------------------------------------------------------------------------------------
subroutine get_env_1d_zpulse_mov_wall(this, tin_sim, x1_sim, env)
  implicit none

  class ( t_zpulse_mov_wall ), intent(inout) :: this
  real(p_double), intent(in) :: x1_sim, tin_sim
  complex(p_double), intent(out) :: env

  env = this%a0 * this % t_envelope(tin_sim, x1_sim)
end subroutine
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
subroutine launch_wall_1d_zpulse_mov_wall( this, b, e, g_space, nx_p_min, t_sim, dt, no_co )
  implicit none

  ! dummy variables

  integer, parameter :: rank = 1

  class( t_zpulse_mov_wall ), intent(inout) :: this
  type( t_vdf ), intent(inout) :: b, e
  type( t_space ), intent(in) :: g_space
  integer, dimension(:), intent(in) :: nx_p_min
  real(p_double), intent(in) :: t_sim, dt
  class( t_node_conf ), intent(in) :: no_co


  ! local variables

  real(p_double), dimension(2,rank) :: g_x_range
  real(p_double), dimension(rank) :: ldx

  real(p_double) :: zmin, zmax, z, z_2

  real(p_double) :: cos_pol, sin_pol
  real(p_double) :: amp,  gpos
  complex(p_double) :: lenv, lenv_2

  integer :: gipos, ipos, i1

  integer :: lbnd, rbnd, i1start, i1finish
  logical :: lbnd_in_domain, rbnd_in_domain, domain_enclosed, lbnd_right_of_domain, rbnd_left_of_domain

  real(p_double) :: phase0_temp, pol_temp

  ! Global box sizes
  call get_x_bnd( g_space, g_x_range )

  ldx(1) = real( b%dx_(1), p_double )
  zmin = real( g_x_range(p_lower,1), p_double )
  zmax = real( g_x_range(p_upper,1), p_double )

  select case( this%interpolation )
  case( p_linear, p_cubic )
    ! Nothing to change here
  case( p_quadratic, p_quartic )
    zmin = zmin + ldx(1)*0.5_p_double
  end select

  cos_pol = cos( this%pol )
  sin_pol = sin( this%pol )
  amp = this%omega0

  ! start near domain boundary 
  if( this%if_launch_from_wall ) then 
    
    if( this%propagation == p_forward ) then
      gipos = 1  
    else
      gipos = floor( ( zmax - zmin ) / ldx(1) ) - 1
    endif
  
  ! else compute the wall position in global coordinates, stored as gpos 
  else 
    gpos = this%wall_pos + this%wall_vel*t_sim
    gipos = int((gpos - zmin)/ldx(1)) + 1
  endif

  ipos = gipos - nx_p_min(1) + 1

  ! left and right index bounds of the calculation in local coordinates 
  if( this%propagation == p_forward ) then 
    lbnd = ipos 
    rbnd = ipos + this%ncells
  else
    lbnd = ipos - this%ncells 
    rbnd = ipos 
  endif

  lbnd_right_of_domain = ( lbnd >= e%nx_(1) )
  rbnd_left_of_domain = (rbnd <= 1)
  lbnd_in_domain  = (lbnd >= 1) .and. (lbnd  <= e%nx_(1))
  rbnd_in_domain = (rbnd >= 1) .and. (rbnd <= e%nx_(1))
  domain_enclosed  = (lbnd < 1)  .and. (rbnd >  e%nx_(1))


  ! coordinates for local field setting in this domain 
  lbnd = max( lbnd, 1 )
  rbnd = min( rbnd, e%nx_(1) )

  if( lbnd_in_domain .or. rbnd_in_domain .or. domain_enclosed ) then  

    i1start = lbnd
    i1finish = rbnd

    if( .not. this%if_add_fields ) then 
      e%f1(:,i1start:i1finish) = 0 
      b%f1(:,i1start:i1finish) = 0 
    endif

    do i1 = i1start, i1finish 

      z = zmin + real( i1 - 1 + nx_p_min(1), p_double ) * ldx(1) 
      z_2  = z + 0.5_p_double * ldx(1) 
      
      call this % get_env_1d(t_sim, z, lenv)
      call this % get_env_1d(t_sim, z_2, lenv_2)

      lenv = amp * lenv
      lenv_2 = amp * lenv_2
      ! e%f2(1, i1, i2) = 0.0_p_double
      e%f1(2, i1) = e%f1(2, i1) &
              + real( this%boost_per_fac * lenv * this % fast_phase( t_sim, z ) * cos_pol, p_k_fld )
      
      e%f1(3, i1) = e%f1(3, i1) &
              + real( this%boost_per_fac * lenv * this % fast_phase( t_sim, z ) * sin_pol, p_k_fld )

      ! b%f2(1, i1, i2) = 0.0_p_double
      b%f1(2, i1) = b%f1(2, i1) &
              +  real( - this%boost_per_fac * lenv_2 * this % fast_phase( t_sim, z_2 ) * sin_pol * this%prop_sign * this%k0/this%omega0, p_k_fld )
      
      b%f1(3, i1) = b%f1(3, i1) &
              + real( this%boost_per_fac * lenv_2 * this % fast_phase( t_sim, z_2 ) * cos_pol * this%prop_sign * this%k0/this%omega0, p_k_fld )

    enddo 

    if (this%pol_type /= 0) then

      phase0_temp = this%phase0 
      pol_temp    = this%pol

      this%phase0 = this%phase0 + real( sign( pi_2, real( this%pol_type, p_double )), p_double )
      this%pol    = this%pol    + real( pi_2, p_double )

      cos_pol = cos( this%pol )
      sin_pol = sin( this%pol )
      
      do i1 = i1start, i1finish 

        z = zmin + real( i1 - 1 + nx_p_min(1), p_double ) * ldx(1) 
        call this % get_env_1d(t_sim, z, lenv)
        call this % get_env_1d(t_sim, z_2, lenv_2)

        lenv = amp * lenv
        lenv_2 = amp * lenv_2

        e%f1(2, i1) = e%f1(2, i1) &
              + real( this%boost_per_fac * lenv * this % fast_phase( t_sim, z ) * cos_pol, p_k_fld )
      
        e%f1(3, i1) = e%f1(3, i1) &
                + real( this%boost_per_fac * lenv * this % fast_phase( t_sim, z ) * sin_pol, p_k_fld )

        ! b%f2(1, i1, i2) = 0.0_p_double
        b%f1(2, i1) = b%f1(2, i1) &
                +  real( - this%boost_per_fac * lenv_2 * this % fast_phase( t_sim, z_2 ) * sin_pol * this%prop_sign * this%k0/this%omega0, p_k_fld )
        
        b%f1(3, i1) = b%f1(3, i1) &
                + real( this%boost_per_fac * lenv_2 * this % fast_phase( t_sim, z_2 ) * cos_pol * this%prop_sign * this%k0/this%omega0, p_k_fld )

      enddo

      this%phase0 = phase0_temp
      this%pol    = pol_temp       

    endif

  endif 


  if( this%if_clean_fields ) then 
    if (this%propagation==p_forward) then
      if( lbnd_right_of_domain ) then 
        e%f1(:,:) = 0 
        b%f1(:,:) = 0 
      else if( lbnd_in_domain ) then 
        e%f1(:,:i1start-1) = 0 
        b%f1(:,:i1start-1) = 0  
      endif
    else if (this%propagation==p_backward) then
      if( rbnd_left_of_domain ) then 
        e%f1(:,:) = 0 
        b%f1(:,:) = 0 
      else if( rbnd_in_domain ) then 
        e%f1(:,i1finish-1:) = 0 
        b%f1(:,i1finish-1:) = 0  
      endif

    endif 

  endif

end subroutine launch_wall_1d_zpulse_mov_wall
!-----------------------------------------------------------------------------------------



!-----------------------------------------------------------------------------------------
subroutine launch_wall_2d_zpulse_mov_wall( this, b, e, g_space, nx_p_min, t_sim, dt, no_co )

  implicit none

  ! dummy variables

  integer, parameter :: rank = 2

  class( t_zpulse_mov_wall ), intent(inout) :: this
  type( t_vdf ), intent(inout) :: b, e
  type( t_space ), intent(in) :: g_space
  integer, dimension(:), intent(in) :: nx_p_min
  real(p_double), intent(in) :: t_sim, dt
  class( t_node_conf ), intent(in) :: no_co


  ! local variables

  real(p_double), dimension(2,rank) :: g_x_range
  real(p_double), dimension(rank) :: ldx

  real(p_double) :: zmin, zmax, rmin, r, r_2, z_2

  real(p_double) :: cos_pol, sin_pol
  real(p_double) :: amp, gpos
  complex(p_double) :: env_00, env_01, env_10, env_11

  integer :: gipos, ipos, i1, i2
  real(p_double) :: z

  integer :: lbnd, rbnd, i1start, i1finish
  logical :: lbnd_in_domain, rbnd_in_domain, domain_enclosed, lbnd_right_of_domain

  real(p_double) :: phase0_temp, pol_temp
  
  ! Global box sizes
  call get_x_bnd( g_space, g_x_range )

  ldx(1:2) = real( b%dx_(1:2), p_double )
  zmin = real( g_x_range(p_lower,1), p_double )
  zmax = real( g_x_range(p_upper,1), p_double )
  rmin = real( g_x_range(p_lower,2), p_double )

  select case( this%interpolation )
  case( p_linear, p_cubic )
    ! Nothing to change here
  case( p_quadratic, p_quartic )
    zmin = zmin + ldx(1)*0.5_p_double
    rmin = rmin + ldx(2)*0.5_p_double
  end select

  cos_pol = cos( this%pol )
  sin_pol = sin( this%pol )
  amp = this%omega0

  ! start near domain boundary 
  if( this%if_launch_from_wall ) then 
    
    if( this%propagation == p_forward ) then
      gipos = 1  
    else
      gipos = floor( ( zmax - zmin ) / ldx(1) ) - 1 
    endif
  
  ! else compute the wall position in global coordinates, stored as gpos 
  else 
    gpos = this%wall_pos + this%wall_vel*t_sim
    gipos = int((gpos - zmin)/ldx(1)) + 1
  endif

  ipos = gipos - nx_p_min(1) + 1

  ! left and right index bounds of the calculation in local coordinates 
  if( this%propagation == p_forward ) then 
    lbnd = ipos 
    rbnd = ipos + this%ncells
  else
    lbnd = ipos - this%ncells 
    rbnd = ipos 
  endif

  lbnd_right_of_domain = ( lbnd >= e%nx_(1) )
  lbnd_in_domain  = (lbnd >= 1) .and. (lbnd  <= e%nx_(1))
  rbnd_in_domain = (rbnd >= 1) .and. (rbnd <= e%nx_(1))
  domain_enclosed  = (lbnd < 1)  .and. (rbnd >  e%nx_(1))

  ! coordinates for local field setting in this domain 
  lbnd = max( lbnd, 1 )
  rbnd = min( rbnd, e%nx_(1) )

  if( lbnd_in_domain .or. rbnd_in_domain .or. domain_enclosed ) then  

    i1start = lbnd
    i1finish = rbnd

    if( .not. this%if_add_fields ) then 
      e%f2(:,i1start:i1finish,:) = 0 
      b%f2(:,i1start:i1finish,:) = 0 
    endif

    do i2 = lbound(e%f2, 3), ubound(e%f2, 3)
       
      r = rmin + (i2 - 1 + nx_p_min(2)) * ldx(2)
      r_2 = r + 0.5_p_double * ldx(2)

      do i1 = i1start, i1finish 

        z = zmin + real( i1 - 1 + nx_p_min(1), p_double ) * ldx(1) 
        z_2  = z + 0.5_p_double * ldx(1) 
        
        call this % get_env_2d(t_sim, z, r, env_00)
        call this % get_env_2d(t_sim, z, r_2, env_01)
        call this % get_env_2d(t_sim,z_2, r, env_10)
        call this % get_env_2d(t_sim,z_2, r_2, env_11)

        ! e%f2(1, i1, i2) = 0.0_p_double
        e%f2(2, i1, i2) = e%f2(2, i1, i2) &
                + real( amp* this%boost_per_fac * env_01*this%fast_phase(t_sim, z) * cos_pol, p_k_fld )
        
        e%f2(3, i1, i2) = e%f2(3, i1, i2) &
                + real( amp* this%boost_per_fac * env_00*this%fast_phase(t_sim, z) * sin_pol, p_k_fld )

        ! b%f2(1, i1, i2) = 0.0_p_double
        b%f2(2, i1, i2) = b%f2(2, i1, i2) &
                +  real( - amp* this%boost_per_fac * env_10 * this%fast_phase(t_sim, z_2) * sin_pol * this%prop_sign * this%k0/this%omega0, p_k_fld )
        
        b%f2(3, i1, i2) = b%f2(3, i1, i2) &
                + real( amp* this%boost_per_fac * env_11 *this%fast_phase(t_sim, z_2) * cos_pol * this%prop_sign * this%k0/this%omega0, p_k_fld )

      enddo 
    enddo 

    if (this%pol_type /= 0) then

      phase0_temp = this%phase0 
      pol_temp    = this%pol

      this%phase0 = this%phase0 + real( sign( pi_2, real( this%pol_type, p_double )), p_double )
      this%pol    = this%pol    + real( pi_2, p_double )

      cos_pol = cos( this%pol )
      sin_pol = sin( this%pol )
      
      do i2 = lbound(e%f2, 3), ubound(e%f2, 3)
         
        r = rmin + (i2 - 1 + nx_p_min(2)) * ldx(2)
        r_2 = r + 0.5_p_double * ldx(2)

        do i1 = i1start, i1finish 

          z = zmin + real( i1 - 1 + nx_p_min(1), p_double ) * ldx(1) 
          call this % get_env_2d(t_sim, z, r, env_00)
          call this % get_env_2d(t_sim, z, r_2, env_01)
          call this % get_env_2d(t_sim,z_2, r, env_10)
          call this % get_env_2d(t_sim,z_2, r_2, env_11)
        
          ! e%f2(1, i1, i2) = 0.0_p_double
          e%f2(2, i1, i2) = e%f2(2, i1, i2) &
                  + real( amp* this%boost_per_fac * env_01*this%fast_phase(t_sim, z) * cos_pol, p_k_fld )
          
          e%f2(3, i1, i2) = e%f2(3, i1, i2) &
                  + real( amp* this%boost_per_fac * env_00*this%fast_phase(t_sim, z) * sin_pol, p_k_fld )

          ! b%f2(1, i1, i2) = 0.0_p_double
          b%f2(2, i1, i2) = b%f2(2, i1, i2) &
                  +  real( - amp* this%boost_per_fac * env_10 * this%fast_phase(t_sim, z_2) * sin_pol * this%prop_sign * this%k0/this%omega0, p_k_fld )
          
          b%f2(3, i1, i2) = b%f2(3, i1, i2) &
                  + real( amp* this%boost_per_fac * env_11 *this%fast_phase(t_sim, z_2) * cos_pol * this%prop_sign * this%k0/this%omega0, p_k_fld )

        enddo 
      enddo 

      this%phase0 = phase0_temp
      this%pol    = pol_temp       

    endif


  endif 


  if( this%if_clean_fields ) then 

    if( lbnd_right_of_domain ) then 
      e%f2(:,:,:) = 0 
      b%f2(:,:,:) = 0 
    else if( lbnd_in_domain ) then 
      e%f2(:,:i1start-1,:) = 0 
      b%f2(:,:i1start-1,:) = 0  
    endif

  endif

end subroutine launch_wall_2d_zpulse_mov_wall
!-----------------------------------------------------------------------------------------



!-----------------------------------------------------------------------------------------
subroutine launch_wall_3d_zpulse_mov_wall( this, b, e, g_space, nx_p_min, t_sim, dt, no_co )

  implicit none

  ! dummy variables

  integer, parameter :: rank = 3

  class( t_zpulse_mov_wall ), intent(inout) :: this
  type( t_vdf ), intent(inout) :: b, e
  type( t_space ), intent(in) :: g_space
  integer, dimension(:), intent(in) :: nx_p_min
  real(p_double), intent(in) :: t_sim, dt
  class( t_node_conf ), intent(in) :: no_co


  ! local variables

  real(p_double), dimension(2,rank) :: g_x_range
  real(p_double), dimension(rank) :: ldx

  real(p_double) :: zmin, zmax, r1min, r2min, r1, r2, r1_2, r2_2, z_2

  real(p_double) :: cos_pol, sin_pol
  real(p_double) :: amp, gpos
  complex(p_double) :: env

  ! delete 
  complex(p_double) :: q
  real(p_double) :: zr, z0 

  integer :: gipos, ipos, i1, i2, i3
  real(p_double) :: z

  integer :: lbnd, rbnd, i1start, i1finish
  logical :: lbnd_in_domain, rbnd_in_domain, domain_enclosed, lbnd_right_of_domain

  real(p_double) :: phase0_temp, pol_temp

  ! Global box sizes
  call get_x_bnd( g_space, g_x_range )

  ldx(1:3) = real( b%dx_(1:3), p_double )
  zmin = real( g_x_range(p_lower,1), p_double )
  zmax = real( g_x_range(p_upper,1), p_double )
  r1min = real( g_x_range(p_lower,2), p_double )
  r2min = real( g_x_range(p_lower,3), p_double )

  select case( this%interpolation )
  case( p_linear, p_cubic )
    ! Nothing to change here
  case( p_quadratic, p_quartic )
    zmin = zmin + ldx(1)*0.5_p_double
    zmax = zmax + ldx(1)*0.5_p_double
    r1min = r1min + ldx(2)*0.5_p_double
    r2min = r2min + ldx(3)*0.5_p_double
  end select

  cos_pol = cos( this%pol )
  sin_pol = sin( this%pol )
  amp = this%omega0

  ! start near domain boundary 
  if( this%if_launch_from_wall ) then 
    
    if( this%propagation == p_forward ) then
      gipos = 1  
    else
      gipos = floor( ( zmax - zmin ) / ldx(1) ) - 1 
    endif
  
  ! else compute the wall position in global coordinates, stored as gpos 
  else 
    gpos = this%wall_pos + this%wall_vel*t_sim
    gipos = int((gpos - zmin)/ldx(1)) + 1
  endif

  ipos = gipos - nx_p_min(1) + 1

  ! left and right index bounds of the calculation in local coordinates 
  if( this%propagation == p_forward ) then 
    lbnd = ipos 
    rbnd = ipos + this%ncells
  else
    lbnd = ipos - this%ncells 
    rbnd = ipos 
  endif

  lbnd_in_domain  = (lbnd >= 1) .and. (lbnd  <= e%nx_(1))
  rbnd_in_domain = (rbnd >= 1) .and. (rbnd <= e%nx_(1))
  domain_enclosed  = (lbnd < 1)  .and. (rbnd >  e%nx_(1))

  ! coordinates for local field setting in this domain 
  lbnd = max( lbnd, 1 )
  rbnd = min( rbnd, e%nx_(1) )

  if( lbnd_in_domain .or. rbnd_in_domain .or. domain_enclosed ) then  

    i1start = lbnd
    i1finish = rbnd
    
    if( .not. this%if_add_fields ) then 
      e%f3(:,i1start:i1finish,:,:) = 0 
      b%f3(:,i1start:i1finish,:,:) = 0 
    endif

    do i2 = lbound(e%f3, 3), ubound(e%f3, 3)
       
      r1 = r1min + (i2 - 1 + nx_p_min(2)) * ldx(2)
      r1_2 = r1 + 0.5_p_double * ldx(2)

      do i3 = lbound(e%f3, 4), ubound(e%f3, 4)
         
        r2 = r2min + (i3 - 1 + nx_p_min(3)) * ldx(3)
        r2_2 = r2 + 0.5_p_double * ldx(3)

        do i1 = i1start, i1finish 

          z = zmin + real( i1 - 1 + nx_p_min(1), p_double ) * ldx(1) 
          z_2  = z + 0.5_p_double * ldx(1) 
          
          call this % get_env_3d(t_sim, z, r1_2, r2, env)
          e%f3(2, i1, i2, i3) = e%f3(2, i1, i2, i3) + real(amp *  this%boost_per_fac &
                        * env * this % fast_phase( t_sim, z) * cos_pol, p_k_fld )

          call this % get_env_3d(t_sim, z, r1, r2_2, env)
          e%f3(3, i1, i2, i3) = e%f3(3, i1, i2, i3) + real(amp *  this%boost_per_fac &
                        * env * this % fast_phase( t_sim, z ) * sin_pol, p_k_fld )

          call this % get_env_3d(t_sim, z_2, r1, r2_2, env)
          b%f3(2, i1, i2, i3) = b%f3(2, i1, i2, i3) + real( - amp * this%boost_per_fac &
                        *env * this % fast_phase( t_sim, z_2 ) &
                        * sin_pol * this%prop_sign * this%k0/this%omega0, p_k_fld )

          call this % get_env_3d(t_sim, z_2, r1_2, r2 , env)
          b%f3(3, i1, i2, i3) = b%f3(3, i1, i2, i3) + real( amp * this%boost_per_fac &
                        * env * this % fast_phase( t_sim, z_2 ) &
                        * cos_pol * this%prop_sign * this%k0/this%omega0, p_k_fld )

        enddo 
      enddo 
    enddo

   if (this%pol_type /= 0) then

      phase0_temp = this%phase0 
      pol_temp    = this%pol

      this%phase0 = this%phase0 + real( sign( pi_2, real( this%pol_type, p_double )), p_double )
      this%pol    = this%pol    + real( pi_2, p_double )

      cos_pol = cos( this%pol )
      sin_pol = sin( this%pol )
      
      do i2 = lbound(e%f3, 3), ubound(e%f3, 3)
         
        r1 = r1min + (i2 - 1 + nx_p_min(2)) * ldx(2)
        r1_2 = r1 + 0.5_p_double * ldx(2)

        do i3 = lbound(e%f3, 3), ubound(e%f3, 3)
           
          r2 = r2min + (i2 - 1 + nx_p_min(3)) * ldx(3)
          r2_2 = r2 + 0.5_p_double * ldx(3)

          do i1 = i1start, i1finish 

            z = zmin + real( i1 - 1 + nx_p_min(1), p_double ) * ldx(1) 
            call this % get_env_3d(t_sim, z, r1_2, r2, env)
            
            e%f3(2, i1, i2, i3) = e%f3(2, i1, i2, i3) + real(amp *  this%boost_per_fac &
                          * env * this % fast_phase( t_sim, z ) * cos_pol, p_k_fld )

            call this % get_env_3d(t_sim, z, r1, r2_2, env)
            e%f3(3, i1, i2, i3) = e%f3(3, i1, i2, i3) + real(amp *  this%boost_per_fac &
                          * env * this % fast_phase( t_sim, z ) * sin_pol, p_k_fld )

            call this % get_env_3d(t_sim, z_2, r1, r2_2, env)
            b%f3(2, i1, i2, i3) = b%f3(2, i1, i2, i3) + real( - amp * this%boost_per_fac &
                          *env * this % fast_phase( t_sim, z_2 ) &
                          * sin_pol * this%prop_sign * this%k0/this%omega0, p_k_fld )

            call this % get_env_3d(t_sim, z_2, r1_2, r2 , env)
            b%f3(3, i1, i2, i3) = b%f3(3, i1, i2, i3) + real( amp * this%boost_per_fac &
                          * env * this % fast_phase( t_sim, z_2 ) &
                          * cos_pol * this%prop_sign * this%k0/this%omega0, p_k_fld )

          enddo 
        enddo 
      enddo

      this%phase0 = phase0_temp
      this%pol    = pol_temp       

    endif

  endif 

end subroutine launch_wall_3d_zpulse_mov_wall
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
subroutine apply_divcorr_2d_zpulse_mov_wall( this, b, e, g_space, nx_p_min, t, dt, no_co )

  implicit none 
  integer, parameter :: rank = 2

  class( t_zpulse_mov_wall ), intent(inout) :: this
  type( t_vdf ), intent(inout) :: b, e
  type( t_space ), intent(in) :: g_space
  integer, dimension(:), intent(in) :: nx_p_min
  real(p_double), intent(in) :: t, dt
  class( t_node_conf ), intent(in) :: no_co

  real(p_double) :: e2p, e2m, e1, b2p, b2m, b1
  integer :: i1, i2

  integer :: lbnd, rbnd
  logical :: lbnd_in_domain, rbnd_in_domain, domain_enclosed

  integer :: i1start, i1finish, direction, i1_initial_data
  integer :: recv_neighbor, send_neighbor
  logical :: start_in_domain, finish_in_domain

  real(p_double), dimension(2,rank) :: g_x_range
  real(p_double), dimension(rank) :: ldx
  real(p_double) :: zmin, zmax, rmin, z_2
  real(p_double) :: gpos, gipos, ipos 

  integer :: size 
  integer :: ierr
  integer, dimension( MPI_STATUS_SIZE ) :: status
  real(p_double), dimension(:,:), allocatable :: divcorr_buf 

  ! start and finish are the start and finish values of the calculation
  ! for forward pulse, start at the end of the wall in the forward direction
  ! and finish at the start of the pulse 
  ! the field outside the fields we just set (up to ipos + ncells ) is assumed to be correct  

  call get_x_bnd( g_space, g_x_range )

  ldx(1:2) = real( b%dx_(1:2), p_double )
  zmin = real( g_x_range(p_lower,1), p_double )
  zmax = real( g_x_range(p_upper,1), p_double )
  rmin = real( g_x_range(p_lower,2), p_double )

  select case( this%interpolation )
  case( p_linear, p_cubic )
    ! Nothing to change here
  case( p_quadratic, p_quartic )
    zmin = zmin + ldx(1)*0.5_p_double
    rmin = rmin + ldx(2)*0.5_p_double
  end select

  ! left and right bounds of the calculation 
  if( this%if_launch_t0 ) then 

      lbnd = 1 - nx_p_min(1)
      rbnd =  int( (zmax - zmin)/ldx(1) ) - nx_p_min(1)

  else

    ! Current wall position
    gpos = this%wall_pos + this%wall_vel*t
    gipos = int((gpos - zmin)/ldx(1)) + 1 
    ipos = gipos - nx_p_min(1) + 1

    if( (gpos - ldx(1) < zmin) .or. (gpos + (this%ncells + 1) * ldx(1) > zmax ) ) then 
      SCR_ROOT("WARNING: zpulse_mov_wall hit the boundary before pulse injection was complete. Stopping injection.")
      this%if_launch = .false.
      return
    endif

    if( this%propagation == p_forward ) then 
      lbnd = ipos 
      rbnd = ipos + this%ncells + 1 
    else
      lbnd = ipos - this%ncells - 1 
      rbnd = ipos 
    endif

  endif

  lbnd_in_domain  = ( lbnd >= 1) .and. (lbnd  <= e%nx_(1))
  rbnd_in_domain = (rbnd >= 1) .and. (rbnd <= e%nx_(1))
  domain_enclosed  = (lbnd < 1)  .and. (rbnd >  e%nx_(1))

  ! wall outside is completely outside this domain: do nothing 
  if( .not. this%if_launch_t0 .and. .not.( lbnd_in_domain .or. rbnd_in_domain .or. domain_enclosed ) ) return 

  if( this%propagation == p_forward ) then 
    direction = -1
    recv_neighbor = p_upper 
    send_neighbor = p_lower  
    i1start = min( rbnd - 1, e%nx_(1) ) 
    i1finish = max( lbnd, 1 )  
    i1_initial_data = rbnd
    start_in_domain = rbnd_in_domain
    finish_in_domain = lbnd_in_domain
  else 
    direction = 1
    recv_neighbor = p_lower 
    send_neighbor = p_upper
    i1start = max( 1, lbnd + 1 ) 
    i1finish = min( rbnd, e%nx_(1) )  
    i1_initial_data = lbnd
    start_in_domain = lbnd_in_domain
    finish_in_domain = rbnd_in_domain 
  endif
  
  ! print *, i1start, i1finish, e%nx_(1), lbnd_in_domain, rbnd_in_domain, domain_enclosed

  size = e%nx_(2)
  allocate( divcorr_buf( 2, size ) )

  ! check if integral start is in domain; else we need calculated values 
  ! from nearest neighbor  
  if( .not. start_in_domain ) then 
    ! print *, "calling recv."
    call mpi_recv( divcorr_buf, 2 * size, MPI_DOUBLE, &
                   no_co%neighbor( recv_neighbor, 1 ) - 1, &
                   0, no_co%comm, status, ierr )
    
    ! print *, "recv"

    CHECK_ERROR( ierr, "zpulse wall: recv fail", p_err_mpi)

  endif 

  do i2 = 1, e%nx_(2) 

    ! starting value: either use value from previous node if the 
    ! divergence integral is not starting locally 
    if( start_in_domain ) then 
      e1 = e%f2( 1, i1_initial_data, i2 )
      b1 = b%f2( 1, i1_initial_data, i2 )
    else
      e1 = divcorr_buf(1,i2)
      b1 = divcorr_buf(2,i2)
    endif 

    do i1 = i1start, i1finish, direction   

      ! old method: 
      ! e2p = e%f2(2,i1+1,i2)
      ! e2m = e%f2(2,i1+1,i2-1)

      e2p = e%f2(2,i1,i2)
      e2m = e%f2(2,i1,i2-1)

      e1 = e1 - direction * ldx(1)*(e2p - e2m)/ldx(2)
      e%f2(1, i1, i2) = real( e1, p_k_fld )

      b2p = b%f2(2,i1,i2+1)
      b2m = b%f2(2,i1,i2)

      b1 = b1 - direction * ldx(1) * (b2p - b2m) / ldx(2) 
      b%f2(1, i1, i2) = real(b1, p_k_fld )

    enddo ! i1

  enddo ! i2

  ! stop 

  ! if integral finish is not in domain, we need to send the results of the integral 
  ! to be finished in the next domain  
  if( .not. finish_in_domain ) then 
  
    divcorr_buf(1,1:size) = e%f2(1, i1finish, 1:size)
    divcorr_buf(2,1:size) = b%f2(1, i1finish, 1:size)

    call mpi_send( divcorr_buf, 2 * size, MPI_DOUBLE, &
                   no_co%neighbor( send_neighbor, 1 ) - 1, &
                   0, no_co%comm, ierr )
  
    CHECK_ERROR( ierr, "zpulse wall: send fail", p_err_mpi)

  endif 
 
 deallocate( divcorr_buf )

end subroutine apply_divcorr_2d_zpulse_mov_wall
!-----------------------------------------------------------------------------------------



!-----------------------------------------------------------------------------------------
subroutine apply_divcorr_3d_zpulse_mov_wall( this, b, e, g_space, nx_p_min, t, dt, no_co )

  implicit none 
  integer, parameter :: rank = 2

  class( t_zpulse_mov_wall ), intent(inout) :: this
  type( t_vdf ), intent(inout) :: b, e
  type( t_space ), intent(in) :: g_space
  integer, dimension(:), intent(in) :: nx_p_min
  real(p_double), intent(in) :: t, dt
  class( t_node_conf ), intent(in) :: no_co

 !  real(p_double) :: e2p, e2m, e1, b2p, b2m, b1
 !  integer :: i1, i2

 !  integer :: lbnd, rbnd
 !  logical :: lbnd_in_domain, rbnd_in_domain, domain_enclosed

 !  integer :: i1start, i1finish, direction, i1_initial_data
 !  integer :: recv_neighbor, send_neighbor
 !  logical :: start_in_domain, finish_in_domain

 !  real(p_double), dimension(2,rank) :: g_x_range
 !  real(p_double), dimension(rank) :: ldx
 !  real(p_double) :: zmin, zmax, rmin
 !  real(p_double) :: gpos, gipos, ipos 

 !  integer :: size 
 !  integer :: ierr
 !  integer, dimension( MPI_STATUS_SIZE ) :: status
 !  real(p_double), dimension(:,:), allocatable :: divcorr_buf 

 !  ! start and finish are the start and finish values of the calculation
 !  ! for forward pulse, start at the end of the wall in the forward direction
 !  ! and finish at the start of the pulse 
 !  ! the field outside the fields we just set (up to ipos + ncells ) is assumed to be correct  

 !  ! Global box sizes
 !  call get_x_bnd( g_space, g_x_range )

 !  ldx(1:2) = real( b%dx_(1:2), p_double )
 !  zmin = real( g_x_range(p_lower,1), p_double )
 !  zmax = real( g_x_range(p_upper,1), p_double )
 !  rmin = real( g_x_range(p_lower,2), p_double )

 !  select case( this%interpolation )
 !  case( p_linear, p_cubic )
 !    ! Nothing to change here
 !  case( p_quadratic, p_quartic )
 !    zmin = zmin + ldx(1)*0.5_p_double
 !    rmin = rmin + ldx(2)*0.5_p_double
 !  end select

 !  ! start near domain boundary 
 !  if( this%if_launch_from_wall ) then 
    
 !    if( this%propagation == p_forward ) then
 !      gipos = 1  
 !    else
 !      gipos = floor( ( zmax - zmin ) / ldx(1) ) - 1 
 !    endif
  
 !  ! else compute the wall position in global coordinates, stored as gpos 
 !  else 
 !    gpos = this%wall_pos + this%wall_vel*t
 !    gipos = int((gpos - zmin)/ldx(1)) + 1
  
 !    if( (gpos - ldx(1) < zmin) .or. (gpos + (this%ncells + 1) * ldx(1) > zmax ) ) then 
 !      SCR_ROOT("WARNING: zpulse_mov_wall hit the boundary before pulse injection was complete. Stopping injection.")
 !      this%if_launch = .false.
 !      return
 !    endif
  
 !  endif

 !  ipos = gipos - nx_p_min(1) + 1

 !  ! left and right bounds of the calculation 
 ! if( this%propagation == p_forward ) then 
 !   lbnd  = ipos 
 !   rbnd = ipos + this%ncells + 1 
 !  else
 !    lbnd = ipos - this%ncells - 1 
 !    rbnd = ipos 
 !  endif

 !  lbnd_in_domain  = ( lbnd >= 1) .and. (lbnd  <= e%nx_(1))
 !  rbnd_in_domain = (rbnd >= 1) .and. (rbnd <= e%nx_(1))
 !  domain_enclosed  = (lbnd < 1)  .and. (rbnd >  e%nx_(1))


 !  ! wall outside is completely outside this domain: do nothing 
 !  if( .not.( lbnd_in_domain .or. rbnd_in_domain .or. domain_enclosed ) ) return 

 !  if( this%propagation == p_forward ) then 
 !    direction = -1
 !    recv_neighbor = p_upper 
 !    send_neighbor = p_lower  
 !    i1start = min( rbnd - 1, e%nx_(1) ) 
 !    i1finish = max( lbnd, 1 )  
 !    i1_initial_data = rbnd
 !    start_in_domain = rbnd_in_domain
 !    finish_in_domain = lbnd_in_domain
 !  else 
 !    direction = 1
 !    recv_neighbor = p_lower 
 !    send_neighbor = p_upper
 !    i1start = max( 1, lbnd + 1 ) 
 !    i1finish = min( rbnd, e%nx_(1) )  
 !    i1_initial_data = lbnd
 !    start_in_domain = lbnd_in_domain
 !    finish_in_domain = rbnd_in_domain 
 !  endif
  
  ! print *, i1start, i1finish, e%nx_(1), lbnd_in_domain, rbnd_in_domain, domain_enclosed

  ! size = e%nx_(2)
  ! allocate( divcorr_buf( 2, size ) )

  ! ! check if integral start is in domain; else we need calculated values 
  ! ! from nearest neighbor  
  ! if( .not. start_in_domain ) then 
  !   ! print *, "calling recv."
  !   call mpi_recv( divcorr_buf, 2 * size, MPI_DOUBLE, &
  !                  no_co%neighbor( recv_neighbor, 1 ) - 1, &
  !                  0, no_co%comm, status, ierr )
   
  !   CHECK_ERROR( ierr, "zpulse wall: recv fail", p_err_mpi)

  ! endif 

  ! do i3 = 1, e%nx_(3) 
  
  !   do i2 = 1, e%nx_(2) 

  !     ! starting value: either use value from previous node if the 
  !     ! divergence integral is not starting locally 
  !     if( start_in_domain ) then 
  !       e1 = e%f3( 1, i1_initial_data, i2 )
  !       b1 = b%f3( 1, i1_initial_data, i2 )
  !     else
  !       e1 = divcorr_buf(1,i2,i3)
  !       b1 = divcorr_buf(2,i2,i3)
  !     endif 

  !     do i1 = i1start, i1finish, direction   

  !       ! old method: 
  !       ! e2p = e%f2(2,i1+1,i2)
  !       ! e2m = e%f2(2,i1+1,i2-1)

  !       e2p = e%f2(2,i1,i2)
  !       e2m = e%f2(2,i1,i2-1)

  !       e1 = e1 - direction * ldx(1)*(e2p - e2m)/ldx(2)
  !       e%f2(1, i1, i2) = real( e1, p_k_fld )

  !       b2p = b%f2(2,i1,i2+1)
  !       b2m = b%f2(2,i1,i2)

  !       b1 = b1 - direction * ldx(1) * (b2p - b2m) / ldx(2) 
  !       b%f2(1, i1, i2) = real(b1, p_k_fld )

  !     enddo ! i1

  !   enddo ! i2
  ! enddo ! i3

  ! ! stop 

  ! ! if integral finish is not in domain, we need to send the results of the integral 
  ! ! to be finished in the next domain  
  ! if( .not. finish_in_domain ) then 
  
  !   divcorr_buf(1,:) = e%f2(1, i1finish, 1:size)
  !   divcorr_buf(2,:) = b%f2(1, i1finish, 1:size)

  !   call mpi_send( divcorr_buf, 2 * size, MPI_DOUBLE, &
  !                  no_co%neighbor( send_neighbor, 1 ) - 1, &
  !                  0, no_co%comm, ierr )
  
  !   CHECK_ERROR( ierr, "zpulse wall: send fail", p_err_mpi)

  ! endif 
 
 ! deallocate( divcorr_buf )

end subroutine apply_divcorr_3d_zpulse_mov_wall
!-----------------------------------------------------------------------------------------



! !-----------------------------------------------------------------------------------------
! subroutine clean_fields_2d_zpulse_mov_wall( this, b, e, g_space, nx_p_min, t, dt, no_co )

!   implicit none 
!   integer, parameter :: rank = 2

!   class( t_zpulse_mov_wall ), intent(inout) :: this
!   type( t_vdf ), intent(inout) :: b, e
!   type( t_space ), intent(in) :: g_space
!   integer, dimension(:), intent(in) :: nx_p_min
!   real(p_double), intent(in) :: t, dt
!   class( t_node_conf ), intent(in) :: no_co

!   call this % get_wall_vars_2d( g_space, lbnd, rbnd )

!   lbnd_in_domain  = (lbnd >= 1) .and. (lbnd  <= e%nx_(1))
!   rbnd_in_domain = (rbnd >= 1) .and. (rbnd <= e%nx_(1))
!   domain_enclosed  = (lbnd < 1)  .and. (rbnd >  e%nx_(1))

!   ! coordinates for local field setting in this domain 
!   lbnd = max( lbnd, 1 )
!   rbnd = min( rbnd, e%nx_(1) )


! end subroutine clean_fields_2d_zpulse_mov_wall
! !-----------------------------------------------------------------------------------------


! !-----------------------------------------------------------------------------------------
! subroutine get_wall_vars_2d_zpulse_mov_wall( this, g_space, nx_p_min, t, &
!                                               zmin, rmin )

!   implicit none 
!   integer, parameter :: rank = 2

!   class( t_zpulse_mov_wall ), intent(inout) :: this
!   type( t_vdf ), intent(inout) :: b, e
!   type( t_space ), intent(in) :: g_space
!   integer, dimension(:), intent(in) :: nx_p_min
!   real(p_double), intent(in) :: t

!   ! Global box sizes
!   call get_x_bnd( g_space, g_x_range )

!   ldx(1:2) = real( b%dx_(1:2), p_double )
!   zmin = real( g_x_range(p_lower,1), p_double )
!   zmax = real( g_x_range(p_upper,1), p_double )
!   rmin = real( g_x_range(p_lower,2), p_double )

!   select case( this%interpolation )
!   case( p_linear, p_cubic )
!     ! Nothing to change here
!   case( p_quadratic, p_quartic )
!     zmin = zmin + ldx(1)*0.5_p_double
!     rmin = rmin + ldx(2)*0.5_p_double
!   end select

!   ! start near domain boundary 
!   if( this%if_launch_from_wall ) then 
    
!     if( this%propagation == p_forward ) then
!       gipos = 1  
!     else
!       gipos = floor( ( zmax - zmin ) / ldx(1) ) - 1 
!     endif
  
!   ! else compute the wall position in global coordinates, stored as gpos 
!   else 
!     gpos = this%wall_pos + this%wall_vel*t
!     gipos = int((gpos - zmin)/ldx(1)) + 1
!   endif

!   ipos = gipos - nx_p_min(1) + 1

!   ! left and right index bounds of the calculation in local coordinates 
!   if( this%propagation == p_forward ) then 
!     lbnd = ipos 
!     rbnd = ipos + this%ncells
!   else
!     lbnd = ipos - this%ncells 
!     rbnd = ipos 
!   endif

! end subroutine get_wall_vars_2d_zpulse_mov_wall
! !-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
subroutine launch_zpulse_mov_wall( this, emf, g_space, &
                                   nx_p_min, g_nx, t, dt, no_co )

  implicit none

  class( t_zpulse_mov_wall ), intent(inout)     :: this
  class( t_emf ) ,  intent(inout) :: emf

  type( t_space ), intent(in) :: g_space
  integer, intent(in), dimension(:) :: nx_p_min, g_nx
  real(p_double), intent(in) :: t
  real(p_double), intent(in) :: dt
  class( t_node_conf ), intent(in) :: no_co

  real(p_double) :: phase0_temp, pol_temp

  if( .not. this%if_launch ) return 
  
  ! Turn off antenna if pulse has ended
  if ( t > this%launch_time + this % t_duration() ) then
    SCR_ROOT( "duration exceeded. ")
    SCR_ROOT( "this%t_duration(): ", this%t_duration() )
    this%if_launch = .false.
    return
  endif
  
  select case (emf%b%x_dim_)
  case(1)
     call this % launch_wall_1d( emf%b, emf%e, g_space, nx_p_min, t, dt, no_co )
  case(2)
     call this % launch_wall_2d( emf%b, emf%e, g_space, nx_p_min, t, dt, no_co )

     if( this%if_apply_divcorr ) then 
        call this % apply_divcorr_2d( emf%b, emf%e, g_space, nx_p_min, t, dt, no_co ) 
     endif

  case(3)
     call this % launch_wall_3d( emf%b, emf%e, g_space, nx_p_min, t, dt, no_co )

     if( this%if_apply_divcorr ) then 
        call this % apply_divcorr_3d( emf%b, emf%e, g_space, nx_p_min, t, dt, no_co ) 
     endif

  case default
     ERROR('Not implemented yet')
     call abort_program( p_err_notimplemented )

  end select


  ! if ( this%if_clean_fields .and. ( t >= this%clean_fields_time ) ) then 
  !   call this % clean_fields( emf%b, emf%e, g_space, nx_p_min, t, dt )
  !   this % if_clean_fields = .false.
  ! endif 

end subroutine launch_zpulse_mov_wall
!-----------------------------------------------------------------------------------------



! !-----------------------------------------------------------------------------------------
! subroutine clean_fields_zpulse_mov_wall(  this, b, e, g_space, nx_p_min, t, dt )

!   implicit none

!   integer, parameter :: rank = 2

!   class( t_zpulse_mov_wall ), intent(inout) :: this
!   type( t_vdf ), intent(inout) :: b, e
!   type( t_space ), intent(in) :: g_space
!   integer, dimension(:), intent(in) :: nx_p_min
!   real(p_double), intent(in) :: t, dt


!   ! local variables
!   real(p_double), dimension(2,rank) :: g_x_range
!   real(p_double), dimension(rank) :: ldx

!   integer :: gipos, ipos, i1
!   real(p_double) :: zmin, gpos 

!   print *, "INFO: clearing fields."

!   call get_x_bnd( g_space, g_x_range )

!   ldx(1:2) = real( b%dx_(1:2), p_double )
!   zmin = real( g_x_range(p_lower,1), p_double )

!   select case( this%interpolation )
!   case( p_linear, p_cubic )
!     ! Nothing to change here
!   case( p_quadratic, p_quartic )
!     zmin = zmin + ldx(1)*0.5_p_double
!   end select

!   ! original wall position
!   gpos = this%wall_pos ! + this%wall_vel*t
!   gipos = int((gpos - zmin)/ldx(1)) + 1 
!   ipos = gipos - nx_p_min(1) + 1

!   do i1 = 1, ipos + this % ncells 
        
!       e%f2(1, i1, :) = 0.0_p_double
!       e%f2(2, i1, :) = 0.0_p_double
!       e%f2(3, i1, :) = 0.0_p_double

!       b%f2(1, i1, :) = 0.0_p_double
!       b%f2(2, i1, :) = 0.0_p_double
!       b%f2(3, i1, :) = 0.0_p_double

!   enddo

! end subroutine clean_fields_zpulse_mov_wall
! !-----------------------------------------------------------------------------------------





end module m_zpulse_mov_wall
