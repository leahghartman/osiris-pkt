# 1 "cyl_modes/os-zpulse-mov-wall-cyl-modes.f03"
# 1 "<built-in>" 1
# 1 "<built-in>" 3
# 467 "<built-in>" 3
# 1 "<command line>" 1
# 1 "<built-in>" 2
# 1 "cyl_modes/os-zpulse-mov-wall-cyl-modes.f03" 2
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
# 2 "cyl_modes/os-zpulse-mov-wall-cyl-modes.f03" 2
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
# 3 "cyl_modes/os-zpulse-mov-wall-cyl-modes.f03" 2

module m_zpulse_mov_wall_cyl_modes

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
# 7 "cyl_modes/os-zpulse-mov-wall-cyl-modes.f03" 2

use m_parameters
use m_math
use m_fparser, only: p_max_expr_len, p_k_fparse, eval
use m_node_conf, only: t_node_conf
use m_zpulse_mov_wall, only : t_zpulse_mov_wall
use m_input_file, only: t_input_file, get_namelist, p_max_nml_section_name
use m_emf_define, only: t_emf_bound, t_emf
use m_emf_cyl_modes, only: t_emf_cyl_modes
use m_space, only: t_space, x_bnd, get_x_bnd
use m_vdf_define, only: t_vdf
use m_grid_define, only: t_grid
use m_grid_cyl_modes,only: t_grid_cyl_modes
use m_space, only: xmin, xmax
use m_cyl_modes, only: t_cyl_modes, new_copy, add, cleanup

use m_zpulse_std, only: p_zpulse_bnormal
use m_zpulse_std, only: p_max_chirp_order, p_polynomial, p_sin2, p_gaussian, p_func
use m_zpulse_std, only: p_const
use m_zpulse_std, only: p_hermite_gaussian, p_hermite_gaussian_astigmatic
use m_zpulse_std, only: p_laguerre_gaussian, p_plane, p_gaussian_asym, t_zpulse, p_forward, p_plane


!-----------------------------------------------------------------------------------------
type, extends( t_zpulse_mov_wall ) :: t_zpulse_mov_wall_cyl_modes
  logical :: if_correct_axis_fields
contains
  procedure :: read_input => read_input_zpulse_mov_wall_cyl_modes
  procedure :: launch => launch_zpulse_mov_wall_cyl_modes
  procedure :: set_fields => set_fields_zpulse_mov_wall_cyl_modes
  procedure :: apply_divcorr => apply_divcorr_zpulse_mov_wall_cyl_modes
  procedure :: check_dimensionality => check_dimensionality_zpulse_mov_wall_cyl_modes
  procedure :: correct_axis_fields => correct_axis_fields_zpulse_mov_wall_cyl_modes

  procedure :: per_envelope => per_envelope_zpulse_mov_wall_cyl_modes
  procedure :: per_envelope_laguerre => per_envelope_laguerre_zpulse_mov_wall_cyl_modes
  procedure :: per_envelope_gaussian => per_envelope_gaussian_zpulse_mov_wall_cyl_modes
  procedure :: t_envelope => t_env_zpulse_mov_wall_cyl_modes

end type t_zpulse_mov_wall_cyl_modes
!-----------------------------------------------------------------------------------------

contains

!-----------------------------------------------------------------------------------------
subroutine read_input_zpulse_mov_wall_cyl_modes( this, input_file, g_space, bnd_con, periodic, grid, sim_options )

  implicit none

  class( t_zpulse_mov_wall_cyl_modes ), intent(inout) :: this
  class( t_input_file ), intent(inout) :: input_file
  type( t_space ), intent(in) :: g_space
  class (t_emf_bound), intent(in) :: bnd_con
  logical, dimension(:), intent(in) :: periodic
  class( t_grid ), intent(in) :: grid
  type( t_options ), intent(in) :: sim_options

  this%if_correct_axis_fields = .true.

  ! Call superclass method
  call this % t_zpulse_mov_wall % read_input( input_file, g_space, bnd_con, periodic, grid, sim_options )

  ! Only allow propagation in dimension 1
  if ( this%direction > 1 ) then
    if ( mpi_node() == 0 ) then
      write(0,*) ""
      write(0,*) "   Error reading zpulse parameters"
      write(0,*) "   propagation only supported in direction 1"
      write(0,*) "   when running with cylindrical modes"
      write(0,*) "   aborting..."
    endif
    stop
  endif

  ! Only allow normal b_type
  if ( this%b_type /= p_zpulse_bnormal ) then
    if ( mpi_node() == 0 ) then
      write(0,*) ""
      write(0,*) "   Error reading zpulse parameters"
      write(0,*) "   btype must be 'normal' when running with cyl_modes"
      write(0,*) "   aborting..."
    endif
    stop
  endif

  ! Don't allow tilted pulses for now
  if ( this%if_tenv_tilt ) then
    if ( mpi_node() == 0 ) then
      write(0,*) ""
      write(0,*) "   Error reading zpulse parameters"
      write(0,*) "   tilt not implemented when running with cyl_modes"
      write(0,*) "   aborting..."
    endif
    stop
  endif

  ! set per_center to be on axis (otherwise the pulse will be invalid)
  ! if the default per_center was set in regular zpulse read_input, then
  ! reset the default per_center to be 0
  if( this%per_center(1) == 0.5_p_double * ( xmin( g_space, 2 ) + xmax( g_space, 2 ) ) ) then
   ! set per_center to be on axis (otherwise the pulse will be invalid)
    this % per_center(1:2) = 0.0_p_double
  else
    if( this%per_center(1) .ne. 0 ) then
      if (mpi_node()==0) print *,"   ERROR: per_center(1) not equal to 0 was requested."
      if (mpi_node()==0) print *,"   When running in cyl_modes, per_center(1) must be 0,"
      if (mpi_node()==0) print *,"   otherwise the pulse will have an infinite number of non-negligible "
      if (mpi_node()==0) print *,"   angular harmonics (instead of just m=1 when per_center = 0 for "
      if (mpi_node()==0) print *,"   cylindrically symmetric pulses). 0 is the default. "
      stop
    endif
  endif

  if( this%per_type .eq. p_hermite_gaussian ) then
    if( this%per_tem_mode(1) > 0 .or. abs(this%per_tem_mode(2)) > 0 ) then
      if (mpi_node()==0) print *,"   ERROR: hermite transverse profiles do not work in quasi-3D."
      stop
    endif
  endif

  if( this%if_launch_from_wall ) then
    if (mpi_node()==0) print *,"   ERROR: if_launch_from_wall is not currently supported in quasi-3D."
    stop
  endif

end subroutine read_input_zpulse_mov_wall_cyl_modes
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! See comment above check_dimensionality_zpulse_cyl_modes for explanation
!-----------------------------------------------------------------------------------------
subroutine check_dimensionality_zpulse_mov_wall_cyl_modes( this )

  class( t_zpulse_mov_wall_cyl_modes ), intent(inout) :: this

  ! do nothing

end subroutine check_dimensionality_zpulse_mov_wall_cyl_modes
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
subroutine launch_zpulse_mov_wall_cyl_modes( this, emf, g_space, nx_p_min, g_nx, t, dt, &
                                             no_co )

  implicit none

  class( t_zpulse_mov_wall_cyl_modes ), intent(inout) :: this
  class( t_emf ), intent(inout) :: emf
  type( t_space ), intent(in) :: g_space
  integer, intent(in), dimension(:) :: nx_p_min, g_nx
  real(p_double), intent(in) :: t
  real(p_double), intent(in) :: dt
  class( t_node_conf ), intent(in) :: no_co

  integer :: n_modes
  integer :: ell

  select type( emf )
  class is( t_emf_cyl_modes )

  ell = abs( this % per_tem_mode(2) )

  n_modes = ubound( emf%b_cyl_m%pf_re, 1 )

  if( this % per_type .eq. p_laguerre_gaussian ) then
    if( n_modes < ell + 1 ) then
      if ( mpi_node() == 0 ) then
        write(0,*) ""
        write(0,*) "    ERROR: num cyl modes must be >= ell + 1 for LG modes"
        write(0,*) "   n_modes: ", n_modes
        write(0,*) "   ell: ", ell
      endif
      stop
    endif
  endif

  if ( this%if_launch .and. ( t >= this%launch_time ) ) then

    ! Turn off antenna if pulse has ended
    if ( t > this % t_duration() + this%launch_time ) then
      if (mpi_node()==0) print *,"duration exceeded. "
      if (mpi_node()==0) print *,"this%t_duration(): ", this%t_duration()
      this%if_launch = .false.
      return
    endif

    call this % set_fields( emf%e_cyl_m, emf%b_cyl_m, g_space, nx_p_min, g_nx, t, dt, &
                            no_co )

    ! If wall has hit the boundary, then do not perform divergence or axis correction
    if ( .not. this%if_launch ) return

    if( this%per_type .ne. p_laguerre_gaussian .or. ell .eq. 0 ) then

      if( this%if_apply_divcorr ) then
        call this % apply_divcorr( emf%e_cyl_m, emf%b_cyl_m, g_space, emf%gix_pos(2), &
                                   nx_p_min, g_nx, t, dt, no_co, 1 )
      endif

      if( this%if_correct_axis_fields ) then
        call this % correct_axis_fields( emf%e_cyl_m, emf%b_cyl_m, g_space, &
                                         emf%gix_pos(2), nx_p_min, g_nx, t, dt, no_co, 1 )
      endif

    else ! handle LG modes

      if( this%if_apply_divcorr ) then
        call this % apply_divcorr( emf%e_cyl_m, emf%b_cyl_m, g_space, emf%gix_pos(2), &
                                   nx_p_min, g_nx, t, dt, no_co, ell + 1 )
        call this % apply_divcorr( emf%e_cyl_m, emf%b_cyl_m, g_space, emf%gix_pos(2), &
                                   nx_p_min, g_nx, t, dt, no_co, ell - 1 )
      endif

      if( this%if_correct_axis_fields ) then
        call this % correct_axis_fields( emf%e_cyl_m, emf%b_cyl_m, g_space, &
                                         emf%gix_pos(2), nx_p_min, g_nx, t, dt, no_co, &
                                         ell + 1 )
        call this % correct_axis_fields( emf%e_cyl_m, emf%b_cyl_m, g_space, &
                                         emf%gix_pos(2), nx_p_min, g_nx, t, dt, no_co, &
                                         ell - 1 )
      endif

    endif

  endif ! this%if_launch .and. ( t >= this%launch_time )

end select

end subroutine launch_zpulse_mov_wall_cyl_modes
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
! B = zhat \times E ---> B_phi = prop_sign * E_r, B_r = - prop_sign * E_phi
!-----------------------------------------------------------------------------------------
subroutine set_fields_zpulse_mov_wall_cyl_modes( this, e_cyl_m, b_cyl_m, g_space, &
                                                 nx_p_min, g_nx, t_sim, dt, no_co )

  implicit none

  class( t_zpulse_mov_wall_cyl_modes ), intent(inout) :: this
  type( t_cyl_modes ), intent(inout) :: e_cyl_m, b_cyl_m
  type( t_space ), intent(in) :: g_space
  integer, intent(in), dimension(:) :: nx_p_min, g_nx
  real(p_double), intent(in) :: t_sim
  real(p_double), intent(in) :: dt
  class( t_node_conf ), intent(in) :: no_co

  real(p_double) :: amp, lenv, lenv_2, zmin, zmax, rmin, z, z_2, r, r_2, dr_2, dz_2, t, t_2
  real(p_double) :: z_sim, z_sim_2
  real(p_double) :: gpos, prop_sign, fast_phase, fast_phase_2

  integer :: ipos, gipos

  integer :: lbnd, rbnd
  logical :: lbnd_in_domain, rbnd_in_domain, domain_enclosed

  integer :: i1, i2, n
  integer :: i1start, i1finish

  integer, parameter :: rank = 2
  real(p_double), dimension(2,rank) :: g_x_range

  real(p_double), dimension(rank) :: ldx

  integer :: ell, sign_ell, ml, mu
  real(p_double) :: per_env_amp, per_env_phase, total_phase

  type(t_vdf), pointer :: e_re0, e_re_l, e_im_l, e_re_u, e_im_u
  type(t_vdf), pointer :: b_re_l, b_im_l, b_re_u, b_im_u

  e_re0 => e_cyl_m%pf_re(0)
  ldx(1) = e_re0%dx(1)
  ldx(2) = e_re0%dx(2)

  ! get zmin, zmax, rmin in absolute coordinates
  call get_x_bnd( g_space, g_x_range )
  zmin = real( g_x_range(p_lower,1), p_double )
  zmax = real( g_x_range(p_upper,1), p_double )
  rmin = real( g_x_range(p_lower,2), p_double )

  amp = this%omega0 * this%a0 * this%boost_per_fac

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
       gipos = e_re0%lbound(2)
    else
      gipos = floor( ( zmax - zmin ) / ldx(1) ) - 1
    endif

  ! else compute the wall position in global coordinates, stored as gpos
  else
    gpos = this%wall_pos + this%wall_vel * t_sim
    gipos = int((gpos - zmin)/ldx(1)) + 1

    if( (gpos - ldx(1) < zmin) .or. (gpos + (this%ncells + 1) * ldx(1) > zmax ) ) then
      if (mpi_node()==0) print *,"WARNING: zpulse_mov_wall hit the boundary before pulse injection"
      if (mpi_node()==0) print *,"was complete. Stopping injection."
      this%if_launch = .false.
      return
    endif

  endif

  ipos = gipos - nx_p_min(1) + 1

  ! switch to local coordinates
  zmin = zmin + real(nx_p_min(1)-1, p_double)*ldx(1)
  rmin = rmin + real(nx_p_min(2)-1, p_double)*ldx(2)

  ! left and right index bounds of the calculation in local coordinates
  if( this%propagation == p_forward ) then
    lbnd = ipos
    rbnd = ipos + this%ncells
  else
    lbnd = ipos - this%ncells
    rbnd = ipos
  endif

  lbnd_in_domain = (lbnd >= 1) .and. (lbnd <= e_re0%nx_(1))
  rbnd_in_domain = (rbnd >= 1) .and. (rbnd <= e_re0%nx_(1))
  domain_enclosed = (lbnd < 1) .and. (rbnd > e_re0%nx_(1))

  lbnd = max( lbnd, 1 )
  rbnd = min( rbnd, e_re0%nx_(1) )

  if( .not.( lbnd_in_domain .or. rbnd_in_domain .or. domain_enclosed ) ) return

  i1start = lbnd
  i1finish = rbnd

  ell = abs( this % per_tem_mode(2) )
  ! Lower and upper modes to which we add fields. Normally this is just ell - 1 and ell + 1,
  ! but we support negative ell values here
  if ( this % per_tem_mode(2) < 0 ) then
    sign_ell = -1
    ml = ell + 1 ! equivalent to abs( this % per_tem_mode(2) - 1 )
    mu = ell - 1 ! equivalent to abs( this % per_tem_mode(2) + 1 )
  else
    sign_ell = +1
    ml = ell - 1
    mu = ell + 1
  endif

  if( .not. this%if_add_fields ) then

    ! The Ez component is forward from transverse components, so don't set the last one
    ! to zero. We will need this value to read in during div_corr.

    e_cyl_m%pf_re(0)%f2(1,i1start:i1finish-1,:) = 0
    b_cyl_m%pf_re(0)%f2(1,i1start:i1finish-1,:) = 0
    e_cyl_m%pf_re(0)%f2(2:3,i1start:i1finish,:) = 0
    b_cyl_m%pf_re(0)%f2(2:3,i1start:i1finish,:) = 0

    do n = 1, ubound( e_cyl_m%pf_re, 1 )

      e_cyl_m%pf_re(n)%f2(1,i1start:i1finish-1,:) = 0
      e_cyl_m%pf_im(n)%f2(1,i1start:i1finish-1,:) = 0
      b_cyl_m%pf_re(n)%f2(1,i1start:i1finish-1,:) = 0
      b_cyl_m%pf_im(n)%f2(1,i1start:i1finish-1,:) = 0
      e_cyl_m%pf_re(n)%f2(2:3,i1start:i1finish,:) = 0
      e_cyl_m%pf_im(n)%f2(2:3,i1start:i1finish,:) = 0
      b_cyl_m%pf_re(n)%f2(2:3,i1start:i1finish,:) = 0
      b_cyl_m%pf_im(n)%f2(2:3,i1start:i1finish,:) = 0

    enddo

  endif

  ! Associate the pointers for conciseness
  if( this%per_type .ne. p_laguerre_gaussian .or. ell .eq. 0 ) then
    e_re_u => e_cyl_m%pf_re(1)
    e_im_u => e_cyl_m%pf_im(1)
    b_re_u => b_cyl_m%pf_re(1)
    b_im_u => b_cyl_m%pf_im(1)
    e_re_l => null()
    e_im_l => null()
    b_re_l => null()
    b_im_l => null()
  else
    e_re_l => e_cyl_m%pf_re(ml)
    b_re_l => b_cyl_m%pf_re(ml)
    e_re_u => e_cyl_m%pf_re(mu)
    b_re_u => b_cyl_m%pf_re(mu)
    if ( ml > 0 ) then
      e_im_l => e_cyl_m%pf_im(ml)
      b_im_l => b_cyl_m%pf_im(ml)
    else
      e_im_l => null()
      b_im_l => null()
    endif
    if ( mu > 0 ) then
      e_im_u => e_cyl_m%pf_im(mu)
      b_im_u => b_cyl_m%pf_im(mu)
    else
      e_im_u => null()
      b_im_u => null()
    endif
  endif

  do i1 = i1start, i1finish

  ! print *, i1start, i1finish, e_re%nx_(1), lbnd_in_domain, rbnd_in_domain, domain_enclosed

    z_sim = zmin + real( i1 - 1, p_double ) * ldx(1)
    z_sim_2 = z_sim + 0.5_p_double * ldx(1)

    ! convert to lab frame coords
    if( this % if_boost ) then
      call this % apply_lorentz_transformation( z_sim, t_sim, z, t )
      call this % apply_lorentz_transformation( z_sim_2, t_sim, z_2, t_2 )
    else
      z = z_sim
      z_2 = z_sim_2
      t = t_sim
      t_2 = t_sim
    endif

    lenv = amp * this % t_envelope( t, z )
    lenv_2 = amp * this % t_envelope( t_2, z_2 )

    fast_phase = this%omega0 * ( z - prop_sign * t ) + this%phase0
    fast_phase_2 = this%omega0 * ( z_2 - prop_sign * t_2 ) + this%phase0

    do i2 = e_re0%lbound(3), e_re0%ubound(3)
    ! do i2 = 2, e_re%nx_(2)

      r = rmin + real(i2 - 1, p_double ) * ldx(2)
      r_2 = r + 0.5_p_double * ldx(2)

      if( this%per_type .ne. p_laguerre_gaussian .or. ell .eq. 0 ) then

        select case( this%pol_type )

        ! lin pol
        case( 0 )

          ! e_r
          call this % per_envelope( z, r_2, per_env_amp, per_env_phase )
          total_phase = per_env_phase + fast_phase

          e_re_u%f2(2,i1,i2) = e_re_u%f2(2,i1,i2) + lenv * per_env_amp * cos( total_phase ) * cos( this%pol )
          e_im_u%f2(2,i1,i2) = e_im_u%f2(2,i1,i2) + lenv * per_env_amp * cos( total_phase ) * sin( this%pol )

          ! e_phi
          call this % per_envelope( z, r, per_env_amp, per_env_phase )
          total_phase = per_env_phase + fast_phase

          e_re_u%f2(3,i1,i2) = e_re_u%f2(3,i1,i2) + lenv * per_env_amp * cos( total_phase ) * sin( this%pol )
          e_im_u%f2(3,i1,i2) = e_im_u%f2(3,i1,i2) - lenv * per_env_amp * cos( total_phase ) * cos( this%pol )

          ! b_r
          call this % per_envelope( z_2, r, per_env_amp, per_env_phase )
          total_phase = per_env_phase + fast_phase_2

          b_re_u%f2(2,i1,i2) = b_re_u%f2(2,i1,i2) - lenv_2 * per_env_amp * cos( total_phase ) * sin( this%pol ) * prop_sign
          b_im_u%f2(2,i1,i2) = b_im_u%f2(2,i1,i2) + lenv_2 * per_env_amp * cos( total_phase ) * cos( this%pol ) * prop_sign

          ! b_phi
          call this % per_envelope( z_2, r_2, per_env_amp, per_env_phase )
          total_phase = per_env_phase + fast_phase_2

          b_re_u%f2(3,i1,i2) = b_re_u%f2(3,i1,i2) + lenv_2 * per_env_amp * cos( total_phase ) * cos( this%pol ) * prop_sign
          b_im_u%f2(3,i1,i2) = b_im_u%f2(3,i1,i2) + lenv_2 * per_env_amp * cos( total_phase ) * sin( this%pol ) * prop_sign

        ! handle circ pol
        ! same as general ell profile below
        case( 1 )

          ! e_r
          call this % per_envelope( z, r_2, per_env_amp, per_env_phase )
          total_phase = per_env_phase + fast_phase

          e_re_u%f2(2,i1,i2) = e_re_u%f2(2,i1,i2) + lenv * per_env_amp * cos( total_phase - this%pol )
          e_im_u%f2(2,i1,i2) = e_im_u%f2(2,i1,i2) - lenv * per_env_amp * sin( total_phase - this%pol )

          ! e_phi
          call this % per_envelope( z, r, per_env_amp, per_env_phase )
          total_phase = per_env_phase + fast_phase

          e_re_u%f2(3,i1,i2) = e_re_u%f2(3,i1,i2) - lenv * per_env_amp * sin( total_phase - this%pol )
          e_im_u%f2(3,i1,i2) = e_im_u%f2(3,i1,i2) - lenv * per_env_amp * cos( total_phase - this%pol )

          ! b_r
          call this % per_envelope( z_2, r, per_env_amp, per_env_phase )
          total_phase = per_env_phase + fast_phase_2

          b_re_u%f2(2,i1,i2) = b_re_u%f2(2,i1,i2) + lenv_2 * per_env_amp * sin( total_phase - this%pol ) * prop_sign
          b_im_u%f2(2,i1,i2) = b_im_u%f2(2,i1,i2) + lenv_2 * per_env_amp * cos( total_phase - this%pol ) * prop_sign

          ! b_phi
          call this % per_envelope( z_2, r_2, per_env_amp, per_env_phase )
          total_phase = per_env_phase + fast_phase_2

          b_re_u%f2(3,i1,i2) = b_re_u%f2(3,i1,i2) + lenv_2 * per_env_amp * cos( total_phase - this%pol ) * prop_sign
          b_im_u%f2(3,i1,i2) = b_im_u%f2(3,i1,i2) - lenv_2 * per_env_amp * sin( total_phase - this%pol ) * prop_sign

        ! same as general ell except with sign switches for imaginary parts
        case( -1 )

          ! e_r
          call this % per_envelope( z, r_2, per_env_amp, per_env_phase )
          total_phase = per_env_phase + fast_phase

          e_re_u%f2(2,i1,i2) = e_re_u%f2(2,i1,i2) + lenv * per_env_amp * cos( total_phase + this%pol )
          e_im_u%f2(2,i1,i2) = e_im_u%f2(2,i1,i2) + lenv * per_env_amp * sin( total_phase + this%pol )

          ! e_phi
          call this % per_envelope( z, r, per_env_amp, per_env_phase )
          total_phase = per_env_phase + fast_phase

          e_re_u%f2(3,i1,i2) = e_re_u%f2(3,i1,i2) + lenv * per_env_amp * sin( total_phase + this%pol )
          e_im_u%f2(3,i1,i2) = e_im_u%f2(3,i1,i2) - lenv * per_env_amp * cos( total_phase + this%pol )

          ! b_r
          call this % per_envelope( z_2, r, per_env_amp, per_env_phase )
          total_phase = per_env_phase + fast_phase_2

          b_re_u%f2(2,i1,i2) = b_re_u%f2(2,i1,i2) - lenv_2 * per_env_amp * sin( total_phase + this%pol ) * prop_sign
          b_im_u%f2(2,i1,i2) = b_im_u%f2(2,i1,i2) + lenv_2 * per_env_amp * cos( total_phase + this%pol ) * prop_sign

          ! b_phi
          call this % per_envelope( z_2, r_2, per_env_amp, per_env_phase )
          total_phase = per_env_phase + fast_phase_2

          b_re_u%f2(3,i1,i2) = b_re_u%f2(3,i1,i2) + lenv_2 * per_env_amp * cos( total_phase + this%pol ) * prop_sign
          b_im_u%f2(3,i1,i2) = b_im_u%f2(3,i1,i2) + lenv_2 * per_env_amp * sin( total_phase + this%pol ) * prop_sign

        end select


      ! else: is a laguerre pulse with ell /= 0
      else

        select case( this%pol_type )

        ! lin pol
        case( 0 )

          ! e_r
          call this % per_envelope_laguerre( z, r_2, per_env_amp, per_env_phase )
          total_phase = per_env_phase + fast_phase

          e_re_l%f2(2,i1,i2) = e_re_l%f2(2,i1,i2) + 0.5 * lenv * per_env_amp * cos( total_phase + this%pol )

          if( ml > 0 ) then
            e_im_l%f2(2,i1,i2) = e_im_l%f2(2,i1,i2) - 0.5 * lenv * per_env_amp * sin( total_phase + this%pol ) * sign_ell
          endif

          e_re_u%f2(2,i1,i2) = e_re_u%f2(2,i1,i2) + 0.5 * lenv * per_env_amp * cos( total_phase - this%pol )

          if( mu > 0 ) then
            e_im_u%f2(2,i1,i2) = e_im_u%f2(2,i1,i2) - 0.5 * lenv * per_env_amp * sin( total_phase - this%pol ) * sign_ell
          endif

          ! e_phi
          call this % per_envelope_laguerre( z, r, per_env_amp, per_env_phase )
          total_phase = per_env_phase + fast_phase

          e_re_l%f2(3,i1,i2) = e_re_l%f2(3,i1,i2) + 0.5 * lenv * per_env_amp * sin( total_phase + this%pol )

          if( ml > 0 ) then
            e_im_l%f2(3,i1,i2) = e_im_l%f2(3,i1,i2) + 0.5 * lenv * per_env_amp * cos( total_phase + this%pol ) * sign_ell
          endif

          e_re_u%f2(3,i1,i2) = e_re_u%f2(3,i1,i2) - 0.5 * lenv * per_env_amp * sin( total_phase - this%pol )

          if( mu > 0 ) then
            e_im_u%f2(3,i1,i2) = e_im_u%f2(3,i1,i2) - 0.5 * lenv * per_env_amp * cos( total_phase - this%pol ) * sign_ell
          endif

          ! b_r
          call this % per_envelope_laguerre( z_2, r, per_env_amp, per_env_phase )
          total_phase = per_env_phase + fast_phase_2

          b_re_l%f2(2,i1,i2) = b_re_l%f2(2,i1,i2) - 0.5 * lenv_2 * per_env_amp * sin( total_phase + this%pol ) * prop_sign

          if( ml > 0 ) then
            b_im_l%f2(2,i1,i2) = b_im_l%f2(2,i1,i2) - 0.5 * lenv_2 * per_env_amp * cos( total_phase + this%pol ) * prop_sign * sign_ell
          endif

          b_re_u%f2(2,i1,i2) = b_re_u%f2(2,i1,i2) + 0.5 * lenv_2 * per_env_amp * sin( total_phase - this%pol ) * prop_sign

          if( mu > 0 ) then
            b_im_u%f2(2,i1,i2) = b_im_u%f2(2,i1,i2) + 0.5 * lenv_2 * per_env_amp * cos( total_phase - this%pol ) * prop_sign * sign_ell
          endif

          ! b_phi
          call this % per_envelope_laguerre( z_2, r_2, per_env_amp, per_env_phase )
          total_phase = per_env_phase + fast_phase_2

          b_re_l%f2(3,i1,i2) = b_re_l%f2(3,i1,i2) + 0.5 * lenv_2 * per_env_amp * cos( total_phase + this%pol ) * prop_sign

          if( ml > 0 ) then
            b_im_l%f2(3,i1,i2) = b_im_l%f2(3,i1,i2) -0.5 * lenv_2 * per_env_amp * sin( total_phase + this%pol ) * prop_sign * sign_ell
          endif

          b_re_u%f2(3,i1,i2) = b_re_u%f2(3,i1,i2) + 0.5 * lenv_2 * per_env_amp * cos( total_phase - this%pol ) * prop_sign

          if( mu > 0 ) then
            b_im_u%f2(3,i1,i2) = b_im_u%f2(3,i1,i2) -0.5 * lenv_2 * per_env_amp * sin( total_phase - this%pol ) * prop_sign * sign_ell
          endif


        ! handle circ pol
        case( 1 )

          ! e_r
          call this % per_envelope_laguerre( z, r_2, per_env_amp, per_env_phase )
          total_phase = per_env_phase + fast_phase

          e_re_u%f2(2,i1,i2) = e_re_u%f2(2,i1,i2) + lenv * per_env_amp * cos( total_phase - this%pol )

          if( mu > 0 ) then
            e_im_u%f2(2,i1,i2) = e_im_u%f2(2,i1,i2) - lenv * per_env_amp * sin( total_phase - this%pol ) * sign_ell
          endif

          ! e_phi
          call this % per_envelope_laguerre( z, r, per_env_amp, per_env_phase )
          total_phase = per_env_phase + fast_phase

          e_re_u%f2(3,i1,i2) = e_re_u%f2(3,i1,i2) - lenv * per_env_amp * sin( total_phase - this%pol )

          if( mu > 0 ) then
            e_im_u%f2(3,i1,i2) = e_im_u%f2(3,i1,i2) - lenv * per_env_amp * cos( total_phase - this%pol ) * sign_ell
          endif

          ! b_r
          call this % per_envelope_laguerre( z_2, r, per_env_amp, per_env_phase )
          total_phase = per_env_phase + fast_phase_2

          b_re_u%f2(2,i1,i2) = b_re_u%f2(2,i1,i2) + lenv_2 * per_env_amp * sin( total_phase - this%pol ) * prop_sign

          if( mu > 0 ) then
            b_im_u%f2(2,i1,i2) = b_im_u%f2(2,i1,i2) + lenv_2 * per_env_amp * cos( total_phase - this%pol ) * prop_sign * sign_ell
          endif

          ! b_phi
          call this % per_envelope_laguerre( z_2, r_2, per_env_amp, per_env_phase )
          total_phase = per_env_phase + fast_phase_2

          b_re_u%f2(3,i1,i2) = b_re_u%f2(3,i1,i2) + lenv_2 * per_env_amp * cos( total_phase - this%pol ) * prop_sign

          if( mu > 0 ) then
            b_im_u%f2(3,i1,i2) = b_im_u%f2(3,i1,i2) - lenv_2 * per_env_amp * sin( total_phase - this%pol ) * prop_sign * sign_ell
          endif


        case( -1 )

          ! e_r
          call this % per_envelope_laguerre( z, r_2, per_env_amp, per_env_phase )
          total_phase = per_env_phase + fast_phase

          e_re_l%f2(2,i1,i2) = e_re_l%f2(2,i1,i2) + lenv * per_env_amp * cos( total_phase + this%pol )

          if( ml > 0 ) then
            e_im_l%f2(2,i1,i2) = e_im_l%f2(2,i1,i2) - lenv * per_env_amp * sin( total_phase + this%pol ) * sign_ell
          endif

          ! e_phi
          call this % per_envelope_laguerre( z, r, per_env_amp, per_env_phase )
          total_phase = per_env_phase + fast_phase

          e_re_l%f2(3,i1,i2) = e_re_l%f2(3,i1,i2) + lenv * per_env_amp * sin( total_phase + this%pol )

          if( ml > 0 ) then
            e_im_l%f2(3,i1,i2) = e_im_l%f2(3,i1,i2) + lenv * per_env_amp * cos( total_phase + this%pol ) * sign_ell
          endif

          ! b_r
          call this % per_envelope_laguerre( z_2, r, per_env_amp, per_env_phase )
          total_phase = per_env_phase + fast_phase_2

          b_re_l%f2(2,i1,i2) = b_re_l%f2(2,i1,i2) - lenv_2 * per_env_amp * sin( total_phase + this%pol ) * prop_sign

          if( ml > 0 ) then
            b_im_l%f2(2,i1,i2) = b_im_l%f2(2,i1,i2) - lenv_2 * per_env_amp * cos( total_phase + this%pol ) * prop_sign * sign_ell
          endif

          ! b_phi
          call this % per_envelope_laguerre( z_2, r_2, per_env_amp, per_env_phase )
          total_phase = per_env_phase + fast_phase_2

          b_re_l%f2(3,i1,i2) = b_re_l%f2(3,i1,i2) + lenv_2 * per_env_amp * cos( total_phase + this%pol ) * prop_sign

          if( ml > 0 ) then
            b_im_l%f2(3,i1,i2) = b_im_l%f2(3,i1,i2) - lenv_2 * per_env_amp * sin( total_phase + this%pol ) * prop_sign * sign_ell
          endif

        end select

      endif

    enddo ! i2

  enddo ! i1


  if( this%if_clean_fields ) then

    if( this%propagation .eq. p_forward ) then

      e_cyl_m%pf_re(0)%f2(:,:i1start-1,:) = 0.0_p_k_fld
      b_cyl_m%pf_re(0)%f2(:,:i1start-1,:) = 0.0_p_k_fld

      do n = 1, ubound( e_cyl_m%pf_re, 1 )

        e_cyl_m%pf_re(n)%f2(:,:i1start-1,:) = 0.0_p_k_fld
        e_cyl_m%pf_im(n)%f2(:,:i1start-1,:) = 0.0_p_k_fld
        b_cyl_m%pf_re(n)%f2(:,:i1start-1,:) = 0.0_p_k_fld
        b_cyl_m%pf_im(n)%f2(:,:i1start-1,:) = 0.0_p_k_fld

      enddo

    else

      e_cyl_m%pf_re(0)%f2(:,i1finish+1:,:) = 0.0_p_k_fld
      b_cyl_m%pf_re(0)%f2(:,i1finish+1:,:) = 0.0_p_k_fld

      do n = 1, ubound( e_cyl_m%pf_re, 1 )

        e_cyl_m%pf_re(n)%f2(:,i1finish+1:,:) = 0.0_p_k_fld
        e_cyl_m%pf_im(n)%f2(:,i1finish+1:,:) = 0.0_p_k_fld
        b_cyl_m%pf_re(n)%f2(:,i1finish+1:,:) = 0.0_p_k_fld
        b_cyl_m%pf_im(n)%f2(:,i1finish+1:,:) = 0.0_p_k_fld

      enddo

    endif

  endif

end subroutine set_fields_zpulse_mov_wall_cyl_modes
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
subroutine per_envelope_zpulse_mov_wall_cyl_modes( this, x1, r, amp, phase )

  class( t_zpulse_mov_wall_cyl_modes ), intent(inout) :: this
  real(p_double), intent(in) :: x1, r
  real(p_double), intent(out) :: amp, phase

  select case( this%per_type )

  case( p_hermite_gaussian )
    call this % per_envelope_gaussian( x1, r, amp, phase )

  case( p_laguerre_gaussian )
    call this % per_envelope_laguerre( x1, r, amp, phase )

  case default
    print *, 'ERROR in os-zpulse-cyl-modes: per_type is not implemented'
    stop

  end select

end subroutine per_envelope_zpulse_mov_wall_cyl_modes
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
subroutine per_envelope_gaussian_zpulse_mov_wall_cyl_modes( this, x1, r, amp, phase )

  class( t_zpulse_mov_wall_cyl_modes ), intent(inout) :: this
  real(p_double), intent(in) :: x1, r
  real(p_double), intent(out) :: amp, phase

  real(p_double) :: z, z0, rho2
  real(p_double) :: rWl2, rWu, curv, gouy_shift

  rho2 = r ** 2

  z = x1 - this%per_focus(1)

  z0 = this%omega0 * this%per_w0(1)**2 / 2

  ! In some situations (i.e. chirps) z0 may be 0 so z = 0 must be treated
  ! as a special case
  if ( z /= 0 ) then
    rWl2 = z0**2 / (z0**2 + z**2)
    curv = 0.5_p_double*rho2*z/(z**2 + z0**2)
    gouy_shift = atan2( z, z0 )
  else
    rWl2 = 1
    curv = 0
    gouy_shift = 0
  endif

  rWu = sqrt(2 * rWl2) / this%per_w0(1)

  ! for mode(1) = 0 and mode(2) = 0 this reproduces the gaussian beam
  amp = sqrt(rWl2) * exp( - rho2 * rWl2 /this%per_w0(1)**2 )

  phase = this%omega0 * curv - gouy_shift - this%omega0 * this%lon_center()


end subroutine per_envelope_gaussian_zpulse_mov_wall_cyl_modes
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
subroutine per_envelope_laguerre_zpulse_mov_wall_cyl_modes( this, x1, r, amp, phase )

  class( t_zpulse_mov_wall_cyl_modes ), intent(inout) :: this
  real(p_double), intent(in) :: x1, r
  real(p_double), intent(out) :: amp, phase

  real(p_double) :: z, z0, rho2
  real(p_double) :: rWl2, rWu, curv, gouy_shift

  rho2 = r**2

  z = x1 - this%per_focus(1)

  z0 = this%omega0 * this%per_w0(1)**2 / 2

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

  amp = sqrt(rWl2) * &
  ( sqrt(rho2)*sqrt(rWl2)/this%per_w0(1) )**abs(this%per_tem_mode(2)) * &
  laguerre( this%per_tem_mode(1), abs(this%per_tem_mode(2)), &
    2*rho2*rWl2/this%per_w0(1)**2 ) * &
  exp( - rho2 * rWl2/this%per_w0(1)**2 )

  phase = this%omega0 * curv - gouy_shift - this%omega0 * this%lon_center()

end subroutine per_envelope_laguerre_zpulse_mov_wall_cyl_modes
!-----------------------------------------------------------------------------------------



!-----------------------------------------------------------------------------------------
function t_env_zpulse_mov_wall_cyl_modes( this, t_sim, z_sim ) result(t_env)

  implicit none

  class( t_zpulse_mov_wall_cyl_modes ), intent(inout) :: this ! intent must be inout because of
                                                ! the eval( f_parser ) function
  real(p_double), intent(in) :: t_sim, z_sim

  real(p_double) :: t_env, t

  real(p_k_fparse), dimension(1) :: t_dbl

  ! ignore pulse broadening
  if( this%propagation .eq. p_forward ) then
    t = t_sim - z_sim + this % xi0
  else
    t = t_sim + z_sim - this % xi0
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

end function t_env_zpulse_mov_wall_cyl_modes
!-----------------------------------------------------------------------------------------


! TODO: Better support single precision in this function
!-----------------------------------------------------------------------------------------
subroutine apply_divcorr_zpulse_mov_wall_cyl_modes( this, e_cyl_m, b_cyl_m, g_space, &
                                                    gir_pos, nx_p_min, g_nx, t, dt, &
                                                    no_co, ell )

  implicit none

  integer, parameter :: rank = 2

  class( t_zpulse_mov_wall_cyl_modes ), intent(inout) :: this
  type( t_cyl_modes ), intent(inout) :: e_cyl_m, b_cyl_m
  type( t_space ), intent(in) :: g_space
  integer, intent(in) :: gir_pos
  integer, intent(in), dimension(:) :: nx_p_min, g_nx
  real(p_double), intent(in) :: t
  real(p_double), intent(in) :: dt
  class( t_node_conf ), intent(in) :: no_co
  integer, intent(in) :: ell


  type( t_vdf ), pointer :: b_re, b_im, e_re, e_im
  real(p_double), dimension(rank) :: ldx
  real(p_double), dimension(2,rank) :: g_x_range
  real(p_double) :: rmin, zmin, zmax

  real(p_double) :: gpos
  integer :: gipos, ipos

  real(p_double) :: e1_re, b1_re, e1_im, b1_im ! e1, b1 field correction
  real(p_double) :: e2p_re, e2m_re, b2p_re, b2m_re, e2p_im, e2m_im, b2p_im, b2m_im, e3_re, e3_im, b3_re, b3_im
  real(p_double) :: rp, rm, invr, invr2, rp2, rm2
  real(p_double) :: rcent, rcent2
  real(p_double) :: dx1, dx2

  integer :: i1, i2

  integer :: lbnd, rbnd
  logical :: lbnd_in_domain, rbnd_in_domain, domain_enclosed

  integer :: i1start, i1finish, direction, i1_initial_data
  integer :: recv_neighbor, send_neighbor
  logical :: start_in_domain, finish_in_domain

  integer :: size, nx1, nx2
  integer :: ierr
  integer, dimension( MPI_STATUS_SIZE ) :: status
  real(p_double), dimension(:,:), allocatable :: divcorr_buf
  real(p_double), dimension(2) :: divcorr_e2
  integer :: gshift_i2, i2_0

  e_re => e_cyl_m%pf_re(ell)
  b_re => b_cyl_m%pf_re(ell)

  if( ell > 0 ) then
    e_im => e_cyl_m%pf_im(ell)
    b_im => b_cyl_m%pf_im(ell)
  else
    ! To suppress compiler warnings
    e_im => null()
    b_im => null()
  endif

  ldx(1) = b_re%dx(1)
  ldx(2) = b_re%dx(2)

  nx1 = e_re%nx_(1)
  nx2 = e_re%nx_(2)

  call get_x_bnd( g_space, g_x_range )
  zmin = real( g_x_range(p_lower,1), p_double )
  zmax = real( g_x_range(p_upper,1), p_double )
  rmin = real( g_x_range(p_lower,2), p_double )

  select case( this%interpolation )
  case( p_linear, p_cubic )
    ! Nothing to change here
  case( p_quadratic, p_quartic )
    zmin = zmin + ldx(1) * 0.5_p_double
    zmax = zmax + ldx(1) * 0.5_p_double
    ! rmin = rmin + ldx(2) * 0.5_p_double
  end select

  ! start near domain boundary
  if( this%if_launch_from_wall ) then

    if( this%propagation == p_forward ) then
      gipos = 1
    else
      gipos = floor( ( zmax - zmin ) / ldx(1) ) - 1
    endif

  ! else compute the wall position in global coordinates, stored as gpos
  else
    gpos = this%wall_pos + this%wall_vel*t
    gipos = int((gpos - zmin)/ldx(1)) + 1
  endif

  ipos = gipos - nx_p_min(1) + 1

  ! switch to local coordinates
  zmin = zmin + real(nx_p_min(1)-1, p_double)*ldx(1)
  rmin = rmin + real(nx_p_min(2)-1, p_double)*ldx(2)

  gshift_i2 = gir_pos - 2

  if ( gshift_i2 < 0 ) then
    i2_0 = 2
  else
    i2_0 = 1
  endif

  ! left and right bounds of the calculation
  if( this%propagation == p_forward ) then
    lbnd = ipos
    rbnd = ipos + this%ncells
  else
    lbnd = ipos - this%ncells - 1
    rbnd = ipos
  endif

  lbnd_in_domain = ( lbnd >= 1) .and. (lbnd <= nx1)
  rbnd_in_domain = (rbnd >= 1) .and. (rbnd <= nx1)
  domain_enclosed = (lbnd < 1) .and. (rbnd > nx1)

  ! wall outside is completely outside this domain: do nothing
  if( .not.( lbnd_in_domain .or. rbnd_in_domain .or. domain_enclosed ) ) return

  if( this%propagation == p_forward ) then
    direction = -1
    recv_neighbor = p_upper
    send_neighbor = p_lower
    i1start = min( rbnd - 1, nx1 )
    i1finish = max( lbnd, 1 )
    i1_initial_data = rbnd
    start_in_domain = rbnd_in_domain
    finish_in_domain = lbnd_in_domain
  else
    direction = 1
    recv_neighbor = p_lower
    send_neighbor = p_upper
    i1start = max( 1, lbnd + 1 )
    i1finish = min( rbnd, nx1 )
    i1_initial_data = lbnd
    start_in_domain = lbnd_in_domain
    finish_in_domain = rbnd_in_domain
  endif

  ! useful debug data
  ! print *, "debug: ", i1start, i1finish, nx1, lbnd_in_domain, rbnd_in_domain, domain_enclosed

  ! If there is a lower neighbor, then e2 needs to be sent and received for the cell in
  ! front of the wall. This value got changed by the last field solve but has not yet
  ! been communicated via update_boundary.
  if ( start_in_domain ) then
    if ( no_co%neighbor( p_lower, 2 ) /= -1 ) then
      call mpi_recv( divcorr_e2, 2, MPI_DOUBLE, no_co%neighbor( p_lower, 2 ) - 1, 0, &
                     no_co%comm, status, ierr )
      e_re%f2(2,i1_initial_data,0) = divcorr_e2(1)
      if (ell > 0) e_im%f2(2,i1_initial_data,0) = divcorr_e2(2)
    endif
    if ( no_co%neighbor( p_upper, 2 ) /= -1 ) then
      divcorr_e2 = 0.0_p_double
      divcorr_e2(1) = e_re%f2(2,i1_initial_data,nx2)
      if (ell > 0) divcorr_e2(2) = e_im%f2(2,i1_initial_data,nx2)
      call mpi_send( divcorr_e2, 2, MPI_DOUBLE, no_co%neighbor( p_upper, 2 ) - 1, 0, &
                     no_co%comm, ierr )
    endif
  endif

  size = nx2 + 1 ! include the first cell below the axis as well
  allocate( divcorr_buf( 8, 0:nx2 ) )

  ! now apply divergence correction
  dx1 = ldx(1)
  dx2 = ldx(2)

  ! check if integral start is in domain
  if( .not. start_in_domain ) then

    call mpi_recv( divcorr_buf, 8 * size, MPI_DOUBLE, &
                  no_co%neighbor( recv_neighbor, 1 ) - 1, &
                  0, no_co%comm, status, ierr )

    call check_error(ierr,"zpulse wall: recv fail",p_err_mpi,"cyl_modes/os-zpulse-mov-wall-cyl-modes.f03",1108)

  endif

  do i2 = i2_0, e_re%nx_(2)

    ! this is for the e-fields
    rp = rmin + (real(i2,p_double) - 0.5_p_double)*dx2 ! e2
    rcent = rmin + (real(i2,p_double) - 1.0_p_double)*dx2 ! e3, e1
    rm = rmin + (real(i2,p_double) - 1.5_p_double)*dx2
    invr = 1.0_p_double/rcent ! e3, e1

    ! this is for the b-fields
    rp2 = rmin + (real(i2,p_double) )*dx2
    rcent2 = rmin + (real(i2,p_double) - 0.5_p_double)*dx2 ! b1, b3
    rm2 = rmin + (real(i2,p_double) - 1.0_p_double)*dx2 ! b2
    invr2 = 1.0_p_double/rcent2 ! b1, b3

    ! starting value: either use value from previous node if the
    ! divergence integral is not starting locally
    if( start_in_domain ) then

      e1_re = e_re%f2(1,i1_initial_data,i2)
      b1_re = b_re%f2(1,i1_initial_data,i2)

      if( ell .gt. 0 ) then
        e1_im = e_im%f2(1,i1_initial_data,i2)
        b1_im = b_im%f2(1,i1_initial_data,i2)
      endif

    else

      e1_re = divcorr_buf(1,i2)
      b1_re = divcorr_buf(3,i2)

      if( ell .gt. 0 ) then
        e1_im = divcorr_buf(2,i2)
        b1_im = divcorr_buf(4,i2)
      endif

      if ( this%propagation == p_forward ) then
        ! Get field values from the above node as well
        e_re%f2(2,nx1+1,i2) = divcorr_buf(5, i2)
        e_re%f2(2,nx1+1,i2-1) = divcorr_buf(5, i2-1)
        e_re%f2(3,nx1+1,i2) = divcorr_buf(7, i2)

        if( ell .gt. 0 ) then
          e_im%f2(2,nx1+1,i2) = divcorr_buf(6, i2)
          e_im%f2(2,nx1+1,i2-1) = divcorr_buf(6, i2-1)
          e_im%f2(3,nx1+1,i2) = divcorr_buf(8, i2)
        endif
      endif

    endif

    do i1 = i1start, i1finish, direction

      if( ell .gt. 0 ) then

        ! used in real part calculations
        b3_im = b_im%f2(3,i1,i2)
        e3_im = e_im%f2(3,i1+1,i2)

        ! used for imag calculations
        e3_re = e_re%f2(3,i1+1,i2)
        b3_re = b_re%f2(3,i1,i2)

        ! e3 im
        e2p_im = e_im%f2(2,i1+1,i2)
        e2m_im = e_im%f2(2,i1+1,i2-1)
        e1_im = e1_im - direction * invr*(dx1/dx2)*(rp*e2p_im - rm*e2m_im) + direction * ell*invr*e3_re*dx1
        e_im%f2(1, i1, i2) = real( e1_im, p_k_fld )

        ! b3 im
        b2p_im = b_im%f2(2,i1,i2+1)
        b2m_im = b_im%f2(2,i1,i2)
        b1_im = b1_im - direction * invr2*(dx1/dx2)*(rp2*b2p_im - rm2*b2m_im) + direction * ell*invr2*b3_re*dx1
        b_im%f2(1, i1, i2) = real( b1_im, p_k_fld )

        ! e3 real
        e2p_re = e_re%f2(2,i1+1,i2)
        e2m_re = e_re%f2(2,i1+1,i2-1)
        e1_re = e1_re - direction * invr*(dx1/dx2)*(rp*e2p_re - rm*e2m_re) - direction * ell*invr*e3_im*dx1
        e_re%f2(1, i1, i2) = real( e1_re, p_k_fld )

        ! b3 re
        b2p_re = b_re%f2(2,i1,i2+1)
        b2m_re = b_re%f2(2,i1,i2)
        b1_re = b1_re - direction * invr2*(dx1/dx2)*(rp2*b2p_re - rm2*b2m_re) - direction * ell*invr2*b3_im*dx1
        b_re%f2(1, i1, i2) = real( b1_re, p_k_fld )

      else
        ! imaginary parts aren't set for mode 0 

        ! e3_im = 0
        ! b3_im = 0

        e3_re = e_re%f2(3,i1+1,i2)
        b3_re = b_re%f2(3,i1,i2)

        ! e3 real
        e2p_re = e_re%f2(2,i1+1,i2)
        e2m_re = e_re%f2(2,i1+1,i2-1)
        e1_re = e1_re - direction * invr*(dx1/dx2)*(rp*e2p_re - rm*e2m_re)
        e_re%f2(1, i1, i2) = real( e1_re, p_k_fld )

        ! b3 im
        b2p_re = b_re%f2(2,i1,i2+1)
        b2m_re = b_re%f2(2,i1,i2)
        b1_re = b1_re - direction * invr2*(dx1/dx2)*(rp2*b2p_re - rm2*b2m_re)
        b_re%f2(1, i1, i2) = real( b1_re, p_k_fld )

      endif

    enddo ! i1

  enddo ! i2


  ! if integral finish is not in domain, we need to send the results of the integral
  ! to be finished in the next domain
  if( .not. finish_in_domain ) then

    divcorr_buf = 0.0_p_double

    divcorr_buf(1,1:nx2) = e_re%f2(1, i1finish, 1:nx2)
    divcorr_buf(3,1:nx2) = b_re%f2(1, i1finish, 1:nx2)
    divcorr_buf(5,0:nx2) = e_re%f2(2, i1finish, 0:nx2)
    divcorr_buf(7,1:nx2) = e_re%f2(3, i1finish, 1:nx2)

    if( ell .gt. 0 ) then
      divcorr_buf(2,1:nx2) = e_im%f2(1, i1finish, 1:nx2)
      divcorr_buf(4,1:nx2) = b_im%f2(1, i1finish, 1:nx2)
      divcorr_buf(6,0:nx2) = e_im%f2(2, i1finish, 0:nx2)
      divcorr_buf(8,1:nx2) = e_im%f2(3, i1finish, 1:nx2)
    endif

    call mpi_send( divcorr_buf, 8 * size, MPI_DOUBLE, &
                  no_co%neighbor( send_neighbor, 1 ) - 1, &
                  0, no_co%comm, ierr )

    call check_error(ierr,"zpulse wall: send fail",p_err_mpi,"cyl_modes/os-zpulse-mov-wall-cyl-modes.f03",1249)

  endif

  deallocate( divcorr_buf )

end subroutine apply_divcorr_zpulse_mov_wall_cyl_modes
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
subroutine correct_axis_fields_zpulse_mov_wall_cyl_modes( this, e_cyl_m, b_cyl_m, &
                                                          g_space, gir_pos, nx_p_min, &
                                                          g_nx, t, dt, no_co, ell )

  implicit none

  integer, parameter :: rank = 2

  class( t_zpulse_mov_wall_cyl_modes ), intent(inout) :: this
  type( t_cyl_modes ), intent(inout) :: e_cyl_m, b_cyl_m
  type( t_space ), intent(in) :: g_space
  integer, intent(in) :: gir_pos
  integer, intent(in), dimension(:) :: nx_p_min, g_nx
  real(p_double), intent(in) :: t
  real(p_double), intent(in) :: dt
  class( t_node_conf ), intent(in) :: no_co
  integer, intent(in) :: ell


  type(t_vdf), pointer :: b_re, b_im, e_re, e_im
  real(p_double) :: zmin
  real(p_double), dimension(rank) :: ldx
  real(p_double), dimension(2,rank) :: g_x_range
  real(p_k_fld) :: dtdr, f_t, o_t, n_e, o_e, phase_factor
  ! real(p_double) :: rmin, zmin, zmax

  real(p_double) :: gpos, zmax
  integer :: gipos, ipos
  integer :: lbnd, rbnd
  logical :: lbnd_in_domain, rbnd_in_domain, domain_enclosed

  integer :: i1, i1start, i1finish, i2, mm

  integer :: gshift_i2, i2_0

  e_re => e_cyl_m%pf_re(ell)
  b_re => b_cyl_m%pf_re(ell)

  if( ell > 0 ) then
    e_im => e_cyl_m%pf_im(ell)
    b_im => b_cyl_m%pf_im(ell)
  else
    ! To suppress compiler warnings
    e_im => null()
    b_im => null()
  endif

  gshift_i2 = gir_pos - 2

  if ( gshift_i2 < 0 ) then
    i2_0 = 2
  else
    i2_0 = 1
  endif

  call get_x_bnd( g_space, g_x_range )
  zmin = real( g_x_range(p_lower,1), p_double )
  zmax = real( g_x_range(p_upper,1), p_double )
  ldx(1) = b_re%dx(1)
  ldx(2) = b_re%dx(2)
  dtdr = real(dt/ldx(2), p_k_fld)
  f_t = real(4.0_p_k_fld/3.0_p_k_fld, p_k_fld)
  o_t = real(1.0_p_k_fld/3.0_p_k_fld, p_k_fld)
  n_e = real(9.0_p_k_fld/8.0_p_k_fld, p_k_fld)
  o_e = real(1.0_p_k_fld/8.0_p_k_fld, p_k_fld)

  if ( mod(ell,2) == 0 ) then
    phase_factor = 1.0_p_k_fld
  else
    phase_factor = -1.0_p_k_fld
  endif

  select case( this%interpolation )
  case( p_linear, p_cubic )
    ! Nothing to change here
  case( p_quadratic, p_quartic )
    zmin = zmin + ldx(1) * 0.5_p_double
    zmax = zmax + ldx(1) * 0.5_p_double
  end select

  ! start near domain boundary
  if( this%if_launch_from_wall ) then

    if( this%propagation == p_forward ) then
      gipos = 1
    else
      gipos = floor( ( zmax - zmin ) / ldx(1) ) - 1
    endif

  ! else compute the wall position in global coordinates, stored as gpos
  else
    gpos = this%wall_pos + this%wall_vel*t
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

  lbnd_in_domain = (lbnd >= 1) .and. (lbnd <= e_re%nx_(1))
  rbnd_in_domain = (rbnd >= 1) .and. (rbnd <= e_re%nx_(1))
  domain_enclosed = (lbnd < 1) .and. (rbnd > e_re%nx_(1))

  lbnd = max( lbnd, 1 )
  rbnd = min( rbnd, e_re%nx_(1) )

  i1start = lbnd
  i1finish = rbnd

  if( .not.( lbnd_in_domain .or. rbnd_in_domain .or. domain_enclosed ) ) return

  ! See comments in os-emf-solver-cyl-modes-yee.f03 for more explanation
  if( gshift_i2 < 0 ) then
     if (ell == 0) then
        do i1 = i1start, i1finish

          ! E1 can be nonzero on the axis
          e_re%f2( 1, i1, 1 ) = e_re%f2( 1, i1, 2 )
          e_re%f2( 1, i1, 0 ) = e_re%f2( 1, i1, 3 )

          ! E2 must be zero on the axis
          e_re%f2( 2, i1, 1 ) = 0.0_p_k_fld
          e_re%f2( 2, i1, 0 ) = - e_re%f2( 2, i1, 2 )

          ! E3 must be zero on the axis
          e_re%f2( 3, i1, 1 ) = - e_re%f2( 3, i1, 2 )
          e_re%f2( 3, i1, 0 ) = - e_re%f2( 3, i1, 3 )

          ! B1 can be nonzero on the axis
          b_re%f2( 1, i1, 1 ) = b_re%f2( 1, i1, 1 ) - 4 * dtdr * e_re%f2( 3, i1, 2 )
          b_re%f2( 1, i1, 0 ) = b_re%f2( 1, i1, 2 )

          ! B2 must be zero on the axis
          b_re%f2( 2, i1, 1 ) = - b_re%f2( 2, i1, 2 )
          b_re%f2( 2, i1, 0 ) = - b_re%f2( 2, i1, 3 )

          ! B3 must be zero on the axis
          b_re%f2( 3, i1, 1 ) = 0.0_p_k_fld
          b_re%f2( 3, i1, 0 ) = - b_re%f2( 3, i1, 2 )
        enddo
      elseif (ell == 1) then
        do i1 = i1start, i1finish
          ! E1 must be zero on the axis
          e_re%f2( 1, i1, 1 ) = - e_re%f2( 1, i1, 2 )
          e_re%f2( 1, i1, 0 ) = - e_re%f2( 1, i1, 3 )
          e_im%f2( 1, i1, 1 ) = - e_im%f2( 1, i1, 2 )
          e_im%f2( 1, i1, 0 ) = - e_im%f2( 1, i1, 3 )

          ! E2 can be nonzero on the axis
          e_re%f2( 2, i1, 1 ) = f_t * e_re%f2( 2, i1, 2 ) - o_t * e_re%f2( 2, i1, 3 )
          e_re%f2( 2, i1, 0 ) = e_re%f2( 2, i1, 2 )
          e_im%f2( 2, i1, 1 ) = f_t * e_im%f2( 2, i1, 2 ) - o_t * e_im%f2( 2, i1, 3 )
          e_im%f2( 2, i1, 0 ) = e_im%f2( 2, i1, 2 )

          ! E3 can be nonzero on the axis
          e_re%f2( 3, i1, 1 ) = 2.0_p_k_fld * e_im%f2( 2, i1, 1 ) - e_re%f2( 3, i1, 2 )
          e_re%f2( 3, i1, 0 ) = 2.0_p_k_fld * e_im%f2( 2, i1, 1 ) - e_re%f2( 3, i1, 3 )
          e_im%f2( 3, i1, 1 ) =-2.0_p_k_fld * e_re%f2( 2, i1, 1 ) - e_im%f2( 3, i1, 2 )
          e_im%f2( 3, i1, 0 ) =-2.0_p_k_fld * e_re%f2( 2, i1, 1 ) - e_im%f2( 3, i1, 3 )

          ! B1 must be zero on the axis
          b_re%f2( 1, i1, 1 ) = 0.0_p_k_fld
          b_re%f2( 1, i1, 0 ) = - b_re%f2( 1, i1, 2 )
          b_im%f2( 1, i1, 1 ) = 0.0_p_k_fld
          b_im%f2( 1, i1, 0 ) = - b_im%f2( 1, i1, 2 )

          ! B2 can be nonzero on the axis
          b_re%f2( 2, i1, 1 ) = b_re%f2( 2, i1, 2 )
          b_re%f2( 2, i1, 0 ) = b_re%f2( 2, i1, 3 )
          b_im%f2( 2, i1, 1 ) = b_im%f2( 2, i1, 2 )
          b_im%f2( 2, i1, 0 ) = b_im%f2( 2, i1, 3 )

          ! B3 can be nonzero on the axis
          b_re%f2( 3, i1, 1 ) = n_e * b_im%f2( 2, i1, 2 ) - o_e * b_im%f2( 2, i1, 3 )
          b_re%f2( 3, i1, 0 ) = b_re%f2( 3, i1, 2 )
          b_im%f2( 3, i1, 1 ) =-n_e * b_re%f2( 2, i1, 2 ) + o_e * b_re%f2( 2, i1, 3 )
          b_im%f2( 3, i1, 0 ) = b_im%f2( 3, i1, 2 )
        enddo
      else ! ell >= 2

        mm = ell - 1

        do i1 = i1start, i1finish
          do i2 = 1, ell - 2
            ! Zero out first few derivatives for E1
            e_re%f2( 1, i1, i2+1 ) = real((2*i2-1)**mm, p_k_fld) / real((2*mm-1)**mm, p_k_fld) * e_re%f2( 1, i1, mm+1 )
            e_im%f2( 1, i1, i2+1 ) = real((2*i2-1)**mm, p_k_fld) / real((2*mm-1)**mm, p_k_fld) * e_im%f2( 1, i1, mm+1 )
          enddo

          ! E1 must be zero on the axis
          e_re%f2( 1, i1, 1 ) = phase_factor * e_re%f2( 1, i1, 2 )
          e_re%f2( 1, i1, 0 ) = phase_factor * e_re%f2( 1, i1, 3 )
          e_im%f2( 1, i1, 1 ) = phase_factor * e_im%f2( 1, i1, 2 )
          e_im%f2( 1, i1, 0 ) = phase_factor * e_im%f2( 1, i1, 3 )

          ! E2 must be zero on the axis
          e_re%f2( 2, i1, 1 ) = 0.0_p_k_fld
          e_re%f2( 2, i1, 0 ) = -phase_factor * e_re%f2( 2, i1, 2 )
          e_im%f2( 2, i1, 1 ) = 0.0_p_k_fld
          e_im%f2( 2, i1, 0 ) = -phase_factor * e_im%f2( 2, i1, 2 )

          ! E3 must be zero on the axis
          e_re%f2( 3, i1, 1 ) = -phase_factor * e_re%f2( 3, i1, 2 )
          e_re%f2( 3, i1, 0 ) = -phase_factor * e_re%f2( 3, i1, 3 )
          e_im%f2( 3, i1, 1 ) = -phase_factor * e_im%f2( 3, i1, 2 )
          e_im%f2( 3, i1, 0 ) = -phase_factor * e_im%f2( 3, i1, 3 )

          do i2 = 1, ell - 2
            ! Zero out first few derivatives for B1
            b_re%f2( 1, i1, i2+1 ) = real(i2**mm, p_k_fld) / real(mm**mm, p_k_fld) * b_re%f2( 1, i1, mm+1 )
            b_im%f2( 1, i1, i2+1 ) = real(i2**mm, p_k_fld) / real(mm**mm, p_k_fld) * b_im%f2( 1, i1, mm+1 )
          enddo

          ! B1 must be zero on the axis
          b_re%f2( 1, i1, 1 ) = 0.0_p_k_fld
          b_re%f2( 1, i1, 0 ) = phase_factor * b_re%f2( 1, i1, 2 )
          b_im%f2( 1, i1, 1 ) = 0.0_p_k_fld
          b_im%f2( 1, i1, 0 ) = phase_factor * b_im%f2( 1, i1, 2 )

          ! B2 must be zero on the axis
          b_re%f2( 2, i1, 1 ) = -phase_factor * b_re%f2( 2, i1, 2 )
          b_re%f2( 2, i1, 0 ) = -phase_factor * b_re%f2( 2, i1, 3 )
          b_im%f2( 2, i1, 1 ) = -phase_factor * b_im%f2( 2, i1, 2 )
          b_im%f2( 2, i1, 0 ) = -phase_factor * b_im%f2( 2, i1, 3 )

          ! B3 must be zero on the axis
          b_re%f2( 3, i1, 1 ) = 0.0_p_k_fld
          b_re%f2( 3, i1, 0 ) = -phase_factor * b_re%f2( 3, i1, 2 )
          b_im%f2( 3, i1, 1 ) = 0.0_p_k_fld
          b_im%f2( 3, i1, 0 ) = -phase_factor * b_im%f2( 3, i1, 2 )
        enddo
      endif
    endif

end subroutine correct_axis_fields_zpulse_mov_wall_cyl_modes
!-----------------------------------------------------------------------------------------



end module m_zpulse_mov_wall_cyl_modes
