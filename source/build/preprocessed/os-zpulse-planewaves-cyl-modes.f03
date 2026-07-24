# 1 "cyl_modes/os-zpulse-planewaves-cyl-modes.f03"
# 1 "<built-in>" 1
# 1 "<built-in>" 3
# 467 "<built-in>" 3
# 1 "<command line>" 1
# 1 "<built-in>" 2
# 1 "cyl_modes/os-zpulse-planewaves-cyl-modes.f03" 2
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
# 2 "cyl_modes/os-zpulse-planewaves-cyl-modes.f03" 2
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
# 3 "cyl_modes/os-zpulse-planewaves-cyl-modes.f03" 2

module m_zpulse_planewaves_cyl_modes

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
# 7 "cyl_modes/os-zpulse-planewaves-cyl-modes.f03" 2

use m_parameters
use m_math
use m_fparser, only: p_max_expr_len, p_k_fparse, eval
use m_node_conf
use m_zpulse_planewaves, only : t_zpulse_planewaves
use m_input_file, only: t_input_file, get_namelist, p_max_nml_section_name
use m_emf_define, only: t_emf_bound, t_emf
use m_emf_cyl_modes, only: t_emf_cyl_modes
use m_space, only: t_space, x_bnd, get_x_bnd, if_move
use m_vdf_define, only: t_vdf
use m_grid_define, only: t_grid
use m_grid_cyl_modes,only: t_grid_cyl_modes
use m_space, only: xmin, xmax
use m_vdf_math, only: add

use m_zpulse_planewaves, only: type_time, type_space, type_vfocus_const, type_vfocus_accel, type_custom_full, type_custom_lon, type_custom_per
use m_zpulse_planewaves, only: k_lon_plane, k_lon_gaussian, k_lon_sin2, k_lon_polynomial
use m_zpulse_planewaves, only: k_per_plane, k_per_gaussian, k_per_laguerre_gaussian, k_per_hermite_gaussian
use m_zpulse_planewaves, only: n_sig_gaussian
use m_zpulse_planewaves, only: im_unit

! testing, debug and control
logical, parameter :: control = .true.

!-----------------------------------------------------------------------------------------
type, extends( t_zpulse_planewaves ) :: t_zpulse_planewaves_cyl_modes

 ! new variables go here
 integer :: m_pulse
 integer, dimension (:), allocatable :: n_array
 real(p_double), dimension (:,:,:), allocatable :: bessel_cache_left

contains
 procedure :: init => init_zpulse_planewaves_cyl_modes
 procedure :: read_input => read_input_zpulse_planewaves_cyl_modes
 procedure :: launch => launch_zpulse_planewaves_cyl_modes
 procedure :: check_dimensionality => check_dimensionality_zpulse_planewaves_cyl_modes

 procedure :: k_function => k_function_zpulse_planewaves_cyl_modes
 procedure :: k_per_function => k_per_function_zpulse_planewaves_cyl_modes

end type t_zpulse_planewaves_cyl_modes
!-----------------------------------------------------------------------------------------

contains

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
! Initialize zpulse_planewaves
subroutine init_zpulse_planewaves_cyl_modes( this, restart, t , dt, b, g_space, nx_p_min, no_co, interpolation )

 implicit none

 class( t_zpulse_planewaves_cyl_modes ), intent(inout) :: this
 logical, intent(in) :: restart
 real( p_double ), intent(in) :: t, dt
 type( t_vdf ) , intent(in) :: b
 type( t_space ), intent(in) :: g_space
 integer, intent(in), dimension(:) :: nx_p_min
 class( t_node_conf ), intent(in) :: no_co
 integer, intent(in) :: interpolation

 integer :: i_n, n_wave_now, num_waves
 integer :: ikz, ikr, n, in, s
 real(p_double) :: kz, kr, omega, max_amp
 complex(p_double) :: f_k
 real(p_double), dimension(2) :: ldx

 ! executable statements
 ldx(1) = b%dx_(1)
 ldx(2) = b%dx_(2)

 ! Need to preserve if_launch value for this zpulse
 this%was_launched = .false.

 if ( restart ) then
  this%was_launched = .true.
 endif

 if ( .not. this%if_launch ) then
  return
 endif

 if ( this%lon_type == k_lon_plane .or. this%per_type == k_per_plane ) then

  allocate( this%kvec_arr( 1, 3 ) )
  allocate( this%evec_arr( 1, 2 ) )

  kz = this%omega0
  kr = this%per_w0(1)

  f_k = this%a0 * s**n * kr! * sqrt( kr**2 + kz**2 )
  !omega = 2/dt * asin( dt * sqrt( sin(kz*ldx(1)/2)**2/ldx(1)**2 + sin(kr*ldx(2)/2)**2/ldx(2)**2 ) )
  omega = sqrt( kz**2 + kr**2 )

  ! multiply components of e and b by f_k
  this%kvec_arr( 1,: ) = (/ this%per_w0(2), kz, kr /)
  this%evec_arr( 1,1 ) = f_k
  this%evec_arr( 1,2 ) = omega
  this%m_pulse = int( this%per_w0(2) ) + 1

  if (mpi_node()==0) print *,this%kvec_arr
  if (mpi_node()==0) print *,this%evec_arr
  if (mpi_node()==0) print *,this%m_pulse

  return
 endif

 ! get minimum and maximum m
 ! ___THIS IS NO LONGER USED AND CAN BE REMOVED ONCE ALL WALL INJECTIONS ARE UPDATED
 select case ( this%pulse_type )
  case ( type_custom_full )
   allocate ( this%n_array( int(this%per2_int_range(2)-this%per2_int_range(1)+1) ) )
   do i_n = 1, size(this%n_array)
    this%n_array(i_n) = int(this%per2_int_range(1)) - 1 + i_n
   enddo
  case default
   select case ( this%per_type )
    case ( k_per_gaussian )
     allocate ( this%n_array(1) )
     this%n_array(1) = 0
    case ( k_per_laguerre_gaussian )
     allocate ( this%n_array(1) )
     this%n_array(1) = this%per_tem_mode(2)
   end select
 end select

 ! calculate basis waves
 n_wave_now = 1

 num_waves = int( ( this%lon_int_range(2) - this%lon_int_range(1) ) / this%dk_waves(1) + 10 )
 num_waves = num_waves * int( ( this%per1_int_range(2) - this%per1_int_range(1) ) / this%dk_waves(2) + 10 )
 num_waves = num_waves * size( this%n_array )

 allocate( this%kvec_arr( num_waves, 3 ) )
 allocate( this%evec_arr( num_waves, 2 ) )

 ! iterate for every plane wave
 do in = 1, size(this%n_array)

  ! get m mode with respective sign
  n = abs( this%n_array(in) )
  s = sign( 1, this%n_array(in) )

  ! calculate maximum amplitude for tolerance
  ! assuming the largest contribution is the principal wave ...
  max_amp = this%tolerance * abs( this%a0 * this%k_function( this%omega0, 0.0_p_double, real(n, p_double) ) * this%dk_waves(1) * this%dk_waves(2) )
  if ( max_amp == 0 ) max_amp = this%tolerance

 ! iterate for every plane wave
 kr = this%per1_int_range(1) - mod( this%per1_int_range(1), this%dk_waves(2) )
 do while ( kr < this%per1_int_range(2) )
  ! calculate direction of each wave - perpendicular
  kr = kr + this%dk_waves(2)
 kz = this%lon_int_range(1) - mod( this%lon_int_range(1), this%dk_waves(1) )
 do while ( kz < this%lon_int_range(2) )
  ! calculate direction of each wave - longitudinal
  kz = kz + this%dk_waves(1)
 ! --------------------------------------------------

  ! calculate intensity of each wave
  f_k = this%a0 * this%k_function( abs(kz), kr, real(n, p_double) ) * this%dk_waves(1) * this%dk_waves(2)

  ! tolerance for speep up
  if ( abs(f_k) < max_amp ) cycle

  f_k = f_k * s**n * kr !* sqrt( kr**2 + kz**2 )
  omega = 2/dt * asin( dt * sqrt( sin(kz*ldx(1)/2)**2/ldx(1)**2 + sin(kr*ldx(2)/2)**2/ldx(2)**2 ) )
  !omega = sqrt( kz**2 + kr**2 )

  ! multiply components of e and b by f_k
  this%kvec_arr( n_wave_now,: ) = (/ this%n_array(in)*1.0_p_double, kz, kr /)
  this%evec_arr( n_wave_now,1 ) = f_k
  this%evec_arr( n_wave_now,2 ) = omega
  n_wave_now = n_wave_now + 1

 enddo
 enddo
 enddo

 this%kvec_arr => this%kvec_arr( :n_wave_now-1,: )
 this%evec_arr => this%evec_arr( :n_wave_now-1,: )

 if (mpi_node()==0) print *,"   - will be injecting", n_wave_now-1, "waves"
 if (mpi_node()==0) print *,"   - longitudinal periodicity", 2*pi/this%dk_waves(1)
 if (mpi_node()==0) print *,"   - radial periodicity (approx)", 2*pi/this%dk_waves(2)
 if (mpi_node()==0) print *,"   - lon_int_range", this%lon_int_range
 if (mpi_node()==0) print *,"   - per1_int_range", this%per1_int_range
 if (mpi_node()==0) print *,"   - dk_waves", this%dk_waves
 if (mpi_node()==0) print *,"   - tolerance", this%tolerance

end subroutine init_zpulse_planewaves_cyl_modes
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
! Check if parameters are good to use with current dimensionality
subroutine check_dimensionality_zpulse_planewaves_cyl_modes( this )

 implicit none

 class( t_zpulse_planewaves_cyl_modes ), intent(inout) :: this

 ! nothing to do

end subroutine check_dimensionality_zpulse_planewaves_cyl_modes
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
! Reading input
subroutine read_input_zpulse_planewaves_cyl_modes( this, input_file, g_space, bnd_con, periodic, grid, sim_options )

 implicit none

 class( t_zpulse_planewaves_cyl_modes ), intent(inout) :: this
 class( t_input_file ), intent(inout) :: input_file
 type( t_space ), intent(in) :: g_space
 class (t_emf_bound), intent(in) :: bnd_con
 logical, dimension(:), intent(in) :: periodic
 class( t_grid ), intent(in) :: grid
 type( t_options ), intent(in) :: sim_options

 real(p_double) :: aux

 ! Call superclass method
 call this % t_zpulse_planewaves % read_input( input_file, g_space, bnd_con, periodic, grid, sim_options )

 ! Adjust integration range
 this%per1_int_range(1) = 0.0_p_double
 this%per1_int_range(2) = this%per1_int_range(2) * (n_sig_gaussian+1) / n_sig_gaussian
 this%dk_waves(2) = this%dk_waves(2) / 2 / sqrt(2.0)

 ! Only longitudinal propagation
 if ( ( this%angles(1) /= 0 .and. this%angles(1) /= pi ) .or. this%angles(2) /= 0 ) then
  if ( mpi_node() == 0 ) then
   write(0,*) ""
   write(0,*) "   Error reading zpulse_planewaves parameters"
   write(0,*) "   only longitudinal propagation allowed when running with cylindrical modes"
   write(0,*) "   aborting..."
  endif
  stop
 endif
 if ( this%angles(1) == pi ) then
  aux = this%lon_int_range(1)
  this%lon_int_range(1) = - this%lon_int_range(2)
  this%lon_int_range(2) = - aux
 endif

 ! Focus on the axis
 if ( this%focus(2) /= 0 .or. this%focus(3) /= 0 ) then
  if ( mpi_node() == 0 ) then
   write(0,*) ""
   write(0,*) "   Error reading zpulse_planewaves parameters"
   write(0,*) "   pulse focus must be on axis when running with cylindrical modes"
   write(0,*) "   aborting..."
  endif
  stop
 endif

 ! Take care of normalization
 this%normalization = 1.0_p_double

 ! Calculate maximum mode needed
 select case ( this%per_type )
  case (k_per_gaussian)
   this%m_pulse = 1
  case (k_per_laguerre_gaussian)
   this%m_pulse = abs( this%per_tem_mode(2) ) + 1
 end select
 select case ( this%pulse_type )
  case (type_custom_full)
   this%m_pulse = int( this%per2_int_range(2) )
 end select

 if (mpi_node()==0) print *,"this%m_pulse", this%m_pulse

end subroutine read_input_zpulse_planewaves_cyl_modes
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
! Launch quasi3d pulse normal
subroutine launch_zpulse_planewaves_cyl_modes_normal( this, e_re, b_re, e_im, b_im, t, nx_p_min, g_space, no_co )

 implicit none

 integer, parameter :: rank = 2

 class( t_zpulse_planewaves_cyl_modes ), intent(inout) :: this
 type( t_vdf ), dimension(this%m_pulse+1), intent(inout) :: e_re, b_re, e_im, b_im
 real(p_double), intent(in) :: t
 integer, dimension(:), intent(in) :: nx_p_min
 type( t_space ), intent(in) :: g_space
 class( t_node_conf ), intent(in) :: no_co ! node configuration

 ! local variables
 real(p_double), dimension(2,rank) :: g_x_range
 real(p_double), dimension(rank) :: ldx

 real(p_double) :: zmin, dz_2
 real(p_double) :: rmin, r, r_2, dr_2

 ! other variables
 integer :: ir, iz, n, in, m, s, iw, total_n_waves
 real(p_double) :: kz, kr, omega, aux_boost, phi0, max_amp
 complex(p_double) :: f_k, amp, exp0
 real(p_double) :: a0_amp, a0_arg, psi_z, psi_z_2

 real(p_double) :: cpsi, spsi, cpsi_2, spsi_2
 real(p_double) :: jp1, j0_2, j0, jp1_2, jp2_2, jm1, jm1_2, jm2_2
 complex(p_double) :: expphi, exppsip, exppsim, exppsip_2, exppsim_2
 integer :: border_x1_up, border_x1_down, border_x2_up, border_x2_down

 ! executable statements
 call get_x_bnd( g_space, g_x_range )
 ldx(1) = b_re(1)%dx_(1)
 ldx(2) = b_re(1)%dx_(2)

 ! spatial coordinates
 ! XMIN ALREADY INCLUDES FOCUS POSITION
 zmin = g_x_range(p_lower,1) + real(nx_p_min(1)-1, p_double)*ldx(1) - this%focus(1)
 rmin = g_x_range(p_lower,2) + real(nx_p_min(2)-1, p_double)*ldx(2)
 dz_2 = ldx(1)/2.0_p_double
 dr_2 = ldx(2)/2.0_p_double

 ! get border if move window
 border_x1_up = b_re(1)%ubound(2)
 border_x1_down = b_re(1)%lbound(2)
 if ( this%wall_injection(1,1) .and. my_ngp( no_co, 1) == 1 ) then
  border_x1_down = 2
 endif
 if ( this%wall_injection(1,2) .and. my_ngp( no_co, 1) == nx( no_co, 1 ) ) then
  border_x1_up = b_re(1)%ubound(2) - 1
 endif
 border_x2_up = b_re(1)%ubound(3)
 border_x2_down = b_re(1)%lbound(3)
 if ( this%wall_injection(2,2) .and. my_ngp( no_co, 2) == nx( no_co, 2 ) ) then
  border_x2_up = b_re(1)%ubound(3) - 4
 endif

 ! store polarization
 phi0 = this%pol
 expphi = exp( im_unit*phi0 )

 ! get actual total number of waves to inject and iterate over them
 total_n_waves = size( this%kvec_arr, 1 )
 do iw = 1, total_n_waves

  n = abs( int( this%kvec_arr(iw,1) ) )
  s = sign( 1, int( this%kvec_arr(iw,1) ) )
  kz = this%kvec_arr(iw,2)
  kr = this%kvec_arr(iw,3)
  f_k = this%evec_arr(iw,1)

  ! calculate omega
  omega = real( this%evec_arr(iw,2) )

  ! boost pulse
  if (this%if_boost) then
   aux_boost = omega
   omega = this%gamma * omega - sqrt(this%gamma**2-1) * kz
   kz = this%gamma * kz - sqrt(this%gamma**2-1) * aux_boost
  endif

  ! start integration scheme
  a0_amp = abs( f_k )
  a0_arg = atan2( dimag( f_k ), real( f_k ) )

  if ( n == 0 ) then
   ! for n = 0 and m = 1
   m = 1

   exp0 = exp( im_unit*s*kz*ldx(1) )
   ! iterate for every space point
   do ir = border_x2_down, border_x2_up
    ! calculate grid positions and exponentials for r
    r = rmin + real(ir-1, p_double) * ldx(2)
    r_2 = r + dr_2
    if (r_2 == 0) r_2 = 1e-6 ! TODO: not the most correct thing to do
    jp1 = bessel_jn(1,kr*r)
    j0 = bessel_jn(0,kr*r)
    j0_2 = bessel_jn(0,kr*r_2)
    jp1_2 = bessel_jn(1,kr*r_2)
    jp2_2 = bessel_jn(2,kr*r_2)
   exppsim = exp( im_unit * ( s * ( kz * zmin - omega * (t + this%lon_position) + a0_arg + this%phase0 ) ) )
   exppsim_2 = exp( im_unit * ( s * ( kz * (zmin+dz_2) - omega * (t + this%lon_position) + a0_arg + this%phase0 ) ) )
   do iz = border_x1_down, border_x1_up
    ! calculate grid positions and exponentials for z
    exppsim = exppsim * exp0
    exppsim_2 = exppsim * exp0
    cpsi = real( exppsim )
    cpsi_2 = real( exppsim_2 )
    spsi = dimag( exppsim )
    spsi_2 = dimag( exppsim_2 )
   ! --------------------------------------------------

    ! e1 at (z_2, r )
    amp = (a0_amp*kr*jp1*spsi_2)/(expphi)
    e_re(m)%f2(1,iz,ir) = e_re(m)%f2(1,iz,ir) + real( amp )
    e_im(m)%f2(1,iz,ir) = e_im(m)%f2(1,iz,ir) - dimag( amp )

    ! e2 at (z , r_2)
    amp = (a0_amp*kz*j0_2*cpsi)/(expphi)
    e_re(m)%f2(2,iz,ir) = e_re(m)%f2(2,iz,ir) + real( amp )
    e_im(m)%f2(2,iz,ir) = e_im(m)%f2(2,iz,ir) - dimag( amp )

    ! e3 at (z , r )
    amp = (im_unit*a0_amp*kz*j0*cpsi)/(expphi)
    e_re(m)%f2(3,iz,ir) = e_re(m)%f2(3,iz,ir) + real( amp )
    e_im(m)%f2(3,iz,ir) = e_im(m)%f2(3,iz,ir) - dimag( amp )

    ! b1 at (z , r_2)
    amp = (-im_unit*a0_amp*kz*kr*jp1_2*spsi)/(expphi*omega)
    b_re(m)%f2(1,iz,ir) = b_re(m)%f2(1,iz,ir) + real( amp )
    b_im(m)%f2(1,iz,ir) = b_im(m)%f2(1,iz,ir) - dimag( amp )

    ! b2 at (z_2, r )
    amp = (-im_unit*a0_amp*(kz**2*j0 + kr/r*jp1)*cpsi_2)/(expphi*omega)
    b_re(m)%f2(2,iz,ir) = b_re(m)%f2(2,iz,ir) + real( amp )
    b_im(m)%f2(2,iz,ir) = b_im(m)%f2(2,iz,ir) - dimag( amp )

    ! b3 at (z_2, r_2)
    amp = (a0_amp*((2*kz**2 + kr**2)*j0_2 - kr**2*jp2_2)*cpsi_2)/(2.*expphi*omega)
    b_re(m)%f2(3,iz,ir) = b_re(m)%f2(3,iz,ir) + real( amp )
    b_im(m)%f2(3,iz,ir) = b_im(m)%f2(3,iz,ir) - dimag( amp )

   enddo
   enddo

  else

   exp0 = exp( im_unit*s*kz*ldx(1) )
   ! iterate for every space point
   do ir = border_x2_down, border_x2_up
    ! calculate grid positions and exponentials for r
    r = rmin + real(ir-1, p_double) * ldx(2)
    r_2 = r + dr_2
    if (r_2 == 0) r_2 = 1e-6 ! TODO: not the most correct thing to do
    jp1 = bessel_jn( n + 1, kr*r )
    j0_2 = bessel_jn( n , kr*r_2 )
    j0 = bessel_jn( n , kr*r )
    jp1_2 = bessel_jn( n + 1, kr*r_2 )
    jp2_2 = bessel_jn( n + 2, kr*r_2 )
    jm1 = bessel_jn( n - 1, kr*r )
    jm1_2 = bessel_jn( n - 1, kr*r_2 )
    jm2_2 = bessel_jn( n - 2, kr*r_2 )
   exppsip = exp( im_unit * ( phi0 + s * ( kz * zmin - omega * (t + this%lon_position) + a0_arg + this%phase0 ) ) )
   exppsim = exp( im_unit * ( - phi0 + s * ( kz * zmin - omega * (t + this%lon_position) + a0_arg + this%phase0 ) ) )
   exppsip_2 = exp( im_unit * ( phi0 + s * ( kz * (zmin+dz_2) - omega * (t + this%lon_position) + a0_arg + this%phase0 ) ) )
   exppsim_2 = exp( im_unit * ( - phi0 + s * ( kz * (zmin+dz_2) - omega * (t + this%lon_position) + a0_arg + this%phase0 ) ) )
   do iz = border_x1_down, border_x1_up
    ! calculate grid positions and exponentials for z
    exppsip = exppsip * exp0
    exppsim = exppsim * exp0
    exppsip_2 = exppsip_2 * exp0
    exppsim_2 = exppsim_2 * exp0
   ! --------------------------------------------------

   ! for n > 0 and m = n - 1
    m = n - 1
    if ( m == 0 ) m = this%m_pulse + 1
    ! e1 at (z_2, r )
    amp = (s*im_unit/2*a0_amp*exppsip_2*kr*jm1)
    e_re(m)%f2(1,iz,ir) = e_re(m)%f2(1,iz,ir) + real( amp )
    e_im(m)%f2(1,iz,ir) = e_im(m)%f2(1,iz,ir) - dimag( amp )

    ! e2 at (z , r_2)
    amp = (a0_amp*exppsip*kz*j0_2)/(2.)
    e_re(m)%f2(2,iz,ir) = e_re(m)%f2(2,iz,ir) + real( amp )
    e_im(m)%f2(2,iz,ir) = e_im(m)%f2(2,iz,ir) - dimag( amp )

    ! e3 at (z , r )
    amp = (-im_unit/2*a0_amp*exppsip*kz*j0)
    e_re(m)%f2(3,iz,ir) = e_re(m)%f2(3,iz,ir) + real( amp )
    e_im(m)%f2(3,iz,ir) = e_im(m)%f2(3,iz,ir) - dimag( amp )


    ! b1 at (z , r_2)
    amp = -0.5*s*(a0_amp*exppsip*kz*kr*jm1_2)/omega
    b_re(m)%f2(1,iz,ir) = b_re(m)%f2(1,iz,ir) + real( amp )
    b_im(m)%f2(1,iz,ir) = b_im(m)%f2(1,iz,ir) - dimag( amp )

    ! b2 at (z_2, r )
    amp = (im_unit/2*a0_amp*exppsip_2*(kr/r*(-1 + n)*jm1 + kz**2*j0))/(omega)
    b_re(m)%f2(2,iz,ir) = b_re(m)%f2(2,iz,ir) + real( amp )
    b_im(m)%f2(2,iz,ir) = b_im(m)%f2(2,iz,ir) - dimag( amp )

    ! b3 at (z_2, r_2)
    amp = -0.25*(a0_amp*exppsip_2*(kr**2*jm2_2 - (2*kz**2 + kr**2)*j0_2))/omega
    b_re(m)%f2(3,iz,ir) = b_re(m)%f2(3,iz,ir) + real( amp )
    b_im(m)%f2(3,iz,ir) = b_im(m)%f2(3,iz,ir) - dimag( amp )

   ! for n > 0 and m = n + 1
    m = n + 1
    ! e1 at (z_2, r )
    amp = (-0.5*s*im_unit*a0_amp*exppsim_2*kr*jp1)
    e_re(m)%f2(1,iz,ir) = e_re(m)%f2(1,iz,ir) + real( amp )
    e_im(m)%f2(1,iz,ir) = e_im(m)%f2(1,iz,ir) - dimag( amp )

    ! e2 at (z , r_2)
    amp = (a0_amp*exppsim*kz*j0_2)/(2.)
    e_re(m)%f2(2,iz,ir) = e_re(m)%f2(2,iz,ir) + real( amp )
    e_im(m)%f2(2,iz,ir) = e_im(m)%f2(2,iz,ir) - dimag( amp )

    ! e3 at (z , r )
    amp = (0.5*im_unit*a0_amp*exppsim*kz*j0)
    e_re(m)%f2(3,iz,ir) = e_re(m)%f2(3,iz,ir) + real( amp )
    e_im(m)%f2(3,iz,ir) = e_im(m)%f2(3,iz,ir) - dimag( amp )


    ! b1 at (z , r_2)
    amp = -0.5*s*(a0_amp*exppsim*kz*kr*jp1_2)/omega
    b_re(m)%f2(1,iz,ir) = b_re(m)%f2(1,iz,ir) + real( amp )
    b_im(m)%f2(1,iz,ir) = b_im(m)%f2(1,iz,ir) - dimag( amp )

    ! b2 at (z_2, r )
    amp = (-im_unit/2*a0_amp*exppsim_2*(kz**2*j0 + kr/r*(1 + n)*jp1))/(omega)
    b_re(m)%f2(2,iz,ir) = b_re(m)%f2(2,iz,ir) + real( amp )
    b_im(m)%f2(2,iz,ir) = b_im(m)%f2(2,iz,ir) - dimag( amp )

    ! b3 at (z_2, r_2)
    amp = (a0_amp*exppsim_2*((2*kz**2 + kr**2)*j0_2 - kr**2*jp2_2))/(4.*omega)
    b_re(m)%f2(3,iz,ir) = b_re(m)%f2(3,iz,ir) + real( amp )
    b_im(m)%f2(3,iz,ir) = b_im(m)%f2(3,iz,ir) - dimag( amp )

   enddo
   enddo

  endif

  ! testing, debug and control
  if ( control .and. mpi_node() == 0 ) then
   if ( mod( iw, total_n_waves/20) == 0 ) then
    write(0,*) " - injected ", iw * 100.0 /total_n_waves, " % of waves"
   endif
  endif
 enddo

end subroutine launch_zpulse_planewaves_cyl_modes_normal
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
! Launch pulse routine
subroutine launch_zpulse_planewaves_cyl_modes( this, emf, g_space, nx_p_min, g_nx, t, dt, no_co )

 implicit none

 class( t_zpulse_planewaves_cyl_modes ), intent(inout) :: this
 class( t_emf ), intent(inout) :: emf
 type( t_space ), intent(in) :: g_space
 integer, intent(in), dimension(:) :: nx_p_min, g_nx
 real(p_double), intent(in) :: t
 real(p_double), intent(in) :: dt
 class( t_node_conf ), intent(in) :: no_co

 ! field variables
 type( t_vdf ), dimension(this%m_pulse+1) :: e_pulse_re, b_pulse_re, e_pulse_im, b_pulse_im

 ! extra variables
 integer :: m

 ! testing, debug and control
 real(p_double) :: start_time, end_time
 real(p_double), dimension(2,2) :: g_x_range

 select type( emf )
 class is( t_emf_cyl_modes )

  ! Check if simulation has enough m modes for pulse
  if( ubound( emf%b_cyl_m%pf_re, 1 ) < this%m_pulse ) then
   if ( mpi_node() == 0 ) then
    write(0,*) ""
    write(0,*) "   ERROR: num cyl for simulation modes must be >= pulse modes"
    write(0,*) "   n_modes: ", ubound( emf%b_cyl_m%pf_re, 1 )
    write(0,*) "   m_pulse: ", this%m_pulse
   endif
   stop
  endif

  ! Launch
  if ( this%if_launch .and. .not. this%was_launched .and. ( ( t >= this%launch_time .and. .not. any(this%wall_injection) ) .or. ( this%initial_injection .and. any(this%wall_injection) ) ) ) then

   if (mpi_node()==0) print *," - injecting zpulse_planewaves"

   ! testing, debug and control
   if ( control ) then
    call cpu_time(start_time)
   endif

   ! create new vdfs to hold the pulse field
   do m = 1, this%m_pulse+1
    call e_pulse_re(m) % new( emf%e, zero = .true. )
    call b_pulse_re(m) % new( emf%b, zero = .true. )
    call e_pulse_im(m) % new( emf%e, zero = .true. )
    call b_pulse_im(m) % new( emf%b, zero = .true. )
   enddo

   call launch_zpulse_planewaves_cyl_modes_normal( this, e_pulse_re, b_pulse_re, e_pulse_im, b_pulse_im, t, nx_p_min, g_space, no_co )

   ! take care of circular polarisation
   if (this%pol_type /= 0) then

    this%pol = this%pol + real( pi_2, p_double )
    this%phase0 = this%phase0 + real( sign( pi_2, real( this%pol_type, p_double )) , p_double )

    call launch_zpulse_planewaves_cyl_modes_normal( this, e_pulse_re, b_pulse_re, e_pulse_im, b_pulse_im, t, nx_p_min, g_space, no_co )

    this%pol = this%pol - real( pi_2, p_double )
    this%phase0 = this%phase0 - real( sign( pi_2, real( this%pol_type, p_double )) , p_double )

   endif

   ! add laser pulse to fields
   do m = 1, this%m_pulse
    call add( emf%e_cyl_m%pf_re(m), e_pulse_re(m) )
    call add( emf%b_cyl_m%pf_re(m), b_pulse_re(m) )
    call add( emf%e_cyl_m%pf_im(m), e_pulse_im(m) )
    call add( emf%b_cyl_m%pf_im(m), b_pulse_im(m) )
   enddo
   ! add component 0
   call add( emf%e_cyl_m%pf_re(0), e_pulse_re(this%m_pulse+1) )
   call add( emf%b_cyl_m%pf_re(0), b_pulse_re(this%m_pulse+1) )

   ! free b_pulse and e_pulse vdfs
   do m = 1, this%m_pulse+1
    call e_pulse_re(m) % cleanup()
    call b_pulse_re(m) % cleanup()
    call e_pulse_im(m) % cleanup()
    call b_pulse_im(m) % cleanup()
   enddo

   ! pulse has been launched, turn off if_launch
   this%was_launched = .true.

   if (mpi_node()==0) print *," - done injecting zpulse_planewaves"

   ! testing, debug and control
   if ( control ) then
    call cpu_time(end_time)
    call get_x_bnd( g_space, g_x_range )
    write(0,*) "<control-loc> ", "Node:", mpi_node(), end_time-start_time
   endif

  endif

  ! launch zpulse from wall
  if ( this%if_launch .and. any( this%wall_injection ) .and. (t > this%t_wall_start) .and. (t < this%t_wall_stop) ) then

   ! create new vdfs to hold the pulse field
   do m = 1, this%m_pulse+1
    call e_pulse_re(m) % new( emf%e, zero = .true. )
    call b_pulse_re(m) % new( emf%b, zero = .true. )
    call e_pulse_im(m) % new( emf%e, zero = .true. )
    call b_pulse_im(m) % new( emf%b, zero = .true. )
   enddo

   ! inject from the left wall
   if ( my_ngp( no_co, 1 ) == 1 .and. this%wall_injection(1,1) ) then
    call launch_zpulse_planewaves_cyl_modes_wall_left( this, e_pulse_re, b_pulse_re, e_pulse_im, b_pulse_im, g_space, nx_p_min, t, dt, no_co )
    ! take care of circular polarisation
    if (this%pol_type /= 0) then
     this%pol = this%pol + real( pi_2, p_double )
     this%phase0 = this%phase0 + real( sign( pi_2, real( this%pol_type, p_double )) , p_double )
     call launch_zpulse_planewaves_cyl_modes_wall_left( this, e_pulse_re, b_pulse_re, e_pulse_im, b_pulse_im, g_space, nx_p_min, t, dt, no_co )
     this%pol = this%pol - real( pi_2, p_double )
     this%phase0 = this%phase0 - real( sign( pi_2, real( this%pol_type, p_double )) , p_double )
    endif
   endif

   ! inject from the right wall
   if ( my_ngp( no_co, 1 ) == nx( no_co, 1 ) .and. this%wall_injection(1,2) ) then
    call launch_zpulse_planewaves_cyl_modes_wall_right( this, e_pulse_re, b_pulse_re, e_pulse_im, b_pulse_im, g_space, nx_p_min, t, dt, no_co )
    ! take care of circular polarisation
    if (this%pol_type /= 0) then
     this%pol = this%pol + real( pi_2, p_double )
     this%phase0 = this%phase0 + real( sign( pi_2, real( this%pol_type, p_double )) , p_double )
     call launch_zpulse_planewaves_cyl_modes_wall_right( this, e_pulse_re, b_pulse_re, e_pulse_im, b_pulse_im, g_space, nx_p_min, t, dt, no_co )
     this%pol = this%pol - real( pi_2, p_double )
     this%phase0 = this%phase0 - real( sign( pi_2, real( this%pol_type, p_double )) , p_double )
    endif
   endif

   ! inject from the upper wall
   if ( my_ngp( no_co, 2 ) == nx( no_co, 2 ) .and. this%wall_injection(2,2) ) then
    call launch_zpulse_planewaves_cyl_modes_wall_top( this, e_pulse_re, b_pulse_re, e_pulse_im, b_pulse_im, g_space, nx_p_min, t, dt, no_co )
    ! take care of circular polarisation
    if (this%pol_type /= 0) then
     this%pol = this%pol + real( pi_2, p_double )
     this%phase0 = this%phase0 + real( sign( pi_2, real( this%pol_type, p_double )) , p_double )
     call launch_zpulse_planewaves_cyl_modes_wall_top( this, e_pulse_re, b_pulse_re, e_pulse_im, b_pulse_im, g_space, nx_p_min, t, dt, no_co )
     this%pol = this%pol - real( pi_2, p_double )
     this%phase0 = this%phase0 - real( sign( pi_2, real( this%pol_type, p_double )) , p_double )
    endif
   endif

   ! add laser pulse to fields
   do m = 1, this%m_pulse
    call add( emf%e_cyl_m%pf_re(m), e_pulse_re(m) )
    call add( emf%b_cyl_m%pf_re(m), b_pulse_re(m) )
    call add( emf%e_cyl_m%pf_im(m), e_pulse_im(m) )
    call add( emf%b_cyl_m%pf_im(m), b_pulse_im(m) )
   enddo
   ! add component 0
   call add( emf%e_cyl_m%pf_re(0), e_pulse_im(this%m_pulse+1) )
   call add( emf%b_cyl_m%pf_re(0), b_pulse_im(this%m_pulse+1) )

   ! free b_pulse and e_pulse vdfs
   do m = 1, this%m_pulse+1
    call e_pulse_re(m) % cleanup()
    call b_pulse_re(m) % cleanup()
    call e_pulse_im(m) % cleanup()
    call b_pulse_im(m) % cleanup()
   enddo

  endif

 end select

end subroutine launch_zpulse_planewaves_cyl_modes
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
! Wave vector weight function
function k_function_zpulse_planewaves_cyl_modes( this, k_lon, k_per1, k_per2 ) result( f_k )

 implicit none
 class( t_zpulse_planewaves_cyl_modes ), intent(inout) :: this
 real(p_double), intent(in) :: k_lon, k_per1, k_per2
 real(p_double) :: kz, kr, m
 complex(p_double) :: f_k

 ! auxiliar variables
 real(p_k_fparse), dimension(3) :: k_dbl

 ! parse data to new variables
 kz = k_lon
 kr = k_per1
 m = k_per2

 !f_k_lon_arg = sign( sqrt( abs( k_lon**2 + (k_per1**2+k_per2**2) * this%abouraddy_factor ) ), k_lon**2 + (k_per1**2+k_per2**2) * this%abouraddy_factor )
 if ( this%pulse_type == type_time ) then
  ! apply perpendicular profile
  f_k = this%k_per_function( kr, m )
  ! waves centered around \omega
  f_k = f_k * this%k_lon_function( sqrt( kz**2 + kr**2 ) )

 else if ( this%pulse_type == type_space ) then
  ! apply perpendicular profile
  f_k = this%k_per_function( kr, m )
  ! waves centered around klon
  f_k = f_k * this%k_lon_function( kz )

 else if ( this%pulse_type == type_vfocus_const ) then
  ! apply perpendicular profile
  f_k = this%k_per_function( kr, m )
  ! waves centered around isosurface of contant focus velocity
  f_k = f_k * this%k_lon_function( ( -kz * this%vfocus + sqrt( kz**2 + kr**2 ) ) / ( 1 - this%vfocus ) )

 else if ( this%pulse_type == type_vfocus_accel ) then
  if ( mpi_node() == 0 ) then
   write(0,*) ""
   write(0,*) "   Error reading zpulse_planewaves parameters"
   write(0,*) "   type vfocus_accel not yet on cyl_modes"
   write(0,*) "   aborting..."
  endif
  stop

 else if ( this%pulse_type == type_custom_full ) then
  k_dbl(1) = real( k_lon , p_k_fparse )
  k_dbl(2) = real( k_per1, p_k_fparse )
  k_dbl(3) = real( k_per2, p_k_fparse )
       f_k = eval( this%f_k_amp_math_func, k_dbl ) * exp( im_unit * eval( this%f_k_pha_math_func, k_dbl ) )

 else if ( this%pulse_type == type_custom_lon ) then
  if ( mpi_node() == 0 ) then
   write(0,*) ""
   write(0,*) "   Error reading zpulse_planewaves parameters"
   write(0,*) "   type custom_lon not yet on cyl_modes"
   write(0,*) "   aborting..."
  endif
  stop

 else if ( this%pulse_type == type_custom_per ) then
  if ( mpi_node() == 0 ) then
   write(0,*) ""
   write(0,*) "   Error reading zpulse_planewaves parameters"
   write(0,*) "   type custom_per not yet on cyl_modes"
   write(0,*) "   aborting..."
  endif
  stop

 endif


end function k_function_zpulse_planewaves_cyl_modes
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
! Transverse wave vector weight function
function k_per_function_zpulse_planewaves_cyl_modes( this, k_per1, k_per2 ) result( f_k_per )

 implicit none
 class( t_zpulse_planewaves_cyl_modes ), intent(in) :: this
 real(p_double), intent(in) :: k_per1, k_per2
 real(p_double) :: kr
 integer :: m
 complex(p_double) :: f_k_per

 ! auxiliar variables
 real(p_double) :: w0
 integer :: p, l

 ! parse data to new variables
 kr = k_per1
 m = int(k_per2)

 select case ( this%per_type )
  case (k_per_gaussian)

   if ( m == 0 ) then
    f_k_per = this%per_w0(1)**2 / 2.0 * exp( - kr**2*this%per_w0(1)**2 / 4 )
   else
    f_k_per = 0.0_p_double
   endif

  case (k_per_laguerre_gaussian)

   p = this%per_tem_mode(1)
   l = this%per_tem_mode(2)
   w0 = this%per_w0(1)

   if ( m == abs(l) ) then

    f_k_per = w0**2 * (-1)**(p) / (4*pi) * 2
    f_k_per = f_k_per * exp( -kr**2*w0**2/4 )
    f_k_per = f_k_per * (kr*w0/sqrt(2.0_p_double))**abs(l)
    f_k_per = f_k_per * laguerre( p, abs(l), kr**2*w0**2/2 )

    ! take care of normalization
    select case ( p )
     case ( 0 ); f_k_per = f_k_per / ( 2**(1.0 + 0.5*Abs(l))*Gamma(1.0 + 0.5*Abs(l)) )
     case ( 1 )
      select case ( l )
      case ( 0 ); f_k_per = f_k_per
      case default; f_k_per = f_k_per / ( -(2**(Abs(l)/2.)*Abs(l)*Gamma(Abs(l)/2.)) )
      end select
     case ( 2 ); f_k_per = f_k_per / ( 2**(1 + Abs(l)/2.)*Gamma(2 + Abs(l)/2.) )
     case ( 3 ); f_k_per = f_k_per / ( -(2**(1 + Abs(l)/2.)*Gamma(2 + Abs(l)/2.)) )
     case ( 4 ); f_k_per = f_k_per / ( 2**(Abs(l)/2.)*Gamma(3 + Abs(l)/2.) )
     case ( 5 ); f_k_per = f_k_per / ( -(2**(Abs(l)/2.)*Gamma(3 + Abs(l)/2.)) )
     case ( 6 ); f_k_per = f_k_per / ( (2**(Abs(l)/2.)*Gamma(4 + Abs(l)/2.))/3. )
     case ( 7 ); f_k_per = f_k_per / ( -0.3333333333333333*(2**(Abs(l)/2.)*Gamma(4 + Abs(l)/2.)) )
    end select

    if ( l < 0 ) f_k_per = f_k_per - 2 * im_unit * dimag(f_k_per)

   else
    f_k_per = 0.0_p_double
   endif

 end select

end function k_per_function_zpulse_planewaves_cyl_modes
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
! Launch quasi3d pulse from the left wall
subroutine launch_zpulse_planewaves_cyl_modes_wall_left( this, e_re, b_re, e_im, b_im, g_space, nx_p_min, t, dt, no_co )

 implicit none

 integer, parameter :: rank = 2

 class( t_zpulse_planewaves_cyl_modes ), intent(inout) :: this
 type( t_vdf ), dimension(this%m_pulse), intent(inout) :: e_re, b_re, e_im, b_im
 type( t_space ), intent(in) :: g_space ! global space information
 integer, intent(in), dimension(:) :: nx_p_min
 real( p_double ), intent(in) :: t ! simulation time
 real( p_double ), intent(in) :: dt ! time step
 class( t_node_conf ), intent(in) :: no_co ! node configuration

 real(p_double), dimension(2,rank) :: g_x_range
 real(p_double), dimension(rank) :: ldx

 integer :: ip, inj, iw, total_n_waves

 real(p_double) :: zmin, z, z_2, dz_2
 real(p_double) :: rmin, r, r_2, dr_2, rdtdz

 ! other variables
 integer :: ikz, ikr, ir, iz, n, in, m, s
 real(p_double) :: kz, kr, omega, aux_boost, phi0, max_amp
 complex(p_double) :: f_k, e2_amp, e3_amp, b2_amp, b3_amp
 real(p_double) :: a0_amp, a0_arg, psi_z, psi_z_2

 real(p_double) :: cpsi, spsi, cpsi_2, spsi_2
 real(p_double) :: jp1, j0_2, j0, jp1_2, jp2_2, jm1, jm1_2, jm2_2
 complex(p_double) :: expphi, exppsip, exppsim, exppsip_2, exppsim_2
 integer :: border_x1_up, border_x1_down, border_x2_up, border_x2_down

 logical :: save_cache

 ! executable statements
 call get_x_bnd( g_space, g_x_range )
 ldx(1) = e_re(1)%dx_(1)
 ldx(2) = e_re(1)%dx_(2)

 ! get border if move window
 border_x2_up = b_re(1)%ubound(3)
 border_x2_down = b_re(1)%lbound(3)
 border_x1_up = b_re(1)%ubound(2)
 border_x1_down = b_re(1)%lbound(2)

 rdtdz = dt/ldx(1)

 ! spatial coordinates
 ! XMIN ALREADY INCLUDES FOCUS POSITION
 zmin = g_x_range(p_lower,1) + real(nx_p_min(1)-1, p_double)*ldx(1) - this%focus(1)
 rmin = g_x_range(p_lower,2) + real(nx_p_min(2)-1, p_double)*ldx(2)
 dz_2 = ldx(1)/2.0_p_double
 dr_2 = ldx(2)/2.0_p_double

 ! store polarization
 phi0 = this%pol
 expphi = exp( im_unit*phi0 )

 ! iterate for every space point
 iz = 2

 ! calculate grid positions and exponentials for z
 z = zmin + real(iz-1, p_double) * ldx(1)
 z_2 = z + dz_2

 ! get actual total number of waves
 total_n_waves = size( this%kvec_arr, 1 )

 ! get cache
 save_cache = .false.
 if ( .not. allocated( this%bessel_cache_left ) ) then
  allocate( this%bessel_cache_left( total_n_waves, border_x2_down:border_x2_up, 8 ) )
  save_cache = .true.
 endif

 ! inject
 do iw = 1, total_n_waves

  n = abs( int( this%kvec_arr(iw,1) ) )
  s = sign( 1, int( this%kvec_arr(iw,1) ) )
  kz = this%kvec_arr(iw,2)
  kr = this%kvec_arr(iw,3)
  f_k = this%evec_arr(iw,1)

  ! calculate omega
  omega = real( this%evec_arr(iw,2) )

  ! boost pulse
  if (this%if_boost) then
   aux_boost = omega
   omega = this%gamma * omega - sqrt(this%gamma**2-1) * kz
   kz = this%gamma * kz - sqrt(this%gamma**2-1) * aux_boost
  endif

  ! start integration scheme
  a0_amp = abs( f_k )
  a0_arg = atan2( dimag( f_k ), real( f_k ) )

  if ( n == 0 ) then
   ! for n = 0 and m = 1
   m = 1

   expphi = exp( im_unit*phi0 )
   cpsi = cos( s * ( kz * z - omega * (t + this%lon_position - dt) + a0_arg + this%phase0 ) )
   cpsi_2 = cos( s * ( kz * z_2 - omega * (t + this%lon_position - dt) + a0_arg + this%phase0 ) )
   spsi = sin( s * ( kz * z - omega * (t + this%lon_position - dt) + a0_arg + this%phase0 ) )
   spsi_2 = sin( s * ( kz * z_2 - omega * (t + this%lon_position - dt) + a0_arg + this%phase0 ) )
   do ir = border_x2_down, border_x2_up
    ! calculate grid positions and exponentials for r
    r = rmin + real(ir-1, p_double) * ldx(2) + ldx(2)
    r_2 = r + dr_2
    if (r_2 == 0) r_2 = 1e-6 ! TODO: not the most correct thing to do
    if ( save_cache ) then
     jp1 = bessel_jn(1,kr*r)
     j0 = bessel_jn(0,kr*r)
     j0_2 = bessel_jn(0,kr*r_2)
     jp1_2 = bessel_jn(1,kr*r_2)
     jp2_2 = bessel_jn(2,kr*r_2)
     this%bessel_cache_left( iw, ir, 1 ) = jp1
     this%bessel_cache_left( iw, ir, 2 ) = j0
     this%bessel_cache_left( iw, ir, 3 ) = j0_2
     this%bessel_cache_left( iw, ir, 4 ) = jp1_2
     this%bessel_cache_left( iw, ir, 5 ) = jp2_2
    else
     jp1 = this%bessel_cache_left( iw, ir, 1 )
     j0 = this%bessel_cache_left( iw, ir, 2 )
     j0_2 = this%bessel_cache_left( iw, ir, 3 )
     jp1_2 = this%bessel_cache_left( iw, ir, 4 )
     jp2_2 = this%bessel_cache_left( iw, ir, 5 )
    endif

   ! --------------------------------------------------

    ! e2 at (z , r )
    e2_amp = (a0_amp*kz*j0_2*cpsi)/(expphi)
    ! e3 at (z , r )
    e3_amp = (im_unit*a0_amp*kz*j0*cpsi)/(expphi)
    ! b2 at (z_2, r )
    b2_amp = (-im_unit*a0_amp*(kz**2*j0 + kr/r*jp1)*cpsi_2)/(expphi*omega)
    ! b3 at (z_2, r_2)
    b3_amp = (a0_amp*((2*kz**2 + kr**2)*j0_2 - kr**2*jp2_2)*cpsi_2)/(2.*expphi*omega)

    ! e2: -d/dz (b3)
    e_re(m)%f2(2,iz,ir) = e_re(m)%f2(2,iz,ir) + real( rdtdz * b3_amp )
    e_im(m)%f2(2,iz,ir) = e_im(m)%f2(2,iz,ir) - dimag( rdtdz * b3_amp )
    ! e3: d/dz (b2)
    e_re(m)%f2(3,iz,ir) = e_re(m)%f2(3,iz,ir) + real( - rdtdz * b2_amp )
    e_im(m)%f2(3,iz,ir) = e_im(m)%f2(3,iz,ir) - dimag( - rdtdz * b2_amp )
    ! b2: d/dz (e3)
    b_re(m)%f2(2,iz-1,ir) = b_re(m)%f2(2,iz-1,ir) + real( - rdtdz * e3_amp )
    b_im(m)%f2(2,iz-1,ir) = b_im(m)%f2(2,iz-1,ir) - dimag( - rdtdz * e3_amp )
    ! b3: -d/dz (e2)
    b_re(m)%f2(3,iz-1,ir) = b_re(m)%f2(3,iz-1,ir) + real( rdtdz * e2_amp )
    b_im(m)%f2(3,iz-1,ir) = b_im(m)%f2(3,iz-1,ir) - dimag( rdtdz * e2_amp )

   enddo

  else

   exppsip = exp( im_unit * ( phi0 + s * ( kz * z - omega * (t + this%lon_position) + a0_arg + this%phase0 ) ) )
   exppsim = exp( im_unit * ( - phi0 + s * ( kz * z - omega * (t + this%lon_position) + a0_arg + this%phase0 ) ) )
   exppsip_2 = exp( im_unit * ( phi0 + s * ( kz * z_2 - omega * (t + this%lon_position) + a0_arg + this%phase0 ) ) )
   exppsim_2 = exp( im_unit * ( - phi0 + s * ( kz * z_2 - omega * (t + this%lon_position) + a0_arg + this%phase0 ) ) )

   do ir = border_x2_down, border_x2_up
    ! calculate grid positions and exponentials for r
    r = rmin + real(ir-1, p_double) * ldx(2)
    r_2 = r + dr_2
    if (r_2 == 0) r_2 = 1e-6 ! TODO: not the most correct thing to do
    if ( save_cache ) then
     jp1 = bessel_jn( n + 1, kr*r )
     j0_2 = bessel_jn( n , kr*r_2 )
     j0 = bessel_jn( n , kr*r )
     jp1_2 = bessel_jn( n + 1, kr*r_2 )
     jp2_2 = bessel_jn( n + 2, kr*r_2 )
     jm1 = bessel_jn( n - 1, kr*r )
     jm1_2 = bessel_jn( n - 1, kr*r_2 )
     jm2_2 = bessel_jn( n - 2, kr*r_2 )
     this%bessel_cache_left( iw, ir, 1 ) = jp1
     this%bessel_cache_left( iw, ir, 2 ) = j0_2
     this%bessel_cache_left( iw, ir, 3 ) = j0
     this%bessel_cache_left( iw, ir, 4 ) = jp1_2
     this%bessel_cache_left( iw, ir, 5 ) = jp2_2
     this%bessel_cache_left( iw, ir, 6 ) = jm1
     this%bessel_cache_left( iw, ir, 7 ) = jm1_2
     this%bessel_cache_left( iw, ir, 8 ) = jm2_2
    else
     jp1 = this%bessel_cache_left( iw, ir, 1 )
     j0_2 = this%bessel_cache_left( iw, ir, 2 )
     j0 = this%bessel_cache_left( iw, ir, 3 )
     jp1_2 = this%bessel_cache_left( iw, ir, 4 )
     jp2_2 = this%bessel_cache_left( iw, ir, 5 )
     jm1 = this%bessel_cache_left( iw, ir, 6 )
     jm1_2 = this%bessel_cache_left( iw, ir, 7 )
     jm2_2 = this%bessel_cache_left( iw, ir, 8 )
    endif

   ! --------------------------------------------------

   ! for n > 0 and m = n - 1
    m = n - 1
    if ( m == 0 ) m = this%m_pulse + 1

    ! e2 at (z , r )
    e2_amp = (a0_amp*exppsip*kz*j0_2)/(2.)
    ! e3 at (z , r )
    e3_amp = (-im_unit/2*a0_amp*exppsip*kz*j0)
    ! b2 at (z_2, r )
    b2_amp = (im_unit/2*a0_amp*exppsip_2*(kr/r*(-1 + n)*jm1 + kz**2*j0))/(omega)
    ! b3 at (z_2, r_2)
    b3_amp = -0.25*(a0_amp*exppsip_2*(kr**2*jm2_2 - (2*kz**2 + kr**2)*j0_2))/omega

    ! e2: -d/dz (b3)
    e_re(m)%f2(2,iz,ir) = e_re(m)%f2(2,iz,ir) + real( rdtdz * b3_amp )
    e_im(m)%f2(2,iz,ir) = e_im(m)%f2(2,iz,ir) - dimag( rdtdz * b3_amp )
    ! e3: d/dz (b2)
    e_re(m)%f2(3,iz,ir) = e_re(m)%f2(3,iz,ir) + real( - rdtdz * b2_amp )
    e_im(m)%f2(3,iz,ir) = e_im(m)%f2(3,iz,ir) - dimag( - rdtdz * b2_amp )
    ! b2: d/dz (e3)
    b_re(m)%f2(2,iz-1,ir) = b_re(m)%f2(2,iz-1,ir) + real( - rdtdz * e3_amp )
    b_im(m)%f2(2,iz-1,ir) = b_im(m)%f2(2,iz-1,ir) - dimag( - rdtdz * e3_amp )
    ! b3: -d/dz (e2)
    b_re(m)%f2(3,iz-1,ir) = b_re(m)%f2(3,iz-1,ir) + real( rdtdz * e2_amp )
    b_im(m)%f2(3,iz-1,ir) = b_im(m)%f2(3,iz-1,ir) - dimag( rdtdz * e2_amp )


   ! for n > 0 and m = n + 1
    m = n + 1

    ! e2 at (z , r )
    e2_amp = (a0_amp*exppsim*kz*j0_2)/(2.)
    ! e3 at (z , r )
    e3_amp = (0.5*im_unit*a0_amp*exppsim*kz*j0)
    ! b2 at (z_2, r )
    b2_amp = (-im_unit/2*a0_amp*exppsim_2*(kz**2*j0 + kr/r*(1 + n)*jp1))/(omega)
    ! b3 at (z_2, r_2)
    b3_amp = (a0_amp*exppsim_2*((2*kz**2 + kr**2)*j0_2 - kr**2*jp2_2))/(4.*omega)

    ! e2: -d/dz (b3)
    e_re(m)%f2(2,iz,ir) = e_re(m)%f2(2,iz,ir) + real( rdtdz * b3_amp )
    e_im(m)%f2(2,iz,ir) = e_im(m)%f2(2,iz,ir) - dimag( rdtdz * b3_amp )
    ! e3: d/dz (b2)
    e_re(m)%f2(3,iz,ir) = e_re(m)%f2(3,iz,ir) + real( - rdtdz * b2_amp )
    e_im(m)%f2(3,iz,ir) = e_im(m)%f2(3,iz,ir) - dimag( - rdtdz * b2_amp )
    ! b2: d/dz (e3)
    b_re(m)%f2(2,iz-1,ir) = b_re(m)%f2(2,iz-1,ir) + real( - rdtdz * e3_amp )
    b_im(m)%f2(2,iz-1,ir) = b_im(m)%f2(2,iz-1,ir) - dimag( - rdtdz * e3_amp )
    ! b3: -d/dz (e2)
    b_re(m)%f2(3,iz-1,ir) = b_re(m)%f2(3,iz-1,ir) + real( rdtdz * e2_amp )
    b_im(m)%f2(3,iz-1,ir) = b_im(m)%f2(3,iz-1,ir) - dimag( rdtdz * e2_amp )

   enddo

  endif


 enddo



end subroutine launch_zpulse_planewaves_cyl_modes_wall_left
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
! Launch quasi3d pulse from the right wall
subroutine launch_zpulse_planewaves_cyl_modes_wall_right( this, e_re, b_re, e_im, b_im, g_space, nx_p_min, t, dt, no_co )

 implicit none

 integer, parameter :: rank = 2

 class( t_zpulse_planewaves_cyl_modes ), intent(inout) :: this
 type( t_vdf ), dimension(this%m_pulse), intent(inout) :: e_re, b_re, e_im, b_im
 type( t_space ), intent(in) :: g_space ! global space information
 integer, intent(in), dimension(:) :: nx_p_min
 real( p_double ), intent(in) :: t ! simulation time
 real( p_double ), intent(in) :: dt ! time step
 class( t_node_conf ), intent(in) :: no_co ! node configuration

 real(p_double), dimension(2,rank) :: g_x_range
 real(p_double), dimension(rank) :: ldx

 integer :: ip, inj, iw, total_n_waves

 real(p_double) :: zmin, z, z_2, dz_2
 real(p_double) :: rmin, r, r_2, dr_2, rdtdz

 ! other variables
 integer :: ikz, ikr, ir, iz, n, in, m, s
 real(p_double) :: kz, kr, omega, aux_boost, phi0, max_amp
 complex(p_double) :: f_k, e2_amp, e3_amp, b2_amp, b3_amp
 real(p_double) :: a0_amp, a0_arg, psi_z, psi_z_2

 real(p_double) :: cpsi, spsi, cpsi_2, spsi_2
 real(p_double) :: jp1, j0_2, j0, jp1_2, jp2_2, jm1, jm1_2, jm2_2
 complex(p_double) :: expphi, exppsip, exppsim, exppsip_2, exppsim_2
 integer :: border_x1_up, border_x1_down, border_x2_up, border_x2_down

 ! executable statements
 call get_x_bnd( g_space, g_x_range )
 ldx(1) = e_re(1)%dx_(1)
 ldx(2) = e_re(1)%dx_(2)

 ! get border if move window
 border_x2_up = b_re(1)%ubound(3)
 border_x2_down = b_re(1)%lbound(3)
 border_x1_up = b_re(1)%ubound(2)
 border_x1_down = b_re(1)%lbound(2)

 rdtdz = dt/ldx(1)

 ! spatial coordinates
 ! XMIN ALREADY INCLUDES FOCUS POSITION
 zmin = g_x_range(p_lower,1) + real(nx_p_min(1)-1, p_double)*ldx(1) - this%focus(1)
 rmin = g_x_range(p_lower,2) + real(nx_p_min(2)-1, p_double)*ldx(2)
 dz_2 = ldx(1)/2.0_p_double
 dr_2 = ldx(2)/2.0_p_double

 ! store polarization
 phi0 = this%pol
 expphi = exp( im_unit*phi0 )

 ! iterate for every space point
 iz = e_re(1)%nx_(1) - 2
 ! calculate grid positions and exponentials for z
 z = zmin + real(iz-1, p_double) * ldx(1)
 z_2 = z + dz_2

 ! get actual total number of waves to inject and iterate over them
 total_n_waves = size( this%kvec_arr, 1 )
 do iw = 1, total_n_waves

  n = abs( int( this%kvec_arr(iw,1) ) )
  s = sign( 1, int( this%kvec_arr(iw,1) ) )
  kz = this%kvec_arr(iw,2)
  kr = this%kvec_arr(iw,3)
  f_k = this%evec_arr(iw,1)

  ! calculate omega
  omega = real( this%evec_arr(iw,2) )

  ! boost pulse
  if (this%if_boost) then
   aux_boost = omega
   omega = this%gamma * omega - sqrt(this%gamma**2-1) * kz
   kz = this%gamma * kz - sqrt(this%gamma**2-1) * aux_boost
  endif

  ! start integration scheme
  a0_amp = abs( f_k )
  a0_arg = atan2( dimag( f_k ), real( f_k ) )

  if ( n == 0 ) then
   ! for n = 0 and m = 1
   m = 1

   expphi = exp( im_unit*phi0 )
   cpsi = cos( s * ( kz * z - omega * (t + this%lon_position - dt) + a0_arg + this%phase0 ) )
   cpsi_2 = cos( s * ( kz * z_2 - omega * (t + this%lon_position - dt) + a0_arg + this%phase0 ) )
   spsi = sin( s * ( kz * z - omega * (t + this%lon_position - dt) + a0_arg + this%phase0 ) )
   spsi_2 = sin( s * ( kz * z_2 - omega * (t + this%lon_position - dt) + a0_arg + this%phase0 ) )
   do ir = border_x2_down, border_x2_up
    ! calculate grid positions and exponentials for r
    r = rmin + real(ir-1, p_double) * ldx(2) + ldx(2)
    r_2 = r + dr_2
    if (r_2 == 0) r_2 = 1e-6 ! TODO: not the most correct thing to do
    jp1 = bessel_jn(1,kr*r)
    j0 = bessel_jn(0,kr*r)
    j0_2 = bessel_jn(0,kr*r_2)
    jp1_2 = bessel_jn(1,kr*r_2)
    jp2_2 = bessel_jn(2,kr*r_2)
   ! --------------------------------------------------

    ! e2 at (z , r )
    e2_amp = (a0_amp*kz*j0_2*cpsi)/(expphi)
    ! e3 at (z , r )
    e3_amp = (im_unit*a0_amp*kz*j0*cpsi)/(expphi)
    ! b2 at (z_2, r )
    b2_amp = (-im_unit*a0_amp*(kz**2*j0 + kr/r*jp1)*cpsi_2)/(expphi*omega)
    ! b3 at (z_2, r_2)
    b3_amp = (a0_amp*((2*kz**2 + kr**2)*j0_2 - kr**2*jp2_2)*cpsi_2)/(2.*expphi*omega)

    ! e2: -d/dz (b3)
    e_re(m)%f2(2,iz,ir) = e_re(m)%f2(2,iz,ir) + real( rdtdz * b3_amp )
    e_im(m)%f2(2,iz,ir) = e_im(m)%f2(2,iz,ir) - dimag( rdtdz * b3_amp )
    ! e3: d/dz (b2)
    e_re(m)%f2(3,iz,ir) = e_re(m)%f2(3,iz,ir) + real( - rdtdz * b2_amp )
    e_im(m)%f2(3,iz,ir) = e_im(m)%f2(3,iz,ir) - dimag( - rdtdz * b2_amp )
    ! b2: d/dz (e3)
    b_re(m)%f2(2,iz-1,ir) = b_re(m)%f2(2,iz-1,ir) + real( - rdtdz * e3_amp )
    b_im(m)%f2(2,iz-1,ir) = b_im(m)%f2(2,iz-1,ir) - dimag( - rdtdz * e3_amp )
    ! b3: -d/dz (e2)
    b_re(m)%f2(3,iz-1,ir) = b_re(m)%f2(3,iz-1,ir) + real( rdtdz * e2_amp )
    b_im(m)%f2(3,iz-1,ir) = b_im(m)%f2(3,iz-1,ir) - dimag( rdtdz * e2_amp )

   enddo

  else

   exppsip = exp( im_unit * ( phi0 + s * ( kz * z - omega * (t + this%lon_position) + a0_arg + this%phase0 ) ) )
   exppsim = exp( im_unit * ( - phi0 + s * ( kz * z - omega * (t + this%lon_position) + a0_arg + this%phase0 ) ) )
   exppsip_2 = exp( im_unit * ( phi0 + s * ( kz * z_2 - omega * (t + this%lon_position) + a0_arg + this%phase0 ) ) )
   exppsim_2 = exp( im_unit * ( - phi0 + s * ( kz * z_2 - omega * (t + this%lon_position) + a0_arg + this%phase0 ) ) )

   do ir = border_x2_down, border_x2_up
    ! calculate grid positions and exponentials for r
    r = rmin + real(ir-1, p_double) * ldx(2)
    r_2 = r + dr_2
    if (r_2 == 0) r_2 = 1e-6 ! TODO: not the most correct thing to do
    jp1 = bessel_jn( n + 1, kr*r )
    j0_2 = bessel_jn( n , kr*r_2 )
    j0 = bessel_jn( n , kr*r )
    jp1_2 = bessel_jn( n + 1, kr*r_2 )
    jp2_2 = bessel_jn( n + 2, kr*r_2 )
    jm1 = bessel_jn( n - 1, kr*r )
    jm1_2 = bessel_jn( n - 1, kr*r_2 )
    jm2_2 = bessel_jn( n - 2, kr*r_2 )

   ! --------------------------------------------------

   ! for n > 0 and m = n - 1
    m = n - 1
    if ( m == 0 ) m = this%m_pulse + 1

    ! e2 at (z , r )
    e2_amp = (a0_amp*exppsip*kz*j0_2)/(2.)
    ! e3 at (z , r )
    e3_amp = (-im_unit/2*a0_amp*exppsip*kz*j0)
    ! b2 at (z_2, r )
    b2_amp = (im_unit/2*a0_amp*exppsip_2*(kr/r*(-1 + n)*jm1 + kz**2*j0))/(omega)
    ! b3 at (z_2, r_2)
    b3_amp = -0.25*(a0_amp*exppsip_2*(kr**2*jm2_2 - (2*kz**2 + kr**2)*j0_2))/omega

    ! e2: -d/dz (b3)
    e_re(m)%f2(2,iz,ir) = e_re(m)%f2(2,iz,ir) + real( rdtdz * b3_amp )
    e_im(m)%f2(2,iz,ir) = e_im(m)%f2(2,iz,ir) - dimag( rdtdz * b3_amp )
    ! e3: d/dz (b2)
    e_re(m)%f2(3,iz,ir) = e_re(m)%f2(3,iz,ir) + real( - rdtdz * b2_amp )
    e_im(m)%f2(3,iz,ir) = e_im(m)%f2(3,iz,ir) - dimag( - rdtdz * b2_amp )
    ! b2: d/dz (e3)
    b_re(m)%f2(2,iz-1,ir) = b_re(m)%f2(2,iz-1,ir) + real( - rdtdz * e3_amp )
    b_im(m)%f2(2,iz-1,ir) = b_im(m)%f2(2,iz-1,ir) - dimag( - rdtdz * e3_amp )
    ! b3: -d/dz (e2)
    b_re(m)%f2(3,iz-1,ir) = b_re(m)%f2(3,iz-1,ir) + real( rdtdz * e2_amp )
    b_im(m)%f2(3,iz-1,ir) = b_im(m)%f2(3,iz-1,ir) - dimag( rdtdz * e2_amp )


   ! for n > 0 and m = n + 1
    m = n + 1

    ! e2 at (z , r )
    e2_amp = (a0_amp*exppsim*kz*j0_2)/(2.)
    ! e3 at (z , r )
    e3_amp = (0.5*im_unit*a0_amp*exppsim*kz*j0)
    ! b2 at (z_2, r )
    b2_amp = (-im_unit/2*a0_amp*exppsim_2*(kz**2*j0 + kr/r*(1 + n)*jp1))/(omega)
    ! b3 at (z_2, r_2)
    b3_amp = (a0_amp*exppsim_2*((2*kz**2 + kr**2)*j0_2 - kr**2*jp2_2))/(4.*omega)

    ! e2: -d/dz (b3)
    e_re(m)%f2(2,iz,ir) = e_re(m)%f2(2,iz,ir) + real( rdtdz * b3_amp )
    e_im(m)%f2(2,iz,ir) = e_im(m)%f2(2,iz,ir) - dimag( rdtdz * b3_amp )
    ! e3: d/dz (b2)
    e_re(m)%f2(3,iz,ir) = e_re(m)%f2(3,iz,ir) + real( - rdtdz * b2_amp )
    e_im(m)%f2(3,iz,ir) = e_im(m)%f2(3,iz,ir) - dimag( - rdtdz * b2_amp )
    ! b2: d/dz (e3)
    b_re(m)%f2(2,iz-1,ir) = b_re(m)%f2(2,iz-1,ir) + real( - rdtdz * e3_amp )
    b_im(m)%f2(2,iz-1,ir) = b_im(m)%f2(2,iz-1,ir) - dimag( - rdtdz * e3_amp )
    ! b3: -d/dz (e2)
    b_re(m)%f2(3,iz-1,ir) = b_re(m)%f2(3,iz-1,ir) + real( rdtdz * e2_amp )
    b_im(m)%f2(3,iz-1,ir) = b_im(m)%f2(3,iz-1,ir) - dimag( rdtdz * e2_amp )

   enddo

  endif


 enddo



end subroutine launch_zpulse_planewaves_cyl_modes_wall_right
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
! Launch quasi3d pulse from the right wall
subroutine launch_zpulse_planewaves_cyl_modes_wall_top( this, e_re, b_re, e_im, b_im, g_space, nx_p_min, t, dt, no_co )

 implicit none

 integer, parameter :: rank = 2

 class( t_zpulse_planewaves_cyl_modes ), intent(inout) :: this
 type( t_vdf ), dimension(this%m_pulse), intent(inout) :: e_re, b_re, e_im, b_im
 type( t_space ), intent(in) :: g_space ! global space information
 integer, intent(in), dimension(:) :: nx_p_min
 real( p_double ), intent(in) :: t ! simulation time
 real( p_double ), intent(in) :: dt ! time step
 class( t_node_conf ), intent(in) :: no_co ! node configuration

 real(p_double), dimension(2,rank) :: g_x_range
 real(p_double), dimension(rank) :: ldx

 integer :: ip, inj, iw, total_n_waves

 real(p_double) :: zmin, z, z_2, dz_2, rdtdr
 real(p_double) :: rmin, r, r_2, dr_2

 ! other variables
 integer :: ikz, ikr, ir, iz, n, in, m, s
 real(p_double) :: kz, kr, omega, aux_boost, phi0, max_amp
 complex(p_double) :: f_k, e1_amp, e3_amp, b1_amp, b3_amp
 real(p_double) :: a0_amp, a0_arg, psi_z, psi_z_2

 real(p_double) :: cpsi, spsi, cpsi_2, spsi_2
 real(p_double) :: jp1, j0_2, j0, jp1_2, jp2_2, jm1, jm1_2, jm2_2
 complex(p_double) :: expphi, exppsip, exppsim, exppsip_2, exppsim_2, exp0
 integer :: border_x1_up, border_x1_down, border_x2_up, border_x2_down

 ! executable statements
 call get_x_bnd( g_space, g_x_range )
 ldx(1) = e_re(1)%dx_(1)
 ldx(2) = e_re(1)%dx_(2)

 ! get border if move window
 border_x2_up = b_re(1)%ubound(3)
 border_x2_down = b_re(1)%lbound(3)
 border_x1_up = b_re(1)%ubound(2)
 border_x1_down = b_re(1)%lbound(2)

 rdtdr = dt/ldx(2)

 ! spatial coordinates
 ! XMIN ALREADY INCLUDES FOCUS POSITION
 zmin = g_x_range(p_lower,1) + real(nx_p_min(1)-1, p_double)*ldx(1) - this%focus(1)
 rmin = g_x_range(p_lower,2) + real(nx_p_min(2)-1, p_double)*ldx(2)
 dz_2 = ldx(1)/2.0_p_double
 dr_2 = ldx(2)/2.0_p_double

 ! store polarization
 phi0 = this%pol
 expphi = exp( im_unit*phi0 )

 ! iterate for every space point
 ir = b_re(1)%ubound(3) - 4

 ! calculate grid positions and exponentials for r
 r = rmin + real(ir-1, p_double) * ldx(2)
 r_2 = r + dr_2

 ! get actual total number of waves to inject and iterate over them
 total_n_waves = size( this%kvec_arr, 1 )
 do iw = 1, total_n_waves

  n = abs( int( this%kvec_arr(iw,1) ) )
  s = sign( 1, int( this%kvec_arr(iw,1) ) )
  kz = this%kvec_arr(iw,2)
  kr = this%kvec_arr(iw,3)
  f_k = this%evec_arr(iw,1)

  if ( n == 0 ) then
   jp1 = bessel_jn(1,kr*r)
   j0 = bessel_jn(0,kr*r)
   j0_2 = bessel_jn(0,kr*r_2)
   jp1_2 = bessel_jn(1,kr*r_2)
   jp2_2 = bessel_jn(2,kr*r_2)
  else
   jp1 = bessel_jn( n + 1, kr*r )
   j0_2 = bessel_jn( n , kr*r_2 )
   j0 = bessel_jn( n , kr*r )
   jp1_2 = bessel_jn( n + 1, kr*r_2 )
   jp2_2 = bessel_jn( n + 2, kr*r_2 )
   jm1 = bessel_jn( n - 1, kr*r )
   jm1_2 = bessel_jn( n - 1, kr*r_2 )
   jm2_2 = bessel_jn( n - 2, kr*r_2 )
  endif

  ! calculate omega
  omega = real( this%evec_arr(iw,2) )

  ! boost pulse
  if (this%if_boost) then
   aux_boost = omega
   omega = this%gamma * omega - sqrt(this%gamma**2-1) * kz
   kz = this%gamma * kz - sqrt(this%gamma**2-1) * aux_boost
  endif

  ! start integration scheme
  a0_amp = abs( f_k )
  a0_arg = atan2( dimag( f_k ), real( f_k ) )

  if ( n == 0 ) then
   ! for n = 0 and m = 1
   m = 1

   exppsim = exp( im_unit * ( s * ( kz * zmin - omega * (t - dt + this%lon_position) + a0_arg + this%phase0 ) ) )
   exppsim_2 = exp( im_unit * ( s * ( kz * (zmin+dz_2) - omega * (t - dt + this%lon_position) + a0_arg + this%phase0 ) ) )
   exp0 = exp( im_unit*s*kz*ldx(1) )
   do iz = border_x1_down, border_x1_up
    ! calculate grid positions and exponentials for z
    exppsim = exppsim * exp0
    exppsim_2 = exppsim * exp0
    cpsi = real( exppsim )
    cpsi_2 = real( exppsim_2 )
    spsi = dimag( exppsim )
    spsi_2 = dimag( exppsim_2 )
   ! --------------------------------------------------

    ! e1 at (z_2, r )
    e1_amp = (a0_amp*kr*jp1*spsi_2)/(expphi)
    ! e3 at (z , r )
    e3_amp = (im_unit*a0_amp*kz*j0*cpsi)/(expphi)
    ! b1 at (z , r_2)
    b1_amp = (-im_unit*a0_amp*kz*kr*jp1_2*spsi)/(expphi*omega)
    ! b3 at (z_2, r_2)
    b3_amp = (a0_amp*((2*kz**2 + kr**2)*j0_2 - kr**2*jp2_2)*cpsi_2)/(2.*expphi*omega)

    ! e1: 1/r d/dr( r b3 )
    e_re(m)%f2(1,iz,ir) = e_re(m)%f2(1,iz,ir) + real( rdtdr * b3_amp )!* r_2 ) / r
    e_im(m)%f2(1,iz,ir) = e_im(m)%f2(1,iz,ir) - dimag( rdtdr * b3_amp )!* r_2 ) / r
    ! e3: - d/dr (b1)
    e_re(m)%f2(3,iz,ir) = e_re(m)%f2(3,iz,ir) + real( - rdtdr * b1_amp )
    e_im(m)%f2(3,iz,ir) = e_im(m)%f2(3,iz,ir) - dimag( - rdtdr * b1_amp )
    ! b1: - 1/r d/dr( r b3 )
    b_re(m)%f2(1,iz,ir) = b_re(m)%f2(1,iz,ir) + real( - rdtdr * e3_amp)! * r ) / r_2
    b_im(m)%f2(1,iz,ir) = b_im(m)%f2(1,iz,ir) - dimag( - rdtdr * e3_amp)! * r ) / r_2
    ! b3: d/dr (e1)
    b_re(m)%f2(3,iz,ir) = b_re(m)%f2(3,iz,ir) + real( rdtdr * e1_amp )
    b_im(m)%f2(3,iz,ir) = b_im(m)%f2(3,iz,ir) - dimag( rdtdr * e1_amp )

   enddo

  else

   ! iterate for every space point
   exppsip = exp( im_unit * ( phi0 + s * ( kz * zmin - omega * (t + this%lon_position) + a0_arg + this%phase0 ) ) )
   exppsim = exp( im_unit * ( - phi0 + s * ( kz * zmin - omega * (t + this%lon_position) + a0_arg + this%phase0 ) ) )
   exppsip_2 = exp( im_unit * ( phi0 + s * ( kz * (zmin+dz_2) - omega * (t + this%lon_position) + a0_arg + this%phase0 ) ) )
   exppsim_2 = exp( im_unit * ( - phi0 + s * ( kz * (zmin+dz_2) - omega * (t + this%lon_position) + a0_arg + this%phase0 ) ) )
   exp0 = exp( im_unit*s*kz*ldx(1) )
   do iz = border_x1_down, border_x1_up
    ! calculate grid positions and exponentials for z
    exppsip = exppsip * exp0
    exppsim = exppsim * exp0
    exppsip_2 = exppsip_2 * exp0
    exppsim_2 = exppsim_2 * exp0
   ! --------------------------------------------------

   ! for n > 0 and m = n - 1
    m = n - 1
    if ( m == 0 ) m = this%m_pulse + 1

    ! e1 at (z_2, r )
    e1_amp = (s*im_unit/2*a0_amp*exppsip_2*kr*jm1)
    ! e3 at (z , r )
    e3_amp = (-im_unit/2*a0_amp*exppsip*kz*j0)
    ! b1 at (z , r_2)
    b1_amp = -0.5*s*(a0_amp*exppsip*kz*kr*jm1_2)/omega
    ! b3 at (z_2, r_2)
    b3_amp = -0.25*(a0_amp*exppsip_2*(kr**2*jm2_2 - (2*kz**2 + kr**2)*j0_2))/omega

    ! e1: 1/r d/dr( r b3 )
    e_re(m)%f2(1,iz,ir) = e_re(m)%f2(1,iz,ir) + real( rdtdr * b3_amp * r_2 ) / r
    e_im(m)%f2(1,iz,ir) = e_im(m)%f2(1,iz,ir) - dimag( rdtdr * b3_amp * r_2 ) / r
    ! e3: - d/dr (b1)
    e_re(m)%f2(3,iz,ir) = e_re(m)%f2(3,iz,ir) + real( - rdtdr * b1_amp )
    e_im(m)%f2(3,iz,ir) = e_im(m)%f2(3,iz,ir) - dimag( - rdtdr * b1_amp )
    ! b1: - 1/r d/dr( r b3 )
    b_re(m)%f2(1,iz,ir) = b_re(m)%f2(1,iz,ir) + real( - rdtdr * e3_amp * r ) / r_2
    b_im(m)%f2(1,iz,ir) = b_im(m)%f2(1,iz,ir) - dimag( - rdtdr * e3_amp * r ) / r_2
    ! b3: d/dr (e1)
    b_re(m)%f2(3,iz,ir) = b_re(m)%f2(3,iz,ir) + real( rdtdr * e1_amp )
    b_im(m)%f2(3,iz,ir) = b_im(m)%f2(3,iz,ir) - dimag( rdtdr * e1_amp )

   ! for n > 0 and m = n + 1
    m = n + 1

    ! e1 at (z_2, r )
    e1_amp = (-0.5*s*im_unit*a0_amp*exppsim_2*kr*jp1)
    ! e3 at (z , r )
    e3_amp = (0.5*im_unit*a0_amp*exppsim*kz*j0)
    ! b1 at (z , r_2)
    b1_amp = -0.5*s*(a0_amp*exppsim*kz*kr*jp1_2)/omega
    ! b3 at (z_2, r_2)
    b3_amp = (a0_amp*exppsim_2*((2*kz**2 + kr**2)*j0_2 - kr**2*jp2_2))/(4.*omega)

    ! e1: 1/r d/dr( r b3 )
    e_re(m)%f2(1,iz,ir) = e_re(m)%f2(1,iz,ir) + real( rdtdr * b3_amp * r_2 ) / r
    e_im(m)%f2(1,iz,ir) = e_im(m)%f2(1,iz,ir) - dimag( rdtdr * b3_amp * r_2 ) / r
    ! e3: - d/dr (b1)
    e_re(m)%f2(3,iz,ir) = e_re(m)%f2(3,iz,ir) + real( - rdtdr * b1_amp )
    e_im(m)%f2(3,iz,ir) = e_im(m)%f2(3,iz,ir) - dimag( - rdtdr * b1_amp )
    ! b1: - 1/r d/dr( r b3 )
    b_re(m)%f2(1,iz,ir) = b_re(m)%f2(1,iz,ir) + real( - rdtdr * e3_amp * r ) / r_2
    b_im(m)%f2(1,iz,ir) = b_im(m)%f2(1,iz,ir) - dimag( - rdtdr * e3_amp * r ) / r_2
    ! b3: d/dr (e1)
    b_re(m)%f2(3,iz,ir) = b_re(m)%f2(3,iz,ir) + real( rdtdr * e1_amp )
    b_im(m)%f2(3,iz,ir) = b_im(m)%f2(3,iz,ir) - dimag( rdtdr * e1_amp )

   enddo

  endif



 enddo

end subroutine launch_zpulse_planewaves_cyl_modes_wall_top
!-----------------------------------------------------------------------------------------



end module m_zpulse_planewaves_cyl_modes
