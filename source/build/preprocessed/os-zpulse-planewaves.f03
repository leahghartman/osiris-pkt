# 1 "zpulse/os-zpulse-planewaves.f03"
# 1 "<built-in>" 1
# 1 "<built-in>" 3
# 467 "<built-in>" 3
# 1 "<command line>" 1
# 1 "<built-in>" 2
# 1 "zpulse/os-zpulse-planewaves.f03" 2
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
# 2 "zpulse/os-zpulse-planewaves.f03" 2
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
# 3 "zpulse/os-zpulse-planewaves.f03" 2

module m_zpulse_planewaves

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
# 7 "zpulse/os-zpulse-planewaves.f03" 2

use m_zpulse_std

use m_node_conf
use m_restart
use m_parameters
use m_math

use m_vdf_define
use m_vdf_math
use m_vdf

use m_emf_define, only: t_emf_bound, t_emf

use m_space
use m_grid_define
use m_fparser
use m_utilities

implicit none

private

! testing, debug and control
logical, parameter :: control = .true.

! longitudinal profiles
integer, parameter :: k_lon_plane = 0
integer, parameter :: k_lon_gaussian = 1
integer, parameter :: k_lon_sin2 = 2
integer, parameter :: k_lon_polynomial = 3

! perpendicular profiles
integer, parameter :: k_per_plane = 0
integer, parameter :: k_per_gaussian = 1
integer, parameter :: k_per_laguerre_gaussian = 2
integer, parameter :: k_per_hermite_gaussian = 3

! pulse_types
integer, parameter :: type_time = 0
integer, parameter :: type_space = 1
integer, parameter :: type_vfocus_const = 2
integer, parameter :: type_vfocus_accel = 3
integer, parameter :: type_custom_full = 4
integer, parameter :: type_custom_lon = 5
integer, parameter :: type_custom_per = 6

! just the imaginaty unity
complex(p_double), parameter :: im_unit = (0.0_p_double, 1.0_p_double)

! parameters for integration ranges
integer, parameter :: n_sig_gaussian = 20
integer, parameter :: n_sig_sin2 = 20
integer, parameter :: n_sig_polynomial = 20

type, extends( t_zpulse ) :: t_zpulse_planewaves

 ! new variables go here
 real(p_double), dimension(3) :: focus
 real(p_double), dimension(3) :: start
 real(p_double), dimension(2) :: angles
 real(p_double), dimension(3) :: beta_boost
 real(p_double) :: pol_elipse

 real(p_double) :: vfocus
 real(p_double) :: afocus
 real(p_double), dimension(3) :: dk_waves
 logical :: normalize
 real(p_double) :: normalization

 real(p_double) :: lon_position
 real(p_double), dimension(2) :: lon_int_range

 real(p_double), dimension(2) :: per1_int_range
 real(p_double), dimension(2) :: per2_int_range

 real(p_double) :: tolerance

 complex(p_double), dimension (:,:), pointer :: evec_arr, bvec_arr
 real(p_double), dimension (:,:), pointer :: kvec_arr

 logical, dimension(3,2) :: wall_injection
 logical :: initial_injection
 real(p_double) :: t_wall_start, t_wall_stop
 logical :: was_launched

 integer :: pulse_type
 type(t_fparser) :: f_k_amp_math_func, f_k_pha_math_func
 type(t_fparser), dimension(2):: jones_vec_real, jones_vec_imag

 logical :: incoherent

 logical :: print_waves

 contains

 procedure :: init => init_zpulse_planewaves
 procedure :: launch => launch_zpulse_planewaves
 procedure :: read_input => read_input_zpulse_planewaves
 procedure :: check_dimensionality => check_dimensionality_zpulse_planewaves

 procedure :: k_function => k_function_zpulse_planewaves
 procedure :: k_lon_function => k_lon_function_zpulse_planewaves
 procedure :: k_per_function => k_per_function_zpulse_planewaves
 procedure :: get_jones_vector => get_jones_vector_zpulse_planewaves

end type t_zpulse_planewaves

public :: t_zpulse_planewaves

public :: type_time, type_space, type_vfocus_const, type_vfocus_accel, type_custom_full, type_custom_lon, type_custom_per
public :: k_lon_plane, k_lon_gaussian, k_lon_sin2, k_lon_polynomial
public :: k_per_plane, k_per_gaussian, k_per_laguerre_gaussian, k_per_hermite_gaussian
public :: n_sig_gaussian
public :: im_unit

contains

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
! Initialize zpulse_planewaves
subroutine init_zpulse_planewaves( this, restart, t , dt, b, g_space, nx_p_min, no_co, interpolation )

 implicit none

 class( t_zpulse_planewaves ), intent(inout) :: this
 logical, intent(in) :: restart
 real( p_double ), intent(in) :: t, dt
 type( t_vdf ) , intent(in) :: b
 type( t_space ), intent(in) :: g_space
 integer, intent(in), dimension(:) :: nx_p_min
 class( t_node_conf ), intent(in) :: no_co
 integer, intent(in) :: interpolation

 logical :: if_launch
 integer :: ik1, ik2, ik3, n_wave_now, num_waves
 real(p_double), dimension(3,3) :: rot
 real(p_double) :: k_long, k_per1, k_per2, dk_long, dk_per1, dk_per2, k1, k2, k3, omega, max_amp
 complex(p_double) :: f_k
 real(p_double), dimension(3) :: x_per, xper, x_per1, xper1, x_per2, xper2
 complex(p_double), dimension(2) :: jones_vec
 real(p_double) :: beta, beta1, beta2, beta3, w_b, k1_b, k2_b, k3_b
 complex(p_double) :: e1_b, e2_b, e3_b, b1_b, b2_b, b3_b
 complex(p_double) :: e1_fac, e2_fac, e3_fac, b1_fac, b2_fac, b3_fac
 complex(p_double), dimension (:), pointer :: fk_arr


 ! Need to preserve if_launch value for this zpulse
 this%was_launched = .false.

 if ( restart ) then
  this%was_launched = .true.
 endif

 if ( .not. this%if_launch ) then
  return
 endif

 if (mpi_node()==0) print *,"   - calculating plane waves"

 if ( this%incoherent ) then
  call init_zpulse_planewaves_incoherent( this, restart, t , dt, b, g_space, nx_p_min, no_co, interpolation )
  return
 endif

 ! calculate plane waves
 n_wave_now = 1
 if ( p_x_dim == 2 ) then

  if ( this%lon_type == k_lon_plane ) then
   num_waves = 1
  else
   num_waves = int( ( this%lon_int_range(2) - this%lon_int_range(1) ) / this%dk_waves(1) + 10 )
  endif
  if ( this%per_type /= k_per_plane ) then
   num_waves = num_waves * int( ( this%per1_int_range(2) - this%per1_int_range(1) ) / this%dk_waves(2) + 10 )
  endif

  allocate( this%kvec_arr( num_waves, 3 ) )
  allocate( this%evec_arr( num_waves, 3 ) )
  allocate( this%bvec_arr( num_waves, 3 ) )
  if ( this%print_waves ) allocate( fk_arr( num_waves ) )

  ! trignometry
  rot(1,:) = (/ cos( this%angles(1) ), -sin( this%angles(1) ), 0.0_p_double /)
  rot(2,:) = (/ sin( this%angles(1) ), cos( this%angles(1) ), 0.0_p_double /)

  ! calculate maximum amplitude for tolerance
  ! assuming the largest contribution is the principal wave ...
  max_amp = this%tolerance * abs( this%a0 * this%k_function( this%omega0, 0.0_p_double, 0.0_p_double ) * this%dk_waves(1) * this%dk_waves(2) )
  if ( max_amp == 0 ) max_amp = this%tolerance

  ! iterate for every plane wave
  k_long = this%lon_int_range(1) - mod( this%lon_int_range(1), this%dk_waves(1) )
  if ( this%lon_type == k_lon_plane ) k_long = this%lon_int_range(1) - this%dk_waves(1)
  do while ( k_long < this%lon_int_range(2) )
   ! calculate direction of each wave - longitudinal
   k_long = k_long + this%dk_waves(1)
  k_per1 = this%per1_int_range(1) - mod( this%per1_int_range(1), this%dk_waves(2) ) - 2*this%dk_waves(2)
  do while ( k_per1 < this%per1_int_range(2) )
   ! calculate direction of each wave - perpendicular
   k_per1 = k_per1 + this%dk_waves(2)
   if ( this%per_type == k_per_plane ) k_per1 = 0.0_p_double
  ! --------------------------------------------------

   ! calculate intensity of each wave
   f_k = this%a0 * this%k_function( k_long, k_per1, 0.0_p_double ) * this%dk_waves(1) * this%dk_waves(2)

   ! tolerance for speep up
   if ( abs(f_k) < max_amp ) cycle

   ! calculate omega corresponding to this wave
   omega = sqrt( k_long**2 + k_per1**2 )

   ! take care of infinite pulse
   if ( this%lon_type == k_lon_plane ) then
     omega = this%omega0
    k_long = sqrt( this%omega0**2 - k_per1**2 )
   endif

   ! calculate direction of wave in rotated space
   k1 = k_long * rot(1,1) + k_per1 * rot(1,2)
   k2 = k_long * rot(2,1) + k_per1 * rot(2,2)

   ! get perpendicular unitary direction
   x_per = (/ -k_per1/sqrt(k_long**2+k_per1**2), k_long/sqrt(k_long**2+k_per1**2), 0.0_p_double /)

   ! calculate direction of perpendicular unitary direction in rotated space
   xper = (/ x_per(1)*rot(1,1)+x_per(2)*rot(1,2), x_per(1)*rot(2,1)+x_per(2)*rot(2,2), 0.0_p_double /)

   ! get jones vector for polarization
   jones_vec = this%get_jones_vector( k_long, k_per1, 0.0_p_double )

   ! calculate component for each field
   e1_fac = xper(1) * jones_vec(1)
   e2_fac = xper(2) * jones_vec(1)
   e3_fac = 1.0_p_double * jones_vec(2)
   b1_fac = - xper(1) * jones_vec(2)
   b2_fac = - xper(2) * jones_vec(2)
   b3_fac = 1.0_p_double * jones_vec(1)
   f_k = f_k * omega

   ! boost pulse
   if (this%if_boost) then
    ! save betas
    beta1 = this%beta_boost(1)
    beta2 = this%beta_boost(2)
    beta3 = this%beta_boost(3)
    beta = sqrt( beta1**2 + beta2**2 + beta3**2 )
    ! boost k-four-vector
    w_b = omega; k1_b = k1; k2_b = k2
    omega = this%gamma * (w_b - beta1*k1_b - beta2*k2_b)
    k1 = k1 + (this%gamma-1) * (beta1*k1_b + beta2*k2_b) * beta1 / beta**2 - this%gamma * w_b * beta1
    k2 = k2 + (this%gamma-1) * (beta1*k1_b + beta2*k2_b) * beta2 / beta**2 - this%gamma * w_b * beta2
    ! boost field factors
    e1_b = e1_fac; e2_b = e2_fac; e3_b = e3_fac; b1_b = b1_fac; b2_b = b2_fac; b3_b = b3_fac;
    e1_fac = this%gamma * ( e1_b + beta2*b3_b - beta3*b2_b ) - (this%gamma-1) * (e1_b*beta1 + e2_b*beta2 + e3_b*beta3) * beta1 / beta**2
    e2_fac = this%gamma * ( e2_b - beta1*b3_b + beta3*b1_b ) - (this%gamma-1) * (e1_b*beta1 + e2_b*beta2 + e3_b*beta3) * beta2 / beta**2
    e3_fac = this%gamma * ( e3_b + beta1*b2_b - beta2*b1_b ) - (this%gamma-1) * (e1_b*beta1 + e2_b*beta2 + e3_b*beta3) * beta3 / beta**2
    b1_fac = this%gamma * ( b1_b - beta2*e3_b + beta3*e2_b ) - (this%gamma-1) * (b1_b*beta1 + b2_b*beta2 + b3_b*beta3) * beta1 / beta**2
    b2_fac = this%gamma * ( b2_b + beta1*e3_b - beta3*e1_b ) - (this%gamma-1) * (b1_b*beta1 + b2_b*beta2 + b3_b*beta3) * beta2 / beta**2
    b3_fac = this%gamma * ( b3_b - beta1*e2_b + beta2*e1_b ) - (this%gamma-1) * (b1_b*beta1 + b2_b*beta2 + b3_b*beta3) * beta3 / beta**2
   endif

   ! start integration scheme
   f_k = f_k * exp( - im_unit * omega * this%lon_position + im_unit * this%phase0 )

   ! multiply components of e and b by f_k
   this%kvec_arr( n_wave_now,: ) = (/ k1, k2, 0.0_p_double /)
   this%evec_arr( n_wave_now,: ) = (/ e1_fac*f_k, e2_fac*f_k, e3_fac*f_k /)
   this%bvec_arr( n_wave_now,: ) = (/ b1_fac*f_k, b2_fac*f_k, b3_fac*f_k /)
   if ( this%print_waves ) fk_arr( n_wave_now ) = this%k_function( k_long, k_per1, 0.0_p_double )
   n_wave_now = n_wave_now + 1

  enddo
  ! exit if infinite wave (only one k_long)
  if ( this%lon_type == k_lon_plane ) exit
  enddo

  this%kvec_arr => this%kvec_arr( :n_wave_now-1,: )
  this%evec_arr => this%evec_arr( :n_wave_now-1,: )
  this%bvec_arr => this%bvec_arr( :n_wave_now-1,: )

 else if ( p_x_dim == 3 ) then

  num_waves = int( ( this%lon_int_range(2) - this%lon_int_range(1) ) / this%dk_waves(1) + 2 )
  num_waves = num_waves * int( ( this%per1_int_range(2) - this%per1_int_range(1) ) / this%dk_waves(2) + 2 )
  num_waves = num_waves * int( ( this%per2_int_range(2) - this%per2_int_range(1) ) / this%dk_waves(3) + 2 )

  allocate( this%kvec_arr( num_waves, 3 ) )
  allocate( this%evec_arr( num_waves, 3 ) )
  allocate( this%bvec_arr( num_waves, 3 ) )
  if ( this%print_waves ) allocate( fk_arr( num_waves ) )

  ! trignometry
  rot(1,:) = (/ cos(this%angles(1))*cos(this%angles(2)), -sin(this%angles(1)), -cos(this%angles(1))*sin(this%angles(2)) /)
  rot(2,:) = (/ sin(this%angles(1))*cos(this%angles(2)), cos(this%angles(1)), -sin(this%angles(1))*sin(this%angles(2)) /)
  rot(3,:) = (/ sin(this%angles(2)), 0.0_p_double, cos(this%angles(2)) /)

  ! calculate maximum amplitude for tolerance
  ! assuming the largest contribution is the principal wave ...
  max_amp = this%tolerance * abs( this%a0 * this%k_function( this%omega0, 0.0_p_double, 0.0_p_double ) * this%dk_waves(1) * this%dk_waves(2) * this%dk_waves(3) )
  if ( max_amp == 0 ) max_amp = this%tolerance

  ! iterate for every plane wave
  k_long = this%lon_int_range(1) - mod( this%lon_int_range(1), this%dk_waves(1) )
  do while ( k_long < this%lon_int_range(2) )
   ! calculate direction of each wave - longitudinal
   k_long = k_long + this%dk_waves(1)
  k_per1 = this%per1_int_range(1) - mod( this%per1_int_range(1), this%dk_waves(2) ) - 2*this%dk_waves(2)
  do while ( k_per1 < this%per1_int_range(2) )
   ! calculate direction of each wave - perpendicular 1
   k_per1 = k_per1 + this%dk_waves(2)
   if ( this%per_type == k_per_plane ) k_per1 = 0.0_p_double
  k_per2 = this%per2_int_range(1) - mod( this%per2_int_range(1), this%dk_waves(3) ) - 2*this%dk_waves(3)
  do while ( k_per2 < this%per2_int_range(2) )
   ! calculate direction of each wave - perpendicular 2
   k_per2 = k_per2 + this%dk_waves(3)
   if ( this%per_type == k_per_plane ) k_per2 = 0.0_p_double
  ! --------------------------------------------------

   ! calculate intensity of each wave
   f_k = this%a0 * this%k_function( k_long, k_per1, k_per2 ) * this%dk_waves(1) * this%dk_waves(2) * this%dk_waves(3)

   ! tolerance for speep up
   if ( abs(f_k) < max_amp ) cycle

   ! calculate omega corresponding to this wave
   omega = sqrt( k_long**2 + k_per1**2 + k_per2**2 )
   !omega = 2/dt * asin( dt * sqrt( sin(k1*ldx(1)/2)**2/ldx(1)**2 + sin(k2*ldx(2)/2)**2/ldx(2)**2 + sin(k3*ldx(3)/2)**2/ldx(3)**2 ) )

   ! take care of infinite pulse
   if ( this%lon_type == k_lon_plane ) then
     omega = this%omega0
    k_long = sqrt( this%omega0**2 - k_per1**2 - k_per2**2 )
   endif
   if (mpi_node()==0) print *,"   - k_long", k_long, "k_per1", k_per1, "k_per2", k_per2, "omega", omega

   ! calculate direction of wave in rotated space
   k1 = k_long * rot(1,1) + k_per1 * rot(1,2) + k_per2 * rot(1,3)
   k2 = k_long * rot(2,1) + k_per1 * rot(2,2) + k_per2 * rot(2,3)
   k3 = k_long * rot(3,1) + k_per1 * rot(3,2) + k_per2 * rot(3,3)

   ! get perpendicular unitary directions
   x_per1 = (/ -k_per1/sqrt(k_long**2+k_per1**2), k_long/sqrt(k_long**2+k_per1**2), 0.0_p_double /)
   x_per2 = (/ -k_long*k_per2/sqrt(k_long**2+k_per1**2)/omega, -k_per1*k_per2/sqrt(k_long**2+k_per1**2)/omega, sqrt(k_long**2+k_per1**2)/omega /)

   ! calculate direction of perpendicular unitary directions in rotated space
   xper1 = (/ x_per1(1)*rot(1,1)+x_per1(2)*rot(1,2)+x_per1(3)*rot(1,3), x_per1(1)*rot(2,1)+x_per1(2)*rot(2,2)+x_per1(3)*rot(2,3), x_per1(1)*rot(3,1)+x_per1(2)*rot(3,2)+x_per1(3)*rot(3,3) /)
   xper2 = (/ x_per2(1)*rot(1,1)+x_per2(2)*rot(1,2)+x_per2(3)*rot(1,3), x_per2(1)*rot(2,1)+x_per2(2)*rot(2,2)+x_per2(3)*rot(2,3), x_per2(1)*rot(3,1)+x_per2(2)*rot(3,2)+x_per2(3)*rot(3,3) /)

   ! get jones vector for polarization
   jones_vec = this%get_jones_vector( k_long, k_per1, k_per2 )

   ! calculate component for each field
   e1_fac = xper1(1) * jones_vec(1) + xper2(1) * jones_vec(2)
   e2_fac = xper1(2) * jones_vec(1) + xper2(2) * jones_vec(2)
   e3_fac = xper1(3) * jones_vec(1) + xper2(3) * jones_vec(2)
   b1_fac = xper2(1) * jones_vec(1) - xper1(1) * jones_vec(2)
   b2_fac = xper2(2) * jones_vec(1) - xper1(2) * jones_vec(2)
   b3_fac = xper2(3) * jones_vec(1) - xper1(3) * jones_vec(2)
   f_k = f_k * omega

   ! boost pulse
   if (this%if_boost) then
    ! save betas
    beta1 = this%beta_boost(1)
    beta2 = this%beta_boost(2)
    beta3 = this%beta_boost(3)
    beta = sqrt( beta1**2 + beta2**2 + beta3**2 )
    ! boost k-four-vector
    w_b = omega; k1_b = k1; k2_b = k2; k3_b = k3
    omega = this%gamma * (w_b - beta1*k1_b - beta2*k2_b - beta2*k3_b)
    k1 = k1 + (this%gamma-1) * (beta1*k1_b + beta2*k2_b + beta3*k3_b) * beta1 / beta**2 - this%gamma * w_b * beta1
    k2 = k2 + (this%gamma-1) * (beta1*k1_b + beta2*k2_b + beta3*k3_b) * beta2 / beta**2 - this%gamma * w_b * beta2
    k3 = k3 + (this%gamma-1) * (beta1*k1_b + beta2*k2_b + beta3*k3_b) * beta3 / beta**2 - this%gamma * w_b * beta3
    ! boost field factors
    e1_b = e1_fac; e2_b = e2_fac; e3_b = e3_fac; b1_b = b1_fac; b2_b = b2_fac; b3_b = b3_fac;
    e1_fac = this%gamma * ( e1_b + beta2*b3_b - beta3*b2_b ) - (this%gamma-1) * (e1_b*beta1 + e2_b*beta2 + e3_b*beta3) * beta1 / beta**2
    e2_fac = this%gamma * ( e2_b - beta1*b3_b + beta3*b1_b ) - (this%gamma-1) * (e1_b*beta1 + e2_b*beta2 + e3_b*beta3) * beta2 / beta**2
    e3_fac = this%gamma * ( e3_b + beta1*b2_b - beta2*b1_b ) - (this%gamma-1) * (e1_b*beta1 + e2_b*beta2 + e3_b*beta3) * beta3 / beta**2
    b1_fac = this%gamma * ( b1_b - beta2*e3_b + beta3*e2_b ) - (this%gamma-1) * (b1_b*beta1 + b2_b*beta2 + b3_b*beta3) * beta1 / beta**2
    b2_fac = this%gamma * ( b2_b + beta1*e3_b - beta3*e1_b ) - (this%gamma-1) * (b1_b*beta1 + b2_b*beta2 + b3_b*beta3) * beta2 / beta**2
    b3_fac = this%gamma * ( b3_b - beta1*e2_b + beta2*e1_b ) - (this%gamma-1) * (b1_b*beta1 + b2_b*beta2 + b3_b*beta3) * beta3 / beta**2
   endif

   ! start integration scheme
   f_k = f_k * exp( - im_unit * omega * this%lon_position + im_unit * this%phase0 )

   ! multiply components of e and b by f_k
   this%kvec_arr( n_wave_now,: ) = (/ k1, k2, k3 /)
   this%evec_arr( n_wave_now,: ) = (/ e1_fac*f_k, e2_fac*f_k, e3_fac*f_k /)
   this%bvec_arr( n_wave_now,: ) = (/ b1_fac*f_k, b2_fac*f_k, b3_fac*f_k /)
   if ( this%print_waves ) fk_arr( n_wave_now ) = this%k_function( k_long, k_per1, k_per2 )
   n_wave_now = n_wave_now + 1

  enddo
  enddo
  ! exit if infinite wave (only one k_long)
  if ( this%lon_type == k_lon_plane ) exit
  enddo

  this%kvec_arr => this%kvec_arr( :n_wave_now-1,: )
  this%evec_arr => this%evec_arr( :n_wave_now-1,: )
  this%bvec_arr => this%bvec_arr( :n_wave_now-1,: )

 endif

 if (mpi_node()==0) print *,"   - will be injecting", n_wave_now-1, "waves"
 if (mpi_node()==0) print *,"   - lon_int_range", this%lon_int_range
 if (mpi_node()==0) print *,"   - per1_int_range", this%per1_int_range
 if (mpi_node()==0) print *,"   - dk_waves", this%dk_waves
 if (mpi_node()==0) print *,"   - tolerance", this%tolerance
 if (mpi_node()==0) print *,"   - max_amp", max_amp


 if ( this%print_waves .and. mpi_node() == 0 ) then
  if (mpi_node()==0) print *,"   - printing waves to file"
  open( 100, file = "waves.dat", status = "unknown" )
  write( 100, * ) "k1 k2 k3 e1r e1i e2r e2i e3r e3i b1r b1i b2r b2i b3r b3i fkr fki"
  do ik1 = 1, size( this%kvec_arr, 1 )
   write( 100, * ) this%kvec_arr( ik1, 1 ), this%kvec_arr( ik1, 2 ), this%kvec_arr( ik1, 3 ), real(this%evec_arr( ik1, 1 )), dimag(this%evec_arr( ik1, 1 )), real(this%evec_arr( ik1, 2 )), dimag(this%evec_arr( ik1, 2 )), real(this%evec_arr( ik1, 3 )), dimag(this%evec_arr( ik1, 3 )), real(this%bvec_arr( ik1, 1 )), dimag(this%bvec_arr( ik1, 1 )), real(this%bvec_arr( ik1, 2 )), dimag(this%bvec_arr( ik1, 2 )), real(this%bvec_arr( ik1, 3 )), dimag(this%bvec_arr( ik1, 3 )), real(fk_arr( ik1 )), dimag(fk_arr( ik1 ))







  enddo
  close( 100 )
  deallocate( fk_arr )
 endif

end subroutine init_zpulse_planewaves
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
! Initialize zpulse_planewaves
subroutine init_zpulse_planewaves_incoherent( this, restart, t , dt, b, g_space, nx_p_min, no_co, interpolation )

 implicit none

 class( t_zpulse_planewaves ), intent(inout) :: this
 logical, intent(in) :: restart
 real( p_double ), intent(in) :: t, dt
 type( t_vdf ) , intent(in) :: b
 type( t_space ), intent(in) :: g_space
 integer, intent(in), dimension(:) :: nx_p_min
 class( t_node_conf ), intent(in) :: no_co
 integer, intent(in) :: interpolation

 logical :: if_launch
 integer :: ik1, ik2, i1, i2, num_waves, n_wave_now
 real(p_double), dimension(3,3) :: rot
 real(p_double) :: k_long, k_per1, k1, k2, omega, max_amp
 complex(p_double) :: f_k
 real(p_double), dimension(3) :: x_per, xper, x_per1, xper1
 complex(p_double), dimension(2) :: jones_vec
 real(p_double) :: beta, beta1, beta2, beta3, w_b, k1_b, k2_b
 complex(p_double) :: e1_b, e2_b, e3_b, b1_b, b2_b, b3_b
 complex(p_double) :: e1_fac, e2_fac, e3_fac, b1_fac, b2_fac, b3_fac

 complex(p_double), dimension (:,:), pointer :: f_k_mat

 complex(p_double), dimension (:,:), pointer :: filter
 integer :: filter_min, filter_max

 ! create basic filter
 filter_min = -10
 filter_max = 10
 allocate( filter(filter_min:filter_max, filter_min:filter_max) )
 f_k = 0
 do i1 = filter_min, filter_max
 do i2 = filter_min, filter_max
  call random_number( k1 )
  call random_number( k2 )
  k1 = 1.0_p_double + sqrt(-2.0*log(k1))*cos(2*pi*k2)
  call random_number( k2 )
  filter( i1, i2 ) = k1 * exp( im_unit *2*pi*k2 )
  f_k = f_k + k1 * exp( im_unit *2*pi*k2 )
 enddo
 enddo
 filter = filter / f_k

 if (mpi_node()==0) print *,"   - calculating plane waves incoherent"

 allocate( f_k_mat( int( ( this%lon_int_range(2) - this%lon_int_range(1) ) / this%dk_waves(1) ), int( ( this%per1_int_range(2) - this%per1_int_range(1) ) / this%dk_waves(2) ) ) )

 max_amp = 0

 ! iterate for every plane wave
 do ik1 = 1, int( ( this%lon_int_range(2) - this%lon_int_range(1) ) / this%dk_waves(1) )
  ! calculate direction of each wave - longitudinal
  k_long = this%lon_int_range(1) - mod( this%lon_int_range(1), this%dk_waves(1) ) + ( ik1-1 ) * this%dk_waves(1)
 do ik2 = 1, int( ( this%per1_int_range(2) - this%per1_int_range(1) ) / this%dk_waves(2) )
  ! calculate direction of each wave - perpendicular
  k_per1 = this%per1_int_range(1) - mod( this%per1_int_range(1), this%dk_waves(2) ) + ( ik2-1 ) * this%dk_waves(2)
 ! --------------------------------------------------
  ! calculate intensity of each wave
  f_k_mat( ik1, ik2 ) = this%k_function( k_long, k_per1, 0.0_p_double )
  if ( abs(f_k_mat( ik1, ik2 )) > max_amp ) then
   max_amp = abs(f_k_mat( ik1, ik2 ))
  endif
 enddo
 enddo

 ! trignometry
 rot(1,:) = (/ cos( this%angles(1) ), -sin( this%angles(1) ), 0.0_p_double /)
 rot(2,:) = (/ sin( this%angles(1) ), cos( this%angles(1) ), 0.0_p_double /)

 ! calculate maximum amplitude for tolerance
 max_amp = this%tolerance * max_amp

 if ( this%lon_type == k_lon_plane ) then
  num_waves = 1
 else
  num_waves = int( ( this%lon_int_range(2) - this%lon_int_range(1) ) / this%dk_waves(1) + 10 )
 endif
 if ( this%per_type /= k_per_plane ) then
  num_waves = num_waves * int( ( this%per1_int_range(2) - this%per1_int_range(1) ) / this%dk_waves(2) + 10 )
 endif
 allocate( this%kvec_arr( num_waves, 3 ) )
 allocate( this%evec_arr( num_waves, 3 ) )
 allocate( this%bvec_arr( num_waves, 3 ) )

 ! iterate for every plane wave
 n_wave_now = 0
 do ik1 = 1, int( ( this%lon_int_range(2) - this%lon_int_range(1) ) / this%dk_waves(1) )
  ! calculate direction of each wave - longitudinal
  k_long = this%lon_int_range(1) - mod( this%lon_int_range(1), this%dk_waves(1) ) + ( ik1-1 ) * this%dk_waves(1)
 do ik2 = 1, int( ( this%per1_int_range(2) - this%per1_int_range(1) ) / this%dk_waves(2) )
  ! calculate direction of each wave - perpendicular
  k_per1 = this%per1_int_range(1) - mod( this%per1_int_range(1), this%dk_waves(2) ) + ( ik2-1 ) * this%dk_waves(2)
 ! --------------------------------------------------

  ! apply filter
  f_k = 0
  do i1 = filter_min, filter_max
  do i2 = filter_min, filter_max
   if ( ik1 - i1 < 1 ) cycle
   if ( ik2 - i2 < 1 ) cycle
   if ( ik1 - i1 > int( ( this%lon_int_range(2) - this%lon_int_range(1) ) / this%dk_waves(1) ) ) cycle
   if ( ik2 - i2 > int( ( this%per1_int_range(2) - this%per1_int_range(1) ) / this%dk_waves(2) ) ) cycle
   f_k = f_k + filter( i1, i2 ) * f_k_mat( ik1 - i1, ik2 - i2 )
  enddo
  enddo

  if ( abs(f_k) < max_amp ) cycle

  ! calculate omega corresponding to this wave
  omega = sqrt( k_long**2 + k_per1**2 )

  ! calculate intensity of each wave
  f_k = this%a0 * f_k * this%dk_waves(1) * this%dk_waves(2)

  ! tolerance for speep up
  if ( abs(f_k) < max_amp ) cycle

  ! calculate direction of wave in rotated space
  k1 = k_long * rot(1,1) + k_per1 * rot(1,2)
  k2 = k_long * rot(2,1) + k_per1 * rot(2,2)

  ! get perpendicular unitary direction
  x_per = (/ -k_per1/sqrt(k_long**2+k_per1**2), k_long/sqrt(k_long**2+k_per1**2), 0.0_p_double /)

  ! calculate direction of perpendicular unitary direction in rotated space
  xper = (/ x_per(1)*rot(1,1)+x_per(2)*rot(1,2), x_per(1)*rot(2,1)+x_per(2)*rot(2,2), 0.0_p_double /)

  ! get jones vector for polarization
  jones_vec = this%get_jones_vector( k_long, k_per1, 0.0_p_double )

  ! calculate component for each field
  e1_fac = xper(1) * jones_vec(1)
  e2_fac = xper(2) * jones_vec(1)
  e3_fac = 1.0_p_double * jones_vec(2)
  b1_fac = - xper(1) * jones_vec(2)
  b2_fac = - xper(2) * jones_vec(2)
  b3_fac = 1.0_p_double * jones_vec(1)
  f_k = f_k * omega

  ! boost pulse
  if (this%if_boost) then
   ! save betas
   beta1 = this%beta_boost(1)
   beta2 = this%beta_boost(2)
   beta3 = this%beta_boost(3)
   beta = sqrt( beta1**2 + beta2**2 + beta3**2 )
   ! boost k-four-vector
   w_b = omega; k1_b = k1; k2_b = k2
   omega = this%gamma * (w_b - beta1*k1_b - beta2*k2_b)
   k1 = k1 + (this%gamma-1) * (beta1*k1_b + beta2*k2_b) * beta1 / beta**2 - this%gamma * w_b * beta1
   k2 = k2 + (this%gamma-1) * (beta1*k1_b + beta2*k2_b) * beta2 / beta**2 - this%gamma * w_b * beta2
   ! boost field factors
   e1_b = e1_fac; e2_b = e2_fac; e3_b = e3_fac; b1_b = b1_fac; b2_b = b2_fac; b3_b = b3_fac;
   e1_fac = this%gamma * ( e1_b + beta2*b3_b - beta3*b2_b ) - (this%gamma-1) * (e1_b*beta1 + e2_b*beta2 + e3_b*beta3) * beta1 / beta**2
   e2_fac = this%gamma * ( e2_b - beta1*b3_b + beta3*b1_b ) - (this%gamma-1) * (e1_b*beta1 + e2_b*beta2 + e3_b*beta3) * beta2 / beta**2
   e3_fac = this%gamma * ( e3_b + beta1*b2_b - beta2*b1_b ) - (this%gamma-1) * (e1_b*beta1 + e2_b*beta2 + e3_b*beta3) * beta3 / beta**2
   b1_fac = this%gamma * ( b1_b - beta2*e3_b + beta3*e2_b ) - (this%gamma-1) * (b1_b*beta1 + b2_b*beta2 + b3_b*beta3) * beta1 / beta**2
   b2_fac = this%gamma * ( b2_b + beta1*e3_b - beta3*e1_b ) - (this%gamma-1) * (b1_b*beta1 + b2_b*beta2 + b3_b*beta3) * beta2 / beta**2
   b3_fac = this%gamma * ( b3_b - beta1*e2_b + beta2*e1_b ) - (this%gamma-1) * (b1_b*beta1 + b2_b*beta2 + b3_b*beta3) * beta3 / beta**2
  endif

  ! start integration scheme
  f_k = f_k * exp( - im_unit * omega * this%lon_position + im_unit * this%phase0 )

  ! multiply components of e and b by f_k
  this%kvec_arr( n_wave_now,: ) = (/ k1, k2, 0.0_p_double /)
  this%evec_arr( n_wave_now,: ) = (/ e1_fac*f_k, e2_fac*f_k, e3_fac*f_k /)
  this%bvec_arr( n_wave_now,: ) = (/ b1_fac*f_k, b2_fac*f_k, b3_fac*f_k /)
  n_wave_now = n_wave_now + 1

 enddo
 enddo

 this%kvec_arr => this%kvec_arr( :n_wave_now-1,: )
 this%evec_arr => this%evec_arr( :n_wave_now-1,: )
 this%bvec_arr => this%bvec_arr( :n_wave_now-1,: )

 if (mpi_node()==0) print *,"   - will be injecting", n_wave_now-1, "waves"
 if (mpi_node()==0) print *,"   - lon_int_range", this%lon_int_range
 if (mpi_node()==0) print *,"   - per1_int_range", this%per1_int_range
 if (mpi_node()==0) print *,"   - dk_waves", this%dk_waves
 if (mpi_node()==0) print *,"   - tolerance", this%tolerance
 if (mpi_node()==0) print *,"   - max_amp", max_amp


end subroutine init_zpulse_planewaves_incoherent
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
! Reading the input file
subroutine read_input_zpulse_planewaves( this, input_file, g_space, bnd_con, periodic, grid, sim_options )

 use m_input_file

 implicit none

 class( t_zpulse_planewaves ), intent(inout) :: this
 class( t_input_file ), intent(inout) :: input_file
 type( t_space ), intent(in) :: g_space
 class (t_emf_bound), intent(in) :: bnd_con
 logical, dimension(:), intent(in) :: periodic
 class( t_grid ), intent(in) :: grid
 type( t_options ), intent(in) :: sim_options

 ! ----------------------------------------------------------------------
 ! local variables
 real(p_double) :: a0
 real(p_double) :: omega0
 integer :: pol_type
 real(p_double) :: pol
 real(p_double) :: pol_elipse
 real(p_double), dimension(2) :: angles
 real(p_double), dimension(3) :: focus
 real(p_double), dimension(3) :: start
 real(p_double) :: phase
 real(p_double) :: tolerance
 logical, dimension(3,2) :: wall_injection
 real(p_double), dimension(3) :: beta_boost
 real(p_double) :: t_wall_start, t_wall_stop
 logical :: initial_injection

 real(p_double), dimension(3) :: periodicity

 character(len=8) :: preset

 logical :: incoherent
 logical :: print_waves

 character(len=16) :: pulse_type
 real(p_double) :: vfocus
 real(p_double) :: afocus
 character(len=p_max_expr_len):: f_k_amp_math_func, f_k_pha_math_func
 character(len=p_max_expr_len), dimension(2):: jones_vec_real, jones_vec_imag

 character(len=16) :: lon_type
 real(p_double) :: lon_position
 real(p_double) :: lon_duration, lon_fwhm, lon_range ! for gaussian profile
 real(p_double) :: lon_rise, lon_flat, lon_fall ! for sin2 and polynomial profiles
 real(p_double), dimension(2) :: lon_int_range

 character(len=16) :: per_type
 real(p_double), dimension(2) :: per_w0, per_fwhm
 integer, dimension(2) :: per_tem_mode
 real(p_double), dimension(2) :: per1_int_range, per2_int_range

 logical :: if_launch
 real(p_double) :: launch_time

 namelist /nl_zpulse_planewaves/ &
  a0, omega0, pol_type, pol, pol_elipse, jones_vec_real, jones_vec_imag, angles, focus, start, phase, preset, beta_boost, periodicity, &
  lon_type, lon_position, lon_duration, lon_fwhm, lon_range, lon_rise, lon_flat, lon_fall, lon_int_range, &
  per_type, per_w0, per_fwhm, per_tem_mode, per1_int_range, per2_int_range, &
  pulse_type, vfocus, afocus, f_k_amp_math_func, f_k_pha_math_func, tolerance, incoherent, &
  wall_injection, t_wall_start, t_wall_stop, initial_injection, if_launch, print_waves

 ! ----------------------------------------------------------------------
 ! auxiliar variables
 integer :: ierr, i1
 real(p_double) :: sig ! auxiliar to define integration ranges
 real(p_double) :: dper, dlon

 ! ----------------------------------------------------------------------
 ! defining default values
 if_launch = .true.
 launch_time = 0.0_p_double

 a0 = 1.0_p_double
 omega0 = 6.28_p_double
 pol_type = 0
 pol = 0.0_p_double
 pol_elipse = 1.0_p_double
 angles = (/ 0.0_p_double, 0.0_p_double /)
 focus = (/ 0.0_p_double, 0.0_p_double, 0.0_p_double /)
 start = (/ -huge(1.0_p_double), -huge(1.0_p_double), -huge(1.0_p_double) /)
 phase = 0.0_p_double
 beta_boost = (/ 0.0_p_double, 0.0_p_double, 0.0_p_double /)

 periodicity = (/ 0.0_p_double, 0.0_p_double, 0.0_p_double /)

 incoherent = .false.
 print_waves = .false.

 preset = "0"
 jones_vec_real(1) = "NO_FUNCTION_SUPPLIED!"
 jones_vec_real(2) = "NO_FUNCTION_SUPPLIED!"
 jones_vec_imag(1) = "NO_FUNCTION_SUPPLIED!"
 jones_vec_imag(2) = "NO_FUNCTION_SUPPLIED!"

 pulse_type = "time"
 vfocus = 0.0_p_double
 afocus = 0.0_p_double
 f_k_amp_math_func = "NO_FUNCTION_SUPPLIED!"
 f_k_pha_math_func = "NO_FUNCTION_SUPPLIED!"
 tolerance = 0.0_p_double

 wall_injection(:,:) = .false.
 t_wall_start = 0.0_p_double
 t_wall_stop = 0.0_p_double
 initial_injection = .false.

 lon_type = "sin2"
 lon_position = -huge(1.0_p_double)
 lon_duration = 10.0_p_double
 lon_range = -huge(1.0_p_double)
 lon_fwhm = 0.0_p_double
 lon_rise = pi
 lon_flat = 0.0_p_double
 lon_fall = pi
 lon_int_range = (/ 0.0_p_double, 0.0_p_double /)

 per_type = "gaussian"
 per_w0 = 2.0_p_double
 per_fwhm = 0.0_p_double
 per_tem_mode = 0
 per1_int_range = (/ 0.0_p_double, 0.0_p_double /)
 per2_int_range = (/ 0.0_p_double, 0.0_p_double /)

 ! ----------------------------------------------------------------------
 ! Get namelist text from input file

 if (mpi_node()==0) print *," - Reading zpulse_planewaves configuration..."

 call get_namelist( input_file, "nl_zpulse_planewaves", ierr )

 if (ierr /= 0) then
 if ( mpi_node() == 0 ) then
  if ( ierr < 0 ) then
  print *, "Error reading zpulse_planewaves parameters"
  else
  print *, "Error: zpulse_planewaves parameters missing"
  endif
  print *, "aborting..."
 endif
 stop
 endif

 ! ----------------------------------------------------------------------
 ! read input file
 read (input_file%nml_text, nml = nl_zpulse_planewaves, iostat = ierr)
 if (ierr /= 0) then
 if ( mpi_node() == 0 ) then
  write(0,*) ""
  write(0,*) "   Error reading zpulse_planewaves parameters"
  write(0,*) "   aborting..."
 endif
 stop
 endif

 ! ----------------------------------------------------------------------
 ! save fundamental parameters
 this%a0 = a0
 this%omega0 = omega0
 this%focus = focus
 this%start = start
 this%phase0 = real( phase * pi_180, p_double )
 this%tolerance = tolerance

 this%incoherent = incoherent
 this%print_waves = print_waves

 if ( pol_type < -1 .or. pol_type > 4 ) then
  if ( mpi_node() == 0 ) then
   write(0,*) ""
   write(0,*) "   Error in zpulse parameters"
   write(0,*) "   pol must be in the range [-1,0,+1,+2,+3]"
   write(0,*) "    -  0 : linear"
   write(0,*) "    - -1 : circular left"
   write(0,*) "    - +1 : circular right"
   write(0,*) "    - +2 : radial"
   write(0,*) "    - +3 : azimuthal"
   write(0,*) "    - +4 : custom - use jones_vec(1:2) to specify pol"
   write(0,*) "   aborting..."
  endif
  stop
 endif
 this%pol_type = pol_type

 if ( pol_elipse < 0 .or. pol_elipse > 1 ) then
  if ( mpi_node() == 0 ) then
   write(0,*) ""
   write(0,*) "   Error in zpulse parameters"
   write(0,*) "   pol_elipse must be in the range [0,1]"
   write(0,*) "   aborting..."
  endif
  stop
 endif
 this%pol_elipse = pol_elipse

 if ( this%pol_type == 4 ) then
  call setup(this%jones_vec_real(1), trim(jones_vec_real(1)), (/"k1", "k2", "k3"/), ierr)
  if (ierr /= 0) then
   if ( mpi_node() == 0 ) then
    write(0,*) ''
    write(0,*) '   Error in zpulse_planewaves parameters'
    write(0,*) '   Supplied jones_vec_real(1) failed to compile'
    write(0,*) '   aborting...'
   endif
     stop
  endif
  call setup(this%jones_vec_real(2), trim(jones_vec_real(2)), (/"k1", "k2", "k3"/), ierr)
  if (ierr /= 0) then
   if ( mpi_node() == 0 ) then
    write(0,*) ''
    write(0,*) '   Error in zpulse_planewaves parameters'
    write(0,*) '   Supplied jones_vec_real(2) failed to compile'
    write(0,*) '   aborting...'
   endif
     stop
  endif
  call setup(this%jones_vec_imag(1), trim(jones_vec_imag(1)), (/"k1", "k2", "k3"/), ierr)
  if (ierr /= 0) then
   if ( mpi_node() == 0 ) then
    write(0,*) ''
    write(0,*) '   Error in zpulse_planewaves parameters'
    write(0,*) '   Supplied jones_vec_imag(1) failed to compile'
    write(0,*) '   aborting...'
   endif
     stop
  endif
  call setup(this%jones_vec_imag(2), trim(jones_vec_imag(2)), (/"k1", "k2", "k3"/), ierr)
  if (ierr /= 0) then
   if ( mpi_node() == 0 ) then
    write(0,*) ''
    write(0,*) '   Error in zpulse_planewaves parameters'
    write(0,*) '   Supplied jones_vec_imag(2) failed to compile'
    write(0,*) '   aborting...'
   endif
     stop
  endif
 endif


 select case ( trim(preset) )
 case ("+x1x2"); angles(1) = 0; angles(2) = 0; pol = 0
 case ("-x1x2"); angles(1) = 180; angles(2) = 0; pol = 180
 case ("+x1x3"); angles(1) = 0; angles(2) = 0; pol = 90
 case ("-x1x3"); angles(1) = 180; angles(2) = 0; pol = 90
 case ("+x2x1"); angles(1) = 90; angles(2) = 0; pol = 180
 case ("-x2x1"); angles(1) = 270; angles(2) = 0; pol = 0
 case ("+x2x3"); angles(1) = 90; angles(2) = 0; pol = 90
 case ("-x2x3"); angles(1) = 270; angles(2) = 0; pol = 90
 case ("+x3x1"); angles(1) = 0; angles(2) = 90; pol = 270
 case ("-x3x1"); angles(1) = 0; angles(2) = -90; pol = 90
 case ("+x3x2"); angles(1) = 0; angles(2) = 90; pol = 0
 case ("-x3x2"); angles(1) = 0; angles(2) = -90; pol = 0
 case ("0"); continue
 case default
  if ( mpi_node() == 0 ) then
   write(0,*) ''
   write(0,*) '   Error in zpulse_planewaves parameters'
   write(0,*) '   preset input does not exist '
   write(0,*) '   aborting...'
  endif
  stop
 end select

 this%angles(1) = real( angles(1) * pi_180, p_double )
 this%angles(2) = real( angles(2) * pi_180, p_double )
 this%pol = real( pol * pi_180, p_double )

 this%normalize = .false.
 this%normalization = 1.0_p_double

 ! ----------------------------------------------------------------------
 do i1 = 1, p_x_dim
  if ( periodicity(i1) == 0.0_p_double ) then
   this%dk_waves(i1) = 2*pi / ( xmax( g_space, i1 ) - xmin( g_space, i1 ) )
  else
   this%dk_waves(i1) = 2*pi / periodicity(i1)
  endif
 enddo

 ! ----------------------------------------------------------------------
 ! save pulse type
 select case ( trim(pulse_type) )
 case ("time")
  this%pulse_type= type_time

 case ("space")
  this%pulse_type= type_space

 case ("vfocus_const")
  this%pulse_type= type_vfocus_const
  if ( vfocus == 1.0_p_double ) then
   if ( mpi_node() == 0 ) then
    write(0,*) ''
    write(0,*) '   Error in zpulse_planewaves parameters'
    write(0,*) '   vfocus cannot be exactly 1 '
    write(0,*) '   aborting...'
   endif
   stop
  endif
  this%vfocus = vfocus

 case ("vfocus_accel")
  this%pulse_type= type_vfocus_accel
  if ( vfocus == 1.0_p_double ) then
   if ( mpi_node() == 0 ) then
    write(0,*) ''
    write(0,*) '   Error in zpulse_planewaves parameters'
    write(0,*) '   vfocus cannot be exactly 1 '
    write(0,*) '   aborting...'
   endif
   stop
  endif
  if ( afocus == 0.0_p_double ) then
   if ( mpi_node() == 0 ) then
    write(0,*) ''
    write(0,*) '   Error in zpulse_planewaves parameters'
    write(0,*) '   afocus cannot be exactly 0 '
    write(0,*) '   aborting...'
   endif
   stop
  endif
  this%vfocus = vfocus
  this%afocus = afocus

 case ("custom_full")
  this%pulse_type = type_custom_full

 case ("custom_lon")
  this%pulse_type = type_custom_lon

 case ("custom_per")
  this%pulse_type = type_custom_per

 case default
  if ( mpi_node() == 0 ) then
   write(0,*) ''
   write(0,*) '   Error in zpulse_planewaves parameters'
   write(0,*) '   pulse_type must be "time", "space", "vfocus_const", '
   write(0,*) '                      "custom_full", "custom_lon" or "custom_per" '
   write(0,*) '   aborting...'
  endif
  stop
 end select

 ! save function parsers for custom profiles
 if ( this%pulse_type == type_custom_full .or. this%pulse_type == type_custom_lon .or. this%pulse_type == type_custom_per ) then
  call setup(this%f_k_amp_math_func, trim(f_k_amp_math_func), (/"k1", "k2", "k3"/), ierr)
  if (ierr /= 0) then
   if ( mpi_node() == 0 ) then
    write(0,*) ''
    write(0,*) '   Error in zpulse_planewaves parameters'
    write(0,*) '   Supplied f_k_amp_math_func failed to compile'
    write(0,*) '   aborting...'
   endif
     stop
  endif
  call setup(this%f_k_pha_math_func, trim(f_k_pha_math_func), (/"k1", "k2", "k3"/), ierr)
  if (ierr /= 0) then
   if ( mpi_node() == 0 ) then
    write(0,*) ''
    write(0,*) '   Error in zpulse_planewaves parameters'
    write(0,*) '   Supplied f_k_pha_math_func failed to compile'
    write(0,*) '   aborting...'
   endif
     stop
  endif
 endif

 ! ----------------------------------------------------------------------
 ! save longitudinal profile parameters
 select case ( trim(lon_type) )
 case ("plane")
  this%lon_type = k_lon_plane
  this%lon_int_range(1) = this%omega0
  this%lon_int_range(2) = this%omega0
  this%dk_waves(1) = 1

 case ("gaussian")
  this%lon_type = k_lon_gaussian

  if ( lon_fwhm > 0.0_p_double ) then
   this%lon_duration = lon_fwhm/sqrt( 2.0_p_double * log(2.0_p_double))
  else
   this%lon_duration = lon_duration
  endif
  this%lon_range = lon_range

  if ( lon_range == -huge(1.0_p_double) ) then
   if ( mpi_node() == 0 ) then
   write(0,*) ''
   write(0,*) '   Error in zpulse_planewaves parameters'
   write(0,*) '   When using "gaussian" longitudinal envelopes, '
   write(0,*) '   lon_range must also be defined'
   write(0,*) '   aborting...'
   endif
   stop
  endif

  if ( lon_int_range(1) == 0 .and. lon_int_range(2) == 0 ) then
   sig = 2.0_p_double*sqrt( 2.0_p_double * log(2.0_p_double)) / this%lon_duration
   this%lon_int_range = (/ this%omega0 - n_sig_gaussian*sig, this%omega0 + n_sig_gaussian*sig /)
  else
   this%lon_int_range(1) = lon_int_range(1)
   this%lon_int_range(2) = lon_int_range(2)
  endif

 case ("sin2")
  this%lon_type = k_lon_sin2

  if ( lon_fwhm > 0.0_p_double ) then
   this%lon_rise = lon_fwhm/2.0_p_double
   this%lon_flat = 0.0_p_double
   this%lon_fall = lon_fwhm/2.0_p_double
  else
   this%lon_rise = lon_rise
   this%lon_flat = lon_flat
   this%lon_fall = lon_fall
  endif

  if ( lon_int_range(1) == 0 .and. lon_int_range(2) == 0 ) then
   sig = 2.0_p_double*pi / (this%lon_rise + 0.1*this%lon_flat + this%lon_fall)
   this%lon_int_range = (/ this%omega0 - n_sig_sin2*sig, this%omega0 + n_sig_sin2*sig /)
  else
   this%lon_int_range(1) = lon_int_range(1)
   this%lon_int_range(2) = lon_int_range(2)
  endif

 case ("polynomial")
  this%lon_type = k_lon_polynomial

  if ( lon_fwhm > 0.0_p_double ) then
   this%lon_rise = lon_fwhm/2.0_p_double
   this%lon_flat = 0.0_p_double
   this%lon_fall = lon_fwhm/2.0_p_double
  else
   this%lon_rise = lon_rise
   this%lon_flat = lon_flat
   this%lon_fall = lon_fall
  endif

  if ( lon_int_range(1) == 0 .and. lon_int_range(2) == 0 ) then
   sig = 2.0_p_double*pi / (this%lon_rise + 0.1*this%lon_flat + this%lon_fall)
   this%lon_int_range = (/ this%omega0 - n_sig_polynomial*sig, this%omega0 + n_sig_polynomial*sig /)
  else
   this%lon_int_range(1) = lon_int_range(1)
   this%lon_int_range(2) = lon_int_range(2)
  endif

 case default
  if ( mpi_node() == 0 ) then
   write(0,*) ''
   write(0,*) '   Error in zpulse_planewaves parameters'
   write(0,*) '   lon_type must be "gaussian" or "sin2" '
   write(0,*) '   no other type has been implemented yet '
   write(0,*) '   aborting...'
  endif
  stop
 end select

 if ( lon_position /= -huge(1.0_p_double) ) then
  this%lon_position = lon_position
 else if( this%start(1) /= -huge(1.0_p_double) .or. this%start(2) /= -huge(1.0_p_double) .or. this%start(3) /= -huge(1.0_p_double) ) then
  this%lon_position = - ( this%focus(1) - this%start(1) ) * cos( this%angles(1) ) * cos( this%angles(2) ) - ( this%focus(2) - this%start(2) ) * sin( this%angles(1) ) * cos( this%angles(2) ) - ( this%focus(3) - this%start(3) ) * sin( this%angles(2) )


 else if( any( wall_injection ) ) then
  select case ( trim(lon_type) )
  case ("plane")
   dlon = 0
  case ("gaussian")
   dlon = this%lon_range
  case ("sin2", "polynomial")
   dlon = this%lon_rise + this%lon_flat + this%lon_fall
  case default
   if ( mpi_node() == 0 ) then
    write(0,*) ''
    write(0,*) '   Error in zpulse_planewaves parameters'
    write(0,*) '   lon_position must be defined with wall injection '
    write(0,*) '   aborting...'
   endif
   stop
  end select
  this%lon_position = xmin( g_space, 1 ) - t_wall_start - dlon/2 - this%focus(1)
 else
  this%lon_position = 0.0_p_double
 endif

 ! ----------------------------------------------------------------------
 ! save perpendicular profile parameters
 select case ( trim(per_type) )
 case ("plane")
  this%per_type = k_per_plane
  this%per1_int_range(1) = 0.0_p_double
  this%per1_int_range(2) = 0.0_p_double
  this%per2_int_range(1) = 0.0_p_double
  this%per2_int_range(2) = 0.0_p_double
  this%dk_waves(2) = 1
  this%dk_waves(3) = 1

 case ("gaussian")
  this%per_type = k_per_gaussian

  ! get main spot size and focal plane
  this%per_w0 = 0.0_p_double
  if (per_fwhm(1) > 0.0_p_double) then
   this%per_w0 = per_fwhm(1)/sqrt(4.0_p_double * log(2.0_p_double))
  else
   this%per_w0 = per_w0(1)
  endif

  ! set integration ranges
  if ( per1_int_range(1) == 0 .and. per1_int_range(2) == 0 ) then
   sig = sqrt( 2.0_p_double ) / this%per_w0(1)
   this%per1_int_range = (/ - n_sig_gaussian*sig, n_sig_gaussian*sig /)
  else
   this%per1_int_range = per1_int_range
  endif

  if ( per2_int_range(1) == 0 .and. per2_int_range(2) == 0 ) then
   sig = sqrt( 2.0_p_double ) / this%per_w0(1)
   this%per2_int_range = (/ - n_sig_gaussian*sig, n_sig_gaussian*sig /)
  else
   this%per2_int_range = per2_int_range
  endif


 case ("laguerre")
  this%per_type = k_per_laguerre_gaussian
  if ( p_x_dim > 1 ) then
   this%per_tem_mode = 0
   do i1 = 1, 2
    this%per_tem_mode(i1) = per_tem_mode(i1)
    if ( per_tem_mode(i1) > 7 ) then
     if ( mpi_node() == 0 ) then
      write(0,*) ''
      write(0,*) '   Error in zpulse-planewaves parameters'
      write(0,*) '   Only TEM modes up to 7 have been implemented. Please contact the'
      write(0,*) '   development team if you need higher order modes.'
      write(0,*) '   aborting...'
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

  ! set integration ranges
  if ( per1_int_range(1) == 0 .and. per1_int_range(2) == 0 ) then
   sig = ( sqrt( 2.0_p_double ) + 0.2*(this%per_tem_mode(1)+this%per_tem_mode(2)) ) / this%per_w0(1)
   this%per1_int_range = (/ - n_sig_gaussian*sig, n_sig_gaussian*sig /)
  else
   this%per1_int_range = per1_int_range
  endif

  if ( per2_int_range(1) == 0 .and. per2_int_range(2) == 0 ) then
   sig = ( sqrt( 2.0_p_double ) + 0.2*(this%per_tem_mode(1)+this%per_tem_mode(2)) ) / this%per_w0(1)
   this%per2_int_range = (/ - n_sig_gaussian*sig, n_sig_gaussian*sig /)
  else
   this%per2_int_range = per2_int_range
  endif


 case ("hermite")
  this%per_type = k_per_hermite_gaussian
  this%per_tem_mode = 0
  this%per_tem_mode = per_tem_mode
  if ( per_tem_mode(1) > 7 .or. per_tem_mode(2) > 7 ) then
   if ( mpi_node() == 0 ) then
    print *, ''
    print *, '   Error in zpulse-planewaves parameters'
    print *, '   Only Laguerre modes up to 7 have been implemented. Please contact the'
    print *, '   development team if you need higher order modes.'
    print *, '   aborting...'
   endif
   stop
  endif

  ! get main spot size and focal plane
  this%per_w0 = 0.0_p_double
  if (per_fwhm(1) > 0.0_p_double) then
   this%per_w0 = per_fwhm(1)/sqrt(4.0_p_double * log(2.0_p_double))
  else
   this%per_w0 = per_w0(1)
  endif

  ! set integration ranges
  if ( per1_int_range(1) == 0 .and. per1_int_range(2) == 0 ) then
   sig = ( sqrt( 2.0_p_double ) + 0.2*this%per_tem_mode(1) ) / this%per_w0(1)
   this%per1_int_range = (/ - n_sig_gaussian*sig, n_sig_gaussian*sig /)
  else
   this%per1_int_range = per1_int_range
  endif

  if ( per2_int_range(1) == 0 .and. per2_int_range(2) == 0 ) then
   sig = ( sqrt( 2.0_p_double ) + 0.2*this%per_tem_mode(2) ) / this%per_w0(1)
   this%per2_int_range = (/ - n_sig_gaussian*sig, n_sig_gaussian*sig /)
  else
   this%per2_int_range = per2_int_range
  endif

 case default
  if ( mpi_node() == 0 ) then
   write(0,*) ''
   write(0,*) '   Error in zpulse_planewaves parameters'
   write(0,*) '   per_type must be "gaussian" or "laguerre"'
   write(0,*) '   no other type has been implemented yet '
   write(0,*) '   aborting...'
  endif
  stop
 end select

 ! ----------------------------------------------------------------------
 ! extra parameters
 this%if_launch = if_launch
 this%launch_time = launch_time

 ! parameters for wall injection (start and stop times)
 this%wall_injection = wall_injection
 this%initial_injection = initial_injection
 if ( any(this%wall_injection) ) then
  if ( t_wall_start /= t_wall_stop ) then
   this%t_wall_start = t_wall_start
   this%t_wall_stop = t_wall_stop
  else
   select case ( trim(lon_type) )
   case ("gaussian")
    dlon = this%lon_range
   case("sin2", "polynomial")
    dlon = this%lon_rise + this%lon_flat + this%lon_fall
   case default
    if ( mpi_node() == 0 ) then
     write(0,*) ''
     write(0,*) '   Error in zpulse_planewaves parameters'
     write(0,*) '   "t_wall_start" and "t_wall_stop" must be defined'
     write(0,*) '   for wall injection to work '
     write(0,*) '   aborting...'
    endif
    stop
   end select
   this%t_wall_start = 0.0_p_double
   this%t_wall_stop = 1.01*dlon
   this%dk_waves(1) = 2*pi / ( this%t_wall_stop - this%t_wall_start )
  endif
  do i1 = 1, p_x_dim
   if ( periodicity(i1) == 0.0_p_double ) then
    this%dk_waves(i1) = this%dk_waves(i1) / 2
   endif
  enddo
 endif

 ! ----------------------------------------------------------------------
 ! adapt integration range according to pulse type
 !if ( this%per_type == k_per_plane .and. this%lon_type == k_lon_plane ) then do nothing
 !if ( this%pulse_type == type_space ) then do nothing
 if ( this%pulse_type == type_custom_full .or. this%pulse_type == type_custom_lon ) then
  ! longitudinal integration range
  if ( lon_int_range(1) == 0 .and. lon_int_range(2) == 0 ) then
   if ( mpi_node() == 0 ) then
    write(0,*) ''
    write(0,*) '   Error in zpulse_planewaves parameters'
    write(0,*) '   lon_int_range needs to be specified when pulse type is "custom" '
    write(0,*) '   aborting...'
   endif
   stop
  endif
  this%lon_int_range = lon_int_range
 endif
 if ( this%pulse_type == type_custom_full .or. this%pulse_type == type_custom_per ) then
  ! transverse integration range 2D
  if ( per1_int_range(1) == 0 .and. per1_int_range(2) == 0 ) then
   if ( mpi_node() == 0 ) then
    write(0,*) ''
    write(0,*) '   Error in zpulse_planewaves parameters'
    write(0,*) '   per1_int_range needs to be specified when pulse type is "custom" '
    write(0,*) '   aborting...'
   endif
   stop
  endif
  this%per1_int_range = per1_int_range
  ! transverse integration range 3D
  if ( per2_int_range(1) == 0 .and. per2_int_range(2) == 0 .and. p_x_dim == 3 ) then
   if ( mpi_node() == 0 ) then
    write(0,*) ''
    write(0,*) '   Error in zpulse_planewaves parameters'
    write(0,*) '   per2_int_range needs to be specified when pulse type is "custom" '
    write(0,*) '   aborting...'
   endif
   stop
  endif
  this%per2_int_range = per2_int_range
 endif

 if ( this%pulse_type == type_time .or. this%pulse_type == type_vfocus_const .or. this%pulse_type == type_vfocus_accel .or. this%pulse_type == type_custom_per ) then
  if ( per1_int_range(1) == 0 .and. per1_int_range(2) == 0 .and. per2_int_range(1) == 0 .and. per2_int_range(2) == 0 .and.lon_int_range(1) == 0 .and. lon_int_range(2) == 0 ) then

   dper = min( this%per1_int_range(2), this%omega0 )
   dlon = this%lon_int_range(2) - this%omega0

   if ( this%pulse_type == type_time ) then
    this%lon_int_range(1) = sqrt( abs(this%omega0**2 - dper**2) ) - dlon
    this%lon_int_range(2) = this%omega0 + dlon
   else if ( this%vfocus < 0.0_p_double ) then
    !if ( dper**2*(this%vfocus**2-1) + (1-this%vfocus)**2*this%omega0**2 > 0 ) then
    ! this%lon_int_range(1) = ( (1-this%vfocus)*this%omega0*this%vfocus + sqrt( dper**2*(this%vfocus**2-1) + (1-this%vfocus)**2*this%omega0**2 ) ) / (1-this%vfocus**2) - dlon
    !else
    ! this%lon_int_range(1) = 0
    !endif
    this%lon_int_range(1) = 0
    this%lon_int_range(2) = this%omega0 + dlon
   else if ( this%vfocus < 1.0_p_double ) then
    this%lon_int_range(1) = this%vfocus / (1+this%vfocus) * this%omega0
    this%lon_int_range(2) = this%omega0 + dlon
   else
    this%lon_int_range(1) = this%omega0 - dlon
    this%lon_int_range(2) = ( (1-this%vfocus)*this%omega0*this%vfocus - sqrt( dper**2*(this%vfocus**2-1) + (1-this%vfocus)**2*this%omega0**2 ) ) / (1-this%vfocus**2) + dlon
   endif

  endif
 endif
 if ( this%per_type == k_per_plane ) then
  this%per1_int_range(1) = 0
  this%per1_int_range(2) = 0
  this%per2_int_range(1) = 0
  this%per2_int_range(2) = 0
  this%dk_waves(2) = 1
  this%dk_waves(3) = 1
 endif
 if ( this%lon_type == k_lon_plane ) then
  this%lon_int_range(1) = this%omega0
  this%lon_int_range(2) = this%omega0
  this%dk_waves(1) = 1
 endif

 if ( this%lon_int_range(1) < 1e-5 * this%omega0 ) then
  this%lon_int_range(1) = 1e-5 * this%omega0
 endif

 ! ----------------------------------------------------------------------
 ! normalize if needed
 if ( this%normalize ) then
  if ( mpi_node() == 0 ) then
   write(0,*) ''
   write(0,*) '   Error in zpulse_planewaves parameters'
   write(0,*) '   normalization not yet implemented '
   write(0,*) '   peak amplitude may not be accurate'
  endif
  !stop
 endif

 ! ----------------------------------------------------------------------
 ! Boost pulse
 this%if_boost = .false.
 if (sim_options%gamma > 1.0_p_double) then
  this%gamma = sim_options%gamma
  if (mpi_node()==0) print *,'   boosting zpulse_planewaves with gamma = ', this%gamma
  this%beta_boost(1) = sqrt( 1 - 1/this%gamma**2 )
  this%if_boost = .true.
 else if ( beta_boost(1) /= 0.0 .or. beta_boost(2) /= 0.0 .or. beta_boost(3) /= 0.0 ) then
  if ( sqrt( beta_boost(1)**2 + beta_boost(2)**2 + beta_boost(3)**2 ) > 1.0 ) then
   if ( mpi_node() == 0 ) then
    write(0,*) ''
    write(0,*) '   Error in zpulse_planewaves parameters'
    write(0,*) '   boosting to reference frame faster than light '
    write(0,*) '   aborting...'
   endif
   stop
  endif
  if (mpi_node()==0) print *,"   boosting zpulse_planewaves "
  this%beta_boost = beta_boost
  this%gamma = 1.0_p_double / sqrt( 1 - beta_boost(1)**2 - beta_boost(2)**2 - beta_boost(3)**2 )
  this%if_boost = .true.
 endif

end subroutine read_input_zpulse_planewaves
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
! Check if parameters are good to use with current dimensionality
subroutine check_dimensionality_zpulse_planewaves( this )

 implicit none

 class( t_zpulse_planewaves ), intent(inout) :: this

 if ( p_x_dim == 1) then
  if ( mpi_node() == 0 ) then
   write(0,*) ''
   write(0,*) '   Error in zpulse_planewaves dimensionality'
   write(0,*) '   not usable in 1d'
   write(0,*) '   aborting...'
  endif
  stop
 endif

 if ( p_x_dim < 3 .and. this%per_type == k_per_laguerre_gaussian ) then
  if ( mpi_node() == 0 ) then
   write(0,*) ''
   write(0,*) '   Error in zpulse_planewaves dimensionality'
   write(0,*) '   laguerre-gaussian only usable in 3d'
   write(0,*) '   aborting...'
  endif
  stop
 endif

end subroutine check_dimensionality_zpulse_planewaves
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
! Launch 2D pulse routine
subroutine launch_zpulse_planewaves_2d_normal( this, e, b, t, nx_p_min, g_space, no_co )

 implicit none

 integer, parameter :: rank = 2

 class( t_zpulse_planewaves ), intent(inout) :: this
 type( t_vdf ), intent(inout) :: e, b
 real(p_double), intent(in) :: t
 integer, dimension(:), intent(in) :: nx_p_min
 type( t_space ), intent(in) :: g_space
 class( t_node_conf ), intent(in) :: no_co ! node configuration

 ! local variables
 real(p_double), dimension(2,rank) :: g_x_range
 real(p_double), dimension(rank) :: ldx

 real(p_double) :: x1min, dx1_2
 real(p_double) :: x2min, dx2_2

 ! other variables
 real(p_double) :: k1, k2
 complex(p_double) :: e1_fac, e2_fac, e3_fac, b1_fac, b2_fac, b3_fac
 integer :: i1, i2
 complex(p_double) :: exp_x1, exp_x2, exp_x10, exp_x20
 integer :: border_x1_up, border_x1_down, border_x2_up, border_x2_down

 integer :: total_n_waves, iw

 ! executable statements
 call get_x_bnd( g_space, g_x_range )
 ldx(1) = b%dx_(1)
 ldx(2) = b%dx_(2)

 ! get border if wall injection
 border_x1_down = b%lbound(2)
 if ( this%wall_injection(1,1) .and. my_ngp( no_co, 1) == 1 ) border_x1_down = 2

 border_x1_up = b%ubound(2)
 if ( this%wall_injection(1,2) .and. my_ngp( no_co, 1) == nx( no_co, 1 ) ) border_x1_up = b%ubound(2) - 4

 border_x2_down = b%lbound(3)
 if ( this%wall_injection(2,1) .and. my_ngp( no_co, 2) == 1 ) border_x2_down = 2

 border_x2_up = b%ubound(3)
 if ( this%wall_injection(2,2) .and. my_ngp( no_co, 2) == nx( no_co, 2 ) ) border_x2_up = b%ubound(3) - 3

 ! spatial coordinates - xmin already includes focus position
 x1min = g_x_range(p_lower,1) + real(nx_p_min(1)-1, p_double)*ldx(1) - this%focus(1)
 x2min = g_x_range(p_lower,2) + real(nx_p_min(2)-1, p_double)*ldx(2) - this%focus(2)
 dx1_2 = ldx(1)/2.0_p_double
 dx2_2 = ldx(2)/2.0_p_double

 ! get actual total number of waves to inject and iterate over them
 total_n_waves = size( this%kvec_arr, 1 )
 do iw = 1, total_n_waves

  k1 = this%kvec_arr( iw, 1 )
  k2 = this%kvec_arr( iw, 2 )
  e1_fac = this%evec_arr( iw, 1 ) * exp( im_unit * ( - sqrt( k1**2 + k2**2 ) * t + k1 * (x1min+dx1_2) + k2 * x2min ) )
  e2_fac = this%evec_arr( iw, 2 ) * exp( im_unit * ( - sqrt( k1**2 + k2**2 ) * t + k1 * x1min + k2 * (x2min+dx2_2) ) )
  e3_fac = this%evec_arr( iw, 3 ) * exp( im_unit * ( - sqrt( k1**2 + k2**2 ) * t + k1 * x1min + k2 * x2min ) )
  b1_fac = this%bvec_arr( iw, 1 ) * exp( im_unit * ( - sqrt( k1**2 + k2**2 ) * t + k1 * x1min + k2 * (x2min+dx2_2) ) )
  b2_fac = this%bvec_arr( iw, 2 ) * exp( im_unit * ( - sqrt( k1**2 + k2**2 ) * t + k1 * (x1min+dx1_2) + k2 * x2min ) )
  b3_fac = this%bvec_arr( iw, 3 ) * exp( im_unit * ( - sqrt( k1**2 + k2**2 ) * t + k1 * (x1min+dx1_2) + k2 * (x2min+dx2_2) ) )

  exp_x10 = exp( im_unit * k1 * ldx(1) )
  exp_x20 = exp( im_unit * k2 * ldx(2) )
  ! iterate for every space point
  exp_x1 = exp( im_unit * k1 * (border_x1_down-1) * ldx(1) )
  do i1 = border_x1_down, border_x1_up
   ! calculate grid positions and exponentials for x1
   exp_x1 = exp_x1 * exp_x10
  exp_x2 = exp( im_unit * k2 * (border_x2_down-1) * ldx(2) )
  do i2 = border_x2_down, border_x2_up
   ! calculate grid positions and exponentials for x2
   exp_x2 = exp_x2 * exp_x20
  ! --------------------------------------------------

   ! e3 at (x1, x2)
   e%f2(3,i1,i2) = e%f2(3,i1,i2) + real( exp_x1 * exp_x2 * e3_fac, p_double)

   ! b3 at (x1_2, x2_2)
   b%f2(3,i1,i2) = b%f2(3,i1,i2) + real( exp_x1 * exp_x2 * b3_fac, p_double)

   ! e1, b2 at (x1_2, x2)
   e%f2(1,i1,i2) = e%f2(1,i1,i2) + real( exp_x1 * exp_x2 * e1_fac, p_double)
   b%f2(2,i1,i2) = b%f2(2,i1,i2) + real( exp_x1 * exp_x2 * b2_fac, p_double)

   ! e2, b1 at (x1, x2_2)
   e%f2(2,i1,i2) = e%f2(2,i1,i2) + real( exp_x1 * exp_x2 * e2_fac, p_double)
   b%f2(1,i1,i2) = b%f2(1,i1,i2) + real( exp_x1 * exp_x2 * b1_fac, p_double)

  enddo
  enddo

  ! testing, debug and control
  if ( control .and. mpi_node() == 0 ) then
   if ( mod( iw, total_n_waves/20) == 0 ) then
    write(0,*) " - injected ", iw * 100.0 /total_n_waves, " % of waves"
   endif
  endif
 enddo

end subroutine launch_zpulse_planewaves_2d_normal
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
! Launch 3D pulse routine
subroutine launch_zpulse_planewaves_3d_normal( this, e, b, t, nx_p_min, g_space, no_co )

 implicit none

 integer, parameter :: rank = 3

 class( t_zpulse_planewaves ), intent(inout) :: this
 type( t_vdf ), intent(inout) :: e, b
 real(p_double), intent(in) :: t
 integer, dimension(:), intent(in) :: nx_p_min
 type( t_space ), intent(in) :: g_space
 class( t_node_conf ), intent(in) :: no_co ! node configuration

 ! local variables
 real(p_double), dimension(2,rank) :: g_x_range
 real(p_double), dimension(rank) :: ldx

 real(p_double) :: x1min, dx1_2
 real(p_double) :: x2min, dx2_2
 real(p_double) :: x3min, dx3_2

 ! other variables
 real(p_double) :: k1, k2, k3
 integer :: i1, i2, i3, iw, total_n_waves
 complex(p_double) :: exp_x1, exp_x2, exp_x3, exp_x10, exp_x20, exp_x30, exp_tot
 complex(p_double) :: e1_fac, e2_fac, e3_fac, b1_fac, b2_fac, b3_fac
 integer :: border_x1_up, border_x1_down, border_x2_up, border_x2_down, border_x3_up, border_x3_down

 ! executable statements
 call get_x_bnd( g_space, g_x_range )
 ldx(1) = b%dx_(1)
 ldx(2) = b%dx_(2)
 ldx(3) = b%dx_(3)

 ! get border if wall injection
 border_x1_up = b%ubound(2)
 border_x1_down = b%lbound(2)
 if ( this%wall_injection(1,1) .and. my_ngp( no_co, 1) == 1 ) then
  border_x1_down = 2
 endif
 if ( this%wall_injection(1,1) .and. my_ngp( no_co, 1) == nx( no_co, 1 ) ) then
  border_x1_up = b%ubound(2) - 1
 endif
 border_x2_up = b%ubound(3)
 border_x2_down = b%lbound(3)
 if ( this%wall_injection(2,1) .and. my_ngp( no_co, 2) == 1 ) then
  border_x2_down = 2
 endif
 if ( this%wall_injection(2,1) .and. my_ngp( no_co, 2) == nx( no_co, 2 ) ) then
  border_x2_up = b%ubound(3) - 1
 endif
 border_x3_up = b%ubound(4)
 border_x3_down = b%lbound(4)
 if ( this%wall_injection(3,1) .and. my_ngp( no_co, 3) == 1 ) then
  border_x3_down = 2
 endif
 if ( this%wall_injection(3,1) .and. my_ngp( no_co, 3) == nx( no_co, 3 ) ) then
  border_x3_up = b%ubound(4) - 1
 endif

 ! spatial coordinates
 ! XMIN ALREADY INCLUDES FOCUS POSITION
 x1min = g_x_range(p_lower,1) + real(nx_p_min(1)-1, p_double)*ldx(1) - this%focus(1)
 x2min = g_x_range(p_lower,2) + real(nx_p_min(2)-1, p_double)*ldx(2) - this%focus(2)
 x3min = g_x_range(p_lower,3) + real(nx_p_min(3)-1, p_double)*ldx(3) - this%focus(3)
 dx1_2 = ldx(1)/2.0_p_double
 dx2_2 = ldx(2)/2.0_p_double
 dx3_2 = ldx(3)/2.0_p_double

 ! get actual total number of waves to inject and iterate over them
 total_n_waves = size( this%kvec_arr, 1 )
 do iw = 1, total_n_waves

  k1 = this%kvec_arr( iw, 1 )
  k2 = this%kvec_arr( iw, 2 )
  k3 = this%kvec_arr( iw, 3 )
  e1_fac = this%evec_arr( iw, 1 ) * exp( im_unit * ( -sqrt( k1**2 + k2**2 + k3**2 ) * t + k1 * (x1min+dx1_2) + k2 * x2min + k3 * x3min ) )
  e2_fac = this%evec_arr( iw, 2 ) * exp( im_unit * ( -sqrt( k1**2 + k2**2 + k3**2 ) * t + k1 * x1min + k2 * (x2min+dx2_2) + k3 * x3min ) )
  e3_fac = this%evec_arr( iw, 3 ) * exp( im_unit * ( -sqrt( k1**2 + k2**2 + k3**2 ) * t + k1 * x1min + k2 * x2min + k3 * (x3min+dx3_2) ) )
  b1_fac = this%bvec_arr( iw, 1 ) * exp( im_unit * ( -sqrt( k1**2 + k2**2 + k3**2 ) * t + k1 * x1min + k2 * (x2min+dx2_2) + k3 * (x3min+dx3_2) ) )
  b2_fac = this%bvec_arr( iw, 2 ) * exp( im_unit * ( -sqrt( k1**2 + k2**2 + k3**2 ) * t + k1 * (x1min+dx1_2) + k2 * x2min + k3 * (x3min+dx3_2) ) )
  b3_fac = this%bvec_arr( iw, 3 ) * exp( im_unit * ( -sqrt( k1**2 + k2**2 + k3**2 ) * t + k1 * (x1min+dx1_2) + k2 * (x2min+dx2_2) + k3 * x3min ) )

  exp_x10 = exp( im_unit * k1 * ldx(1) )
  exp_x20 = exp( im_unit * k2 * ldx(2) )
  exp_x30 = exp( im_unit * k3 * ldx(3) )

  ! iterate for every space point
  exp_x1 = exp( im_unit * k1 * (border_x1_down-1)*ldx(1) )
  do i1 = border_x1_down, border_x1_up
   ! calculate grid positions and exponentials for x1
   exp_x1 = exp_x1 * exp_x10
  exp_x2 = exp( im_unit * k2 * (border_x2_down-1)*ldx(2) )
  do i2 = border_x2_down, border_x2_up
   ! calculate grid positions and exponentials for x2
   exp_x2 = exp_x2 * exp_x20
  exp_x3 = exp( im_unit * k3 * (border_x3_down-1)*ldx(3) )
  do i3 = border_x3_down, border_x3_up
   ! calculate grid positions and exponentials for x3
   exp_x3 = exp_x3 * exp_x30
   exp_tot = exp_x1 * exp_x2 * exp_x3
  ! --------------------------------------------------

    ! e1 at (x1_2, x2 , x3 )
    e%f3(1,i1,i2,i3) = e%f3(1,i1,i2,i3) + real( exp_tot * e1_fac, p_double)
    ! e2 at (x1 , x2_2, x3 )
    e%f3(2,i1,i2,i3) = e%f3(2,i1,i2,i3) + real( exp_tot * e2_fac, p_double)
    ! e3 at (x1 , x2 , x3_2)
    e%f3(3,i1,i2,i3) = e%f3(3,i1,i2,i3) + real( exp_tot * e3_fac, p_double)
    ! b1 at (x1 , x2_2, x3_2)
    b%f3(1,i1,i2,i3) = b%f3(1,i1,i2,i3) + real( exp_tot * b1_fac, p_double)
    ! b2 at (x1_2, x2 , x3_2)
    b%f3(2,i1,i2,i3) = b%f3(2,i1,i2,i3) + real( exp_tot * b2_fac, p_double)
    ! b3 at (x1_2, x2_2, x3 )
    b%f3(3,i1,i2,i3) = b%f3(3,i1,i2,i3) + real( exp_tot * b3_fac, p_double)

  enddo
  enddo
  enddo
  ! testing, debug and control
  if ( control .and. mpi_node() == 0 ) then
   if ( mod( iw, total_n_waves/20) == 0 ) then
    write(0,*) " - injected ", iw * 100.0 /total_n_waves, " % of waves"
   endif
  endif
 enddo

end subroutine launch_zpulse_planewaves_3d_normal
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
! Launch pulse routine
subroutine launch_zpulse_planewaves( this, emf, g_space, nx_p_min, g_nx, t, dt, no_co )

 class( t_zpulse_planewaves ), intent(inout) :: this
 class( t_emf ) , intent(inout) :: emf

 type( t_space ), intent(in) :: g_space

 integer, intent(in), dimension(:) :: nx_p_min, g_nx
 real(p_double), intent(in) :: t
 real(p_double), intent(in) :: dt
 class( t_node_conf ), intent(in) :: no_co

 ! local variables
 type( t_vdf ) :: e_pulse, b_pulse
 type( t_space ) :: lg_space
 integer, dimension(p_x_dim) :: lnx_p_min, lg_nx

 ! testing, debug and control
 real(p_double) :: start_time, end_time
 real(p_double), dimension(2,p_x_dim) :: g_x_range

 ! launch zpulse
 if ( this%if_launch .and. .not. this%was_launched .and. ( ( t >= this%launch_time .and. .not. any(this%wall_injection) ) .or. ( this%initial_injection .and. any(this%wall_injection) ) ) ) then

  ! testing, debug and control
  if ( control ) then
   call cpu_time(start_time)
  endif

  if (mpi_node()==0) print *," - injecting zpulse_planewaves"

  ! initialize local variables
  lg_space = g_space
  lnx_p_min = nx_p_min(1:p_x_dim)
  lg_nx = g_nx(1:p_x_dim)

  ! create electric and magnetic field vdfs
  call e_pulse % new( emf%e, zero = .true. )
  call b_pulse % new( emf%b, zero = .true. )

  ! choose dimension
  select case (p_x_dim)
  case (2)
   call launch_zpulse_planewaves_2d_normal( this, e_pulse, b_pulse, t, nx_p_min, g_space, no_co )
  case (3)
   call launch_zpulse_planewaves_3d_normal( this, e_pulse, b_pulse, t, nx_p_min, g_space, no_co )
  end select

  ! add the created pulses to the simulation
  call add( emf%e, e_pulse )
  call add( emf%b, b_pulse )

  ! clear memory
  call e_pulse % cleanup()
  call b_pulse % cleanup()

  ! zpulse has been launched
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
  select case (p_x_dim)
  case (2)
   call launch_zpulse_planewaves_wall_2d( this, emf%e, emf%b, g_space, nx_p_min, t, dt, no_co )
  case (3)
   call launch_zpulse_planewaves_wall_3d( this, emf%e, emf%b, g_space, nx_p_min, t, dt, no_co )
  end select
 endif

end subroutine launch_zpulse_planewaves
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
! Wave vector weight function
function k_function_zpulse_planewaves( this, k_lon, k_per1, k_per2 ) result( f_k )

 implicit none
 class( t_zpulse_planewaves ), intent(inout):: this
 real(p_double), intent(in) :: k_lon, k_per1, k_per2
 complex(p_double) :: f_k
 real(p_k_fparse), dimension(3) :: k_dbl

 ! auxiliar
 real(p_double) :: vg, dvg, vg_range, vgmin, vgmax

 !f_k_lon_arg = sign( sqrt( abs( k_lon**2 + (k_per1**2+k_per2**2) * this%abouraddy_factor ) ), k_lon**2 + (k_per1**2+k_per2**2) * this%abouraddy_factor )
 if ( this%pulse_type == type_time ) then
  ! apply perpendicular profile
  f_k = this%k_per_function( k_per1, k_per2 )
  ! waves centered around \omega
  f_k = f_k * this%k_lon_function( sqrt( k_lon**2 + k_per1**2 + k_per2**2 ) )

 else if ( this%pulse_type == type_space ) then
  ! apply perpendicular profile
  f_k = this%k_per_function( k_per1, k_per2 )
  ! waves centered around klon
  f_k = f_k * this%k_lon_function( k_lon )

 else if ( this%pulse_type == type_vfocus_const ) then
  ! apply perpendicular profile
  f_k = this%k_per_function( k_per1, k_per2 )
  ! waves centered around isosurface of contant focus velocity
  f_k = f_k * this%k_lon_function( ( -k_lon * this%vfocus + sqrt( k_lon**2 + (k_per1**2+k_per2**2) ) ) / ( 1 - this%vfocus ) )

 else if ( this%pulse_type == type_vfocus_accel ) then
  ! apply perpendicular profile
  f_k = this%k_per_function( k_per1, k_per2 )

  ! calculate velocity associated with wave
  if ( k_lon == this%omega0 .and. k_per1 == 0.0 .and. k_per2 == 0.0 ) then
   vg = this%vfocus
  else
   vg = ( sqrt(k_lon**2 + k_per1**2 + k_per2**2) - this%omega0 ) / (k_lon - this%omega0)
  endif
  ! apply corresponding phase
  f_k = f_k * exp( -im_unit * (sqrt(k_lon**2 + k_per1**2 + k_per2**2) * this%afocus * (vg-this%vfocus))**2 )
  !f_k = f_k * exp( -im_unit * this%afocus * (vg-this%vfocus) )

  ! maximum integration beta
  dvg = 0.01
  vg_range = 0.03
  vgmin = this%vfocus
  vgmax = this%vfocus + vg_range
  ! apply h(\beta) function
  if ( vg < vgmin - dvg ) then
   f_k = f_k * 0.0
  else if ( vg < vgmin ) then
   f_k = f_k * (1.0 - (vgmin - vg) / dvg)
  else if ( vg < vgmax ) then
   f_k = f_k * 1.0
  else if ( vg < vgmax + dvg ) then
   f_k = f_k * (1 - (vg - vgmax) / dvg)
  else
   f_k = f_k * 0.0
  endif


 else if ( this%pulse_type == type_custom_full ) then
  k_dbl(1) = real( k_lon , p_k_fparse )
  k_dbl(2) = real( k_per1, p_k_fparse )
  k_dbl(3) = real( k_per2, p_k_fparse )
       f_k = eval( this%f_k_amp_math_func, k_dbl ) * exp( im_unit * eval( this%f_k_pha_math_func, k_dbl ) )

 else if ( this%pulse_type == type_custom_lon ) then
  k_dbl(1) = real( k_lon , p_k_fparse )
  k_dbl(2) = real( k_per1, p_k_fparse )
  k_dbl(3) = real( k_per2, p_k_fparse )
       f_k = eval( this%f_k_amp_math_func, k_dbl ) * exp( im_unit * eval( this%f_k_pha_math_func, k_dbl ) )
  f_k = f_k * this%k_per_function( k_per1, k_per2 )

 else if ( this%pulse_type == type_custom_per ) then
  k_dbl(1) = real( k_lon , p_k_fparse )
  k_dbl(2) = real( k_per1, p_k_fparse )
  k_dbl(3) = real( k_per2, p_k_fparse )
       f_k = eval( this%f_k_amp_math_func, k_dbl ) * exp( im_unit * eval( this%f_k_pha_math_func, k_dbl ) )
  f_k = f_k * this%k_lon_function( sqrt( k_lon**2 + k_per1**2 + k_per2**2 ) )

 endif

end function k_function_zpulse_planewaves
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
! Transverse wave vector weight function
function k_per_function_zpulse_planewaves( this, k_per1, k_per2 ) result( f_k_per )

 implicit none
 class( t_zpulse_planewaves ), intent(in) :: this
 real(p_double), intent(in) :: k_per1, k_per2
 complex(p_double) :: f_k_per

 ! auxiliar variables
 real(p_double) :: sig, kr, ktheta, w0
 integer :: p, l

 select case ( p_x_dim )
  case(2)
   select case ( this%per_type )
   case (k_per_plane)
    f_k_per = 1.0_p_double
   case (k_per_gaussian)

    sig = sqrt( 2.0_p_double ) / this%per_w0(1)
    f_k_per = 1.0 / sqrt( 2.0*pi ) / sig * exp( -k_per1**2/2/sig**2 )
    !if (mpi_node()==0) print *,k_per1, f_k_per

   case (k_per_hermite_gaussian)

    sig = sqrt( 2.0_p_double ) / this%per_w0(1)
    f_k_per = 1.0 / sqrt( 2.0*pi ) / sig * exp( -k_per1**2/2/sig**2 ) * hermite( this%per_tem_mode(1), k_per1/sig ) * (-im_unit)**this%per_tem_mode(1)

   end select

  case(3)
   select case ( this%per_type )
   case (k_per_plane)
    f_k_per = 1.0_p_double
   case (k_per_gaussian)

    f_k_per = this%per_w0(1)**2 / ( 4.0*pi ) * exp( -( k_per1**2 + k_per2**2 ) * this%per_w0(1)**2 / 4 )

   case (k_per_hermite_gaussian)

    sig = sqrt( 2.0_p_double ) / this%per_w0(1)
    f_k_per = 1.0 / sqrt( 2.0*pi ) / sig * exp( -k_per1**2/2/sig**2 ) * hermite( this%per_tem_mode(1), k_per1/sig ) * (-im_unit)**this%per_tem_mode(1)
    f_k_per = f_k_per * 1.0 / sqrt( 2.0*pi ) / sig * exp( -k_per2**2/2/sig**2 ) * hermite( this%per_tem_mode(2), k_per2/sig ) * (-im_unit)**this%per_tem_mode(2)

   case (k_per_laguerre_gaussian)

    kr = sqrt( k_per1**2 + k_per2**2 )
    ktheta = atan2( k_per2, k_per1 )
    p = this%per_tem_mode(1)
    l = this%per_tem_mode(2)
    w0 = this%per_w0(1)

    f_k_per = w0**2 * (-im_unit)**l * (-1)**(p-l) / (4*pi) * 2.0
    f_k_per = f_k_per * (kr*w0/sqrt(2.0_p_double))**abs(l)
    f_k_per = f_k_per * laguerre( p, abs(l), kr**2*w0**2/2 )
    f_k_per = f_k_per * exp( -kr**2*w0**2/4 + im_unit*l*ktheta )

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
    end select

 end select

end function k_per_function_zpulse_planewaves
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
! Longitudinal wave vector weight function
function k_lon_function_zpulse_planewaves( this, k ) result( f_k_lon )

 use, intrinsic :: ieee_arithmetic

 implicit none
 class( t_zpulse_planewaves ), intent(in) :: this
 real(p_double), intent(in) :: k
 real(p_double) :: k_lon
 complex(p_double) :: f_k_lon

 ! auxiliar variables
 real(p_double) :: sig, real_part, imag_part

 k_lon = k

 select case ( this%lon_type )
 case (k_lon_plane)
  f_k_lon = 1.0_p_double
 case (k_lon_gaussian)

  sig = 2.0 * sqrt( 2.0 * log(2.0) ) / this%lon_duration
  f_k_lon = 1.0 / sqrt( 2.0*pi ) / sig * exp( -(k_lon-this%omega0)**2/2/sig**2 )

 case (k_lon_sin2)

  real_part = (pi**1.5*((-(this%lon_fall**2*(k_lon-this%omega0)**2) + pi**2)*sin(((k_lon-this%omega0)*(this%lon_fall + this%lon_flat - this%lon_rise))/2.) + (-(this%lon_fall**2*(k_lon-this%omega0)**2) + 2*pi**2 - (k_lon-this%omega0)**2*this%lon_rise**2)*sin(((k_lon-this%omega0)*(this%lon_fall + this%lon_flat + this%lon_rise))/2.) + (-pi**2 + (k_lon-this%omega0)**2*this%lon_rise**2)*sin(this%lon_fall*(k_lon-this%omega0) - ((k_lon-this%omega0)*(this%lon_fall + this%lon_flat + this%lon_rise))/2.)))/(2.*sqrt(2.0_p_double)*(k_lon-this%omega0)*(this%lon_fall*(k_lon-this%omega0) - pi)*(this%lon_fall*(k_lon-this%omega0) + pi)*(-pi + (k_lon-this%omega0)*this%lon_rise)*(pi + (k_lon-this%omega0)*this%lon_rise))
  imag_part = (pi**1.5*((this%lon_fall**2*(k_lon-this%omega0)**2 - pi**2)*cos(((k_lon-this%omega0)*(this%lon_fall + this%lon_flat - this%lon_rise))/2.) + (k_lon-this%omega0)**2*(this%lon_fall**2 - this%lon_rise**2)*cos(((k_lon-this%omega0)*(this%lon_fall + this%lon_flat + this%lon_rise))/2.) + (pi**2 - (k_lon-this%omega0)**2*this%lon_rise**2)*cos(this%lon_fall*(k_lon-this%omega0) - ((k_lon-this%omega0)*(this%lon_fall + this%lon_flat + this%lon_rise))/2.)))/(2.*sqrt(2.0_p_double)*(k_lon-this%omega0)*(this%lon_fall*(k_lon-this%omega0) - pi)*(this%lon_fall*(k_lon-this%omega0) + pi)*(-pi + (k_lon-this%omega0)*this%lon_rise)*(pi + (k_lon-this%omega0)*this%lon_rise))
  f_k_lon = (real_part + im_unit * imag_part)

  if ( k_lon-this%omega0 == 0 ) f_k_lon = (this%lon_fall + 2.*this%lon_flat + this%lon_rise)/(2.*sqrt(2.*pi))

  f_k_lon = f_k_lon / sqrt(2.0_p_double * pi)

 case (k_lon_polynomial)

  f_k_lon = ((0,-1)*this%lon_fall*(k_lon-this%omega0)*(720. + this%lon_fall**4*(k_lon-this%omega0)**4)*this%lon_rise**5*cos((this%lon_fall*(k_lon-this%omega0))/2.) + (0,1)*exp((im_unit*(k_lon-this%omega0)*(this%lon_fall + 2*this%lon_flat + this%lon_rise))/2.)*this%lon_fall**5*(k_lon-this%omega0)*this%lon_rise*(720. + (k_lon-this%omega0)**4*this%lon_rise**4)*cos(((k_lon-this%omega0)*this%lon_rise)/2.) + (0,1440)*this%lon_rise**5*sin((this%lon_fall*(k_lon-this%omega0))/2.) - (0,120.)*this%lon_fall**2*(k_lon-this%omega0)**2*this%lon_rise**5*sin((this%lon_fall*(k_lon-this%omega0))/2.) + this%lon_fall**5*(k_lon-this%omega0)**5*this%lon_rise**5*sin((this%lon_fall*(k_lon-this%omega0))/2.) + 2*exp((im_unit*(this%lon_fall + this%lon_flat)*(k_lon-this%omega0))/2.)*this%lon_fall**5*(k_lon-this%omega0)**5*this%lon_rise**5*sin((this%lon_flat*(k_lon-this%omega0))/2.) - (0,1440.)*exp((im_unit*(k_lon-this%omega0)*(this%lon_fall + 2*this%lon_flat + this%lon_rise))/2.)*this%lon_fall**5*sin(((k_lon-this%omega0)*this%lon_rise)/2.) + (0,120)*exp((im_unit*(k_lon-this%omega0)*(this%lon_fall + 2*this%lon_flat + this%lon_rise))/2.)*this%lon_fall**5*(k_lon-this%omega0)**2*this%lon_rise**2*sin(((k_lon-this%omega0)*this%lon_rise)/2.) + exp((im_unit*(k_lon-this%omega0)*(this%lon_fall + 2*this%lon_flat + this%lon_rise))/2.)*this%lon_fall**5*(k_lon-this%omega0)**5*this%lon_rise**5*sin(((k_lon-this%omega0)*this%lon_rise)/2.)) / (exp((im_unit*(k_lon-this%omega0)*(this%lon_flat + this%lon_rise))/2.)*this%lon_fall**5*(k_lon-this%omega0)**6*sqrt(2*pi)*this%lon_rise**5)
# 2016 "zpulse/os-zpulse-planewaves.f03"
  if ( k_lon-this%omega0 == 0 ) f_k_lon = (this%lon_fall + 2.*this%lon_flat + this%lon_rise)/(2.*sqrt(2.*pi))

  f_k_lon = f_k_lon / sqrt(2.0_p_double * pi)

 end select

end function k_lon_function_zpulse_planewaves
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
! Get Jones Vector for plane wave
function get_jones_vector_zpulse_planewaves( this, k_lon, k_per1, k_per2 ) result( jones_vec )

 implicit none
 class( t_zpulse_planewaves ), intent(inout) :: this
 real(p_double), intent(in) :: k_lon, k_per1, k_per2
 complex(p_double), dimension(2) :: jones_vec
   real(p_k_fparse), dimension(3) :: k_dbl

 if ( abs(this%pol_type) <= 1 ) then ! linear and circular
  jones_vec(1) = cos( this%pol ) - this%pol_type*im_unit * sin( this%pol ) * this%pol_elipse
  jones_vec(2) = sin( this%pol ) + this%pol_type*im_unit * cos( this%pol ) * this%pol_elipse
 else if ( this%pol_type == 2 ) then ! radial
  if ( sqrt( k_per1**2 + k_per2**2 ) /= 0 ) then
   jones_vec(1) = k_per1 / sqrt( k_per1**2 + k_per2**2 )
   jones_vec(2) = k_per2 / sqrt( k_per1**2 + k_per2**2 )
  else
   jones_vec(1) = 0.0
   jones_vec(2) = 0.0
  endif
 else if ( this%pol_type == 3 ) then ! azimuthal
  if ( sqrt( k_per1**2 + k_per2**2 ) /= 0 ) then
   jones_vec(1) = - k_per2 / sqrt( k_per1**2 + k_per2**2 )
   jones_vec(2) = k_per1 / sqrt( k_per1**2 + k_per2**2 )
  else
   jones_vec(1) = 0.0
   jones_vec(2) = 0.0
  endif
 else if ( this%pol_type == 4 ) then ! custom
  k_dbl(1) = real( k_lon , p_k_fparse )
  k_dbl(2) = real( k_per1, p_k_fparse )
  k_dbl(3) = real( k_per2, p_k_fparse )
       jones_vec(1) = eval( this%jones_vec_real(1), k_dbl ) + im_unit * eval( this%jones_vec_imag(1), k_dbl )
       jones_vec(2) = eval( this%jones_vec_real(2), k_dbl ) + im_unit * eval( this%jones_vec_imag(2), k_dbl )
 endif

end function get_jones_vector_zpulse_planewaves
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
! Launch 2D pulse from wall routine
subroutine launch_zpulse_planewaves_wall_2d( this, e, b, g_space, nx_p_min, t, dt, no_co )

 implicit none

 integer, parameter :: rank = 2

 class( t_zpulse_planewaves ), intent(inout) :: this
 type( t_vdf ), intent(inout) :: e ! magnetic field
 type( t_vdf ), intent(inout) :: b ! magnetic field
 type( t_space ), intent(in) :: g_space ! global space information
 integer, intent(in), dimension(:) :: nx_p_min
 real( p_double ), intent(in) :: t ! simulation time
 real( p_double ), intent(in) :: dt ! time step
 class( t_node_conf ), intent(in) :: no_co ! node configuration

 ! local variables
 real(p_double), dimension(2,rank) :: g_x_range
 real(p_double), dimension(rank) :: ldx

 real(p_double) :: x1min, dx1_2
 real(p_double) :: x2min, dx2_2

 ! other variables
 real(p_double) :: k1, k2
 complex(p_double) :: e1_fac, e2_fac, e3_fac, b1_fac, b2_fac, b3_fac
 integer :: i1, i2
 complex(p_double) :: exp_x1, exp_x2, exp_x10, exp_x20
 integer :: border_x1_up, border_x1_down, border_x2_up, border_x2_down

 integer :: total_n_waves, iw

 ! executable statements
 call get_x_bnd( g_space, g_x_range )
 ldx(1) = b%dx_(1)
 ldx(2) = b%dx_(2)

 ! get border if wall injection
 border_x1_up = b%ubound(2)
 border_x1_down = b%lbound(2)
 if ( this%wall_injection(1,1) .and. my_ngp( no_co, 1) == 1 ) then
  border_x1_down = 1
 endif
 if ( this%wall_injection(1,1) .and. my_ngp( no_co, 1) == nx( no_co, 1 ) ) then
  border_x1_up = b%ubound(2) - 3
 endif
 border_x2_up = b%ubound(3)
 border_x2_down = b%lbound(3)
 if ( this%wall_injection(2,1) .and. my_ngp( no_co, 2) == 1 ) then
  border_x2_down = 1
 endif
 if ( this%wall_injection(2,1) .and. my_ngp( no_co, 2) == nx( no_co, 2 ) ) then
  border_x2_up = b%ubound(3) - 3
 endif

 ! spatial coordinates - xmin already includes focus position and partition location
 x1min = g_x_range(p_lower,1) + real(nx_p_min(1)-1, p_double)*ldx(1) - this%focus(1)
 x2min = g_x_range(p_lower,2) + real(nx_p_min(2)-1, p_double)*ldx(2) - this%focus(2)
 dx1_2 = ldx(1)/2.0_p_double
 dx2_2 = ldx(2)/2.0_p_double

 ! get actual total number of waves to inject and iterate over them
 total_n_waves = size( this%kvec_arr, 1 )
 do iw = 1, total_n_waves

  k1 = this%kvec_arr( iw, 1 )
  k2 = this%kvec_arr( iw, 2 )
  e1_fac = this%evec_arr( iw, 1 ) * exp( im_unit * ( - sqrt( k1**2 + k2**2 ) * (t-dt) + k1 * (x1min+dx1_2) + k2 * x2min ) )
  e2_fac = this%evec_arr( iw, 2 ) * exp( im_unit * ( - sqrt( k1**2 + k2**2 ) * (t-dt) + k1 * x1min + k2 * (x2min+dx2_2) ) )
  e3_fac = this%evec_arr( iw, 3 ) * exp( im_unit * ( - sqrt( k1**2 + k2**2 ) * (t-dt) + k1 * x1min + k2 * x2min ) )
  b1_fac = this%bvec_arr( iw, 1 ) * exp( im_unit * ( - sqrt( k1**2 + k2**2 ) * (t-dt) + k1 * x1min + k2 * (x2min+dx2_2) ) )
  b2_fac = this%bvec_arr( iw, 2 ) * exp( im_unit * ( - sqrt( k1**2 + k2**2 ) * (t-dt) + k1 * (x1min+dx1_2) + k2 * x2min ) )
  b3_fac = this%bvec_arr( iw, 3 ) * exp( im_unit * ( - sqrt( k1**2 + k2**2 ) * (t-dt) + k1 * (x1min+dx1_2) + k2 * (x2min+dx2_2) ) )

  exp_x10 = exp( im_unit * k1 * ldx(1) )
  exp_x20 = exp( im_unit * k2 * ldx(2) )

  ! x1 down wall
  if ( this%wall_injection(1,1) .and. my_ngp( no_co, 1 ) == 1 ) then
   ! iterate for every space point
   i1 = 2
   exp_x1 = exp( im_unit * k1 * (border_x1_down-1+i1)*ldx(1) )
   exp_x2 = exp( im_unit * k2 * (border_x2_down-1)*ldx(2) )
   do i2 = border_x2_down, border_x2_up
    ! calculate grid positions and exponentials for x2
    exp_x2 = exp_x2 * exp_x20
   ! --------------------------------------------------
    b%f2(2,i1,i2) = b%f2(2,i1,i2) - dt/ldx(1) * real( exp_x1 * exp_x2 * e3_fac, p_double)
    b%f2(3,i1,i2) = b%f2(3,i1,i2) + dt/ldx(1) * real( exp_x1 * exp_x2 * e2_fac, p_double)
    e%f2(2,i1,i2) = e%f2(2,i1,i2) + dt/ldx(1) * real( exp_x1 * exp_x2 * b3_fac, p_double)
    e%f2(3,i1,i2) = e%f2(3,i1,i2) - dt/ldx(1) * real( exp_x1 * exp_x2 * b2_fac, p_double)
   enddo
  endif

  ! x1 upper wall
  if ( this%wall_injection(1,2) .and. my_ngp( no_co, 1 ) == nx( no_co, 1 ) ) then
   ! iterate for every space point
   i1 = b%ubound(2) - 4
   exp_x1 = exp( im_unit * k1 * (border_x1_down+1+i1)*ldx(1) )
   exp_x2 = exp( im_unit * k2 * (border_x2_down-1)*ldx(2) )
   do i2 = border_x2_down, border_x2_up
    ! calculate grid positions and exponentials for x2
    exp_x2 = exp_x2 * exp_x20
   ! --------------------------------------------------
    b%f2(2,i1,i2) = b%f2(2,i1,i2) + dt/ldx(1) * real( exp_x1 * exp_x2 * e3_fac, p_double)
    b%f2(3,i1,i2) = b%f2(3,i1,i2) - dt/ldx(1) * real( exp_x1 * exp_x2 * e2_fac, p_double)
    e%f2(2,i1,i2) = e%f2(2,i1,i2) - dt/ldx(1) * real( exp_x1 * exp_x2 * b3_fac, p_double)
    e%f2(3,i1,i2) = e%f2(3,i1,i2) + dt/ldx(1) * real( exp_x1 * exp_x2 * b2_fac, p_double)
   enddo
  endif

  ! x2 down wall
  if ( this%wall_injection(2,1) .and. my_ngp( no_co, 2 ) == 1 ) then
   ! iterate for every space point
   i2 = 2
   exp_x1 = exp( im_unit * k1 * (border_x1_down-1)*ldx(1) )
   exp_x2 = exp( im_unit * k2 * (border_x2_down-1+i2)*ldx(2) )
   do i1 = border_x1_down, border_x1_up
    ! calculate grid positions and exponentials for x1
    exp_x1 = exp_x1 * exp_x10
   ! --------------------------------------------------
    b%f2(1,i1,i2) = b%f2(1,i1,i2) + dt/ldx(2) * real( exp_x1 * exp_x2 * e3_fac, p_double)
    b%f2(3,i1,i2) = b%f2(3,i1,i2) - dt/ldx(2) * real( exp_x1 * exp_x2 * e1_fac, p_double)
    e%f2(1,i1,i2) = e%f2(1,i1,i2) - dt/ldx(2) * real( exp_x1 * exp_x2 * b3_fac, p_double)
    e%f2(3,i1,i2) = e%f2(3,i1,i2) + dt/ldx(2) * real( exp_x1 * exp_x2 * b1_fac, p_double)
   enddo
  endif

  ! x2 upper wall
  if ( this%wall_injection(2,2) .and. my_ngp( no_co, 2 ) == nx( no_co, 2 ) ) then
   ! iterate for every space point
   i2 = b%ubound(3) - 4
   exp_x1 = exp( im_unit * k1 * (border_x1_down-1)*ldx(1) )
   exp_x2 = exp( im_unit * k2 * (border_x2_down+1+i2)*ldx(2) )
   do i1 = border_x1_down, border_x1_up
    ! calculate grid positions and exponentials for x1
    exp_x1 = exp_x1 * exp_x10
   ! --------------------------------------------------
    b%f2(1,i1,i2) = b%f2(1,i1,i2) - dt/ldx(2) * real( exp_x1 * exp_x2 * e3_fac, p_double)
    b%f2(3,i1,i2) = b%f2(3,i1,i2) + dt/ldx(2) * real( exp_x1 * exp_x2 * e1_fac, p_double)
    e%f2(1,i1,i2) = e%f2(1,i1,i2) + dt/ldx(2) * real( exp_x1 * exp_x2 * b3_fac, p_double)
    e%f2(3,i1,i2) = e%f2(3,i1,i2) - dt/ldx(2) * real( exp_x1 * exp_x2 * b1_fac, p_double)
   enddo
  endif



 enddo

end subroutine launch_zpulse_planewaves_wall_2d
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
! Launch 3D pulse from wall routine
subroutine launch_zpulse_planewaves_wall_3d( this, e, b, g_space, nx_p_min, t, dt, no_co )

 implicit none

 integer, parameter :: rank = 3

 class( t_zpulse_planewaves ), intent(inout) :: this
 type( t_vdf ), intent(inout) :: e ! magnetic field
 type( t_vdf ), intent(inout) :: b ! magnetic field
 type( t_space ), intent(in) :: g_space ! global space information
 integer, intent(in), dimension(:) :: nx_p_min
 real( p_double ), intent(in) :: t ! simulation time
 real( p_double ), intent(in) :: dt ! time step
 class( t_node_conf ), intent(in) :: no_co ! node configuration

 ! local variables
 real(p_double), dimension(2,rank) :: g_x_range
 real(p_double), dimension(rank) :: ldx

 real(p_double) :: x1min, dx1_2
 real(p_double) :: x2min, dx2_2
 real(p_double) :: x3min, dx3_2

 ! other variables
 real(p_double) :: k1, k2, k3
 integer :: i1, i2, i3, iw, total_n_waves
 complex(p_double) :: exp_x1, exp_x2, exp_x3, exp_x10, exp_x20, exp_x30
 complex(p_double) :: e1_fac, e2_fac, e3_fac, b1_fac, b2_fac, b3_fac
 integer :: border_x1_up, border_x1_down, border_x2_up, border_x2_down, border_x3_up, border_x3_down

 ! executable statements
 call get_x_bnd( g_space, g_x_range )
 ldx(1) = b%dx_(1)
 ldx(2) = b%dx_(2)
 ldx(3) = b%dx_(3)

 ! get border if wall injection
 border_x1_up = b%ubound(2)
 border_x1_down = b%lbound(2)
 border_x2_up = b%ubound(3)
 border_x2_down = b%lbound(3)
 border_x3_up = b%ubound(4)
 border_x3_down = b%lbound(4)

 ! spatial coordinates
 ! XMIN ALREADY INCLUDES FOCUS POSITION
 x1min = g_x_range(p_lower,1) + real(nx_p_min(1)-1, p_double)*ldx(1) - this%focus(1)
 x2min = g_x_range(p_lower,2) + real(nx_p_min(2)-1, p_double)*ldx(2) - this%focus(2)
 x3min = g_x_range(p_lower,3) + real(nx_p_min(3)-1, p_double)*ldx(3) - this%focus(3)
 dx1_2 = ldx(1)/2.0_p_double
 dx2_2 = ldx(2)/2.0_p_double
 dx3_2 = ldx(3)/2.0_p_double

 ! get actual total number of waves to inject and iterate over them
 total_n_waves = size( this%kvec_arr, 1 )
 do iw = 1, total_n_waves

  k1 = this%kvec_arr( iw, 1 )
  k2 = this%kvec_arr( iw, 2 )
  k3 = this%kvec_arr( iw, 3 )
  e1_fac = this%evec_arr( iw, 1 ) * exp( im_unit * ( -sqrt( k1**2 + k2**2 + k3**2 ) * (t-dt) + k1 * (x1min+dx1_2) + k2 * x2min + k3 * x3min ) )
  e2_fac = this%evec_arr( iw, 2 ) * exp( im_unit * ( -sqrt( k1**2 + k2**2 + k3**2 ) * (t-dt) + k1 * x1min + k2 * (x2min+dx2_2) + k3 * x3min ) )
  e3_fac = this%evec_arr( iw, 3 ) * exp( im_unit * ( -sqrt( k1**2 + k2**2 + k3**2 ) * (t-dt) + k1 * x1min + k2 * x2min + k3 * (x3min+dx3_2) ) )
  b1_fac = this%bvec_arr( iw, 1 ) * exp( im_unit * ( -sqrt( k1**2 + k2**2 + k3**2 ) * (t-dt) + k1 * x1min + k2 * (x2min+dx2_2) + k3 * (x3min+dx3_2) ) )
  b2_fac = this%bvec_arr( iw, 2 ) * exp( im_unit * ( -sqrt( k1**2 + k2**2 + k3**2 ) * (t-dt) + k1 * (x1min+dx1_2) + k2 * x2min + k3 * (x3min+dx3_2) ) )
  b3_fac = this%bvec_arr( iw, 3 ) * exp( im_unit * ( -sqrt( k1**2 + k2**2 + k3**2 ) * (t-dt) + k1 * (x1min+dx1_2) + k2 * (x2min+dx2_2) + k3 * x3min ) )

  exp_x10 = exp( im_unit * k1 * ldx(1) )
  exp_x20 = exp( im_unit * k2 * ldx(2) )
  exp_x30 = exp( im_unit * k3 * ldx(3) )

  ! x1 down wall
  if ( this%wall_injection(1,1) .and. my_ngp( no_co, 1 ) == 1 ) then
   ! iterate for every space point
   i1 = 2
   exp_x1 = exp( im_unit * k1 * (border_x1_down-1+i1)*ldx(1) )
   exp_x2 = exp( im_unit * k2 * (border_x2_down-1)*ldx(2) )
   do i2 = border_x2_down, border_x2_up
    ! calculate grid positions and exponentials for x2
    exp_x2 = exp_x2 * exp_x20
   exp_x3 = exp( im_unit * k3 * (border_x3_down-1)*ldx(3) )
   do i3 = border_x3_down, border_x3_up
    ! calculate grid positions and exponentials for x3
    exp_x3 = exp_x3 * exp_x30
   ! --------------------------------------------------
    b%f3(2,i1-1,i2,i3) = b%f3(2,i1-1,i2,i3) - dt/ldx(1) * real( exp_x1 * exp_x2 * exp_x3 * e3_fac, p_double)
    b%f3(3,i1-1,i2,i3) = b%f3(3,i1-1,i2,i3) + dt/ldx(1) * real( exp_x1 * exp_x2 * exp_x3 * e2_fac, p_double)
    e%f3(2,i1,i2,i3) = e%f3(2,i1,i2,i3) + dt/ldx(1) * real( exp_x1 * exp_x2 * exp_x3 * b3_fac, p_double)
    e%f3(3,i1,i2,i3) = e%f3(3,i1,i2,i3) - dt/ldx(1) * real( exp_x1 * exp_x2 * exp_x3 * b2_fac, p_double)
   enddo
   enddo
  endif

  ! x1 upper wall
  ! TODO: not working perfectly fine
  if ( this%wall_injection(1,2) .and. my_ngp( no_co, 1 ) == nx( no_co, 1 ) ) then
   ! iterate for every space point
   i1 = b%nx_(1) - 0
   exp_x1 = exp( im_unit * k1 * (border_x1_down-1+i1)*ldx(1) )
   exp_x2 = exp( im_unit * k2 * (border_x2_down-1)*ldx(2) )
   do i2 = border_x2_down, border_x2_up
    ! calculate grid positions and exponentials for x2
    exp_x2 = exp_x2 * exp_x20
   exp_x3 = exp( im_unit * k3 * (border_x3_down-1)*ldx(3) )
   do i3 = border_x3_down, border_x3_up
    ! calculate grid positions and exponentials for x3
    exp_x3 = exp_x3 * exp_x30
   ! --------------------------------------------------
    b%f3(2,i1,i2,i3) = b%f3(2,i1,i2,i3) + dt/ldx(1) * real( exp_x1 * exp_x2 * exp_x3 * e3_fac, p_double)
    b%f3(3,i1,i2,i3) = b%f3(3,i1,i2,i3) - dt/ldx(1) * real( exp_x1 * exp_x2 * exp_x3 * e2_fac, p_double)
    e%f3(2,i1,i2,i3) = e%f3(2,i1,i2,i3) - dt/ldx(1) * real( exp_x1 * exp_x2 * exp_x3 * b3_fac, p_double)
    e%f3(3,i1,i2,i3) = e%f3(3,i1,i2,i3) + dt/ldx(1) * real( exp_x1 * exp_x2 * exp_x3 * b2_fac, p_double)
   enddo
   enddo
  endif

  ! x2 down wall
  if ( this%wall_injection(2,1) .and. my_ngp( no_co, 2 ) == 1 ) then
   ! iterate for every space point
   i2 = 2
   exp_x2 = exp( im_unit * k2 * (border_x2_down-1+i2)*ldx(2) )
   exp_x1 = exp( im_unit * k1 * (border_x1_down-1)*ldx(1) )
   do i1 = border_x1_down, border_x1_up
    ! calculate grid positions and exponentials for x1
    exp_x1 = exp_x1 * exp_x10
   exp_x3 = exp( im_unit * k3 * (border_x3_down-1)*ldx(3) )
   do i3 = border_x3_down, border_x3_up
    ! calculate grid positions and exponentials for x3
    exp_x3 = exp_x3 * exp_x30
   ! --------------------------------------------------
    b%f3(1,i1,i2-1,i3) = b%f3(1,i1,i2-1,i3) + dt/ldx(2) * real( exp_x1 * exp_x2 * exp_x3 * e3_fac, p_double)
    b%f3(3,i1,i2-1,i3) = b%f3(3,i1,i2-1,i3) - dt/ldx(2) * real( exp_x1 * exp_x2 * exp_x3 * e1_fac, p_double)
    e%f3(1,i1,i2,i3) = e%f3(1,i1,i2,i3) - dt/ldx(2) * real( exp_x1 * exp_x2 * exp_x3 * b3_fac, p_double)
    e%f3(3,i1,i2,i3) = e%f3(3,i1,i2,i3) + dt/ldx(2) * real( exp_x1 * exp_x2 * exp_x3 * b1_fac, p_double)
   enddo
   enddo
  endif

  ! x2 upper wall
  ! TODO: not working perfectly fine
  if ( this%wall_injection(2,2) .and. my_ngp( no_co, 2 ) == nx( no_co, 2 ) ) then
   ! iterate for every space point
   i2 = b%nx_(2) - 0
   exp_x2 = exp( im_unit * k2 * (border_x2_down-1+i2)*ldx(2) )
   exp_x1 = exp( im_unit * k1 * (border_x1_down-1)*ldx(1) )
   do i1 = border_x1_down, border_x1_up
    ! calculate grid positions and exponentials for x1
    exp_x1 = exp_x1 * exp_x10
   exp_x3 = exp( im_unit * k3 * (border_x3_down-1)*ldx(3) )
   do i3 = border_x3_down, border_x3_up
    ! calculate grid positions and exponentials for x3
    exp_x3 = exp_x3 * exp_x30
   ! --------------------------------------------------
    b%f3(1,i1,i2,i3) = b%f3(1,i1,i2,i3) - dt/ldx(2) * real( exp_x1 * exp_x2 * exp_x3 * e3_fac, p_double)
    b%f3(3,i1,i2,i3) = b%f3(3,i1,i2,i3) + dt/ldx(2) * real( exp_x1 * exp_x2 * exp_x3 * e1_fac, p_double)
    e%f3(1,i1,i2,i3) = e%f3(1,i1,i2,i3) + dt/ldx(2) * real( exp_x1 * exp_x2 * exp_x3 * b3_fac, p_double)
    e%f3(3,i1,i2,i3) = e%f3(3,i1,i2,i3) - dt/ldx(2) * real( exp_x1 * exp_x2 * exp_x3 * b1_fac, p_double)
   enddo
   enddo
  endif

  ! x3 down wall
  if ( this%wall_injection(3,1) .and. my_ngp( no_co, 3 ) == 1 ) then
   ! iterate for every space point
   i3 = 2
   exp_x3 = exp( im_unit * k3 * (border_x3_down-1+i3)*ldx(3) )
   exp_x1 = exp( im_unit * k1 * (border_x1_down-1)*ldx(1) )
   do i1 = border_x1_down, border_x1_up
    ! calculate grid positions and exponentials for x1
    exp_x1 = exp_x1 * exp_x10
   exp_x2 = exp( im_unit * k2 * (border_x2_down-1)*ldx(2) )
   do i2 = border_x2_down, border_x2_up
    ! calculate grid positions and exponentials for x2
    exp_x2 = exp_x2 * exp_x20
   ! --------------------------------------------------
    b%f3(1,i1,i2,i3-1) = b%f3(1,i1,i2,i3-1) - dt/ldx(3) * real( exp_x1 * exp_x2 * exp_x3 * e2_fac, p_double)
    b%f3(2,i1,i2,i3-1) = b%f3(2,i1,i2,i3-1) + dt/ldx(3) * real( exp_x1 * exp_x2 * exp_x3 * e1_fac, p_double)
    e%f3(1,i1,i2,i3) = e%f3(1,i1,i2,i3) + dt/ldx(3) * real( exp_x1 * exp_x2 * exp_x3 * b2_fac, p_double)
    e%f3(2,i1,i2,i3) = e%f3(2,i1,i2,i3) - dt/ldx(3) * real( exp_x1 * exp_x2 * exp_x3 * b1_fac, p_double)
   enddo
   enddo
  endif

  ! x3 upper wall
  ! TODO: not working perfectly fine
  if ( this%wall_injection(3,2) .and. my_ngp( no_co, 3 ) == nx( no_co, 3 ) ) then
   ! iterate for every space point
   i3 = b%nx_(3) - 0
   exp_x3 = exp( im_unit * k3 * (border_x3_down-1+i3)*ldx(3) )
   exp_x1 = exp( im_unit * k1 * (border_x1_down-1)*ldx(1) )
   do i1 = border_x1_down, border_x1_up
    ! calculate grid positions and exponentials for x1
    exp_x1 = exp_x1 * exp_x10
   exp_x2 = exp( im_unit * k2 * (border_x2_down-1)*ldx(2) )
   do i2 = border_x2_down, border_x2_up
    ! calculate grid positions and exponentials for x2
    exp_x2 = exp_x2 * exp_x20
   ! --------------------------------------------------
    b%f3(1,i1,i2,i3) = b%f3(1,i1,i2,i3) + dt/ldx(3) * real( exp_x1 * exp_x2 * exp_x3 * e2_fac, p_double)
    b%f3(2,i1,i2,i3) = b%f3(2,i1,i2,i3) - dt/ldx(3) * real( exp_x1 * exp_x2 * exp_x3 * e1_fac, p_double)
    e%f3(1,i1,i2,i3) = e%f3(1,i1,i2,i3) - dt/ldx(3) * real( exp_x1 * exp_x2 * exp_x3 * b2_fac, p_double)
    e%f3(2,i1,i2,i3) = e%f3(2,i1,i2,i3) + dt/ldx(3) * real( exp_x1 * exp_x2 * exp_x3 * b1_fac, p_double)
   enddo
   enddo
  endif


 enddo

end subroutine launch_zpulse_planewaves_wall_3d
!-----------------------------------------------------------------------------------------



end module m_zpulse_planewaves
