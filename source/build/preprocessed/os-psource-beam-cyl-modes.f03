# 1 "cyl_modes/psource/os-psource-beam-cyl-modes.f03"
# 1 "<built-in>" 1
# 1 "<built-in>" 3
# 467 "<built-in>" 3
# 1 "<command line>" 1
# 1 "<built-in>" 2
# 1 "cyl_modes/psource/os-psource-beam-cyl-modes.f03" 2
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
# 2 "cyl_modes/psource/os-psource-beam-cyl-modes.f03" 2
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
# 3 "cyl_modes/psource/os-psource-beam-cyl-modes.f03" 2

module m_psource_beam_cyl_modes

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
# 7 "cyl_modes/psource/os-psource-beam-cyl-modes.f03" 2

use m_species_define
use m_psource_std_cyl_modes
use m_parameters
use m_fparser
use m_psource_std
use m_current_define
use m_node_conf

implicit none

private

interface transform_focal
  module procedure transform_focal_beam_cyl_modes
end interface

integer, parameter :: p_gaussian = 2 ! gaussian

type, extends(t_psource_std_cyl_modes) :: t_psource_beam_cyl_modes

  ! Focal plane distance after acceleration
  real(p_k_part), dimension(p_p_dim) :: focal_dist

  ! Params to rotate p2x2 and px3x3 phase space
  real(p_k_part), dimension(p_p_dim) :: alpha

  ! Beam uth parameter at Focal Plane
  real(p_k_part), dimension(p_p_dim) :: uth

  real(p_k_part) :: gamma

contains

  ! Methods from abstract class
  procedure :: read_input => read_input_beam_cyl_modes
  procedure :: inject => inject_beam_cyl_modes

end type t_psource_beam_cyl_modes


public :: t_psource_beam_cyl_modes


contains

!-------------------------------------------------------------------------------
! Read information from input file
!-------------------------------------------------------------------------------
subroutine read_input_beam_cyl_modes( this, input_file, coordinates )
  use m_input_file

  implicit none

  class( t_psource_beam_cyl_modes ), intent(inout) :: this
  class( t_input_file ), intent(inout) :: input_file
  integer, intent(in) :: coordinates


  ! minimum density to inject particles
  real(p_k_part) :: den_min

  ! global density (multiplies the selected density profile)
  real(p_k_part) :: density

  ! ratio of sigma_y/sigma_x (only use in quasi-3D)
  real(p_k_part) :: aspect_ratio

  ! variables needed for Gaussian profiles
  real(p_k_part), dimension(p_x_dim) :: gauss_center
  real(p_k_part), dimension(p_x_dim) :: gauss_sigma
  integer, dimension( p_x_dim ) :: gauss_n_sigma
  real(p_k_part), dimension(2,p_x_dim) :: gauss_range

  ! This is only used for t_psource_constq but is included here for simplicity
  !integer, dimension(p_max_dim) :: sample_rate

  real(p_double), dimension(p_p_dim) :: focal_dist, alpha
  real(p_double), dimension(p_p_dim) :: uth
  real(p_double) :: gamma


  ! namelists of input variables
  namelist /nl_profile/ den_min, density, gauss_n_sigma, gauss_range, &
  gauss_center, gauss_sigma, focal_dist, alpha, &
  uth, gamma, aspect_ratio

  integer :: j
  integer :: ierr

  aspect_ratio = 1.0_p_k_part
  den_min = 0
  density = 1.0
  gauss_center = 0.0_p_k_part
  gauss_sigma = huge( 1.0_p_k_part )
  gauss_n_sigma = 0
  gauss_range(p_lower,:) = -gauss_sigma
  gauss_range(p_upper,:) = gauss_sigma

  focal_dist = 0.0
  alpha = 0.0
  uth = 0
  gamma = 1

  ! Get namelist text from input file
  call get_namelist( input_file, "nl_profile", ierr )
  if ( ierr == 0 ) then
    read (input_file%nml_text, nml = nl_profile, iostat = ierr)
    if (ierr /= 0) then
      if ( mpi_node() == 0 ) then
        write(0,*) ""
        write(0,*) "   Error reading profile information"
        write(0,*) "   aborting..."
      endif
      stop
    endif
  else
    if (mpi_node()==0) print *,"   - profile parameters missing, using uniform density"
  endif
  ! minimum density for injection
  this%den_min = den_min
  this%aspect_ratio = aspect_ratio
  this%density = density
  this % if_mult = .true.
  do j = 1, p_x_dim
    this%type(j) = p_gaussian

    this%gauss_center(j) = gauss_center(j)
    this%gauss_sigma(j) = abs( gauss_sigma(j) )
    if ( gauss_n_sigma(j) > 0 ) then
      this%gauss_range(p_lower,j) = gauss_center(j)-abs(gauss_n_sigma(j)*gauss_sigma(j))
      this%gauss_range(p_upper,j) = gauss_center(j)+abs(gauss_n_sigma(j)*gauss_sigma(j))
    else
      this%gauss_range(:,j) = gauss_range(:,j)
    endif
  enddo

  this%uth = uth
  this%alpha = alpha
  this%focal_dist = focal_dist
  this%gamma = gamma
  call transform_focal(this)
end subroutine read_input_beam_cyl_modes
!---------------------------------------------------------------------------------------------------




!-----------------------------------------------------------------------------------------
! Transform the particle beam (3D/Quasi-3D/rz version)
!-----------------------------------------------------------------------------------------
subroutine transform_focal_beam_cyl_modes( this )

  use m_psource_std
  use m_psource_std_cyl_modes
  implicit none

  class( t_psource_beam_cyl_modes ), intent( inout ) :: this

  real(p_double) :: sigma_mul, p1, dens_mul, temp
  integer :: j

  p1 = sqrt( this % gamma ** 2 - 1 )
  dens_mul = 1.0_p_k_part

  do j = 2, p_p_dim
    ! focal plane distance defined after n_acceleration
    this%focal_dist(j) = this%focal_dist(j)
    !transform params

    if(j == p_p_dim) then
      temp = (this%gauss_sigma(p_p_dim-1)*this%aspect_ratio/(sigma_mul ))
      sigma_mul = sqrt(1 + this%focal_dist(j)**2 *this%uth(j)**2 &
        /(p1**2 * (this%gauss_sigma(p_p_dim-1)*this%aspect_ratio/(sigma_mul) )**2 ))
    else
      sigma_mul = sqrt(1 + this%focal_dist(j)**2 *this%uth(j)**2 /(p1**2 * this%gauss_sigma(j)**2 ))
      this%gauss_sigma(j) = sigma_mul * this%gauss_sigma(j)
      this%gauss_range(:,j) = sigma_mul * this%gauss_range(:,j)
    endif
    dens_mul = dens_mul * sigma_mul
    this%uth(j) = this%uth(j)/sigma_mul
    if(this%focal_dist(j) /= 0.0) then
      this%alpha(j) = p1 / this%focal_dist(j) * (1.0_p_k_part- 1.0_p_k_part/sigma_mul**2)
    endif
  enddo
  ! transform density
  this%density = this%density /dens_mul
  this%density = this%density * this%aspect_ratio

  !recalculate aspect ratio
  this%aspect_ratio = temp * sigma_mul/this%gauss_sigma(p_p_dim-1)

  !transform density back with new aspect ratio
  this%density = this%density /this%aspect_ratio

  !adjust range if needed
  if(this%aspect_ratio > 1.0_p_double) then
    this%gauss_range(:,p_p_dim-1) = this%gauss_range(:,p_p_dim-1) * this%aspect_ratio
  endif


end subroutine transform_focal_beam_cyl_modes
!-----------------------------------------------------------------------------------------



function inject_beam_cyl_modes( this, species, ig_xbnd_inj, jay, no_co, bnd_cross, &
                                node_cross, send_msg, recv_msg ) result(num_inj)

  use m_species_udist
  use m_random
  implicit none

  class(t_psource_beam_cyl_modes), intent(inout) :: this
  class(t_species), intent(inout), target :: species
  integer, dimension(:, :), intent(in) :: ig_xbnd_inj
  class( t_current ), intent(inout) :: jay
  class( t_node_conf ), intent(in) :: no_co
  type( t_part_idx ), dimension(2), intent(inout) :: bnd_cross
  type( t_part_idx ), intent(inout) :: node_cross
  type( t_spec_msg ), dimension(2), intent(inout) :: send_msg, recv_msg

  integer :: num_inj, i, j
  real(p_double), dimension( p_p_dim + 1 ) :: x
  real(p_double) :: u1

  ! Call superclass to inject particles
  num_inj = this % t_psource_std_cyl_modes % inject( species, ig_xbnd_inj, jay, no_co, &
                                              bnd_cross, node_cross, send_msg, recv_msg )
  do i = species%num_par+1, species%num_par + num_inj
    species % p(1,i) = this % uth(1) * rng % genrand_gaussian( )
    species % p(2,i) = this % uth(2) * rng % genrand_gaussian( )
    species % p(3,i) = this % uth(3) * rng % genrand_gaussian( )
  enddo

  if ( any( this % focal_dist /= 0 ) ) then
    do i = species%num_par+1, species%num_par + num_inj
      ! 3D/rz geometry
      call species % get_position( i, x )
      do j = 2, p_p_dim
        species%p(j,i) = species%p(j,i) - x(j+1) * this%alpha(j)
      enddo
    enddo
  endif

  u1 = sqrt( this % gamma ** 2 - 1 )
  species%udist%ufl(1) = u1
  do i = species%num_par+1, species%num_par + num_inj
    species%p(1,i) = species%p(1,i) + u1
  enddo


end function inject_beam_cyl_modes

end module m_psource_beam_cyl_modes
