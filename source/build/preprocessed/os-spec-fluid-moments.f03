# 1 "overdense/os-spec-fluid-moments.f03"
# 1 "<built-in>" 1
# 1 "<built-in>" 3
# 467 "<built-in>" 3
# 1 "<command line>" 1
# 1 "<built-in>" 2
# 1 "overdense/os-spec-fluid-moments.f03" 2
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
# 2 "overdense/os-spec-fluid-moments.f03" 2
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
# 3 "overdense/os-spec-fluid-moments.f03" 2


module m_spec_fluid_moments

use m_parameters
use m_node_conf
use m_grid
use m_grid_define
use m_species_define
use m_vdf
use m_vdf_define

implicit none

private

!-----------------------------------------------------------------------------------------
! moments 'class'
!-----------------------------------------------------------------------------------------
type t_fluid_moments

  type(t_vdf) :: count
  integer :: n_count = -1

  type(t_vdf) :: density
  integer :: n_density = -1

  type(t_vdf) :: current
  integer :: n_current = -1

  !type(t_vdf) :: fluid_velocity

  type(t_vdf) :: abs_momentum
  integer :: n_abs_momentum = -1

  !type(t_vdf) :: abs_proper_velocity

  type(t_vdf) :: fourth_root_abs_momentum
  integer :: n_fourth_root_abs_momentum = -1

  !type(t_vdf) :: fourth_root_abs_proper_velocity

  type(t_vdf) :: pressure
  integer :: n_pressure = -1

  type(t_vdf) :: fourth_root_pressure
  integer :: n_fourth_root_pressure = -1

  type(t_vdf) :: rest_frame_pressure
  integer :: n_rest_frame_pressure = -1

  !type(t_vdf) :: temperature

  type(t_vdf) :: pressure_tensor
  integer :: n_pressure_tensor = -1

  type(t_vdf) :: heat_flux
  integer :: n_heat_flux = -1

  ! hard-coded to 0 for now, could be higher with future implementations
  integer :: interpolation = 0

end type t_fluid_moments

interface setup
  module procedure setup_fluid_moment
end interface setup

interface cleanup
  module procedure cleanup_fluid_moments
end interface cleanup

interface calculate
  module procedure calculate_fluid_moment
end interface calculate

interface read_fluid_moment
  module procedure read_fluid_moment_single
  module procedure read_fluid_moment_single_scalar
end interface read_fluid_moment

public :: t_fluid_moments
public :: setup, cleanup, calculate, read_fluid_moment

contains

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
recursive subroutine setup_fluid_moment(this, moment_name, grid, no_co, dx)

  type(t_fluid_moments), intent(inout) :: this
  character(len=*), intent(in) :: moment_name
  class(t_grid), intent(in) :: grid
  class( t_node_conf ), intent(in) :: no_co
  real(p_double), dimension(:), intent(in) :: dx

  integer, dimension(2,p_x_dim) :: gc_num

  gc_num = this%interpolation ! this isn't actually quite right (except for zeroth level interpolation)

  select case (moment_name)

  case('density')

    if( this%density%f_dim() < 1 ) &
      call this%density%new( p_x_dim, 1, grid%my_nx(3,:), gc_num, dx, .true. )

  case('current')

    if( this%current%f_dim() < 1 ) &
      call this%current%new( p_x_dim, p_p_dim, grid%my_nx(3,:), gc_num, dx, .true. )

  case('fluid_velocity')

    if( this%density%f_dim() < 1 ) &
      call setup(this, 'density', grid, no_co, dx)

    if( this%current%f_dim() < 1 ) &
      call setup(this, 'current', grid, no_co, dx)

  case('abs_momentum')

    if( this%abs_momentum%f_dim() < 1 ) &
      call this%abs_momentum%new( p_x_dim, p_p_dim, grid%my_nx(3,:), gc_num, dx, .true. )

  case('abs_proper_velocity')

    if( this%density%f_dim() < 1 ) &
      call setup(this, 'density', grid, no_co, dx)
    if( this%abs_momentum%f_dim() < 1 ) &
      call setup(this, 'abs_momentum', grid, no_co, dx)

  case('fourth_root_abs_momentum')

    if( this%abs_momentum%f_dim() < 1 ) &
      call this%fourth_root_abs_momentum%new( p_x_dim, p_p_dim, grid%my_nx(3,:), gc_num, dx, .true. )

  case('fourth_root_abs_proper_velocity')

    if( this%density%f_dim() < 1 ) &
      call setup(this, 'density', grid, no_co, dx)
    if( this%fourth_root_abs_momentum%f_dim() < 1 ) &
      call setup(this, 'fourth_root_abs_momentum', grid, no_co, dx)

  case('pressure')

    if( this%pressure%f_dim() < 1 ) &
      call this%pressure%new( p_x_dim, 1, grid%my_nx(3,:), gc_num, dx, .true. )

  case('fourth_root_pressure')

    if( this%fourth_root_pressure%f_dim() < 1 ) &
      call this%fourth_root_pressure%new( p_x_dim, 1, grid%my_nx(3,:), gc_num, dx, .true.)

  case('rest_frame_pressure')

    call setup(this, 'fluid_velocity', grid, no_co, dx)
    if( this%rest_frame_pressure%f_dim() < 1 ) &
      call this%rest_frame_pressure%new( p_x_dim, 1, grid%my_nx(3,:), gc_num, dx, .true. )

  case('temperature')

    if( this%density%f_dim() < 1 ) &
      call setup(this, 'density', grid, no_co, dx)
    if( this%pressure%f_dim() < 1 ) &
      call setup(this, 'pressure', grid, no_co, dx)

  case('fourth_root_temperature')

    if( this%density%f_dim() < 1 ) &
      call setup(this, 'density', grid, no_co, dx)
    if( this%fourth_root_pressure%f_dim() < 1 ) &
      call setup(this, 'fourth_root_pressure', grid, no_co, dx)

  case('pressure_tensor')

    if( this%pressure_tensor%f_dim() < 1 ) &
      call this%pressure_tensor%new( p_x_dim, p_p_dim*p_p_dim, grid%my_nx(3,:), gc_num, &
                                      dx, .true. )

  case('heat_flux')

    if( this%heat_flux%f_dim() < 1 ) &
      call this%heat_flux%new( p_x_dim, p_p_dim, grid%my_nx(3,:), gc_num, dx, .true. )

  case default

    write(err_buf__,*) 'moment ', moment_name, ' not implemented';call err__("overdense/os-spec-fluid-moments.f03",191)
    write(err_buf__,*) 'aborting setup_fluid_moment';call err__("overdense/os-spec-fluid-moments.f03",192)
    call abort_program( p_err_notimplemented )

  end select


end subroutine setup_fluid_moment
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
recursive subroutine calculate_fluid_moment(this, species, moment_name, n)

  ! Dummy variables
  type(t_fluid_moments), intent(inout) :: this
  class(t_species), intent(in) :: species
  character(len=*), intent(in) :: moment_name
  integer, intent(in) :: n

  ! Local variables
  integer :: i, i0, i1
  real(p_k_part), dimension(p_cache_size) :: particle_moments
  integer :: component

  ! Executable statements

  select case (moment_name)

  case('density')

    if( this%n_density .eq. n) return

    this%density = 0.0_p_k_part
    ! operator overloading, nice
    do i = 1, species%num_par, p_cache_size

      i0 = i
      i1 = i + p_cache_size - 1
      if( i1 .gt. species%num_par) i1 = species%num_par
      call get_particle_moments( species, i0, i1, particle_moments, 'density' )
      call deposit_particle_moments( species, this%density, 1, particle_moments, i0, i1, &
                                     this%interpolation )

    enddo

    this%n_density = n

  case('current')

    if( this%n_current .eq. n) return

    this%current = 0.0_p_k_part

    do component = 1, p_p_dim

      do i = 1, species%num_par, p_cache_size

        i0 = i
        i1 = i + p_cache_size - 1
        if( i1 .gt. species%num_par) i1 = species%num_par
        call get_particle_moments(species, i0, i1, particle_moments, 'current', component)
        call deposit_particle_moments(species, this%current, component, particle_moments,&
                                      i0, i1, this%interpolation )

      enddo

    enddo

    this%n_current = n

  case('fluid_velocity')

    call calculate(this, species, 'density', n)
    call calculate(this, species, 'current', n)

  case('abs_momentum')

    if( this%n_abs_momentum .eq. n) return

    this%abs_momentum = 0.0_p_k_part

    do component = 1, p_p_dim

      do i = 1, species%num_par, p_cache_size

        i0 = i
        i1 = i + p_cache_size - 1
        if( i1 .gt. species%num_par) i1 = species%num_par
        call get_particle_moments(species, i0, i1, particle_moments, 'abs_momentum', component)
        call deposit_particle_moments(species, this%abs_momentum, component, particle_moments,&
                                      i0, i1, this%interpolation )

      enddo

    enddo

    this%n_abs_momentum = n

  case('abs_proper_velocity')

    call calculate(this, species, 'density', n)
    call calculate(this, species, 'abs_momentum', n)

  case('fourth_root_abs_momentum')

    if( this%n_fourth_root_abs_momentum .eq. n) return

    this%fourth_root_abs_momentum = 0.0_p_k_part

    do component = 1, p_p_dim

      do i = 1, species%num_par, p_cache_size

        i0 = i
        i1 = i + p_cache_size - 1
        if( i1 .gt. species%num_par) i1 = species%num_par
        call get_particle_moments(species, i0, i1, particle_moments, 'fourth_root_abs_momentum', component)
        call deposit_particle_moments(species, this%fourth_root_abs_momentum, component, particle_moments,&
                                      i0, i1, this%interpolation )

      enddo

    enddo

    this%n_fourth_root_abs_momentum = n

  case('fourth_root_abs_proper_velocity')

    call calculate(this, species, 'density', n)
    call calculate(this, species, 'fourth_root_abs_momentum', n)

  case('pressure')

    if( this%n_pressure .eq. n) return

    this%pressure = 0.0_p_k_part

    do i = 1, species%num_par, p_cache_size

      i0 = i
      i1 = i + p_cache_size - 1
      if( i1 .gt. species%num_par) i1 = species%num_par
      call get_particle_moments( species, i0, i1, particle_moments, 'pressure' )
      call deposit_particle_moments( species, this%pressure, 1, particle_moments, i0, i1,&
                                     this%interpolation )

    enddo

    this%n_pressure = n

  case('fourth_root_pressure')

    if( this%n_fourth_root_pressure .eq. n) return

    this%fourth_root_pressure = 0.0_p_k_part

    do i = 1, species%num_par, p_cache_size

      i0 = i
      i1 = i + p_cache_size - 1
      if( i1 .gt. species%num_par) i1 = species%num_par
      call get_particle_moments(species, i0, i1, particle_moments, 'fourth_root_pressure')
      call deposit_particle_moments( species, this%fourth_root_pressure, 1, &
                                      particle_moments, i0, i1, this%interpolation )

    enddo

    this%n_fourth_root_pressure = n

  case('rest_frame_pressure')

    call calculate(this, species, 'fluid_velocity', n)
    if( this%n_rest_frame_pressure .eq. n) return

    continue

  case('temperature')

    call calculate(this, species, 'density', n)
    call calculate(this, species, 'pressure', n)

  case('fourth_root_temperature')

    call calculate(this, species, 'density', n)
    call calculate(this, species, 'fourth_root_pressure', n)

  case('pressure_tensor')

    if( this%n_pressure_tensor .eq. n) return

    continue

  case('heat_flux')

    if( this%n_heat_flux .eq. n) return

    continue

  case default

    write(err_buf__,*) 'moment ', moment_name, ' not implemented';call err__("overdense/os-spec-fluid-moments.f03",393)
    write(err_buf__,*) 'aborting calculate_fluid_moment';call err__("overdense/os-spec-fluid-moments.f03",394)
    call abort_program( p_err_notimplemented )

  end select

end subroutine calculate_fluid_moment
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
subroutine cleanup_fluid_moments( this )

  type( t_fluid_moments), intent( inout ) :: this

  call this % count % cleanup()
  call this % density % cleanup()
  call this % current % cleanup()
  call this % abs_momentum % cleanup()
  call this % fourth_root_abs_momentum % cleanup()
  call this % pressure % cleanup()
  call this % fourth_root_pressure % cleanup()
  call this % rest_frame_pressure % cleanup()
  call this % pressure_tensor % cleanup()
  call this % heat_flux % cleanup()

end subroutine cleanup_fluid_moments
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
subroutine get_particle_moments( species, i0, i1, particle_moments, moment_name, &
                                 component, component2)

  ! Dummy variables
  class(t_species), intent(in) :: species
  integer, intent(in) :: i0, i1
  real(p_k_part), dimension(:), intent(inout) :: particle_moments
  character(len=*), intent(in) :: moment_name
  integer, intent(in), optional :: component
  integer, intent(in), optional :: component2

  ! Local variables
  integer :: npar
  integer :: i,j,k
  real(p_k_part) :: p_squared, gamma

  ! Executable statements

  particle_moments(:) = 0.0_p_k_part

  npar = i1 - i0 + 1

  select case (moment_name)

    case('count')

      particle_moments = 1.0_p_k_part

    case('density')

      do i = 1, npar

        j = i + i0 - 1
        particle_moments(i) = species%q(j) * species%rqm

      enddo

    case('current')

      if( .not. present(component) ) then
        write(err_buf__,*) "Current moment requested but no component specified";call err__("overdense/os-spec-fluid-moments.f03",466)
        call abort_program( -1 )
      endif

      do i = 1, npar

        j = i + i0 - 1
        p_squared = 0.0_p_k_part
        do k = 1, p_p_dim
          p_squared = p_squared + species%p(k,j)**2
        enddo
        gamma = sqrt(p_squared + 1)
        particle_moments(i) = species%q(j) * species%p(component, j) / gamma

      enddo

    case('abs_momentum')

      if( .not. present(component) ) then
        write(err_buf__,*) "abs_momentum moment requested but no component specified";call err__("overdense/os-spec-fluid-moments.f03",485)
        call abort_program( -1 )
      endif

      do i = 1, npar

        j = i + i0 - 1
        particle_moments(i) = species%q(j) * species%rqm * abs(species%p(component, j))

      enddo

    case('fourth_root_abs_momentum')

      if( .not. present(component) ) then
        write(err_buf__,*) "fourth_root_abs_momentum moment requested but no component specified";call err__("overdense/os-spec-fluid-moments.f03",499)
        call abort_program( -1 )
      endif

      do i = 1, npar

        j = i + i0 - 1
        particle_moments(i) = species%q(j) * species%rqm * &
                              sqrt(sqrt(abs(species%p(component, j))))

      enddo

    case('pressure')

      do i = 1, npar

        j = i + i0 - 1
        p_squared = 0.0_p_k_part
        do k = 1, p_p_dim
          p_squared = p_squared + species%p(k,j)**2
        enddo
        gamma = sqrt(p_squared + 1)
        particle_moments(i) = species%q(j) * species%rqm * p_squared/(gamma + 1)

      enddo

    case('fourth_root_pressure')

      do i = 1, npar

        j = i + i0 - 1
        p_squared = 0.0_p_k_part
        do k = 1, p_p_dim
          p_squared = p_squared + species%p(k,j)**2
        enddo
        gamma = sqrt(p_squared + 1)
        particle_moments(i) = species%q(j) * species%rqm * sqrt(sqrt(p_squared/(gamma + 1)))

      enddo

    case('rest_frame_pressure')

      ! I have to think more abou the best way to do this
      continue

    case('pressure_tensor')

      ! if( .not. present(component) ) then
      ! write(err_buf__,*) "Pressure tensor moment requested but no 1st component specified";call err__("overdense/os-spec-fluid-moments.f03",547)
      ! call abort_program( -1 )
      ! endif

      ! if( .not. present(component2) ) then
      ! write(err_buf__,*) "Pressure tensor moment requested but no 2nd component specified";call err__("overdense/os-spec-fluid-moments.f03",552)
      ! call abort_program( -1 )
      ! endif

      ! do i = 1, npar

      ! j = i + i0 - 1
      ! p_squared = 0.0_p_k_part
      ! do k = 1, p_p_dim
      ! p_squared = p_squared + species%p(k,j)**2
      ! enddo
      ! gamma = sqrt(p_squared + 1)
      ! particle_moments(i) = species%q(j) * species%rqm * species%p(component, j) &
      ! * species%p(component2, j) / ( gamma + 1 )

      ! enddo
      ! I'm not 100% sure that's relativisticly correct (or correct at all)

      continue

    case('heat_flux')

      if( .not. present(component) ) then
        write(err_buf__,*) "Heat flux moment requested but no component specified";call err__("overdense/os-spec-fluid-moments.f03",575)
        call abort_program( -1 )
      endif

      do i = 1, npar

        j = i + i0 - 1
        p_squared = 0.0_p_k_part
        do k = 1, p_p_dim
          p_squared = p_squared + species%p(k,j)**2
        enddo
        gamma = sqrt(p_squared + 1)
        particle_moments(i) = species%q(j) * species%rqm * p_squared/(gamma + 1) &
                               * species%p(component, j) / gamma

      enddo

    case default

      write(err_buf__,*) 'moment ', moment_name, ' not implemented';call err__("overdense/os-spec-fluid-moments.f03",594)
      write(err_buf__,*) 'aborting get_particle_moments';call err__("overdense/os-spec-fluid-moments.f03",595)
      call abort_program( p_err_notimplemented )

  end select

end subroutine get_particle_moments
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
subroutine deposit_particle_moments( species, vdf, component, particle_moments, i0, i1, &
                                     interpolation )

  ! Dummy variables
  class( t_species ), intent(in) :: species
  type( t_vdf ), intent(inout) :: vdf
  integer, intent(in) :: component
  real(p_k_part), dimension(:), intent(in) :: particle_moments
  integer, intent(in) :: i0, i1
  integer, intent(in) :: interpolation

  ! Local variables
  integer :: npar
  integer :: i, j
  integer :: ix1, ix2, ix3

  ! Executable statements

  npar = i1 - i0 + 1

  select case (interpolation)

  case(0)

    select case (p_x_dim)

    case(1)

      do i = 1, npar

        j = i + i0 - 1
        ix1 = species%ix(1,j)
        vdf%f1(component, ix1 ) = vdf%f1(component, ix1 ) &
                                  + particle_moments(i)

      enddo

    case(2)

      do i = 1, npar

        j = i + i0 - 1
        ix1 = species%ix(1,j)
        ix2 = species%ix(2,j)
        vdf%f2(component, ix1, ix2 ) = vdf%f2(component, ix1, ix2 ) &
                                       + particle_moments(i)

      enddo

    case(3)

      do i = 1, npar

        j = i + i0 - 1
        ix1 = species%ix(1,j)
        ix2 = species%ix(2,j)
        ix3 = species%ix(3,j)
        vdf%f3(component, ix1, ix2, ix3 ) = vdf%f3(component, ix1, ix2, ix3 ) &
                                            + particle_moments(i)

      enddo

    end select

  case default

      write(err_buf__,*) interpolation, 'th level fluid moment interpolation not implemented';call err__("overdense/os-spec-fluid-moments.f03",672)
      write(err_buf__,*) 'aborting deposit_particle_moments';call err__("overdense/os-spec-fluid-moments.f03",673)
      call abort_program( p_err_notimplemented )

  end select

end subroutine deposit_particle_moments
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
subroutine read_fluid_moment_single_scalar( this, species, index, moment_name, result )

  ! Dummy variables
  type( t_fluid_moments ), intent(in) :: this
  class( t_species ), intent(in) :: species
  integer, intent(in) :: index
  character(len=*), intent(in) :: moment_name
  real(p_k_part), intent(out) :: result

  ! Local variables
  real(p_k_part), dimension(1) :: result_array

  ! Executable statements

  call read_fluid_moment_single( this, species, index, moment_name, result_array )
  result = result_array(1)

end subroutine read_fluid_moment_single_scalar
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
subroutine read_fluid_moment_single( this, species, index, moment_name, result )

  ! Dummy variables
  type( t_fluid_moments ), intent(in) :: this
  class( t_species ), intent(in) :: species
  integer, intent(in) :: index
  character(len=*), intent(in) :: moment_name
  real(p_k_part), dimension(:), intent(out) :: result

  ! Local variables
  integer :: i
  real(p_k_part) :: density, pressure, fourth_root_pressure
  real(p_k_part), dimension(p_p_dim) :: abs_momentum, fourth_root_abs_momentum

  ! Executable statements

  select case (moment_name)

    case('density')

      call interpolate_moment_single( this%density, 1, species, index, this%interpolation, result(1) )

    case('current')

      do i = 1, p_p_dim
        call interpolate_moment_single( this%current, i, species, index, this%interpolation, result(i) )
      enddo

    case('abs_momentum')

      do i = 1, p_p_dim
        call interpolate_moment_single( this%abs_momentum, i, species, index, this%interpolation, result(i) )
      enddo

    case('abs_proper_velocity')

      call interpolate_moment_single( this%density, 1, species, index, this%interpolation, density )
      do i = 1, p_p_dim
        call interpolate_moment_single( this%abs_momentum, i, species, index, this%interpolation, abs_momentum(i) )
      enddo
      if ( density > 0.0_p_k_part ) then
        do i = 1, p_p_dim
          result(i) = abs_momentum(i) / density
        enddo
      else
        result = 0.0_p_k_part
      endif

    case('fourth_root_abs_momentum')

      do i = 1, p_p_dim
        call interpolate_moment_single( this%fourth_root_abs_momentum, i, species, index, this%interpolation, result(i) )
      enddo

    case('fourth_root_abs_proper_velocity')

      call interpolate_moment_single( this%density, 1, species, index, this%interpolation, density )
      do i = 1, p_p_dim
        call interpolate_moment_single( this%fourth_root_abs_momentum, i, species, index, this%interpolation, fourth_root_abs_momentum(i) )
      enddo
      if ( density > 0.0_p_k_part ) then
        do i = 1, p_p_dim
          result(i) = fourth_root_abs_momentum(i) / density
        enddo
      else
        result = 0.0_p_k_part
      endif

    case('pressure')

      call interpolate_moment_single( this%pressure, 1, species, index, this%interpolation, result(1) )

    case('fourth_root_pressure')

      call interpolate_moment_single( this%fourth_root_pressure, 1, species, index, this%interpolation, result(1) )

    case('temperature')

      call interpolate_moment_single( this%density, 1, species, index, this%interpolation, density )
      call interpolate_moment_single( this%pressure, 1, species, index, this%interpolation, pressure )

      if( density .gt. 0.0_p_k_part) then
        result(1) = pressure / density
      else
        result(1) = 0.0_p_k_part
      endif

    case('fourth_root_temperature')

      call interpolate_moment_single( this%density, 1, species, index, this%interpolation, density )
      call interpolate_moment_single( this%fourth_root_pressure, 1, species, index, this%interpolation, fourth_root_pressure )

      if( density .gt. 0.0_p_k_part) then
        result(1) = fourth_root_pressure / density
      else
        result(1) = 0.0_p_k_part
      endif

    case default

      write(err_buf__,*) 'moment', moment_name, 'not implemented';call err__("overdense/os-spec-fluid-moments.f03",807)
      write(err_buf__,*) 'aborting read_fluid_moment_single';call err__("overdense/os-spec-fluid-moments.f03",808)
      call abort_program( p_err_notimplemented )

  end select

end subroutine read_fluid_moment_single
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
subroutine interpolate_moment_single( vdf, component, species, index, interpolation, result )

  ! Dummy variables
  type(t_vdf), intent(in) :: vdf
  integer, intent(in) :: component
  class(t_species), intent(in) :: species
  integer, intent(in) :: index
  integer, intent(in) :: interpolation
  real(p_k_part), intent(out) :: result

  ! Local variables
  integer :: ix1, ix2, ix3

  ! Executable statements

  select case (interpolation)

  case(0)

    select case (p_x_dim)

    case(1)

      ix1 = species%ix( 1, index )
      result = vdf%f1( component, ix1 )

    case(2)

      ix1 = species%ix( 1, index )
      ix2 = species%ix( 2, index )
      result = vdf%f2( component, ix1, ix2)

    case(3)

      ix1 = species%ix( 1, index )
      ix2 = species%ix( 2, index )
      ix3 = species%ix( 3, index )
      result = vdf%f3( component, ix1, ix2, ix3)

    end select

  case default

    write(err_buf__,*) interpolation, 'th level fluid moment interpolation not implemented';call err__("overdense/os-spec-fluid-moments.f03",862)
    write(err_buf__,*) 'aborting interpolate_moment_single';call err__("overdense/os-spec-fluid-moments.f03",863)
    call abort_program( p_err_notimplemented )

  end select

end subroutine interpolate_moment_single
!-------------------------------------------------------------------------------

end module m_spec_fluid_moments
