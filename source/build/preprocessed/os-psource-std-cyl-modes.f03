# 1 "cyl_modes/psource/os-psource-std-cyl-modes.f03"
# 1 "<built-in>" 1
# 1 "<built-in>" 3
# 467 "<built-in>" 3
# 1 "<command line>" 1
# 1 "<built-in>" 2
# 1 "cyl_modes/psource/os-psource-std-cyl-modes.f03" 2
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
# 2 "cyl_modes/psource/os-psource-std-cyl-modes.f03" 2
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
# 3 "cyl_modes/psource/os-psource-std-cyl-modes.f03" 2

module m_psource_std_cyl_modes

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
# 7 "cyl_modes/psource/os-psource-std-cyl-modes.f03" 2

use m_parameters
use m_fparser
use m_current_define
use m_node_conf
use m_psource_std
use m_species_cyl_modes_define
use m_species_define, only: t_species, t_psource, t_part_idx, t_spec_msg
use m_species_udist



use m_input_file

implicit none

private


! density psource types

! Already defined in os-spec-define.f03
! NOTE: p_gaussian is also re-defined in inject_std_cyl_modes at the end of the file
integer, parameter :: p_none = -1 ! no injection (always return 0)
integer, parameter :: p_uniform = 0 ! uniform
integer, parameter :: p_pw_linear = 1 ! piecewise-linear
integer, parameter :: p_gaussian = 2 ! gaussian
! integer, parameter :: p_channel = 3 ! parabolic channel
! integer, parameter :: p_sphere = 4 ! sphere
integer, parameter :: p_func = 5 ! math function

! used in cathode only
! integer, parameter :: p_step = 10
! integer, parameter :: p_thruster = 11

! maximum number of points for the piecewise-linear profiles
integer, parameter :: p_num_x_max = 100

type, extends(t_psource_std) :: t_psource_std_cyl_modes

  ! vars if needed
  ! sigma_y/sigma_x (Quasi-3D parameter)
  real(p_k_part) :: aspect_ratio

contains

  ! Methods from abstract class to override
  procedure :: read_input => read_input_std_cyl_modes
  procedure :: inject => inject_std_cyl_modes

  ! Additional methods
  ! Allow access to density values from outside
  procedure :: get_den_value_cyl => get_den_value_std_cyl_modes


end type t_psource_std_cyl_modes

! Add interface to superclass inject function to avoid type casting
interface
  function inject_std_cyl_modes( this, species, ig_xbnd_inj, jay, no_co, bnd_cross, &
                                 node_cross, send_msg, recv_msg ) result(num_inj)

  import t_psource_std_cyl_modes, t_species, t_current, t_node_conf, t_part_idx, &
         t_spec_msg

  class( t_psource_std_cyl_modes), intent(inout) :: this
  class( t_species ), intent(inout), target :: species
  integer, dimension(:, :), intent(in) :: ig_xbnd_inj
  class( t_current ), intent(inout) :: jay
  class( t_node_conf ), intent(in) :: no_co
  type( t_part_idx ), dimension(2), intent(inout) :: bnd_cross
  type( t_part_idx ), intent(inout) :: node_cross
  type( t_spec_msg ), dimension(2), intent(inout) :: send_msg, recv_msg

  integer :: num_inj

  end function
end interface

public :: t_psource_std_cyl_modes

public :: p_channel, p_sphere

contains

!-------------------------------------------------------------------------------
! Read information from input file
!-------------------------------------------------------------------------------
subroutine read_input_std_cyl_modes( this, input_file, coordinates )

  implicit none

  class( t_psource_std_cyl_modes ), intent(inout) :: this
  class( t_input_file ), intent(inout) :: input_file
  integer, intent(in) :: coordinates

  ! number of points in piecewise-linear profiles
  integer :: num_x
  ! components needed for piecewise linear profiles
  real(p_k_part), dimension(p_num_x_max,p_x_dim) :: x, fx

  ! variable used to choose different types of profiles
  character(20), dimension(p_x_dim) :: profile_type

  ! minimum density to inject particles
  real(p_k_part) :: den_min

  ! ratio of sigma_y/sigma_x (only use in quasi-3D)
  real(p_k_part) :: aspect_ratio

  ! global density (multiplies the selected density profile)
  real(p_k_part) :: density

  ! variables needed for Gaussian profiles
  real(p_k_part), dimension(p_x_dim) :: gauss_center
  real(p_k_part), dimension(p_x_dim) :: gauss_sigma
  integer, dimension( p_x_dim ) :: gauss_n_sigma
  real(p_k_part), dimension(2,p_x_dim) :: gauss_range

  ! variables needed for parabolic channel profiles
  integer :: channel_dir
  real(p_k_part) :: channel_bottom
  real(p_k_part) :: channel_r0
  real(p_k_part) :: channel_depth
  real(p_k_part) :: channel_size
  real(p_k_part), dimension(p_x_dim) :: channel_center
  real(p_k_part) :: channel_wall
  real(p_k_part), dimension(2) :: channel_pos

  ! variables needed for sphere profiles
  real(p_k_part), dimension(p_x_dim) :: sphere_center
  real(p_k_part) :: sphere_radius

  ! variables needed for math func profiles
  character(len = p_max_expr_len) :: math_func_expr

  ! This is only used for t_psource_constq but is included here for simplicity
  integer, dimension(p_max_dim) :: sample_rate

  ! namelists of input variables
  namelist /nl_profile/ profile_type, den_min, &
                        density, num_x, x, fx, gauss_n_sigma, gauss_range, &
                        gauss_center, gauss_sigma, &
                        channel_dir, channel_r0, channel_depth, channel_size, &
                        channel_center, channel_wall, channel_pos, channel_bottom, &
                        sphere_center, sphere_radius, math_func_expr, &
                        sample_rate, aspect_ratio

  integer :: i, j
  integer :: ierr


  ! Default profile type (same as uniform)
  profile_type = "default"

  ! default aspect ratio sigma_x = sigma_y
  aspect_ratio = 1.0_p_k_part
  ! Minimum density for injection
  den_min = 0

  !default values for piecewise-linear profile
  num_x = -1
  density = 1.0

  x = - huge( 1.0_p_k_part )
  fx = 0.0_p_k_part


  ! default values for gaussian psources
  gauss_center = 0.0_p_k_part
  gauss_sigma = huge( 1.0_p_k_part )
  gauss_n_sigma = 0
  gauss_range(p_lower,:) = -gauss_sigma
  gauss_range(p_upper,:) = gauss_sigma

  ! default values for channel parameters
  channel_dir = 1 ! default along x1
  channel_r0 = -huge( 1.0_p_k_part ) ! no default r0
  channel_depth = -huge( 1.0_p_k_part ) ! no default depth
  channel_size = -huge( 1.0_p_k_part ) ! no default size
  channel_center = -huge( 1.0_p_k_part ) ! no default size
  channel_wall = 0.0 ! default wall size
  channel_pos = -huge( 1.0_p_k_part ) ! no default pos
  channel_bottom = 1.0 ! default bottom density

  ! default values for sphere parameters
  sphere_center = 0.0
  sphere_radius = 0.0

  ! default values for math function
  math_func_expr = "NO_FUNCTION_SUPPLIED!"

  ! Variables for fixed charge injection - not used in this psource
  sample_rate = 32

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

  ! set the default psource type
  if ( profile_type(1) == "default" ) then
    if (num_x > 0) then
       profile_type(1) = "piecewise-linear"
    else
      profile_type(1) = "uniform"
    endif
  endif

  ! set driver aspect ratio
  this%aspect_ratio = aspect_ratio

  ! minimum density for injection
  this%den_min = den_min
  this%density = density

  ! This variable is only used by t_psource_constq but it was simpler to include it here
  ! When using constant charge injection this parameter will be checked in a subclass
  do i = 1, p_x_dim
    this%sample_rate(i) = sample_rate(i)
  enddo


  ! parse non mult options first

  this % if_mult = .false.

  select case ( trim( profile_type(1) ) )
  case ( "piecewise-linear" )
    this % if_mult = .true.
  case ( "gaussian" )
    this % if_mult = .true.
  case ( "channel" )
    this%type = p_channel
    if (p_x_dim == 1) then
      if ( mpi_node() == 0 ) then
        write(0,*) "(*error*) Channel psource is not available in 1D"
        write(0,*) "(*error*) aborting..."
      endif
      stop
    endif

    ! set channel variables
    this%channel_dir = channel_dir
    if ((this%channel_dir < 1) .or. (this%channel_dir > p_x_dim)) then
      if ( mpi_node() == 0 ) then
        write(0,*) "(*error*) Invalid channel direction"
        write(0,*) "(*error*) aborting..."
      endif
      stop
    endif

    this%channel_r0 = channel_r0
    if ( this%channel_r0 <= 0.0 ) then
      if ( mpi_node() == 0 ) then
        write(0,*) "(*error*) Invalid channel_r0 value, must be > 0.0"
        write(0,*) "(*error*) aborting..."
      endif
      stop
    endif

    this%channel_depth = channel_depth
    if ( this%channel_depth < 0.0 ) then
      if ( mpi_node() == 0 ) then
        write(0,*) "(*error*) Invalid channel_depth value, must be >= 0.0"
        write(0,*) "(*error*) aborting..."
      endif
      stop
    endif

    this%channel_size = channel_size
    if ( this%channel_size <= 0.0 ) then
      if ( mpi_node() == 0 ) then
        write(0,*) "(*error*) Invalid channel diameter (channel_size) value &
        &must be > 0.0"
        write(0,*) "(*error*) aborting..."
      endif
      stop
    endif

    do j = 1, p_x_dim-1
      this%channel_center(j) = channel_center(j)
      if (this%channel_center(j) == -huge( 1.0_p_k_part )) then
        if ( mpi_node() == 0 ) then
          write(0,*) "(*error*) channel_center(",j,") was not specified"
          write(0,*) "(*error*) aborting..."
        endif
        stop
      endif
    enddo

    this%channel_wall = channel_wall
    do j = 1, 2
      this%channel_pos(j) = channel_pos(j)
      if (this%channel_pos(j) == -huge( 1.0_p_k_part )) then
        if ( mpi_node() == 0 ) then
          write(0,*) "(*error*) channel_pos(",j,") was not specified"
          write(0,*) "(*error*) aborting..."
        endif
        stop
      endif
    enddo

    this%channel_bottom = channel_bottom

  case ( "sphere" )
    this%type = p_sphere

    ! set sphere variables
    this%sphere_center = sphere_center
    this%sphere_radius = sphere_radius

  case ( "math func" )
    this%type = p_func

    ! set math func variables
    this%math_func_expr = trim(math_func_expr)

  case ( "uniform" )

    this%type = p_uniform

! if not all directions set to uniform use a separable variable function
    do j = 2, p_x_dim
      if ( (trim( profile_type(j) ) /= "uniform") ) this%if_mult = .true.
    enddo

  case default
    if ( mpi_node() == 0 ) then
      write(0,*) '(*error*) Invalid profile_type(1), "', trim(profile_type(1)), '"'
      write(0,*) '(*error*) aborting...'
    endif
    stop
  end select

  ! allocate x, fx if using piecewise-linear
  this%num_x = num_x
  if ((this%if_mult) .and. ( num_x > 0 )) then
    call alloc(this%x, (/ num_x, p_x_dim /),"cyl_modes/psource/os-psource-std-cyl-modes.f03",357)
    call alloc(this%fx, (/ num_x, p_x_dim /),"cyl_modes/psource/os-psource-std-cyl-modes.f03",358)
  endif

  if ( this%if_mult ) then
    ! parse mult options
    do j=1, p_x_dim
      if ( profile_type(j) == "default" ) then
        if (num_x > 0) then
          profile_type(j) = "piecewise-linear"
        else
          profile_type(j) = "uniform"
        endif
      endif


      select case ( trim( profile_type(j) ) )
      case ( "uniform" )
        this%type(j) = p_uniform

      case ( "piecewise-linear" )
        if (this%num_x <= 0) then
          if ( mpi_node() == 0 ) then
            write(0,*) "(*error*) Invalid num_x. num_x must be > 0 when using &
                        &'piecewise-linear' density psources"
            write(0,*) "(*error*) aborting..."
          endif
          stop
        endif
        this%type(j) = p_pw_linear

        do i = 1, num_x-1
          if ( x(i,j) >= x(i+1,j) ) then
            if ( mpi_node() == 0 ) then
              write(0,*) "(*error*) Invalid x value. All values must be larger than predecessor."
              write(0,*) "(*error*) aborting..."
            endif
            stop
          endif
        enddo

        ! store the piecewise linear parameters for this direction
        do i=1, num_x
          this%x(i,j) = x(i,j)
          this%fx(i,j) = fx(i,j)
        enddo

      case ( "gaussian" )
        this%type(j) = p_gaussian

        this%gauss_center(j) = gauss_center(j)
        this%gauss_sigma(j) = abs( gauss_sigma(j) )
        if ( gauss_n_sigma(j) > 0 ) then
          if ( this%aspect_ratio > 1.0_p_k_part .and. j == p_r_dim ) then
            this%gauss_range(p_lower,j) = gauss_center(j)-&
                                          abs(gauss_n_sigma(j)*gauss_sigma(j)*this%aspect_ratio)
            this%gauss_range(p_upper,j) = gauss_center(j)+&
                                          abs(gauss_n_sigma(j)*gauss_sigma(j)*this%aspect_ratio)
          else
            this%gauss_range(p_lower,j) = gauss_center(j)-abs(gauss_n_sigma(j)*gauss_sigma(j))
            this%gauss_range(p_upper,j) = gauss_center(j)+abs(gauss_n_sigma(j)*gauss_sigma(j))
          endif
        else
          this%gauss_range(:,j) = gauss_range(:,j)
        endif

      case default
        if ( mpi_node() == 0 ) then
          write(0,*) "(*error*) Invalid profile_type(",j,") : '", trim(profile_type(j)), "'"
          write(0,*) "(*error*) You can only use 'piecewise-linear', 'gaussian', and 'uniform' &
                      &if you are using a separable function psource."
          write(0,*) "(*error*) aborting..."
        endif
        stop
      end select
    enddo

  else
    ! check that no other options were given
    do j = 2, p_x_dim
      if (profile_type(j) /= "default") then
        if ( mpi_node() == 0 ) then
          write(0,*) "(*error*) Invalid profile_type(",j,") : '", trim(profile_type(j)), "'"
          write(0,*) "(*error*) When using : '",trim(profile_type(1)),"' you cannot &
                      &specify any other options for the profile_type"
          write(0,*) "(*error*) aborting..."
        endif
        stop
      endif
    enddo

    ! Check math function
    if ( this%type(1) == p_func ) then
      select case (p_x_dim)
      case (1)
        call setup(this%math_func, trim(this%math_func_expr), (/'x1'/), ierr)
      case (2)
        call setup(this%math_func, trim(this%math_func_expr), (/'x1','x2'/), ierr)
      case (3)
        call setup(this%math_func, trim(this%math_func_expr), (/'x1','x2','x3'/), ierr)
      end select

      if (ierr /= 0) then
        if ( mpi_node() == 0 ) then
          write(0,*) "(*error*) Invalid function supplied : '", &
                      trim(this%math_func_expr), "'"
        endif
        stop
      endif

    endif
  endif


end subroutine read_input_std_cyl_modes
!---------------------------------------------------------------------------------------------------




!---------------------------------------------------------------------------------------------------
! Gets density values for a set of npxx points (Quasi-3D)
!---------------------------------------------------------------------------------------------------
subroutine get_den_value_std_cyl_modes( this, pxx, npxx, pptheta, den_value )

  implicit none

  ! this needs to be inout because of the eval function
  class( t_psource_std_cyl_modes ), intent(inout) :: this
  real(p_k_part), dimension(:,:), intent(in) :: pxx
  integer, intent(in) :: npxx, pptheta
  real(p_k_part), dimension(:,:), intent(out) :: den_value

  ! local variables

  integer :: i, j, k, l
  real(p_k_part) :: r, theta
  integer :: i1, i2
  integer :: x_dim
  real(p_k_part), parameter :: pi = 3.14159265358979323846264_p_k_part
  ! executable statements
  if (this%type(1) == p_none) then
    den_value = 0.0_p_k_part
    return
  endif
  ! iterate over theta
  do l =1, pptheta
    if (this%if_mult) then ! psource is a product of p_x_dim functions

      den_value(:,l) = 1.0d0


      do i = 1, p_x_dim

        select case ( this%type(i) )

        case(p_pw_linear) ! piecewise linear -----------------------

          do k = 1, npxx
            if (pxx(i,k) <= this%x(1,i) ) then
              den_value(k,l) = den_value(k,l) * this%fx(1,i)
            elseif (pxx(i,k) >= this%x(this%num_x,i) ) then
              den_value(k,l) = den_value(k,l) * this%fx(this%num_x,i)
            else
              j = 2
              do while (pxx(i,k) > this%x(j,i))
                j = j+1
              enddo

              den_value(k,l) = den_value(k,l) * ( this%fx(j-1,i) &
                                + ( this%fx(j,i) - this%fx(j-1,i) ) &
                                / ( this%x( j,i) - this%x( j-1,i) ) &
                                * ( pxx(i,k) - this%x( j-1,i) ) )

            endif
          enddo


        case(p_gaussian) ! gaussian --------------------------------
          do k = 1, npxx

            if ((pxx(i,k) < this%gauss_range(p_lower,i)) .or. &
                (pxx(i,k) > this%gauss_range(p_upper,i)) ) then
              den_value(k,l) = 0.0
            else
              if(i < p_x_dim) then
                den_value(k,l) = den_value(k,l) * &
                exp(-(pxx(i,k)-this%gauss_center(i))**2/(2*this%gauss_sigma(i)**2 ))
              else
                if ( this%aspect_ratio == 1.0_p_k_part ) then
                  den_value(k,l) = den_value(k,l) * &
                  exp(-(pxx(i,k)-this%gauss_center(i))**2/(2*this%gauss_sigma(i)**2 ))
                else
                  ! need to do both x and y here
                  theta = (l - 1) * ( 2.0 * pi / pptheta ) + pxx(p_t_dim,k)
                  den_value(k,l) = den_value(k,l) * &
                  exp(-(pxx(i,k)-this%gauss_center(i))**2/(2*this%gauss_sigma(i)**2 )) * &
                  exp(-(pxx(i,k)-this%gauss_center(i))**2/(2*this%gauss_sigma(i)**2 ) * &
                    sin(theta)**2 * (1.0_p_k_part/this%aspect_ratio**2 -1.0_p_k_part) )
                endif
              endif
            endif

          enddo

        case( p_uniform ) ! uniform -----

          ! nothing required

        end select !case ( this%type(i) )

      enddo !i = 1, p_x_dim



    else ! (this%if_mult)

      ! psource is an arbitrary function of pxx

      select case ( this%type(1) )
      case(p_uniform) ! uniform ----------------------------------

        den_value(:,l) = 1.0d0

      case(p_func) ! math func -----------------------------------

        do k = 1, npxx
          den_value(k,l) = real( eval( this%math_func, real( pxx(:,k), p_k_fparse ) ), p_k_part )
        enddo

      case (p_sphere) ! sphere

        do k = 1, npxx
          r = ( pxx(1,k) - this%sphere_center(1) )**2
          do i = 2, p_x_dim
            r = r + ( pxx(i,k) - this%sphere_center(i) )**2
          enddo
          r = sqrt(r)

          if (r <= this%sphere_radius) then
            den_value(k,l) = 1.0_p_k_part
          else
            den_value(k,l) = 0.0_p_k_part
          endif
        enddo

      case (p_channel) ! parabolic channel

        if ( p_x_dim == 3 ) then
          i1 = 3 - this%channel_dir
          i2 = 3
          if (this%channel_dir == 3) then
            i1 = 1
            i2 = 2
          endif
        endif

        do k = 1, npxx

          if ((pxx(this%channel_dir,k) >= this%channel_pos(1)) .and. &
            (pxx(this%channel_dir,k) <= this%channel_pos(2))) then

            ! get distance to the center of the channel

            select case (p_x_dim)
            case (1)
              r = 0.0
            case (2)
              r = abs( pxx(3 - this%channel_dir,k) - this%channel_center(1))
            case (3)
              x_dim = p_x_dim-1
              r = sqrt((pxx(i1,k) - this%channel_center(1))**2 + &
                        (pxx(i2,k) - this%channel_center(x_dim))**2)

            end select

            ! get parabolic psource

            if (r <= this%channel_size/2) then ! inside the channel
              den_value(k,l) = this%channel_bottom + &
                                this%channel_depth * &
                                (r/this%channel_r0)**2
            else ! outside the channel

              if (this%channel_wall <= 0.0) then ! finite channel
                den_value(k,l) = this%channel_bottom + &
                                  this%channel_depth * &
                                  ((this%channel_size/2)/this%channel_r0)**2

              else ! leaky channel

                if (r < (this%channel_size/2 + this%channel_wall )) then
                  den_value(k,l) = (-r + this%channel_size/2 + this%channel_wall) / &
                                    this%channel_wall * (this%channel_bottom + &
                                    this%channel_depth * &
                                    ((this%channel_size/2)/this%channel_r0)**2)
                else
                  den_value(k,l) = 0.0d0
                endif

              endif

            endif

          else
            den_value(k,l) = 0.0d0
          endif

        enddo


      case default ! invalid psource or psource not specified
        den_value(1: npxx,l) = 1.0d0

      end select


    endif
  enddo
  den_value(1:npxx,1:pptheta) = den_value(1:npxx,1:pptheta) * this % density

end subroutine get_den_value_std_cyl_modes
!---------------------------------------------------------------------------------------------------


end module m_psource_std_cyl_modes


!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
function inject_std_cyl_modes( this, species, ig_xbnd_inj, jay, no_co, bnd_cross, node_cross, &
                                        send_msg, recv_msg ) result(num_inj)

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
# 691 "cyl_modes/psource/os-psource-std-cyl-modes.f03" 2

  use m_parameters
  use m_current_define, only : t_current
  use m_node_conf, only : t_node_conf
  use m_psource_std
  use m_psource_std_cyl_modes
  use m_species_cyl_modes_define
  use m_species_define, only: t_species, t_psource, t_part_idx, t_spec_msg, &
                              p_spec_buf_block
  use m_species_udist



  use m_random, only: rng

  implicit none

  ! dummy variables

  class( t_psource_std_cyl_modes), intent(inout) :: this
  class( t_species ), intent(inout), target :: species
  integer, dimension(:, :), intent(in) :: ig_xbnd_inj
  class( t_current ), intent(inout) :: jay
  class( t_node_conf ), intent(in) :: no_co
  type( t_part_idx ), dimension(2), intent(inout) :: bnd_cross
  type( t_part_idx ), intent(inout) :: node_cross
  type( t_spec_msg ), dimension(2), intent(inout) :: send_msg, recv_msg

  ! local variables
  integer, parameter :: p_max_x_dim = 3
  integer, parameter :: p_gaussian = 2

  ! number of particles to inject
  integer :: num_inj

  ! cell size (p_x_dim)
  real(p_double), dimension(p_max_x_dim) :: dx
  real(p_double), dimension(p_max_x_dim) :: g_xmin
  real(p_k_part) :: theta

  ! random theta offset
  real(p_k_part), dimension(:,:), pointer :: rand_theta
  real(p_k_part) :: spoke_shift
  integer, dimension(2) :: r_lb, r_ub ! boundaries for above array
  logical :: store_theta ! whether or not we need to store rand_theta values

  ! half distance between particles
  real(p_k_part), dimension(p_max_x_dim) :: dxp_2

  integer :: i, i1, i2, ipart, j

  real(p_k_part) :: x1

  ! volume that each particle occupies
  real(p_k_part) :: pvol

  ! number of particles per cell
  integer :: ppcell, ppc2
  real(p_k_part) :: ppc2_r

  ! number of particles per ring (high order cyl modes only )
  integer :: pptheta

  ! particle positions (global / inside cell )
  real(p_k_part), dimension(:,:), pointer :: ppos, ppos_cell

  ! particle charges
  real(p_k_part), dimension(:,:), pointer :: pcharge

  ! charge normalization coefficients for each particle (ppcell)
  real(p_k_part), dimension(:), pointer :: alpha

  class( t_species_cyl_modes ), pointer :: spec_cyl
  real(p_k_part), parameter :: pi = 3.14159265358979323846264_p_k_part

  ! Set total number of injected particles to 0
  num_inj = 0

  ! Get a pointer of the correct kind to the species object
  select type(species)
  class is (t_species_cyl_modes)
    spec_cyl => species
  class default
    write(err_buf__,*) 'inject_cyl_modes was called with a species of the wrong class';call err__("cyl_modes/psource/os-psource-std-cyl-modes.f03",774)
    call abort_program(p_err_invalid)
  end select

  ! check if num_par_x is > 0 for all directions
  ! if not return silently
  do i=1, p_x_dim
    if (spec_cyl%num_par_x(i) <= 0) return
  enddo

  ! Get minimum density to inject

  do i = 1, p_x_dim
    ! get cell size
    dx(i) = spec_cyl%dx(i)

    ! get global mininum (this is shifted by +0.5 cells from global simulation values)
    g_xmin( i ) = spec_cyl%g_box( p_lower , i )

    ! get half distance between particles
    dxp_2(i) = 0.5_p_k_part/spec_cyl%num_par_x(i)
  enddo

  ! find total number particles per cell
  ppcell = spec_cyl%num_par_x(1)
  do i = 2, p_x_dim
    ppcell = ppcell * spec_cyl%num_par_x(i)
  enddo
  ppc2 = spec_cyl%num_par_x(2)
  ppc2_r = real( ppc2, p_k_part )
  pptheta = spec_cyl%num_par_theta

  ! initialize temp buffers
  call alloc(ppos, (/p_x_dim+1, ppcell/),"cyl_modes/psource/os-psource-std-cyl-modes.f03",807) ! extra position for theta of first particle
  call alloc(ppos_cell, (/p_x_dim, ppcell/),"cyl_modes/psource/os-psource-std-cyl-modes.f03",808)
  call alloc(pcharge, (/ppcell,pptheta/),"cyl_modes/psource/os-psource-std-cyl-modes.f03",809)
  call alloc(alpha, (/ppcell/),"cyl_modes/psource/os-psource-std-cyl-modes.f03",810)

  ! If using an asymmetric gaussian profile in r, then store random theta values (costly).
  ! This will actually be necessary for any profile which uses theta values in
  ! get_den_value_cyl. Currently this is only the below configuration.
  store_theta = this%if_mult .and. this%type(p_r_dim) == p_gaussian &
                .and. this%aspect_ratio /= 1.0_p_k_part
  if (store_theta .and. spec_cyl%rand_theta) then
    r_lb = (/ ig_xbnd_inj(p_lower,1), ig_xbnd_inj(p_lower,2) /)
    r_ub = (/ ig_xbnd_inj(p_upper,1), ig_xbnd_inj(p_upper,2) /)
    call alloc(rand_theta, r_lb, r_ub,"cyl_modes/psource/os-psource-std-cyl-modes.f03",820)
  endif

  ! If we are not shifting spokes in a cell, then multiply the offset by 0.
  if (spec_cyl%shift_cell_spokes) then
    spoke_shift = 2.0_p_k_part * pi / pptheta / ppcell
  else
    spoke_shift = 0.0_p_k_part
  endif

  ! find normalization factor
  pvol = sign( 1.0_p_k_part/ppcell, spec_cyl%rqm )

  ! If using high order cylindrical modes account for additional particles along theta
  pvol = pvol / real(pptheta)


  ! Get position of particles inside the cell
  ! The position will always be in the range [-0.5, +0.5] regardless of interpolation type

  i = 0
  ! Only handle p_x_dim = 2
  do i1 = 1, spec_cyl%num_par_x(1)
    x1 = (2*i1 - 1 - spec_cyl%num_par_x(1)) * dxp_2(1)
    do i2 = 1, spec_cyl%num_par_x(2)
      i = i + 1
      ppos_cell(1,i) = x1
      ppos_cell(2,i) = (2*i2 - 1 - spec_cyl%num_par_x(2)) * dxp_2(2)
    enddo
  enddo

  ! we could write a recursive version of the following
  ! that would work for any number of dimensions...

  ipart = spec_cyl%num_par + 1
  num_inj = 0

  ! Inject particles
  ! Only handle p_x_dim = 2

  ! loop through all cells and count number of particles to inject
  ! don't inject on axis

  do i1 = ig_xbnd_inj(p_lower,1), ig_xbnd_inj(p_upper,1)
    ppos(1,:) = real( g_xmin(1) + (ppos_cell(1,:) + (i1-1))*dx(1), p_k_part )

    do i2 = ig_xbnd_inj(p_lower,2), ig_xbnd_inj(p_upper,2)
      ppos(2,:) = real( g_xmin(2) + (ppos_cell(2,:) + (i2-1))*dx(2), p_k_part )

      ! Calculate theta if needed for get_den_value_cyl
      ! Here we calculate one random number per cell. Each particle ring is then offset
      ! from that random number by an increment of 2*pi/pptheta/ppcell. Each particle in
      ! a given ring is spaced out by increments of 2*pi/pptheta. This ensures we sample
      ! the greatest space in theta within each cell.
      if (store_theta) then
        if (spec_cyl%rand_theta) then
          call rng % harvest_real2( rand_theta(i1,i2), &
            (/ i1 - spec_cyl%my_nx_p( p_lower, 1 ) + 1, &
               i2 - spec_cyl%my_nx_p( p_lower, 2 ) + 1 /) )
          rand_theta(i1,i2) = rand_theta(i1,i2) * 2.0_p_k_part * pi
          do i = 1, ppcell
            ppos(p_t_dim,i) = rand_theta(i1,i2) + (i-1) * spoke_shift
          enddo
        else
          do i = 1, ppcell
            ppos(p_t_dim,i) = (i-1) * spoke_shift
          enddo
        endif
      endif

      call this % get_den_value_cyl( ppos, ppcell, pptheta, pcharge )
      do i = 1, ppcell
        do j = 1, pptheta
          if ((pcharge(i,j) > this%den_min) .and. &
            (ppos(p_r_dim,i) > 0.0_p_k_part)) then

            num_inj = num_inj + 1
          endif
        enddo
      enddo

    enddo
  enddo

  ! When using high order modes account for multiple particles per ring
  num_inj = num_inj

  if ( num_inj > 0 ) then

    ! check if the buffer size is sufficient and grow it if necessary
    if ( num_inj > spec_cyl%num_par_max - spec_cyl%num_par ) then
      call spec_cyl%grow_buffer(spec_cyl%num_par_max + num_inj + p_spec_buf_block )
    endif

    ! loop through all the injection cells and
    ! inject particles, normalizing charge
    do i2 = ig_xbnd_inj(p_lower,2), ig_xbnd_inj(p_upper,2)
      ppos(2,:) = real( g_xmin(2) + (ppos_cell(2,:) + (i2-1))*dx(2), p_k_part )

      ! set the correct alpha coefficients for the radial position
      do i = 1, ppcell
        alpha(i) = spec_cyl % norm_charge_axis( i2, ppos(p_r_dim,i) )
      enddo

      do i1 = ig_xbnd_inj(p_lower,1), ig_xbnd_inj(p_upper,1)
        ppos(1,:) = real( g_xmin(1) + (ppos_cell(1,:) + (i1-1))*dx(1), p_k_part )


        ! See above comment for how we calculate theta positions of particles
        if (store_theta) then
          if (spec_cyl%rand_theta) then
            do i = 1, ppcell
              ppos(p_t_dim,i) = rand_theta(i1,i2) + (i-1) * spoke_shift
            enddo
          else
            do i = 1, ppcell
              ppos(p_t_dim,i) = (i-1) * spoke_shift
            enddo
          endif
        else
          if (spec_cyl%rand_theta) then
            call rng % harvest_real2( ppos(p_t_dim,1), &
              (/ i1 - spec_cyl%my_nx_p( p_lower, 1 ) + 1, &
                 i2 - spec_cyl%my_nx_p( p_lower, 2 ) + 1 /) )
            ppos(p_t_dim,1) = ppos(p_t_dim,1) * 2.0_p_k_part * pi
            do i = 2, ppcell
              ppos(p_t_dim,i) = ppos(p_t_dim,1) + (i-1) * spoke_shift
            enddo
          else
            do i = 1, ppcell
              ppos(p_t_dim,i) = (i-1) * spoke_shift
            enddo
          endif
        endif

        call this % get_den_value_cyl( ppos, ppcell, pptheta, pcharge )

        do i=1, ppcell

          do j = 1, pptheta
            if ((pcharge(i,j) > this%den_min) .and. &
              (ppos(p_r_dim,i) > 0.0_p_k_part)) then
              ! add particle
              spec_cyl%x(1, ipart) = ppos_cell(1,i)
              spec_cyl%x(2, ipart) = ppos_cell(2,i)
              spec_cyl%ix(1, ipart) = i1 - spec_cyl%my_nx_p( p_lower, 1 ) + 1
              spec_cyl%ix(2, ipart) = i2 - spec_cyl%my_nx_p( p_lower, 2 ) + 1
              spec_cyl%q(ipart) = pcharge(i,j) * pvol * ppos( p_r_dim, i ) * alpha(i)
              ! add cartesian data for higher-order modes
              theta = (j - 1) * ( 2.0_p_k_part * pi / pptheta ) + ppos(p_t_dim,i)
              spec_cyl%x(3, ipart) = ppos( p_r_dim, i )*cos( theta )
              spec_cyl%x(4, ipart) = ppos( p_r_dim, i )*sin( theta )

              ipart = ipart + 1
            endif
          enddo

        enddo
      enddo
    enddo

  endif

  call freemem(ppos,"cyl_modes/psource/os-psource-std-cyl-modes.f03",983)
  call freemem(ppos_cell,"cyl_modes/psource/os-psource-std-cyl-modes.f03",984)
  call freemem(pcharge,"cyl_modes/psource/os-psource-std-cyl-modes.f03",985)
  call freemem(alpha,"cyl_modes/psource/os-psource-std-cyl-modes.f03",986)
  if (store_theta .and. spec_cyl%rand_theta) then
    call freemem(rand_theta,"cyl_modes/psource/os-psource-std-cyl-modes.f03",988)
  endif

  if ( num_inj > 0 ) then

    ! set momentum of injected particles
    call set_momentum( spec_cyl, spec_cyl%num_par+1, spec_cyl%num_par + num_inj )




  endif

end function inject_std_cyl_modes
!-----------------------------------------------------------------------------------------
