# 1 "spec/psource/os-psource-python.f03"
# 1 "<built-in>" 1
# 1 "<built-in>" 3
# 467 "<built-in>" 3
# 1 "<command line>" 1
# 1 "<built-in>" 2
# 1 "spec/psource/os-psource-python.f03" 2
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
# 2 "spec/psource/os-psource-python.f03" 2
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
# 3 "spec/psource/os-psource-python.f03" 2

module m_psource_python

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
# 7 "spec/psource/os-psource-python.f03" 2

use m_species_define, only : t_species, t_part_idx, t_spec_msg
use m_parameters
use m_psource_std, only : t_psource_std
use m_fparser, only : p_max_expr_len
use m_current_define, only : t_current
use m_node_conf, only : t_node_conf
use m_vdf_define, only : t_vdf
use m_input_file, only : t_input_file, get_namelist
use m_callpy

implicit none

private

type, extends( t_psource_std ) :: t_psource_python

  ! variables for python function fields
  character(len = p_max_expr_len) :: py_mod = "", py_func = ""

  ! grid to hold the python density
  type( t_vdf ) :: den_grid

  ! parameters for the grid
  real(p_k_part), dimension(p_max_dim) :: xmin, xmax
  real(p_double), dimension(p_max_dim) :: dx

contains

  procedure :: read_input => read_input_python
  procedure :: if_inject => if_inject_python
  procedure :: inject => inject_python
  procedure :: cleanup => cleanup_python
  procedure :: get_den_value => get_den_value_python

  procedure :: load_density => load_density_python

end type t_psource_python

! Add interface to superclass inject function to avoid type casting
interface
  function inject_std( this, species, ig_xbnd_inj, jay, no_co, bnd_cross, node_cross, &
                         send_msg, recv_msg ) result(num_inj)

  import t_psource_std, t_species, t_current, t_node_conf, t_part_idx, t_spec_msg

  class(t_psource_std), intent(inout) :: this
  class(t_species), intent(inout), target :: species
  integer, dimension(:,:), intent(in) :: ig_xbnd_inj
  class( t_current ), intent(inout) :: jay
  class( t_node_conf ), intent(in) :: no_co
  type( t_part_idx ), dimension(2), intent(inout) :: bnd_cross
  type( t_part_idx ), intent(inout) :: node_cross
  type( t_spec_msg ), dimension(2), intent(inout) :: send_msg, recv_msg

  integer :: num_inj

  end function
end interface

interface interpolate
  module procedure interpolate_vdf_simple
end interface

public :: t_psource_python, interpolate

contains

!-------------------------------------------------------------------------------
! Read information from input file
!-------------------------------------------------------------------------------
subroutine read_input_python( this, input_file, coordinates )

  implicit none

  class( t_psource_python ), intent(inout) :: this
  class( t_input_file ), intent(inout) :: input_file
  integer, intent(in) :: coordinates

  ! minimum density to inject particles
  real(p_k_part) :: den_min

  ! global density (multiplies the selected density profile)
  real(p_k_part) :: density

  ! python
  character(len = p_max_expr_len) :: py_mod, py_func

  ! namelists of input variables
  namelist /nl_profile/ den_min, density, py_mod, py_func

  integer :: ierr

  if ( .not. if_py_util() ) then
    if ( mpi_node() == 0 ) then
      write(0,*) "   Error: requested ""python"" for species init_type, but code was"
      write(0,*) "   not compiled with python support. Please see the"
      write(0,*) "   PY_FCOMPILEFLAGS and PY_FLINKFLAGS flags in your config file."
      write(0,*) "   aborting..."
    endif
    stop
  endif

  ! Minimum density for injection
  den_min = 0

  ! default density
  density = 1.0

  ! python defaults
  py_mod = "NO_PYTHON_MODULE_SUPPLIED!"
  py_func = "NO_PYTHON_FUNCTION_SUPPLIED!"

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
    if ( mpi_node() == 0 ) then
      write(0,*) ""
      write(0,*) "   Error, profile parameters missing for python initialization"
      write(0,*) "   aborting..."
    endif
    stop
  endif

  ! minimum density for injection
  this%den_min = den_min
  this%density = density

  this%py_mod = trim(py_mod)
  this%py_func = trim(py_func)

end subroutine read_input_python
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
function if_inject_python( this )

  implicit none

  logical :: if_inject_python
  class( t_psource_python ), intent( in ) :: this

  if_inject_python = .true.

end function if_inject_python
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
function inject_python( this, species, ig_xbnd_inj, jay, no_co, bnd_cross, node_cross, &
                        send_msg, recv_msg ) result(num_inj)
!-----------------------------------------------------------------------------------------
! Inject the particles from a python function
!-----------------------------------------------------------------------------------------

  implicit none

  class(t_psource_python), intent(inout) :: this
  class(t_species), intent(inout), target :: species
  integer, dimension(:,:), intent(in) :: ig_xbnd_inj
  class( t_current ), intent(inout) :: jay
  class( t_node_conf ), intent(in) :: no_co
  type( t_part_idx ), dimension(2), intent(inout) :: bnd_cross
  type( t_part_idx ), intent(inout) :: node_cross
  type( t_spec_msg ), dimension(2), intent(inout) :: send_msg, recv_msg

  integer :: num_inj

  if ( this%den_grid%f_dim() < 1 ) then

    call this % load_density()

  endif

  num_inj = inject_std( this, species, ig_xbnd_inj, jay, no_co, bnd_cross, node_cross, &
                        send_msg, recv_msg )

end function inject_python
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Cleanup t_psource_python object
!-----------------------------------------------------------------------------------------
subroutine cleanup_python( this )

  implicit none

  class( t_psource_python ), intent( inout ) :: this

  call this%den_grid%cleanup()

  call this%t_psource_std%cleanup()

end subroutine cleanup_python
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine load_density_python( this )

  implicit none

  class( t_psource_python ), intent( inout ) :: this

  integer, dimension(p_x_dim) :: nx
  integer, dimension(2, p_x_dim) :: gc_num
  integer :: i

  ! Call python interface to get data
  call call_function( trim(this%py_mod), trim(this%py_func) )

  call get_state( "nx", nx )
  call get_state( "xmin", this%xmin(1:p_x_dim) )
  call get_state( "xmax", this%xmax(1:p_x_dim) )

  do i = 1, p_x_dim
    if (nx(i) <= 1) then
      write(err_buf__,*) 'In load_density_python, nx must be > 1 in each dimension.';call err__("spec/psource/os-psource-python.f03",233)
      write(err_buf__,*) 'Incompatible size of density grid.';call err__("spec/psource/os-psource-python.f03",234)
      call abort_program(p_err_rstrd)
    endif
  enddo

  this%dx(1:p_x_dim) = real((this%xmax(1:p_x_dim) - this%xmin(1:p_x_dim)) / (nx - 1), &
                            p_double)
  gc_num = 0

  ! Initialize a vdf to fit the python density data
  call this%den_grid%new( p_x_dim, 1, nx, gc_num, this%dx(1:p_x_dim), .false. )

  ! Populate the initialized vdf with the python data (assume to match)
  select case( p_x_dim )
  case(1)
    call get_state( "data", this%den_grid%f1(1, :) )
  case(2)
    call get_state( "data", this%den_grid%f2(1, :, :) )
  case(3)
    call get_state( "data", this%den_grid%f3(1, :, :, :) )
  end select

end subroutine load_density_python
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Gets density values for a set of npxx points
!-----------------------------------------------------------------------------------------
subroutine get_den_value_python( this, pxx, npxx, den_value )

  implicit none

  ! this needs to be inout because of the eval function
  class( t_psource_python ), intent(inout) :: this
  real(p_k_part), dimension(:,:), intent(in) :: pxx
  integer, intent(in) :: npxx
  real(p_k_part), dimension(:), intent(out) :: den_value

  integer :: i, k
  real(p_k_part), dimension(p_max_dim) :: xmin, xmax
  real(p_double), dimension(p_max_dim) :: dx
  integer, dimension(p_max_dim) :: nx, ix, ix_a, ix_b
  real(p_k_fld), dimension(p_max_dim) :: x, x_a, x_b

  ! Get local copies
  xmin = this%xmin
  xmax = this%xmax
  dx = this%dx
  nx(1:p_x_dim) = this%den_grid%nx_(1:p_x_dim)

  ! Loop over particles
  do k = 1, npxx

    ! Get the ix and x values in each dimension
    do i = 1, p_x_dim

      ! The if statements are done separately to be faster, see ntrim in os-spec-push.f03

      ! Check if out of bounds lower
      if (pxx(i,k) <= xmin(i)) then
        ix_a(i) = 1
      else
        ix_a(i) = 0
      endif
      x_a(i) = 0.0_p_k_fld

      ! Check if out of bounds upper
      if (pxx(i,k) >= xmax(i)) then
        ix_b(i) = nx(i) - 1
        x_b(i) = 1.0_p_k_fld
      else
        ix_b(i) = 0
        x_b(i) = 0.0_p_k_fld
      endif

      ix(i) = ix_a(i) + ix_b(i)
      x(i) = x_a(i) + x_b(i)

      if (ix(i) == 0) then

        ! Find the coordinate of this particle within the grid
        x(i) = (pxx(i,k) - xmin(i)) / dx(i)
        ix(i) = int(x(i))
        x(i) = x(i) - ix(i)
        ix(i) = ix(i) + 1 ! For Fortran indexing

      endif

    enddo

    ! Interpolate to find the density value
    den_value(k) = interpolate_vdf_simple( this%den_grid, ix(1:p_x_dim), x(1:p_x_dim) )

  enddo

  den_value(1: npxx) = den_value(1: npxx) * this % density

end subroutine get_den_value_python
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Perform 1D, 2D, or 3D interpolation to find the density.
! The ix values must be within the vdf indices, and the x values range between 0 and 1.
! Thus, a point at [0.25, 0.75] on a grid with 2 cells going from 0 to 1 in each dimension
! will have values ix = [0, 0] and x = [0.25, 0.75].
!-----------------------------------------------------------------------------------------
function interpolate_vdf_simple( vdf, ix, x )

  implicit none

  type( t_vdf ), intent(in) :: vdf
  integer, dimension(:), intent(in) :: ix
  real(p_k_fld), dimension(:), intent(in) :: x

  real(p_k_part) :: interpolate_vdf_simple
  real(p_k_fld) :: den

  real(p_k_fld), dimension(p_max_dim) :: x1
  integer, dimension(p_max_dim) :: ix1

  x1(1:p_x_dim) = 1.0_p_k_fld - x(1:p_x_dim)
  ix1(1:p_x_dim) = ix(1:p_x_dim) + 1

  select case( p_x_dim )
  case(1)

    den = vdf%f1(1,ix(1)) * x1(1) + vdf%f1(1,ix1(1)) * x(1)

  case(2)

    den = vdf%f2(1, ix(1), ix(2)) * x1(1) * x1(2) &
        + vdf%f2(1,ix1(1), ix(2)) * x(1) * x1(2) &
        + vdf%f2(1, ix(1),ix1(2)) * x1(1) * x(2) &
        + vdf%f2(1,ix1(1),ix1(2)) * x(1) * x(2)

  case(3)

    den = vdf%f3(1, ix(1), ix(2), ix(3)) * x1(1) * x1(2) * x1(3) &
        + vdf%f3(1,ix1(1), ix(2), ix(3)) * x(1) * x1(2) * x1(3) &
        + vdf%f3(1, ix(1),ix1(2), ix(3)) * x1(1) * x(2) * x1(3) &
        + vdf%f3(1,ix1(1),ix1(2), ix(3)) * x(1) * x(2) * x1(3) &
        + vdf%f3(1, ix(1), ix(2),ix1(3)) * x1(1) * x1(2) * x(3) &
        + vdf%f3(1,ix1(1), ix(2),ix1(3)) * x(1) * x1(2) * x(3) &
        + vdf%f3(1, ix(1),ix1(2),ix1(3)) * x1(1) * x(2) * x(3) &
        + vdf%f3(1,ix1(1),ix1(2),ix1(3)) * x(1) * x(2) * x(3)

  end select

  interpolate_vdf_simple = real( den, p_k_part )

end function interpolate_vdf_simple
!-----------------------------------------------------------------------------------------

end module m_psource_python
