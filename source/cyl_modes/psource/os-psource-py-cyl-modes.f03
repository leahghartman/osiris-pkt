#include "os-config.h"
#include "os-preprocess.fpp"

module m_psource_py_cyl_modes

#include "memory/memory.h"

use m_species_define, only : t_species, t_part_idx, t_spec_msg
use m_parameters
use m_psource_std_cyl_modes, only : t_psource_std_cyl_modes
use m_psource_python, only : interpolate
use m_fparser, only : p_max_expr_len
use m_current_define, only : t_current
use m_node_conf, only : t_node_conf
use m_vdf_define, only : t_vdf
use m_input_file, only : t_input_file, get_namelist
use m_callpy

implicit none

private

type, extends( t_psource_std_cyl_modes ) :: t_psource_py_cyl_modes

  ! variables for python function fields
  character(len = p_max_expr_len) :: py_mod = "", py_func = ""

  ! grid to hold the python density
  type( t_vdf ) :: den_grid

  ! parameters for the grid
  real(p_k_part), dimension(p_x_dim) :: xmin, xmax
  real(p_double), dimension(p_x_dim) :: dx

contains

  procedure :: read_input => read_input_py_cyl_modes
  procedure :: if_inject => if_inject_py_cyl_modes
  procedure :: inject => inject_py_cyl_modes
  procedure :: cleanup => cleanup_py_cyl_modes
  procedure :: get_den_value_cyl => get_den_value_py_cyl_modes

  procedure :: load_density => load_density_py_cyl_modes

end type t_psource_py_cyl_modes

! Add interface to superclass inject function to avoid type casting
interface
  function inject_std_cyl_modes( this, species, ig_xbnd_inj, jay, no_co, bnd_cross, &
                                 node_cross, send_msg, recv_msg ) result(num_inj)

  import t_psource_std_cyl_modes, t_species, t_current, t_node_conf, t_part_idx, &
         t_spec_msg

  class( t_psource_std_cyl_modes), intent(inout) :: this
  class( t_species ), intent(inout), target :: species
  integer, dimension(:, :), intent(in) :: ig_xbnd_inj
  class( t_current ), intent(inout)   :: jay
  class( t_node_conf ), intent(in)     :: no_co
  type( t_part_idx ), dimension(2), intent(inout) :: bnd_cross
  type( t_part_idx ), intent(inout) :: node_cross
  type( t_spec_msg ), dimension(2), intent(inout) :: send_msg, recv_msg

  integer :: num_inj

  end function
end interface

public :: t_psource_py_cyl_modes

contains

!-------------------------------------------------------------------------------
! Read information from input file
!-------------------------------------------------------------------------------
subroutine read_input_py_cyl_modes( this, input_file, coordinates )

  implicit none

  class( t_psource_py_cyl_modes ), intent(inout) :: this
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

end subroutine read_input_py_cyl_modes
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
function if_inject_py_cyl_modes( this )

  implicit none

  logical :: if_inject_py_cyl_modes
  class( t_psource_py_cyl_modes ), intent( in ) :: this

  if_inject_py_cyl_modes = .true.

end function if_inject_py_cyl_modes
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
function inject_py_cyl_modes( this, species, ig_xbnd_inj, jay, no_co, bnd_cross, &
                              node_cross, send_msg, recv_msg ) result(num_inj)
!-----------------------------------------------------------------------------------------
! Inject the particles from a python function
!-----------------------------------------------------------------------------------------

  implicit none

  class(t_psource_py_cyl_modes), intent(inout) :: this
  class(t_species), intent(inout), target :: species
  integer, dimension(:,:), intent(in) :: ig_xbnd_inj
  class( t_current ), intent(inout)   :: jay
  class( t_node_conf ), intent(in)     :: no_co
  type( t_part_idx ), dimension(2), intent(inout) :: bnd_cross
  type( t_part_idx ), intent(inout) :: node_cross
  type( t_spec_msg ), dimension(2), intent(inout) :: send_msg, recv_msg

  integer :: num_inj

  if ( this%den_grid%f_dim() < 1 ) then

    call this % load_density()

  endif

  num_inj = inject_std_cyl_modes( this, species, ig_xbnd_inj, jay, no_co, bnd_cross, &
                                  node_cross, send_msg, recv_msg )

end function inject_py_cyl_modes
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Cleanup t_psource_py_cyl_modes object
!-----------------------------------------------------------------------------------------
subroutine cleanup_py_cyl_modes( this )

  implicit none

  class( t_psource_py_cyl_modes ), intent( inout ) :: this

  call this%den_grid%cleanup()

  call this%t_psource_std%cleanup()

end subroutine cleanup_py_cyl_modes
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine load_density_py_cyl_modes( this )

  implicit none

  class( t_psource_py_cyl_modes ), intent( inout ) :: this

  integer, dimension(p_x_dim) :: nx
  integer, dimension(2, p_x_dim) :: gc_num
  integer :: i

  ! Call python interface to get data
  call call_function( trim(this%py_mod), trim(this%py_func) )

  call get_state( "nx", nx )
  call get_state( "xmin", this%xmin )
  call get_state( "xmax", this%xmax )

  do i = 1, p_x_dim
    if (nx(i) <= 1) then
      ERROR('In load_density_py_cyl_modes, nx must be > 1 in each dimension.')
      ERROR('Incompatible size of density grid.')
      call abort_program(p_err_rstrd)
    endif
  enddo

  this%dx = real((this%xmax - this%xmin) / (nx - 1), p_double)
  gc_num = 0

  ! Initialize a vdf to fit the python density data
  call this%den_grid%new( p_x_dim, 1, nx, gc_num, this%dx, .false. )

  ! Populate the initialized vdf with the python data (assume to match)
  call get_state( "data", this%den_grid%f2(1, :, :) )

end subroutine load_density_py_cyl_modes
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Gets density values for a set of npxx points
!-----------------------------------------------------------------------------------------
subroutine get_den_value_py_cyl_modes( this, pxx, npxx, pptheta, den_value )

  implicit none

  ! this needs to be inout because of the eval function
  class( t_psource_py_cyl_modes ), intent(inout) :: this
  real(p_k_part), dimension(:,:), intent(in) :: pxx
  integer, intent(in) :: npxx, pptheta
  real(p_k_part), dimension(:,:), intent(out) :: den_value

  integer :: i, k, l
  real(p_k_part), dimension(p_x_dim) :: xmin, xmax
  real(p_double), dimension(p_x_dim) :: dx
  integer, dimension(p_x_dim) :: nx, ix, ix_a, ix_b
  real(p_k_fld), dimension(p_x_dim) :: x, x_a, x_b

  ! Get local copies
  xmin = this%xmin
  xmax = this%xmax
  dx = this%dx
  nx = this%den_grid%nx_(1:p_x_dim)

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
      x_a(i)  = 0.0_p_k_fld

      ! Check if out of bounds upper
      if (pxx(i,k) >= xmax(i)) then
        ix_b(i) = nx(i) - 1
        x_b(i)  = 1.0_p_k_fld
      else
        ix_b(i) = 0
        x_b(i)  = 0.0_p_k_fld
      endif

      ix(i) = ix_a(i) + ix_b(i)
      x(i)  =  x_a(i) +  x_b(i)

      if (ix(i) == 0) then

        ! Find the coordinate of this particle within the grid
        x(i) = (pxx(i,k) - xmin(i)) / dx(i)
        ix(i) = int(x(i))
        x(i) = x(i) - ix(i)
        ix(i) = ix(i) + 1 ! For Fortran indexing

      endif

    enddo

    ! Interpolate to find the density value
    den_value(k, 1) = interpolate( this%den_grid, ix, x )

  enddo

  den_value(1:npxx, 1) = den_value(1:npxx, 1) * this % density

  ! Copy the density values to the other particles in theta
  do l = 2, pptheta
    do k = 1, npxx
      den_value(k, l) = den_value(k, 1)
    enddo
  enddo

end subroutine get_den_value_py_cyl_modes
!-----------------------------------------------------------------------------------------

end module m_psource_py_cyl_modes
