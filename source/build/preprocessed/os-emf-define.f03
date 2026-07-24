# 1 "emf/os-emf-define.f03"
# 1 "<built-in>" 1
# 1 "<built-in>" 3
# 467 "<built-in>" 3
# 1 "<command line>" 1
# 1 "<built-in>" 2
# 1 "emf/os-emf-define.f03" 2
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
# 2 "emf/os-emf-define.f03" 2
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
# 3 "emf/os-emf-define.f03" 2

module m_emf_define

use m_system
use m_parameters

use m_vdf
use m_vdf_define
use m_vdf_smooth
use m_vdf_report
use m_vdf_comm

use m_wall_define

use m_diagnostic_utilities

use m_fparser

use m_input_file

use m_space, only : t_space
use m_grid_define, only : t_grid
use m_node_conf, only : t_node_conf
use m_time_step, only : t_time_step
use m_restart, only : t_restart_handle
use m_current_define, only : t_current

use m_emf_diag
use m_logprof, only : begin_event, end_event

implicit none

private

! string to id restart data
character(len=*), parameter, public :: p_emf_rst_id = "emf rst data - 0x0006"

!-------------------------------------------------------------------------------
! Module constants and variables
!-------------------------------------------------------------------------------

! type of smooth
integer, parameter :: p_emfsmooth_none = 0, & ! No smooth
                      p_emfsmooth_stand = 1, & ! Smooth E&B for particle push only
                      p_emfsmooth_nci = 2, & ! NCI field smooth (for particle push only)
                      p_emfsmooth_local = 3 ! Direct E&B smooth

! Type of external field to use
integer, parameter :: p_extfld_none = 0, &
                      p_extfld_static = 1, &
                      p_extfld_dynamic = 2

! constants for external/initial fields
integer, parameter :: p_emf_none = 0, &
                      p_emf_uniform = 1, &
                      p_emf_math = 2, &
                      p_emf_dipole = 3, &
                      p_emf_python = 4

integer, parameter, dimension(3) :: p_emf_hc_e = (/ p_hc_x, p_hc_y, p_hc_z /)

integer, parameter, dimension(3) :: p_emf_hc_b = (/ p_hc_y + p_hc_z, &
                                                    p_hc_x + p_hc_z, &
                                                    p_hc_x + p_hc_y /)

! timing events
integer :: fsolverev = 0, fsmoothev, getpsi_ev

!-------------------------------------------------------------------------------
! Psi calculation class definition
!-------------------------------------------------------------------------------

type :: t_emf_psi

  ! vdf holding psi values
  type( t_vdf ) :: data

  ! Iteration when psi was last calculated
  integer :: n = -1

  ! Communicator used for parallel psi calculations
  integer :: comm = mpi_comm_null

end type t_emf_psi

!-------------------------------------------------------------------------------
! Boundary conditions
!-------------------------------------------------------------------------------

# 1 "emf/boundary/os-emf-bnd-define.f03" 1

! string to id restart data
character(len=*), parameter :: p_bcemf_rst_id = "bcemf rst data - 0x0002"

!-------------------------------------------------------------------------------
! t_lindman class definition
!-------------------------------------------------------------------------------

type:: t_lindman

    integer :: idir = -1 ! direction of lindman bc (1->p_x_dim)
    integer :: orientation ! upper or lower wall
    integer, dimension(p_f_dim-1) :: perp ! perpendicular field directions
    integer, dimension(2) :: range ! range of buffers

    type (t_wall) :: e_buffer,b_buffer ! e and b buffers for storing the
                                                    ! fields at previous timesteps

    real (p_k_fld), dimension(p_max_dim) :: dtdxi ! (dt/dx/2)
    real (p_k_fld) :: disp_corr ! dispersion correction
                                                    ! (1-dx/(dt*c))/(1+dx/(dt*c))

    integer , dimension(p_max_dim) :: imin, imax

    integer , dimension(2) :: bound

end type t_lindman

!-------------------------------------------------------------------------------
! t_vpml class definition
!-------------------------------------------------------------------------------

type:: t_vpml

    ! Simulation coordinates (cartesian/r-z cylindrical)
    integer :: coordinates

    ! Position of lower axial boundary for r-z cylindrical simulations
    integer :: gir_pos

    ! Boundaries are sorted as follows
    ! x1_lower, x1_upper, x2_lower, x2_upper, x3_lower, x3_upper

    ! Number of vpml boundaries
    integer :: n_vpml = 0

    ! Current vpml when looping over update boundary with tiles
    integer :: c

    ! Perp. direction of the boundary
    integer, dimension(:), pointer :: dir_vpml => null()

    ! Boundary location (p_lower, p_upper)
    integer, dimension(:), pointer :: loc_vpml => null()

    ! Number of cells in vpml boundaries
    integer :: bnd_size

    logical, dimension(2, p_x_dim) :: if_diffuse

    ! dt/dx
    real(p_k_fld), dimension(p_x_dim) :: dtdx

    ! Contact walls for each corner
    integer, dimension(:,:,:), pointer :: pos_corner => null()

    ! Coeficients of E and B in eq (27) of Vay's paper: JCP 165, 511 (2000)
    ! [
    ! direction => wall direction (p_x_dim elements)
    ! coeff. idx => index of coefficients in eq. (27) of Vay's paper (3 elements)
    ! array pos. => position in the vpml boundary (bnd_size elements)
    ! ]

    real(p_k_fld), dimension(:,:,:), pointer :: coef_e_low => null()

    real(p_k_fld), dimension(:,:,:), pointer :: coef_b_low => null()

    real(p_k_fld), dimension(:,:,:), pointer :: coef_e_up => null()

    real(p_k_fld), dimension(:,:,:), pointer :: coef_b_up => null()

    ! 2 walls per vpml boundary
    type (t_wall), pointer, dimension(:) :: wall_array_e => null()

    type (t_wall), pointer, dimension(:) :: wall_array_b => null()

end type t_vpml

!-------------------------------------------------------------------------------
! t_emf_bound class definition
!-------------------------------------------------------------------------------

type :: t_emf_bound

    ! type of boundary conditions for each boundary

    !indexes: (lower-upper bound, direction)
    integer, dimension(2,p_max_dim) :: type

    integer, dimension(2,p_max_dim) :: global_type

    ! lindman boundary data
    type (t_lindman), dimension(2):: lindman

    ! pml boundary data (all boundaries)
    type (t_vpml) :: vpml_all
    integer :: vpml_bnd_size
    logical, dimension(2,p_x_dim) :: vpml_diffuse

contains

    procedure :: read_input => read_input_emf_bound
    procedure :: init => init_emf_bound
    procedure :: cleanup => cleanup_emf_bound
    procedure :: write_checkpoint => write_checkpoint_emf_bound
    procedure :: restart_read => restart_read_emf_bound
    procedure :: update_boundary => update_boundary_emf_bound
    procedure :: update_boundary_b => update_boundary_bfld_emf_bound
    procedure :: update_boundary_e => update_boundary_efld_emf_bound
    procedure :: move_window => move_window_emf_bound
    procedure :: list_algorithm => list_algorithm_emf_bound
    procedure :: reshape => reshape_emf_bound

end type t_emf_bound

interface
    subroutine read_input_emf_bound( this, input_file, periodic, if_move )
        import t_emf_bound, t_input_file
        class (t_emf_bound), intent( inout ) :: this
        class( t_input_file ), intent(inout) :: input_file
        logical, dimension(:), intent(in) :: periodic, if_move
    end subroutine
end interface

interface
    subroutine init_emf_bound( this, b, e, dt, part_interpolation, &
                               space, no_co, grid, restart, restart_handle )
        import t_emf_bound, t_vdf, t_space, t_node_conf, t_grid, t_restart_handle, p_double
        class (t_emf_bound), intent( inout ) :: this
        type( t_vdf ), intent(in) :: b, e
        real( p_double ), intent(in) :: dt
        integer, intent(in) :: part_interpolation
        type( t_space ), intent(in) :: space
        class( t_node_conf ), intent(in) :: no_co
        class( t_grid ), intent(in) :: grid
        logical, intent(in) :: restart
        type( t_restart_handle ), intent(in) :: restart_handle
    end subroutine
end interface

interface
    subroutine cleanup_emf_bound( this )
        import t_emf_bound
        class (t_emf_bound), intent( inout ) :: this
    end subroutine
end interface

interface
    subroutine write_checkpoint_emf_bound( this, restart_handle )
        import t_emf_bound, t_restart_handle
        class (t_emf_bound), intent(in) :: this
        type( t_restart_handle ), intent(inout) :: restart_handle
    end subroutine
end interface

interface
    subroutine restart_read_emf_bound( this, restart_handle )
        import t_emf_bound, t_restart_handle
        class (t_emf_bound), intent(inout) :: this
        type( t_restart_handle ), intent(in) :: restart_handle
    end subroutine
end interface

interface
    subroutine update_boundary_emf_bound( this, b, e, g_space, no_co, send_msg, recv_msg )
        import t_emf_bound, t_vdf, t_space, t_node_conf, t_vdf_msg
        class( t_emf_bound ), intent(inout) :: this
        type( t_vdf ), intent(inout) :: b, e
        type( t_space ), intent(in) :: g_space
        class( t_node_conf ), intent(in) :: no_co
        type(t_vdf_msg), dimension(2), intent(inout) :: send_msg, recv_msg
    end subroutine
end interface

interface
    subroutine update_boundary_bfld_emf_bound( this, e, b, step )
        import t_emf_bound, t_vdf
        class( t_emf_bound ), intent(inout) :: this
        class( t_vdf ), intent(inout) :: e, b
        integer, intent(in) :: step
    end subroutine
end interface

interface
    subroutine update_boundary_efld_emf_bound( this, e, b )
        import t_emf_bound, t_vdf
        class( t_emf_bound ), intent(inout) :: this
        class( t_vdf ), intent(inout) :: e, b
    end subroutine
end interface

interface
    subroutine move_window_emf_bound( this, space )
        import t_emf_bound, t_space
        class( t_emf_bound ), intent(inout) :: this
        type( t_space ), intent(in) :: space
    end subroutine
end interface

interface
    subroutine list_algorithm_emf_bound( this )
        import t_emf_bound
        class( t_emf_bound ), intent(in) :: this
    end subroutine list_algorithm_emf_bound
end interface

interface
    subroutine reshape_emf_bound( this, old_lb, new_lb, no_co, send_msg, recv_msg )
        import t_emf_bound, t_grid, t_node_conf, t_vdf_msg
        class( t_emf_bound ), intent(inout) :: this
        class(t_grid), intent(in) :: old_lb, new_lb
        class( t_node_conf ), intent(in) :: no_co
        type(t_vdf_msg), dimension(2), intent(inout) :: send_msg, recv_msg
    end subroutine reshape_emf_bound
end interface

public :: p_bcemf_rst_id, t_lindman, t_vpml, t_emf_bound
# 93 "emf/os-emf-define.f03" 2

!-------------------------------------------------------------------------------
! t_emf_solver class definition
!-------------------------------------------------------------------------------

# 1 "emf/solver/os-emf-solver-define.f03" 1
! This file is meant to be included by os-emf-define.f03

integer, parameter, public :: p_maxlen_emf_solver_name = 32

type, abstract, public :: t_emf_solver

  integer :: gir_pos

contains
  procedure(init_emf_solver), deferred :: init
  procedure(advance_emf_solver), deferred :: advance
  procedure(min_gc_emf_solver), deferred :: min_gc
  procedure(name_emf_solver), deferred :: name
  procedure(read_input_emf_solver), deferred :: read_input
  procedure(test_stability_emf_solver), deferred :: test_stability
  procedure :: cleanup => cleanup_emf_solver
end type

abstract interface
  subroutine init_emf_solver( this, gix_pos, coords )
    import t_emf_solver
    class( t_emf_solver ), intent(inout) :: this
    integer, intent(in), dimension(:) :: gix_pos
    integer, intent(in) :: coords
  end subroutine
end interface

abstract interface
  subroutine advance_emf_solver( this, e, b, jay, dt, bnd_con )
    import t_emf_solver, t_vdf, p_double, t_emf_bound
    class( t_emf_solver ), intent(inout) :: this
    class( t_vdf ), intent(inout) :: e, b, jay
    real(p_double), intent(in) :: dt
    class( t_emf_bound ), intent(inout) :: bnd_con
  end subroutine
end interface

abstract interface
  subroutine min_gc_emf_solver( this, min_gc )
    import t_emf_solver, p_x_dim
    class( t_emf_solver ), intent(in) :: this
    integer, dimension(:,:), intent(inout) :: min_gc
  end subroutine
end interface

abstract interface
  function name_emf_solver( this )
    import t_emf_solver, p_maxlen_emf_solver_name
    class( t_emf_solver ), intent(in) :: this
    character( len = p_maxlen_emf_solver_name ) :: name_emf_solver
  end function name_emf_solver
end interface

abstract interface
  subroutine read_input_emf_solver( this, input_file, dx, dt )
    import t_emf_solver, t_input_file, p_double
    class( t_emf_solver ), intent(inout) :: this
    class( t_input_file ), intent(inout) :: input_file
    real(p_double), dimension(:), intent(in) :: dx
    real(p_double), intent(in) :: dt
  end subroutine read_input_emf_solver
end interface

abstract interface
  subroutine test_stability_emf_solver( this, dt, dx )
    import t_emf_solver, p_double
    class( t_emf_solver ), intent(in) :: this
    real(p_double), intent(in) :: dt
    real(p_double), dimension(:), intent(in) :: dx
  end subroutine test_stability_emf_solver
end interface

! Most solvers don't use the cleanup feature, so a simple subroutine is defined
! in os-emf-solver that doesn't do anything
interface
  subroutine cleanup_emf_solver( this )
    import t_emf_solver
    class( t_emf_solver ), intent(inout) :: this
  end subroutine cleanup_emf_solver
end interface
# 99 "emf/os-emf-define.f03" 2

!-------------------------------------------------------------------------------
! t_emf_gridval class definition
!-------------------------------------------------------------------------------
! Parameters for definining EMF grid field values (used by initial and external
! fields)
!-------------------------------------------------------------------------------

type :: t_emf_gridval

  ! Type of grid values to set (none, uniform, math or dipole)
  integer, dimension(p_f_dim) :: type_b
  integer, dimension(p_f_dim) :: type_e

  ! variables for uniform fields
  real(p_k_fld), dimension(p_f_dim) :: uniform_b0 = 0.0_p_k_fld
  real(p_k_fld), dimension(p_f_dim) :: uniform_e0 = 0.0_p_k_fld

  ! variables for math function fields

  logical :: dynamic = .false. ! the field values vary over time
  character(len = p_max_expr_len), dimension(p_f_dim) :: mfunc_expr_b = " "
  type(t_fparser), dimension(p_f_dim) :: mfunc_b

  character(len = p_max_expr_len), dimension(p_f_dim) :: mfunc_expr_e = " "
  type(t_fparser), dimension(p_f_dim) :: mfunc_e

  ! variables for dipole fields
  real(p_k_fld), dimension(p_f_dim) :: dipole_b_m
  real(p_k_fld), dimension(p_x_dim) :: dipole_b_x0
  real(p_k_fld) :: dipole_b_r0

  real(p_k_fld), dimension(p_f_dim) :: dipole_e_p
  real(p_k_fld), dimension(p_x_dim) :: dipole_e_x0
  real(p_k_fld) :: dipole_e_r0

  ! variables for python function fields
  character(len = p_max_expr_len) :: py_mod = "", py_func = ""

  integer :: interpolation

  ! If flds should be initialized when performing the moving window
  logical :: init_move_window = .true.

  contains

  procedure :: set_fld_values => set_fld_values_emf_gridval
  procedure :: write_checkpoint => write_checkpoint_emf_gridval
  procedure :: restart_read => restart_read_emf_gridval
  procedure :: setup => setup_emf_gridval
  procedure :: cleanup => cleanup_emf_gridval

end type t_emf_gridval

interface
  subroutine set_fld_values_emf_gridval( this, e, b, g_space, nx_p_min, t, lbin, ubin )
    import t_emf_gridval, t_vdf, t_space, p_double
    class(t_emf_gridval), intent( inout ) :: this
    type(t_vdf), intent(inout) :: e, b
    type( t_space ), intent(in) :: g_space
    integer, dimension(:), intent(in) :: nx_p_min
    real( p_double ), intent(in) :: t
    integer, dimension(:), intent(in), optional :: lbin
    integer, dimension(:), intent(in), optional :: ubin
  end subroutine
end interface

interface
  subroutine write_checkpoint_emf_gridval( this, restart_handle )
    import t_emf_gridval, t_restart_handle
    class( t_emf_gridval ), intent(in) :: this
    type( t_restart_handle ), intent(inout) :: restart_handle
  end subroutine
end interface

interface
  subroutine restart_read_emf_gridval( this, restart_handle )
    import t_emf_gridval, t_restart_handle
    class( t_emf_gridval ), intent(inout) :: this
    type( t_restart_handle ), intent(in) :: restart_handle
  end subroutine
end interface

interface
  subroutine setup_emf_gridval( this, dynamic, interpolation, restart, restart_handle )
    import t_emf_gridval, t_restart_handle
    class( t_emf_gridval ), intent(inout) :: this
    logical, intent(in) :: dynamic
    integer, intent(in) :: interpolation
    logical, intent(in) :: restart
    type( t_restart_handle ), intent(in) :: restart_handle
  end subroutine
end interface

interface
  subroutine cleanup_emf_gridval( this )
    import t_emf_gridval
    class( t_emf_gridval ), intent(inout) :: this
  end subroutine
end interface

!-------------------------------------------------------------------------------
! t_emf class definition
!-------------------------------------------------------------------------------

type :: t_emf

  ! time step that the field data held within this t_emf was calculated.
  ! Used by external calculation accelerators (e.g. GPUs)
  ! to signal the need for the field data to be copied from the GPU back to Osiris
  ! for processing or diagnostic output.
  integer :: n_current_emf = 0

  ! electromagnetic field values on the grid:
  type( t_vdf ) :: b, e

  ! center fields on the corner of the cell for particle interpolation
  logical :: part_grid_center

  ! Initial values for EMF
  class( t_emf_gridval ), pointer :: init_emf => null()

  ! em field values for particle interpolation
  type( t_vdf ), pointer :: b_part => null(), e_part => null()
  logical :: part_fld_alloc

  ! type of field solver to use
  class( t_emf_solver ), pointer :: solver => null()
  procedure(new_solver_emf_int), pointer :: new_solver => null()

  ! switches that turns on external fields
  integer :: ext_fld = p_extfld_none

  ! external electromagnetic field values
  type( t_vdf ) :: ext_b, ext_e

  ! external magnetic and electric fields
  class( t_emf_gridval ), pointer :: ext_emf => null()

  ! boundary conditions for the electromagnetic field
  class( t_emf_bound ), pointer :: bnd_con => null()

  ! spatial smoothing for electric and magnetic fields
  ! to be applied before advancing particles
  integer :: smooth_type, smooth_niter, smooth_nmax
  type( t_smooth ) :: smooth

  ! diagnostic set up for fields
  class( t_diag_emf ), pointer :: diag => null()

  ! psi diagnostic
  type( t_emf_psi ) :: psi

  ! coordinates type
  integer :: coordinates

  ! Position on global grid
  integer, dimension(p_max_dim) :: gix_pos

contains

  procedure :: allocate_objs => allocate_objs_emf
  procedure :: read_input => read_input_emf
  procedure :: init => init_emf
  procedure :: advance => advance_emf
  procedure :: report => report_diag_emf
  procedure :: report_energy => report_energy_emf
  procedure :: fill_data => fill_data_emf

  procedure :: solver_type_class_emf
  procedure :: solver_type_name_emf
  generic :: f_solver_type => solver_type_class_emf, solver_type_name_emf

  procedure :: cleanup => cleanup_emf
  procedure :: write_checkpoint => write_checkpoint_emf
  procedure :: restart_read => restart_read_emf
  procedure :: move_window => move_window_emf

  procedure :: update_boundary => update_boundary_emf
  procedure :: update_particle_fld => update_particle_fld_emf

  procedure :: list_algorithm => list_algorithm_emf
  procedure :: min_gc => min_gc_emf
  procedure :: test_stability => test_stability_emf
  procedure :: dx => dx_emf
  procedure :: if_charge_cons => if_charge_cons_emf
  procedure :: bc_type => bc_type_emf
  procedure :: get_diag_buffer_size => get_diag_buffer_size_emf

  procedure :: reshape => reshape_emf

end type t_emf

interface
  subroutine new_solver_emf_int( this, name )
  import t_emf
  class( t_emf ), intent(inout) :: this
  character( len = * ), intent(in) :: name
  end subroutine
end interface

interface
  subroutine allocate_objs_emf( this )
    import t_emf
    class( t_emf ), intent(inout) :: this
  end subroutine
end interface

interface
  subroutine read_input_emf( this, input_file, periodic, if_move, grid, dx, dt, gamma )

    import t_emf, t_input_file, t_grid, p_double

    class( t_emf ), intent(inout) :: this
    class( t_input_file ), intent(inout) :: input_file
    logical, dimension(:), intent(in) :: periodic, if_move
    class( t_grid ), intent(in) :: grid
    real(p_double), dimension(:), intent(in) :: dx
    real(p_double), intent(in) :: dt
    real(p_double), intent(in) :: gamma
  end subroutine
end interface

interface
  subroutine init_emf( this, part_grid_center, part_interpolation, &
                        g_space, grid, gc_min, dx, tstep, tmin, tmax, &
                        no_co, send_msg, recv_msg, restart, restart_handle, sim_options )

    import t_emf, t_space, p_double, t_node_conf, t_vdf_msg
    import t_restart_handle, t_options, t_grid, t_time_step

    class( t_emf ), intent( inout ), target :: this

    logical, intent(in) :: part_grid_center
    integer, intent(in) :: part_interpolation
    type( t_space ), intent(in) :: g_space
    class( t_grid ), intent(in) :: grid
    integer, dimension(:,:), intent(in) :: gc_min
    type( t_time_step ), intent(in) :: tstep
    real(p_double), intent(in) :: tmin, tmax
    real( p_double ), dimension(:), intent(in) :: dx

    class( t_node_conf ), intent(in), target :: no_co
    type(t_vdf_msg), dimension(2), intent(inout) :: send_msg, recv_msg

    logical, intent(in) :: restart
    type( t_restart_handle ), intent(in) :: restart_handle

    type( t_options ), intent(in) :: sim_options

  end subroutine
end interface

interface

  subroutine report_diag_emf( this, g_space, grid, no_co, tstep, t, send_msg, recv_msg )
  import t_emf, t_space, t_grid, t_node_conf, t_time_step, p_double, t_vdf_msg
  class( t_emf ), intent(inout) :: this
  type( t_space ), intent(in) :: g_space
  class( t_grid ), intent(in) :: grid
  class( t_node_conf ), intent(in) :: no_co
  type( t_time_step ), intent(in) :: tstep
  real(p_double), intent(in) :: t
  type(t_vdf_msg), dimension(2), intent(inout) :: send_msg, recv_msg

  end subroutine
end interface

interface
  subroutine report_energy_emf( this, no_co, tstep, t, tmin )
  import t_emf, t_node_conf, t_time_step, p_double
  class( t_emf ), intent(inout) :: this
  class( t_node_conf ), intent(in) :: no_co
  type( t_time_step ), intent(in) :: tstep
  real(p_double), intent(in) :: t, tmin

  end subroutine
end interface

interface
  subroutine fill_data_emf(this, tstep)
  import t_emf, t_time_step
  class( t_emf ), intent(inout) :: this
  type( t_time_step ), intent(in) :: tstep
  end subroutine
end interface

interface
  function solver_type_class_emf( this )
  import t_emf
  integer :: solver_type_class_emf
  class( t_emf ), intent(in) :: this
  end function
end interface

interface
  function solver_type_name_emf( this, name )
  import t_emf
  integer :: solver_type_name_emf
  class( t_emf ), intent(in) :: this
  character( len = * ), intent(in) :: name
  end function
end interface

interface
  subroutine cleanup_emf( this )
  import t_emf
  class( t_emf ), intent(inout) :: this
  end subroutine
end interface

interface
  subroutine write_checkpoint_emf( this, restart_handle )
  import t_emf, t_restart_handle
  class( t_emf ), intent(in) :: this
  type( t_restart_handle ), intent(inout) :: restart_handle
  end subroutine
end interface

interface
  subroutine restart_read_emf( this, restart_handle )
  import t_emf, t_restart_handle
  class( t_emf ), intent(inout) :: this
  type( t_restart_handle ), intent(in) :: restart_handle
  end subroutine
end interface

interface
  subroutine move_window_emf( this, g_space, nx_p_min, need_fld_val )
  import t_emf, t_space
  class( t_emf ), intent( inout ) :: this
  type( t_space ), intent(in) :: g_space
  integer, dimension(:), intent(in) :: nx_p_min
  logical, intent(in) :: need_fld_val
  end subroutine
end interface

interface
  subroutine update_boundary_emf( this, g_space, no_co, send_msg, recv_msg )
    import t_emf, t_space, t_node_conf, t_vdf_msg
    class( t_emf ), intent( inout ) :: this
    type( t_space ), intent(in) :: g_space
    class( t_node_conf ), intent(in) :: no_co
    type(t_vdf_msg), dimension(2), intent(inout) :: send_msg, recv_msg
  end subroutine update_boundary_emf
end interface

interface
  subroutine update_particle_fld_emf( this, n, t, dt, space, nx_p_min )
    import t_emf, t_space, p_double
    class( t_emf ), intent( inout ) :: this
    integer, intent(in) :: n
    real( p_double ), intent(in) :: t, dt
    type( t_space ), intent(in) :: space
    integer, dimension(:), intent(in) :: nx_p_min
  end subroutine update_particle_fld_emf
end interface

interface
  subroutine reshape_emf( this, old_lb, new_lb, no_co, send_msg, recv_msg )
    import t_emf, t_grid, t_node_conf, t_vdf_msg
    class(t_emf), intent(inout), target :: this
    class( t_grid ), intent(in) :: old_lb, new_lb
    class( t_node_conf ), intent(in) :: no_co
    type(t_vdf_msg), dimension(2), intent(inout) :: send_msg, recv_msg
  end subroutine reshape_emf
end interface


! -------------------------------------------------------------------------------------------------

public :: t_diag_emf, t_emf_gridval, t_emf, t_emf_psi

public :: fsolverev, fsmoothev, getpsi_ev

public :: p_emfsmooth_none, p_emfsmooth_stand, p_emfsmooth_nci, p_emfsmooth_local

public :: p_emf_none, p_emf_uniform, p_emf_math, p_emf_dipole, p_emf_python
public :: p_emf_hc_e, p_emf_hc_b

public :: p_extfld_none, p_extfld_static, p_extfld_dynamic

contains

!-----------------------------------------------------------------------------------------
! Advance the EM fields 1 iteration
!-----------------------------------------------------------------------------------------
subroutine advance_emf( this, jay, dt, send_msg, recv_msg )

  implicit none

  class( t_emf ), intent( inout ) :: this
  ! We must use the t_current object and not a vdf because some subclasses of t_emf
  ! use subclasses of t_current
  class( t_current ), intent(inout) :: jay
  real(p_double), intent(in) :: dt
  type(t_vdf_msg), dimension(2), intent(inout) :: send_msg, recv_msg

  call begin_event(fsolverev)
  call this % solver % advance( this % e, this % b, jay % pf(1), dt, this % bnd_con )
  call end_event(fsolverev)

end subroutine

!-----------------------------------------------------------------------------------------
! returns the number of guard cells required by this
! object. Does not require setup
!-----------------------------------------------------------------------------------------
function min_gc_emf( this )

  implicit none

  class( t_emf ), intent(in) :: this
  integer, dimension(2,p_x_dim) :: min_gc_emf

  integer, dimension(2,p_x_dim) :: gc_num_temp
  integer :: i

  ! get number of guard cells required by field solver
  call this % solver % min_gc(min_gc_emf)

  ! Get number of guard cells required by smoothing and choose the maximum
  gc_num_temp(p_lower,:) = smooth_order(this%smooth)
  gc_num_temp(p_upper,:) = gc_num_temp(p_lower,:)+1

  do i = 1, p_x_dim
    min_gc_emf(p_lower,i)=max(min_gc_emf(p_lower,i),gc_num_temp(p_lower,i))
    min_gc_emf(p_upper,i)=max(min_gc_emf(p_upper,i),gc_num_temp(p_upper,i))
  enddo

end function min_gc_emf

!-----------------------------------------------------------------------------------------
! Test resolution and time step for solver algorithm stability
!-----------------------------------------------------------------------------------------
subroutine test_stability_emf( this, dt, dx )

  implicit none

  class( t_emf ), intent(in) :: this
  real(p_double), intent(in) :: dt
  real(p_double), dimension(:), intent(in) :: dx

  ! each solver has a different stability criteria
  call this % solver % test_stability( dt, dx )

end subroutine test_stability_emf

!-----------------------------------------------------------------------------------------
! Printout the algorithm for the EMF object
!-----------------------------------------------------------------------------------------
subroutine list_algorithm_emf( this )

  implicit none

  class( t_emf ), intent(in) :: this

  integer :: i

  print *, ' '
  print *, 'Field Solver:'
  print *, '- ', trim( this % solver % name() ), ' solver'

  print *, '- Guard Cells:'
  do i = 1, p_x_dim
    print '(A,I0,A,I0,A,I0,A)', '     x', i, ' : [', this%b%gc_num_( p_lower, i ), &
                                                ', ', this%b%gc_num_( p_upper, i ) , ']'
  enddo

  ! external fields
  select case ( this%ext_fld )
    case (p_extfld_none)
      print *, '- Using self generated fields only'
    case (p_extfld_static)
      print *, '- Using static external field source'
    case (p_extfld_dynamic)
      print *, '- Using dynamic external field source'
  end select

  ! field smoothing
  select case ( this%smooth_type )
    case ( p_emfsmooth_none )
      print *, '- No field smoothing'
    case ( p_emfsmooth_stand )
      print *, '- Normal field smoothing, particle interpolation fields only'
    case ( p_emfsmooth_nci )
      print *, '- NCI filtering, particle interpolation fields only'
    case ( p_emfsmooth_local )
      print *, '- Direct filtering of EM fields'
  end select

  ! Boundary conditions
  call this % bnd_con % list_algorithm( )

end subroutine list_algorithm_emf

function dx_emf( this )

  implicit none

  class(t_emf), intent(in) :: this
  real(p_double) , dimension(p_x_dim) :: dx_emf

  integer :: i

  ! executable statements
  do i = 1, this%e%x_dim_
    dx_emf(i) = this%e%dx_(i)
  enddo

end function dx_emf

subroutine get_diag_buffer_size_emf( this, gnx, diag_buffer_size )

  implicit none

  class( t_emf ), intent(in) :: this
  integer, dimension(:), intent(in) :: gnx
  integer, intent(inout) :: diag_buffer_size

  integer, dimension(3) :: ln_avg

  integer :: i, bsize

  ln_avg = n_avg( this%diag%reports, p_x_dim )

  if (ln_avg(1) > 0) then
    bsize = gnx(1) / ln_avg(1)
    do i = 2, p_x_dim
      bsize = bsize * ( gnx(i) / ln_avg(i) )
    enddo

    ! we need 2 buffers (why?)
    bsize = 2*bsize
    if ( bsize > diag_buffer_size ) diag_buffer_size = bsize
  endif

end subroutine get_diag_buffer_size_emf

function if_charge_cons_emf( this, tstep )

  use m_time_step

  implicit none

  class( t_emf ), intent(in) :: this
  type( t_time_step ), intent(in) :: tstep

  logical :: if_charge_cons_emf

  if_charge_cons_emf = test_if_report( tstep, this%diag%ndump_fac_charge_cons )

end function if_charge_cons_emf

!-----------------------------------------------------------------------------------------
! gives the boundary type description variable for this emf
!-----------------------------------------------------------------------------------------
function bc_type_emf( this )

  implicit none

  class(t_emf), intent(in) :: this
  integer, dimension(2,p_x_dim) :: bc_type_emf

  bc_type_emf = this % bnd_con % type( 1:2, 1:p_x_dim )

end function bc_type_emf


end module
