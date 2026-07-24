#include "os-config.h"
#include "os-preprocess.fpp"

module m_emf_cyl_modes

#include "memory/memory.h"

use m_system
use m_parameters
use m_cyl_modes
use m_emf_define
use m_vdf_define,  only: t_vdf, t_vdf_report
use m_vdf_report,  only: if_report, report_vdf
use m_vdf_comm,    only: irecv, isend, irecv_wait, isend_wait, t_vdf_msg, update_periodic_vdf
use m_emf_diag,    only: t_diag_emf, diag_emf_ev, diag_emf_boost_ev, diag_emf_boost_extract_ev
use m_space,       only: t_space, x_dim, nx_move, if_move, xmin, xmax
use m_node_conf,   only: t_node_conf, neighbor, periodic, p_sum, reduce_array, root
use m_grid_define, only: t_grid
use m_restart,     only: t_restart_handle, restart_io_write, restart_io_read
use m_current_define, only: t_current
use m_logprof,     only: begin_event, end_event, create_event
use m_grid_cyl_modes, only: t_grid_cyl_modes
use m_current_cyl_modes, only: t_current_cyl_modes
use m_emf_solver_cyl_modes_yee, only: t_emf_solver_cyl_modes_yee
use m_emf_solver_cyl_modes_fei, only: t_emf_solver_cyl_modes_fei
use m_time_step,   only: t_time_step, test_if_report, n, dt, ndump
use m_emf_bound,   only: read_nml, type
use m_emf_bound_cyl_modes
use m_input_file,  only: t_input_file, get_namelist, disp_out
use m_emf_solver_cyl_modes
use m_vdf_math,    only: total, add
use m_math,        only: pi
use m_fparser,     only: p_max_expr_len, t_fparser, p_k_fparse, eval, cleanup, setup
use m_vdf_smooth,  only: read_nml, smooth
use m_emf_ncifilter, only: filter_nci
use m_emf_gridcenter, only: grid_center
use m_boosted_diag_cyl_modes, only : t_boosted_diag_cyl_modes
use m_callpy

use m_vdf_comm
use m_node_conf


use m_logprof


implicit none

private

! string to id restart data
character(len=*), parameter :: p_emf_cyl_modes_rst_id = "emf_cyl_modes rst data - 0x0001"
character(len=*), parameter :: p_emf_gridval_cyl_modes_rst_id = "emf_gridval_cyl_modes rst data - 0x0000"

!-------------------------------------------------------------------------------
! t_emf_psi_cyl_modes class definition
!-------------------------------------------------------------------------------
type :: t_emf_psi_cyl_modes

  ! vdf holding psi values
  type( t_cyl_modes ) :: data

  ! Iteration when psi was last calculated
  integer :: n = -1

  ! Communicator used for parallel psi calculations
  integer :: comm = mpi_comm_null

end type t_emf_psi_cyl_modes



!-------------------------------------------------------------------------------
! t_emf_cyl_modes class definition
!-------------------------------------------------------------------------------

type, extends( t_emf ) :: t_emf_cyl_modes

  integer :: n_cyl_modes = -1

  type(t_emf_psi_cyl_modes) :: psi_cyl_m

  type( t_cyl_modes ) :: b_cyl_m
  type( t_cyl_modes ) :: e_cyl_m

  type( t_cyl_modes ) :: ext_b_cyl_m
  type( t_cyl_modes ) :: ext_e_cyl_m

  type( t_cyl_modes ), pointer :: b_part_cyl_m => null()
  type( t_cyl_modes ), pointer :: e_part_cyl_m => null()

contains

  procedure :: allocate_objs => allocate_objs_emf_cyl_modes
  procedure :: read_input => read_input_emf_cyl_modes
  procedure :: init => init_emf_cyl_modes
  procedure :: advance => advance_emf_cyl_modes
  procedure :: report => report_diag_emf_cyl_modes
  procedure :: report_energy => report_energy_emf_cyl_modes
  procedure :: cleanup => cleanup_emf_cyl_modes
  procedure :: write_checkpoint => write_checkpoint_emf_cyl_modes
  procedure :: restart_read_cm => restart_read_cm_emf_cyl_modes
  procedure :: move_window => move_window_emf_cyl_modes
  procedure :: update_boundary => update_boundary_emf_cyl_modes

  procedure :: solver_type_class_emf => solver_type_class_emf_cyl_modes
  procedure :: solver_type_name_emf => solver_type_name_emf_cyl_modes

  procedure :: update_particle_fld => update_particle_fld_emf_cyl_modes

  procedure :: reshape => reshape_emf_cyl_modes

end type t_emf_cyl_modes

type, extends( t_diag_emf ) :: t_diag_emf_cyl_modes

  integer :: n_cyl_modes = -1

contains

  procedure :: allocate_objs       => allocate_objs_diag_emf_cyl_modes
  procedure :: init                => init_diag_emf_cyl_modes
  procedure :: avail_report_quants => avail_report_quants_emf_cyl_modes
  procedure :: init_report_quants  => init_report_quants_emf_cyl_modes

end type t_diag_emf_cyl_modes

type, extends( t_emf_gridval ) :: t_emf_gridval_cyl_modes

  integer :: n_cyl_modes = -1

  ! variables for uniform fields in cyl mode decomposition
  real(p_k_fld), dimension(p_max_cyl_modes,p_f_dim) :: uniform_e0_re
  real(p_k_fld), dimension(p_max_cyl_modes,p_f_dim) :: uniform_e0_im
  real(p_k_fld), dimension(p_max_cyl_modes,p_f_dim) :: uniform_b0_re
  real(p_k_fld), dimension(p_max_cyl_modes,p_f_dim) :: uniform_b0_im

  ! variables for math function fields in cyl mode decomposition
  character(len = p_max_expr_len), dimension(p_max_cyl_modes,p_f_dim) :: mfunc_expr_b_re = " "
  type(t_fparser),                 dimension(p_max_cyl_modes,p_f_dim) :: mfunc_b_re
  character(len = p_max_expr_len), dimension(p_max_cyl_modes,p_f_dim) :: mfunc_expr_b_im = " "
  type(t_fparser),                 dimension(p_max_cyl_modes,p_f_dim) :: mfunc_b_im

  character(len = p_max_expr_len), dimension(p_max_cyl_modes,p_f_dim) :: mfunc_expr_e_re = " "
  type(t_fparser),                 dimension(p_max_cyl_modes,p_f_dim) :: mfunc_e_re
  character(len = p_max_expr_len), dimension(p_max_cyl_modes,p_f_dim) :: mfunc_expr_e_im = " "
  type(t_fparser),                 dimension(p_max_cyl_modes,p_f_dim) :: mfunc_e_im

contains

  procedure :: set_fld_values_cm => set_fld_values_cm_emf_gridval_cyl_modes
  procedure :: write_checkpoint  => write_checkpoint_emf_gridval_cyl_modes
  procedure :: restart_read      => restart_read_emf_gridval_cyl_modes
  procedure :: setup             => setup_emf_gridval_cyl_modes
  procedure :: cleanup           => cleanup_emf_gridval_cyl_modes

end type t_emf_gridval_cyl_modes

! diagnostic quantities specific to the quasi-3D geometry
character(len=13), dimension(33), parameter :: p_report_quants = &
  (/ 'e1_cyl_m     ', 'e2_cyl_m     ', 'e3_cyl_m     ', &
     'b1_cyl_m     ', 'b2_cyl_m     ', 'b3_cyl_m     ', &
     'ext_e1_cyl_m ', 'ext_e2_cyl_m ', 'ext_e3_cyl_m ', &
     'ext_b1_cyl_m ', 'ext_b2_cyl_m ', 'ext_b3_cyl_m ', &
     'part_e1_cyl_m', 'part_e2_cyl_m', 'part_e3_cyl_m', &
     'part_b1_cyl_m', 'part_b2_cyl_m', 'part_b3_cyl_m', &
     'ene_e1_cyl_m ', 'ene_e2_cyl_m ', 'ene_e3_cyl_m ', &
     'ene_b1_cyl_m ', 'ene_b2_cyl_m ', 'ene_b3_cyl_m ', &
     'ene_e_cyl_m  ', 'ene_b_cyl_m  ', 'ene_emf_cyl_m', &
     'div_e_cyl_m  ', 'div_b_cyl_m  ', 'psi_cyl_m    ', &
     's1_int_cyl_m ', 's2_int_cyl_m ', 's3_int_cyl_m ' /)

character(len=14), dimension(6), parameter :: p_boosted_report_quants = &
   (/ 'e1_boost_cyl_m', 'e2_boost_cyl_m', 'e3_boost_cyl_m', &
      'b1_boost_cyl_m', 'b2_boost_cyl_m', 'b3_boost_cyl_m' /)


! cyl_modes specific diagnostics
integer, parameter :: p_e1_cyl_m = 1, p_e2_cyl_m = 2, p_e3_cyl_m = 3
integer, parameter :: p_b1_cyl_m = 4, p_b2_cyl_m = 5, p_b3_cyl_m = 6
integer, parameter :: p_ext_e1_cyl_m = 7, p_ext_e2_cyl_m = 8, p_ext_e3_cyl_m = 9
integer, parameter :: p_ext_b1_cyl_m = 10, p_ext_b2_cyl_m = 11, p_ext_b3_cyl_m = 12
integer, parameter :: p_part_e1_cyl_m = 13, p_part_e2_cyl_m = 14, p_part_e3_cyl_m = 15
integer, parameter :: p_part_b1_cyl_m = 16, p_part_b2_cyl_m = 17, p_part_b3_cyl_m = 18
integer, parameter :: p_ene_e1_cyl_m = 19, p_ene_e2_cyl_m = 20, p_ene_e3_cyl_m = 21
integer, parameter :: p_ene_b1_cyl_m = 22, p_ene_b2_cyl_m = 23, p_ene_b3_cyl_m = 24
integer, parameter :: p_ene_e_cyl_m = 25, p_ene_b_cyl_m = 26, p_ene_emf_cyl_m = 27
integer, parameter :: p_div_e_cyl_m = 28, p_div_b_cyl_m = 29, p_psi_cyl_m = 30
integer, parameter :: p_s1_int_cyl_m = 31, p_s2_int_cyl_m = 32, p_s3_int_cyl_m = 33

! exported symbols
public :: t_emf_cyl_modes, t_emf_psi_cyl_modes
public :: t_diag_emf_cyl_modes

!-----------------------------------------------------------------------------------------
contains

!-----------------------------------------------------------------------------------------
!       Allocate any objects contained within the EMF (t_emf) object
!-----------------------------------------------------------------------------------------
subroutine allocate_objs_emf_cyl_modes( this )

  implicit none

  class( t_emf_cyl_modes ), intent(inout) :: this

  if ( .not. associated( this%new_solver ) ) then
    this%new_solver => new_solver_emf_cyl_modes
  endif

  if ( .not. associated( this%bnd_con ) ) then
    allocate( t_emf_bound_cyl_modes :: this%bnd_con )
  endif

  if ( .not. associated( this%diag ) ) then
    allocate( t_diag_emf_cyl_modes :: this%diag )
  endif

  if ( .not. associated( this%init_emf ) ) then
    allocate( t_emf_gridval_cyl_modes :: this%init_emf )
  endif

  if ( .not. associated( this%ext_emf ) ) then
    allocate( t_emf_gridval_cyl_modes :: this%ext_emf )
  endif

  ! call superclass
  call this % t_emf % allocate_objs()

end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine allocate_objs_diag_emf_cyl_modes( this )

  implicit none

  class( t_diag_emf_cyl_modes ), intent( inout )  ::  this

  if ( .not. associated( this%boosted_diag ) ) then
    allocate( t_boosted_diag_cyl_modes :: this%boosted_diag )
  endif

  ! call superclass
  call this % t_diag_emf % allocate_objs()

end subroutine allocate_objs_diag_emf_cyl_modes
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Read input for emf_cyl_modes object
!   - this routine mainly calls the superclass (t_emf) routine
!-----------------------------------------------------------------------------------------
subroutine read_input_emf_cyl_modes( this, input_file, periodic, if_move, grid, dx, dt, gamma )

  implicit none

  class( t_emf_cyl_modes ), intent(inout) :: this
  class( t_input_file ), intent(inout) :: input_file
  logical, dimension(:), intent(in) :: periodic, if_move
  class( t_grid ), intent(in)  :: grid
  real(p_double), dimension(:), intent(in) :: dx
  real(p_double), intent(in) :: dt
  real(p_double), intent(in) :: gamma

  character(len=20) :: solver

  ! smooth
  character(len=20) :: smooth_type
  integer :: smooth_niter, smooth_nmax

  character(len=20) :: ext_fld

  ! Initial fields
  character(len=20), dimension(p_f_dim) :: type_init_b ! type of initial b-field
  character(len=20), dimension(p_f_dim) :: type_init_e ! type of initial e-field

  ! magnitude of initial e and b field cylindrical modes
  real(p_k_fld), dimension(0:p_max_cyl_modes, p_f_dim) :: init_b0_re
  real(p_k_fld), dimension(0:p_max_cyl_modes, p_f_dim) :: init_e0_re
  real(p_k_fld), dimension(1:p_max_cyl_modes, p_f_dim) :: init_b0_im
  real(p_k_fld), dimension(1:p_max_cyl_modes, p_f_dim) :: init_e0_im

  character(len = p_max_expr_len), dimension(0:p_max_cyl_modes,p_f_dim) :: init_b_re_mfunc
  character(len = p_max_expr_len), dimension(0:p_max_cyl_modes,p_f_dim) :: init_e_re_mfunc
  character(len = p_max_expr_len), dimension(1:p_max_cyl_modes,p_f_dim) :: init_b_im_mfunc
  character(len = p_max_expr_len), dimension(1:p_max_cyl_modes,p_f_dim) :: init_e_im_mfunc

  character(len = p_max_expr_len) :: init_py_mod, init_py_func

  logical :: init_move_window

  ! external fields
  character(len=20), dimension(p_f_dim) :: type_ext_b ! type of external b-field
  character(len=20), dimension(p_f_dim) :: type_ext_e ! type of external e-field

  ! magnitude of external e and b field cylindrical modes
  real(p_k_fld), dimension(0:p_max_cyl_modes, p_f_dim) :: ext_b0_re
  real(p_k_fld), dimension(0:p_max_cyl_modes, p_f_dim) :: ext_e0_re
  real(p_k_fld), dimension(1:p_max_cyl_modes, p_f_dim) :: ext_b0_im
  real(p_k_fld), dimension(1:p_max_cyl_modes, p_f_dim) :: ext_e0_im

  character(len = p_max_expr_len), dimension(0:p_max_cyl_modes,p_f_dim) :: ext_b_re_mfunc
  character(len = p_max_expr_len), dimension(0:p_max_cyl_modes,p_f_dim) :: ext_e_re_mfunc
  character(len = p_max_expr_len), dimension(1:p_max_cyl_modes,p_f_dim) :: ext_b_im_mfunc
  character(len = p_max_expr_len), dimension(1:p_max_cyl_modes,p_f_dim) :: ext_e_im_mfunc

  character(len = p_max_expr_len) :: ext_py_mod, ext_py_func

  ! dipole fields not yet supported in cyl_modes

  namelist /nl_el_mag_fld/ solver, smooth_type, smooth_niter, smooth_nmax, &
                           type_init_b, type_init_e, &
                           init_b0_re, init_e0_re, init_b0_im, init_e0_im, &
                           init_b_re_mfunc, init_e_re_mfunc, &
                           init_b_im_mfunc, init_e_im_mfunc, &
                           init_py_mod, init_py_func, init_move_window, &
                           ext_fld, type_ext_b, type_ext_e, &
                           ext_b0_re, ext_e0_re, ext_b0_im, ext_e0_im, &
                           ext_b_re_mfunc, ext_e_re_mfunc, &
                           ext_b_im_mfunc, ext_e_im_mfunc, &
                           ext_py_mod, ext_py_func

  integer :: ierr, i, j

  solver = "yee"

  smooth_type = "none"
  smooth_niter = 1
  smooth_nmax = -1

  ! Initial fields
  type_init_b = "uniform"
  type_init_e = "uniform"
  init_b0_re  = 0.0_p_k_fld
  init_b0_im  = 0.0_p_k_fld
  init_e0_re  = 0.0_p_k_fld
  init_e0_im  = 0.0_p_k_fld

  init_b_re_mfunc = "NO_FUNCTION_SUPPLIED!"
  init_b_im_mfunc = "NO_FUNCTION_SUPPLIED!"
  init_e_re_mfunc = "NO_FUNCTION_SUPPLIED!"
  init_e_im_mfunc = "NO_FUNCTION_SUPPLIED!"

  init_py_mod = "NO_PYTHON_MODULE_SUPPLIED!"
  init_py_mod = "NO_PYTHON_FUNCTION_SUPPLIED!"

  init_move_window = .true.

  ! external fields
  ext_fld    = "none"
  type_ext_b = "uniform"
  type_ext_e = "uniform"
  ext_b0_re  = 0.0_p_k_fld
  ext_b0_im  = 0.0_p_k_fld
  ext_e0_re  = 0.0_p_k_fld
  ext_e0_im  = 0.0_p_k_fld

  ext_b_re_mfunc = "NO_FUNCTION_SUPPLIED!"
  ext_b_im_mfunc = "NO_FUNCTION_SUPPLIED!"
  ext_e_re_mfunc = "NO_FUNCTION_SUPPLIED!"
  ext_e_im_mfunc = "NO_FUNCTION_SUPPLIED!"

  ext_py_mod = "NO_PYTHON_MODULE_SUPPLIED!"
  ext_py_mod = "NO_PYTHON_FUNCTION_SUPPLIED!"

  ! Get namelist text from input file
  call get_namelist( input_file, "nl_el_mag_fld", ierr )

  if (ierr == 0) then
    read (input_file%nml_text, nml = nl_el_mag_fld, iostat = ierr)
    if (ierr /= 0) then
      print *, "Error reading emf parameters"
      print *, "aborting..."
      stop
    endif
  else
    if (ierr < 0) then
      print *, "Error reading emf parameters"
      print *, "aborting..."
      stop
    else
      if (disp_out(input_file)) then
        SCR_ROOT(" - emf parameters missing, using default")
      endif
    endif
  endif

  ! allocate objects emf cyl modes (this calls superclass allocate_objs as well)
  call this % allocate_objs()

  ! setup cyl_modes specific data
  select type( grid ); class is ( t_grid_cyl_modes )
    this%n_cyl_modes = grid%n_cyl_modes
  end select

  ! Create new solver
  call this % new_solver( solver )
  if (.not. associated(this%solver)) then
    print *, "Error reading emf parameters"
    print *, "aborting..."
    stop
  endif

  ! smooth type
  select case(trim(smooth_type))
  case ( "none" )
    this%smooth_type = p_emfsmooth_none
  case ( "stand" )
    this%smooth_type = p_emfsmooth_stand
  case ( "nci" )
    this%smooth_type = p_emfsmooth_nci
  case ( "local" )

    ! Check zero niter
    if (smooth_niter == 0) then
      print *, "Error reading emf parameters"
      print *, "'smooth_niter' must be defined and > 0 for 'local' smooth"
      print *, "aborting..."
    stop

    else

      this%smooth_type = p_emfsmooth_local
      this%smooth_niter = smooth_niter
      this%smooth_nmax = smooth_nmax

    endif

  case default
    print *, "Error reading emf parameters"
    print *, "Invalid value for the smooth parameter"
    print *, "Available smooth types are 'none', 'stand', 'nci' and 'local'"
    print *, "aborting..."
    stop
  end select

  ! Initial fields
   do i=1, p_f_dim
    select case(trim(type_init_b(i)))
    case( "uniform" )
      this%init_emf%type_b(i) = p_emf_uniform
    case( "math func" )
      this%init_emf%type_b(i) = p_emf_math
    case( "python" )

      if ( .not. if_py_util() ) then
        if ( mpi_node() == 0 ) then
          write(0,*) "   Error: requested ""python"" for type_init_b, but code was"
          write(0,*) "   not compiled with python support. Please see the"
          write(0,*) "   PY_FCOMPILEFLAGS and PY_FLINKFLAGS flags in your config file."
          write(0,*) "   aborting..."
        endif
        stop
      endif

      ! Check if simulation is moving
      if (init_move_window .and. if_move(1)) then
        if ( mpi_node() == 0 ) then
          write(0,*) "(*warning*) Python routines for initializing fields will be"
          write(0,*) "            called each time the simulation domain moves."
          write(0,*) "            Consider setting init_move_window = .false."
        endif
      endif

      this%init_emf%type_b(i) = p_emf_python
      this%init_emf%py_mod = trim(init_py_mod)
      this%init_emf%py_func = trim(init_py_func)
    case default
      print *, "Error reading emf parameters"
      print *, "Invalid value for the type_init_b parameter"
      print *, "Available initial B-field types in cylindrical coords are"
      print *, "'uniform', 'math func' and 'python'"
      print *, "aborting..."
      stop
    end select

    select case(trim(type_init_e(i)))
    case( "uniform" )
      this%init_emf%type_e(i) = p_emf_uniform
    case( "math func" )
      this%init_emf%type_e(i) = p_emf_math
    case( "python" )

      if ( .not. if_py_util() ) then
        if ( mpi_node() == 0 ) then
          write(0,*) "   Error: requested ""python"" for type_init_e, but code was"
          write(0,*) "   not compiled with python support. Please see the"
          write(0,*) "   PY_FCOMPILEFLAGS and PY_FLINKFLAGS flags in your config file."
          write(0,*) "   aborting..."
        endif
        stop
      endif

      ! Check if simulation is moving
      if (init_move_window .and. if_move(1)) then
        if ( mpi_node() == 0 ) then
          write(0,*) "(*warning*) Python routines for initializing fields will be"
          write(0,*) "            called each time the simulation domain moves."
          write(0,*) "            Consider setting init_move_window = .false."
        endif
      endif

      this%init_emf%type_e(i) = p_emf_python
      this%init_emf%py_mod = trim(init_py_mod)
      this%init_emf%py_func = trim(init_py_func)
    case default
      print *, "Error reading emf parameters"
      print *, "Invalid value for the type_init_e parameter"
      print *, "Available initial E-field types in cylindrical coords are"
      print *, "'uniform', 'math func' and 'python'"
      print *, "aborting..."
      stop
    end select
  end do

  ! set initial e and b fields modal decomposition
  select type( init_emf => this%init_emf )
  class is( t_emf_gridval_cyl_modes )

    init_emf%n_cyl_modes = this%n_cyl_modes

    do j = 1, p_f_dim
      init_emf%uniform_e0(j) = init_e0_re(0,j)
      init_emf%uniform_b0(j) = init_b0_re(0,j)

      init_emf%mfunc_expr_e(j) = trim(init_e_re_mfunc(0,j))
      init_emf%mfunc_expr_b(j) = trim(init_b_re_mfunc(0,j))
    enddo

    ! store values for rest of modes in the special cyl modes arrays
    do j = 1, p_f_dim
      do i = 1, this%n_cyl_modes
        init_emf%uniform_e0_re(i,j) = init_e0_re(i,j)
        init_emf%uniform_e0_im(i,j) = init_e0_im(i,j)
        init_emf%uniform_b0_re(i,j) = init_b0_re(i,j)
        init_emf%uniform_b0_im(i,j) = init_b0_im(i,j)

        init_emf%mfunc_expr_e_re(i,j) = trim(init_e_re_mfunc(i,j))
        init_emf%mfunc_expr_e_im(i,j) = trim(init_e_im_mfunc(i,j))
        init_emf%mfunc_expr_b_re(i,j) = trim(init_b_re_mfunc(i,j))
        init_emf%mfunc_expr_b_im(i,j) = trim(init_b_im_mfunc(i,j))
      enddo
    enddo
  end select

  this%init_emf%init_move_window = init_move_window

  ! external fields
  select case(trim(ext_fld))
  case ( "none" )
    this%ext_fld = p_extfld_none
  case ( "static" )
    this%ext_fld = p_extfld_static
  case ( "dynamic" )
    this%ext_fld = p_extfld_dynamic
  case default
    print *, "Error reading emf parameters"
    print *, "Invalid value for the ext_fld parameter"
    print *, "Available external field types are 'none', 'static' and 'dynamic'"
    print *, "aborting..."
    stop
  end select

  do i=1, p_f_dim
    select case(trim(type_ext_b(i)))
    case( "none" )
      this%ext_emf%type_b(i) = p_emf_none
    case( "uniform" )
      this%ext_emf%type_b(i) = p_emf_uniform
    case( "math func" )
      this%ext_emf%type_b(i) = p_emf_math
    case( "python" )

      if ( .not. if_py_util() ) then
        if ( mpi_node() == 0 ) then
          write(0,*) "   Error: requested ""python"" for type_ext_b, but code was"
          write(0,*) "   not compiled with python support. Please see the"
          write(0,*) "   PY_FCOMPILEFLAGS and PY_FLINKFLAGS flags in your config file."
          write(0,*) "   aborting..."
        endif
        stop
      endif

      this%ext_emf%type_b(i) = p_emf_python
      this%ext_emf%py_mod = trim(ext_py_mod)
      this%ext_emf%py_func = trim(ext_py_func)
    case default
      print *, "Error reading emf parameters"
      print *, "Invalid value for the type_ext_b parameter"
      print *, "Available external B-field types in cylindrical coords are"
      print *, "'none', 'uniform', 'math func' and 'python'"
      print *, "aborting..."
      stop
    end select

    select case(trim(type_ext_e(i)))
    case( "none" )
      this%ext_emf%type_e(i) = p_emf_none
    case( "uniform" )
      this%ext_emf%type_e(i) = p_emf_uniform
    case( "math func" )
      this%ext_emf%type_e(i) = p_emf_math
    case( "python" )

      if ( .not. if_py_util() ) then
        if ( mpi_node() == 0 ) then
          write(0,*) "   Error: requested ""python"" for type_ext_e, but code was"
          write(0,*) "   not compiled with python support. Please see the"
          write(0,*) "   PY_FCOMPILEFLAGS and PY_FLINKFLAGS flags in your config file."
          write(0,*) "   aborting..."
        endif
        stop
      endif

      this%ext_emf%type_e(i) = p_emf_python
      this%ext_emf%py_mod = trim(ext_py_mod)
      this%ext_emf%py_func = trim(ext_py_func)
    case default
      print *, "Error reading emf parameters"
      print *, "Invalid value for the type_ext_e parameter"
      print *, "Available external E-field types in cylindrical coords are"
      print *, "'none', 'uniform', 'math func' and 'python'"
      print *, "aborting..."
      stop
    end select
  end do

  ! set external e and b fields modal decomposition
  select type( ext_emf => this%ext_emf )
  class is( t_emf_gridval_cyl_modes )

    ext_emf%n_cyl_modes = this%n_cyl_modes

    do j = 1, p_f_dim
      ext_emf%uniform_e0(j) = ext_e0_re(0,j)
      ext_emf%uniform_b0(j) = ext_b0_re(0,j)

      ext_emf%mfunc_expr_e(j) = trim(ext_e_re_mfunc(0,j))
      ext_emf%mfunc_expr_b(j) = trim(ext_b_re_mfunc(0,j))
    enddo

    ! store values for rest of modes in the special cyl modes arrays
    do j = 1, p_f_dim
      do i = 1, this%n_cyl_modes
        ext_emf%uniform_e0_re(i,j) = ext_e0_re(i,j)
        ext_emf%uniform_e0_im(i,j) = ext_e0_im(i,j)
        ext_emf%uniform_b0_re(i,j) = ext_b0_re(i,j)
        ext_emf%uniform_b0_im(i,j) = ext_b0_im(i,j)

        ext_emf%mfunc_expr_e_re(i,j) = trim(ext_e_re_mfunc(i,j))
        ext_emf%mfunc_expr_e_im(i,j) = trim(ext_e_im_mfunc(i,j))
        ext_emf%mfunc_expr_b_re(i,j) = trim(ext_b_re_mfunc(i,j))
        ext_emf%mfunc_expr_b_im(i,j) = trim(ext_b_im_mfunc(i,j))
      enddo
    enddo
  end select

  ! read boundary condtion information
  call this % bnd_con % read_input( input_file, periodic, if_move )

  ! Do not allow PML with local smoothing, for now...
  if (this%smooth_type == p_emfsmooth_local) then
    do i=1, p_x_dim
      if (.not. periodic(i)) then
        if ( (type(this%bnd_con,p_lower,i) == p_bc_vpml) &
              .or. (type(this%bnd_con,p_upper,i) == p_bc_vpml) ) then
          print *, "PML boundary conditions cannot be used with 'local' EMF smoothing."
          print *, "Aborting..."
          stop
        endif
      endif
    end do
  endif

  ! vpml currently not working for cylindrical modes; needs to be re-worked!
  if (.not. periodic(1)) then
    if ( (type(this%bnd_con,p_lower,1) == p_bc_vpml) &
          .or. (type(this%bnd_con,p_upper,1) == p_bc_vpml) ) then
      print *, "PML boundary conditions currently not working in x1 for cyl_modes."
      print *, "Consider using lindman instead."
      print *, "Aborting..."
      stop
    endif
  endif

  ! Read field solver parameters, if necessary
  call this % solver % read_input( input_file, dx, dt )

  ! read smoothing information
  call read_nml( this%smooth, input_file )

  ! read diagnostics information
  select type( diag => this%diag ); class is ( t_diag_emf_cyl_modes )
    diag % n_cyl_modes = this % n_cyl_modes
  end select
  select type( boosted_diag => this%diag%boosted_diag ); class is ( t_boosted_diag_cyl_modes )
    boosted_diag%n_cyl_modes = this%n_cyl_modes
  end select

  call this % diag % read_input( input_file, gamma )

  ! Must enumerate the solvers here to access the advance_cyl_modes function
  select type( solver => this%solver )
  class is( t_emf_solver_cyl_modes_yee )

    solver%n_cyl_modes = this%n_cyl_modes

  class is( t_emf_solver_cyl_modes_fei )

    solver%n_cyl_modes = this%n_cyl_modes

  end select

end subroutine read_input_emf_cyl_modes
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Initialize emf_cyl_modes object
!   - this routine also calls the superclass (t_emf) setup routine first
!-----------------------------------------------------------------------------------------
subroutine init_emf_cyl_modes( this, part_grid_center, part_interpolation, &
                               g_space, grid, gc_min, dx, tstep, tmin, tmax, &
                               no_co, send_msg, recv_msg, &
                               restart, restart_handle, sim_options )

  implicit none

  class( t_emf_cyl_modes ), intent( inout ), target  ::  this

  logical, intent(in) :: part_grid_center
  integer, intent(in) :: part_interpolation
  type( t_space ),     intent(in) :: g_space
  class( t_grid ), intent(in) :: grid
  integer, dimension(:,:), intent(in) :: gc_min
  real( p_double ),   intent(in) :: tmin, tmax
  type( t_time_step ), intent(in) :: tstep
  real( p_double ), dimension(:), intent(in) :: dx
  class( t_node_conf ), intent(in), target :: no_co
  type(t_vdf_msg), dimension(2), intent(inout) :: send_msg, recv_msg
  logical, intent(in) :: restart
  type( t_restart_handle ), intent(in) :: restart_handle
  type( t_options ), intent(in) :: sim_options

  ! setup superclass data
  call this % t_emf % init( part_grid_center, part_interpolation, &
                            g_space, grid, gc_min, dx, tstep, tmin, tmax, &
                            no_co, send_msg, recv_msg, &
                            restart, restart_handle, sim_options )

  if ( restart ) then
    ! The t_cyl_modes objects must be read in after init based on how restarts
    ! are handled with the t_emf object
    call this%restart_read_cm( restart_handle )
  else
    ! setup main fields
    call setup( this%e_cyl_m, this%n_cyl_modes, this%e )
    call setup( this%b_cyl_m, this%n_cyl_modes, this%b )


    ! set initial values for e and b
    select type( init_emf => this%init_emf ); class is( t_emf_gridval_cyl_modes )
      call init_emf%set_fld_values_cm( this%e_cyl_m, this%b_cyl_m, g_space, &
                                       grid%my_nx( p_lower, : ), 0.0d0 )
    end select
  endif

  ! setup external fields in same way regardless of whether or not its a restart
  ! no need to read/write that info, just get it from input deck
  if ( this%ext_fld /= p_extfld_none ) then
    ! setup external fields
    call setup( this%ext_e_cyl_m, this%n_cyl_modes, this%ext_e )
    call setup( this%ext_b_cyl_m, this%n_cyl_modes, this%ext_b )

    ! setup the external field object has already been called

    ! setup the values of static external fields
    ! values of dynamic external fields will be updated at every time step in update_particle_fld
    if ( this%ext_fld == p_extfld_static ) then
      select type( ext_emf => this%ext_emf ); class is( t_emf_gridval_cyl_modes )
        call ext_emf%set_fld_values_cm( this%ext_e_cyl_m, this%ext_b_cyl_m, g_space, &
                                        grid%my_nx( p_lower, : ), 0.0d0 )
      end select
    endif

  endif

  if ( this%part_fld_alloc ) then
    ! allocate and setup b_part/e_part
    allocate( this%b_part_cyl_m )
    allocate( this%e_part_cyl_m )

    call setup( this%b_part_cyl_m, this%n_cyl_modes, this%b_part )
    call setup( this%e_part_cyl_m, this%n_cyl_modes, this%e_part )
  else
    ! just point to main emf data
    this%b_part_cyl_m => this%b_cyl_m
    this%e_part_cyl_m => this%e_cyl_m
  endif

end subroutine init_emf_cyl_modes
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Cleanup emf_cyl_modes object
!-----------------------------------------------------------------------------------------
subroutine cleanup_emf_cyl_modes( this )

  implicit none

  class(t_emf_cyl_modes), intent(inout) :: this

  ! cleanup superclass
  call this % t_emf % cleanup()

  ! cleanup local data structures
  if ( this%part_fld_alloc ) then
    call cleanup( this%e_part_cyl_m )
    call cleanup( this%b_part_cyl_m )

    deallocate( this%e_part_cyl_m )
    deallocate( this%b_part_cyl_m )
  endif

  call cleanup( this%e_cyl_m )
  call cleanup( this%b_cyl_m )

  ! cleanup psi
  call cleanup_psi_cyl_modes( this%psi_cyl_m )

  call cleanup( this%ext_e_cyl_m )
  call cleanup( this%ext_b_cyl_m )

end subroutine cleanup_emf_cyl_modes
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Advance EMF fields and cyl_modes laser
!-----------------------------------------------------------------------------------------
subroutine advance_emf_cyl_modes( this, jay, dt, send_msg, recv_msg )

  implicit none

  class( t_emf_cyl_modes ), intent( inout )  ::  this

  class( t_current ), intent(inout) :: jay
  real(p_double), intent(in) :: dt

  type(t_vdf_msg), dimension(2), intent(inout) :: send_msg, recv_msg

  call begin_event(fsolverev)

  select type( bnd_con => this%bnd_con )
  class is( t_emf_bound_cyl_modes )

    select type( jay )
    class is( t_current_cyl_modes )

      ! Must enumerate the solvers here to access the advance_cyl_modes function
      select type( solver => this%solver )
      class is( t_emf_solver_cyl_modes_yee )

        call solver % advance_cyl_modes( this%e_cyl_m, this%b_cyl_m, jay%jay_cyl_m, &
                                          dt, bnd_con )

      class is( t_emf_solver_cyl_modes_fei )

        call solver % advance_cyl_modes( this%e_cyl_m, this%b_cyl_m, jay%jay_cyl_m, &
                                          dt, bnd_con )

      end select
    end select
  end select

  call end_event( fsolverev )

end subroutine advance_emf_cyl_modes
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine write_checkpoint_emf_cyl_modes( this, restart_handle )
!-----------------------------------------------------------------------------------------
!       write object information into a restart file
!-----------------------------------------------------------------------------------------

  implicit none

  class( t_emf_cyl_modes ), intent(in) :: this
  type( t_restart_handle ), intent(inout) :: restart_handle

  character(len=*), parameter :: err_msg = 'error writing restart data for emf_cyl_modes object.'
  integer :: ierr

  ! write superclass checkpoint data first
  call this % t_emf % write_checkpoint( restart_handle )

  ! write cyl_modes checkpoint id
  restart_io_wr( p_emf_cyl_modes_rst_id, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  ! write restart data for cylindrical modes
  call restart_write( this%b_cyl_m, restart_handle )
  call restart_write( this%e_cyl_m, restart_handle )

end subroutine write_checkpoint_emf_cyl_modes
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine restart_read_cm_emf_cyl_modes( this, restart_handle )
!-----------------------------------------------------------------------------------------
!       read object information from a restart file after init is complete
!-----------------------------------------------------------------------------------------

  implicit none

  class( t_emf_cyl_modes ), intent(inout) :: this
  type( t_restart_handle ), intent(in) :: restart_handle

  character(len=*), parameter :: err_msg = 'error reading restart data for emf_cyl_modes object.'
  character(len=len(p_emf_cyl_modes_rst_id)) :: rst_id
  integer :: ierr

  ! read checkpoint id
  restart_io_rd( rst_id, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  ! check if restart file is compatible
  if ( rst_id /= p_emf_cyl_modes_rst_id ) then
    ERROR('Corrupted restart file, or restart file ')
    ERROR('from incompatible binary (emf_cyl_modes)')
    ERROR('rst_id = "', rst_id, '"' )
    call abort_program(p_err_rstrd)
  endif

  ! read checkpoint data for cylindrical modes
  call restart_read( this%b_cyl_m, this%b, restart_handle )
  call restart_read( this%e_cyl_m, this%e, restart_handle )

end subroutine restart_read_cm_emf_cyl_modes
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine move_window_emf_cyl_modes( this, g_space, nx_p_min, need_fld_val )
!-----------------------------------------------------------------------------------------
!       move boundaries of the electro-magnetic field
!-----------------------------------------------------------------------------------------

  implicit none

  class( t_emf_cyl_modes ), intent( inout )  ::  this
  type( t_space ),     intent(in) :: g_space
  integer, dimension(:), intent(in) :: nx_p_min
  logical, intent(in) :: need_fld_val

  integer, dimension( p_max_dim ) :: lb

  ! first call superclass's move_window
  call this % t_emf % move_window(g_space, nx_p_min, need_fld_val)

  ! then do cylindrical modes part

  if (nx_move( g_space, 1 ) > 0) then

    ! move window for e and b field
    call move_window( this%e_cyl_m, g_space )
    call move_window( this%b_cyl_m, g_space )

    if ( need_fld_val .and. this%init_emf%init_move_window ) then
      lb(1) = this%b%nx_(1) - nx_move( g_space, 1 )
      lb(2:p_x_dim ) = 1 - this%b%gc_num_( p_lower, 2:p_x_dim )

      select type( init_emf => this%init_emf ); class is( t_emf_gridval_cyl_modes )
        call init_emf%set_fld_values_cm( this%e_cyl_m, this%b_cyl_m, g_space, &
                                         nx_p_min, 0.0d0, lbin = lb )
      end select
    endif

    ! move window for external fields: this is only required for static external fields,
    ! dynamic fields are recalculated at every iteration
    if ( this%ext_fld == p_extfld_static ) then
      ! shift local data
      call move_window( this%ext_b_cyl_m, g_space )
      call move_window( this%ext_e_cyl_m, g_space )

      lb(1) = this%b%nx_(1) + this%b%gc_num_( p_upper, 1 ) - nx_move( g_space, 1 )
      lb(2:p_x_dim ) = 1 - this%b%gc_num_( p_lower, 2:p_x_dim )

      select type( ext_emf => this%ext_emf ); class is( t_emf_gridval_cyl_modes )
        call ext_emf%set_fld_values_cm( this%ext_e_cyl_m, this%ext_b_cyl_m, g_space, &
                                        nx_p_min, 0.0d0, lbin = lb )
      end select
    endif

  endif

end subroutine move_window_emf_cyl_modes
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine update_boundary_emf_cyl_modes( this, g_space, no_co, send_msg, recv_msg )

  class( t_emf_cyl_modes ), intent( inout )  ::  this
  type( t_space ),     intent(in) :: g_space
  class( t_node_conf ), intent(in) :: no_co
  type(t_vdf_msg), dimension(2), intent(inout) :: send_msg, recv_msg

  integer :: i

#ifdef __CYL_MODES_DEBUG__
  call check_nan( this%b_cyl_m, "B field before update_boundary" )
  call check_nan( this%e_cyl_m, "E field before update_boundary" )
#endif

  ! update_boundary for the superclass
  call this % t_emf % update_boundary( g_space, no_co, send_msg, recv_msg )

  select type( bnd_con => this%bnd_con )
  class is( t_emf_bound_cyl_modes )

    ! update boundary for cylindrical modes
    do i = 1, this%n_cyl_modes
      call bnd_con%update_boundary_cm( this%b_cyl_m%pf_re(i), this%e_cyl_m%pf_re(i), &
                                        g_space, no_co , i, p_real, send_msg, recv_msg)
      call bnd_con%update_boundary_cm( this%b_cyl_m%pf_im(i), this%e_cyl_m%pf_im(i), &
                                        g_space, no_co, i, p_imag, send_msg, recv_msg)
    enddo

  end select

#ifdef __CYL_MODES_DEBUG__
  call check_nan( this%b_cyl_m, "B field after update_boundary" )
  call check_nan( this%e_cyl_m, "E field after update_boundary" )
#endif

end subroutine update_boundary_emf_cyl_modes
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
function solver_type_class_emf_cyl_modes( this )

    implicit none

    integer :: solver_type_class_emf_cyl_modes
    class( t_emf_cyl_modes ), intent(in) :: this

    select type( solver => this%solver )
    class is ( t_emf_solver_cyl_modes_yee )
        solver_type_class_emf_cyl_modes = p_emf_cyl_modes_yee
    class is ( t_emf_solver_cyl_modes_fei )
        solver_type_class_emf_cyl_modes = p_emf_cyl_modes_fei
    class default
        ! Unknown class
        solver_type_class_emf_cyl_modes = -1
    end select

end function solver_type_class_emf_cyl_modes
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
function solver_type_name_emf_cyl_modes( this, name )

    implicit none

    integer :: solver_type_name_emf_cyl_modes
    class( t_emf_cyl_modes ), intent(in) :: this
    character( len = * ), intent(in) :: name

    select case( trim(name) )
    case( "yee" )
        solver_type_name_emf_cyl_modes = p_emf_cyl_modes_yee
    case( "custom", "fei" )
        solver_type_name_emf_cyl_modes = p_emf_cyl_modes_fei
    case default
        ! Unknown class
        solver_type_name_emf_cyl_modes = -1
    end select

end function solver_type_name_emf_cyl_modes
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
! These allow the list of report quantities to be expanded / replaced by subclasses
!-----------------------------------------------------------------------------------------
function avail_report_quants_emf_cyl_modes( this )

  class( t_diag_emf_cyl_modes ), intent(in) :: this

  integer :: avail_report_quants_emf_cyl_modes

  avail_report_quants_emf_cyl_modes = size( p_report_quants )

end function avail_report_quants_emf_cyl_modes
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine init_report_quants_emf_cyl_modes( this )

  class( t_diag_emf_cyl_modes ), intent(inout) :: this

  if ( .not. associated( this % report_quants ) ) then
    allocate( this % report_quants( this % avail_report_quants( ) ) )
  endif

  ! These reports, and only these reports
  this % report_quants( 1 : size(p_report_quants) ) = p_report_quants


  ! Boosted diagnostic reports
  if ( .not. associated( this % boosted_report_quants ) ) then
    allocate( this % boosted_report_quants( size(p_boosted_report_quants) ) )
  endif

  ! A subclass may have a larger set of diagnostics
  this%boosted_report_quants(1:size(p_boosted_report_quants)) = p_boosted_report_quants

end subroutine init_report_quants_emf_cyl_modes
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Initialize individual diagnostics
!-----------------------------------------------------------------------------------------
subroutine init_diag_emf_cyl_modes( this, no_ext_fld, part_fld_alloc, interpolation )

  implicit none

  class( t_diag_emf_cyl_modes ), intent(inout) :: this
  logical, intent(in) :: no_ext_fld
  logical, intent(in) :: part_fld_alloc
  integer, intent(in) :: interpolation

  type(t_vdf_report), pointer :: report
  integer, parameter :: izero = ichar('0')
  real(p_double) :: interp_offset

  ! Field quantities only move in the x1 direction based on interpolation
  select case( interpolation )
  case( p_linear, p_cubic )
    interp_offset = 0.0_p_double
  case( p_quadratic, p_quartic )
    interp_offset = 0.5_p_double
  case default
    interp_offset = 0.0_p_double
    ERROR('Interpolation value not supported')
    call abort_program(p_err_invalid)
  end select

  ! setup cyl_modes specific diagnostics
  report => this%reports
  do
    if ( .not. associated( report ) ) exit

    ! Set default units for all diagnostics
    report%xname  = (/'x1', 'x2', 'x3'/)
    report%xlabel = (/'x_1', 'x_2', 'x_3'/)
    report%xunits = (/'c / \omega_p', 'c / \omega_p', 'c / \omega_p'/)

    ! these are just dummy values for now
    report%time_units = '1 / \omega_p'
    report%dt         = 1.0

    ! all field quantities are aligned in time
    report%offset_t = 0.0_p_double

    report%fileLabel = ''
    report%basePath  = trim(path_mass) // 'FLD' // p_dir_sep

    ! Set specific diagnostic parameters
    select case ( report%quant )
    case ( p_e1_cyl_m, p_e2_cyl_m, p_e3_cyl_m )
      report%label = 'E_'//char(izero + 1 + report%quant - p_e1_cyl_m)
      report%units = 'm_e c \omega_p e^{-1}'
      report%offset_x = (/ 0.0_p_double+interp_offset, 0.0_p_double, 0.0_p_double /)
      report%offset_x(report%quant-p_e1_cyl_m+1) = report%offset_x(report%quant-p_e1_cyl_m+1) + 0.5_p_double

    case ( p_b1_cyl_m, p_b2_cyl_m, p_b3_cyl_m )
      report%label = 'B_'//char(izero + 1 + report%quant - p_b1_cyl_m)
      report%units = 'm_e c \omega_p e^{-1}'
      report%offset_x = (/ 0.5_p_double+interp_offset, 0.5_p_double, 0.5_p_double /)
      report%offset_x(report%quant-p_b1_cyl_m+1) = report%offset_x(report%quant-p_b1_cyl_m+1) - 0.5_p_double

    case ( p_psi_cyl_m )
        report%label     = '\Psi'
        report%units     = 'm c^2 e^{-1}'
        ! aligned with e1
        report%offset_x = (/ 0.5_p_double+interp_offset, 0.0_p_double, 0.0_p_double /)

    case ( p_ext_e1_cyl_m, p_ext_e2_cyl_m, p_ext_e3_cyl_m )
      if ( no_ext_fld ) then
        SCR_ROOT('(*warning*) External field diagnostic requested, but external fields')
        SCR_ROOT('            not in use.')
        report%ndump = 0
      else
        report%label = 'E_'//char(izero + 1 + report%quant - p_ext_e1_cyl_m)//'^{ext}'
        report%units = 'm_e c \omega_p e^{-1}'
        report%offset_x = (/ 0.0_p_double+interp_offset, 0.0_p_double, 0.0_p_double /)
        report%offset_x(report%quant-p_ext_e1_cyl_m+1) = report%offset_x(report%quant-p_ext_e1_cyl_m+1) + 0.5_p_double
      endif

    case ( p_ext_b1_cyl_m, p_ext_b2_cyl_m, p_ext_b3_cyl_m )
      if ( no_ext_fld ) then
        SCR_ROOT('(*warning*) External field diagnostic requested, but external fields')
        SCR_ROOT('            not in use.')
        report%ndump = 0
      else
        report%label = 'B_'//char(izero + 1 + report%quant - p_ext_b1_cyl_m)//'^{ext}'
        report%units = 'm_e c \omega_p e^{-1}'
        report%offset_x = (/ 0.5_p_double+interp_offset, 0.5_p_double, 0.5_p_double /)
        report%offset_x(report%quant-p_ext_b1_cyl_m+1) = report%offset_x(report%quant-p_ext_b1_cyl_m+1) - 0.5_p_double
      endif

    case ( p_part_e1_cyl_m, p_part_e2_cyl_m, p_part_e3_cyl_m )
      if ( .not. part_fld_alloc ) then
        SCR_ROOT('(*warning*) Particle fields diagnostic requested, but external/smoothed')
        SCR_ROOT('            fields are not in use. Use main field diagnostics instead.')
        report%ndump = 0
      else
        report%label = 'E_'//char(izero + 1 + report%quant - p_part_e1_cyl_m)//'^{part}'
        report%units = 'm_e c \omega_p e^{-1}'
        report%offset_x = (/ 0.0_p_double+interp_offset, 0.0_p_double, 0.0_p_double /)
        report%offset_x(report%quant-p_part_e1_cyl_m+1) = report%offset_x(report%quant-p_part_e1_cyl_m+1) + 0.5_p_double
      endif

    case ( p_part_b1_cyl_m, p_part_b2_cyl_m, p_part_b3_cyl_m )
      if ( .not. part_fld_alloc ) then
        SCR_ROOT('(*warning*) Particle fields diagnostic requested, but external/smoothed')
        SCR_ROOT('            fields are not in use. Use main field diagnostics instead.')
        report%ndump = 0
      else
        report%label = 'B_'//char(izero + 1 + report%quant - p_part_b1_cyl_m)//'^{part}'
        report%units = 'm_e c \omega_p e^{-1}'
        report%offset_x = (/ 0.5_p_double+interp_offset, 0.5_p_double, 0.5_p_double /)
        report%offset_x(report%quant-p_part_b1_cyl_m+1) = report%offset_x(report%quant-p_part_b1_cyl_m+1) - 0.5_p_double
      endif

    case ( p_ene_e1_cyl_m, p_ene_e2_cyl_m, p_ene_e3_cyl_m )
      report%label = 'E_'//char(izero + 1 + report%quant - p_ene_e1_cyl_m)//'^2'
      report%units = 'm_e^2 c^2 \omega_p^2 e^{-2}'
      report%offset_x = (/ 0.0_p_double+interp_offset, 0.0_p_double, 0.0_p_double /)
      report%offset_x(report%quant-p_ene_e1_cyl_m+1) = report%offset_x(report%quant-p_ene_e1_cyl_m+1) + 0.5_p_double

    case ( p_ene_b1_cyl_m, p_ene_b2_cyl_m, p_ene_b3_cyl_m )
      report%label = 'B_'//char(izero + 1 + report%quant - p_ene_b1_cyl_m)//'^2'
      report%units = 'm_e^2 c^2 \omega_p^2 e^{-2}'
      report%offset_x = (/ 0.5_p_double+interp_offset, 0.5_p_double, 0.5_p_double /)
      report%offset_x(report%quant-p_ene_b1_cyl_m+1) = report%offset_x(report%quant-p_ene_b1_cyl_m+1) - 0.5_p_double

    case ( p_ene_e_cyl_m )
      report%label = 'E^2'
      report%units = 'm_e^2 c^2 \omega_p^2 e^{-2}'

    case ( p_ene_b_cyl_m )
      report%label = 'B^2'
      report%units = 'm_e^2 c^2 \omega_p^2 e^{-2}'

    case ( p_ene_emf_cyl_m )
      report%label = 'E^2 + B^2'
      report%units = 'm_e^2 c^2 \omega_p^2 e^{-2}'

    case ( p_div_e_cyl_m )
      report%label = '\bf{\nabla}\cdot\bf{E}'
      report%units = 'm_e c \omega_p e^{-1}'
      report%offset_x = (/ 0.0_p_double+interp_offset, 0.0_p_double, 0.0_p_double /)

    case ( p_div_b_cyl_m )
      report%label = '\bf{\nabla}\cdot\bf{B}'
      report%units = 'm_e c \omega_p e^{-1}'
      report%offset_x = (/ 0.5_p_double+interp_offset, 0.5_p_double, 0.5_p_double /)

    case ( p_s1_int_cyl_m, p_s2_int_cyl_m, p_s3_int_cyl_m )
      report%label = 'S_'//char(izero + 1 + report%quant - p_s1_int_cyl_m)
      report%units = 'm_e^2 c^2 \omega_p^2 e^{-2}'
      report%offset_x = (/ 0.0_p_double+interp_offset, 0.0_p_double, 0.0_p_double /)

    case default
      ! must be a superclass diagnostic, ignore
      ! print *, 'In setup_diag_cyl_modes, unknown report = ', report % name
      continue

    end select

    ! process next report
    report => report % next
  enddo

  if (diag_emf_ev==0) then
    diag_emf_ev = create_event('EMF diagnostics')
    if( this%if_use_boosted_diag ) then
      diag_emf_boost_extract_ev = create_event('boosted EMF diag extract portion')
      diag_emf_boost_ev = create_event('boosted EMF diag write portion')
    endif
  endif

end subroutine init_diag_emf_cyl_modes
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!       report on electro-magnetic field - diagnostic
!-----------------------------------------------------------------------------------------
subroutine report_diag_emf_cyl_modes( this, g_space, grid, no_co, tstep, t, send_msg, recv_msg )

  implicit none

  class( t_emf_cyl_modes ), intent(inout) :: this
  type( t_space ), intent(in) :: g_space
  class( t_grid ), intent(in) :: grid
  class( t_node_conf ), intent(in) :: no_co
  type( t_time_step ), intent(in) :: tstep
  real(p_double), intent(in) :: t
  type(t_vdf_msg), dimension(2), intent(inout) :: send_msg, recv_msg

  ! temporary objects for diagnostics
  type( t_cyl_modes ) :: cyl_m_a

  type( t_vdf_report ), pointer :: rep
  integer :: i

  type( t_cyl_modes ), pointer :: psi_ptr => null()
  call begin_event( diag_emf_ev )

  ! run cyl_modes specific diagnostics
  rep => this%diag%reports
  do
    if ( .not. associated( rep ) ) exit

    if ( if_report( rep, tstep ) ) then

      select case ( rep%quant )
      case ( p_e1_cyl_m, p_e2_cyl_m, p_e3_cyl_m )
        call report_cyl_modes( rep, this%e_cyl_m, &
          rep%quant - p_e1_cyl_m + 1, g_space, grid, no_co, tstep, t )
      case ( p_psi_cyl_m )
         ! calculate psi diagnostic
         call get_psi_cyl_modes( this, n(tstep), no_co, psi_ptr, send_msg, recv_msg )
         ! report it

         call report_cyl_modes( rep, psi_ptr, 1, g_space, grid, no_co, tstep, t )

         psi_ptr => null()
      case ( p_b1_cyl_m, p_b2_cyl_m, p_b3_cyl_m )
        call report_cyl_modes( rep, this%b_cyl_m, &
          rep%quant - p_b1_cyl_m + 1, g_space, grid, no_co, tstep, t )

      case ( p_ext_e1_cyl_m, p_ext_e2_cyl_m, p_ext_e3_cyl_m )
        call report_cyl_modes( rep, this%ext_e_cyl_m, &
          rep%quant - p_ext_e1_cyl_m + 1, g_space, grid, no_co, tstep, t )

      case ( p_ext_b1_cyl_m, p_ext_b2_cyl_m, p_ext_b3_cyl_m )
        call report_cyl_modes( rep, this%ext_b_cyl_m, &
          rep%quant - p_ext_b1_cyl_m + 1, g_space, grid, no_co, tstep, t )

      case ( p_part_e1_cyl_m, p_part_e2_cyl_m, p_part_e3_cyl_m )
        call report_cyl_modes( rep, this%e_part_cyl_m, &
          rep%quant - p_part_e1_cyl_m + 1, g_space, grid, no_co, tstep, t )

      case ( p_part_b1_cyl_m, p_part_b2_cyl_m, p_part_b3_cyl_m )
        call report_cyl_modes( rep, this%b_part_cyl_m, &
          rep%quant - p_part_b1_cyl_m + 1, g_space, grid, no_co, tstep, t )

      case ( p_ene_e1_cyl_m, p_ene_e2_cyl_m, p_ene_e3_cyl_m )
        ! calculate energy in field component
        call new_copy( cyl_m_a, this%e_cyl_m, f_dim=1, copy=.true., &
                       fc=rep%quant - p_ene_e1_cyl_m + 1 )
        call pow( cyl_m_a, 2 )

        ! report it
        call report_cyl_modes( rep, cyl_m_a, 1, g_space, grid, no_co, tstep, t )

        ! free temporary memory
        call cleanup( cyl_m_a, cleanup_0=.true. )

      case ( p_ene_b1_cyl_m, p_ene_b2_cyl_m, p_ene_b3_cyl_m )
        ! calculate energy in field component
        call new_copy( cyl_m_a, this%b_cyl_m, f_dim=1, copy=.true., &
                       fc=rep%quant - p_ene_b1_cyl_m + 1 )
        call pow( cyl_m_a, 2 )

        ! report it
        call report_cyl_modes( rep, cyl_m_a, 1, g_space, grid, no_co, tstep, t )

        ! free temporary memory
        call cleanup( cyl_m_a, cleanup_0=.true. )

      case ( p_ene_e_cyl_m )
        ! calculate energy in electric field
        call new_copy( cyl_m_a, this%e_cyl_m, f_dim=1, copy=.true., fc=1 )
        call pow( cyl_m_a, 2 )
        do i = 2, 3
          call add_pow( cyl_m_a, this%e_cyl_m, 2, fc=i )
        enddo

        ! report it
        call report_cyl_modes( rep, cyl_m_a, 1, g_space, grid, no_co, tstep, t )

        ! free temporary memory
        call cleanup( cyl_m_a, cleanup_0=.true. )

      case ( p_ene_b_cyl_m )
        ! calculate energy in magnetic field
        call new_copy( cyl_m_a, this%b_cyl_m, f_dim=1, copy=.true., fc=1 )
        call pow( cyl_m_a, 2 )
        do i = 2, 3
          call add_pow( cyl_m_a, this%b_cyl_m, 2, fc=i )
        enddo

        ! report it
        call report_cyl_modes( rep, cyl_m_a, 1, g_space, grid, no_co, tstep, t )

        ! free temporary memory
        call cleanup( cyl_m_a, cleanup_0=.true. )

      case ( p_ene_emf_cyl_m )
        ! calculate energy in electro-magnetic field
        call new_copy( cyl_m_a, this%e_cyl_m, f_dim=1, copy=.true., fc=1 )
        call pow( cyl_m_a, 2 )
        do i = 2, 3
          call add_pow( cyl_m_a, this%e_cyl_m, 2, fc=i )
        enddo
        do i = 1, 3
          call add_pow( cyl_m_a, this%b_cyl_m, 2, fc=i )
        enddo

        ! report it
        call report_cyl_modes( rep, cyl_m_a, 1, g_space, grid, no_co, tstep, t )

        ! free temporary memory
        call cleanup( cyl_m_a, cleanup_0=.true. )

      case ( p_div_e_cyl_m )
        ! calculate electric field divergence
        call new_copy( cyl_m_a, this%e_cyl_m, f_dim=1, zero=.true. )
        call div( this%e_cyl_m, cyl_m_a, this%gix_pos(2), .false. )

        ! report it
        call report_cyl_modes( rep, cyl_m_a, 1, g_space, grid, no_co, tstep, t )

        ! free temporary memory
        call cleanup( cyl_m_a, cleanup_0=.true. )

      case ( p_div_b_cyl_m )
        ! calculate magnetic field divergence
        call new_copy( cyl_m_a, this%b_cyl_m, f_dim=1, zero=.true. )
        call div( this%b_cyl_m, cyl_m_a, this%gix_pos(2), .true. )

        ! report it
        call report_cyl_modes( rep, cyl_m_a, 1, g_space, grid, no_co, tstep, t )

        ! free temporary memory
        call cleanup( cyl_m_a, cleanup_0=.true. )

      case ( p_s1_int_cyl_m, p_s2_int_cyl_m, p_s3_int_cyl_m )
        ! calculate component of Poynting flux
        call new_copy( cyl_m_a, this%e_cyl_m, f_dim=1, zero=.false. )
        call poynting_emf_cyl_modes( this, rep%quant - p_s1_int_cyl_m + 1, cyl_m_a )

        ! report it
        call report_cyl_modes( rep, cyl_m_a, 1, g_space, grid, no_co, tstep, t )

        ! free temporary memory
        call cleanup( cyl_m_a, cleanup_0=.true. )

      end select

    endif

    rep => rep%next
  enddo

  if( this%diag%if_use_boosted_diag ) then

    call begin_event( diag_emf_boost_extract_ev )

    select type( boosted_diag => this%diag%boosted_diag ); class is ( t_boosted_diag_cyl_modes )
      call boosted_diag%extract_cyl_modes( tstep, t, grid, g_space, this%e_cyl_m, this%b_cyl_m )
    end select

    call end_event( diag_emf_boost_extract_ev )

    call begin_event( diag_emf_boost_ev )

    call this%diag%boosted_diag%report( tstep, grid )

    call end_event( diag_emf_boost_ev )

  endif

  call end_event( diag_emf_ev )

end subroutine report_diag_emf_cyl_modes
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Poynting Flux (S) Calculation for cyl modes
! Here we are dumping poynting flux data averaged over the angular coordinate theta:
!     S = integral_{0}^{2*\pi} E x B d\theta
! However, we dump the contribution from each mode separately.
! It is calculated at the lower corner of the cell (in the same position as the charge)
!-----------------------------------------------------------------------------------------
subroutine poynting_emf_cyl_modes( this, comp, s )

  implicit none

  class( t_emf_cyl_modes ), intent(in) :: this
  integer, intent(in) :: comp
  type( t_cyl_modes ), intent(inout) :: s

  ! integer :: n_modes
  integer :: mode

  ! normalization factor associated with angular integration
  real(p_k_fld) :: norm

  ! add mode 0 contribution to poynting vector field
  norm = 1.0_p_k_fld
  call poynting_2d_cyl_modes(this%e_cyl_m%pf_re(0), this%b_cyl_m%pf_re(0), comp, s%pf_re(0), norm)

  ! add higher modes real and imag components contribution to poynting vector field
  do mode = 1, this%n_cyl_modes
    norm = 0.5_p_k_fld
    call poynting_2d_cyl_modes(this%e_cyl_m%pf_re(mode), this%b_cyl_m%pf_re(mode), comp, s%pf_re(mode), norm)
    call poynting_2d_cyl_modes(this%e_cyl_m%pf_im(mode), this%b_cyl_m%pf_im(mode), comp, s%pf_im(mode), norm)
  enddo

end subroutine poynting_emf_cyl_modes
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine poynting_2d_cyl_modes( e_vdf_data, b_vdf_data, comp, s, norm )

  implicit none

  type( t_vdf ), intent(in)    :: e_vdf_data
  type( t_vdf ), intent(in)    :: b_vdf_data
  integer, intent(in)          :: comp
  type( t_vdf ), intent(inout) :: s
  real(p_k_fld), intent(in)    :: norm ! normalization factor associated with angular integration

  integer :: i1, i2
  real( p_k_fld ) :: e1, e2, e3, b1, b2, b3

  select case ( comp )

  case (1)

    do i2 = 1, e_vdf_data%nx_(2)
      do i1 = 1, e_vdf_data%nx_(1)

        e2 = 0.5_p_k_fld * ( e_vdf_data%f2(2, i1, i2) + e_vdf_data%f2(2, i1, i2-1) )
        e3 = e_vdf_data%f2(3, i1, i2)

        b2 = 0.5_p_k_fld * ( b_vdf_data%f2( 2, i1, i2  ) + b_vdf_data%f2( 2, i1-1, i2) )
        b3 = 0.25_p_k_fld * ( b_vdf_data%f2( 3, i1, i2  ) + b_vdf_data%f2( 3, i1-1, i2) +&
                              b_vdf_data%f2( 3, i1, i2-1) + b_vdf_data%f2( 3, i1-1, i2-1))

        S%f2(1, i1, i2) = norm * (e2*b3 - b2*e3)   ! S1

      enddo
    enddo

  case (2)

    do i2 = 1, e_vdf_data%nx_(2)
      do i1 = 1, e_vdf_data%nx_(1)

        e1 = 0.5_p_k_fld * ( e_vdf_data%f2(1, i1, i2) + e_vdf_data%f2(1, i1-1, i2) )
        e3 = e_vdf_data%f2(3, i1, i2)

        b1 = 0.5_p_k_fld * ( b_vdf_data%f2( 1, i1, i2  ) + b_vdf_data%f2( 1, i1, i2-1)  )
        b3 = 0.25_p_k_fld * ( b_vdf_data%f2( 3, i1, i2  ) + b_vdf_data%f2( 3, i1-1, i2) +&
                              b_vdf_data%f2( 3, i1, i2-1) + b_vdf_data%f2( 3, i1-1, i2-1))

        S%f2(1, i1, i2) = norm * (e3*b1 - e1*b3)   ! S2

      enddo
    enddo

  case (3)

    do i2 = 1, e_vdf_data%nx_(2)
      do i1 = 1, e_vdf_data%nx_(1)

        e1 = 0.5_p_k_fld * ( e_vdf_data%f2(1, i1, i2) + e_vdf_data%f2(1, i1-1, i2) )
        e2 = 0.5_p_k_fld * ( e_vdf_data%f2(2, i1, i2) + e_vdf_data%f2(2, i1, i2-1) )

        b1 = 0.5_p_k_fld * ( b_vdf_data%f2( 1, i1, i2  ) + b_vdf_data%f2( 1, i1, i2-1)  )
        b2 = 0.5_p_k_fld * ( b_vdf_data%f2( 2, i1, i2  ) + b_vdf_data%f2( 2, i1-1, i2) )

        S%f2(1, i1, i2) = norm * (e1*b2 - b1*e2)   ! S3

      enddo
    enddo

  end select

end subroutine poynting_2d_cyl_modes
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Update fields to be used by particles
!-----------------------------------------------------------------------------------------
subroutine update_particle_fld_emf_cyl_modes( this, n, t, dt, space, nx_p_min )

  implicit none

  class( t_emf_cyl_modes ), intent( inout )  ::  this
  integer, intent(in) :: n
  real( p_double ), intent(in) :: t, dt
  type( t_space ), intent(in) :: space
  integer, dimension(:), intent(in) :: nx_p_min

  integer :: m

  ! Smoothing and external fields

  if ( this%ext_fld == p_extfld_dynamic ) then
    call this%ext_emf%set_fld_values( this%ext_e, this%ext_b, space, nx_p_min, t )
    select type( ext_emf => this%ext_emf ); class is( t_emf_gridval_cyl_modes )
      call ext_emf%set_fld_values_cm( this%ext_e_cyl_m, this%ext_b_cyl_m, space, &
                                      nx_p_min, t )
    end select
  endif

  ! do local smooth of fields if necessary
  if ( this%smooth_type == p_emfsmooth_local ) then
    ! Stop smoothing after given iteration
    if (this%smooth_nmax > 0 .and. n >= this%smooth_nmax) then
      this%smooth_type = p_emfsmooth_none
    else if ( mod( n, this%smooth_niter ) == 0) then

      ! smooth e and b fields
      call begin_event( fsmoothev )

      call smooth( this%e_cyl_m, this%smooth )
      call smooth( this%b_cyl_m, this%smooth )

      call end_event( fsmoothev )
    endif
  endif


  ! Update fields for particle interpolation
  if ( this%part_fld_alloc ) then

    ! Copy current emf field values to e_part and b_part
    ! mode 0 is contained in the main arrays
    this%e_part = this%e
    this%b_part = this%b

    ! copy rest of modes
    do m = 1, this%n_cyl_modes
      this % e_part_cyl_m % pf_re(m) = this % e_cyl_m % pf_re(m)
      this % e_part_cyl_m % pf_im(m) = this % e_cyl_m % pf_im(m)
      this % b_part_cyl_m % pf_re(m) = this % b_cyl_m % pf_re(m)
      this % b_part_cyl_m % pf_im(m) = this % b_cyl_m % pf_im(m)
    enddo

    ! smooth e and b fields using standard smoothing
    select case ( this%smooth_type )
    case( p_emfsmooth_stand )
      call begin_event( fsmoothev )

      call smooth( this%e_part_cyl_m, this%smooth )
      call smooth( this%b_part_cyl_m, this%smooth )

      call end_event( fsmoothev )

    case( p_emfsmooth_nci )
      call begin_event( fsmoothev )
      call filter_nci( this%e_part, this%b_part, dt, this%part_grid_center )
      do m = 1, this%n_cyl_modes
        call filter_nci( this%e_part_cyl_m%pf_re(m), this%b_part_cyl_m%pf_re(m), dt, &
                         this%part_grid_center )
        call filter_nci( this%e_part_cyl_m%pf_im(m), this%b_part_cyl_m%pf_im(m), dt, &
                         this%part_grid_center )
      enddo
      call end_event( fsmoothev )

    end select

    ! Add external fields (this could be optimized in the situations where the external fields
    ! are constant - avoid reading ext_e0, ext_b0)
    if ( this%ext_fld /= p_extfld_none ) then
      call add( this%e_part_cyl_m, this%ext_e_cyl_m )
      call add( this%b_part_cyl_m, this%ext_b_cyl_m )
    endif

    ! Grid center is always false

  endif

end subroutine update_particle_fld_emf_cyl_modes
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Update fields to be used by particles
!-----------------------------------------------------------------------------------------
subroutine reshape_emf_cyl_modes( this, old_lb, new_lb, no_co, send_msg, recv_msg )

  implicit none

  class( t_emf_cyl_modes ), intent( inout ), target  ::  this
  class( t_grid ), intent(in) :: old_lb, new_lb
  class( t_node_conf ), intent(in) :: no_co
  type(t_vdf_msg), dimension(2), intent(inout) :: send_msg, recv_msg

  call this % t_emf % reshape( old_lb, new_lb, no_co, send_msg, recv_msg )

  ! reshape the vdf objects
  call reshape_copy( this%b_cyl_m, old_lb, new_lb, this%b, no_co, send_msg, recv_msg )
  call reshape_copy( this%e_cyl_m, old_lb, new_lb, this%e, no_co, send_msg, recv_msg )

  ! if necessary reshape external fields
  if ( this%ext_fld /= p_extfld_none ) then

    ! Reshape the vdfs
    call reshape_nocopy( this%ext_b_cyl_m, new_lb, this%ext_b )
    call reshape_nocopy( this%ext_e_cyl_m, new_lb, this%ext_e )

    ! Recalculate the fields for static fields (dynamic fields will be recalculated at next
    ! iteration
    if ( this%ext_fld /= p_extfld_static ) then

      ! I need to pass the necessary parameters down here
      ERROR( 'Not implemented yet' )
      call abort_program( p_err_notimplemented )
    endif
  endif

  ! if necessary reshape particle interpolation fields
  ! these do not require copying since they are updated from e and b at each timestep
  if ( this%part_fld_alloc ) then
    call reshape_nocopy( this%b_part_cyl_m, new_lb, this%b_part )
    call reshape_nocopy( this%e_part_cyl_m, new_lb, this%e_part )
  endif

  ! reshape vdf objects if needed

  ! Reshape psi object
  if ( associated( this%psi_cyl_m%data%pf_re ) ) then
    call reshape_nocopy( this%psi_cyl_m%data, new_lb )
  endif

end subroutine reshape_emf_cyl_modes
!-----------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
! Initialize psi object
!---------------------------------------------------------------------------------------------------
subroutine init_psi_cyl_modes( emf, no_co )

  implicit none

  class( t_emf_cyl_modes ), intent(inout) :: emf
  class( t_node_conf ), intent(in) :: no_co

  integer :: nodes_x, ngp_x
  integer :: color, key, ierr

  call new_copy( emf%psi_cyl_m%data, emf%e_cyl_m, f_dim=1, zero=.true.)

  ! Initialize communicator if necessary
  nodes_x = nx( no_co, 1 )

  if ( nodes_x > 1 ) then
    ngp_x = my_ngp( no_co, 1 )

    ! Communicator for MPI_SCAN use

    ! Reverse the order of the individual ranks so that the scan function adds the values
    ! right to left
    key = nodes_x - ngp_x

    ! split communicators along x2 and x3 directions so that nodes in the same communicator
    ! share the same x2 and x3 coordinates
    color = my_ngp( no_co, 2 )

    call MPI_COMM_SPLIT( comm( no_co ), color, key, emf%psi_cyl_m%comm, ierr )

  else

    emf%psi_cyl_m%comm = mpi_comm_null

  endif

end subroutine init_psi_cyl_modes
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! Cleanup psi object
!---------------------------------------------------------------------------------------------------
subroutine cleanup_psi_cyl_modes( psi_cyl_m )

  implicit none

  type( t_emf_psi_cyl_modes ), intent(inout) :: psi_cyl_m

  integer :: ierr

  call cleanup(psi_cyl_m%data, cleanup_0 = .true.)
  psi_cyl_m%n = -1

  if ( psi_cyl_m%comm /= mpi_comm_null ) then
    call MPI_COMM_FREE( psi_cyl_m%comm, ierr )
    psi_cyl_m%comm = mpi_comm_null
  endif

end subroutine cleanup_psi_cyl_modes
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------
subroutine get_psi_cyl_modes( emf, n, no_co, psi_cyl_m, send_msg, recv_msg )

  implicit none

  class( t_emf_cyl_modes ), intent(inout), target :: emf
  integer, intent(in) :: n
  class( t_node_conf ), intent(in) :: no_co

  type( t_cyl_modes ), pointer :: psi_cyl_m
  type(t_vdf_msg), dimension(2), intent(inout) :: send_msg, recv_msg


  ! Initialize the Psi object if required
  if ( emf%psi_cyl_m%n < 0 ) then
    call init_psi_cyl_modes( emf, no_co )
  endif


  ! If data for this iteration is not available calculate it
  if ( n /= emf%psi_cyl_m%n ) then
    call update_psi_cyl_modes( emf, no_co, send_msg, recv_msg )
    emf%psi_cyl_m%n = n
  endif

  ! return pointer to psi data
  psi_cyl_m => emf%psi_cyl_m%data

end subroutine get_psi_cyl_modes
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! Calculate the Psi diagnostic in parallel, MPI_SCAN version
!---------------------------------------------------------------------------------------------------
subroutine update_psi_cyl_modes( emf, no_co, send_msg, recv_msg )

  implicit none

  class( t_emf_cyl_modes ), intent(inout) :: emf
  class( t_node_conf ), intent(in) :: no_co
  type(t_vdf_msg), dimension(2), intent(inout) :: send_msg, recv_msg


  integer :: i1, i2, nx1, nx2, mode, total_num_modes
  real( p_k_fld ) :: dx

  ! real/imaginary vdf pointers
  type(t_vdf), pointer :: psi_re, psi_im, e1_re, e1_im


  real( p_k_fld ), dimension(:,:), pointer   :: accQ3D, acc_outQ3D

  integer :: nodes_x, ngp_x
  integer :: mpi_type, count, ierr

  !---

  ! sychronize nodes (this is mostly for benchmarking purposes)
  call no_co % barrier( )

  call begin_event(getpsi_ev)

  ! calculate local psi
  dx = real( emf%e%dx_(1), p_k_fld )
  nx1 = emf%psi_cyl_m%data%pf_re(0)%nx_( 1 )
  nx2 = emf%psi_cyl_m%data%pf_re(0)%nx_( 2 )

  do mode = 0, emf%n_cyl_modes

    psi_re => emf%psi_cyl_m%data%pf_re(mode)
    e1_re => emf%e_cyl_m%pf_re(mode)

    do i2 = 1, psi_re%nx_( 2 )
      psi_re%f2( 1, nx1, i2 ) = e1_re%f2( 1, nx1, i2 ) * dx
      do i1 = nx1-1, 1, -1
        psi_re%f2( 1, i1, i2 ) = psi_re%f2( 1, i1+1, i2 ) + e1_re%f2( 1, i1, i2 ) * dx
      enddo
    enddo

    if(mode > 0) then
      psi_im => emf%psi_cyl_m%data%pf_im(mode)
      e1_im => emf%e_cyl_m%pf_im(mode)

      do i2 = 1, psi_re%nx_( 2 )
        psi_im%f2( 1, nx1, i2 ) = e1_im%f2( 1, nx1, i2 ) * dx
        do i1 = nx1-1, 1, -1
          psi_im%f2( 1, i1, i2 ) = psi_im%f2( 1, i1+1, i2 ) + e1_im%f2( 1, i1, i2 ) * dx
        enddo
      enddo
    endif



  enddo



  ! Add contribution from other nodes to the right
  if ( no_num( no_co ) > 1 ) then

    nodes_x = nx( no_co, 1 )

    if ( nodes_x > 1 ) then

       ngp_x = my_ngp( no_co, 1 )

       select case ( p_k_fld )
         case( p_single )
           mpi_type = MPI_REAL
         case( p_double )
           mpi_type = MPI_DOUBLE_PRECISION
       end select

       total_num_modes = 1 + 2 * emf%n_cyl_modes
       ! create temp buffer with leftmost values
       call alloc( accQ3D, (/ nx2, total_num_modes /) )
       call alloc( acc_outQ3D, (/ nx2, total_num_modes /) )

       ! copy mode 0 - real component
       do i2 = 1, nx2
          accQ3D( i2, 1 ) = emf%psi_cyl_m%data%pf_re(0)%f2( 1, 1, i2 )
        enddo

      !copy rest of modes - real + im
       do mode = 1, emf%n_cyl_modes
          psi_re => emf%psi_cyl_m%data%pf_re(mode)
          psi_im => emf%psi_cyl_m%data%pf_im(mode)
          do i2 = 1, nx2
          accQ3D( i2, 2 * mode ) = psi_re%f2( 1, 1, i2 )
          accQ3D( i2, 2 * mode + 1 ) = psi_im%f2( 1, 1, i2 )
          enddo
       enddo

       ! add contributions from all nodes in communicator
       ! (the openmpi-1.2.9 implementation of MPI_EXSCAN + MPI_IN_PLACE is broken )
       count = nx2 * total_num_modes
       call MPI_EXSCAN( accQ3D, acc_outQ3D, count, mpi_type, MPI_SUM, emf%psi_cyl_m%comm, ierr )

       ! add contribution to local data
       if ( ngp_x < nodes_x ) then
           ! real mode m = 0 first
           psi_re => emf%psi_cyl_m%data%pf_re(0)
           do i2 = 1, nx2
             do i1 = 1, nx1
               psi_re%f2( 1, i1, i2 ) = psi_re%f2( 1, i1, i2 ) + acc_outQ3D( i2, 1 )
             enddo
           enddo

           ! do other modes
           do mode = 1, emf%n_cyl_modes
           psi_re => emf%psi_cyl_m%data%pf_re(mode)
           psi_im => emf%psi_cyl_m%data%pf_im(mode)
            do i2 = 1, nx2
               do i1 = 1, nx1
                 psi_re%f2( 1, i1, i2 ) = psi_re%f2( 1, i1, i2 ) + acc_outQ3D( i2, 2 * mode )
                 psi_im%f2( 1, i1, i2 ) = psi_im%f2( 1, i1, i2 ) + acc_outQ3D( i2, 2 * mode+1 )
            enddo
           enddo
           enddo

       endif

       ! delete temp buffer with leftmost values
       call freemem( accQ3D )
       call freemem( acc_outQ3D )


    endif
    ! fill guard cells with values from another node
    ! (this should be made optional)
    call update_boundary( emf%psi_cyl_m%data%pf_re(0), p_vdf_replace, no_co, send_msg, recv_msg )
    do mode = 1, emf%n_cyl_modes
    call update_boundary( emf%psi_cyl_m%data%pf_re(mode), p_vdf_replace, no_co, send_msg, recv_msg )
    call update_boundary( emf%psi_cyl_m%data%pf_im(mode), p_vdf_replace, no_co, send_msg, recv_msg )
    enddo
  endif

  call end_event(getpsi_ev)


end subroutine update_psi_cyl_modes
!---------------------------------------------------------------------------------------------------





!-----------------------------------------------------------------------------------------
!       report electro-magnetic field energy
!-----------------------------------------------------------------------------------------
subroutine report_energy_emf_cyl_modes( this, no_co, tstep, t, tmin )
!-----------------------------------------------------------------------------------------

  implicit none

  class( t_emf_cyl_modes ), intent(inout) :: this
  class( t_node_conf ), intent(in) :: no_co
  type( t_time_step ), intent(in) :: tstep
  real(p_double), intent(in) :: t, tmin

  real(p_double), dimension(4*p_f_dim*this%n_cyl_modes+2*p_f_dim) :: temp_int
  real(p_double), dimension(2*p_f_dim*this%n_cyl_modes+  p_f_dim) :: temp_b, temp_e
  real(p_double), dimension(p_f_dim) :: offset
  integer :: i, j, lb, ierr
  character(len=8), dimension(4*p_f_dim*this%n_cyl_modes+2*p_f_dim) :: fld_lbls

  character(len=256) :: full_name, path, formatter

  call begin_event( diag_emf_ev )

  ! reports on integrated field energy
  if (test_if_report( tstep, this%diag%ndump_fac_ene_int ) ) then

    call this%fill_data( tstep )
    temp_int = 0.0_p_double

    ! get local e and b field integrals
    if ( p_f_dim /= 3 ) then
      ERROR('report_energy_emf_cyl_modes not supported for ',p_f_dim,' field dimensions')
      call abort_program( p_err_notimplemented )
    endif

    offset = 0.0_p_double
    offset(1) = 0.5_p_double
    offset(3) = 0.5_p_double
    call total( this%b_cyl_m, temp_b, this%gix_pos(2), pow = 2, off = offset )

    offset = 0.0_p_double
    offset(2) = 0.5_p_double
    call total( this%e_cyl_m, temp_e, this%gix_pos(2), pow = 2, off = offset )

    ! interleave temp_b and temp_e arrays into temp_int
    do i = 1, 2*this%n_cyl_modes + 1
      lb = p_f_dim * (i-1)
      temp_int(2*lb+1:2*lb+3) = temp_b(lb+1:lb+3)
      temp_int(2*lb+4:2*lb+6) = temp_e(lb+1:lb+3)
    enddo

    ! fill in fld_lbls for mode 0
    do j = 1, p_f_dim
      write( fld_lbls(j        ), '(A,I0,A)' ) 'B', j, '-0-RE'
      write( fld_lbls(j+p_f_dim), '(A,I0,A)' ) 'E', j, '-0-RE'
      fld_lbls(j        ) = adjustr(fld_lbls(j        ))
      fld_lbls(j+p_f_dim) = adjustr(fld_lbls(j+p_f_dim))
    enddo

    ! fill in fld_lbls for other modes
    do i = 1, this%n_cyl_modes
      lb = 2*p_f_dim + 4*(i-1)*p_f_dim
      do j = 1, p_f_dim
        write( fld_lbls(lb+j          ), '(A,I0,A,I0,A)' ) 'B', j, '-', i, '-RE'
        write( fld_lbls(lb+j+  p_f_dim), '(A,I0,A,I0,A)' ) 'E', j, '-', i, '-RE'
        write( fld_lbls(lb+j+2*p_f_dim), '(A,I0,A,I0,A)' ) 'B', j, '-', i, '-IM'
        write( fld_lbls(lb+j+3*p_f_dim), '(A,I0,A,I0,A)' ) 'E', j, '-', i, '-IM'
        fld_lbls(lb+j          ) = adjustr(fld_lbls(lb+j          ))
        fld_lbls(lb+j+  p_f_dim) = adjustr(fld_lbls(lb+j+  p_f_dim))
        fld_lbls(lb+j+2*p_f_dim) = adjustr(fld_lbls(lb+j+2*p_f_dim))
        fld_lbls(lb+j+3*p_f_dim) = adjustr(fld_lbls(lb+j+3*p_f_dim))
      enddo
    enddo

    ! sum up results from all nodes
    call reduce_array( no_co, temp_int, operation = p_sum )

    ! save the data
    if ( root(no_co) ) then
      ! normalize energies
      ! 1/2*dv in cyl coordinates equates to 1/2*2*pi*dv = pi*dv
      temp_int = temp_int*this%e%dvol()*pi

      ! setup path and file names
      path  =  trim(path_hist)
      full_name = trim(path) // 'fld_ene'

      ! Open file and position at the last record
      if ( t == tmin ) then

        call mkdir( path, ierr )

        open (unit=file_id_fldene, file=full_name, status = 'REPLACE' , &
              form='formatted')

        ! Write Header
        write(file_id_fldene, '(A)' ) '! EM field energy per field component'
        write(formatter, '(A,I0,A)' ) '( A6, 1X,A15,', size(temp_int), '(1X,A23) )'
        write(file_id_fldene,formatter) 'Iter','Time    ',fld_lbls
      else

        open (unit=file_id_fldene, file=full_name, position = 'append', &
              form='formatted')
      endif

      write(formatter, '(A,I0,A)' ) '( I6, 1X,g15.8, ', size(temp_int), '(1X,es23.16) )'
      write(file_id_fldene,formatter) n(tstep), t, temp_int
      close(file_id_fldene)

    endif

  endif

  call end_event( diag_emf_ev )

end subroutine report_energy_emf_cyl_modes
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! set field values for the cyl_modes emf_gridval object
!-----------------------------------------------------------------------------------------
subroutine set_fld_values_cm_emf_gridval_cyl_modes( this, ext_e, ext_b, g_space, &
                                                    nx_p_min, t, lbin, ubin )

  implicit none

  class(t_emf_gridval_cyl_modes), intent(inout) :: this
  type( t_cyl_modes ), intent(inout) :: ext_e, ext_b
  type( t_space ), intent(in) :: g_space
  integer, dimension(:), intent(in) :: nx_p_min
  real( p_double ), intent(in) :: t
  integer, dimension(:), intent(in), optional :: lbin
  integer, dimension(:), intent(in), optional :: ubin

  real(p_double), dimension(p_max_dim) ::  g_xmin, ldx
  real(p_k_fparse), dimension(p_max_dim+1) :: x_eval

  real(p_double), dimension(p_max_dim,p_f_dim) :: odx_b
  real(p_double), dimension(p_max_dim,p_f_dim) :: odx_e
  integer :: i, j1, j2, m, mm, x_dim, gshift_i2

  integer, dimension(p_max_dim) :: lb, ub
  integer, dimension(p_max_dim) :: shift
  real(p_k_fld) :: f_t, o_t, n_e, o_e, phase_factor

  real(p_double), dimension(2,p_max_dim) :: x_bnd
  character(2) :: fld, mode, cmplx

  g_xmin(1:p_x_dim) = xmin(g_space)
  ldx(1:p_x_dim) = ext_b%pf_re(0)%dx()
  shift(1:p_x_dim) = nx_p_min(1:p_x_dim) - 1

  gshift_i2 = nx_p_min(p_r_dim) - 2

  ! only p_x_dim = 2 allowed w/ cyl modes
  x_dim = 2

  ! get boundaries to inject
  if (present(lbin)) then
    lb(1:2) = lbin(1:2)
  else
    lb(1) = lbound( ext_b%pf_re(0)%f2, 2)
    lb(2) = lbound( ext_b%pf_re(0)%f2, 3)
  endif
  if (present(ubin)) then
    ub(1:2) = ubin(1:2)
  else
    ub(1) = ubound( ext_b%pf_re(0)%f2, 2)
    ub(2) = ubound( ext_b%pf_re(0)%f2, 3)
  endif

  ! set the grid offsets

  odx_b(1,1) = 0.0_p_double*ldx(1)
  odx_b(2,1) = 0.5_p_double*ldx(2)

  odx_b(1,2) = 0.5_p_double*ldx(1)
  odx_b(2,2) = 0.0_p_double*ldx(2)

  odx_b(1,3) = 0.5_p_double*ldx(1)
  odx_b(2,3) = 0.5_p_double*ldx(2)

  odx_e(1,1) = 0.5_p_double*ldx(1)
  odx_e(2,1) = 0.0_p_double*ldx(2)

  odx_e(1,2) = 0.0_p_double*ldx(1)
  odx_e(2,2) = 0.5_p_double*ldx(2)

  odx_e(1,3) = 0.0_p_double*ldx(1)
  odx_e(2,3) = 0.0_p_double*ldx(2)

  select case( this%interpolation )
  case( p_linear, p_cubic )
    ! Nothing to change here
  case( p_quadratic, p_quartic )
    ! Only add shift in z, not r
    do j2 = 1, p_f_dim
      odx_b(1,j2) = odx_b(1,j2) + 0.5_p_double*ldx(1)
      odx_e(1,j2) = odx_e(1,j2) + 0.5_p_double*ldx(1)
    enddo
  end select

  do m = 1, this%n_cyl_modes
    do i = 1, p_f_dim

      ! set the external b field
      select case (this%type_b(i))
      case(p_emf_uniform)
        ! set the external field to the supplied uniform value
        do j2 = lb(2), ub(2)
          do j1 = lb(1), ub(1)
            ext_b % pf_re(m) % f2(i,j1,j2) = this % uniform_b0_re(m,i)
            ext_b % pf_im(m) % f2(i,j1,j2) = this % uniform_b0_im(m,i)
          enddo
        enddo

      case(p_emf_math)

        ! set the external field to the function value
        ! also set the value in the guard cells
        x_eval(3) = t
        do j2 = lb(2), ub(2)
          do j1 = lb(1), ub(1)

            ! set the proper correction to account for field positions
            x_eval(1) = g_xmin(1) + (j1-1+shift(1))*ldx(1)     + odx_b(1, i)
            x_eval(2) = g_xmin(2) + (j2-1+shift(2))*ldx(x_dim) + odx_b(2, i)

            ext_b%pf_re(m)%f2(i,j1,j2) = real( eval( this%mfunc_b_re(m,i), x_eval ), p_k_fld )
            ext_b%pf_im(m)%f2(i,j1,j2) = real( eval( this%mfunc_b_im(m,i), x_eval ), p_k_fld )

          enddo
        enddo

      case(p_emf_python)

        x_bnd(1,1) = g_xmin(1) + (lb(1)-1+shift(1))*ldx(1) + odx_b(1,i)
        x_bnd(1,2) = g_xmin(2) + (lb(2)-1+shift(2))*ldx(2) + odx_b(2,i)

        x_bnd(2,1) = g_xmin(1) + (ub(1)-1+shift(1))*ldx(1) + odx_b(1,i)
        x_bnd(2,2) = g_xmin(2) + (ub(2)-1+shift(2))*ldx(2) + odx_b(2,i)

        fld = ""
        write(fld,'(A,I1)') "b", i
        mode = ""
        write(mode,'(I2)') m

        cmplx = "re"
        call set_state( "data", ext_b%pf_re(m)%f2(i,lb(1):ub(1),lb(2):ub(2)) )
        call set_state( "fld", trim(fld) )
        call set_state( "mode", trim(mode) )
        call set_state( "cmplx", trim(cmplx) )
        call set_state( "x_bnd", x_bnd(:,1:p_x_dim) )
        call set_state( "g_xmin", g_xmin(1:p_x_dim) )
        call set_state( "g_xmax", xmax(g_space) )
        call set_state( "nx_p_min", nx_p_min(1:p_x_dim) )
        call set_state( "t", (/t/) )

        call call_function( trim(this%py_mod), trim(this%py_func) )

        call get_state( "data", ext_b%pf_re(m)%f2(i,lb(1):ub(1),lb(2):ub(2)) )

        cmplx = "im"
        call set_state( "data", ext_b%pf_im(m)%f2(i,lb(1):ub(1),lb(2):ub(2)) )
        call set_state( "cmplx", trim(cmplx) )

        call call_function( trim(this%py_mod), trim(this%py_func) )

        call get_state( "data", ext_b%pf_im(m)%f2(i,lb(1):ub(1),lb(2):ub(2)) )

      end select

      ! set the external e field
      select case (this%type_e(i))
      case(p_emf_uniform)
        ! set the external field to the supplied uniform value
        do j2 = lb(2), ub(2)
          do j1 = lb(1), ub(1)
            ext_e % pf_re(m) % f2(i,j1,j2) = this % uniform_e0_re(m,i)
            ext_e % pf_im(m) % f2(i,j1,j2) = this % uniform_e0_im(m,i)
          enddo
        enddo

      case(p_emf_math)

        ! set the external field to the function value
        ! also set the value in the guard cells
        x_eval(3) = t
        do j2 = lb(2), ub(2)
          do j1 = lb(1), ub(1)

            ! set the proper correction to account for field positions
            x_eval(1) = g_xmin(1) + (j1-1+shift(1))*ldx(1) + odx_e(1, i)
            x_eval(2) = g_xmin(2) + (j2-1+shift(2))*ldx(2) + odx_e(2, i)

            ext_e%pf_re(m)%f2(i,j1,j2) = real( eval( this%mfunc_e_re(m,i), x_eval ), p_k_fld )
            ext_e%pf_im(m)%f2(i,j1,j2) = real( eval( this%mfunc_e_im(m,i), x_eval ), p_k_fld )
          enddo
        enddo

      case(p_emf_python)

        x_bnd(1,1) = g_xmin(1) + (lb(1)-1+shift(1))*ldx(1) + odx_e(1,i)
        x_bnd(1,2) = g_xmin(2) + (lb(2)-1+shift(2))*ldx(2) + odx_e(2,i)

        x_bnd(2,1) = g_xmin(1) + (ub(1)-1+shift(1))*ldx(1) + odx_e(1,i)
        x_bnd(2,2) = g_xmin(2) + (ub(2)-1+shift(2))*ldx(2) + odx_e(2,i)

        fld = ""
        write(fld,'(A,I1)') "e", i
        mode = ""
        write(mode,'(I2)') m

        cmplx = "re"
        call set_state( "data", ext_e%pf_re(m)%f2(i,lb(1):ub(1),lb(2):ub(2)) )
        call set_state( "fld", trim(fld) )
        call set_state( "mode", trim(mode) )
        call set_state( "cmplx", trim(cmplx) )
        call set_state( "x_bnd", x_bnd(:,1:p_x_dim) )
        call set_state( "g_xmin", g_xmin(1:p_x_dim) )
        call set_state( "g_xmax", xmax(g_space) )
        call set_state( "nx_p_min", nx_p_min(1:p_x_dim) )
        call set_state( "t", (/t/) )

        call call_function( trim(this%py_mod), trim(this%py_func) )

        call get_state( "data", ext_e%pf_re(m)%f2(i,lb(1):ub(1),lb(2):ub(2)) )

        cmplx = "im"
        call set_state( "data", ext_e%pf_im(m)%f2(i,lb(1):ub(1),lb(2):ub(2)) )
        call set_state( "cmplx", trim(cmplx) )

        call call_function( trim(this%py_mod), trim(this%py_func) )

        call get_state( "data", ext_e%pf_im(m)%f2(i,lb(1):ub(1),lb(2):ub(2)) )

      end select

    enddo ! i = 1, p_f_dim
  enddo ! m = 1, this%n_cyl_modes

  ! fix near-axis field values
  ! note, the guard cells that reside below the axis do not hold true field values,
  ! but hold the field values appropriate for particle interpolation.
  ! i.e., we must have the following for E and B (not crucial):
  !     E_z(-dr/2)   =  (-1)^m E_z(dr/2)
  !     E_r(-dr)     = -(-1)^m E_r(dr)
  !     E_phi(-dr/2) = -(-1)^m E_phi(dr/2)
  ! if quartic interpolation is ever implemented, also handle cell -1 in r
  if ( gshift_i2 < 0 ) then

    f_t = real(4.0_p_k_fld/3.0_p_k_fld, p_k_fld)
    o_t = real(1.0_p_k_fld/3.0_p_k_fld, p_k_fld)
    n_e = real(9.0_p_k_fld/8.0_p_k_fld, p_k_fld)
    o_e = real(1.0_p_k_fld/8.0_p_k_fld, p_k_fld)

    do m = 0, this%n_cyl_modes

      if (m == 0) then
        do j1 = lb(1), ub(1)
          ! E1 can be nonzero on the axis
          ext_e%pf_re(m)%f2( 1, j1, 1 ) = ext_e%pf_re(m)%f2( 1, j1, 2 )
          ext_e%pf_re(m)%f2( 1, j1, 0 ) = ext_e%pf_re(m)%f2( 1, j1, 3 )

          ! E2 must be zero on the axis
          ext_e%pf_re(m)%f2( 2, j1, 1 ) = 0.0_p_k_fld
          ext_e%pf_re(m)%f2( 2, j1, 0 ) = - ext_e%pf_re(m)%f2( 2, j1, 2 )

          ! E3 must be zero on the axis
          ext_e%pf_re(m)%f2( 3, j1, 1 ) = - ext_e%pf_re(m)%f2( 3, j1, 2 )
          ext_e%pf_re(m)%f2( 3, j1, 0 ) = - ext_e%pf_re(m)%f2( 3, j1, 3 )

          ! B1 can be nonzero on the axis
          ! Don't modify ext_b%pf_re(m)%f2( 1, j1, 1 )
          ext_b%pf_re(m)%f2( 1, j1, 0 ) = ext_b%pf_re(m)%f2( 1, j1, 2 )

          ! B2 must be zero on the axis
          ext_b%pf_re(m)%f2( 2, j1, 1 ) = - ext_b%pf_re(m)%f2( 2, j1, 2 )
          ext_b%pf_re(m)%f2( 2, j1, 0 ) = - ext_b%pf_re(m)%f2( 2, j1, 3 )

          ! B3 must be zero on the axis
          ext_b%pf_re(m)%f2( 3, j1, 1 ) = 0.0_p_k_fld
          ext_b%pf_re(m)%f2( 3, j1, 0 ) = - ext_b%pf_re(m)%f2( 3, j1, 2 )
        enddo
      elseif (m == 1) then
        do j1 = lb(1), ub(1)
          ! E1 must be zero on the axis
          ext_e%pf_re(m)%f2( 1, j1, 1 ) = - ext_e%pf_re(m)%f2( 1, j1, 2 )
          ext_e%pf_re(m)%f2( 1, j1, 0 ) = - ext_e%pf_re(m)%f2( 1, j1, 3 )
          ext_e%pf_im(m)%f2( 1, j1, 1 ) = - ext_e%pf_im(m)%f2( 1, j1, 2 )
          ext_e%pf_im(m)%f2( 1, j1, 0 ) = - ext_e%pf_im(m)%f2( 1, j1, 3 )

          ! E2 can be nonzero on the axis
          ext_e%pf_re(m)%f2( 2, j1, 1 ) = f_t * ext_e%pf_re(m)%f2( 2, j1, 2 ) - o_t * ext_e%pf_re(m)%f2( 2, j1, 3 )
          ext_e%pf_re(m)%f2( 2, j1, 0 ) = ext_e%pf_re(m)%f2( 2, j1, 2 )
          ext_e%pf_im(m)%f2( 2, j1, 1 ) = f_t * ext_e%pf_im(m)%f2( 2, j1, 2 ) - o_t * ext_e%pf_im(m)%f2( 2, j1, 3 )
          ext_e%pf_im(m)%f2( 2, j1, 0 ) = ext_e%pf_im(m)%f2( 2, j1, 2 )

          ! E3 can be nonzero on the axis
          ext_e%pf_re(m)%f2( 3, j1, 1 ) = 2.0_p_k_fld * ext_e%pf_im(m)%f2( 2, j1, 1 ) - ext_e%pf_re(m)%f2( 3, j1, 2 )
          ext_e%pf_re(m)%f2( 3, j1, 0 ) = 2.0_p_k_fld * ext_e%pf_im(m)%f2( 2, j1, 1 ) - ext_e%pf_re(m)%f2( 3, j1, 3 )
          ext_e%pf_im(m)%f2( 3, j1, 1 ) =-2.0_p_k_fld * ext_e%pf_re(m)%f2( 2, j1, 1 ) - ext_e%pf_im(m)%f2( 3, j1, 2 )
          ext_e%pf_im(m)%f2( 3, j1, 0 ) =-2.0_p_k_fld * ext_e%pf_re(m)%f2( 2, j1, 1 ) - ext_e%pf_im(m)%f2( 3, j1, 3 )

          ! B1 must be zero on the axis
          ext_b%pf_re(m)%f2( 1, j1, 1 ) = 0.0_p_k_fld
          ext_b%pf_re(m)%f2( 1, j1, 0 ) = - ext_b%pf_re(m)%f2( 1, j1, 2 )
          ext_b%pf_im(m)%f2( 1, j1, 1 ) = 0.0_p_k_fld
          ext_b%pf_im(m)%f2( 1, j1, 0 ) = - ext_b%pf_im(m)%f2( 1, j1, 2 )

          ! B2 can be nonzero on the axis
          ext_b%pf_re(m)%f2( 2, j1, 1 ) = ext_b%pf_re(m)%f2( 2, j1, 2 )
          ext_b%pf_re(m)%f2( 2, j1, 0 ) = ext_b%pf_re(m)%f2( 2, j1, 3 )
          ext_b%pf_im(m)%f2( 2, j1, 1 ) = ext_b%pf_im(m)%f2( 2, j1, 2 )
          ext_b%pf_im(m)%f2( 2, j1, 0 ) = ext_b%pf_im(m)%f2( 2, j1, 3 )

          ! B3 can be nonzero on the axis
          ext_b%pf_re(m)%f2( 3, j1, 1 ) = n_e * ext_b%pf_im(m)%f2( 2, j1, 2 ) - o_e * ext_b%pf_im(m)%f2( 2, j1, 3 )
          ext_b%pf_re(m)%f2( 3, j1, 0 ) = ext_b%pf_re(m)%f2( 3, j1, 2 )
          ext_b%pf_im(m)%f2( 3, j1, 1 ) =-n_e * ext_b%pf_re(m)%f2( 2, j1, 2 ) + o_e * ext_b%pf_re(m)%f2( 2, j1, 3 )
          ext_b%pf_im(m)%f2( 3, j1, 0 ) = ext_b%pf_im(m)%f2( 3, j1, 2 )
        enddo
      else ! m >= 2

        if ( mod(m,2) == 0 ) then
          phase_factor = 1.0_p_k_fld
        else
          phase_factor = -1.0_p_k_fld
        endif

        mm = m - 1

        do j1 = lb(1), ub(1)
          do j2 = 1, m - 2
            ! Zero out first few derivatives for E1
            ext_e%pf_re(m)%f2( 1, j1, j2+1 ) = real((2*j2-1)**mm, p_k_fld) / real((2*mm-1)**mm, p_k_fld) * ext_e%pf_re(m)%f2( 1, j1, mm+1 )
            ext_e%pf_im(m)%f2( 1, j1, j2+1 ) = real((2*j2-1)**mm, p_k_fld) / real((2*mm-1)**mm, p_k_fld) * ext_e%pf_im(m)%f2( 1, j1, mm+1 )
          enddo

          ! E1 must be zero on the axis
          ext_e%pf_re(m)%f2( 1, j1, 1 ) = phase_factor * ext_e%pf_re(m)%f2( 1, j1, 2 )
          ext_e%pf_re(m)%f2( 1, j1, 0 ) = phase_factor * ext_e%pf_re(m)%f2( 1, j1, 3 )
          ext_e%pf_im(m)%f2( 1, j1, 1 ) = phase_factor * ext_e%pf_im(m)%f2( 1, j1, 2 )
          ext_e%pf_im(m)%f2( 1, j1, 0 ) = phase_factor * ext_e%pf_im(m)%f2( 1, j1, 3 )

          ! E2 must be zero on the axis
          ext_e%pf_re(m)%f2( 2, j1, 1 ) = 0.0_p_k_fld
          ext_e%pf_re(m)%f2( 2, j1, 0 ) = -phase_factor * ext_e%pf_re(m)%f2( 2, j1, 2 )
          ext_e%pf_im(m)%f2( 2, j1, 1 ) = 0.0_p_k_fld
          ext_e%pf_im(m)%f2( 2, j1, 0 ) = -phase_factor * ext_e%pf_im(m)%f2( 2, j1, 2 )

          ! E3 must be zero on the axis
          ext_e%pf_re(m)%f2( 3, j1, 1 ) = -phase_factor * ext_e%pf_re(m)%f2( 3, j1, 2 )
          ext_e%pf_re(m)%f2( 3, j1, 0 ) = -phase_factor * ext_e%pf_re(m)%f2( 3, j1, 3 )
          ext_e%pf_im(m)%f2( 3, j1, 1 ) = -phase_factor * ext_e%pf_im(m)%f2( 3, j1, 2 )
          ext_e%pf_im(m)%f2( 3, j1, 0 ) = -phase_factor * ext_e%pf_im(m)%f2( 3, j1, 3 )

          do j2 = 1, m - 2
            ! Zero out first few derivatives for B1
            ext_b%pf_re(m)%f2( 1, j1, j2+1 ) = real(j2**mm, p_k_fld) / real(mm**mm, p_k_fld) * ext_b%pf_re(m)%f2( 1, j1, mm+1 )
            ext_b%pf_im(m)%f2( 1, j1, j2+1 ) = real(j2**mm, p_k_fld) / real(mm**mm, p_k_fld) * ext_b%pf_im(m)%f2( 1, j1, mm+1 )
          enddo

          ! B1 must be zero on the axis
          ext_b%pf_re(m)%f2( 1, j1, 1 ) = 0.0_p_k_fld
          ext_b%pf_re(m)%f2( 1, j1, 0 ) = phase_factor * ext_b%pf_re(m)%f2( 1, j1, 2 )
          ext_b%pf_im(m)%f2( 1, j1, 1 ) = 0.0_p_k_fld
          ext_b%pf_im(m)%f2( 1, j1, 0 ) = phase_factor * ext_b%pf_im(m)%f2( 1, j1, 2 )

          ! B2 must be zero on the axis
          ext_b%pf_re(m)%f2( 2, j1, 1 ) = -phase_factor * ext_b%pf_re(m)%f2( 2, j1, 2 )
          ext_b%pf_re(m)%f2( 2, j1, 0 ) = -phase_factor * ext_b%pf_re(m)%f2( 2, j1, 3 )
          ext_b%pf_im(m)%f2( 2, j1, 1 ) = -phase_factor * ext_b%pf_im(m)%f2( 2, j1, 2 )
          ext_b%pf_im(m)%f2( 2, j1, 0 ) = -phase_factor * ext_b%pf_im(m)%f2( 2, j1, 3 )

          ! B3 must be zero on the axis
          ext_b%pf_re(m)%f2( 3, j1, 1 ) = 0.0_p_k_fld
          ext_b%pf_re(m)%f2( 3, j1, 0 ) = -phase_factor * ext_b%pf_re(m)%f2( 3, j1, 2 )
          ext_b%pf_im(m)%f2( 3, j1, 1 ) = 0.0_p_k_fld
          ext_b%pf_im(m)%f2( 3, j1, 0 ) = -phase_factor * ext_b%pf_im(m)%f2( 3, j1, 2 )
        enddo
      endif

    enddo ! m = 0, this%n_cyl_modes

  endif ! ( gshift_i2 < 0 )

end subroutine set_fld_values_cm_emf_gridval_cyl_modes
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine write_checkpoint_emf_gridval_cyl_modes( this, restart_handle )

  implicit none

  class( t_emf_gridval_cyl_modes ), intent(in) :: this
  type( t_restart_handle ), intent(inout) :: restart_handle

  character(len=*), parameter :: err_msg = 'error writing restart data for emf_gridval_cyl_modes object.'
  integer :: ierr

  call this % t_emf_gridval % write_checkpoint( restart_handle )

  restart_io_wr( p_emf_gridval_cyl_modes_rst_id, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  restart_io_wr( this%uniform_e0_re, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  restart_io_wr( this%uniform_e0_im, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  restart_io_wr( this%uniform_b0_re, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  restart_io_wr( this%uniform_b0_im, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  restart_io_wr( this%mfunc_expr_b_re, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  restart_io_wr( this%mfunc_expr_b_im, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  restart_io_wr( this%mfunc_expr_e_re, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  restart_io_wr( this%mfunc_expr_e_im, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

end subroutine write_checkpoint_emf_gridval_cyl_modes
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine restart_read_emf_gridval_cyl_modes( this, restart_handle )

  implicit none

  class( t_emf_gridval_cyl_modes ), intent(inout) :: this
  type( t_restart_handle ), intent(in) :: restart_handle

  character(len=*), parameter :: err_msg = 'error reading restart data for emf_gridval_cyl_modes object.'
  character(len=len(p_emf_gridval_cyl_modes_rst_id)) :: rst_id
  integer :: ierr

  call this % t_emf_gridval % restart_read( restart_handle )

  restart_io_rd( rst_id, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  ! check if restart file is compatible
  if ( rst_id /= p_emf_gridval_cyl_modes_rst_id) then
    ERROR('Corrupted restart file, or restart file')
    ERROR('from incompatible binary (emf_gridval_cyl_modes)')
    ERROR('rst_id = "', rst_id,'"')
    call abort_program(p_err_rstrd)
  endif

  restart_io_rd( this%uniform_e0_re, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd )

  restart_io_rd( this%uniform_e0_im, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd )

  restart_io_rd( this%uniform_b0_re, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd )

  restart_io_rd( this%uniform_b0_im, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd )

  restart_io_rd( this%mfunc_expr_b_re, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd )

  restart_io_rd( this%mfunc_expr_b_im, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd )

  restart_io_rd( this%mfunc_expr_e_re, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd )

  restart_io_rd( this%mfunc_expr_e_im, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd )


end subroutine restart_read_emf_gridval_cyl_modes
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine setup_emf_gridval_cyl_modes( this, dynamic, interpolation, restart, restart_handle )

  implicit none

  class( t_emf_gridval_cyl_modes ), intent(inout) :: this
  logical, intent(in) :: dynamic
  integer, intent(in) :: interpolation
  logical, intent(in) :: restart
  type( t_restart_handle ), intent(in) :: restart_handle

  integer :: i, j, nvars, ierr
  character(len=2), dimension(p_x_dim+1) :: vars

  ! if restarting read information from checkpoint file
  if ( restart ) then
    call this%restart_read( restart_handle )
  endif

  this%dynamic = dynamic
  this%interpolation = interpolation

  ! compile math functions if required
  do i=1, p_x_dim
    vars(i)='x'//char(ichar('0')+i)
  enddo

  if ( this%dynamic ) then
    nvars = p_x_dim + 1
    vars(nvars) = 't'
  else
    nvars = p_x_dim
  endif

  do i= 1, p_f_dim
    if (this%type_b(i) == p_emf_math) then
      call setup(this%mfunc_b(i), trim(this%mfunc_expr_b(i)), vars(1:nvars), ierr)
      call disp_error( trim(this%mfunc_expr_b(i)), "B", ierr )

      do j = 1, this%n_cyl_modes
        call setup(this%mfunc_b_re(j,i), trim(this%mfunc_expr_b_re(j,i)), vars(1:nvars), ierr)
        call disp_error( trim(this%mfunc_expr_b_re(j,i)), "B", ierr )

        call setup(this%mfunc_b_im(j,i), trim(this%mfunc_expr_b_im(j,i)), vars(1:nvars), ierr)
        call disp_error( trim(this%mfunc_expr_b_im(j,i)), "B", ierr )
      enddo
    endif

    if (this%type_e(i) == p_emf_math) then
      call setup(this%mfunc_e(i), trim(this%mfunc_expr_e(i)), vars(1:nvars), ierr)
      call disp_error( trim(this%mfunc_expr_e(i)), "E", ierr )

      do j = 1, this%n_cyl_modes
        call setup(this%mfunc_e_re(j,i), trim(this%mfunc_expr_e_re(j,i)), vars(1:nvars), ierr)
        call disp_error( trim(this%mfunc_expr_e_re(j,i)), "E", ierr )

        call setup(this%mfunc_e_im(j,i), trim(this%mfunc_expr_e_im(j,i)), vars(1:nvars), ierr)
        call disp_error( trim(this%mfunc_expr_e_im(j,i)), "E", ierr )
      enddo
    endif
  enddo

  contains

  subroutine disp_error( mfunc, comp, ierr )

    implicit none

    character(len=*), intent(in) :: mfunc
    character(len=1), intent(in) :: comp
    integer, intent(in) :: ierr

    if ( ierr /= 0) then
      print *, "(*error*) Error compiling supplied function :"
      print *, trim(mfunc)
      print *, "(*error*) for external/initial ",comp," field values."
      print *, "(*error*) bailing out..."
      call abort_program()
    endif

  end subroutine disp_error

end subroutine setup_emf_gridval_cyl_modes
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine cleanup_emf_gridval_cyl_modes( this )

  implicit none

  class( t_emf_gridval_cyl_modes ), intent(inout) :: this
  integer :: i, j

  call this % t_emf_gridval % cleanup()

  do i= 1, p_f_dim
    if (this%type_b(i) == p_emf_math) then
      do j = 1, this%n_cyl_modes
        call cleanup(this%mfunc_b_re(j,i))
        call cleanup(this%mfunc_b_im(j,i))
      enddo
    endif

    if (this%type_e(i) == p_emf_math) then
      do j = 1, this%n_cyl_modes
        call cleanup(this%mfunc_e_re(j,i))
        call cleanup(this%mfunc_e_im(j,i))
      enddo
    endif
  enddo


end subroutine cleanup_emf_gridval_cyl_modes
!-----------------------------------------------------------------------------------------

end module m_emf_cyl_modes
