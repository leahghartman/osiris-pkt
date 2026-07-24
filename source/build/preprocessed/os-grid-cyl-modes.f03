# 1 "cyl_modes/os-grid-cyl-modes.f03"
# 1 "<built-in>" 1
# 1 "<built-in>" 3
# 467 "<built-in>" 3
# 1 "<command line>" 1
# 1 "<built-in>" 2
# 1 "cyl_modes/os-grid-cyl-modes.f03" 2
! os-grid-cyl-modes.f03

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
# 4 "cyl_modes/os-grid-cyl-modes.f03" 2
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
# 5 "cyl_modes/os-grid-cyl-modes.f03" 2

module m_grid_cyl_modes

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
# 9 "cyl_modes/os-grid-cyl-modes.f03" 2

use m_system
use m_parameters
use m_grid_define
use m_input_file
use m_fparser, only : setup, p_max_expr_len
use m_cyl_modes, only : p_max_cyl_modes


implicit none

private

!-------------------------------------------------------------------------------
! t_grid_cyl_modes class definition
!-------------------------------------------------------------------------------
type, extends( t_grid ) :: t_grid_cyl_modes

  integer :: n_cyl_modes = -1 ! number of modes in quasi-3D geometry

  contains

  procedure :: read_input => read_input_grid_cyl_modes

end type t_grid_cyl_modes


public :: t_grid_cyl_modes

contains

!---------------------------------------------------
subroutine read_input_grid_cyl_modes( this, input_file, partition )
!---------------------------------------------------
! read necessary information from input
!---------------------------------------------------

  implicit none

  class( t_grid_cyl_modes ), intent(inout) :: this
  class( t_input_file ), intent(inout) :: input_file
  integer, intent(in), dimension(:) :: partition

  integer, dimension(p_x_dim) :: nx_p
  integer, dimension(p_x_dim) :: io_nmerge
  character(15) :: coordinates
  integer :: n_cyl_modes

  logical, dimension(p_x_dim) :: load_balance

  logical :: any_lb

  character(len=16) :: lb_type, lb_gather
  integer :: n_dynamic, start_load_balance
  logical :: balance_on_start
  real( p_single ) :: max_imbalance, cell_weight

  integer :: ndump_global_load, ndump_node_load, ndump_grid_load
  !define a string
  character (len=p_max_expr_len) :: spatial_loaddensity

  namelist /nl_grid/ nx_p, coordinates, n_cyl_modes, io_nmerge, &
                             load_balance, lb_type, lb_gather, n_dynamic, start_load_balance, &
                             balance_on_start, max_imbalance, cell_weight, &
                             ndump_global_load, ndump_node_load, ndump_grid_load, spatial_loaddensity


  integer :: i, ierr

  ! this will be replaced when dimension can be set at run time
  this%x_dim = p_x_dim


  nx_p = 0
  coordinates = "cylindrical" ! default is cylindrical in quasi-3D
  io_nmerge = 1
  n_cyl_modes = 0

  load_balance = .false.
  lb_type = "none"
  lb_gather = "sum"
  n_dynamic = -1
  spatial_loaddensity = "NO FUNCTION DEFINED"

  ! Default is to start load balance at timestep 0 (initialization)
  start_load_balance = -1
  balance_on_start = .false.

  max_imbalance = 0.0
  cell_weight = 0.0

  ndump_global_load = 0
  ndump_node_load = 0
  ndump_grid_load = 0

  ! Get namelist text from input file
  call get_namelist( input_file, "nl_grid", ierr )

  if (ierr /= 0) then
    if (ierr < 0) then
      print *, "Error reading grid parameters"
    else
      print *, "Error grid parameters missing"
    endif
    print *, "aborting..."
    stop
  endif

  read (input_file%nml_text, nml = nl_grid, iostat = ierr)
  if (ierr /= 0) then
    print *, "Error reading grid parameters"
    print *, "aborting..."
    stop
  endif

  ! validate global grid parameters
  do i = 1, this%x_dim
    if (nx_p(i) <= 0) then
      print *, "Invalid grid parameters, nx_p must be >= 0 ",&
                "for all dimensions."
      print *, "nx_p(",p_x_dim,") = ", nx_p
      print *, "(also check the number of dimensions of the binary)"
      print *, "aborting..."
      stop
    endif
  enddo

  ! store global grid size
  this%g_nx(1:this%x_dim) = nx_p(1:this%x_dim)

  select case ( trim( coordinates ) )
  case ( "cylindrical" )
    if (disp_out(input_file)) then
      if (mpi_node()==0) print *,'cylindrical coordinates - B1 on axis'
    endif
    this%coordinates = p_cylindrical_b
    if ( p_x_dim /= 2 ) then
      print *, "Invalid coordinates, cylindrical coordinates "
      print *, "are only available in 2D."
      stop
    endif
    this%n_cyl_modes = n_cyl_modes
    if ( n_cyl_modes < 0 .or. n_cyl_modes > p_max_cyl_modes ) then
      print *, "Invalid number of modes, n_cyl_modes must be between:"
      print '(A,I0,A)', " 0 and ", p_max_cyl_modes, " inclusive."
      stop
    endif
  case default
    print *, "Invalid coordinates, coordinates must be:"
    print *, "cylindrical' (quasi-3D)"
    stop
  end select

  this%start_load_balance = start_load_balance
  this%balance_on_start = balance_on_start

  select case (trim(lb_type))
  case("none")
    this%lb_type = p_no_load_balance
    this%load_balance = .false.

  case("test")
    this%lb_type = p_test_load_balance
    this%load_balance = .true.
    this%n_dynamic = n_dynamic

  case("static")
    this%lb_type = p_static_load_balance
    this%cell_weight = cell_weight

    ! The start load balance parameter is ignored, since this occurs a timestep 0
    this%start_load_balance = -1

  case("dynamic")
    this%lb_type = p_dynamic_load_balance
    if (n_dynamic < 0) then
      print *, ""
      print *, "Error in grid parameters"
      print *, "n_dynamic must be >= 1 when using dynamic load balancing"
      print *, "aborting..."
      stop
    endif
    this%n_dynamic = n_dynamic

    this%max_imbalance = max_imbalance
    this%cell_weight = cell_weight

    select case (trim( lb_gather ))
    case ("sum")
      this%lb_gather_max = .false.
    case ("max")
      this%lb_gather_max = .true.
    case default
      print *, ""
      print *, "Error in grid parameters"
      print *, "lb_gather must be 'sum' or 'max' when using dynamic load balancing"
      print *, "aborting..."
      stop
    end select
  case("expression")
    this%lb_type = p_expr_load_balance
    this%cell_weight = cell_weight

    select case (p_x_dim)
    case (1)
      call setup(this%spatial_loaddensity, trim(spatial_loaddensity), (/'x1'/), ierr)
    case (2)
      call setup(this%spatial_loaddensity, trim(spatial_loaddensity), (/'x1','x2'/), ierr)
    case (3)
      call setup(this%spatial_loaddensity, trim(spatial_loaddensity), (/'x1','x2','x3'/), ierr)
    end select

    ! check if function compiled ok
    if (ierr /= 0) then
      print *, ""
      print *, "Error compiling spatial_loaddensity"
      print *, "aborting..."
      stop
    endif

  case default
    print *, ""
    print *, "Error in grid parameters"
    print *, "Invalid lb_type selected."
    print *, "aborting..."
    stop
  end select

  if ( this%lb_type /= p_no_load_balance ) then

    ! Store load balance directions
    any_lb = .false.
    do i = 1, this%x_dim
      this%load_balance(i) = load_balance(i)
      if ( load_balance(i) ) any_lb = .true.
    enddo

    if ( .not. any_lb ) then
      write(0,*) ""
      write(0,*) "Error in grid parameters"
      write(0,*) "Load balance type '", trim(lb_type), "' selected, but no load balance"
      write(0,*) "direction specified"
      write(0,*) "aborting..."
      stop
    endif

    ! store cell weight
    this%cell_weight = cell_weight
    if ( cell_weight < 0.0 ) then
      write(0,*) ""
      write(0,*) "Error in grid parameters"
      write(0,*) "cell_weight must be >= 0."
      write(0,*) "aborting..."
      stop
    endif


  endif


  ! grid reports
  this%ndump_global_load = ndump_global_load
  this%ndump_node_load = ndump_node_load
  this%ndump_grid_load = ndump_grid_load

  ! Data merging for file output
  do i = 1, this%x_dim

    if ( io_nmerge(i) < 1 ) then
      write(0,*) "Error in grid parameters"
      write(0,*) "io_nmerge(:) must be >= 1 in all directions."
      write(0,*) "aborting..."
      stop
    endif

    if ( mod(partition(i), io_nmerge(i) ) /= 0 ) then
      write(0,*) "Error in grid parameters"
      write(0,*) "io_nmerge(:) must divide number of nodes in parallel partition exactly."
      write(0,*) "aborting..."
      stop
    endif

    this % io % n_merge(i) = io_nmerge(i)
  enddo

end subroutine read_input_grid_cyl_modes
!---------------------------------------------------

end module m_grid_cyl_modes
