# 1 "pgc/os-tridiag.f03"
# 1 "<built-in>" 1
# 1 "<built-in>" 3
# 467 "<built-in>" 3
# 1 "<command line>" 1
# 1 "<built-in>" 2
# 1 "pgc/os-tridiag.f03" 2
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
# 2 "pgc/os-tridiag.f03" 2
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
# 3 "pgc/os-tridiag.f03" 2

! use division optimization `PGC_DIVISION_OPTIMIZATION` to optimize
! complex division in terms of precision and roundoff errors
!
!#define PGC_DIVISION_OPTIMIZATION

! use LAPACK libraries for tridiagonal solver if you can finde compiler flag
!
!#define PGC_WITH_LAPACK

module m_tridiag

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
# 16 "pgc/os-tridiag.f03" 2

use m_system
use m_parameters

implicit none

private

! parameters for solver type
integer, parameter :: p_thomas_alg = 1
integer, parameter :: p_lu_alg = 2
integer, parameter :: p_lapack_alg = 3

! parameters for tridiagonal bands
integer, parameter :: p_tridiag_lower = -1
integer, parameter :: p_tridiag_middle = 0
integer, parameter :: p_tridiag_upper = 1

! parameters for grid informations
integer, parameter :: p_id_min = 1
integer, parameter :: p_id_max = 2
integer, parameter :: p_grid_size = 3

type t_line_nconf
  ! node conf for line communications

  ! mpi communicator
  integer :: comm = mpi_comm_null

  ! number of nodes
  integer :: no_num

  ! id of process
  integer :: my_aid

  ! flag for periodic node grid
  logical :: ifpr
end type t_line_nconf

type t_msg
  ! message type for packing and sending

  ! recipient/sender node
  integer :: node = -1

  ! indecies for each dimension
  integer :: l_dim1, u_dim1
  integer :: l_dim2, u_dim2

  ! message size
  integer :: msg_size

  ! message type
  integer :: msg_type

  ! message buffer
  complex( p_k_fld ), allocatable, dimension(:,:) :: buffer

  ! message tag (optional)
  integer :: tag = 0

  ! message request id
  integer :: handle = MPI_REQUEST_NULL
end type t_msg

type, abstract :: t_base_tridiag_solver
  ! base type for tridiagonl matrix solver
  ! function to calculate the coefficients
  procedure(coefficient_function), pointer :: coeff_func
  ! storage for coffiecients (-1: lower, 0: middle, 1: upper)
  complex( p_k_fld ), pointer, dimension(:,:) :: coeff
  integer :: lb, ub
  ! factor for periodic systems
  complex(p_k_fld) :: gam
  !
  contains
  procedure(base_tridiag_setup), deferred :: setup
  procedure(base_tridiag_solve), deferred :: solve
  procedure(base_tridiag_decompose), deferred :: decompose
  procedure(base_tridiag_cleanup), deferred :: cleanup
end type
abstract interface
  subroutine base_tridiag_setup( this, lb, ub )
    import t_base_tridiag_solver
    implicit none
    class(t_base_tridiag_solver), intent(inout) :: this
    integer, intent(in) :: lb
    integer, intent(in) :: ub
  end subroutine
end interface
abstract interface
  subroutine base_tridiag_solve( this, rhs )
    import t_base_tridiag_solver, p_k_fld
    implicit none
    class(t_base_tridiag_solver), intent(in) :: this
    complex(p_k_fld), intent(inout), dimension(this%lb:this%ub) :: rhs
  end subroutine
end interface
abstract interface
  subroutine base_tridiag_decompose( this )
    import t_base_tridiag_solver
    implicit none
    class(t_base_tridiag_solver), intent(inout) :: this
  end subroutine
end interface
abstract interface
  subroutine base_tridiag_cleanup( this )
    import t_base_tridiag_solver
    implicit none
    class(t_base_tridiag_solver), intent(inout) :: this
  end subroutine
end interface

! Tridiagonal solver based on the Thomas algorithm
type, extends(t_base_tridiag_solver) :: t_thomas_solver
  ! tridiagonal matrix solver based on Thomas algorithm
  ! gamma stores intermediate factors for forward and
  ! backward elimination - reduces amount of operations
  complex( p_k_fld ), pointer, dimension(:) :: gamma
  contains
  procedure :: setup => thomas_setup
  procedure :: solve => thomas_solve
  procedure :: decompose => thomas_decompose
  procedure :: cleanup => thomas_cleanup
end type

! Tridiagonal solver based on the LU decomposition
type, extends(t_base_tridiag_solver) :: t_lu_solver
  ! tridiagonal matrix solver based on LU decomposition
  complex( p_k_fld ), pointer, dimension(:) :: u
  complex( p_k_fld ), pointer, dimension(:) :: ipiv
  contains
  procedure :: setup => lu_setup
  procedure :: solve => lu_solve
  procedure :: decompose => lu_decompose
  procedure :: cleanup => lu_cleanup
end type
# 166 "pgc/os-tridiag.f03"
type t_tridiag
  ! tridiagonal class which handles finding solution and communication

  ! solver for solving tridiagonal system
  class(t_base_tridiag_solver), pointer :: solver

  ! storage for right hand side
  complex( p_k_fld ), pointer, dimension(:,:) :: rhs => null()

  ! for periodic boundaries - based on Sherman-Morrison Formula
  ! u_vec - non-tridiag elements / independent solution
  complex( p_k_fld ), pointer, dimension(:) :: u_vec

  ! trigger for iterative improvement of a solution to linear equations
  logical :: error_cleaning

  ! right-hand-side smoothing
  real( p_k_fld ), pointer, dimension(:) :: smoothing_kernel
  complex( p_k_fld ), pointer, dimension(:) :: raw_buffer
  integer :: smoothing_order
  logical :: pre_smoothing
  logical :: post_smoothing

  ! information about algebraic system - in transversal direction
  integer :: local_size = -1
  integer :: global_size = -1
  integer :: id_min, id_max
  integer :: ext_min, ext_max
  integer :: lb, ub

  ! in longitudinal direction
  integer :: stripes
  integer :: length

  ! information of node configuration
  type( t_line_nconf ) :: nconf

  ! message buffer for transpose communications
  type( t_msg ), private, allocatable, dimension(:) :: src_msg
  type( t_msg ), private, allocatable, dimension(:) :: rhs_msg

  ! guard cell information - 1: left cells, 2: right, 3: total number
  integer, dimension(1:3) :: gc_num

  contains

  procedure :: setup => setup_tridiag
  procedure :: solve => solve_tridiag
  procedure :: cleanup => cleanup_tridiag
  procedure :: setup_smoothing => setup_smoothing_tridiag
  procedure :: smooth => smooth_tridiag
  procedure :: periodic_preCalculation => periodic_preCalculation_tridiag
end type t_tridiag

abstract interface
  ! interface for coefficent function
  subroutine coefficient_function(this , omega , dt , dxi, if_periodic )
    import t_base_tridiag_solver, p_k_fld, p_double
    implicit none
    class(t_base_tridiag_solver), intent(inout) :: this
    real(p_k_fld), intent(in) :: omega
    real(p_double), intent(in) :: dt , dxi
    logical, intent(in) :: if_periodic
  end subroutine
end interface

! export type for discritization matrix
public :: t_line_nconf
public :: t_tridiag, t_base_tridiag_solver
public :: p_lapack_alg, p_lu_alg, p_thomas_alg
public :: p_tridiag_lower, p_tridiag_middle, p_tridiag_upper

! symbols for communication
public :: root, last

contains

!-----------------------------------------------------------------------------
subroutine setup_tridiag(this, coeff_func, solver_type, no_co, &
                         grid, gc_min, dim, l_shift, u_shift)
  ! used for allocating arrays and deriving communiction and grid objects

  use omp_lib


  use m_node_conf, only : t_node_conf, n_threads!, my_aid
  use m_grid_define, only : t_grid
  use m_grid

  implicit none

  ! input parameters
  class( t_tridiag ), intent(inout) :: this
  procedure(coefficient_function) :: coeff_func
  integer, intent(in) :: solver_type
  class( t_grid ), intent(in) :: grid
  class( t_node_conf ), intent(in) :: no_co
  integer, dimension(1:, 1:), intent(in) :: gc_min
  integer, intent(in) :: dim

  integer, intent(in), optional :: l_shift ! shift of lower index grid
  integer, intent(in), optional :: u_shift ! shift of upper index grid

  ! aux parameters
  integer :: n_nodes, i, nt
  integer, allocatable, dimension(:,:) :: grid_info, nodes_stripes
  integer, dimension(1:2) :: grid_shift

  ! allocates appropriate solver type
  select case( solver_type )
  case( p_thomas_alg )
    allocate( t_thomas_solver :: this%solver )
  case( p_lu_alg )
    allocate( t_lu_solver :: this%solver )
  case( p_lapack_alg )



    write(err_buf__,*) "(* error *) compiled without lapack support!";call err__("pgc/os-tridiag.f03",284)
    write(err_buf__,*) "(* error *) recompile with 'PGC_WITH_LAPACK'";call err__("pgc/os-tridiag.f03",285)
    call abort_program( p_err_invalid )

  case default
    write(err_buf__,*) "(* error *) wrong solver type for tridiagonal solver";call err__("pgc/os-tridiag.f03",289)
    call abort_program( p_err_invalid )
  end select

  ! repack shifts before gathering grid info
  grid_shift = 0

  if(present(l_shift)) then
    grid_shift(1) = l_shift
  endif

  if(present(u_shift)) then
    grid_shift(2) = u_shift
  endif

  ! obtain information of line node configuration
  call allocate_line_nconf( this%nconf, no_co, dim )

  ! obtain informations of grid
  call gather_grid_info( this, no_co, grid, grid_shift, gc_min, dim )

  ! store finite difference coeffcients
  this%solver%coeff_func => coeff_func
  ! allocate discritization matrix objects
  call this%solver%setup( this%lb, this%ub )

  ! number of nodes
  n_nodes = this%nconf%no_num
  this%stripes = get_stripes( this, this%nconf )

  ! right-hand-side
  this%rhs => null()

  if( n_nodes > 1 ) then
    ! allocate messages
    allocate( this%src_msg(1:n_nodes) )
    allocate( this%rhs_msg(1:n_nodes) )

    ! allocate/share grid information in transversal direction
    allocate(grid_info( 1:n_nodes, 1:2 ))
    grid_info = global_grids_info( this, this%nconf )

    ! allocate/share grid information in longitudinal direction
    allocate(nodes_stripes( 1:n_nodes, 1:2 ))
    nodes_stripes = global_stripes_info( this, this%nconf )

    do i = 1, n_nodes
      call allocate_msg( this%src_msg(i), i, nodes_stripes(i, :), &
                         (/ this%ext_min, this%ext_max /))
      call allocate_msg( this%rhs_msg(i), i, (/ 1, this%stripes /), &
                         grid_info(i, :) )
    enddo

    call alloc(this%rhs, (/ 1, this%lb /), (/ this%stripes, this%ub /),"pgc/os-tridiag.f03",342)

    deallocate( grid_info )
    deallocate( nodes_stripes )
  endif

  ! allocate memory for periodic systems - use subsystem
  if( this%nconf%ifpr ) then

    nt = n_threads( no_co )



    call alloc(this%u_vec, (/ this%lb /), (/ this%ub /),"pgc/os-tridiag.f03",355)
  endif
end subroutine setup_tridiag
!-----------------------------------------------------------------------------
subroutine thomas_setup( this, lb, ub )
  ! setup routine for thomas algorithm
  implicit none

  ! input parameters
  class(t_thomas_solver), intent(inout) :: this
  integer, intent(in) :: lb
  integer, intent(in) :: ub

  this%lb = lb
  this%ub = ub
  call alloc(this%gamma, (/ lb /), (/ ub /),"pgc/os-tridiag.f03",370)
  call alloc(this%coeff, (/ p_tridiag_lower, lb /), (/ p_tridiag_upper, ub /),"pgc/os-tridiag.f03",371)
end subroutine thomas_setup
!-----------------------------------------------------------------------------
subroutine lu_setup( this, lb, ub )
  ! setup routine for thomas algorithm
  implicit none

  ! input parameters
  class(t_lu_solver), intent(inout) :: this
  integer, intent(in) :: lb
  integer, intent(in) :: ub

  ! base members for tridiagonal system
  this%lb = lb
  this%ub = ub
  call alloc(this%coeff, (/ p_tridiag_lower, lb /), (/ p_tridiag_upper, ub /),"pgc/os-tridiag.f03",386)

  ! specific members for LU solver
  call alloc(this%u, (/ lb /), (/ ub-2 /),"pgc/os-tridiag.f03",389)
  call alloc(this%ipiv, (/ lb /), (/ ub /),"pgc/os-tridiag.f03",390)

end subroutine lu_setup
!-----------------------------------------------------------------------------
# 411 "pgc/os-tridiag.f03"
!-----------------------------------------------------------------------------
subroutine setup_smoothing_tridiag( this, kernel, pre, post )
  ! setup routine for tridiagonal solver smoothing
  ! TODO: add a general formulation for an arbitary sized window
  ! works at the moment with a kernel with a window of 3

  implicit none

  ! input parameters
  class( t_tridiag ), intent(inout) :: this
  real( p_k_fld ), intent(in), dimension(1:, 1:) :: kernel
  logical, intent(in) :: pre
  logical, intent(in) :: post

  ! auxiliary parameters
  integer :: order, window
  integer :: i, k
  real( p_k_fld ) :: f1, f2, sum

  ! store smoothing information
  this%pre_smoothing = pre
  this%post_smoothing = post

  if( pre .or. post ) then

    ! first dimension is the smoothing order
    order = size(kernel, dim=2)
    window = size(kernel, dim=1)

    this%smoothing_order = order

    ! reduce smoothing kernel
    call alloc(this%smoothing_kernel, (/ -order /), (/ order /),"pgc/os-tridiag.f03",443)
    this%smoothing_kernel = 0.0_p_k_fld
    this%smoothing_kernel(0) = 1.0_p_k_fld

    do i = 1, order
      f1 = 0.0_p_k_fld
      do k = -(i-1), i
        this%smoothing_kernel(k-1) = &
          this%smoothing_kernel(k-1) + kernel(3, i) * this%smoothing_kernel(k)
        f2 = kernel(1, i) * this%smoothing_kernel(k)
        this%smoothing_kernel(k) = kernel(2,i) * this%smoothing_kernel(k) + f1
        f1 = f2
      enddo
    enddo

    ! normalize
    sum = 0.0_p_k_fld
    do i = -order, order
      sum = sum + this%smoothing_kernel(i)
    enddo
    do i = -order, order
      this%smoothing_kernel(i) = this%smoothing_kernel(i) / sum
    enddo

    if( this%global_size <= 0 ) then
      print *, '(*error*) before setting up the smoothing routine for'
      print *, '(*error*) tridiagonal solver, use its setup routine'
      stop
    endif

    call alloc(this%raw_buffer, (/ this%lb - order /), (/ this%ub + order /),"pgc/os-tridiag.f03",473)

  ! check that smothing kernel is symmetric - not active by default
  ! set to 1 if you want to run the test
# 487 "pgc/os-tridiag.f03"
  endif
end subroutine setup_smoothing_tridiag
!-----------------------------------------------------------------------------
subroutine smooth_tridiag( this )
  ! smooth routine for tridiagonal solver
  implicit none

  class( t_tridiag ), intent(inout) :: this

  ! auxiliary parameters
  integer :: jmin, jmax, i, j, k
  complex( p_k_fld ) :: f

  ! length of rhs
  jmin = lbound(this%rhs, dim=2)
  jmax = ubound(this%rhs, dim=2)

  ! for every stipe
  do i = 1, this%stripes
    ! copy inner buffer
    do j = jmin, jmax
      this%raw_buffer(j) = this%rhs(i, j)
    enddo
    ! mirror boundaries
    do k = 1, this%smoothing_order
      this%raw_buffer(jmin - k) = -1.0_p_k_fld * this%rhs(i, jmin + k)
      this%raw_buffer(jmax + k) = -1.0_p_k_fld * this%rhs(i, jmax - k)
    enddo

    ! apply smoothing kernel
    do j = jmin, jmax
      f = 0.0_p_k_fld
      do k = -this%smoothing_order, this%smoothing_order
        f = f + this%raw_buffer(j+k) * this%smoothing_kernel(k)
      enddo
      this%rhs(i, j) = f
    enddo
  enddo
end subroutine smooth_tridiag
!-----------------------------------------------------------------------------
subroutine solve_tridiag( this, source )
  ! solving routine for tridiagonal solver

  implicit none

  ! input parameters
  class( t_tridiag ), intent(inout) :: this
  complex( p_k_fld ), dimension(this%length, this%ext_min:this%ext_max), &
    intent(inout), target :: source

  ! transpose: source -> this%rhs if communication needed else take source
  if (this%nconf%no_num /= 1 ) then
    call pack_msgs( this%src_msg, source, this%length, &
                    this%ext_min, this%ext_max )
    call transpose( this%nconf, this%src_msg, this%rhs_msg )
    call unpack_msgs( this%rhs_msg, this%rhs, this%stripes, &
                      this%lb, this%ub )
  else
    this%rhs => source
  endif

  if (this%pre_smoothing) then
    call this%smooth()
  endif

  ! find solution depending if periodic or non-periodic system
  if( this%nconf%ifpr ) then
    call periodic_solve( this )
  else
    call non_periodic_solve( this )
  endif

  if (this%post_smoothing) then
    call this%smooth()
  endif

  ! retranspose source and this%rhs
  if (this%nconf%no_num /= 1 ) then
    call pack_msgs( this%rhs_msg, this%rhs, this%stripes, &
                    this%lb, this%ub )
    call transpose( this%nconf, this%rhs_msg, this%src_msg )
    call unpack_msgs( this%src_msg, source, this%length, &
                      this%ext_min, this%ext_max )
  else
    this%rhs => null()
  endif
end subroutine solve_tridiag
!-----------------------------------------------------------------------------
subroutine thomas_cleanup( this )
  ! cleanup routine for thomas solver
  implicit none
  ! input parameters
  class( t_thomas_solver ), intent(inout) :: this

  ! tridiagonal system arrays
  call freemem(this%gamma,"pgc/os-tridiag.f03",582)
  call freemem(this%coeff,"pgc/os-tridiag.f03",583)
end subroutine thomas_cleanup
!-----------------------------------------------------------------------------
subroutine lu_cleanup( this )
  ! cleanup routine for lu solver
  implicit none
  ! input parameters
  class(t_lu_solver), intent(inout) :: this

  ! tridiagonal system arrays
  call freemem(this%coeff,"pgc/os-tridiag.f03",593)
  call freemem(this%u,"pgc/os-tridiag.f03",594)
  call freemem(this%ipiv,"pgc/os-tridiag.f03",595)
end subroutine lu_cleanup
!-----------------------------------------------------------------------------
# 609 "pgc/os-tridiag.f03"
!-----------------------------------------------------------------------------
subroutine cleanup_tridiag( this )
  ! cleanup routine for tridiagonal solver

  implicit none

  ! input parameters
  class( t_tridiag ), intent(inout) :: this

  ! deallcoate solver
  call this%solver%cleanup()
  deallocate( this%solver )

  ! cleanup messages
  if( this%nconf%no_num /= 1 ) then
    call cleanup_msg( this%src_msg )
    call cleanup_msg( this%rhs_msg )
  endif

  ! if right-hand-side pointer is set, unset it
  if( this%nconf%no_num /= 1 ) then
    call freemem(this%rhs,"pgc/os-tridiag.f03",630)
    nullify(this%rhs)
  else
    nullify(this%rhs)
  endif

  if( this%pre_smoothing .or. this%post_smoothing ) then
    call freemem(this%smoothing_kernel,"pgc/os-tridiag.f03",637)
    call freemem(this%raw_buffer,"pgc/os-tridiag.f03",638)
    nullify( this%smoothing_kernel, this%raw_buffer )
  endif

  ! periodic helpers
  if( this%nconf%ifpr ) then
    call freemem(this%u_vec,"pgc/os-tridiag.f03",644)
  endif
end subroutine cleanup_tridiag
!-----------------------------------------------------------------------------
subroutine thomas_decompose( this )
  ! decomposition routine for thomas algorithm

  implicit none

  ! input parameters
  class(t_thomas_solver), intent(inout) :: this

  ! auxiliary
  integer :: i

  ! calculate gamma
# 688 "pgc/os-tridiag.f03"
  this%gamma(this%lb) = 1.0_p_k_fld / this%coeff(p_tridiag_middle, this%lb)
  do i = this%lb+1, this%ub
    this%gamma(i) = 1.0_p_k_fld / &
      (this%coeff(p_tridiag_middle, i) - this%coeff(p_tridiag_lower, i) * &
       this%gamma(i-1) * this%coeff( p_tridiag_upper, i ) )
  end do


end subroutine thomas_decompose
!-----------------------------------------------------------------------------
subroutine lu_decompose( this )
  ! decomposition routine for thomas algorithm

  implicit none

  ! input parameters
  class(t_lu_solver), intent(inout) :: this

  ! auxiliary
  integer :: i
  complex(p_k_fld) :: fact, temp

  do i = this%lb , this%ub-1
    this%ipiv(i) = i
  enddo

  do i = this%lb , this%ub - 2
    this%u(i) = 0.0_p_k_fld
  enddo

  do i = this%lb, this%ub-2
    if(abs_val(this%coeff(p_tridiag_middle, i)) >= &
       abs_val(this%coeff(p_tridiag_lower, i))) then
      ! no row interchange required, eliminate dl(i)
      if(abs_val(this%coeff(p_tridiag_middle, i)) /= 0.0_p_k_fld) then
        fact = this%coeff(p_tridiag_lower, i) / this%coeff(p_tridiag_middle, i)
        this%coeff(p_tridiag_lower, i) = fact
        this%coeff(p_tridiag_middle, i+1) = this%coeff(p_tridiag_middle, i+1) - &
          fact * this%coeff(p_tridiag_upper, i)
      endif
    else
      ! interchange rows i and i+1, eliminate dl(i)
      fact = this%coeff(p_tridiag_middle, i) / this%coeff(p_tridiag_lower, i)
      this%coeff(p_tridiag_middle, i) = this%coeff(p_tridiag_lower, i)
      this%coeff(p_tridiag_lower, i) = fact
      temp = this%coeff(p_tridiag_upper, i)
      this%coeff(p_tridiag_upper, i) = this%coeff(p_tridiag_middle, i+1)
      this%coeff(p_tridiag_middle, i+1) = temp - fact*this%coeff(p_tridiag_middle, i+1)
      this%u(i) = this%coeff(p_tridiag_upper, i+1)
      this%coeff(p_tridiag_upper, i+1) = -fact*this%coeff(p_tridiag_upper, i+1)
      this%ipiv(i) = i + 1
    endif
  enddo

  if(this%ub > 1) then
    i = this%ub - 1

    if(abs_val(this%coeff(p_tridiag_middle, i)) >= &
       abs_val(this%coeff(p_tridiag_lower, i))) then
      if(abs_val(this%coeff(p_tridiag_middle, i)) /= 0.0) then
        fact = this%coeff(p_tridiag_lower, i) / this%coeff(p_tridiag_middle, i)
        this%coeff(p_tridiag_lower, i) = fact
        this%coeff(p_tridiag_middle, i+1) = this%coeff(p_tridiag_middle, i+1) - &
          fact * this%coeff(p_tridiag_upper, i)
      endif
    else
      fact = this%coeff(p_tridiag_middle, i) / this%coeff(p_tridiag_lower, i)
      this%coeff(p_tridiag_middle, i) = this%coeff(p_tridiag_lower, i)
      this%coeff(p_tridiag_lower, i) = fact
      temp = this%coeff(p_tridiag_upper, i)
      this%coeff(p_tridiag_upper, i) = this%coeff(p_tridiag_middle, i+1)
      this%coeff(p_tridiag_middle, i+1) = temp - fact*this%coeff(p_tridiag_middle, i+1)
      this%ipiv(i) = i + 1
   endif
  endif

  ! check for a zero on the diagonal of u.
  do i = this%lb, this%ub
    if(abs_val(this%coeff(p_tridiag_middle, i)) == 0.0) then
      print *, '(* error *) diagonal element ', i, 'is zero'
      exit
    endif
  enddo

  ! check diagonal dominance
  do i = this%lb, this%ub
    if(abs(this%coeff(p_tridiag_middle, i)) < &
        (abs(this%coeff(p_tridiag_lower, i) + &
         abs(this%coeff(p_tridiag_upper, i))))) then
      print *, "diagonal dominance is broken at i = ", i
      stop
    endif
  enddo

end subroutine lu_decompose
!-----------------------------------------------------------------------------
# 792 "pgc/os-tridiag.f03"
!-----------------------------------------------------------------------------
function abs_val( zdum )

  implicit none
  complex(p_k_fld), intent(in) :: zdum
  real(p_k_fld) :: abs_val

  abs_val = abs( real( zdum ) ) + abs( aimag( zdum ) )

end function abs_val
!-----------------------------------------------------------------------------
subroutine non_periodic_solve( this )
  ! TODO: modify to consider solver
  implicit none

  ! input parameters
  class( t_tridiag ), intent(inout) :: this

  ! auxiliary parameters
  integer :: i, j, jmin, jmax, stripes
  complex( p_k_fld ), pointer :: rhs(:,:)
  complex( p_k_fld ), allocatable :: rhs_correction(:,:), residual(:)
  complex( p_k_fld ), pointer :: coeff(:,:)

  ! extract info about system - required for shared decleration of openMP
  rhs => this%rhs
  coeff => this%solver%coeff
  stripes = this%stripes

  if(this%error_cleaning) then
    jmin = this%lb
    jmax = this%ub

    ! allocate rhs_correction
    allocate(rhs_correction(1:stripes, jmin:jmax))
    allocate(residual(jmin:jmax))

    ! copy previous rhs
    !$omp parallel do private(j), shared(rhs, rhs_correction, jmin, jmax)
    do i = 1, stripes
      do j = jmin, jmax
        rhs_correction(i,j) = rhs(i,j)
      enddo
    enddo
    !$omp end parallel do

    ! solve with an error
    !$omp parallel do shared(stripes, rhs)
    do i = 1, stripes
      call this%solver%solve(rhs(i, :))
    enddo
    !$omp end parallel do

    ! construct rhs for error equation, find the error and substract it
    !$omp parallel do private(j, residual) &
    !$omp& shared(stripes, rhs, rhs_correction, coeff)
    do i = 1, stripes
      call apply_mat(residual, coeff, rhs(i, :))

      do j = jmin, jmax
        rhs_correction(i,j) = residual(j) - rhs_correction(i,j)
      enddo

      call this%solver%solve(rhs_correction(i, :))

      do j = jmin, jmax
        rhs(i,j) = rhs(i,j) - rhs_correction(i,j)
      enddo
    enddo
    !$omp end parallel do

    deallocate(rhs_correction, residual)
  else
    !$omp parallel do shared(stripes, rhs)
    do i = 1, stripes
      call this%solver%solve(rhs(i, :))
    enddo
    !$omp end parallel do
  endif
end subroutine non_periodic_solve
!-----------------------------------------------------------------------------
subroutine periodic_preCalculation_tridiag( this )
  ! calculate indepent solution for periodic case

  implicit none

  class( t_tridiag ), intent(inout) :: this

  complex( p_k_fld ), pointer :: u_vec(:)
  integer :: jmin, jmax

  u_vec => this%u_vec

  ! access indecies
  jmin = this%lb
  jmax = this%ub

  u_vec = 0.0_p_k_fld
  u_vec( jmin ) = this%solver%gam
  u_vec( jmax ) = this%solver%coeff(p_tridiag_upper, jmax)

  call this%solver%solve(u_vec)
end subroutine periodic_preCalculation_tridiag
!-----------------------------------------------------------------------------
subroutine periodic_solve( this )
  ! TODO: modify to consider solver
  ! Solving a periodc System based on Sherman-Morrison formula


  use omp_lib


  implicit none

  ! input parameters
  class( t_tridiag ), intent(inout) :: this

  complex( p_k_fld ) :: factor
  ! corresponds elements in semi-tridiag matrix
  integer :: i, j, jmin, jmax, stripes, tid, lower
  complex( p_k_fld ), pointer :: rhs(:,:)
  complex( p_k_fld ), pointer :: u_vec(:)
  complex( p_k_fld ), pointer :: coeff(:,:)
  complex( p_k_fld ) :: gam

  rhs => this%rhs
  u_vec => this%u_vec
  coeff => this%solver%coeff
  stripes = this%stripes
  gam = this%solver%gam
  lower = p_tridiag_lower

  ! access indecies
  jmin = this%lb
  jmax = this%ub

  ! slice along x-direction

  tid = omp_get_thread_num()




  !$omp parallel do &
  !$omp& private(factor, j) &
  !$omp& shared(rhs, coeff, lower) &
  !$omp& shared(jmin, jmax, gam, u_vec)
  do i = 1, stripes
    call this%solver%solve(rhs(i, :))
    factor = &
      (rhs(i, jmin) + (coeff(lower, jmin) * rhs(i, jmax)) / gam) / &
      (1.0_p_k_fld + u_vec(jmin) + (coeff(lower, jmin) * u_vec(jmax)) / gam)

    do j = jmin, jmax
      rhs(i, j) = rhs(i, j) - factor * u_vec(j)
    enddo
  enddo
  !$omp end parallel do

  nullify( rhs, u_vec )
end subroutine periodic_solve
!-----------------------------------------------------------------------------
subroutine thomas_solve( this, rhs )
  implicit none

  ! input parameters
  class(t_thomas_solver), intent(in) :: this
  complex(p_k_fld), intent(inout), dimension(this%lb:this%ub) :: rhs

  ! auxiliary variables
  complex( p_k_fld ), dimension(this%lb:this%ub):: tmp
  integer :: i, imin, imax

  imin = this%lb
  imax = this%ub

  ! forward elimination
  tmp(imin) = this%gamma(imin) * rhs(imin)
  do i = imin+1, imax
    tmp(i) = this%gamma(i) * (rhs(i) - this%coeff(p_tridiag_lower, i) * tmp(i-1))
  enddo

  ! backward elimination
  rhs(imax) = tmp(imax)
  do i = imax-1, imin, -1
    rhs(i) = tmp(i) - this%gamma(i) * this%coeff(p_tridiag_upper, i) * rhs(i+1)
  enddo
end subroutine thomas_solve
!-----------------------------------------------------------------------------
subroutine lu_solve( this, rhs )
  implicit none

  ! input parameters
  class(t_lu_solver), intent(in) :: this
  complex(p_k_fld), intent(inout), dimension(this%lb:this%ub) :: rhs

  ! auxiliary variables
  complex(p_k_fld) :: tmp
  integer :: i

  ! Solve L*y = b
  do i = this%lb, this%ub-1
    if(this%ipiv(i) == i) then
      rhs(i+1) = rhs(i+1) - this%coeff(p_tridiag_lower, i) * rhs(i)
    else
      tmp = rhs(i)
      rhs(i) = rhs(i+1)
      rhs(i+1) = tmp - this%coeff(p_tridiag_lower, i) * rhs(i)
    endif
  enddo

  ! Solve U*x = c
  rhs(this%ub) = rhs(this%ub) / this%coeff(p_tridiag_middle, this%ub)

  rhs(this%ub-1) = &
    (rhs(this%ub-1) - this%coeff(p_tridiag_upper, this%ub-1) * rhs(this%ub)) / &
    this%coeff(p_tridiag_middle, this%ub-1)
  do i = this%ub-2, this%lb, -1
    rhs(i) = ( &
        rhs(i) - (this%coeff(p_tridiag_upper, i) * rhs(i+1)) - &
        (this%u(i) * rhs(i+2)) &
      ) / this%coeff(p_tridiag_middle, i)
  enddo

end subroutine lu_solve
!-----------------------------------------------------------------------------
# 1048 "pgc/os-tridiag.f03"
!-----------------------------------------------------------------------------
subroutine apply_mat(b, coeff, x)
  implicit none

  ! input parameters
  complex( p_k_fld ), intent(out), dimension(:) :: b
  complex( p_k_fld ), intent(in), &
    dimension(p_tridiag_lower:p_tridiag_upper, size(b)) :: coeff
  complex( p_k_fld ), intent(inout), dimension(size(b)) :: x

  ! auxiliary variables
  integer :: i, imin, imax

  imin = lbound(b, dim=1)
  imax = ubound(b, dim=1)

  b(imin) = coeff(p_tridiag_middle, imin) * x(imin) + &
            coeff(p_tridiag_upper, imin) * x(imin+1)

  do i = imin+1, imax-1
    b(i) = coeff(p_tridiag_lower, i) * x(i-1) + &
           coeff(p_tridiag_middle, i) * x(i) + &
           coeff(p_tridiag_upper, i) * x(i+1)
  enddo

  b(imax) = coeff(p_tridiag_lower, imax) * x(imax-1) + &
            coeff(p_tridiag_middle, imax) * x(imax)
end subroutine apply_mat
!-----------------------------------------------------------------------------
subroutine transpose( nconf, send_msg, recv_msg )
  implicit none

  ! input parameters
  type( t_line_nconf ), intent(in) :: nconf

  type( t_msg ), dimension(nconf%no_num), intent(inout) :: send_msg
  type( t_msg ), dimension(nconf%no_num), intent(inout) :: recv_msg

  ! auxiliary variables
  integer :: i

  ! communicate buffers else copy buffer from send to recv
  do i = 1, nconf%no_num
    if( i /= my_aid(nconf)) then
      call irecv_msg( recv_msg(i), nconf )
      call isend_msg( send_msg(i), nconf )
    else
      recv_msg(i)%buffer = send_msg(i)%buffer
    endif
  enddo

  call waitall_msg( recv_msg )
  call waitall_msg( send_msg )
end subroutine transpose
!-----------------------------------------------------------------------------
subroutine allocate_line_nconf( nconf, no_co, dim )
  use m_node_conf, only : t_node_conf, comm, periodic, no_num

  implicit none

  ! input/output parameters
  type( t_line_nconf ), intent(inout) :: nconf
  class( t_node_conf ) , intent(in) :: no_co
  integer, intent(in) :: dim

  ! auxiliary
  logical :: if_periodic(p_x_dim)
  logical :: remain_dims(p_x_dim)
  integer :: comm_new = mpi_comm_null
  integer :: new_rank
  integer :: comm_size

  integer :: ierr

  ! construct subgrid array
  remain_dims = .false.
  remain_dims(dim) = .true.

  if ( no_num(no_co) == 1 ) then
    nconf%comm = comm(no_co)
    nconf%no_num = 1
    nconf%my_aid = 1
  else
    ! devide comm world
    call mpi_cart_sub( comm(no_co), remain_dims, comm_new, ierr )
    call check_mpi_error( ierr, "line_conf :: mpi_cart_sub" )
    ! gather communicator informations
    call mpi_comm_size( comm_new, comm_size, ierr )
    call check_mpi_error( ierr, "line_conf :: mpi_comm_size" )
    call mpi_comm_rank( comm_new, new_rank, ierr )
    call check_mpi_error( ierr, "line_conf :: mpi_comm_rank" )

    nconf%comm = comm_new
    nconf%no_num = comm_size
    nconf%my_aid = new_rank + 1
  endif

  ! get if periodic
  if_periodic = periodic( no_co )
  nconf%ifpr = if_periodic( dim )
end subroutine allocate_line_nconf
!-----------------------------------------------------------------------------
subroutine check_mpi_error( ierr, msg )
  implicit none

  integer, intent(in) :: ierr
  character(len=*), intent(in) :: msg

  if ( ierr /= MPI_SUCCESS ) then
    write(err_buf__,*) msg;call err__("pgc/os-tridiag.f03",1157)
    call abort_program( p_err_mpi )
  endif
end subroutine check_mpi_error
!-----------------------------------------------------------------------------
subroutine gather_grid_info( this, no_co, grid, grid_shift, gc_min, dim )
  use m_node_conf, only : t_node_conf!, my_aid
  use m_grid_define, only : t_grid
  use m_grid

  implicit none

  class( t_tridiag ), intent(inout) :: this
  class( t_node_conf ), intent(in) :: no_co
  class( t_grid ), intent(in) :: grid
  integer, dimension(1:2) :: grid_shift
  integer, dimension(1:, 1:), intent(in) :: gc_min
  integer, intent(in) :: dim

  ! auxiliary parameters
  integer, dimension(3,grid%x_dim) :: nx_p_grid

  this%global_size = grid%g_nx(dim) - grid_shift(1) - grid_shift(2)
  ! nx_p_grid = nx_p( grid, no_co, no_co%my_aid() )
  nx_p_grid = nx_p( grid, no_co, no_co%my_aid() )
  this%id_min = nx_p_grid(1, dim)
  this%id_max = nx_p_grid(2, dim)
  this%local_size = nx_p_grid(3, dim)
  this%length = nx_p_grid(3, 1)

  ! indices for external grid - right-hand-side of solution
  this%ext_min = 1
  this%ext_max = this%local_size

  this%gc_num(1) = gc_min(1, dim)
  this%gc_num(2) = gc_min(2, dim)
  this%gc_num(3) = gc_min(1, dim) + gc_min(2, dim)

  ! grid indices
  this%lb = 1 + grid_shift(1)
  this%ub = grid%g_nx(dim) + grid_shift(2)

  if(root(this%nconf)) then
    this%id_min = this%id_min + grid_shift(1)
    this%local_size = this%local_size - grid_shift(1)
    this%ext_min = this%ext_min + grid_shift(1)
  endif

  if (last( this%nconf) ) then
    this%id_max = this%id_max + grid_shift(2)
    this%local_size = this%local_size - grid_shift(2)
    this%ext_max = this%ext_max + grid_shift(2)
  endif

end subroutine gather_grid_info
!-----------------------------------------------------------------------------
function global_grids_info( this, nconf )
  implicit none

  class( t_tridiag ), intent(inout) :: this
  type( t_line_nconf ), intent(in) :: nconf
  integer, dimension(1:nconf%no_num, 1:2) :: global_grids_info

  ! auxiliary parameters
  integer, dimension(1:nconf%no_num, 1:2) :: aux_global_grid_info
  integer :: buffer_size
  integer :: ierr

  ! construct local grid information - first
  global_grids_info = 0
  global_grids_info(my_aid(nconf), p_id_min) = this%id_min
  global_grids_info(my_aid(nconf), p_id_max) = this%id_max

  ! pre-communication setup
  aux_global_grid_info = 0
  buffer_size = 2 * nconf%no_num

  ! reduce global grid information on root
  call mpi_reduce( global_grids_info, aux_global_grid_info, buffer_size, &
                   mpi_integer, mpi_sum, 0, nconf%comm, ierr )

  if( root(nconf) ) then
    global_grids_info = aux_global_grid_info
  endif

  ! broadcast global grid informations
  call mpi_bcast( global_grids_info, buffer_size, mpi_integer, 0, &
                  comm(nconf), ierr )
end function global_grids_info
!-----------------------------------------------------------------------------
function global_stripes_info( this, nconf )
  implicit none

  class( t_tridiag ), intent(inout) :: this
  type( t_line_nconf ), intent(in) :: nconf
  integer, dimension(1:nconf%no_num, 1:2) :: global_stripes_info

  ! auxiliary parameters
  integer, dimension(1:nconf%no_num) :: aux_stripes
  integer :: buffer_size
  integer :: ierr
  integer :: i

  call mpi_gather( this%stripes, 1, mpi_integer, aux_stripes, 1, &
                   mpi_integer, 0, comm(nconf), ierr )

  ! if root calculate stripes index
  if( root(nconf) ) then
    do i = 1, nconf%no_num
      if( i == 1 ) then
        global_stripes_info(i, 1) = 1
        global_stripes_info(i, 2) = aux_stripes(i)
      else
        global_stripes_info(i, 1) = (i - 1) * aux_stripes(i-1) + 1
        global_stripes_info(i, 2) = &
          global_stripes_info(i, 1) + aux_stripes(i) - 1
      endif
    enddo
  endif

  ! broadcast global grid informations
  buffer_size = 2 * nconf%no_num

  call mpi_bcast( global_stripes_info, buffer_size, mpi_integer, 0, &
                  comm(nconf), ierr )
end function global_stripes_info
!-----------------------------------------------------------------------------
function get_stripes( this, line_conf )
  implicit none

  class( t_tridiag ), intent(inout) :: this
  type( t_line_nconf ), intent(in) :: line_conf

  integer :: get_stripes

  ! auxiliary parameters
  integer :: local_stripes, num_stripes

  num_stripes = this%length
  local_stripes = num_stripes / line_conf%no_num

  if (last(line_conf)) then
    get_stripes = num_stripes - (local_stripes * (line_conf%no_num - 1))
  else
    get_stripes = local_stripes
  endif
end function get_stripes
!-----------------------------------------------------------------------------
subroutine allocate_msg( msg, node, dim1, dim2)
  implicit none

  type( t_msg ), intent(inout) :: msg
  integer, intent(in) :: node

  integer, intent(in), dimension(2) :: dim1
  integer, intent(in), dimension(2) :: dim2

  msg%node = node

  select case( p_k_fld )
  case( p_single )
    msg%msg_type = MPI_COMPLEX
  case( p_double )
    msg%msg_type = MPI_DOUBLE_COMPLEX
  end select

  msg%l_dim1 = dim1(1)
  msg%u_dim1 = dim1(2)

  msg%l_dim2 = dim2(1)
  msg%u_dim2 = dim2(2)

  msg%msg_size = (dim1(2) - dim1(1) + 1) * (dim2(2) - dim2(1) + 1)

  allocate( msg%buffer( dim1(1):dim1(2), dim2(1):dim2(2) ) )
end subroutine allocate_msg
!-----------------------------------------------------------------------------
subroutine cleanup_msg( msg )
  implicit none

  type( t_msg ), dimension(:), intent(inout) :: msg

  !auxiliary
  integer :: i

  do i = lbound(msg, 1), ubound(msg, 1)
    deallocate( msg(i)%buffer )
  enddo
end subroutine cleanup_msg
!-----------------------------------------------------------------------------
subroutine pack_msgs( msgs, content, content_length, l_bound, u_bound )
  implicit none

  type( t_msg ), intent(inout), dimension(:) :: msgs
  integer, intent(in) :: content_length, l_bound, u_bound
  complex( p_k_fld ), intent(in), &
    dimension(content_length, l_bound:u_bound) :: content

  ! auxiliary
  integer :: num_messages, l_dim1, u_dim1, l_dim2, u_dim2, i

  num_messages = size( msgs )

  do i = 1, num_messages
    ! extract msg information for accessing right regions of content
    l_dim1 = msgs(i)%l_dim1
    u_dim1 = msgs(i)%u_dim1

    l_dim2 = msgs(i)%l_dim2
    u_dim2 = msgs(i)%u_dim2

    ! use subregion of content for
    msgs(i)%buffer = content( l_dim1:u_dim1, l_dim2:u_dim2 )
  enddo
end subroutine pack_msgs
!-----------------------------------------------------------------------------
subroutine unpack_msgs( msgs, content, content_length, l_bound, u_bound )
  implicit none

  type( t_msg ), intent(in), dimension(:) :: msgs
  integer, intent(in) :: content_length, l_bound, u_bound
  complex( p_k_fld ), intent(inout), &
    dimension(content_length, l_bound:u_bound) :: content

  ! auxiliary
  integer :: num_messages, l_dim1, u_dim1, l_dim2, u_dim2, i

  num_messages = size( msgs )

  do i = 1, num_messages
    ! extract msg information for accessing right regions of content
    l_dim1 = msgs(i)%l_dim1
    u_dim1 = msgs(i)%u_dim1

    l_dim2 = msgs(i)%l_dim2
    u_dim2 = msgs(i)%u_dim2

    ! use subregion of content for
    content( l_dim1:u_dim1, l_dim2:u_dim2 ) = msgs(i)%buffer
  enddo
end subroutine unpack_msgs
!-----------------------------------------------------------------------------
subroutine irecv_msg( recv, nconf )
  implicit none

  type( t_line_nconf ), intent(in) :: nconf
  type( t_msg ), intent(inout) :: recv

  ! auxiliary
  integer :: ierr

  call mpi_irecv( recv%buffer, recv%msg_size, recv%msg_type, recv%node-1, &
                  recv%tag, comm(nconf), recv%handle, ierr )
end subroutine irecv_msg
!-----------------------------------------------------------------------------
subroutine isend_msg( send, nconf )
  implicit none

  type( t_line_nconf ), intent(in) :: nconf
  type( t_msg ), intent(inout) :: send

  ! auxiliary
  integer :: ierr

  call mpi_isend( send%buffer, send%msg_size, send%msg_type, send%node-1, &
                  send%tag, comm(nconf), send%handle, ierr )
end subroutine isend_msg
!-----------------------------------------------------------------------------
subroutine waitall_msg( msgs )
  implicit none

  type( t_msg ), dimension(:), intent(inout) :: msgs

  ! auxiliary
  integer :: i, j, ierr, count
  integer, allocatable, dimension(:,:) :: status
  integer, allocatable, dimension(:) :: request

  count = size( msgs )

  allocate( status( MPI_STATUS_SIZE, count ) )
  allocate( request( count ) )

  j = 0
  do i = 1, count
    if( msgs(i)%handle /= MPI_REQUEST_NULL ) then
      j = j + 1
      request(j) = msgs(i)%handle
    endif
  enddo

  ! check if waiting is required
  if( j > 0 ) then
    call mpi_waitall( j, request, status, ierr )
  endif

  ! clear message handles
  do i = 1, count
    msgs(i)%handle = MPI_REQUEST_NULL
  enddo

  deallocate( status )
  deallocate( request )
end subroutine waitall_msg
!-----------------------------------------------------------------------------
function my_aid( nconf )
  implicit none

  type( t_line_nconf ), intent(in) :: nconf
  integer :: my_aid

  my_aid = nconf%my_aid
end function my_aid
!-----------------------------------------------------------------------------
function root( nconf )
  implicit none

  type( t_line_nconf ), intent(in) :: nconf
  logical :: root

  root = ( nconf%my_aid == 1 )
end function root
!-----------------------------------------------------------------------------
function comm( nconf )
  implicit none

  type( t_line_nconf ), intent(in) :: nconf
  integer :: comm

  comm = nconf%comm
end function comm
!-----------------------------------------------------------------------------
function inner( nconf )
  implicit none

  type( t_line_nconf ), intent(in) :: nconf
  logical :: inner

  inner = (( nconf%my_aid /= 1 ) .and. ( nconf%my_aid /= nconf%no_num ) )
  inner = ( inner .and. (nconf%no_num /= 1 ) )
end function inner
!-----------------------------------------------------------------------------
function last( nconf )
  implicit none

  type( t_line_nconf ), intent(in) :: nconf
  logical :: last

  last = ( nconf%my_aid == nconf%no_num )
end function last
!-----------------------------------------------------------------------------
function reduce_stripes( arr )
  implicit none

  integer, dimension(:), intent(in) :: arr
  integer :: reduce_stripes

  ! auxiliary parameters
  integer :: i, tmp

  tmp = 0

  do i = lbound(arr, 1), ubound(arr, 1) - 1
    tmp = tmp + arr(i)
  enddo

  reduce_stripes = tmp
end function reduce_stripes
end module
