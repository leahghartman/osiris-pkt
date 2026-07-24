# 1 "random/os-random-class-grid.f03"
# 1 "<built-in>" 1
# 1 "<built-in>" 3
# 467 "<built-in>" 3
# 1 "<command line>" 1
# 1 "<built-in>" 2
# 1 "random/os-random-class-grid.f03" 2
!-------------------------------------------------------------------------------
! An Array of Random number generators
!
! Used to specifing a distinct random number generator for each
! spacial cell. This is usefull for making results that are indepedent of
! spatial (e.g. MPI) decomposition. Used mainly for testing.
!
! This class is a 'facade' object where each method calls a corresponding method
! of the 't_random' class.
!-------------------------------------------------------------------------------

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
# 13 "random/os-random-class-grid.f03" 2
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
# 14 "random/os-random-class-grid.f03" 2

module m_random_class_grid

use m_parameters
use m_restart
use m_random_class
use m_grid_define
use m_node_conf

private

! string to id restart data
character(len=*), parameter :: p_rand_rst_id = "grid of random num gens rst data - 0x0000"

! needed in Fortran to cconstruct 'an array of pointers'
type :: t_random_grid_holder
  class( t_random ), pointer :: ptr
end type

type, extends( t_random ) :: t_random_grid

  class(t_random), public, pointer :: base_rng => null()

  integer :: nx = 0
  integer :: ny = 0
  integer :: nz = 0
  integer :: global_min_ix = 0
  integer :: global_min_iy = 0
  integer :: global_min_iz = 0
  integer :: global_nx = 0
  integer :: global_ny = 0
  integer :: global_nz = 0


    type( t_random_grid_holder ), dimension(:), public, pointer :: rngs_ => NULL()






contains

  procedure :: genrand_int32 => int32_random_grid
  procedure :: init_genrand_scalar => init_genrand_scalar
  procedure :: init_genrand_scalar_grid => init_scalar_random_grid

  procedure :: genrand_gaussian => genrand_gaussian_grid

  procedure :: write_checkpoint => write_checkpoint_random_grid
  procedure :: read_checkpoint => read_checkpoint_random_grid

  procedure :: class_name => class_name_random_grid

end type t_random_grid


public :: t_random_grid, t_random_grid_holder

contains

function class_name_random_grid( this )

  implicit none

  class( t_random_grid ), intent(in) :: this
  character( len = p_maxlen_random_class_name ) :: class_name_random_grid

  class_name_random_grid = "Grid o' "//this%base_rng%class_name()

end function class_name_random_grid

!
! Some Macros that compute the offset added to the seed at a cell in the grid
!
# 98 "random/os-random-class-grid.f03"
!-------------------------------------------------------------------------------
! Needed to satisfy the fulfill our abtract parent class (even though we don't use it) 
!-------------------------------------------------------------------------------
subroutine init_genrand_scalar(this, s )
  implicit none
  class( t_random_grid ), intent(inout) :: this
  integer, intent(in) :: s

end subroutine init_genrand_scalar

!-------------------------------------------------------------------------------
! generates a integer random number on the [-2147483648, 2147483647] interval
! (all possible 32 bit integers) -(2**31), (2**31)-1
!-------------------------------------------------------------------------------
function int32_random_grid( this, ix )
  use m_system
  implicit none

  class( t_random_grid ), intent(inout) :: this
  integer, dimension(:), optional :: ix
  integer :: int32_random_grid


    int32_random_grid = this%rngs_( ix(1) ) %ptr% genrand_int32()






end function int32_random_grid

!-------------------------------------------------------------------------------
! generates a random number with a normally distributed deviate with 0 mean
! and unit variance. This routine follows the method found in section 7.2
! (pg. 289) of Numerical Recipes in C, 2nd edition.
!
! This routine uses genrand_real3 that always returns numbers in the (0,1)
! interval.
!
! The module variables iset and gset need to saved at restart.
!
! Note: Watch out this is tricky.. we must intercept this call then call the underlying
! the corresponding method on the rng state t_random objects in the grid.
! BUT THAT'S NOT THE REASON WHY WE MAKE OUR OWN 'genrand_gaussian_grid' for the grid.
! All the mapping to the grid could have been taken care of by our 'int32_random_grid'
! above.. this works for all the other t_random proceedures.. BUT NOT 'genrand_gaussian_grid'.
! It almost works.. but we need one other thing... because 'genrand_gaussian_grid'
! keeps extra state in between calls (it produces two numbers at a time, so it returns one
! and then returns the other without any additional calcualtion the next time called).
! And this extra resides in the t_random objects in the gird. Ao in order to properly restore that
! state, we must call 'genrand_gaussian_grid' on the objects in the grid directly.
!-------------------------------------------------------------------------------
function genrand_gaussian_grid( this, ix )

  implicit none

  class( t_random_grid ), intent(inout) :: this
  integer, dimension(:), optional :: ix

  real(p_double) :: genrand_gaussian_grid


    genrand_gaussian_grid = this%rngs_( ix(1) ) %ptr% genrand_gaussian( ix)






end function

!-------------------------------------------------------------------------------
! Initializes random number generator seed with a scalar integer
!-------------------------------------------------------------------------------
subroutine init_scalar_random_grid(this, s, grid, class_idx )
  use m_grid
  use m_restart

  implicit none

  class( t_random_grid ), intent(inout) :: this
  integer, intent(in) :: s
  class( t_grid ), intent(in) :: grid
  integer, intent(inout) :: class_idx

  integer, dimension(1) :: temp_nx
  integer :: x, y, z
  integer :: seed_offset

  temp_nx = my_nx( grid )

  !
  ! Notice that the 'grid o rngs' is 2 larger then needed (one extra cell at the lower index,
  ! one extra at the upper index). These are 'boundry cells'.
  ! Such 'boundry cells' are 'reflective' edges so that if a an index is 1 cell ouside the bounds of
  ! the edges, then the RNG used is the one that is at the appropaite edge within proper bounds.
  ! This is used, for example, in thermal boundry conditions. Technically, at this time, such
  ! boundry cells are only really needed on nodes that haves edges along the simulation's 
  ! 'physical boundries' But the boundry cells are only pointers and the memory usage is minimal
  ! and worth the cost of extra code complexity for now.
  !


    ! gather configuration info
    this%nx = temp_nx(1)
    this%global_min_ix = my_nx_p_min( grid,1 )
    this%global_nx = g_nx(grid, 1)

    ! allocate the grid.
    allocate(this%rngs_(0:this%nx+1))

    ! allocate and init the rng at each cell.
    do x= 1, this%nx
      seed_offset = ( x+this%global_min_ix-1)
      allocate(this%rngs_(x)%ptr, source=this%base_rng )
      call this%rngs_(x)%ptr%init_genrand_scalar( s + seed_offset )
    end do

    !
    ! map the 'boundry cells'
    !
    this%rngs_(0 )%ptr => this%rngs_(1 )%ptr
    this%rngs_(this%nx+1)%ptr => this%rngs_(this%nx)%ptr
# 323 "random/os-random-class-grid.f03"
end subroutine init_scalar_random_grid

!-------------------------------------------------------------------------------
! Write random number generator state to a restart file
!-------------------------------------------------------------------------------
subroutine write_checkpoint_random_grid( this, restart_handle )

    implicit none

    class( t_random_grid ), intent(in) :: this
    class( t_restart_handle ), intent(inout) :: restart_handle

    character(len=*), parameter :: err_msg = 'error writing restart data for random number generator.'
    integer :: ierr, x, y, z

    call restart_io_write("p_rand_rst_id", p_rand_rst_id, restart_handle, ierr)
    call check_error(ierr,err_msg,p_err_rstwrt,"random/os-random-class-grid.f03",339)

    ! now loop trough aour grid of generators and serialize their restart state.

      do x=1, this%nx
        call this%rngs_(x) %ptr%write_checkpoint( restart_handle )
      enddo
# 362 "random/os-random-class-grid.f03"
end subroutine write_checkpoint_random_grid

!-------------------------------------------------------------------------------
! Write random number generator state to a restart file
!-------------------------------------------------------------------------------
subroutine read_checkpoint_random_grid( this, restart_handle )
    implicit none

    class( t_random_grid ), intent(inout) :: this
    class( t_restart_handle ), intent(in) :: restart_handle

    integer :: ierr, x, y, z
    character(len=len(p_rand_rst_id)) :: rst_id
    character(len=*), parameter :: err_msg = 'error reading restart data for random number generator grid'

    call restart_io_read(rst_id, restart_handle, ierr)
    call check_error(ierr,err_msg,p_err_rstrd,"random/os-random-class-grid.f03",378)

    ! check if restart file is compatible
    if ( rst_id /= p_rand_rst_id) then
        write(err_buf__,*) 'Corrupted restart file, or restart file ';call err__("random/os-random-class-grid.f03",382)
        write(err_buf__,*) 'from incompatible binary (random number generator)';call err__("random/os-random-class-grid.f03",383)
        call abort_program(p_err_rstrd)
    endif



      do x=1, this%nx
        call this%rngs_(x) %ptr%read_checkpoint( restart_handle )
      enddo
# 408 "random/os-random-class-grid.f03"
end subroutine read_checkpoint_random_grid



end module m_random_class_grid
