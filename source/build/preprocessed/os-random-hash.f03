# 1 "random/os-random-hash.f03"
# 1 "<built-in>" 1
# 1 "<built-in>" 3
# 467 "<built-in>" 3
# 1 "<command line>" 1
# 1 "<built-in>" 2
# 1 "random/os-random-hash.f03" 2
!-------------------------------------------------------------------------------
! Implements the Hashing RNG from Numerical Recipes Chapter 7.5
! Honestly Josh doesn't think it's very good, but it's standard and cheap
!
! Josh, 2016
!
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
# 10 "random/os-random-hash.f03" 2
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
# 11 "random/os-random-hash.f03" 2

module m_random_hash

use m_parameters
use m_restart
use m_random_class

private

! string to id restart data
character(len=*), parameter :: p_rand_rst_id = "random_hash rst data - 0x0000"

type, extends( t_random ) :: t_random_hash

    ! state data
    integer :: seed
    integer :: index

contains

    procedure :: genrand_int32 => hash_int32
    procedure :: init_genrand_scalar => init_hash_scalar

    procedure :: write_checkpoint => write_checkpoint_hash
    procedure :: read_checkpoint => read_checkpoint_hash

    procedure :: class_name => class_name_hash

    procedure :: has_cycled

end type t_random_hash

public :: t_random_hash

contains

function has_cycled( this )

  implicit none

  class( t_random_hash ), intent(in) :: this
  logical :: has_cycled

  ! Technically this actually only means it has half-cycled, but that's hackish to begin with;
  ! and this way the test can be done only occasionally with a decent safety factor
  has_cycled = (this%index .lt. 0)

end function has_cycled


function class_name_hash( this )

  implicit none

  class( t_random_hash ), intent(in) :: this
  character( len = p_maxlen_random_class_name ) :: class_name_hash

  class_name_hash = "Hashing RNG"

end function class_name_hash


!-------------------------------------------------------------------------------
! generates a integer random number on the [-2147483648, 2147483647] interval
! (all possible 32 bit integers) -(2**31), (2**31)-1
!-------------------------------------------------------------------------------
function hash_int32( this, ix )

  implicit none

  class( t_random_hash ), intent(inout) :: this
  integer, dimension(:), optional :: ix
  integer :: hash_int32

  integer, parameter :: NITER = 4
  !integer, dimension(NITER), parameter :: c1 = (/ -z'45569779', z'1e17d32c', z'03bcdc3c', z'0f33d1b2' /)
  !integer, dimension(NITER), parameter :: c2 = (/ z'4b0f3b58', -z'178b0f3d', z'6955c5a6', z'55a7ca46' /)

  integer, dimension(NITER), parameter :: c1 = (/ -1163302777, 504877868, 62708796, 255054258 /)
  integer, dimension(NITER), parameter :: c2 = (/ 1259289432, -394989373, 1767228838, 1437059654 /)

  integer, parameter :: mask = int(z'ffff')


  integer :: i, ia, ib, iswap, itmph, itmpl, lword, irword

  lword = this%seed
  irword = this%index

  do i = 1, NITER
    iswap = irword
    ia = ieor(irword, c1(i))
    itmpl = iand(ia, mask)
    itmph = ishft(ia, -16)
    ib = itmpl*itmpl + not(itmph*itmph)
    ia = ishft(ib, -16)
    ib = ishft(ib, 16)
    irword = ieor(ieor((ia + ib), c2(i)) + itmpl*itmph, lword)
    lword = iswap
  enddo

  this%index = this%index + 1

  hash_int32 = irword


end function hash_int32

!-------------------------------------------------------------------------------
! Initializes random number generator seed with a scalar integer
!-------------------------------------------------------------------------------
subroutine init_hash_scalar(this, s)

  implicit none

  class( t_random_hash ), intent(inout) :: this
  integer, intent(in) :: s

  this%seed = s
  this%index = 1

end subroutine init_hash_scalar


!-------------------------------------------------------------------------------
! Write random number generator state to a restart file
!-------------------------------------------------------------------------------
subroutine write_checkpoint_hash( this, restart_handle )

    implicit none

    class( t_random_hash ), intent(in) :: this
    class( t_restart_handle ), intent(inout) :: restart_handle

    character(len=*), parameter :: err_msg = 'error writing restart data for random number generator.'
    integer :: ierr

    call restart_io_write("p_rand_rst_id", p_rand_rst_id, restart_handle, ierr)
    call check_error(ierr,err_msg,p_err_rstwrt,"random/os-random-hash.f03",149)

    call restart_io_write("this%seed", this%seed, restart_handle, ierr)
    call check_error(ierr,err_msg,p_err_rstwrt,"random/os-random-hash.f03",152)

    call restart_io_write("this%index", this%index, restart_handle, ierr)
    call check_error(ierr,err_msg,p_err_rstwrt,"random/os-random-hash.f03",155)

    ! Write gaussian distribution data. This is done here for simplicity
    call restart_io_write("this%gaussian_set", this%gaussian_set, restart_handle, ierr)
    call check_error(ierr,err_msg,p_err_rstwrt,"random/os-random-hash.f03",159)

    call restart_io_write("this%new_gaussian", this%new_gaussian, restart_handle, ierr)
    call check_error(ierr,err_msg,p_err_rstwrt,"random/os-random-hash.f03",162)

end subroutine write_checkpoint_hash

!-------------------------------------------------------------------------------
! Write random number generator state to a restart file
!-------------------------------------------------------------------------------
subroutine read_checkpoint_hash( this, restart_handle )
    implicit none

    class( t_random_hash ), intent(inout) :: this
    class( t_restart_handle ), intent(in) :: restart_handle

    integer :: ierr
    character(len=len(p_rand_rst_id)) :: rst_id
    character(len=*), parameter :: err_msg = 'error reading restart data for random number generator'

    call restart_io_read(rst_id, restart_handle, ierr)
    call check_error(ierr,err_msg,p_err_rstrd,"random/os-random-hash.f03",180)

    ! check if restart file is compatible
    if ( rst_id /= p_rand_rst_id) then
        write(err_buf__,*) 'Corrupted restart file, or restart file ';call err__("random/os-random-hash.f03",184)
        write(err_buf__,*) 'from incompatible binary (random number generator)';call err__("random/os-random-hash.f03",185)
        call abort_program(p_err_rstrd)
    endif

    call restart_io_read(this % seed, restart_handle, ierr)
    call check_error(ierr,err_msg,p_err_rstrd,"random/os-random-hash.f03",190)

    call restart_io_read(this % index, restart_handle, ierr)
    call check_error(ierr,err_msg,p_err_rstrd,"random/os-random-hash.f03",193)

    ! Read gaussian distribution data. This is done here for simplicity
    call restart_io_read(this%gaussian_set, restart_handle, ierr)
    call check_error(ierr,err_msg,p_err_rstrd,"random/os-random-hash.f03",197)

    call restart_io_read(this%new_gaussian, restart_handle, ierr)
    call check_error(ierr,err_msg,p_err_rstrd,"random/os-random-hash.f03",200)

end subroutine read_checkpoint_hash


end module m_random_hash
