# 1 "random/os-random-kiss.f03"
# 1 "<built-in>" 1
# 1 "<built-in>" 3
# 467 "<built-in>" 3
# 1 "<command line>" 1
# 1 "<built-in>" 2
# 1 "random/os-random-kiss.f03" 2
!-------------------------------------------------------------------------------
! KISS Pseudo Random Number Generator
! Marsaglia, George (2003) "Random Number Generators," Journal of Modern Applied
! Statistical Methods: Vol. 2: Iss. 1, Article 2.
!
! Period > 2^124 ~ 10^37
!
! Memory: 4 x 32 bit integers (16 bytes)
! Speed: ~ 520 M numbers/s
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
# 13 "random/os-random-kiss.f03" 2
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
# 14 "random/os-random-kiss.f03" 2

module m_random_kiss

use m_parameters
use m_restart
use m_random_class

private

! string to id restart data
character(len=*), parameter :: p_rand_rst_id = "random_kiss rst data - 0x0000"

type, extends( t_random ) :: t_random_kiss

  ! state data
  integer :: x = 123456789
  integer :: y = 362436000
  integer :: z = 521288629
  integer :: c = 7654321

  contains

  procedure :: genrand_int32 => kiss_int32
  procedure :: init_genrand_scalar => init_kiss_scalar

  procedure :: write_checkpoint => write_checkpoint_kiss
  procedure :: read_checkpoint => read_checkpoint_kiss

  procedure :: class_name => class_name_kiss

end type t_random_kiss

public :: t_random_kiss

contains

function class_name_kiss( this )

  implicit none

  class( t_random_kiss ), intent(in) :: this
  character( len = p_maxlen_random_class_name ) :: class_name_kiss

  class_name_kiss = "KISS"

end function class_name_kiss

!-------------------------------------------------------------------------------
! generates a integer random number on the [-2147483648, 2147483647] interval
! (all possible 32 bit integers) -(2**31), (2**31)-1
!-------------------------------------------------------------------------------
function kiss_int32( this, ix )

  implicit none

  class( t_random_kiss ), intent(inout) :: this
  integer, dimension(:), optional :: ix
  integer :: kiss_int32

  integer, parameter :: p_int64 = selected_int_kind(10)

  integer(p_int64) :: t

  this%x = 69069 * this%x + 12345
  this%y = ieor( this%y, ishft(this%y, +13))
  this%y = ieor( this%y, ishft(this%y, -17))
  this%y = ieor( this%y, ishft(this%y, +5))

  t = 698769069_p_int64 * this%z + this%c
  this%z = t
  this%c = ishft( t, +32 )

  kiss_int32 = this%x + this%y + this%z

end function kiss_int32

!-------------------------------------------------------------------------------
! Initializes random number generator seed with a scalar integer
!-------------------------------------------------------------------------------
subroutine init_kiss_scalar(this, s)

  implicit none

  class( t_random_kiss ), intent(inout) :: this
  integer, intent(in) :: s

  this%x = ieor( int(z'80000000'), s )
  this%y = 362436000
  this%z = 521288629
  this%c = 1812433253 * ieor( this%x, ishft(this%x, -30) )

end subroutine init_kiss_scalar


!-------------------------------------------------------------------------------
! Write random number generator state to a restart file
!-------------------------------------------------------------------------------
subroutine write_checkpoint_kiss( this, restart_handle )

    implicit none

    class( t_random_kiss ), intent(in) :: this
    class( t_restart_handle ), intent(inout) :: restart_handle

    character(len=*), parameter :: err_msg = 'error writing restart data for random number generator.'
    integer :: ierr

    ! Restart tag
    call restart_io_write("p_rand_rst_id", p_rand_rst_id, restart_handle, ierr)
    call check_error(ierr,err_msg,p_err_rstwrt,"random/os-random-kiss.f03",123)

    ! State data
    call restart_io_write("this%x", this%x, restart_handle, ierr)
    call check_error(ierr,err_msg,p_err_rstwrt,"random/os-random-kiss.f03",127)

    call restart_io_write("this%y", this%y, restart_handle, ierr)
    call check_error(ierr,err_msg,p_err_rstwrt,"random/os-random-kiss.f03",130)

    call restart_io_write("this%z", this%z, restart_handle, ierr)
    call check_error(ierr,err_msg,p_err_rstwrt,"random/os-random-kiss.f03",133)

    call restart_io_write("this%c", this%c, restart_handle, ierr)
    call check_error(ierr,err_msg,p_err_rstwrt,"random/os-random-kiss.f03",136)

    ! Write gaussian distribution data. This is done here for simplicity
    call restart_io_write("this%gaussian_set", this%gaussian_set, restart_handle, ierr)
    call check_error(ierr,err_msg,p_err_rstwrt,"random/os-random-kiss.f03",140)

    call restart_io_write("this%new_gaussian", this%new_gaussian, restart_handle, ierr)
    call check_error(ierr,err_msg,p_err_rstwrt,"random/os-random-kiss.f03",143)

end subroutine write_checkpoint_kiss

!-------------------------------------------------------------------------------
! Write random number generator state to a restart file
!-------------------------------------------------------------------------------
subroutine read_checkpoint_kiss( this, restart_handle )
    implicit none

    class( t_random_kiss ), intent(inout) :: this
    class( t_restart_handle ), intent(in) :: restart_handle

    integer :: ierr
    character(len=len(p_rand_rst_id)) :: rst_id
    character(len=*), parameter :: err_msg = 'error reading restart data for random number generator'

    call restart_io_read(rst_id, restart_handle, ierr)
    call check_error(ierr,err_msg,p_err_rstrd,"random/os-random-kiss.f03",161)

    ! check if restart file is compatible
    if ( rst_id /= p_rand_rst_id) then
        write(err_buf__,*) 'Corrupted restart file, or restart file ';call err__("random/os-random-kiss.f03",165)
        write(err_buf__,*) 'from incompatible binary (random number generator)';call err__("random/os-random-kiss.f03",166)
        call abort_program(p_err_rstrd)
    endif

    ! State data
    call restart_io_read(this % x, restart_handle, ierr)
    call check_error(ierr,err_msg,p_err_rstrd,"random/os-random-kiss.f03",172)

    call restart_io_read(this % y, restart_handle, ierr)
    call check_error(ierr,err_msg,p_err_rstrd,"random/os-random-kiss.f03",175)

    call restart_io_read(this % z, restart_handle, ierr)
    call check_error(ierr,err_msg,p_err_rstrd,"random/os-random-kiss.f03",178)

    call restart_io_read(this % c, restart_handle, ierr)
    call check_error(ierr,err_msg,p_err_rstrd,"random/os-random-kiss.f03",181)

    ! Read gaussian distribution data. This is done here for simplicity
    call restart_io_read(this%gaussian_set, restart_handle, ierr)
    call check_error(ierr,err_msg,p_err_rstrd,"random/os-random-kiss.f03",185)

    call restart_io_read(this%new_gaussian, restart_handle, ierr)
    call check_error(ierr,err_msg,p_err_rstrd,"random/os-random-kiss.f03",188)

end subroutine read_checkpoint_kiss


end module m_random_kiss
