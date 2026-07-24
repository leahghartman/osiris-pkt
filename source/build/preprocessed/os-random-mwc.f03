# 1 "random/os-random-mwc.f03"
# 1 "<built-in>" 1
# 1 "<built-in>" 3
# 467 "<built-in>" 3
# 1 "<command line>" 1
# 1 "<built-in>" 2
# 1 "random/os-random-mwc.f03" 2
!-------------------------------------------------------------------------------
! Marsaglia MWC (multiply-with-carry) PRNG
!
! Period : ?
!
! Memory: 2 x 32 bit integers (8 bytes)
! Speed: ~ 620 M numbers/s
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
# 11 "random/os-random-mwc.f03" 2
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
# 12 "random/os-random-mwc.f03" 2

module m_random_mwc

use m_parameters
use m_restart
use m_random_class

private

! string to id restart data
character(len=*), parameter :: p_rand_rst_id = "random_mwc rst data - 0x0000"

type, extends( t_random ) :: t_random_mwc

    ! state data
    integer :: m_w = 12345
    integer :: m_z = 67890

contains

    procedure :: genrand_int32 => mwc_int32
    procedure :: init_genrand_scalar => init_mwc_scalar

    procedure :: write_checkpoint => write_checkpoint_mwc
    procedure :: read_checkpoint => read_checkpoint_mwc

    procedure :: class_name => class_name_mwc

end type t_random_mwc

public :: t_random_mwc

contains

function class_name_mwc( this )

  implicit none

  class( t_random_mwc ), intent(in) :: this
  character( len = p_maxlen_random_class_name ) :: class_name_mwc

  class_name_mwc = "Multiply With Carry Random Number Generator (MWC)"

end function class_name_mwc

!-------------------------------------------------------------------------------
! generates a integer random number on the [-2147483648, 2147483647] interval
! (all possible 32 bit integers) -(2**31), (2**31)-1
!-------------------------------------------------------------------------------
function mwc_int32( this, ix )

  implicit none

  class( t_random_mwc ), intent(inout) :: this
  integer, dimension(:), optional :: ix
  integer :: mwc_int32

  integer, parameter :: MASK = int(z'ffff') ! least significant 16 bits

  ! Reference C implementation
  ! m_z = 36969 * (m_z & 65535) + (m_z >> 16);
  ! m_w = 18000 * (m_w & 65535) + (m_w >> 16);
  ! return (m_z << 16) + m_w; /* 32-bit result */

  this % m_z = 36969 * iand( this % m_z, MASK ) + ishft( this % m_z, -16 )
  this % m_w = 18000 * iand( this % m_w, MASK ) + ishft( this % m_w, -16 )

  mwc_int32 = ishft( this % m_z, 16 ) + this % m_w

end function mwc_int32

!-------------------------------------------------------------------------------
! Initializes random number generator seed with a scalar integer
!-------------------------------------------------------------------------------
subroutine init_mwc_scalar(this, s)

  implicit none

  class( t_random_mwc ), intent(inout) :: this
  integer, intent(in) :: s

  this % m_z = ieor( int(z'80000000'), s )

  ! In the unlikely event that a user uses s = 0x80000000 set a default seed
  ! to prevent a bad sequence
  if( this % m_z == 0 ) this % m_z = 123456789

  this % m_w = 1812433253 * ieor( this % m_z, ishft(this % m_z, -30) )

end subroutine init_mwc_scalar


!-------------------------------------------------------------------------------
! Write random number generator state to a restart file
!-------------------------------------------------------------------------------
subroutine write_checkpoint_mwc( this, restart_handle )

    implicit none

    class( t_random_mwc ), intent(in) :: this
    class( t_restart_handle ), intent(inout) :: restart_handle

    character(len=*), parameter :: err_msg = 'error writing restart data for random number generator.'
    integer :: ierr

    call restart_io_write("p_rand_rst_id", p_rand_rst_id, restart_handle, ierr)
    call check_error(ierr,err_msg,p_err_rstwrt,"random/os-random-mwc.f03",118)

    ! State data
    call restart_io_write("this%m_w", this%m_w, restart_handle, ierr)
    call check_error(ierr,err_msg,p_err_rstwrt,"random/os-random-mwc.f03",122)

    call restart_io_write("this%m_z", this%m_z, restart_handle, ierr)
    call check_error(ierr,err_msg,p_err_rstwrt,"random/os-random-mwc.f03",125)

    ! Write gaussian distribution data. This is done here for simplicity
    call restart_io_write("this%gaussian_set", this%gaussian_set, restart_handle, ierr)
    call check_error(ierr,err_msg,p_err_rstwrt,"random/os-random-mwc.f03",129)

    call restart_io_write("this%new_gaussian", this%new_gaussian, restart_handle, ierr)
    call check_error(ierr,err_msg,p_err_rstwrt,"random/os-random-mwc.f03",132)

end subroutine write_checkpoint_mwc

!-------------------------------------------------------------------------------
! Write random number generator state to a restart file
!-------------------------------------------------------------------------------
subroutine read_checkpoint_mwc( this, restart_handle )
    implicit none

    class( t_random_mwc ), intent(inout) :: this
    class( t_restart_handle ), intent(in) :: restart_handle

    integer :: ierr
    character(len=len(p_rand_rst_id)) :: rst_id
    character(len=*), parameter :: err_msg = 'error reading restart data for random number generator'

    call restart_io_read(rst_id, restart_handle, ierr)
    call check_error(ierr,err_msg,p_err_rstrd,"random/os-random-mwc.f03",150)

    ! check if restart file is compatible
    if ( rst_id /= p_rand_rst_id) then
        write(err_buf__,*) 'Corrupted restart file, or restart file ';call err__("random/os-random-mwc.f03",154)
        write(err_buf__,*) 'from incompatible binary (random number generator)';call err__("random/os-random-mwc.f03",155)
        call abort_program(p_err_rstrd)
    endif

    ! State data
    call restart_io_read(this % m_w, restart_handle, ierr)
    call check_error(ierr,err_msg,p_err_rstrd,"random/os-random-mwc.f03",161)

    call restart_io_read(this % m_z, restart_handle, ierr)
    call check_error(ierr,err_msg,p_err_rstrd,"random/os-random-mwc.f03",164)

    ! Read gaussian distribution data. This is done here for simplicity
    call restart_io_read(this%gaussian_set, restart_handle, ierr)
    call check_error(ierr,err_msg,p_err_rstrd,"random/os-random-mwc.f03",168)

    call restart_io_read(this%new_gaussian, restart_handle, ierr)
    call check_error(ierr,err_msg,p_err_rstrd,"random/os-random-mwc.f03",171)

end subroutine read_checkpoint_mwc


end module m_random_mwc
