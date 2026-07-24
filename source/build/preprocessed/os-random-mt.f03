# 1 "random/os-random-mt.f03"
# 1 "<built-in>" 1
# 1 "<built-in>" 3
# 467 "<built-in>" 3
# 1 "<command line>" 1
# 1 "<built-in>" 2
# 1 "random/os-random-mt.f03" 2
!-------------------------------------------------------------------------------
! Fortran implementation of the Mersenne Twister random number generator
! (period 2**19937-1). For details on the random number generator see:
!
! http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html
!
! and references therein.
!
! Ported to fortran from the distribution mt19937ar.c
!
! Period = 2^19937-1
!
! Memory: 625 x 32 bit integers (2.5 kb)
! Speed: ~ 1078 M numbers/s
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
# 18 "random/os-random-mt.f03" 2
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
# 19 "random/os-random-mt.f03" 2

module m_random_mt

use m_parameters
use m_restart
use m_random_class

private

! period parameters
integer, parameter :: N = 624
integer, parameter :: M = 397

integer, parameter :: MATRIX_A = int(z'9908b0df') ! constant vector a
integer, parameter :: UPPER_MASK = int(z'80000000') ! most significant w-r bits
integer, parameter :: LOWER_MASK = int(z'7fffffff') ! least significant r bits


! default seed
integer, parameter :: DEFAULT_SEED = 5489

! mag01(x) = x * MATRIX_A for x=0,1
integer, parameter, dimension(0:1) :: mag01 = (/ 0, MATRIX_A /)


! string to id restart data
character(len=*), parameter :: p_rand_rst_id = "random_mt rst data - 0x0000"

type, extends( t_random ) :: t_random_mt

    ! state data
    integer, dimension(N) :: mt ! the array for the state vector
    integer :: mti

contains

    procedure :: genrand_int32 => mt_int32
    procedure :: init_genrand_scalar => init_mt_scalar

    procedure :: write_checkpoint => write_checkpoint_mt
    procedure :: read_checkpoint => read_checkpoint_mt

    procedure :: class_name => class_name_mt

end type t_random_mt

public :: t_random_mt

contains

function class_name_mt( this )

  implicit none

  class( t_random_mt ), intent(in) :: this
  character( len = p_maxlen_random_class_name ) :: class_name_mt

  class_name_mt = "Mersenne Twister"

end function class_name_mt

!-------------------------------------------------------------------------------
! generates a integer random number on the [-2147483648, 2147483647] interval
! (all possible 32 bit integers) -(2**31), (2**31)-1
!-------------------------------------------------------------------------------
function mt_int32( this, ix )

  implicit none

  class( t_random_mt ), intent(inout) :: this
  integer, dimension(:), optional :: ix
  integer :: mt_int32
  integer :: y

  integer :: kk

  if (this%mti > N) then ! generate N words at one time

    do kk = 1, N-M
      y = ior(iand(this%mt(kk),UPPER_MASK),iand(this%mt(kk+1),LOWER_MASK))
      this%mt(kk) = ieor(ieor(this%mt(kk + M),ishft( y, -1)), mag01(iand(y,1)))
    enddo

    do kk = N-M+1, N-1
      y = ior(iand(this%mt(kk),UPPER_MASK),iand(this%mt(kk+1),LOWER_MASK))
      this%mt(kk) = ieor(ieor(this%mt(kk+(M-N)),ishft(y,-1)), mag01(iand(y,1)))
    enddo

    y = ior(iand(this%mt(N),UPPER_MASK),iand(this%mt(1),LOWER_MASK))
    this%mt(N) = ieor(ieor(this%mt(M),ishft(y,-1)), mag01(iand(y,1)))

    this % mti = 1
  endif

  y = this%mt( this%mti )
  this%mti = this%mti + 1

  ! Tempering
  y = ieor(y, ishft( y, -11))
  y = ieor(y, iand( ishft( y, 7 ), int(z'9d2c5680')))
  y = ieor(y, iand( ishft( y, 15 ), int(z'efc60000')))
  y = ieor(y, ishft( y, -18))

  mt_int32 = y

end function mt_int32

!-------------------------------------------------------------------------------
! Initializes random number generator seed with a scalar integer
!-------------------------------------------------------------------------------
subroutine init_mt_scalar(this, s)

  implicit none

  class( t_random_mt ), intent(inout) :: this
  integer, intent(in) :: s
  integer :: i

  this%mt(1) = ieor( s, int(z'80000000') )
  do i = 2, N

    this%mt(i) = (1812433253 * ieor( this%mt(i-1), ishft(this%mt(i-1), -30)) + i-1)
    ! See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier.
    ! In the previous versions, MSBs of the seed affect
    !only MSBs of the array mt[].

  enddo

  this%mti = N

end subroutine init_mt_scalar

!-------------------------------------------------------------------------------
! Write random number generator state to a restart file
!-------------------------------------------------------------------------------
subroutine write_checkpoint_mt( this, restart_handle )

    implicit none

    class( t_random_mt ), intent(in) :: this
    class( t_restart_handle ), intent(inout) :: restart_handle

    character(len=*), parameter :: err_msg = 'error writing restart data for random number generator.'
    integer :: ierr

    call restart_io_write("p_rand_rst_id", p_rand_rst_id, restart_handle, ierr)
    call check_error(ierr,err_msg,p_err_rstwrt,"random/os-random-mt.f03",165)

    ! State data
    call restart_io_write("this%mti", this%mti, restart_handle, ierr)
    call check_error(ierr,err_msg,p_err_rstwrt,"random/os-random-mt.f03",169)

    call restart_io_write("this%mt", this%mt, restart_handle, ierr)
    call check_error(ierr,err_msg,p_err_rstwrt,"random/os-random-mt.f03",172)

    ! Write gaussian distribution data. This is done here for simplicity
    call restart_io_write("this%gaussian_set", this%gaussian_set, restart_handle, ierr)
    call check_error(ierr,err_msg,p_err_rstwrt,"random/os-random-mt.f03",176)

    call restart_io_write("this%new_gaussian", this%new_gaussian, restart_handle, ierr)
    call check_error(ierr,err_msg,p_err_rstwrt,"random/os-random-mt.f03",179)

end subroutine write_checkpoint_mt

!-------------------------------------------------------------------------------
! Write random number generator state to a restart file
!-------------------------------------------------------------------------------
subroutine read_checkpoint_mt( this, restart_handle )
    implicit none

    class( t_random_mt ), intent(inout) :: this
    class( t_restart_handle ), intent(in) :: restart_handle

    integer :: ierr
    character(len=len(p_rand_rst_id)) :: rst_id
    character(len=*), parameter :: err_msg = 'error reading restart data for random number generator'

    call restart_io_read(rst_id, restart_handle, ierr)
    call check_error(ierr,err_msg,p_err_rstrd,"random/os-random-mt.f03",197)

    ! check if restart file is compatible
    if ( rst_id /= p_rand_rst_id) then
        write(err_buf__,*) 'Corrupted restart file, or restart file ';call err__("random/os-random-mt.f03",201)
        write(err_buf__,*) 'from incompatible binary (random number generator)';call err__("random/os-random-mt.f03",202)
        call abort_program(p_err_rstrd)
    endif

    ! State data
    call restart_io_read(this % mti, restart_handle, ierr)
    call check_error(ierr,err_msg,p_err_rstrd,"random/os-random-mt.f03",208)

    call restart_io_read(this % mt, restart_handle, ierr)
    call check_error(ierr,err_msg,p_err_rstrd,"random/os-random-mt.f03",211)

    ! Read gaussian distribution data. This is done here for simplicity
    call restart_io_read(this%gaussian_set, restart_handle, ierr)
    call check_error(ierr,err_msg,p_err_rstrd,"random/os-random-mt.f03",215)

    call restart_io_read(this%new_gaussian, restart_handle, ierr)
    call check_error(ierr,err_msg,p_err_rstrd,"random/os-random-mt.f03",218)

end subroutine read_checkpoint_mt


end module m_random_mt
