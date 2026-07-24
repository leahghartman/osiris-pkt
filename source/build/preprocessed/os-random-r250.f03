# 1 "random/os-random-r250.f03"
# 1 "<built-in>" 1
# 1 "<built-in>" 3
# 467 "<built-in>" 3
# 1 "<command line>" 1
# 1 "<built-in>" 2
# 1 "random/os-random-r250.f03" 2
!-------------------------------------------------------------------------------
! F03 Implementation of the R250 Pseudo-random number generator
!
! algorithm from:
! Kirkpatrick, S., and E. Stoll, A Very Fast Shift-Register Sequence Random Number
! Generator, Journal of Computational Physics (1981), v. 40. p. 517
!
! Period = 2^250 - 1
! Speed: ~ 1100 M numbers/s
! Memory: 252 x 32 bit integers (1 kbyte)
!
! see also:
! Maier, W.L., 1991; A Fast Pseudo Random Number Generator,
! Dr. Dobb's Journal, May, pp. 152 - 157
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
# 18 "random/os-random-r250.f03" 2
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
# 19 "random/os-random-r250.f03" 2

module m_random_r250

use m_parameters
use m_restart
use m_random_class

private

! string to id restart data
character(len=*), parameter :: p_rand_rst_id = "random_r250 rst data - 0x0000"

type, extends( t_random ) :: t_random_r250

  ! state data
  integer :: index
  integer, dimension(250) :: buffer

contains

  procedure :: genrand_int32 => r250_int32
  procedure :: init_genrand_scalar => init_r250_scalar

  procedure :: write_checkpoint => write_checkpoint_r250
  procedure :: read_checkpoint => read_checkpoint_r250

  procedure :: class_name => class_name_r250

end type t_random_r250

public :: t_random_r250

contains

function class_name_r250( this )

  implicit none

  class( t_random_r250 ), intent(in) :: this
  character( len = p_maxlen_random_class_name ) :: class_name_r250

  class_name_r250 = "R250 Pseudo-random number generator"

end function class_name_r250

!-------------------------------------------------------------------------------
! generates a integer random number on the [-2147483648, 2147483647] interval
! (all possible 32 bit integers) -(2**31), (2**31)-1
!-------------------------------------------------------------------------------
function r250_int32( this, ix )

  implicit none

  class( t_random_r250 ), intent(inout) :: this
  integer, dimension(:), optional :: ix
  integer :: r250_int32

  integer :: i, j

  ! find indexes in r250 buffer
  i = this%index
  j = i + 103
  if ( j > 250 ) j = j - 250

  ! generate random number
  r250_int32 = ieor( this%buffer(i), this%buffer(j) )
  this%buffer(i) = r250_int32

  ! increment pointer for next time
  this%index = this%index + 1
  if (this%index > 250) this%index = 1

end function r250_int32

!-------------------------------------------------------------------------------
! Initializes random number generator seed with a scalar integer
!-------------------------------------------------------------------------------
subroutine init_r250_scalar(this, s)

  implicit none

  class( t_random_r250 ), intent(inout) :: this
  integer, intent(in) :: s

  integer, parameter :: step = 7

  ! Parameters for 32 bit random number generation
  integer(p_int64) :: ms_bit = int(z'80000000',kind=p_int64)
  integer(p_int64) :: half_range = int(z'40000000',kind=p_int64)
  integer(p_int64) :: all_bits = int(z'FFFFFFFF',kind=p_int64)

  integer :: lcrm_seed, i, k, mask, msb

  ! Fill r250 buffer with BITS-1 bit values
  lcrm_seed = s
  do i = 1, 250
    this%buffer(i) = lcmr( lcrm_seed )
  enddo

  ! Set some MSBs to 1
  do i = 1, 250
    if ( lcmr( lcrm_seed ) > half_range ) then
        this%buffer(i) = ior( int(this%buffer(i), kind = p_int64), ms_bit )
    endif
  enddo

  ! Turn on diagonal bit
  msb = ms_bit

  ! Turn off the leftmost bits
  mask = all_bits

  do i = 1, 32
    ! Select word to operate on
    k = step * (i-1) + 4
    ! Turn off bits left of the diagonal
    this%buffer(k) = iand( this%buffer(k), mask )
    ! Tor on the diagonal bit
    this%buffer(k) = ior( this%buffer(k), msb )
    msb = ishft( msb, -1 )
    mask = ishft( mask, -1 )
  enddo

  this%index = 1

contains

  ! The Linear Congruential Method, the "minimal standard generator"
  ! is used for the initialization of the R250 PRNG
  ! Park & Miller, 1988, Comm of the ACM, 31(10), pp. 1192-1201

  ! this generates numbers in the interval 0 .. 2^31-1 (positive 32 bit ints)

  function lcmr( seed )
    implicit none

    integer :: lcmr
    integer, intent(inout) :: seed

    integer, parameter :: a = 16807
    integer, parameter :: int32_max = 2147483647
    integer, parameter :: q = int32_max / a
    integer, parameter :: r = mod( int32_max, a )
    integer :: hi, lo, test

    hi = seed / q
    lo = mod( seed, q )

    test = a * lo - r * hi

    if ( test > 0 ) then
      seed = test
    else
      seed = test + int32_max
    endif

    lcmr = seed

  end function lcmr


end subroutine init_r250_scalar

!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
! Write random number generator state to a restart file
!-------------------------------------------------------------------------------
subroutine write_checkpoint_r250( this, restart_handle )

    implicit none

    class( t_random_r250 ), intent(in) :: this
    class( t_restart_handle ), intent(inout) :: restart_handle

    character(len=*), parameter :: err_msg = 'error writing restart data for random number generator.'
    integer :: ierr

    call restart_io_write("p_rand_rst_id", p_rand_rst_id, restart_handle, ierr)
    call check_error(ierr,err_msg,p_err_rstwrt,"random/os-random-r250.f03",199)

    ! State data
    call restart_io_write("this%index", this%index, restart_handle, ierr)
    call check_error(ierr,err_msg,p_err_rstwrt,"random/os-random-r250.f03",203)

    call restart_io_write("this%buffer", this%buffer, restart_handle, ierr)
    call check_error(ierr,err_msg,p_err_rstwrt,"random/os-random-r250.f03",206)

    ! Write gaussian distribution data. This is done here for simplicity
    call restart_io_write("this%gaussian_set", this%gaussian_set, restart_handle, ierr)
    call check_error(ierr,err_msg,p_err_rstwrt,"random/os-random-r250.f03",210)

    call restart_io_write("this%new_gaussian", this%new_gaussian, restart_handle, ierr)
    call check_error(ierr,err_msg,p_err_rstwrt,"random/os-random-r250.f03",213)

end subroutine write_checkpoint_r250

!-------------------------------------------------------------------------------
! Write random number generator state to a restart file
!-------------------------------------------------------------------------------
subroutine read_checkpoint_r250( this, restart_handle )
    implicit none

    class( t_random_r250 ), intent(inout) :: this
    class( t_restart_handle ), intent(in) :: restart_handle

    integer :: ierr
    character(len=len(p_rand_rst_id)) :: rst_id
    character(len=*), parameter :: err_msg = 'error reading restart data for random number generator'

    call restart_io_read(rst_id, restart_handle, ierr)
    call check_error(ierr,err_msg,p_err_rstrd,"random/os-random-r250.f03",231)

    ! check if restart file is compatible
    if ( rst_id /= p_rand_rst_id) then
        write(err_buf__,*) 'Corrupted restart file, or restart file ';call err__("random/os-random-r250.f03",235)
        write(err_buf__,*) 'from incompatible binary (random number generator)';call err__("random/os-random-r250.f03",236)
        call abort_program(p_err_rstrd)
    endif

    ! State data
    call restart_io_read(this % index, restart_handle, ierr)
    call check_error(ierr,err_msg,p_err_rstrd,"random/os-random-r250.f03",242)

    call restart_io_read(this % buffer, restart_handle, ierr)
    call check_error(ierr,err_msg,p_err_rstrd,"random/os-random-r250.f03",245)

    ! Read gaussian distribution data. This is done here for simplicity
    call restart_io_read(this%gaussian_set, restart_handle, ierr)
    call check_error(ierr,err_msg,p_err_rstrd,"random/os-random-r250.f03",249)

    call restart_io_read(this%new_gaussian, restart_handle, ierr)
    call check_error(ierr,err_msg,p_err_rstrd,"random/os-random-r250.f03",252)

end subroutine read_checkpoint_r250


end module m_random_r250
