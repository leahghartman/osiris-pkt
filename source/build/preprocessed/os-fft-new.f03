# 1 "fft/os-fft-new.f03"
# 1 "<built-in>" 1
# 1 "<built-in>" 3
# 467 "<built-in>" 3
# 1 "<command line>" 1
# 1 "<built-in>" 2
# 1 "fft/os-fft-new.f03" 2
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
# 2 "fft/os-fft-new.f03" 2
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
# 3 "fft/os-fft-new.f03" 2

! Only one-dimensional fft wrappers are currently implemented in this module
module m_fft_new

use, intrinsic :: iso_c_binding
use m_system

! use m_math, only : pi

implicit none

real(p_double), parameter, private :: pi = 3.14159265358979323846264_p_double


! define data type of fft/ifft wrappers
# 41 "fft/os-fft-new.f03"
integer, parameter :: type_real = C_DOUBLE
integer, parameter :: type_complex = C_DOUBLE_COMPLEX
# 62 "fft/os-fft-new.f03"
type :: t_array_ptr
  real(p_double), dimension(:), allocatable :: arr
end type


!-------------------------------------------------------------------------------
type :: t_fft_manager

  integer :: rank
  integer, dimension(:), allocatable :: shape
  integer :: size

  type(C_PTR) :: rbuf_cptr, cbuf1_cptr, cbuf2_cptr
  type(C_PTR) :: plan_fft_r2c, plan_ifft_c2r, plan_fft_c2c, plan_ifft_c2c

  real(type_real), dimension(:), pointer :: rbuf
  complex(type_complex), dimension(:), pointer :: cbuf1
  complex(type_complex), dimension(:), pointer :: cbuf2

  real(type_real), dimension(:,:), pointer :: rbuf_2d
  complex(type_complex), dimension(:,:), pointer :: cbuf1_2d
  complex(type_complex), dimension(:,:), pointer :: cbuf2_2d

  type(t_array_ptr), dimension(:), allocatable :: fft_freqs

contains

  procedure :: init => init_fft_manager
  procedure :: cleanup => cleanup_fft_manager
  procedure :: check_fftw_support => check_fftw_support_fft_manager

  procedure :: fft_c2c => fft_c2c_fft_manager
  procedure :: ifft_c2c => ifft_c2c_fft_manager

  procedure :: fft_2d_c2c => fft_2d_c2c_fft_manager
  procedure :: ifft_2d_c2c => ifft_2d_c2c_fft_manager

  ! to do
  ! procedure :: fft_r2c => fft_r2c_fft_manager
  ! procedure :: ifft_c2r => fft_c2r_fft_manager

  procedure :: compute_fft_freqs => compute_fft_freqs_fft_manager

end type t_fft_manager
!-------------------------------------------------------------------------------


public :: t_fft_manager

contains

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! only compile fftw functions if we linked fftw. otherwise
! compile stubs so that the code still compiles (but will not run.)
# 446 "fft/os-fft-new.f03"
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
subroutine init_fft_manager( this,shape_in, rank )

  implicit none

  class( t_fft_manager ), intent(inout) :: this
  integer, dimension(:), intent(in) :: shape_in
  integer, intent(in) :: rank

  call this % check_fftw_support()

end subroutine init_fft_manager
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
subroutine cleanup_fft_manager( this )

  implicit none
  class( t_fft_manager ), intent(inout) :: this

  call this % check_fftw_support()

end subroutine cleanup_fft_manager
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
subroutine compute_fft_freqs_fft_manager( this, spacing )

  implicit none
  class( t_fft_manager ), intent(inout) :: this
  real(p_double), dimension(:), intent(in) :: spacing

  call this % check_fftw_support()

end subroutine compute_fft_freqs_fft_manager
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
subroutine check_fftw_support_fft_manager( this )

  implicit none
  class( t_fft_manager ), intent(in) :: this

  if (mpi_node()==0) print *,"ERROR: attempted to call t_fft_manager manager subroutine, "
  if (mpi_node()==0) print *,"but the code was not compiled with FFTW."

  stop

end subroutine check_fftw_support_fft_manager
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
subroutine fft_c2c_fft_manager( this, x )

  implicit none
  class( t_fft_manager ), intent(inout) :: this
  complex(p_double), intent(inout), dimension(:) :: x

  call this % check_fftw_support()

end subroutine fft_c2c_fft_manager
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
subroutine ifft_c2c_fft_manager( this, x )

  implicit none
  class( t_fft_manager ), intent(inout) :: this
  complex(p_double), intent(inout), dimension(:) :: x

  call this % check_fftw_support()

end subroutine ifft_c2c_fft_manager
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
subroutine fft_2d_c2c_fft_manager( this, x )

  implicit none
  class( t_fft_manager ), intent(inout) :: this
  complex(p_double), intent(inout), dimension(:,:) :: x

  call this % check_fftw_support()

end subroutine fft_2d_c2c_fft_manager
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
subroutine ifft_2d_c2c_fft_manager( this, x )

  implicit none
  class( t_fft_manager ), intent(inout) :: this
  complex(p_double), intent(inout), dimension(:,:) :: x

  call this % check_fftw_support()

end subroutine ifft_2d_c2c_fft_manager
!-------------------------------------------------------------------------------




!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! end check for FFTW_ENABLED

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------


end module m_fft_new
