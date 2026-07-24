# 1 "fft/os-fft.f03"
# 1 "<built-in>" 1
# 1 "<built-in>" 3
# 467 "<built-in>" 3
# 1 "<command line>" 1
# 1 "<built-in>" 2
# 1 "fft/os-fft.f03" 2
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
# 2 "fft/os-fft.f03" 2
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
# 3 "fft/os-fft.f03" 2

! Only one-dimensional fft wrappers are currently implemented in this module
module m_fft

use, intrinsic :: iso_c_binding
use m_system

implicit none

! define data type of fft/ifft wrappers
# 30 "fft/os-fft.f03"
integer, parameter :: type_real = C_DOUBLE
integer, parameter :: type_complex = C_DOUBLE_COMPLEX
# 60 "fft/os-fft.f03"
!-------------------------------------------------------------------------------
! fft object with the same routines as the original module
!-------------------------------------------------------------------------------
type :: t_fft

  type(C_PTR) :: plan_fft_x, plan_fft_y, plan_fft_z, &
                  plan_ifft_x, plan_ifft_y, plan_ifft_z, &
                  ptr_real_x, ptr_real_y, ptr_real_z, &
                  ptr_complex_x, ptr_complex_y, ptr_complex_z

  real(type_real), dimension(:), pointer :: fftw_real_x => null(), fftw_real_y => null(),&
                                            fftw_real_z => null()

  complex(type_complex), dimension(:), pointer :: fftw_complex_x => null(), &
                                                  fftw_complex_y => null(), &
                                                  fftw_complex_z => null()

contains

  procedure :: fftw_init_1d
  procedure :: fftw_init2_1d
  generic :: init => fftw_init_1d, fftw_init2_1d

  procedure :: fftw_r2c_1d2
  procedure :: fftw_r2c_1d3
  generic :: fft => fftw_r2c_1d2, fftw_r2c_1d3

  procedure :: fftw_idft_1d0
  procedure :: fftw_idft_1d1
  procedure :: fftw_idft_1d2
  procedure :: fftw_c2r_1d2
  procedure :: fftw_c2r_1d3
  generic :: ifft => fftw_idft_1d0, fftw_idft_1d1, fftw_idft_1d2, fftw_c2r_1d2, &
                       fftw_c2r_1d3

  procedure :: cleanup => fftw_cleanup_1d

end type t_fft
!-------------------------------------------------------------------------------

logical :: initialized = .false.

private
public :: t_fft, type_real, type_complex

contains
# 583 "fft/os-fft.f03"
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! NOTE: The following 6 routines are a stopgap fix to allow this code to compile
! on systems without FFTW. It will be removed once a permanent FFT
! library in introduced.
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

subroutine err_out( name )
  implicit none
  character(*) :: name

  write(0,*) "An FFT routine was called but the Osiris was not compiled with FFT "
  write(0,*) "support via the FFTW library. "
  write(0,*) "Please see the FFTW_DIR parameter in your system's configuration file. "
  write(0,*) "The routine called was: ", trim(name)
  call abort_program( )

end subroutine
!-------------------------------------------------------------------------------
subroutine fftw_init_1d( me, n )
!-------------------------------------------------------------------------------
  use, intrinsic :: iso_c_binding
  implicit none

  ! dummy variables
  class( t_fft ), intent(inout) :: me
  integer, intent(in) :: n

  call err_out("fftw_init_1d")
end subroutine
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine fftw_init2_1d( me, n, is_backward )
!-------------------------------------------------------------------------------
  use, intrinsic :: iso_c_binding
  implicit none

  ! dummy variables
  class( t_fft ), intent(inout) :: me
  integer, intent(in) :: n
  logical, intent(in) :: is_backward

  call err_out("fftw_init2_1d")
end subroutine fftw_init2_1d
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine fftw_idft_1d0( me, f, n )
!-------------------------------------------------------------------------------
  use, intrinsic :: iso_c_binding
  implicit none

  class( t_fft ), intent(inout) :: me
  real(p_double), intent(inout), dimension(:, :) :: f
  integer, intent(in) :: n

  call err_out("fftw_idft_1d0")
end subroutine fftw_idft_1d0
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine fftw_idft_1d1( me, f, n )
!-------------------------------------------------------------------------------
  use, intrinsic :: iso_c_binding
  implicit none

  class( t_fft ), intent(inout) :: me
  real(p_double), intent(inout), dimension(:, :, :) :: f
  integer, dimension(2), intent(in) :: n

  call err_out("fftw_idft_1d1")
end subroutine fftw_idft_1d1
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine fftw_idft_1d2( me, f, n )
!-------------------------------------------------------------------------------
  use, intrinsic :: iso_c_binding
  implicit none

  class( t_fft ), intent(inout) :: me
  real(p_double), intent(inout), dimension(:, :, :, :) :: f
  integer, dimension(2), intent(in) :: n

  call err_out("fftw_idft_1d2")
end subroutine fftw_idft_1d2
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine fftw_r2c_1d2( me, f, n1, n2 )
!-------------------------------------------------------------------------------
  use, intrinsic :: iso_c_binding
  implicit none

  ! dummy variables
  class( t_fft ), intent(inout) :: me
  real(type_real), intent(inout), dimension(:,:,:), pointer :: f
  integer, intent(in) :: n1, n2

  call err_out("fftw_r2c_1d2")
end subroutine
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine fftw_c2r_1d2( me, f, n1, n2 )
!-------------------------------------------------------------------------------
  use, intrinsic :: iso_c_binding
  implicit none

  ! dummy variables
  class( t_fft ), intent(inout) :: me
  real(type_real), intent(inout), dimension(:,:,:), pointer :: f
  integer, intent(in) :: n1, n2

  call err_out("fftw_c2r_1d2")
end subroutine
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine fftw_r2c_1d3( me, f, n1, n2, n3 )
!-------------------------------------------------------------------------------
  use, intrinsic :: iso_c_binding
  implicit none

  ! dummy variables
  class( t_fft ), intent(inout) :: me
  real(type_real), intent(inout), dimension(:,:,:, :), pointer :: f
  integer, intent(in) :: n1, n2, n3

  call err_out("fftw_r2c_1d3")
end subroutine
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine fftw_c2r_1d3( me, f, n1, n2, n3 )
!-------------------------------------------------------------------------------
  use, intrinsic :: iso_c_binding
  implicit none

  ! dummy variables
  class( t_fft ), intent(inout) :: me
  real(type_real), intent(inout), dimension(:,:,:, :), pointer :: f
  integer, intent(in) :: n1, n2, n3

  call err_out("fftw_c2r_1d3")
end subroutine
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine fftw_cleanup_1d( me )
!-------------------------------------------------------------------------------
  use, intrinsic :: iso_c_binding
  implicit none

  class( t_fft ), intent(inout) :: me

  call err_out("fftw_cleanup_1d")
end subroutine
!-------------------------------------------------------------------------------



end module m_fft
