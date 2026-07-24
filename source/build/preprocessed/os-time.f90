# 1 "os-time.f90"
# 1 "<built-in>" 1
# 1 "<built-in>" 3
# 467 "<built-in>" 3
# 1 "<command line>" 1
# 1 "<built-in>" 2
# 1 "os-time.f90" 2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! time class
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
# 6 "os-time.f90" 2
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
# 7 "os-time.f90" 2


module m_time

use m_restart
use m_utilities
use m_parameters

implicit none

private

! string to id restart data
character(len=*), parameter :: p_time_rst_id = "time rst data - 0x0002"


type :: t_time

  private

  real(p_double) :: tmin, tmax, t

end type t_time


interface read_nml
  module procedure read_nml_time
end interface

interface setup
  module procedure setup_time
end interface

interface restart_write
  module procedure restart_write_time
end interface

interface restart_read
  module procedure restart_read_time
end interface

interface update
  module procedure update_time
end interface

interface test
  module procedure test_time
end interface

interface tmin
  module procedure tmin
end interface

interface tmax
  module procedure tmax
end interface

interface t
  module procedure t
end interface

! declare things that should be public
public :: t_time, read_nml, setup
public :: restart_write
public :: update, test
public :: tmin, tmax, t


contains


!---------------------------------------------------
!---------------------------------------------------
subroutine read_nml_time( this, input_file )

  use m_input_file

  implicit none

  type( t_time ), intent(out) :: this
  class( t_input_file ), intent(inout) :: input_file

  real(p_double) :: tmin, tmax

  namelist /nl_time/ tmin, tmax

  integer :: ierr

! executable statements

  tmin = 0.0_p_double
  tmax = 0.0_p_double

  ! Get namelist text from input file
  call get_namelist( input_file, "nl_time", ierr )

  if (ierr /= 0) then
 if (ierr < 0) then
   print *, "Error reading time parameters"
 else
   print *, "Error: time parameters missing"
 endif
 print *, "aborting..."
 stop
  endif


  read (input_file%nml_text, nml = nl_time, iostat = ierr)
  if (ierr /= 0) then
 print *, "Error reading time parameters"
 print *, "aborting..."
 stop
  endif

  this%tmin = tmin
  this%tmax = tmax

end subroutine read_nml_time
!---------------------------------------------------

!---------------------------------------------------
subroutine setup_time( this, initial, dt, restart, restart_handle )
!---------------------------------------------------
! sets up this data structure from the given information
!---------------------------------------------------

! <use-statements for this subroutine>

  implicit none

! dummy variables

  type( t_time ), intent( inout ) :: this
  real(p_double), intent( in ) :: initial, dt
  logical, intent(in) :: restart
  type( t_restart_handle ), intent(in) :: restart_handle

! local variables

  ! set initial time
  if ( restart ) then
 call restart_read( this, restart_handle )
  else
 this%t = initial
  endif

! Check to assure that this calculation will run for at least one cycle

  if ( this%tmax < ( this%t + 2 * dt ) ) then
  this%tmax = this%t + 2 * dt
  endif

end subroutine setup_time
!---------------------------------------------------


!---------------------------------------------------
subroutine restart_write_time( this, restart_handle )
!---------------------------------------------------
! write object information into a restart file
!---------------------------------------------------

! <use-statements for this subroutine>

  implicit none

! dummy variables

  type( t_time ), intent(in) :: this
  type( t_restart_handle ), intent(inout) :: restart_handle

! local variables

  integer :: ierr

! executable statements

  call restart_io_write("p_time_rst_id", p_time_rst_id, restart_handle, ierr)
  if ( ierr/=0 ) then
 write(err_buf__,*) "error writing restart data for time object.";call err__("os-time.f90",186)
 call abort_program(p_err_rstwrt)
  endif

  call restart_io_write("this%t", this%t, restart_handle, ierr)
  if ( ierr/=0 ) then
 write(err_buf__,*) 'error writing restart data for time object.';call err__("os-time.f90",192)
 call abort_program(p_err_rstwrt)
  endif
  call restart_io_write("this%tmin", this%tmin, restart_handle, ierr)
  if ( ierr/=0 ) then
 write(err_buf__,*) 'error writing restart data for time object.';call err__("os-time.f90",197)
 call abort_program(p_err_rstwrt)
  endif


end subroutine restart_write_time
!---------------------------------------------------


!---------------------------------------------------
subroutine restart_read_time( this, restart_handle )
!---------------------------------------------------
! read object information from a restart file
!---------------------------------------------------

! <use-statements for this subroutine>

  implicit none

! dummy variables

  type( t_time ), intent(inout) :: this
  type( t_restart_handle ), intent(in) :: restart_handle

! local variables

  integer :: ierr
  character(len=len(p_time_rst_id)) :: rst_id


! executable statements

! check if restart file is compatible
  call restart_io_read(rst_id, restart_handle, ierr)
  if ( ierr/=0 ) then
  write(err_buf__,*) 'error reading restart data for time object.';call err__("os-time.f90",232)
  call abort_program(p_err_rstrd)
  endif

  if ( rst_id /= p_time_rst_id) then
    write(err_buf__,*) 'Corrupted restart file, or restart file';call err__("os-time.f90",237)
    write(err_buf__,*) 'from incompatible binary (time)';call err__("os-time.f90",238)
 call abort_program(p_err_rstrd)
  endif

  call restart_io_read(this%t, restart_handle, ierr)
  if ( ierr/=0 ) then
 write(err_buf__,*) 'error reading restart data for time object.';call err__("os-time.f90",244)
 call abort_program(p_err_rstrd)
  endif
  call restart_io_read(this%tmin, restart_handle, ierr)
  if ( ierr/=0 ) then
 write(err_buf__,*) 'error reading restart data for time object.';call err__("os-time.f90",249)
 call abort_program(p_err_rstrd)
  endif

end subroutine restart_read_time
!---------------------------------------------------


!---------------------------------------------------
subroutine update_time( this, dt, n )
!---------------------------------------------------
! update time using new iteration value
!---------------------------------------------------

implicit none

! dummy variables

type( t_time ), intent( inout ) :: this

real(p_double), intent( in ) :: dt
integer, intent(in) :: n


this%t = this%tmin + n*dt


end subroutine update_time
!---------------------------------------------------



!---------------------------------------------------
function test_time( this )
!---------------------------------------------------
! use "this" to determine new value of "test"
!---------------------------------------------------

! <use-statements for this subroutine>

  implicit none

! dummy variables

  logical :: test_time

  type( t_time ), intent( in ) :: this

! local variables - none

! executable statements

! get value of result
  test_time = ( this%t < this%tmax )

end function test_time
!---------------------------------------------------


!---------------------------------------------------
function tmin( this )
!---------------------------------------------------
! gives the minimum (starting) time
!---------------------------------------------------

  implicit none

! dummy variables

  real(p_double) :: tmin

  type( t_time ), intent(in) :: this

! local variables - none

! executable statements

  tmin = this%tmin

end function tmin
!---------------------------------------------------


!---------------------------------------------------
function tmax( this )
!---------------------------------------------------
! gives the maximum (final) time
!---------------------------------------------------

  implicit none

! dummy variables

  real(p_double) :: tmax

  type( t_time ), intent(in) :: this

! local variables - none

! executable statements

  tmax = this%tmax

end function tmax
!---------------------------------------------------


!---------------------------------------------------
function t( this )
!---------------------------------------------------
! gives the current time
!---------------------------------------------------

  implicit none

! dummy variables

  real(p_double) :: t

  type( t_time ), intent(in) :: this

! local variables - none

! executable statements

  t = this%t

end function t
!---------------------------------------------------


end module m_time
