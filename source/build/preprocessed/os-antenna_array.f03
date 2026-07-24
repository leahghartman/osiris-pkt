# 1 "os-antenna_array.f03"
# 1 "<built-in>" 1
# 1 "<built-in>" 3
# 467 "<built-in>" 3
# 1 "<command line>" 1
# 1 "<built-in>" 2
# 1 "os-antenna_array.f03" 2

!#define DEBUG_FILE 1

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
# 5 "os-antenna_array.f03" 2
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
# 6 "os-antenna_array.f03" 2

module m_antenna_array

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
# 10 "os-antenna_array.f03" 2

use m_restart
use m_space, only : t_space
use m_emf_define, only : t_emf
use m_parameters
use m_vdf_define, only : t_vdf
use m_antenna
use m_input_file, only : t_input_file, disp_out, get_namelist

implicit none
private

real(p_double),parameter:: VER_NO_ANT_ARRAY=1.0

type t_antenna_array

  integer :: n_antenna
  class(t_antenna), pointer, dimension(:) :: ant_array => null()

contains

  procedure :: read_input => read_input_antenna_array
  procedure :: init => init_antenna_array
  procedure :: cleanup => cleanup_antenna_array
  procedure :: launch => launch_ant_array
  procedure :: restart_read => restart_read_antenna_array
  procedure :: write_checkpoint => write_checkpoint_antenna_array
  procedure :: allocate_objs => allocate_objs_antenna_array

end type t_antenna_array

interface antenna_on
  module procedure antenna_array_on
end interface

public :: t_antenna_array, antenna_on

contains

! **************************************************************
! **************************************************************

!---------------------------------------------------
subroutine allocate_objs_antenna_array( this )
!---------------------------------------------------

  implicit none

  class( t_antenna_array ), intent(inout) :: this

  allocate( t_antenna :: this%ant_array(this%n_antenna) )

end subroutine allocate_objs_antenna_array
!---------------------------------------------------

!---------------------------------------------------
subroutine restart_read_antenna_array(this, restart_handle)
!---------------------------------------------------

  implicit none

  class( t_antenna_array ), intent(inout) :: this
  type( t_restart_handle ), intent(in) :: restart_handle

  real (p_double) :: ver
  integer:: n_antenna, i, ierr

  call restart_io_read(ver, restart_handle, ierr)
  if ( ierr/=0 ) then
    write(err_buf__,*) 'error reading restart data for antenna_array object.';call err__("os-antenna_array.f03",79)
    call abort_program(p_err_rstrd)
  endif
  if(ver > VER_NO_ANT_ARRAY) then
    write(err_buf__,*) 'ANT_ARRAY: invalid restart file, ver #=',ver;call err__("os-antenna_array.f03",83)
    call abort_program(p_err_invalid)
  end if

  call restart_io_read(n_antenna, restart_handle, ierr)
  if ( ierr/=0 ) then
    write(err_buf__,*) 'error reading restart data for antenna_array object.';call err__("os-antenna_array.f03",89)
    call abort_program(p_err_rstrd)
  endif
  if ( this%n_antenna /= n_antenna ) then
    write(err_buf__,*) 'error reading restart data for antenna_array object.';call err__("os-antenna_array.f03",93)
    write(err_buf__,*) 'The number of antennas in the input file does not match the restart file';call err__("os-antenna_array.f03",94)
    call abort_program(p_err_rstrd)
  endif

  do i=1,this%n_antenna
    call this % ant_array(i) % restart_read( restart_handle )
  end do

end subroutine restart_read_antenna_array
!---------------------------------------------------

!---------------------------------------------------
subroutine write_checkpoint_antenna_array(this,restart_handle)
!---------------------------------------------------

  implicit none

  class(t_antenna_array), intent(in) :: this

  type( t_restart_handle ), intent(inout) :: restart_handle

  integer :: i, ierr

  call restart_io_write("VER_NO_ANT_ARRAY", VER_NO_ANT_ARRAY, restart_handle, ierr)
  if ( ierr/=0 ) then
    write(err_buf__,*) 'error writing restart data for antenna_array object.';call err__("os-antenna_array.f03",119)
    call abort_program(p_err_rstwrt)
  endif
  call restart_io_write("this%n_antenna", this%n_antenna, restart_handle, ierr)
  if ( ierr/=0 ) then
    write(err_buf__,*) 'error writing restart data for antenna_array object.';call err__("os-antenna_array.f03",124)
    call abort_program(p_err_rstwrt)
  endif

  do i=1,this%n_antenna
    call this % ant_array(i) % write_checkpoint( restart_handle )
  end do

end subroutine write_checkpoint_antenna_array
!---------------------------------------------------

!---------------------------------------------------
subroutine read_input_antenna_array(this, input_file)
!---------------------------------------------------

  implicit none

  class(t_antenna_array), intent(inout) :: this
  class( t_input_file ), intent(inout) :: input_file

  integer :: i
  integer :: n_antenna
  integer :: ierr

  namelist /nl_antenna_array/ n_antenna

  n_antenna=0

  ! Get namelist text from input file
  call get_namelist( input_file, "nl_antenna_array", ierr )

  if (ierr == 0) then
    read (input_file%nml_text, nml = nl_antenna_array, iostat = ierr)
    if (ierr /= 0) then
      print *, ""
      print *, "Error reading ant_array parameters "
      print *, "aborting..."
      stop
    endif
  else
    if (disp_out(input_file)) then
      if (mpi_node()==0) print *," - no antenna array specified"
    endif
    this%n_antenna=0
    return
  endif

  this%n_antenna=n_antenna
  if (n_antenna <= 0) return

  call this % allocate_objs()

  do i=1,n_antenna
    if (mpi_node()==0 .and. disp_out(input_file)) then
      print ('(A,I0,A)'),"  - antenna (",i,") configuration..."
    endif
    call this % ant_array(i) % read_input( input_file )
  end do

  end subroutine read_input_antenna_array
!---------------------------------------------------


!---------------------------------------------------
subroutine launch_ant_array( this, emf, dt, t, g_space, nx_p_min )
!---------------------------------------------------

  implicit none

  class(t_antenna_array), intent(in) :: this
  class(t_emf), intent(inout) :: emf
  real(p_double), intent(in) :: dt
  real(p_double), intent(in) :: t
  type(t_space), intent(in) :: g_space
  integer, dimension(:), intent(in) :: nx_p_min

  integer :: i

  do i=1,this%n_antenna
    if (antenna_on(this%ant_array(i))) then
        call this%ant_array(i)%antenna( emf, dt, t, g_space, nx_p_min )
    else
       print *, 'the i-th antenna is NOT on',i
    endif
  enddo

end subroutine launch_ant_array
!---------------------------------------------------

!---------------------------------------------------
logical function antenna_array_on(this)
!---------------------------------------------------

  implicit none
  class(t_antenna_array), intent(in) :: this

  antenna_array_on=(this%n_antenna > 0)

end function antenna_array_on
!---------------------------------------------------


!-------------------------------------------------------------------------------
! sets up this data structure from the given information
!-------------------------------------------------------------------------------
subroutine init_antenna_array( this, restart, restart_handle )

  implicit none
  class(t_antenna_array), intent(inout) :: this
  logical, intent(in) :: restart
  type( t_restart_handle ), intent(in) :: restart_handle

  ! no setup is required unless restarting

  if ( restart ) then
    call this % restart_read( restart_handle )
  endif

end subroutine init_antenna_array
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! cleanup dynamically allocated memory
!-------------------------------------------------------------------------------
subroutine cleanup_antenna_array( this )

  implicit none

  class(t_antenna_array), intent(inout) :: this

  if ( this%n_antenna > 0 ) then
    deallocate( this%ant_array )
  endif

end subroutine cleanup_antenna_array
!-------------------------------------------------------------------------------

end module m_antenna_array
