# 1 "os-logprof.f90"
# 1 "<built-in>" 1
# 1 "<built-in>" 3
# 467 "<built-in>" 3
# 1 "<command line>" 1
# 1 "<built-in>" 2
# 1 "os-logprof.f90" 2
!-------------------------------------------------------------------------------------------------
! m_logprof module
!
! This module handles logging and profiling of events. It can optionally be used with the MPE
! library and/or with the PAPI library, just set the compiler macros __USE_MPE__ and or
! __USE_PAPI__ and link with the appropriate libs.
!
!-------------------------------------------------------------------------------------------------

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
# 11 "os-logprof.f90" 2
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
# 12 "os-logprof.f90" 2

module m_logprof

use m_parameters
use m_system
use m_diagfile, only : p_time_length

implicit none

! restrict access to things explicitly declared public
private

integer, parameter :: p_max_events = 64
integer, parameter :: p_max_event_name = 128

!---------------------------------------------------------------------------------------------------
! PAPI event data
!---------------------------------------------------------------------------------------------------
# 52 "os-logprof.f90"
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! MPE event data
!---------------------------------------------------------------------------------------------------
# 69 "os-logprof.f90"
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! t_event class definition
!---------------------------------------------------------------------------------------------------
type :: t_event

  character(len = p_max_event_name) :: name
# 86 "os-logprof.f90"
  integer( p_int64 ) :: timer_begin
  real( p_double ) :: total_time


end type t_event
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! Module variables
!---------------------------------------------------------------------------------------------------

! Total number of events created
integer, save :: num_events = 0

! Events array
type( t_event ), dimension(p_max_events), save :: events

!---------------------------------------------------------------------------------------------------

interface create_event
  module procedure create_event
end interface

interface begin_event
  module procedure begin_event
end interface

interface end_event
  module procedure end_event
end interface

interface total_event_time
  module procedure total_event_time
end interface

interface list_total_event_times
  module procedure list_total_event_times
end interface



public :: create_event, begin_event, end_event
public :: total_event_time, list_total_event_times
# 141 "os-logprof.f90"
 contains

!---------------------------------------------------------------------------------------------------
! Create an event for profiling
!---------------------------------------------------------------------------------------------------
function create_event( name )

  implicit none

  integer :: create_event
  character(len=*), intent(in) :: name





  num_events = num_events+1

  if ( num_events > p_max_events ) then
    if (mpi_node() == 0 ) then
      write(0, "('(*error*) Unable to create event ', A,', too many events.')" ) trim(name)
      write(0, "('(*error*) The maximum number of events is currently set at ',I0)" ) p_max_events
    endif
    call abort_program(p_err_invalid)
  endif

  create_event = num_events

  events(num_events)%name = name
# 183 "os-logprof.f90"
  events(num_events)%total_time = 0.0D0


end function create_event
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
! Start an event
!---------------------------------------------------------------------------------------------------
subroutine begin_event( event_num )

  implicit none

  integer, intent(in) :: event_num







  if( if_profiling() ) then






   ! Read the current timer ticks
   events(event_num)%timer_begin = timer_ticks()

  endif

end subroutine begin_event
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! Finish an event
!---------------------------------------------------------------------------------------------------
subroutine end_event( event_num )

  implicit none

  integer, intent(in) :: event_num





  integer( p_int64 ) :: timer_end

  if( if_profiling() ) then
! Terminate the event



  timer_end = timer_ticks()






! Accumulate time



  events(event_num)%total_time = events(event_num)%total_time + &
                                 timer_interval_seconds(events(event_num)%timer_begin, timer_end)

  endif

end subroutine end_event
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! Return the total event time for a given event
!---------------------------------------------------------------------------------------------------
function total_event_time( event_num )

  implicit none

  real( p_double ) :: total_event_time

  integer, intent(in) :: event_num




  total_event_time = events(event_num)%total_time


end function total_event_time
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! List the avg, min and max times for all events
!---------------------------------------------------------------------------------------------------
subroutine list_total_event_times( n, comm, label )

  !use mpi

  implicit none

  ! dummy variables
  integer, intent(in) :: n, comm
  character(len=*), optional, intent(in) :: label

  ! local variables
  integer :: i, rank, size, ierr
  real( p_double ) , dimension( num_events ) :: time, time_min, time_max, time_avg
  character(len=128) :: fname
  character(len=2) :: tim_len
# 306 "os-logprof.f90"
  fname = trim(path_time)//p_dir_sep//'timings-'
  if ( present( label ) ) then
    fname = trim(fname) // trim(label)
  else
    write( tim_len, '(i2)' ) p_time_length
    write( fname, '(A,i'//trim(tim_len)//'.'//trim(tim_len)//')' ) trim(fname), n
  endif

  if ( comm == MPI_COMM_NULL ) then


 rank = 0

 ! Create output directory if necessary
 call mkdir( trim(path_time), ierr )
 if ( ierr /= 0 ) then
   write( 0, * ) 'Error creating TIMINGS ', strerror( ierr )
   return
 endif

    open ( unit = file_id_prof , file = trim(fname), form = 'formatted' )

 write( file_id_prof, '(A39,1X,A19)' ) &
        ' Event', 'Total [s]'
    ! Serial run
 write( file_id_prof, '(A)' ) &
     '-----------------------------------------------------------'






 do i = 1, num_events
   write( file_id_prof, '(A39,1X,E19.3)' ) trim( events(i)%name ), events(i)%total_time
 enddo


    close(file_id_prof)

  else

    call mpi_comm_rank( comm, rank, ierr )
    call mpi_comm_size( comm, size, ierr )

    ! Parallel run





    do i = 1, num_events
      time(i) = events(i)%total_time
    enddo


    ! Get min, max and avg times
    call mpi_reduce( time, time_min, num_events, MPI_DOUBLE_PRECISION, MPI_MIN, 0, comm, ierr )
    call mpi_reduce( time, time_max, num_events, MPI_DOUBLE_PRECISION, MPI_MAX, 0, comm, ierr )

    call mpi_reduce( time, time_avg, num_events, MPI_DOUBLE_PRECISION, MPI_SUM, 0, comm, ierr )

    if ( rank == 0 ) then

       do i = 1, num_events
         time_avg(i) = time_avg(i) / size
       enddo

    ! Create output directory if necessary
    call mkdir( trim(path_time), ierr )
    if ( ierr /= 0 ) then
   write( 0, * ) 'Error creating TIMINGS ', strerror( ierr )
   return
    endif

    open ( unit = file_id_prof , file = trim(fname), form = 'formatted' )

       write( file_id_prof, * ) ' Iterations = ', n
       write( file_id_prof, * ) ' '
       write( file_id_prof, '(A39,3(1X,A19))' ) &
                    'Event', 'Avg [s]', 'Min [s]', 'Max [s]'
       write( file_id_prof, '(A)' ) &
                    '------------------------------------------------------------------------------'

       do i = 1, num_events
         write( file_id_prof, '(A39,3(1X,E19.6))' ) trim( events(i)%name ), &
             time_avg(i), time_min(i), time_max(i)
       enddo

       close(file_id_prof)

    endif

  endif
# 437 "os-logprof.f90"
end subroutine list_total_event_times
!---------------------------------------------------------------------------------------------------
# 670 "os-logprof.f90"
end module m_logprof
