# 1 "os-timeavg.f90"
# 1 "<built-in>" 1
# 1 "<built-in>" 3
# 467 "<built-in>" 3
# 1 "<command line>" 1
# 1 "<built-in>" 2
# 1 "os-timeavg.f90" 2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! time average class
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This is a support class for doing time averages of grid
! quantities. This class can take care of two types of grid
! data: VDF objects, and standard F90 single precision arrays

! (1,2 and 3D). The VDF objects exist on the parallel universe
! and the grids are node local. When doing a diagnostic on grid
! quantities they are summed from all the nodes and written by
! node 0.
!
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
# 17 "os-timeavg.f90" 2
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
# 18 "os-timeavg.f90" 2

module m_time_avg

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
# 22 "os-timeavg.f90" 2

use m_parameters
use m_utilities

use m_node_conf
use m_grid_define
use m_space

use m_utilities

use m_diagnostic_utilities

use m_vdf_define
use m_vdf_math
use m_vdf_report
use m_vdf_average
use m_vdf_comm
use m_vdf

implicit none

! restrict access to things explicitly declared public
private

integer, parameter :: p_tavg_none = -1
integer, parameter :: p_tavg_vdf = 0
integer, parameter :: p_tavg_grid = 1

type :: t_time_avg

  ! type of time average object
  integer :: tavg_type = p_tavg_none

  ! vdf data for time averaged diagnostics
  type( t_vdf ) :: vdf_data

  ! grid data for time averaged diagnostics
  real(p_diag_prec), dimension(:), pointer :: buffer => null()
  integer :: rank
  integer, dimension(3) :: n

  ! number of datsets already added for time averaging
  ! -1 denotes memory not initialized
  integer :: n_tavg = -1

end type t_time_avg

interface setup
  module procedure setup_grid_time_avg
end interface setup

interface cleanup
  module procedure cleanup_time_avg
end interface cleanup

interface add
  module procedure add_vdf_time_avg
end interface add

interface get_buffer
  module procedure get_buffer_time_avg
end interface get_buffer

interface reset
  module procedure reset_time_avg
end interface reset

interface report
! module procedure report_time_avg_vdf
  module procedure report_time_avg_grid
end interface report

interface if_time_avg_add
  module procedure if_time_avg_add_single
  module procedure if_time_avg_add_multi
end interface if_time_avg_add

interface initialized
  module procedure initialized_time_avg
end interface initialized

interface reshape_obj
  module procedure reshape_time_avg
end interface

interface do_average
  module procedure do_average_time_avg
end interface

public :: t_time_avg, setup, cleanup, reshape_obj
public :: add, reset, do_average, get_buffer
public :: report
! public :: report_ave
public :: if_time_avg_add, initialized

 contains



!---------------------------------------------------
  subroutine setup_grid_time_avg( this, rank, n )
!---------------------------------------------------
! sets up this data structure from the given information
!---------------------------------------------------

! <use-statements for this subroutine>

    implicit none

! dummy variables

    type( t_time_avg ), intent( inout ) :: this
    integer, intent(in) :: rank
    integer, dimension(:), intent(in) :: n

! executable statements

    ! check if object is not already setup
    call asrt__(this%tavg_type == p_tavg_none,"this%tavg_type == p_tavg_none","os-timeavg.f90",140)

    call asrt__(rank >= 1 .and. rank <= 3,"rank >= 1 .and. rank <= 3","os-timeavg.f90",142)

    ! set the type
    this%tavg_type = p_tavg_grid

    this%rank = rank
    this%n(1:rank) = n(1:rank)

    ! allocate the grid
    call alloc(this%buffer, (/ product(n(1:rank)) /),"os-timeavg.f90",151)

    ! reset the data structure
    call reset(this)

  end subroutine setup_grid_time_avg
!---------------------------------------------------

!---------------------------------------------------
subroutine cleanup_time_avg( this )
!---------------------------------------------------
! cleanup this object
!---------------------------------------------------

   implicit none

   type( t_time_avg ), intent( inout ) :: this

   select case (this%tavg_type)
     case ( p_tavg_vdf )
      call this%vdf_data%cleanup()

     case ( p_tavg_grid )
      call freemem(this%buffer,"os-timeavg.f90",174)

     case default ! clearing a non-initialized object is ok
      continue
   end select

   ! clear the type
   this%tavg_type = p_tavg_none

   ! clear the number of datasets added
   this%n_tavg = -1

end subroutine cleanup_time_avg
!---------------------------------------------------

!---------------------------------------------------
  subroutine reset_time_avg( this )
!---------------------------------------------------
! reset the average
!---------------------------------------------------

    implicit none

! dummy variables

    type( t_time_avg ), intent( inout ) :: this

! local variables

! executable statements

    ! Check if we have a valid object
    call asrt__(this%tavg_type /= p_tavg_none,"this%tavg_type /= p_tavg_none","os-timeavg.f90",206)

    ! reset data
    select case ( this%tavg_type )
      case ( p_tavg_vdf )

         this%vdf_data = 0.0

      case ( p_tavg_grid )
         this%buffer = 0

    end select

    ! clear the number of datasets added
    this%n_tavg = 0

  end subroutine reset_time_avg
!---------------------------------------------------

!---------------------------------------------------
  subroutine add_vdf_time_avg( this, vdf_a, fcomp )
!---------------------------------------------------
! adds a time frame for the time averaged

!---------------------------------------------------

    implicit none

! dummy variables

    type( t_time_avg ), intent(inout) :: this
    type( t_vdf ), intent(in) :: vdf_a
    integer, intent(in) :: fcomp

! local variables

! executable statements

     ! Check if we have a valid object
     call asrt__(this%tavg_type == p_tavg_vdf,"this%tavg_type == p_tavg_vdf","os-timeavg.f90",245)

     this%n_tavg = this%n_tavg + 1
     call add( this%vdf_data, 1, vdf_a, fcomp )

  end subroutine add_vdf_time_avg
!---------------------------------------------------

!---------------------------------------------------
! Gets a pointer to the local data buffer to add a
! time frame externally. This also advances n_tavg
!---------------------------------------------------
subroutine get_buffer_time_avg( this, buffer )

  implicit none

  type( t_time_avg ), intent(inout) :: this
  real(p_diag_prec), dimension(:), pointer :: buffer

  call asrt__(this%tavg_type == p_tavg_grid,"this%tavg_type == p_tavg_grid","os-timeavg.f90",264)

  this%n_tavg = this%n_tavg + 1
  buffer => this%buffer

end subroutine get_buffer_time_avg
!---------------------------------------------------

!---------------------------------------------------
subroutine do_average_time_avg( this )
!---------------------------------------------------
! Divides the stored values by n_tavg
!---------------------------------------------------

   implicit none

   type( t_time_avg ), intent(inout) :: this

    select case ( this%tavg_type )
      case ( p_tavg_vdf )
         call mult( this%vdf_data, 1.0_p_k_fld / this%n_tavg )

      case ( p_tavg_grid )
         this%buffer = this%buffer / this%n_tavg

    end select

end subroutine do_average_time_avg
!---------------------------------------------------

!---------------------------------------------------
subroutine report_time_avg_grid( this, no_co, diagFile, name )
!---------------------------------------------------
! save diagnostic information regarding this tavg_diag
!---------------------------------------------------

  implicit none

! dummy variables and functions

  type( t_time_avg ), intent(inout) :: this ! tavg_diag object

  class( t_node_conf ), intent(in) :: no_co ! node configuration of the object

  class(t_diag_file), allocatable :: diagFile
  character(len=*), intent(in) :: name

  real(p_diag_prec), dimension(:,:), pointer :: buffer2
  real(p_diag_prec), dimension(:,:,:), pointer :: buffer3
  ! local variables
  integer :: bsize, ierr

! executable statements

  ! Check if we have a valid object
  call asrt__(this%tavg_type == p_tavg_grid,"this%tavg_type == p_tavg_grid","os-timeavg.f90",319)

  if ( no_num(no_co) > 1) then
    bsize = size( this%buffer )

    if ( root(no_co) ) then
       call MPI_REDUCE(MPI_IN_PLACE, this%buffer, bsize, mpi_real_type( p_diag_prec ), MPI_SUM, 0, comm(no_co), ierr)
    else
       call MPI_REDUCE(this%buffer, 0, bsize, mpi_real_type( p_diag_prec ), MPI_SUM, 0, comm(no_co), ierr)
    endif
  endif

  if ( root(no_co) ) then
    ! create the file
    call diagFile % open( p_diag_create )

    ! Write the data
    select case ( diagFile % grid % ndims )
      case ( 1 )
        call diagFile % add_dataset( name, this%buffer )

      case ( 2 )
        buffer2(1:this%n(1), 1:this%n(2)) => this%buffer
        call diagFile % add_dataset( name, buffer2 )

      case ( 3 )
        buffer3(1:this%n(1), 1:this%n(2), 1:this%n(3)) => this%buffer
        call diagFile % add_dataset( name, buffer3 )
    end select
    !call diagFile % add_dataset( name, reshape( this%buffer, this%n ) )

    ! Close the file
    call diagFile % close( )

    deallocate( diagFile )
  endif

end subroutine report_time_avg_grid
!---------------------------------------------------

!---------------------------------------------------
  function if_time_avg_add_single( n, ndump, ntavg )
!---------------------------------------------------
! returns true if it is time to add a time
! frame to the time average
!---------------------------------------------------

    implicit none

! dummy variables
    integer, intent(in) :: n, ndump, ntavg
    logical :: if_time_avg_add_single

! local variables

! executable statements

    if_time_avg_add_single = (ndump*ntavg > 0)
    if (if_time_avg_add_single) then
         if_time_avg_add_single = ( (mod(n, ndump) > ndump - ntavg) .or. &
                                     (mod(n, ndump) == 0) )
    endif

  end function if_time_avg_add_single
!---------------------------------------------------

!---------------------------------------------------
  function if_time_avg_add_multi( n, ndump, ntavg )
!---------------------------------------------------
! returns true if it is time to add a time
! frame to the time average
!---------------------------------------------------

    implicit none

    ! dummy variables
    integer, intent(in) :: n, ntavg
    integer, dimension(:), intent(in) :: ndump
    logical :: if_time_avg_add_multi

    integer :: i

    if_time_avg_add_multi = .false.
    do i = 1, size( ndump )
      if_time_avg_add_multi = (ndump(i)*ntavg > 0)
      if (if_time_avg_add_multi) then
           if_time_avg_add_multi = ( (mod(n, ndump(i)) > ndump(i) - ntavg) .or. &
                                      (mod(n, ndump(i)) == 0) )
      endif
      if ( if_time_avg_add_multi ) exit

    enddo

  end function if_time_avg_add_multi
!---------------------------------------------------

!---------------------------------------------------
  function initialized_time_avg( this )
!---------------------------------------------------
! returns true if it is time to add a time
! frame to the time average
!---------------------------------------------------

    implicit none

! dummy variables
    logical :: initialized_time_avg

    type( t_time_avg ), intent(in) :: this

! local variables

! executable statements

    initialized_time_avg = (this%tavg_type /= p_tavg_none)

  end function initialized_time_avg
!---------------------------------------------------

!---------------------------------------------------
subroutine reshape_time_avg( this, old_lb, new_lb, no_co, send_msg, recv_msg )
!---------------------------------------------------
   implicit none

   ! dummy vars
   type(t_time_avg), intent(inout) :: this
   class(t_grid), intent(in) :: old_lb, new_lb
   class( t_node_conf ), intent(in) :: no_co
   type(t_vdf_msg), dimension(2), intent(inout) :: send_msg, recv_msg

   ! only vdf time averages need reshaping
   if ( this%tavg_type == p_tavg_vdf ) then
     call reshape_copy( this%vdf_data, old_lb, new_lb, no_co, send_msg, recv_msg )
   endif

end subroutine reshape_time_avg
!---------------------------------------------------

end module m_time_avg
