# 1 "os-cross.f90"
# 1 "<built-in>" 1
# 1 "<built-in>" 3
# 467 "<built-in>" 3
# 1 "<command line>" 1
# 1 "<built-in>" 2
# 1 "os-cross.f90" 2
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
# 2 "os-cross.f90" 2
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
# 3 "os-cross.f90" 2



module m_cross

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
# 9 "os-cross.f90" 2

  use m_parameters
  use m_restart

  use m_utilities

  use stringutil

  implicit none

  private

  real(p_k_fld), parameter, dimension(17) :: &
  Li_ene= (/ .9e-5, 1.0e-5, 1.3125e-5, 1.75e-5, 2e-5, 2.5e-5, &
    6.125e-5, 1e-4, 1e-3, 1e-2, 1e-1, 5e-1, &
     8e-1, 1e0, 1e1, 3e2, 1e10/)

  real(p_k_fld), parameter, dimension(17) :: &
  Li_cross=(/ tiny(1.0e0), 1.5433e-16, 1.975e-16, 2.2091e-16, 2.2439e-16, &
     2.2149e-16, 1.5648e-16,1.1973e-16, 2.0805e-17, 2.7214e-18, &
     4.0665e-19, 1.8192e-19,1.6436e-19, 1.5948e-19, 1.5465e-19, &
     1.5546e-19, 1.5546e-19 /)

  type :: t_cross

! allow access to type components only to module procedures
 private

! number of points in profile
 integer :: num_ene


! allocate later: ene,cross( num_ene)
 real(p_k_fld), dimension(:), pointer :: ene => null()
 real(p_k_fld), dimension(:), pointer :: cross => null()


  end type t_cross

  interface read_nml
 module procedure read_nml_cross
  end interface

  interface setup
 module procedure setup_cross
  end interface

  interface cleanup
    module procedure cleanup_cross
  end interface

  interface restart_write
 module procedure restart_write_cross
  end interface

  interface restart_read
 module procedure restart_read_cross
  end interface

  interface value
 module procedure value_of_cross
  end interface

! declare things that should be public
  public :: t_cross, cleanup, read_nml, setup
  public :: restart_write, restart_read
  public :: value


 contains

!! add gas_type as input so that for known gases, no need to input cross-section
!---------------------------------------------------
subroutine read_nml_cross( this, input_file, gas_type )

  use m_input_file

  implicit none

  type( t_cross ), intent(out) :: this
  class( t_input_file ), intent(inout) :: input_file
  character(len=*), intent(in) :: gas_type

! maximum number of points for the function profiles used here
! necessary because of the limitations in the use of namelists
! namelist can not include variable-size data structures
  integer, parameter :: p_num_ene_max = 100

! variables that sets the dimensionality of the profile
! and the maximum number of points
! in case of piecewise-linear profiles
  integer :: num_ene
 ! character(20) :: gas_type

! ene is in the unit of Mev, cross is in unit of cm-2
! components needed for piecewise linear profiles
! Here I changed the ene and cross section to be log format to avoid repeating log computation later
  real(p_k_fld), dimension(p_num_ene_max) :: ene, cross

   integer :: i

   namelist /nl_num_ene/ num_ene
   namelist /nl_cross/ cross, ene

  integer::ierr

! executable statements

   if((lowercase(gas_type)) == "custom") then

   num_ene = 0

      call get_namelist( input_file, "nl_num_ene", ierr )
   read (input_file%nml_text, nml = nl_num_ene)

   if ( num_ene > 0 ) then
  call alloc(this%ene, (/ num_ene /),"os-cross.f90",125)
  this%ene = 0.0_p_k_fld

  call alloc(this%cross, (/ num_ene /),"os-cross.f90",128)
  this%cross = 0.0_p_k_fld

  ene = - huge( 1.0_p_k_fld )
  cross = 0.0_p_k_fld
   endif

      call get_namelist( input_file, "nl_cross", ierr )
   read (input_file%nml_text, nml = nl_cross)

   else
  num_ene=17

  call alloc(this%ene, (/ num_ene /),"os-cross.f90",141)
  this%ene = 0.0_p_k_fld

  call alloc(this%cross, (/ num_ene /),"os-cross.f90",144)
  this%cross = 0.0_p_k_fld

  ene = - huge( 1.0_p_k_fld )
  cross = 0.0_p_k_fld
   endif

   select case ( gas_type )

  case ("h")

    print *, "   Using hydrogen as neutral gas.finish later"


  case ("he")

    print *, "   Using helium as neutral gas not finished yet."


  case ("li")

    print *, "   Using lithium as neutral gas."
    ene(1:17)= Li_ene
    cross(1:17)=Li_cross


  case ("ar")

    print *, "   Using argon as neutral gas. not finished yet"


  case ("custom")

    print *, "   Using a customized neutral gas."

  case default

    print *, "(*error*) Not a valid neutral gas -> ", gas_type
    print *, "Available gases:"
    print *, "    - H"
    print *, "    - He"
    print *, "    - Li"
    print *, "    - Ar"
    print *, "    - Custom"
    call abort_program()
   end select

   this%num_ene=num_ene
   do i=1, this%num_ene
  this%ene(i) = log10( ene(i) )
  this%cross(i) = log10( cross(i) )
   enddo

   do i=1, num_ene
  print *, "ene",this%ene(i)
  print *, "cross",this%cross(i)
   enddo

end subroutine read_nml_cross
!---------------------------------------------------


!---------------------------------------------------
subroutine setup_cross(this)
!---------------------------------------------------
! setup profile for later use
!---------------------------------------------------

  implicit none

  type( t_cross ), intent( inout ) :: this

  integer :: i

  do i =1, this%num_ene-1
     if ( this%ene(i+1) <= this%ene(i) ) then
        this%ene( i+1) = this%ene( i)
        this%cross(i+1) = this%cross(i)
     endif
  enddo

end subroutine setup_cross
!---------------------------------------------------


subroutine cleanup_cross( this )

  implicit none
  type( t_cross ), intent( inout ) :: this

  call freemem(this%ene,"os-cross.f90",234)
  call freemem(this%cross,"os-cross.f90",235)

end subroutine cleanup_cross

!---------------------------------------------------
  real(p_k_fld) function value_of_cross( this, ex )
!---------------------------------------------------
! Computes value at a point from a piecewise linear profile.
! Here it is the log value of ene and cross can be treated linearly
! the input ex is in the scale of mc^2
!---------------------------------------------------
 implicit none

! dummy variables

 type( t_cross ), intent(in) :: this
 real(p_k_fld), intent(in) :: ex

! local variables

 integer :: li
 real(p_k_fld) :: e_log, cross_log ! log10 of ex

! executable statements

! Locate exx between two values of this%ene, and obtain
! value of density by linear interpolation.
! If exx lies outside range, density=limit-density.

! exx position at which cross desired
! this%ene array of energy at which functions are specified
! this%cross array of corresponding values

    value_of_cross =10 ** this%cross(1)
    e_log = log10(ex*.511)
    if ( e_log <= this%ene(1) ) return
    do li=2,this%num_ene
    if ( e_log <= this%ene(li) ) then
    cross_log = this%cross(li-1) + ( this%cross(li) - this%cross(li-1) ) / &
                     (this%ene( li) - this%ene( li-1) ) * ( e_log - this%ene( li-1) )

    value_of_cross = 10**cross_log

    return
    endif
    enddo

    value_of_cross = 10**this%cross(this%num_ene)

 end function value_of_cross
!---------------------------------------------------


!---------------------------------------------------
  subroutine restart_write_cross( this, restart_handle )
!---------------------------------------------------
! write object information into a restart file
!---------------------------------------------------

! <use-statements for this subroutine>

 implicit none

! dummy variables

 type( t_cross ), intent(in) :: this
    type( t_restart_handle ), intent(inout) :: restart_handle

! local variables

 logical :: if_x, if_cross
 integer :: ierr

! executable statements

 if_x = associated( this%ene )
 if_cross = associated( this%cross )

    call restart_io_write("this%num_ene", this%num_ene, restart_handle, ierr)
    if ( ierr/=0 ) then
      write(err_buf__,*) 'error writing restart data for cross-section object.';call err__("os-cross.f90",315)
      call abort_program( p_err_rstwrt)
    endif
    call restart_io_write("if_x", if_x, restart_handle, ierr)
    if ( ierr/=0 ) then
      write(err_buf__,*) 'error writing restart data for cross-section object.';call err__("os-cross.f90",320)
      call abort_program( p_err_rstwrt)
    endif
    call restart_io_write("if_cross", if_cross, restart_handle, ierr)
    if ( ierr/=0 ) then
      write(err_buf__,*) 'error writing restart data for cross-section object.';call err__("os-cross.f90",325)
      call abort_program( p_err_rstwrt)
    endif

 if ( if_x ) then
      call restart_io_write("size(this%ene)", size(this%ene), restart_handle, ierr)
      if ( ierr/=0 ) then
        write(err_buf__,*) 'error writing restart data for cross-section object.';call err__("os-cross.f90",332)
        call abort_program( p_err_rstwrt)
      endif
      call restart_io_write("this%ene", this%ene, restart_handle, ierr)
      if ( ierr/=0 ) then
        write(err_buf__,*) 'error writing restart data for cross-section object.';call err__("os-cross.f90",337)
        call abort_program( p_err_rstwrt)
      endif
 endif
 if ( if_cross ) then
      call restart_io_write("size(this%cross)", size(this%cross), restart_handle, ierr)
      if ( ierr/=0 ) then
        write(err_buf__,*) 'error writing restart data for cross-section object.';call err__("os-cross.f90",344)
        call abort_program( p_err_rstwrt)
      endif
      call restart_io_write("this%cross", this%cross, restart_handle, ierr)
      if ( ierr/=0 ) then
        write(err_buf__,*) 'error writing restart data for cross-section object.';call err__("os-cross.f90",349)
        call abort_program( p_err_rstwrt)
      endif
 endif

  end subroutine restart_write_cross
!---------------------------------------------------

!---------------------------------------------------
  subroutine restart_read_cross( this, restart_handle )
!---------------------------------------------------
! read object information from a restart file
!---------------------------------------------------

! <use-statements for this subroutine>

 implicit none

! dummy variables

 type( t_cross ), intent(inout) :: this
    type( t_restart_handle ), intent(in) :: restart_handle

! local variables

 integer :: s1, ierr
 logical :: if_x, if_cross

! executable statements

 call cleanup(this)

 call restart_io_read(this%num_ene, restart_handle, ierr)
 if ( ierr/=0 ) then
      write(err_buf__,*) 'error reading restart data for cross-section object.';call err__("os-cross.f90",383)
      call abort_program( p_err_rstwrt)
    endif
 call restart_io_read(if_x, restart_handle, ierr)
 if ( ierr/=0 ) then
      write(err_buf__,*) 'error reading restart data for cross-section object.';call err__("os-cross.f90",388)
      call abort_program( p_err_rstwrt)
    endif
 call restart_io_read(if_cross, restart_handle, ierr)
 if ( ierr/=0 ) then
      write(err_buf__,*) 'error reading restart data for cross-section object.';call err__("os-cross.f90",393)
      call abort_program( p_err_rstwrt)
    endif

 if ( if_x ) then
       call restart_io_read(s1, restart_handle, ierr)
       if ( ierr/=0 ) then
         write(err_buf__,*) 'error reading restart data for cross-section object.';call err__("os-cross.f90",400)
         call abort_program( p_err_rstwrt)
       endif

    call alloc(this%ene, (/s1/),"os-cross.f90",404)
       call restart_io_read(this%ene, restart_handle, ierr)
       if ( ierr/=0 ) then
         write(err_buf__,*) 'error reading restart data for cross-section object.';call err__("os-cross.f90",407)
         call abort_program( p_err_rstwrt)
       endif
 endif
 if ( if_cross ) then
       call restart_io_read(s1, restart_handle, ierr)
       if ( ierr/=0 ) then
         write(err_buf__,*) 'error reading restart data for cross-section object.';call err__("os-cross.f90",414)
         call abort_program( p_err_rstwrt)
       endif

    call alloc(this%cross, (/s1/),"os-cross.f90",418)
       call restart_io_read(this%cross, restart_handle, ierr)
       if ( ierr/=0 ) then
         write(err_buf__,*) 'error reading restart data for cross-section object.';call err__("os-cross.f90",421)
         call abort_program( p_err_rstwrt)
       endif
 endif

  end subroutine restart_read_cross
!---------------------------------------------------


end module m_cross
