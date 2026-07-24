# 1 "cyl_modes/os-antenna-cyl-modes.f03"
# 1 "<built-in>" 1
# 1 "<built-in>" 3
# 467 "<built-in>" 3
# 1 "<command line>" 1
# 1 "<built-in>" 2
# 1 "cyl_modes/os-antenna-cyl-modes.f03" 2
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
# 4 "cyl_modes/os-antenna-cyl-modes.f03" 2
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
# 5 "cyl_modes/os-antenna-cyl-modes.f03" 2

module m_antenna_cyl_modes

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
# 9 "cyl_modes/os-antenna-cyl-modes.f03" 2

use m_space, only : t_space, xmin
use m_emf_define, only : t_emf
use m_emf_bound, only : is_open
use m_emf_cyl_modes, only : t_emf_cyl_modes
use m_parameters
use m_math, only : pi_180
use m_vdf_define, only : t_vdf
use m_antenna
use m_antenna_array, only : t_antenna_array
use m_input_file, only : t_input_file

implicit none

private

type, extends( t_antenna ) :: t_antenna_cyl_modes

  ! no class-specific member data

contains

  procedure :: read_input => read_input_antenna_cyl_modes
  procedure :: antenna => antenna_2d_cyl_modes

end type t_antenna_cyl_modes

! Define here instead of in separate file since only one method is overridden
type, extends( t_antenna_array ) :: t_antenna_array_cyl_modes

  ! no class-specific member data

contains

  procedure :: allocate_objs => allocate_objs_antenna_array_cyl_modes

end type t_antenna_array_cyl_modes

interface l1
  module procedure l1_1
end interface

interface l2
  module procedure l2_1
end interface

interface l3
  module procedure l3_1
end interface

interface l4
  module procedure l4_1
end interface

interface l5
  module procedure l5_1
end interface

interface l6
  module procedure l6_1
end interface

public :: t_antenna_cyl_modes, t_antenna_array_cyl_modes

contains

!---------------------------------------------------
subroutine allocate_objs_antenna_array_cyl_modes( this )
!---------------------------------------------------

  implicit none

  class( t_antenna_array_cyl_modes ), intent(inout) :: this

  allocate( t_antenna_cyl_modes :: this%ant_array(this%n_antenna) )

end subroutine allocate_objs_antenna_array_cyl_modes
!---------------------------------------------------

!----------------------------------------------------
subroutine read_input_antenna_cyl_modes( this, input_file )
!----------------------------------------------------
! read input values for antenna
!----------------------------------------------------

  implicit none
  class( t_antenna_cyl_modes ), intent(inout) :: this
  class( t_input_file ), intent(inout) :: input_file

  real (p_k_fld) :: fac1

  call this % t_antenna % read_input( input_file )

  ! For now asymmetric beams are effectively ignored, only rad_x, etc. processed
  ! TODO: put in asymmetric beams

  ! Handle quasi-3d-specific modifications
  if (this%focus.ne.0.0) then !! renormalize some quantities

    fac1=sqrt(1.0+(this%focus*this%focus)/(this%z_r(1)*this%z_r(1)))

    ! this%phase2(1)=-0.5*atan(this%focus/this%z_r(1))
    this%phase2(1)=-atan(this%focus/this%z_r(1))

    ! We want this%a0 = a0/fac1, but a factor of 1/sqrt(fac1) is
    ! applied in read_input_antenna
    this%a0=this%a0/sqrt(fac1)

  endif

  if ( this%direction > 1 ) then
    if ( mpi_node() == 0 ) then
      write(0,*) ""
      write(0,*) "   Error reading antenna parameters"
      write(0,*) "   antenna only supported in direction 1"
      write(0,*) "   when running with cylindrical modes"
      write(0,*) "   aborting..."
    endif
    stop
  endif

end subroutine read_input_antenna_cyl_modes
!----------------------------------------------------


!-----------------------------------------------------------------------------------------
! Antenna module for hybrid cylindrical modes
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine antenna_2d_cyl_modes(this, emf, dt, t, g_space, nx_p_min )
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------

  implicit none
  class(t_antenna_cyl_modes), intent(in) :: this
  class(t_emf), intent(inout) :: emf
  real (p_double),intent(in):: t, dt
  type (t_space),intent(in)::g_space
  integer, dimension(:), intent(in) :: nx_p_min

  !! local variables
  integer ind1
  real (p_k_fld),dimension(p_x_dim):: delta_x,lower_bound
  integer ,dimension(2,2):: array_bound
  real (p_k_fld) :: a_t,a_tx,phase

  real (p_k_fld),dimension(2)::r,coord
  real (p_k_fld) :: local_time,true_time,local_omega
  integer :: ind_wall
  integer :: i
  type (t_vdf), pointer :: efield_re, efield_im

  ! if not on an injecting wall return silently
  if ( .not. is_open(emf%bnd_con, this%direction, this%side) ) return

  select type (emf)
  class is (t_emf_cyl_modes)
    efield_re => emf%b_cyl_m%pf_re(1)
    efield_im => emf%b_cyl_m%pf_im(1)
  end select

  delta_x = real( efield_re%dx(), p_k_fld)
  do i = 1,p_x_dim
    lower_bound(i) = real(xmin(g_space,i) + (nx_p_min(i) - 1)*efield_re%dx( i ), p_k_fld)
  end do

  array_bound=efield_re%gc_num()
  array_bound(p_lower,1)=array_bound(p_lower,1)+lbound(efield_re%f2,dim=2)
  array_bound(p_lower,2)=array_bound(p_lower,2)+lbound(efield_re%f2,dim=3)

  array_bound(p_upper,1)=ubound(efield_re%f2,dim=2)-array_bound(p_upper,1)
  array_bound(p_upper,2)=ubound(efield_re%f2,dim=3)-array_bound(p_upper,2)


  true_time = real(t, p_k_fld)
  local_time= true_time - this%delay
  local_omega = (this%omega0+this%chirp/2*local_time)

  if(this%jitteromg.ne.0) then
    phase=(local_omega*local_time)+this%phase+this%phase2(1)*2.0 &
          -this%jitterw/this%jitteromg*cos(this%jitteromg*local_time)
  else
    phase=(local_omega*local_time)+this%phase+this%phase2(1)*2.0
  end if


  if(this%ant_type.eq.p_env_linear) then
    if(this%jitteromg.ne.0) then
      phase=local_omega*local_time+this%phase &
            -this%jitterw/this%jitteromg*cos(this%jitteromg*local_time)
    else
      phase=local_omega*local_time+this%phase
    end if
    a_t=local_omega*envelop_env_linear(this,true_time)*sin(phase)
  else if(this%ant_type.eq.p_env_gaussian) then
    if(this%jitteromg.ne.0) then
      phase=local_omega*local_time+this%phase &
            -this%jitterw/this%jitteromg*cos(this%jitteromg*local_time)
    else
      phase=local_omega*local_time+this%phase
    end if
    a_t=this%a0*local_omega*envelop_env_gaussian(this,true_time)*sin(phase)
  else if(this%ant_type.eq.p_env_supergaussian) then
    if(this%jitteromg.ne.0) then
      phase=local_omega*local_time+this%phase &
            -this%jitterw/this%jitteromg*cos(this%jitteromg*local_time)
    else
      phase=local_omega*local_time+this%phase
    end if
    a_t=this%a0*local_omega*envelop_env_supergaussian(this,true_time)*sin(phase)
  else if(this%ant_type.eq.p_env_alfven) then
    if(local_time > 0.0) then
      phase=this%omega0*local_time+this%phase
      a_t=this%a0 * sin(phase) * envelop_env_gaussian(this,true_time)
    else
      a_t=0.0
    end if
  end if

  !write(*,*)'a_t=',a_t,phase,true_time


  ! --- Take care of spatial dependence here ---

  if(this%direction==1) then

    if(this%side.eq.p_lower) then
      ind_wall=array_bound(p_lower,1)
    else
      ind_wall=array_bound(p_upper,1)
    end if

    do ind1=lbound(efield_re%f2,dim=3),array_bound(p_upper,2)
      coord(1)=lower_bound(2)+delta_x(2)*&
              & (ind1-array_bound(p_lower,2))-this%x0
      r(1)=abs(coord(1)-this%x0)/this%rad_x
      a_tx=a_t*dt/delta_x(1)*profile(r(1))*&
          & cos(this%curv_const(1)*(coord(1))*(coord(1)))

      efield_re%f2(2,ind_wall,ind1)= efield_re%f2(2,ind_wall,ind1)+&
                                   & a_tx*sin(this%pol)
      efield_re%f2(3,ind_wall,ind1)= efield_re%f2(3,ind_wall,ind1)+&
                                   & a_tx*cos(this%pol)

      efield_im%f2(3,ind_wall,ind1)= efield_im%f2(3,ind_wall,ind1)+&
                                   & a_tx*sin(this%pol)
      efield_im%f2(2,ind_wall,ind1)= efield_im%f2(2,ind_wall,ind1)+&
                                   & a_tx*cos(this%pol)

    end do

  end if

end subroutine antenna_2d_cyl_modes
!-----------------------------------------------------------------------------------------

! ************************************************
! ************************************************
! various hermite polynomials
! ************************************************
! ************************************************
!-----------------------------------------------------------------------------------------
function l1_1(x)
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
  real (p_k_fld), intent(in) :: x
  real (p_k_fld) :: l1_1

  l1_1=x

end function l1_1

! ************************************************
! ************************************************
!-----------------------------------------------------------------------------------------
function l2_1(x)
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
  implicit none
  real (p_k_fld), intent(in) :: x
  real (p_k_fld) :: l2_1

  l2_1=0.5*(3.0*x*x-1.0)

end function l2_1
!-----------------------------------------------------------------------------------------

! ************************************************
! ************************************************
!-----------------------------------------------------------------------------------------
function l3_1(x)
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
  implicit none
  real (p_k_fld), intent(in) :: x
  real (p_k_fld) :: l3_1

  l3_1= 0.5*x*(5.0*x*x-3.0)

end function l3_1
!-----------------------------------------------------------------------------------------

! ************************************************
! ************************************************
!-----------------------------------------------------------------------------------------
function l4_1(x)
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
  real (p_k_fld), intent(in) :: x
  real (p_k_fld) :: l4_1
  real (p_k_fld) :: x2

  x2=x*x
  l4_1=0.125*x2*(35.0*x2-30.0)+0.375

end function l4_1
!-----------------------------------------------------------------------------------------

! ************************************************
! ************************************************
!-----------------------------------------------------------------------------------------
function l5_1(x)
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
  implicit none

  real (p_k_fld), intent(in) :: x
  real (p_k_fld) :: l5_1
  real (p_k_fld) :: x2

  x2=x*x
  l5_1=0.125*x*(63.0*x2*x2-70.0*x2+15)

end function l5_1
!-----------------------------------------------------------------------------------------

! ************************************************
! ************************************************
!-----------------------------------------------------------------------------------------
pure function l6_1(x)
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
  implicit none

  real (p_k_fld), intent(in) :: x
  real (p_k_fld) :: l6_1
  real (p_k_fld) :: x2

  x2=x*x
  l6_1=0.0625*(x2*(105.0+x2*(-315.0+231.0*x2))) - 0.3125

end function l6_1
!-----------------------------------------------------------------------------------------

!! ************************************************
!! ************************************************
!-----------------------------------------------------------------------------------------
function envelop_env_linear(this,time)
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
  implicit none
  class (t_antenna),intent(in):: this
  real (p_k_fld),intent(in):: time

  real (p_k_fld) :: envelop_env_linear

  real (p_k_fld) :: t_normal,local_time

  local_time=time-this%delay
  if(local_time < 0) then
    envelop_env_linear=0.0
  else if(local_time < this%t_rise) then
    t_normal=(local_time/this%t_rise)
    envelop_env_linear=t_normal
  else if(local_time < (this%t_rise+this%t_flat)) then
    envelop_env_linear=1.0
  else if(local_time < (this%t_rise+this%t_flat+this%t_fall)) then
    t_normal=1.0-(local_time-(this%t_rise+this%t_flat))/this%t_fall
    envelop_env_linear=t_normal
  else
    envelop_env_linear=0.0
  end if
  envelop_env_linear=2.0*envelop_env_linear

end function envelop_env_linear
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
function envelop_env_gaussian(this,time)
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
  implicit none
  class (t_antenna),intent(in):: this
  real (p_k_fld),intent(in):: time

  real (p_k_fld) :: envelop_env_gaussian

  real (p_k_fld) :: t_normal,local_time

  local_time=time-this%delay
  if(local_time < 0) then
    envelop_env_gaussian=0.0
  else if(local_time < this%t_rise) then
    t_normal=(local_time/this%t_rise)
    envelop_env_gaussian=t_normal*t_normal*t_normal*&
                          (10.0-t_normal*(15.0-6.0*t_normal))
  else if(local_time < (this%t_rise+this%t_flat)) then
    envelop_env_gaussian=1.0
  else if(local_time < (this%t_rise+this%t_flat+this%t_fall)) then
    t_normal=1.0-(local_time-(this%t_rise+this%t_flat))/this%t_fall
    envelop_env_gaussian=t_normal*t_normal*t_normal*&
                          (10.0-t_normal*(15.0-6.0*t_normal))
  else
    envelop_env_gaussian=0.0
  end if
  envelop_env_gaussian=2.0*envelop_env_gaussian

end function envelop_env_gaussian
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
function envelop_env_supergaussian(this,time)
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
  implicit none
  class (t_antenna),intent(in):: this
  real (p_k_fld),intent(in):: time

  real (p_k_fld) :: envelop_env_supergaussian

  real (p_k_fld) :: t_normal,local_time

  local_time=time-this%delay
  if(local_time < 0) then
    envelop_env_supergaussian=0.0
  else if(local_time < this%t_rise) then
    t_normal=(1-(local_time/this%t_rise))*3.0
    envelop_env_supergaussian=exp(-(t_normal)**2/2-this%eps*(t_normal)**this%power)
  else if(local_time < (this%t_rise+this%t_flat)) then
    envelop_env_supergaussian=1.0
  else if(local_time < (this%t_rise+this%t_flat+this%t_fall)) then
    t_normal=((local_time-(this%t_rise+this%t_flat))/this%t_fall)*3.0
    envelop_env_supergaussian=exp(-(t_normal)**2/2-this%eps*(t_normal)**this%power)

  else
    envelop_env_supergaussian=0.0
  end if
  envelop_env_supergaussian=2.0*envelop_env_supergaussian

end function envelop_env_supergaussian
!-----------------------------------------------------------------------------------------

! memory allocation routines omitted here

end module m_antenna_cyl_modes
