#include "os-config.h"
#include "os-preprocess.fpp"

module m_detlist

use m_system
use m_parameters
use detector_cls
use class_SpheDet
use class_CartDet
use m_input_file, only : t_input_file
use m_space, only : t_space
use m_node_conf, only : t_node_conf
use m_time_step, only : t_time_step
use m_restart, only : t_restart_handle


type :: t_detector_list

  class( detector ), pointer :: list => null()

contains

  procedure :: read_input => read_input_detector_list
  procedure :: init       => init_detector_list
  procedure :: report     => report_detector_list
  procedure :: cleanup    => cleanup_detector_list
  procedure :: fill_emf   => fill_emf_detector_list
  procedure :: restart_write => restart_write_detector_list
  procedure :: restart_read => restart_read_detector_list

 end type t_detector_list

contains

!-----------------------------------------------------------------------------------------
subroutine read_input_detector_list( this, input_file, dt )

  implicit none

  class( t_detector_list ), intent(inout) :: this
  class( t_input_file ), intent(inout) :: input_file
  real(p_double),intent(in) :: dt

  class( Detector ), pointer :: tail, item
  integer:: place

  place = 0

  tail => null()
  do
    select case (trim(input_file % get_section_name()))
    case('nl_spherical')
      allocate( sphe_Detector :: item )
    case('nl_cartesian')
      allocate( cart_Detector :: item )
    case default
      ! Not a detector section, exit
      exit
    end select

    ! Read input for detector
    call item % read_input( input_file, dt )
    item % next => null()
    place=place+1
    item%place=place
    if (associated(tail)) then
      tail % next => item
    else
      this % list => item
    endif

    tail => item

  enddo

end subroutine read_input_detector_list
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine init_detector_list( this, sp_name, no_co, restart, restart_handle )

  implicit none

  class( t_detector_list ), intent(inout) :: this
  character(len = 64) ,intent(in):: sp_name
  class( t_node_conf ), intent(in) :: no_co
  logical, intent(in) :: restart
  type( t_restart_handle ), intent(in) :: restart_handle

  class( Detector ), pointer :: det

  det => this%list
  do
    if ( .not. associated(det)) exit
    call det % init(sp_name,no_co,restart,restart_handle)
    det => det % next
  enddo

end subroutine init_detector_list
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine cleanup_detector_list( this )

  implicit none

  class( t_detector_list ), intent(inout) :: this
  class( Detector ), pointer :: det, next

  det => this%list
  do
    if ( .not. associated(det)) exit

    next => det % next
    call det % destruct()
    deallocate( det )
    det => next
  enddo

end subroutine cleanup_detector_list
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine report_detector_list( this,t,tstep,no_co )

  implicit none

  class( t_detector_list ), intent(inout) :: this
  real(p_double), intent(in) :: t
  type( t_time_step ), intent(in) :: tstep
  class( t_node_conf ), intent(in) :: no_co

  class( Detector ), pointer :: det

  det => this%list
  do
    if ( .not. associated(det)) exit
    call det % report(t,tstep,no_co)
    det => det % next
  enddo

end subroutine report_detector_list
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine fill_emf_detector_list(this,x,beta,time,xprev,betaprev,charge,np,tid)

  implicit none

  class( t_detector_list ), intent(inout) :: this
  real(p_k_part), dimension(3,p_cache_size), intent(in) :: x, xprev
  real(p_k_part), dimension(p_p_dim,p_cache_size), intent(in) :: beta, betaprev
  real(p_k_part), dimension(p_cache_size), intent(in) :: charge
  real(p_double), intent(in) :: time
  integer, intent(in) :: np, tid

  class( Detector ), pointer :: det

  det => this%list
  do
    if ( .not. associated(det)) exit
    call det % fill_emf(x,beta,time,xprev,betaprev,charge,np,tid)
    det => det % next
  enddo

end subroutine fill_emf_detector_list
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine restart_write_detector_list( this, restart_handle )

  implicit none

  class( t_detector_list ), intent(in) :: this
  type( t_restart_handle ), intent(inout) :: restart_handle

  class( Detector ), pointer :: det

  det => this%list
  do
    if ( .not. associated(det) ) exit
    call det % restart_write( restart_handle )
    det => det % next
  enddo

end subroutine restart_write_detector_list
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine restart_read_detector_list( this, restart_handle )

  implicit none

  class( t_detector_list ), intent(inout) :: this
  type( t_restart_handle ), intent(in) :: restart_handle

  class( Detector ), pointer :: det

  det => this%list
  do
    if ( .not. associated(det) ) exit
    call det % restart_read( restart_handle )
    det => det % next
  enddo

end subroutine restart_read_detector_list
!-----------------------------------------------------------------------------------------

end module m_detlist
