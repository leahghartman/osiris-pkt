module m_simulation_pkt

#include "os-config.h"
#include "os-preprocess.fpp"

use m_simulation
use m_system

! this is just meant to prevent re-exporting m_simulation symbols, everything is public
private

type, extends( t_simulation ) :: t_simulation_pkt

contains

   procedure :: allocate_objs      => allocate_objs_pkt

end type t_simulation_pkt

public :: t_simulation_pkt

contains

!-----------------------------------------------------------------------------------------
! Initialize simulation objects
!-----------------------------------------------------------------------------------------
subroutine allocate_objs_pkt( sim )

  use m_particles_pkt

  implicit none

  class( t_simulation_pkt ), intent(inout) :: sim

  SCR_ROOT('Allocating pkt particle objects.')

  allocate( t_particles_pkt :: sim%part )

  ! allocate any remaining objects in superclass
  call sim % t_simulation % allocate_objs()

end subroutine allocate_objs_pkt

end module m_simulation_pkt
