module m_simulation_overdense

#include "os-config.h"
#include "os-preprocess.fpp"

use m_simulation
use m_system

private

type, extends( t_simulation ) :: t_simulation_overdense

contains

  procedure :: allocate_objs => allocate_objs_overdense

end type t_simulation_overdense

public :: t_simulation_overdense

contains

!-----------------------------------------------------------------------------------------
! Initialize simulation objects
!-----------------------------------------------------------------------------------------
subroutine allocate_objs_overdense( sim )

  use m_particles_overdense

  implicit none

  class( t_simulation_overdense ), intent(inout) :: sim

  SCR_ROOT('Allocating damping particle objects.')

  allocate( t_particles_overdense :: sim%part )

  ! allocate any remaining objects in superclass
  call sim % t_simulation % allocate_objs()

end subroutine allocate_objs_overdense


end module m_simulation_overdense