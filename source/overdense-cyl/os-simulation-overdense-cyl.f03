module m_simulation_overdense_cyl

#include "os-config.h"
#include "os-preprocess.fpp"

use m_system
use m_simulation_cyl_modes, only: t_simulation_cyl_modes
use m_particles_overdense_cyl, only: t_particles_overdense_cyl

implicit none

private

type, extends( t_simulation_cyl_modes ) :: t_simulation_overdense_cyl

contains
  procedure :: allocate_objs => allocate_objs_overdense_cyl

end type t_simulation_overdense_cyl

public :: t_simulation_overdense_cyl

contains

!-----------------------------------------------------------------------------------------
! Initialize simulation objects
!-----------------------------------------------------------------------------------------
subroutine allocate_objs_overdense_cyl( sim )

  implicit none

  class( t_simulation_overdense_cyl ), intent(inout) :: sim

  SCR_ROOT('Allocating overdense quasi-3D objects.')

  allocate( t_particles_overdense_cyl :: sim%part )

  ! allocate any remaining objects in superclass
  call sim % t_simulation_cyl_modes % allocate_objs()

end subroutine allocate_objs_overdense_cyl




end module m_simulation_overdense_cyl
