module m_simulation_cyl_modes

#include "os-config.h"
#include "os-preprocess.fpp"

use m_system
use m_simulation, only : t_simulation

implicit none

private

type, extends( t_simulation ) :: t_simulation_cyl_modes

contains
  procedure :: allocate_objs => allocate_objs_cyl_modes

end type t_simulation_cyl_modes

public :: t_simulation_cyl_modes

contains

!-----------------------------------------------------------------------------------------
! Initialize simulation objects
!-----------------------------------------------------------------------------------------
subroutine allocate_objs_cyl_modes( sim )

  use m_grid_cyl_modes, only : t_grid_cyl_modes
  use m_zpulse_cyl_modes, only : t_zpulse_list_cyl_modes
  use m_emf_cyl_modes, only : t_emf_cyl_modes
  use m_particles_cyl_modes, only : t_particles_cyl_modes
  use m_current_cyl_modes, only : t_current_cyl_modes
  use m_antenna_cyl_modes, only : t_antenna_array_cyl_modes

  implicit none

  class( t_simulation_cyl_modes ), intent(inout) :: sim

  SCR_ROOT('Allocating cyl_modes zpulse, emf, and particle objects.')

  allocate( t_grid_cyl_modes :: sim%grid )
  allocate( t_zpulse_list_cyl_modes :: sim%zpulse_list )
  allocate( t_antenna_array_cyl_modes :: sim%antenna_array )
  allocate( t_emf_cyl_modes :: sim%emf )
  if ( .not. associated( sim%part ) ) then
    allocate( t_particles_cyl_modes :: sim%part )
  end if
  allocate( t_current_cyl_modes :: sim%jay )

  ! allocate any remaining objects in superclass
  call sim % t_simulation % allocate_objs()

end subroutine allocate_objs_cyl_modes




end module m_simulation_cyl_modes
