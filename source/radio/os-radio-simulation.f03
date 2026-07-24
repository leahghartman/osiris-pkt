module m_simulation_rad

#include "os-config.h"
#include "os-preprocess.fpp"

use m_simulation, only : t_simulation
use m_system

private

type, extends( t_simulation ) :: t_simulation_rad

contains

  procedure :: allocate_objs  => allocate_objs_rad
  procedure :: list_algorithm => list_algorithm_rad

end type t_simulation_rad

public :: t_simulation_rad

contains

!-----------------------------------------------------------------------------------------
! Initialize simulation objects
!-----------------------------------------------------------------------------------------
subroutine allocate_objs_rad( sim )

  use m_particles_rad

  implicit none

  class( t_simulation_rad ), intent(inout) :: sim

  SCR_ROOT('Allocating Radiative particle objects.')

  allocate( t_particles_rad :: sim%part )

  ! allocate any remaining objects in superclass
  call sim % t_simulation % allocate_objs()

end subroutine allocate_objs_rad
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine list_algorithm_rad( sim )

  implicit none

  class( t_simulation_rad ), intent(in) :: sim

  print *, ""
  print *, "Simulation mode: RaDiO"
  call sim % t_simulation % list_algorithm()

end subroutine list_algorithm_rad
!-----------------------------------------------------------------------------------------

end module m_simulation_rad