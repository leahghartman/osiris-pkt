module m_simulation_rad_cyl

#include "os-config.h"
#include "os-preprocess.fpp"

!use m_simulation, only : t_simulation
use m_simulation_cyl_modes, only: t_simulation_cyl_modes

use m_system

private

type, extends( t_simulation_cyl_modes ) :: t_simulation_rad_cyl

contains

  procedure :: allocate_objs  => allocate_objs_rad_cyl
  procedure :: list_algorithm => list_algorithm_rad_cyl

end type t_simulation_rad_cyl

public :: t_simulation_rad_cyl

contains

!-----------------------------------------------------------------------------------------
! Initialize simulation objects
!-----------------------------------------------------------------------------------------
subroutine allocate_objs_rad_cyl( sim )

  use m_particles_rad_cyl

  implicit none

  class( t_simulation_rad_cyl ), intent(inout) :: sim

  SCR_ROOT('Allocating Radiative particle objects.')

  allocate( t_particles_radiat_cyl :: sim%part )

  ! allocate any remaining objects in superclass
  call sim % t_simulation_cyl_modes % allocate_objs()

end subroutine allocate_objs_rad_cyl
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine list_algorithm_rad_cyl( sim )

  implicit none

  class( t_simulation_rad_cyl ), intent(in) :: sim

  print *, ""
  print *, "Simulation mode: Cylindrical RaDiO"
  call sim % t_simulation % list_algorithm()

end subroutine list_algorithm_rad_cyl
!-----------------------------------------------------------------------------------------

end module m_simulation_rad_cyl
