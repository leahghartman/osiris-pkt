module m_simulation_neutral_spin

#include "os-config.h"
#include "os-preprocess.fpp"

use m_parameters
use m_simulation

implicit none

private

type, extends( t_simulation ) :: t_simulation_neutral_spin

contains
  procedure :: allocate_objs      => allocate_objs_neutral_spin

end type t_simulation_neutral_spin

public :: t_simulation_neutral_spin

contains

!-----------------------------------------------------------------------------------------
! Initialize simulation objects
!-----------------------------------------------------------------------------------------
subroutine allocate_objs_neutral_spin( sim )
  
  use m_particles_neutral_spin, only : t_particles_neutral_spin
  
  implicit none
  
  class( t_simulation_neutral_spin ), intent(inout) :: sim
  
  SCR_ROOT('Allocating particle objects for TDSE ionization.')

  allocate( t_particles_neutral_spin :: sim%part )
  
  ! allocate any remaining objects in superclass
  call sim % t_simulation % allocate_objs()

end subroutine allocate_objs_neutral_spin

end module m_simulation_neutral_spin