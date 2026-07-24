module m_simulation_shear

#include "os-config.h"
#include "os-preprocess.fpp"

use m_system
use m_parameters
use m_simulation

implicit none

private

type, extends( t_simulation ) :: t_simulation_shear

contains

   procedure :: allocate_objs      => allocate_objs_shear

end type t_simulation_shear

public :: t_simulation_shear

contains

!-----------------------------------------------------------------------------------------
! Initialize simulation objects
!-----------------------------------------------------------------------------------------
subroutine allocate_objs_shear( sim )

  use m_emf_shear
  use m_particles_shear

  implicit none

  class( t_simulation_shear ), intent(inout) :: sim

  SCR_ROOT('Allocating SHEAR emf and particle objects.')

  allocate( t_emf_shear :: sim%emf )
  allocate( t_particles_shear :: sim%part )

  ! allocate any remaining objects in superclass
  call sim % t_simulation % allocate_objs()

end subroutine allocate_objs_shear




end module m_simulation_shear
