module m_simulation_pgc

#include "os-config.h"
#include "os-preprocess.fpp"

use m_system
use m_parameters
use m_simulation

implicit none

private

type, extends( t_simulation ) :: t_simulation_pgc

contains

   procedure :: allocate_objs      => allocate_objs_pgc

end type t_simulation_pgc

public :: t_simulation_pgc

contains

!-----------------------------------------------------------------------------------------
! Initialize simulation objects
!-----------------------------------------------------------------------------------------
subroutine allocate_objs_pgc( sim )

  use m_emf_pgc_define
  use m_particles_pgc

  implicit none

  class( t_simulation_pgc ), intent(inout) :: sim

  SCR_ROOT('Allocating PGC emf and particle objects.')

  allocate( t_emf_pgc :: sim%emf )
  allocate( t_particles_pgc :: sim%part )

  ! allocate any remaining objects in superclass
  call sim % t_simulation % allocate_objs()

end subroutine allocate_objs_pgc




end module m_simulation_pgc
