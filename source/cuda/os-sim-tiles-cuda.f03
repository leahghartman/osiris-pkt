#include "os-config.h"
#include "os-preprocess.fpp"

module m_sim_tiles_cuda

use m_simulation_tiles
use m_node_conf_cuda, only : t_node_conf_cuda
use m_system

#include "memory/memory.h"

implicit none

private

type, extends( t_simulation_tiles ) :: t_sim_tiles_cuda

contains

  procedure :: allocate_objs => allocate_objs_tiles_cuda

end type t_sim_tiles_cuda

public :: t_sim_tiles_cuda

contains

!-----------------------------------------------------------------------------------------
! Initialize simulation objects
!-----------------------------------------------------------------------------------------
subroutine allocate_objs_tiles_cuda( sim )

  use m_particles_cuda
  use m_current_cuda

  implicit none

  class( t_sim_tiles_cuda ), intent(inout) :: sim

  if ( .not. associated( sim%no_co ) ) then
    allocate( t_node_conf_cuda :: sim%no_co )
  endif

  if ( .not. associated( sim%jay ) ) then
    allocate( t_current_cuda :: sim%jay )
  endif

  if ( .not. associated( sim%part ) ) then
    allocate( t_particles_cuda :: sim%part )
  endif

  ! call superclass
  call sim % t_simulation_tiles % allocate_objs()

end subroutine allocate_objs_tiles_cuda
!-----------------------------------------------------------------------------------------


end module m_sim_tiles_cuda
