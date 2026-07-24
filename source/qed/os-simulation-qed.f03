module m_simulation_qed

#include "os-config.h"
#include "os-preprocess.fpp"

use m_simulation
use m_system

! this is just meant to prevent re-exporting m_simulation symbols, everything is public
private

type, extends( t_simulation ) :: t_simulation_qed

contains

   procedure :: allocate_objs      => allocate_objs_qed
   procedure :: report_load        => report_load_sim_qed

end type t_simulation_qed

public :: t_simulation_qed

contains

!-----------------------------------------------------------------------------------------
! Initialize simulation objects
!-----------------------------------------------------------------------------------------
subroutine allocate_objs_qed( sim )

  use m_particles_qed

  implicit none

  class( t_simulation_qed ), intent(inout) :: sim

  SCR_ROOT('Allocating QED particle objects.')

  allocate( t_particles_qed :: sim%part )

  ! allocate any remaining objects in superclass
  call sim % t_simulation % allocate_objs()

end subroutine allocate_objs_qed

!-----------------------------------------------------------------------------------------
! Report on the load of the simulation, including particles + photons.
!-----------------------------------------------------------------------------------------
subroutine report_load_sim_qed( this )
    
  use m_grid_define
  use m_particles
  use m_dynamic_loadbalance
  use m_simulation
  use m_time_step
  use m_time
  use m_logprof
  use m_particles_qed
  
  implicit none
  
  class ( t_simulation_qed ), intent(inout) ::  this
  
  call begin_event(report_load_ev)
  
  ! global load : total, max, min, avg parts per node
  if ( this%grid%if_report( n(this%tstep), p_global ) ) then
      call report_global_load( this%part, n(this%tstep), this%no_co )
  endif
  
  ! node load : number of particles per node for all nodes
  if ( this%grid%if_report( n(this%tstep), p_node ) ) then
      call report_node_load_part_qed( this%part, n(this%tstep), this%grid%ndump( p_node ), t(this%time), &
      this%no_co )
  endif
  
  ! grid load : number of particles per cell for all grid points
  if ( this%grid%if_report( n(this%tstep), p_grid ) ) then
      call report_grid_load( this%part, n(this%tstep), this%grid%ndump( p_grid ), t(this%time), &
      this%g_space, this%grid, this%no_co )
  endif
  
  call end_event(report_load_ev)
  
end subroutine report_load_sim_qed
!-----------------------------------------------------------------------------------------

end module m_simulation_qed