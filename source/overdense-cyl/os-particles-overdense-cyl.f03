
#include "os-config.h"
#include "os-preprocess.fpp"

module m_particles_overdense_cyl

#include "memory/memory.h"

use m_system
use m_parameters
use m_particles_cyl_modes, only: t_particles_cyl_modes
use m_particles_overdense
use m_spec_overdense_cyl
use m_particles_define, only : pushev
use m_spec_overdense
use m_spec_splitting, only : if_split, t_recombine
use m_species_define, only : t_species
use m_species_memory
use m_emf_define, only : t_emf
use m_time_step, only : t_time_step, n
use m_current_define, only : t_current
use m_logprof, only : begin_event, end_event, create_event
use m_node_conf, only : t_node_conf
use m_grid_define, only : t_grid
use m_restart, only : t_restart_handle
use m_space, only : t_space
use m_bnd, only : t_bnd
use m_neutral_cyl_modes,        only: t_neutral_cyl_modes
use m_zpulse, only : t_zpulse_list

implicit none

private

! particle type, extends the cyl_modes type and add arguments/methods from overdense
type, extends( t_particles_cyl_modes ) :: t_particles_overdense_cyl

contains

  procedure :: allocate_objs   => allocate_objs_part_overdense_cyl
  procedure :: init            => init_part_overdense_cyl
  procedure :: advance_deposit => advance_deposit_part_overdense_cyl

end type t_particles_overdense_cyl

public :: t_particles_overdense_cyl

contains

!-----------------------------------------------------------------------------------------
! Overrides allocate_objs method and allocates t_species_overdense_cyl objects instead
!-----------------------------------------------------------------------------------------
subroutine allocate_objs_part_overdense_cyl( this )

  implicit none

  class( t_particles_overdense_cyl ), intent(inout) :: this

  integer :: i
  class(t_species_overdense_cyl), pointer :: spec, tail

  if ( this%num_species > 0 ) then

    allocate( t_species_overdense_cyl :: spec )
    this % species => spec
    do i = 2, this%num_species
      tail => spec
      allocate( t_species_overdense_cyl :: spec )
      tail % next => spec
    enddo

    if ( this%num_cathode > 0 ) then
      allocate( this%cathode( this%num_cathode ) )
    endif

#ifdef __HAS_IONIZATION__
    if ( this%num_neutral > 0 ) then
    allocate(t_neutral_cyl_modes :: this%neutral( this%num_neutral ))
    endif
#endif

  endif

end subroutine allocate_objs_part_overdense_cyl
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Sets up this data structure from the given information
!-----------------------------------------------------------------------------------------
subroutine init_part_overdense_cyl( this, g_space, jay, emf, grid, no_co, bnd, &
                                    zpulse_list, ndump_fac, restart, restart_handle, t, &
                                    tstep, tmin, tmax, sim_options )

  implicit none

  class( t_particles_overdense_cyl ), intent(inout) :: this

  type( t_space ),     intent(in) :: g_space
  class( t_emf ) ,intent(inout) :: emf
  class( t_current ) ,intent(inout) :: jay
  class( t_grid ), intent(in) :: grid
  class( t_node_conf ), intent(in) :: no_co
  class( t_bnd ), intent(inout) :: bnd
  class( t_zpulse_list ), intent(in) :: zpulse_list
  integer, intent(in) :: ndump_fac
  logical, intent(in) :: restart
  type( t_restart_handle ), intent(in) :: restart_handle
  real(p_double), intent(in) :: t
  type(t_time_step), intent(in) :: tstep
  real(p_double), intent(in) :: tmin, tmax
  type( t_options ), intent(in) :: sim_options

  call this % t_particles_cyl_modes%init( g_space, jay, emf, grid, no_co, bnd, &
                                          zpulse_list, ndump_fac, restart,restart_handle,&
                                          t, tstep, tmin, tmax, sim_options )

  ! setup events in this file
  if (split_ev==0) then
    split_ev       = create_event('particle splitting')
    recombine_ev   = create_event('particle recombination')
    prep_damper_ev = create_event('particle damping preparation')
  endif

end subroutine init_part_overdense_cyl
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Prepare fluid moment arrays if damping
!-----------------------------------------------------------------------------------------
subroutine advance_deposit_part_overdense_cyl( this, emf, jay, tstep, t, no_co, options )

  implicit none

  class( t_particles_overdense_cyl ), intent(inout) :: this
  class( t_emf ), intent( inout )  ::  emf
  class( t_current ), intent(inout) :: jay

  real(p_double), intent(in) :: t
  type( t_time_step ), intent(in) :: tstep
  class( t_node_conf ), intent(in) :: no_co
  type( t_options ), intent(in) :: options

  class(t_species), pointer :: species
  type(t_recombine), pointer :: recombine

  ! prepare the damper if required
  call begin_event( pushev )
  species => this % species
  do
    if (.not. associated(species)) exit
    select type(species)
    class is (t_species_overdense_cyl)
      if ( t >= species%push_start_time ) then

        call begin_event( prep_damper_ev )
        call species % prep_damper( n(tstep), t )
        call end_event( prep_damper_ev )

      endif
    end select
    species => species % next
  enddo
  call end_event( pushev )

  call this % t_particles_cyl_modes % advance_deposit( emf, jay, tstep, t, no_co, options)

  ! split and recombine particles if required
  call begin_event( pushev )
  species => this % species
  do
    if (.not. associated(species)) exit
    select type(species)
    class is (t_species_overdense_cyl)
      if ( t >= species%push_start_time ) then

        call begin_event( split_ev )
        call species % split( n(tstep) )
        call end_event( split_ev )

        call begin_event( recombine_ev )
        recombine => species%splitter%recombine
        do
          if (.not. associated(recombine)) exit
          call species % recombine( recombine, n(tstep), t )
          recombine => recombine%next
        enddo
        call end_event( recombine_ev )

      endif
    end select
    species => species % next
  enddo
  call end_event( pushev )

end subroutine advance_deposit_part_overdense_cyl
!-----------------------------------------------------------------------------------------


end module m_particles_overdense_cyl
