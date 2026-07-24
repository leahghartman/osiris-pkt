#include "os-config.h"
#include "os-preprocess.fpp"

module m_particles_shear

#include "memory/memory.h"

use m_system
use m_space
use m_grid_define
use m_node_conf
use m_restart
use m_time_step
use m_input_file

use m_particles_define
use m_species_define
use m_current_define

implicit none

private

type, extends( t_particles ) :: t_particles_shear

  ! no additional class members
  ! all species must be of type species_shear
  ! all neutrals must be of type neutral_shear

contains

  procedure :: advance_deposit => advance_deposit_shear
  procedure :: allocate_objs => allocate_objs_part_shear

end type t_particles_shear

public :: t_particles_shear

contains

!-----------------------------------------------------------------------------------------
! Overrides allocate_objs method and allocates t_species_shear and t_neutral_shear objects
! instead
!-----------------------------------------------------------------------------------------
subroutine allocate_objs_part_shear( this )

  use m_species_define
  use m_species_shear
  use m_cathode

  implicit none

  class( t_particles_shear ), intent(inout) :: this

  integer :: i
  class(t_species_shear), pointer :: spec, tail

  if ( this%num_species > 0 ) then
    allocate( t_species_shear :: spec )
    this % species => spec
    do i = 2, this%num_species
      tail => spec
      allocate( t_species_shear :: spec )
      tail % next => spec
    enddo

    if ( this%num_cathode > 0 ) then
      allocate( this%cathode( this%num_cathode ) )
    endif
  endif

  ! currently ionization is not supported
  if ( this%num_neutral > 0 ) then
    SCR_ROOT('(*warning*) Ionization is currently disabled in shear mode')
  endif

end subroutine allocate_objs_part_shear
!-----------------------------------------------------------------------------------------

subroutine advance_deposit_shear( this, emf, jay, tstep, t, no_co, options )

  use m_emf_define
  use m_emf_shear
  use m_vdf_define
  use m_vdf_math
  use m_time_step
  use m_node_conf
  use m_current_define

  implicit none

  class( t_particles_shear ), intent(inout) :: this
  class( t_emf ), intent( inout )  ::  emf
  class( t_current ), intent(inout) :: jay

  real(p_double), intent(in) :: t
  type( t_time_step ), intent(in) :: tstep

  class( t_node_conf ), intent(in) :: no_co
  type( t_options ), intent(in) :: options

  ! call superclass pusher
  call this % t_particles % advance_deposit( emf, jay, tstep, t, no_co, options )

end subroutine advance_deposit_shear

end module m_particles_shear
