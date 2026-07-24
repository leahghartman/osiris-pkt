!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Particles module for TDSE model
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "os-config.h"
#include "os-preprocess.fpp"

module m_particles_neutral_spin

use m_system
use m_parameters
use m_particles_define, only : t_particles
use m_species_define,   only : t_species
use m_neutral_spin,     only : t_neutral_spin

implicit none

private

type, extends( t_particles ) :: t_particles_neutral_spin

contains

  procedure :: allocate_objs => allocate_objs_part_neutral_spin
  ! procedure :: init => init_part_neutral_spin
  ! procedure :: advance_deposit => advance_deposit_part_neutral_spin

end type t_particles_neutral_spin

public :: t_particles_neutral_spin


contains

!-----------------------------------------------------------------------------------------
! Allocate particle species, cathodes, etc.
!-----------------------------------------------------------------------------------------
subroutine allocate_objs_part_neutral_spin( this )

  implicit none

  class( t_particles_neutral_spin ), intent(inout) :: this

  integer :: i
  class(t_species), pointer :: spec, tail


  if ( this%num_species > 0 ) then

    allocate( spec )
    this % species => spec
    do i = 2, this%num_species
      tail => spec
      allocate( spec )
      tail % next => spec
    enddo

    if ( this%num_cathode > 0 ) then
      allocate( this%cathode( this%num_cathode ) )
    endif

#ifdef __HAS_IONIZATION__
    if ( this%num_neutral > 0 ) then
      allocate( t_neutral_spin :: this%neutral( this%num_neutral ) )
    endif
#endif

  endif

end subroutine allocate_objs_part_neutral_spin
!-----------------------------------------------------------------------------------------

end module m_particles_neutral_spin