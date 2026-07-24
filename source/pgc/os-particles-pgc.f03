#include "os-config.h"
#include "os-preprocess.fpp"

module m_particles_pgc

#include "memory/memory.h"

use m_system
use m_parameters
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

type, extends( t_particles ) :: t_particles_pgc

  ! no additional class members
  ! all species must be of type species_pgc
  ! all neutrals must be of type neutral_pgc

contains

  procedure :: advance_deposit => advance_deposit_pgc
  procedure :: allocate_objs => allocate_objs_part_pgc

end type t_particles_pgc

public :: t_particles_pgc

contains

!-----------------------------------------------------------------------------------------
! Overrides allocate_objs method and allocates t_species_pgc and t_neutral_pgc objects
! instead
!-----------------------------------------------------------------------------------------
subroutine allocate_objs_part_pgc( this )

  use m_species_define
  use m_species_pgc
  use m_neutral_pgc
  use m_cathode

  implicit none

  class( t_particles_pgc ), intent(inout) :: this

  integer :: i
  class( t_species_pgc ), pointer :: spec, tail

  if ( this % num_species > 0 ) then

    allocate(spec)
    this % species => spec
    do i = 2, this % num_species
      tail => spec
      allocate(spec)
      tail % next => spec
    enddo

    if ( this % num_cathode > 0 ) then
      allocate( this%cathode( this%num_cathode ) )
    endif

#ifdef __HAS_IONIZATION__
	  if ( this%num_neutral > 0 ) then
      allocate( t_neutral_pgc :: this%neutral( this%num_neutral ) )
    endif
#endif

  endif

end subroutine allocate_objs_part_pgc
!-----------------------------------------------------------------------------------------

subroutine advance_deposit_pgc( this, emf, jay, tstep, t, no_co, options )

  use m_emf_define,     only : t_emf
  use m_emf_pgc_define, only : t_emf_pgc
  use m_vdf_define
  use m_vdf_math
  use m_time_step
  use m_node_conf
  use m_current_define

  implicit none

  class( t_particles_pgc ), intent(inout) :: this
  class( t_emf ), intent( inout )  ::  emf
  class( t_current ), intent(inout) :: jay

  real(p_double), intent(in) :: t
  type( t_time_step ), intent(in) :: tstep

  class( t_node_conf ), intent(in) :: no_co
  type( t_options ), intent(in) :: options

  !print *, '[pgc] In advance_deposit_pgc'

  ! before calling the pusher we must reset chi

  select type( emf )

  class is ( t_emf_pgc )

	! reset_chi() only works with t_emf_pgc objects, and since emf was declared as
	! t_emf we must use a select type construct to check that the object is indeed a
	! t_emf_pgc instance

    call emf % reset_chi()

  class default

    ! This must never happen, advance_deposit_pgc must always be called with a
    ! t_emf_pgc object
    call abort_program( p_err_invalid )

  end select

  ! call superclass pusher
  call this % t_particles % advance_deposit( emf, jay, tstep, t, no_co, options )

  ! reduce chi
  select type( emf )
  class is ( t_emf_pgc )

    call reduce( emf%chi, 1 )

  class default

    ! This must never happen, advance_deposit_pgc must always be called with a
    ! t_emf_pgc object
    call abort_program( p_err_invalid )

  end select

end subroutine advance_deposit_pgc

end module m_particles_pgc
