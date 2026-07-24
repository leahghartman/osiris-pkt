#include "os-config.h"
#include "os-preprocess.fpp"

module m_particles_cuda

use m_particles_define

implicit none

private

type, extends( t_particles ) :: t_particles_cuda

contains

  procedure :: allocate_objs => allocate_objs_part_cuda

end type t_particles_cuda

public :: t_particles_cuda

contains

!-----------------------------------------------------------------------------------------
! Overrides allocate_objs method and allocates t_species_cuda instead
!-----------------------------------------------------------------------------------------
subroutine allocate_objs_part_cuda( this )

  use m_species_define
  use m_cathode
  use m_species_cuda

  implicit none

  class( t_particles_cuda ), intent(inout) :: this

  integer :: i
  class( t_species_cuda ), pointer :: spec, tail

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
      allocate( this%neutral( this%num_neutral ) )
    endif
#endif

  endif

end subroutine allocate_objs_part_cuda
!-----------------------------------------------------------------------------------------

end module m_particles_cuda
