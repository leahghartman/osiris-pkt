#include "os-config.h"
#include "os-preprocess.fpp"

module m_profile_cyl_modes

#include "memory/memory.h"

use m_psource_std_cyl_modes
use m_psource_beam_cyl_modes
use m_psource_file
use m_psource_py_cyl_modes
use m_species_define

implicit none

private

interface new_source_cyl_modes
    module procedure new_source_cyl_modes
end interface new_source_cyl_modes

public :: new_source_cyl_modes

contains

function new_source_cyl_modes( type, input_file, coordinates )

  use m_input_file

  implicit none

  class( t_psource ), pointer :: new_source_cyl_modes
  character(len=*), intent(in) :: type
  class( t_input_file ), intent(inout) :: input_file
  integer, intent(in) :: coordinates
  ! create t_profile object of the selected kind
  select case ( trim(type) )
  case('standard','profile')
    allocate( t_psource_std_cyl_modes :: new_source_cyl_modes )
  ! case('constq')
  !   allocate( t_psource_constq :: new_source )
  case('beamfocus')
    allocate( t_psource_beam_cyl_modes :: new_source_cyl_modes )
  ! case('betatronfel')
  !   allocate( t_psource_betatronfel :: new_source )
  case('file')
    allocate( t_psource_file :: new_source_cyl_modes )
  case('python')
    allocate( t_psource_py_cyl_modes :: new_source_cyl_modes )
  case default
    new_source_cyl_modes => null()
  end select

  ! If successfull read input parameters and finish initalization
  if ( associated(new_source_cyl_modes) ) then
    call new_source_cyl_modes % read_input( input_file, coordinates )
  endif

end function new_source_cyl_modes

end module m_profile_cyl_modes