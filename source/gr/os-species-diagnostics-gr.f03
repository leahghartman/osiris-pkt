#include "os-config.h"
#include "os-preprocess.fpp"

module m_species_diagnostics_gr

#include "memory/memory.h"

use m_system
use m_parameters

use m_species_define_gr, only : t_species_gr

implicit none

private

interface report_nemit_gr
    module procedure report_nemit_gr
end interface report_nemit_gr

public :: report_nemit_gr

contains

subroutine report_nemit_gr( this, g_space, grid, no_co, tstep, t )
  
  use stringutil
  use m_vdf_define, only : t_vdf_report
  use m_space, only : t_space
  use m_grid_define, only : t_grid
  use m_node_conf, only : t_node_conf
  use m_time_step

  implicit none

  class( t_species_gr ), intent(inout) :: this
  type( t_space ),      intent(in) :: g_space
  class( t_grid ),      intent(in) :: grid
  class( t_node_conf ),  intent(in) :: no_co
  type( t_time_step ),  intent(in) :: tstep
  real( p_double ) :: t

  type( t_vdf_report ) :: info

  ! Prepare the metadata
  info%name = 'nemit'
  info%ndump = ndump(tstep)

  info%xname  = (/'x1', 'x2', 'x3'/)

  info%xlabel = (/'x_1', 'x_2', 'x_3'/)
  info%xunits = (/'c / \omega_p', 'c / \omega_p', 'c / \omega_p'/)

  info%time_units = '1 / \omega_p'
  info%dt         = dt(tstep)

  info%label = '# emissions'
  info%units = ''

  info%n = n(tstep)
  info%t = t

  info%fileLabel = ''
  info%path = trim(path_mass) // &
              'GROUP' // p_dir_sep // &
              trim(this%name) // p_dir_sep // &
              'nemit' // p_dir_sep

  info%filename = 'nemit-' // trim(this%name) // '-' // idx_string( n(tstep)/ndump(tstep), 6 )

  ! write the file
  call this % nemit % write( info, 1, g_space, grid, no_co )
  
  ! reset vdf values
  call this % nemit % zero()
  
end subroutine report_nemit_gr

end module m_species_diagnostics_gr