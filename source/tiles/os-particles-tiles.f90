!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     particle functions for tile load balance
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "os-config.h"
#include "os-preprocess.fpp"

module m_particles_tiles

  use m_system
  use m_parameters
  use m_particles_define, only : t_particles
  use m_species_define, only : t_species
  use m_species_loadbalance_tiles, only : add_particle_load, add_density_load
  use m_grid_define, only : t_grid
  use m_grid_tiles, only : t_grid_group
  use m_node_conf, only : t_node_conf
  use m_node_conf_tiles, only : t_node_conf_tiles
  use m_tile_conf, only : p_topology_viktor

  implicit none

  private

  interface add_load
    module procedure add_load_particles_tiles
  end interface

  public :: add_load

 contains

!-----------------------------------------------------------------------------------------
! Sets the particle load.
! For iteration n = 0 the load is set for the whole simulation volume (since this is
! called when no node partitions exist yet), for n > 0 the load is set for the
! local volume/particles only
!-----------------------------------------------------------------------------------------
subroutine add_load_particles_tiles( this, grid, no_co, n )

  implicit none

  class( t_particles ), intent(inout) :: this
  class( t_grid ), intent(inout) :: grid
  class( t_node_conf ), intent(in) :: no_co
  integer, intent(in) :: n

  class( t_species ), pointer :: species

  ! loop through all species and add to int_load array
  species => this % species
  if ( n > 0 ) then

    select type(grid); class is(t_grid_group)
    select type(no_co); class is(t_node_conf_tiles)
    do
      if (.not. associated(species)) exit
      select case(no_co%ti_co%topology_type)
      case(p_topology_viktor)
        select case(p_x_dim)
        case(2)
          call add_particle_load( species, no_co, til_load_2=grid%til_load_2 )
        case(3)
          call add_particle_load( species, no_co, til_load_3=grid%til_load_3 )
        end select
      case default
        call add_particle_load( species, no_co, int_load=grid%int_load )
      end select
      species => species % next
    enddo
    end select
    end select

  else

    do
      if (.not. associated(species)) exit
      call add_density_load( species, grid, no_co )
      species => species % next
    enddo

  endif

end subroutine add_load_particles_tiles
!-----------------------------------------------------------------------------------------

end module m_particles_tiles
