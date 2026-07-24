#include "os-config.h"
#include "os-preprocess.fpp"

module m_species_loadbalance_tiles

#include "memory/memory.h"

use m_parameters
use m_species_define, only : t_species, t_psource, p_cell_near
use m_node_conf, only : t_node_conf
use m_node_conf_tiles, only : t_node_conf_tiles
use m_tile_conf, only : p_topology_viktor
use m_grid_define, only : t_grid
use m_grid_tiles, only : t_grid_group
use m_psource_std, only : t_psource_std

implicit none

private

interface add_density_load
  module procedure add_density_load_tiles
end interface

interface add_particle_load
  module procedure add_particle_load_tiles
end interface

public :: add_density_load, add_particle_load

contains

!-------------------------------------------------------------------------------
! Determines the number of particles per tile based on
! the density profile for the species
!-------------------------------------------------------------------------------
subroutine add_density_load_tiles( this, grid, no_co )
!-------------------------------------------------------------------------------

  implicit none

  ! this must be inout because of the den_value calls
  class( t_species ), intent(inout) :: this
  class( t_grid ), intent(inout) :: grid
  class( t_node_conf ), intent(in) :: no_co

  integer, dimension(p_max_dim) :: g_nx
  integer, dimension(2,p_max_dim) :: my_nx
  integer :: npar_cell

  integer :: i

  real(p_k_part), dimension(p_max_dim) :: dx
  real(p_k_part), dimension(p_max_dim ) :: g_xmin
  real(p_k_part) :: x1, x2

  ! Global tile number in each direction
  integer, dimension( p_max_dim ) :: g_til

  ! Grids/node (l_nx) and remainder of the same (mod_res), and l_nx1 = l_nx + 1
  ! Patterned after even_dist in os-grid-parallel
  integer, dimension(p_max_dim) :: l_nx, l_nx1, mod_res

  ! offset past region where we need an extra cell?
  ! 0 if cell is in partition region with an extra cell (partitions have l_nx1)
  ! 1 if cell is past partition region with an extra cell (partitions have l_nx)
  integer, dimension(p_max_dim) :: ext_off

  integer :: i1, i2, i3
  integer :: j1, j2, j3
  integer :: til

  ! particle positions (global / inside cell )
  real(p_k_part), dimension(:,:), pointer :: ppos, ppos_cell

  ! particle charges
  real(p_k_part), dimension(:), pointer :: pcharge

  ! increment in position for the particles
  real(p_k_part), dimension(3) :: dxpart

  real( p_double ) :: cshift

  class(t_psource), pointer :: profile


  profile => this % source

  ! This will only work if using a t_profile_std (or subclass) profile
  select type (profile)

  class is (t_psource_std)

  select type(no_co); class is(t_node_conf_tiles)
    g_til(1:p_x_dim) = no_co % ti_co % g_til
  end select

  ! get grid sizes and distance between particles
  do i = 1, p_x_dim
    g_nx(i) = grid%g_nx(i)
    my_nx(p_lower:p_upper, i) = grid%my_nx(p_lower:p_upper, i)
    g_xmin( i ) = grid%g_box( p_lower, i )

    ! the species%dx variable has not been defined yet
    dx(i) = ( grid%g_box( p_upper, i ) - grid%g_box( p_lower, i ) ) / g_nx(i)

    dxpart(i) = 1.0_p_k_part / this%num_par_x(i)

    l_nx(i) = g_nx(i)/g_til(i)
    mod_res(i) = mod(g_nx(i),g_til(i))
  enddo

  l_nx1 = l_nx + 1

  ! get total number of particles per cell
  npar_cell = this%num_par_x(1)
  do i = 2, p_x_dim
    npar_cell = npar_cell*this%num_par_x(i)
  enddo

  ! initialize temp buffers
  call alloc( ppos, (/ p_x_dim, npar_cell /) )
  call alloc( ppos_cell, (/ p_x_dim, npar_cell /) )
  call alloc( pcharge, (/ npar_cell /))


  if ( this%pos_type == p_cell_near ) then

    ! positions inside cell are from -0.5 to +0.5
    cshift = 0.5_p_double

    if ( this%coordinates /= p_cylindrical_b ) then

      ! The particle global position is considered to be calculated from the center of the cell, so
      ! that a particle with gix = 1, x = -0.5 will be at the left edge of the simulation

      ! To simplify absolute position calculations we can just shift the global min by half a cell
      do i = 1, p_x_dim
        g_xmin(i) = g_xmin(i) + 0.5_p_double * this%dx(i)
      enddo

      ! all quantities depending on global positions must (should, it's just half a cell) be
      ! adjusted

    else
      ! In cylindrical coordinates we cannot do this for the radial direction, because the algorithm
      ! puts the radial boundary at position r = 0, which is placed at the center of cell 1.

      ! In this case, for even interpolation, particles with gix = 2, x = -0.5 will be in the axis,
      ! and for odd interpolation, this would be gix = 1, x = +0.5

      g_xmin(1) = g_xmin(1) + 0.5_p_double * this%dx(1)
    endif

  else

    cshift = 0.0_p_double
  endif

  select type(grid); class is(t_grid_group)
  select type(no_co); class is(t_node_conf_tiles)

  ! count number of particles to inject in each direction
  select case( p_x_dim )

  case (1)
    i = 0
    do i1 = 1, this%num_par_x(1)
      i = i + 1
      ppos_cell(1,i) = real( dxpart(1) * ( i1 - 0.5_p_double ) - cshift, p_k_part )
    enddo

    do i1 = my_nx(p_lower,1), my_nx(p_upper,1)
      ppos(1,:) = (ppos_cell(1,:) + (i1-1))*dx(1) + g_xmin(1)
      ! Are we past region with 1 extra cell? 0 if not, 1 if yes
      if(mod_res(1)==0) then
        ext_off(1) = 1
      else
        ext_off(1) = min( 1, (i1-1)/(l_nx1(1)*mod_res(1)) )
      endif
      ! Tile number in x1
      j1 = ext_off(1)*mod_res(1) + &
      ( i1 - ext_off(1) * l_nx1(1) * mod_res(1) - 1 )/( l_nx1(1) - ext_off(1) ) + 1
      call profile % get_den_value( ppos, npar_cell, pcharge )
      til = no_co%ti_co%hi1( j1 )

      do i = 1, npar_cell
        if (pcharge(i) >= profile%den_min) then
          grid%int_load(til) = grid%int_load(til) + 1
        endif
      enddo
    enddo

  case (2)
    i = 0
    do i1 = 1, this%num_par_x(1)
      x1 = real( dxpart(1) * ( i1 - 0.5_p_double ) - cshift, p_k_part )
      do i2 = 1, this%num_par_x(2)
        i = i + 1
        ppos_cell(1,i) = x1
        ppos_cell(2,i) = real( dxpart(2) * ( i2 - 0.5_p_double ) - cshift, p_k_part )
      enddo
    enddo

    do i1 = my_nx(p_lower,1), my_nx(p_upper,1)
      ppos(1,:) = (ppos_cell(1,:) + (i1-1))*dx(1) + g_xmin(1)
      ! Are we past region with 1 extra cell? 0 if not, 1 if yes
      if(mod_res(1)==0) then
        ext_off(1) = 1
      else
        ext_off(1) = min( 1, (i1-1)/(l_nx1(1)*mod_res(1)) )
      endif
      ! Tile number in x1
      j1 = ext_off(1)*mod_res(1) + &
      ( i1 - ext_off(1) * l_nx1(1) * mod_res(1) - 1 )/( l_nx1(1) - ext_off(1) ) + 1

      do i2 = my_nx(p_lower,2), my_nx(p_upper,2)
        ppos(2,:) = (ppos_cell(2,:) + (i2-1))*dx(2) + g_xmin(2)
        ! Are we past region with 1 extra cell? 0 if not, 1 if yes
        if(mod_res(2)==0) then
          ext_off(2) = 1
        else
          ext_off(2) = min( 1, (i2-1)/(l_nx1(2)*mod_res(2)) )
        endif
        ! Tile number in x2
        j2 = ext_off(2)*mod_res(2) + &
        ( i2 - ext_off(2) * l_nx1(2) * mod_res(2) - 1 )/( l_nx1(2) - ext_off(2) ) + 1

        call profile % get_den_value( ppos, npar_cell, pcharge )

        select case( no_co%ti_co%topology_type )
        case(p_topology_viktor)

          do i = 1, npar_cell
            if (pcharge(i) >= profile%den_min) then
              grid%til_load_2(j1,j2) = grid%til_load_2(j1,j2) + 1
            endif
          enddo

        case default

          til = no_co%ti_co%hi2( j1, j2 )

          do i = 1, npar_cell
            if (pcharge(i) >= profile%den_min) then
              grid%int_load(til) = grid%int_load(til) + 1
            endif
          enddo

        end select

      enddo

    enddo

  case (3)
    i = 0
    do i1 = 1, this%num_par_x(1)
      x1 = real( dxpart(1) * ( i1 - 0.5_p_double) - cshift, p_k_part )
      do i2 = 1, this%num_par_x(2)
        x2 = real( dxpart(2) * ( i2 - 0.5_p_double ) - cshift, p_k_part )
        do i3 = 1, this%num_par_x(3)
          i = i + 1
          ppos_cell(1,i) = x1
          ppos_cell(2,i) = x2
          ppos_cell(3,i) = real( dxpart(3) * ( i3 - 0.5_p_double ) - cshift, p_k_part )
        enddo
      enddo
    enddo


    do i1 = my_nx(p_lower,1), my_nx(p_upper,1)
      ppos(1,:) = (ppos_cell(1,:) + (i1-1))*dx(1) + g_xmin(1)
      ! Are we past region with 1 extra cell? 0 if not, 1 if yes
      if(mod_res(1)==0) then
        ext_off(1) = 1
      else
        ext_off(1) = min( 1, (i1-1)/(l_nx1(1)*mod_res(1)) )
      endif
      ! Tile number in x1
      j1 = ext_off(1)*mod_res(1) + &
      ( i1 - ext_off(1) * l_nx1(1) * mod_res(1) - 1 )/( l_nx1(1) - ext_off(1) ) + 1

      do i2 = my_nx(p_lower,2), my_nx(p_upper,2)
        ppos(2,:) = (ppos_cell(2,:) + (i2-1))*dx(2) + g_xmin(2)
        ! Are we past region with 1 extra cell? 0 if not, 1 if yes
        if(mod_res(2)==0) then
          ext_off(2) = 1
        else
          ext_off(2) = min( 1, (i2-1)/(l_nx1(2)*mod_res(2)) )
        endif
        ! Tile number in x2
        j2 = ext_off(2)*mod_res(2) + &
        ( i2 - ext_off(2) * l_nx1(2) * mod_res(2) - 1 )/( l_nx1(2) - ext_off(2) ) + 1

        do i3 = my_nx(p_lower,3), my_nx(p_upper,3)
          ppos(3,:) = (ppos_cell(3,:) + (i3-1))*dx(3) + g_xmin(3)
          ! Are we past region with 1 extra cell? 0 if not, 1 if yes
          if(mod_res(3)==0) then
            ext_off(3) = 1
          else
            ext_off(3) = min( 1, (i3-1)/(l_nx1(3)*mod_res(3)) )
          endif
          ! Tile number in x3
          j3 = ext_off(3)*mod_res(3) + &
          ( i3 - ext_off(3) * l_nx1(3) * mod_res(3) - 1 )/( l_nx1(3) - ext_off(3) ) + 1

          call profile % get_den_value( ppos, npar_cell, pcharge )

          select case( no_co%ti_co%topology_type )
          case(p_topology_viktor)

            do i = 1, npar_cell
              if (pcharge(i) >= profile%den_min) then
                grid%til_load_3(j1,j2,j3) = grid%til_load_3(j1,j2,j3) + 1
              endif
            enddo

          case default

            til = no_co%ti_co%hi3( j1, j2, j3 )

            do i = 1, npar_cell
              if (pcharge(i) >= profile%den_min) then
                grid%int_load(til) = grid%int_load(til) + 1
              endif
            enddo

          end select

        enddo

      enddo

    enddo

  end select

  end select
  end select

  call freemem( ppos )
  call freemem( ppos_cell )
  call freemem( pcharge )

  class default

    if ( mpi_node() == 0 ) then
      write(0,*) "(*warning*) add_density_load_tiles is not supported by the selected &
                 &species initialization type"
    endif

  end select

end subroutine add_density_load_tiles
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Determines the number of particles per tile
!-------------------------------------------------------------------------------
subroutine add_particle_load_tiles( this, no_co, int_load, til_load_2, til_load_3 )

  implicit none

  class( t_species ), intent(in) :: this
  class( t_node_conf ), intent(in) :: no_co
  real( p_double ), dimension(:), pointer, optional :: int_load
  real( p_double ), dimension(:,:), pointer, optional :: til_load_2
  real( p_double ), dimension(:,:,:), pointer, optional :: til_load_3

  select type(no_co); class is(t_node_conf_tiles)

  select case( no_co%ti_co%topology_type )
  case(p_topology_viktor)

    select case(p_x_dim)
    case(2)
      til_load_2(no_co%ti_co%g_tgp(1),no_co%ti_co%g_tgp(2)) = &
        til_load_2(no_co%ti_co%g_tgp(1),no_co%ti_co%g_tgp(2)) + this%num_par
    case(3)
      til_load_3(no_co%ti_co%g_tgp(1),no_co%ti_co%g_tgp(2),no_co%ti_co%g_tgp(3)) = &
        til_load_3(no_co%ti_co%g_tgp(1),no_co%ti_co%g_tgp(2),no_co%ti_co%g_tgp(3)) + this%num_par
    end select

  case default

    int_load(no_co%ti_co%g_til_aid) = int_load(no_co%ti_co%g_til_aid) + this%num_par

  end select

  end select

end subroutine add_particle_load_tiles
!-------------------------------------------------------------------------------


end module m_species_loadbalance_tiles