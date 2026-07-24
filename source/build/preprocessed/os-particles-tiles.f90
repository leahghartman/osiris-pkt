# 1 "tiles/os-particles-tiles.f90"
# 1 "<built-in>" 1
# 1 "<built-in>" 3
# 467 "<built-in>" 3
# 1 "<command line>" 1
# 1 "<built-in>" 2
# 1 "tiles/os-particles-tiles.f90" 2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! particle functions for tile load balance
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# 1 "./os-config.h" 1
! Configuration file for osiris

! ----------------------------------------------------------------------------------------
! Algorithm options
! ----------------------------------------------------------------------------------------

! ----------------------------------------------------------------------------------------
! System options
! ----------------------------------------------------------------------------------------

! MPI supports MPI_IN_PLACE in global operations
!#define 1

! Use OpenMP
!#define 1

! Use SIMD optimized code
!#define SIMD
!#define SIMD_SSE
!#define SIMD_AVX
!#define SIMD_BGQ
!#define SIMD_MIC

! Use SION for checkpointing (not available on all systems)
!#define __RST_IO__ = __RST_SION__

! Use log files
!#define __USE_LOG__

! Use MPE for logging and profiling
!#define __USE_MPE__

! Compiler does not fully support the sizeof intrinsic which is a Fortran 2003 feature
!#define __NO_SIZEOF__

! Use the PAPI library for profiling
!#define __USE_PAPI__


! ----------------------------------------------------------------------------------------
! Distribution Options
! ----------------------------------------------------------------------------------------
!
! Optional modules to be removed for distribution
! These are all turned on by default unless the __DISTRO__ preprocessor macro is defined
!
! ----------------------------------------------------------------------------------------



! Include Ionization module


! Include binary collisions module


! Include particle tracking module


! Include perfectly matched layers boundary conditions for EMF


! Include spin advance module which is disabled by default
! #define __HAS_SPIN__

! Include a debug flag for cylindrical modes simulations
! #define __CYL_MODES_DEBUG__

! Compile functions for reporting B fields in the RaDiO module
! #define __HAS_RAD_BFLD__
# 6 "tiles/os-particles-tiles.f90" 2
# 1 "./os-preprocess.fpp" 1
!
! File: os-preprocess.fpp
!
! A set of preprocessing macros for osiris, and the fpp/cpp preprocessor
!
!



! Macros for the IBM CPP/GNU preprocessor




! Assertion macro



! Debug functions
# 33 "./os-preprocess.fpp"
! SCR hostname functions







! LOG functions




! ERROR functions






! WARNING functions




! functions for restart io
# 7 "tiles/os-particles-tiles.f90" 2

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
