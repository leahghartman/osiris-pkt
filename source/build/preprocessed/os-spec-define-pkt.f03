# 1 "pkt/os-spec-define-pkt.f03"
# 1 "<built-in>" 1
# 1 "<built-in>" 3
# 467 "<built-in>" 3
# 1 "<command line>" 1
# 1 "<built-in>" 2
# 1 "pkt/os-spec-define-pkt.f03" 2
!-----------------------------------------------------------------------------------------

! pkt species definition module

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
# 6 "pkt/os-spec-define-pkt.f03" 2
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
# 7 "pkt/os-spec-define-pkt.f03" 2

module m_species_define_pkt

# 1 "./memory/memory.h" 1
! Include file for the memory module
! Files must include this file using #include "memory/memory.h" rather than just using the module
! the #include must be placed where the module use statement would usually be i.e.
!
! module module2
!
! use module1
! #include "memory/memory.h"
!
# 19 "./memory/memory.h"
use memory
# 11 "pkt/os-spec-define-pkt.f03" 2

use m_system
use m_parameters
use m_space
use m_grid_define
use m_node_conf
use m_restart
use m_time_step
use m_input_file

use m_emf_define
use m_vdf_define

use m_species_define

use m_vdf_comm , only : t_vdf_msg
use m_current_define , only : t_current

implicit none

private

! string to id restart data
character(len=*), parameter :: p_elecposi_rst_id = "elec/posi rst data - 0x0000"

!-----------------------------------------------------------------------------------------
! t_species_pkt
! - species with pkt pusher and merging
!-----------------------------------------------------------------------------------------

type, extends( t_species ) :: t_species_pkt

  real(p_k_part) :: n0 !reference density [cm^-3]
  real(p_k_part) :: omega_p0 !reference electron plasma frequency
  real(p_double) :: pkt_g_cutoff, p_emit_cutoff ! cutoff parameters for pkt emissions

  integer :: ndump_fac_rad

  ! emission flag
  logical :: if_damp_classical
  logical :: if_damp_pkt

  ! radiated energy
  real(p_double), dimension(5) :: energy_rad = 0.0_p_double

  ! pointer to pkt group photons
  class( t_species ), pointer :: photons ! photons, always associated

  ! pointer to reciprocal pkt species
  class( t_species ), pointer :: recip_spec => null()

  ! pair production model
  integer :: pairprod

  real(p_double) :: pkt_zeta
  real(p_double) :: pairprod_gthr
  real(p_double) :: pairprod_gpair

  ! radiation and chi spectrum parameters
  real(p_double) :: radspect_emin, chi_emit_min
  real(p_double) :: radspect_emax, chi_emit_max
  integer :: radspect_bins, chi_emit_nbins

  ! radiation spectrum arrays
  real(p_double), dimension(:), pointer :: radspect => null()
  real(p_double), dimension(:), pointer :: radspect_oor => null()
  real(p_double), dimension(:), pointer :: chispect_photons => null()
  real(p_double), dimension(:), pointer :: chispect_leptons => null()

  ! spherical radiation spectrum parameters
  real(p_double) :: radanglespect_thresh
  integer :: radanglespect_bins

  ! spherical radiation spectrum arrays
  real(p_double), dimension(:,:,:), pointer :: radanglespect => null()

  ! spherical radiation spectrum parameters
  real(p_double) :: radsphere_emin, radsphere_emax
  integer :: radsphere_nbins
  logical :: if_add_classical_radsphere

  ! spherical radiation spectrum arrays
  real(p_double), dimension(:,:,:), pointer :: radsphere => null()

  ! angular radiation detector ndump
  integer :: ndump_fac_raddetector

  ! radiation detector parameters
  real(p_double) :: raddetector_emin, raddetector_emax, raddetector_phimin, raddetector_phimax, &
                    raddetector_thetamin, raddetector_thetamax
  integer :: raddetector_nbins, raddetector_dir

  ! radiation detector array
  real(p_double), dimension(:,:), pointer :: raddetector => null()

  ! load balance parameters, obtained from particles object
  real(p_double) :: dlb_fld2_thresh
  real(p_double) :: dlb_part_weight

  type(t_vdf), pointer :: photon_density_grid => null()
  type(t_vdf), pointer :: photon_gradient_grid => null()

contains
  procedure :: allocate_objs => allocate_objs_spec_pkt
  procedure :: read_input => read_input_spec_pkt
  procedure :: set_dudt => set_dudt_spec_pkt
  procedure :: init => init_spec_pkt
  procedure :: list_algorithm => list_algorithm_spec_pkt
  procedure :: cleanup => cleanup_spec_pkt
  procedure :: reshape => reshape_spec_pkt
  procedure :: add_particle_load => add_particle_load_spec_pkt

end type t_species_pkt

! type bound procedure (method) interfaces

interface
  subroutine set_dudt_spec_pkt( this )

    import t_species_pkt

    class( t_species_pkt ), intent(inout) :: this

  end subroutine
end interface

interface add_particle_load_pkt
  module procedure add_particle_load_spec_pkt_func
end interface

public :: t_species_pkt, add_particle_load_pkt

contains

!-----------------------------------------------------------------------------------------
! Allocates algorithm specific objects
!-----------------------------------------------------------------------------------------
subroutine allocate_objs_spec_pkt( this )

  implicit none

  class( t_species_pkt ), intent(inout) :: this

end subroutine allocate_objs_spec_pkt
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Read pkt species information from input file
!-----------------------------------------------------------------------------------------
subroutine read_input_spec_pkt( this, input_file, def_name, periodic, if_move, grid, &
  dt, read_prof, sim_options )

    class( t_species_pkt ), intent(inout) :: this
    class( t_input_file ), intent(inout) :: input_file

    character(len = *), intent(in) :: def_name
    logical, dimension(:), intent(in) :: periodic, if_move
    class( t_grid ), intent(in) :: grid

    real(p_double), intent(in) :: dt
    logical, intent(in) :: read_prof

    type(t_options), intent(in) :: sim_options

    ! call superclass read_input
    call this % t_species % read_input( input_file, def_name, periodic, if_move, grid, &
                            dt, read_prof, sim_options )

    ! allocate algorithm specific objects
    call this % allocate_objs()

    ! set the omega_p0 and n0 variables
    this % omega_p0 = real( sim_options % omega_p0, p_k_part )
    this % n0 = real( sim_options % n0, p_k_part )

end subroutine read_input_spec_pkt
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Printout the algorithm used by the pusher
!-----------------------------------------------------------------------------------------
subroutine list_algorithm_spec_pkt( this )

  implicit none
  class( t_species_pkt ), intent(in) :: this

  if (mpi_node()==0) print *,"Species : ", trim(this%name)

  ! Push type
  select case ( this%push_type )
  case (p_std)
      print *, '   - pkt standard (Boris) pusher'
  case(p_exact)
      print *, '   - pkt exact pusher'
  case(p_exact_rr)
      print *, '   - pkt exact pusher with radiation reaction'
  case(p_std_feedback)
      print *, '   - pkt standard pusher with photon feedback'
  end select

  if ( this%free_stream ) then
      print *, '- Free streaming particles (no dudt)'
      write(0,*) '- (*warning*) ', trim(this%name), ' are free streaming!'
  endif

end subroutine list_algorithm_spec_pkt
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Cleans up algorithm specific objects
!-----------------------------------------------------------------------------------------
subroutine cleanup_spec_pkt( this )

  implicit none

  class(t_species_pkt) , intent(inout) :: this

  ! cleanup local data structures
  ! call this % emit % cleanup()

  ! cleanup superclass
  call this % t_species % cleanup()

end subroutine cleanup_spec_pkt
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Initializes pkt species
!-----------------------------------------------------------------------------------------
subroutine init_spec_pkt( this, sp_id, interpolation, grid_center, grid, g_space, emf, jay, &
  no_co, send_vdf, recv_vdf, bnd_cross, node_cross, send_spec, &
  recv_spec, ndump_fac, restart, restart_handle, sim_options, tstep, tmin, tmax )
!-----------------------------------------------------------------------------------------
! sets up this data structure from the given information
!-----------------------------------------------------------------------------------------

  implicit none

  class( t_species_pkt ), intent(inout) :: this

  integer, intent(in) :: sp_id
  integer, intent(in) :: interpolation
  logical, intent(in) :: grid_center
  class( t_grid ), intent(in) :: grid
  type( t_space ), intent(in) :: g_space
  class( t_emf ), intent(inout) :: emf
  class( t_current ), intent(inout) :: jay
  class( t_node_conf ), intent(in) :: no_co
  type(t_vdf_msg), dimension(2), intent(inout) :: send_vdf, recv_vdf
  type( t_part_idx ), dimension(2), intent(inout) :: bnd_cross
  type( t_part_idx ), intent(inout) :: node_cross
  type( t_spec_msg ), dimension(2), intent(inout) :: send_spec, recv_spec
  integer, intent(in) :: ndump_fac
  logical, intent(in) :: restart
  type( t_restart_handle ), intent(in) :: restart_handle
  type(t_options), intent(in) :: sim_options
  type( t_time_step ), intent(in) :: tstep
  real(p_double), intent(in) :: tmin, tmax

  ! local variables
  integer, dimension(2, p_x_dim) :: gc_num
  real( p_double ), dimension(p_x_dim) :: dx

  ! call superclass init
  call this % t_species % init( sp_id, interpolation, grid_center, grid, g_space, emf, jay, &
      no_co, send_vdf, recv_vdf, bnd_cross, node_cross, send_spec, &
      recv_spec, ndump_fac, restart, restart_handle, sim_options, tstep, tmin, tmax )

  ! overwrite pusher option
  call this % set_dudt()

end subroutine init_spec_pkt
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! redistribute the species particle data through all nodes when node grids change and
! store new grid boundaries
!-----------------------------------------------------------------------------------------
subroutine reshape_spec_pkt( this, old_grid, new_grid, msg_patt, no_co, send_msg, recv_msg )
!-----------------------------------------------------------------------------------------

  use m_vdf_comm, only : reshape_copy

  implicit none

  class( t_species_pkt ), intent(inout) :: this
  type( t_msg_patt ), intent(in) :: msg_patt
  class( t_grid ), intent(in) :: old_grid, new_grid
  class( t_node_conf ), intent(in) :: no_co
  type(t_vdf_msg), dimension(2), intent(inout) :: send_msg, recv_msg

  call this % t_species % reshape( old_grid, new_grid, msg_patt, no_co, send_msg, recv_msg )

end subroutine reshape_spec_pkt
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Wrapper for load function
!-----------------------------------------------------------------------------------------
subroutine add_particle_load_spec_pkt( this, grid, emf )

  implicit none

  class( t_species_pkt ), intent(in) :: this
  class( t_grid ), intent(inout) :: grid
  class( t_emf ), intent(in) :: emf

  call add_particle_load_pkt( this, grid, emf, this%omega_p0, this%dlb_fld2_thresh, &
                              this%dlb_part_weight )

end subroutine add_particle_load_spec_pkt
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Determines the number of particles per cell along all directions based on
! the particle positions for the species
!-----------------------------------------------------------------------------------------
subroutine add_particle_load_spec_pkt_func( this, grid, emf, omega_p0, dlb_fld2_thresh, &
                                            dlb_part_weight )

  implicit none

  class( t_species ), intent(in) :: this
  class( t_grid ), intent(inout) :: grid
  class( t_emf ), intent(in) :: emf
  real(p_k_part), intent(in) :: omega_p0
  real(p_double), intent(in) :: dlb_fld2_thresh, dlb_part_weight

  integer, dimension( p_max_dim ) :: dim_off, ix
  real( p_double ), dimension(:), pointer :: int_load
  integer :: i, j, np, pp, ptrcur

  real(p_double), parameter :: hbar = 1.054571800d-34 ! Planck constant over 2 pi
  real(p_double), parameter :: cl = 299792458.0d0 ! speed of light in vacuum [m/s]
  real(p_double), parameter :: me = 9.10938356d-31 ! electron rest mass [kg]
  real(p_double) :: norm_schw2, weight, e2, b2

  real(p_k_part), dimension(p_p_dim,p_cache_size) :: bp, ep

  ! Get the offset of the int_load array for each direction
  dim_off(1) = 0
  do i = 2, grid%x_dim
    dim_off(i) = dim_off(i-1) + grid%g_nx(i-1)
  enddo

  ! Add the correction for global cell positions to dim_off
  do i = 1, grid%x_dim
    dim_off(i) = dim_off(i) + ( this%my_nx_p( p_lower, i ) - 1 )
  enddo

  ! just for clarity
  int_load => grid%int_load

  ! normalization for field squared
  norm_schw2 = ( hbar*omega_p0/(me*cl**2.d0) )**2

  do ptrcur = 1, this%num_par, p_cache_size

    ! check if last copy of table and set np
    if( ptrcur + p_cache_size > this%num_par ) then
      np = this%num_par - ptrcur + 1
    else
      np = p_cache_size
    endif

    ! modify bp & ep to include timestep and charge-to-mass ratio
    ! and perform half the electric field acceleration.
    ! Result is stored in UTEMP.

    call this % get_emf( emf, bp, ep, np, ptrcur )

    pp = ptrcur
    do i = 1, np

      do j = 1, p_x_dim
        ix(j) = dim_off(j) + this%ix(j,pp)
      enddo

      e2 = ep(1,i)**2 + ep(2,i)**2 + ep(3,i)**2
      b2 = bp(1,i)**2 + bp(2,i)**2 + bp(3,i)**2

      if ( e2 + b2 > dlb_fld2_thresh ) then
        weight = dlb_part_weight
      else
        weight = 1 + (e2 + b2) * norm_schw2
      endif

      do j = 1, p_x_dim
        int_load(ix(j)) = int_load(ix(j)) + weight
      enddo

      pp = pp + 1

    enddo

  enddo

end subroutine add_particle_load_spec_pkt_func
!-----------------------------------------------------------------------------------------

end module m_species_define_pkt
