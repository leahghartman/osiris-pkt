!-----------------------------------------------------------------------------------------

! QED species definition module

#include "os-preprocess.fpp"
#include "os-config.h"

module m_species_define_qed

#include "memory/memory.h"

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
use m_mergedel

use m_vdf_comm , only : t_vdf_msg
use m_current_define , only : t_current
use m_bremsstrahlung
use m_trident_coul

implicit none

private

! string to id restart data
character(len=*), parameter :: p_elecposi_rst_id = "elec/posi rst data - 0x0000"

!-----------------------------------------------------------------------------------------
! t_species_qed
!  - species with qed pusher and merging
!-----------------------------------------------------------------------------------------

type, extends( t_species ) :: t_species_qed

  real(p_k_part) :: n0 !reference density [cm^-3]
  real(p_k_part) :: omega_p0 !reference electron plasma frequency
  real(p_double) :: qed_g_cutoff, p_emit_cutoff ! cutoff parameters for qed emissions

  integer :: ndump_fac_rad

  ! emission flag
  logical :: if_damp_classical
  logical :: if_damp_qed

  ! radiated energy
  real(p_double), dimension(5) :: energy_rad = 0.0_p_double

  ! pointer to qed group photons
  class( t_species ), pointer :: photons ! photons, always associated

  ! pointer to reciprocal qed species
  class( t_species ), pointer :: recip_spec => null()

  ! pair production model
  integer :: pairprod 

  real(p_double) :: qed_zeta
  real(p_double) :: pairprod_gthr
  real(p_double) :: pairprod_gpair

  ! particle merging / deleting info
  type( t_merge_info ) :: merge_info

  ! emission events diagnostic
  integer :: ndump_fac_n_emit
  type( t_vdf ) :: n_emit

  ! chi parameter diagnostic
  integer :: ndump_fac_chi
  type( t_vdf ) :: weight, weight_chi

  ! radiation spectrum diagnostic ndump
  integer :: ndump_fac_radspect, ndump_fac_chi_emit

  ! radiation  and chi spectrum parameters
  real(p_double) :: radspect_emin, chi_emit_min
  real(p_double) :: radspect_emax, chi_emit_max
  integer :: radspect_bins, chi_emit_nbins

  ! radiation spectrum arrays
  real(p_double), dimension(:), pointer :: radspect => null()  
  real(p_double), dimension(:), pointer :: radspect_oor => null() 
  real(p_double), dimension(:), pointer :: chispect_photons => null() 
  real(p_double), dimension(:), pointer :: chispect_leptons => null() 

  ! spherical radiation spectrum diagnostic ndump
  integer :: ndump_fac_radanglespect

  ! spherical radiation spectrum parameters
  real(p_double) :: radanglespect_thresh
  integer :: radanglespect_bins

  ! spherical radiation spectrum arrays
  real(p_double), dimension(:,:,:), pointer :: radanglespect => null()

  ! spherical radiation spectrum diagnostic ndump
  integer :: ndump_fac_radsphere

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

  ! bremsstrahlung photon emission
  type( t_bremsstrahlung_info ) :: bremsstrahlung_info

  ! Trident Coulomb pair creation
  type( t_trident_coul_info ) :: trident_coul_info

  ! Diagnostic pour trident_coul pair production
  integer :: ndump_fac_pairs_from_elec
  integer :: num_pairs_ct = 0
  real( p_double ) :: weight_pairs_ct = 0.0_p_double
  real( p_double ) :: energy_pairs_ct = 0.0_p_double

  ! load balance parameters, obtained from particles object
  real(p_double) :: dlb_fld2_thresh
  real(p_double) :: dlb_part_weight

contains
  procedure :: allocate_objs     => allocate_objs_spec_qed
  procedure :: read_input        => read_input_spec_qed
  procedure :: set_dudt          => set_dudt_spec_qed
  procedure :: init              => init_spec_qed
  procedure :: list_algorithm    => list_algorithm_spec_qed
  procedure :: report            => report_spec_qed
  procedure :: cleanup           => cleanup_spec_qed
  procedure :: mergedel          => mergedel_spec_qed
  procedure :: report_energy_rad => report_energy_rad_spec_qed
  procedure :: report_pairs_ct   => report_pairs_ct_spec_qed
  procedure :: get_energy_pairs_ct  => get_energy_pairs_ct
  procedure :: reshape           => reshape_spec_qed
  procedure :: add_particle_load => add_particle_load_spec_qed

end type t_species_qed

type, extends( t_diag_species ) :: t_diag_species_qed

contains

  !procedure :: avail_report_quants => avail_report_quants_species_qed
  !procedure :: init_report_quants  => init_report_quants_species_qed
  procedure :: read_input          => read_input_diag_spec_qed

end type t_diag_species_qed

! diagnostic quantities specific to the spherical geometry
character(len=12), dimension(2), parameter :: p_report_quants_species_qed = &
                                            (/ 'emit_n  ', &
                                              'emit_ene' /)

! type bound procedure (method) interfaces

interface 
  subroutine report_spec_qed( this, emf, g_space, grid, no_co, tstep, t, tmin, send_msg, &
                              recv_msg )
  
    import t_species_qed, t_emf, t_space, t_grid, t_node_conf, t_time_step, p_double, t_vdf_msg

    class( t_species_qed ),  intent(inout) :: this
    
    class( t_emf ),     intent(inout) :: emf
    type( t_space ),       intent(in) :: g_space
    class( t_grid ),       intent(in) :: grid
    class( t_node_conf ),  intent(in) :: no_co
    type( t_time_step ),   intent(in) :: tstep
    real(p_double),        intent(in) :: t, tmin
    type(t_vdf_msg), dimension(2), intent(inout) :: send_msg, recv_msg

  end subroutine
end interface

interface
  subroutine set_dudt_spec_qed( this  )

    import t_species_qed

    class( t_species_qed ), intent(inout) :: this

  end subroutine
end interface

interface 
  subroutine report_energy_rad_spec_qed( this, no_co, tstep, t, tmin )
  
    import t_species_qed, t_node_conf, t_time_step, p_double

    class( t_species_qed ), intent(in) :: this
    
    class( t_node_conf ),  intent(in) :: no_co
    type( t_time_step ),   intent(in) :: tstep
    real(p_double),        intent(in) :: t, tmin

  end subroutine
end interface

interface 
  subroutine report_pairs_ct_spec_qed( this, no_co, tstep, t, tmin )
  
    import t_species_qed, t_node_conf, t_time_step, p_double

    class( t_species_qed ), intent(in) :: this
    
    class( t_node_conf ),  intent(in) :: no_co
    type( t_time_step ),   intent(in) :: tstep
    real(p_double),        intent(in) :: t, tmin

  end subroutine
end interface

interface add_particle_load_qed
  module procedure add_particle_load_spec_qed_func
end interface

public :: t_species_qed, add_particle_load_qed

contains

!-----------------------------------------------------------------------------------------
! Allocates algorithm specific objects
!-----------------------------------------------------------------------------------------
subroutine allocate_objs_spec_qed( this )

  implicit none

  class( t_species_qed ), intent(inout) :: this
 
  ! allocate algorithm specific t_diag_species object
  if ( .not. associated( this%diag ) ) then
      allocate( t_diag_species_qed :: this%diag )
  endif

  ! setup radiated energy spectrum diagnostic
  if (this%ndump_fac_radspect > 0) then
      ! allocate arrays
      call alloc( this%radspect, (/ this%radspect_bins /) )
      call alloc( this%radspect_oor, (/ 2 /) )
      
      ! set all values to zero
      this%radspect = 0.0_p_double
      this%radspect_oor = 0
  endif

  ! setup chi spectrum diagnostic
  if (this%ndump_fac_chi_emit > 0) then
      ! allocate arrays
      call alloc( this%chispect_leptons, (/ this%chi_emit_nbins /) )
      call alloc( this%chispect_photons, (/ this%chi_emit_nbins /) )

      ! set all values to zero
      this%chispect_leptons = 0.0_p_double
      this%chispect_photons = 0.0_p_double

  endif

  ! setup a spectrum differential in energy and angle
  if (this%ndump_fac_radanglespect > 0) then
      ! allocate arrays
      call alloc( this%radanglespect, (/ this%radanglespect_bins, this%radanglespect_bins, 6 /) )
          
      ! set all values to zero
      this%radspect = 0.0_p_double
  endif

  ! setup spherical radiated energy spectrum diagnostic
  if (this%ndump_fac_radsphere > 0) then
      ! allocate arrays
      call alloc( this%radsphere, (/ this%radsphere_nbins, this%radsphere_nbins, 6 /) )
          
      ! set all values to zero
      this%radsphere = 0.0_p_double
  endif

    ! setup the detector diagnostic
  if (this%ndump_fac_raddetector > 0) then
      ! allocate arrays
      call alloc( this%raddetector, (/ this%raddetector_nbins, this%raddetector_nbins /) )
          
      ! set all values to zero
      this%raddetector = 0.0_p_double
  endif

  ! no need to allocate superclass objects, this is done on t_species % read_input()

end subroutine allocate_objs_spec_qed
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Read qed species information from input file
!-----------------------------------------------------------------------------------------
subroutine read_input_spec_qed( this, input_file, def_name, periodic, if_move, grid, &
  dt, read_prof, sim_options )

    class( t_species_qed ),  intent(inout) :: this
    class( t_input_file ), intent(inout) :: input_file

    character(len = *),    intent(in) :: def_name
    logical, dimension(:), intent(in) :: periodic, if_move
    class( t_grid ),               intent(in) :: grid

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
    this % n0       = real( sim_options % n0, p_k_part )

    ! read additional merge information
    call this % merge_info % read_input( input_file )

end subroutine read_input_spec_qed
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Printout the algorithm used by the pusher
!-----------------------------------------------------------------------------------------
subroutine list_algorithm_spec_qed( this )

  implicit none
  class( t_species_qed ), intent(in) :: this

  SCR_ROOT("Species : ", trim(this%name))

  ! Push type
  select case ( this%push_type )
  case (p_std)
      print *, '   - QED standard (Boris) pusher'
  case(p_exact)
      print *, '   - QED exact pusher'
  case(p_exact_rr)
      print *, '   - QED exact pusher with radiation reaction'
  end select

  if ( this%free_stream ) then
      print *, '- Free streaming particles (no dudt)'
      write(0,*) '- (*warning*) ', trim(this%name), ' are free streaming!'
  endif

end subroutine list_algorithm_spec_qed
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Cleans up algorithm specific objects
!-----------------------------------------------------------------------------------------
subroutine cleanup_spec_qed( this )

  implicit none

  class(t_species_qed) , intent(inout) :: this

  ! cleanup local data structures
  ! call this % emit % cleanup()

  ! cleanup superclass
  call this % t_species % cleanup()

  ! cleanup radiated energy spectrum diagnostic arrays
  if (this%ndump_fac_radspect > 0) then
    ! allocate arrays
    call freemem( this%radspect )
    call freemem( this%radspect_oor )
  endif

  ! cleanup chi spectrum diagnostic arrays
  if (this%ndump_fac_chi_emit > 0) then
    ! allocate arrays
    call freemem( this%chispect_leptons )
    call freemem( this%chispect_photons )
  endif

  ! cleanup spherical radiated energy spectrum diagnostic arrays
  if (this%ndump_fac_radanglespect > 0) then
    ! allocate arrays
    call freemem( this%radanglespect )
  endif

  ! cleanup spherical radiated energy spectrum diagnostic arrays
  if (this%ndump_fac_radsphere > 0) then
    ! deallocate arrays
    call freemem( this%radsphere )
  endif

  ! cleanup detector diagnostic array
  if (this%ndump_fac_raddetector > 0) then
    ! deallocate arrays
    call freemem( this%raddetector )
  endif

  ! cleanup emission events diagnostic vdf
  if (this%ndump_fac_n_emit > 0) then
    call this % n_emit % cleanup()
  endif

  ! cleanup chi parameter diagnostic vdf
  if (this%ndump_fac_chi > 0) then
    call this % weight_chi % cleanup()
    call this % weight % cleanup()
  endif

end subroutine cleanup_spec_qed
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Initializes qed species
!-----------------------------------------------------------------------------------------
subroutine init_spec_qed( this, sp_id, interpolation, grid_center, grid, g_space, emf, jay, &
  no_co, send_vdf, recv_vdf, bnd_cross, node_cross, send_spec, &
  recv_spec, ndump_fac, restart, restart_handle, sim_options, tstep, tmin, tmax )
!-----------------------------------------------------------------------------------------
!       sets up this data structure from the given information
!-----------------------------------------------------------------------------------------

  implicit none

  class( t_species_qed ), intent(inout) :: this

  integer,     intent(in) :: sp_id
  integer, intent(in) :: interpolation
  logical, intent(in) :: grid_center
  class( t_grid ), intent(in)     :: grid
  type( t_space ),     intent(in) :: g_space
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

  ! reset radiated energy diagnostic
  this % energy_rad = 0.0_p_double

  ! Trident_coul diagnostic
  this%num_pairs_ct = 0
  this%weight_pairs_ct = 0.0
  this%energy_pairs_ct = 0.0

  ! setup emission events diagnostic
  if (this%ndump_fac_n_emit > 0) then
      ! We set guard cells to one so that the dlb communication doesn't error out
      gc_num = 1
      dx = emf%dx()
      call this % n_emit % new( p_x_dim, 1, this%my_nx_p(3,1:p_x_dim), gc_num, dx, .true. )
  endif

  ! setup chi parameter diagnostic
  if (this%ndump_fac_chi > 0) then
      ! We set guard cells to one so that the dlb communication doesn't error out
      gc_num = 1
      dx = emf%dx()
      call this % weight % new( p_x_dim, 1, this%my_nx_p(3,1:p_x_dim), gc_num, dx, .true. )
      call this % weight_chi % new( p_x_dim, 1, this%my_nx_p(3,1:p_x_dim), gc_num, dx, .true. )
  endif

  ! setup vdf with ion density for bremsstrahlung
  if (this%bremsstrahlung_info%if_bremsstrahlung) then
      ! The number of guard cells required depends on the interpolation level
      gc_num(p_lower,:) = this%interpolation
      gc_num(p_upper,:) = this%interpolation + 1
      dx = 0.
      call this % bremsstrahlung_info % rho_ion % new( p_x_dim, 1, this%my_nx_p(3,1:p_x_dim), gc_num, dx, .true. )
  endif

    ! setup vdf with ion density for trident_coul
  if (this%trident_coul_info%if_trident_coul) then
      ! The number of guard cells required depends on the interpolation level
      gc_num(p_lower,:) = this%interpolation
      gc_num(p_upper,:) = this%interpolation + 1
      dx = 0.
      call this % trident_coul_info % rho_ion % new( p_x_dim, 1, this%my_nx_p(3,1:p_x_dim), gc_num, dx, .true. )
  endif

end subroutine init_spec_qed
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Merge species
!-----------------------------------------------------------------------------------------
subroutine mergedel_spec_qed( this, n, nxc, t, no_co )

  use m_system
  use m_node_conf
  use m_mergedel

  implicit none

  class( t_species_qed ), intent(inout) :: this
  integer, intent(in) :: n
  integer,        dimension(:), intent(in) :: nxc
  real(p_double),               intent(in) :: t
  class( t_node_conf ), intent(in) :: no_co

  ! merge charged particles
  call mergedel( this, this%merge_info, n, nxc, t, no_co, p_charge )

end subroutine mergedel_spec_qed
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Initialize QED species specific diagnostics
!-----------------------------------------------------------------------------------------
subroutine read_input_diag_spec_qed( this, input_file, gamma )
  
  use m_species_diagnostics

  class ( t_diag_species_qed ), intent(inout) :: this
  class( t_input_file ), intent(inout) :: input_file
  real(p_double), intent(in) :: gamma

  call this % t_diag_species % read_input ( input_file, gamma )

end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Return pair energy on local node
!-----------------------------------------------------------------------------------------
function get_energy_pairs_ct( this )

  implicit none

  class( t_species_qed ), intent(in) :: this
  real( p_double ) :: get_energy_pairs_ct

  integer :: i
  real( p_double ) :: energy, cell_volume

  energy = this % energy_pairs_ct

  ! Normalize to cell size and charge over mass ratio
  cell_volume = this % dx(1)
  do i = 2, p_x_dim
    cell_volume = cell_volume * this % dx(i)
  enddo
  
  get_energy_pairs_ct = energy * cell_volume

end function get_energy_pairs_ct
!---------------------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
! redistribute the species particle data through all nodes when node grids change and
! store new grid boundaries
!-----------------------------------------------------------------------------------------
subroutine reshape_spec_qed( this, old_grid, new_grid, msg_patt, no_co, send_msg, recv_msg )
!-----------------------------------------------------------------------------------------

  use m_vdf_comm, only : reshape_copy

  implicit none

  class( t_species_qed ), intent(inout) :: this
  type( t_msg_patt ), intent(in) :: msg_patt
  class( t_grid ), intent(in) :: old_grid, new_grid
  class( t_node_conf ), intent(in) :: no_co
  type(t_vdf_msg), dimension(2), intent(inout) :: send_msg, recv_msg

  call this % t_species % reshape( old_grid, new_grid, msg_patt, no_co, send_msg, recv_msg )

  if ( this%n_emit%f_dim() > 0 ) then
    call reshape_copy( this%n_emit, old_grid, new_grid, no_co, send_msg, recv_msg )
  endif

  if ( this%weight%f_dim() > 0 ) then
    call reshape_copy( this%weight, old_grid, new_grid, no_co, send_msg, recv_msg )
    call reshape_copy( this%weight_chi, old_grid, new_grid, no_co, send_msg, recv_msg )
  endif

end subroutine reshape_spec_qed
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Wrapper for load function
!-----------------------------------------------------------------------------------------
subroutine add_particle_load_spec_qed( this, grid, emf )

  implicit none

  class( t_species_qed ), intent(in) :: this
  class( t_grid ), intent(inout) :: grid
  class( t_emf ), intent(in) :: emf

  call add_particle_load_qed( this, grid, emf, this%omega_p0, this%dlb_fld2_thresh, &
                              this%dlb_part_weight )

end subroutine add_particle_load_spec_qed
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Determines the number of particles per cell along all directions based on
! the particle positions for the species
!-----------------------------------------------------------------------------------------
subroutine add_particle_load_spec_qed_func( this, grid, emf, omega_p0, dlb_fld2_thresh, &
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
  real(p_double), parameter :: cl   = 299792458.0d0   ! speed of light in vacuum [m/s]
  real(p_double), parameter :: me   = 9.10938356d-31  ! electron rest mass [kg]
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

end subroutine add_particle_load_spec_qed_func
!-----------------------------------------------------------------------------------------

end module m_species_define_qed
