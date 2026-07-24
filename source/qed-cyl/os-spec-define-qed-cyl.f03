!-----------------------------------------------------------------------------------------

! QED species definition module

#include "os-preprocess.fpp"
#include "os-config.h"

module m_species_define_qedcyl

#include "memory/memory.h"

use m_system
use m_parameters
use m_species_define
use m_species_define_qed,        only: add_particle_load_qed
use m_particles_cyl_modes,       only: t_particles_cyl_modes
use m_species_cyl_modes_define,  only: t_species_cyl_modes
use m_input_file,                only: t_input_file, get_namelist
use m_vdf_report,                only: new
use m_species_collisions,        only: read_nml
use m_grid_define,               only: t_grid, t_msg_patt
use m_space,                     only: t_space
use m_emf_define,                only: t_emf
use m_node_conf,                 only: t_node_conf
use m_current_define,            only: t_current
use m_bnd,                       only: t_bnd
use m_species,                   only: move_window
use m_vdf_comm,                  only : t_vdf_msg
use m_vdf_define
use m_time_step
use m_restart
use m_emf_cyl_modes,             only : t_emf_cyl_modes
use m_bremsstrahlung_qedcyl,     only : t_bremsstrahlung_info_qedcyl
use m_trident_coul_qedcyl,       only : t_trident_coul_info_qedcyl
use m_cyl_modes,                 only: t_cyl_modes, setup

implicit none

private

!-----------------------------------------------------------------------------------------
! t_species_qed
!  - species with qed pusher and merging
!-----------------------------------------------------------------------------------------

type, extends( t_species_cyl_modes ) :: t_species_qedcyl

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
  class( t_species_cyl_modes ), pointer :: photons ! photons, always associated

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

  ! Bremsstrahlung photon creation
  type( t_bremsstrahlung_info_qedcyl ) :: bremsstrahlung_info

  ! Trident Coulomb pair creation
  type( t_trident_coul_info_qedcyl ) :: trident_coul_info

  ! Diagnostic pour trident_coul pair production
  integer :: ndump_fac_pairs_from_elec
  integer :: num_pairs_ct = 0
  real( p_double ) :: weight_pairs_ct = 0.0_p_double
  real( p_double ) :: energy_pairs_ct = 0.0_p_double

  ! pointer to reciprocal qed species
  class( t_species_cyl_modes ), pointer :: recip_spec => null()

  ! load balance parameters, obtained from particles object
  real(p_double) :: dlb_fld2_thresh
  real(p_double) :: dlb_part_weight

contains

  procedure :: allocate_objs     => allocate_objs_spec_qedcyl
  procedure :: read_input        => read_input_spec_qedcyl
  procedure :: set_dudt          => set_dudt_spec_qedcyl
  procedure :: init              => init_spec_qedcyl
  procedure :: list_algorithm    => list_algorithm_spec_qedcyl
  procedure :: report            => report_spec_qedcyl
  procedure :: report_energy_rad => report_energy_rad_spec_qedcyl
  procedure :: cleanup           => cleanup_spec_qedcyl
  procedure :: report_pairs_ct   => report_pairs_ct_spec_qedcyl
  procedure :: get_energy_pairs_ct  => get_energy_pairs_ct_qedcyl
  procedure :: reshape           => reshape_spec_qedcyl
  procedure :: add_particle_load => add_particle_load_spec_qedcyl

end type t_species_qedcyl

interface 
  subroutine report_spec_qedcyl( this, emf, g_space, grid, no_co, tstep, t, tmin, &
                                 send_msg, recv_msg )
  
    import t_species_qedcyl, t_emf, t_space, t_grid, t_node_conf, t_time_step, p_double, t_vdf_msg

    class( t_species_qedcyl ),  intent(inout) :: this
    
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
  subroutine set_dudt_spec_qedcyl( this  )

    import t_species_qedcyl

    class( t_species_qedcyl ), intent(inout) :: this

  end subroutine
end interface

interface 
  subroutine report_energy_rad_spec_qedcyl( this, no_co, tstep, t, tmin )
  
    import t_species_qedcyl, t_node_conf, t_time_step, p_double

    class( t_species_qedcyl ), intent(in) :: this
    
    class( t_node_conf ),  intent(in) :: no_co
    type( t_time_step ),   intent(in) :: tstep
    real(p_double),        intent(in) :: t, tmin

  end subroutine
end interface

interface 
  subroutine report_pairs_ct_spec_qedcyl( this, no_co, tstep, t, tmin )
  
    import t_species_qedcyl, t_node_conf, t_time_step, p_double

    class( t_species_qedcyl ), intent(in) :: this
    
    class( t_node_conf ),  intent(in) :: no_co
    type( t_time_step ),   intent(in) :: tstep
    real(p_double),        intent(in) :: t, tmin

  end subroutine
end interface

public t_species_qedcyl

! type diad QED
type, extends( t_diag_species ) :: t_diag_species_qedcyl

contains

  !procedure :: avail_report_quants => avail_report_quants_species_qed
  !procedure :: init_report_quants  => init_report_quants_species_qed
  procedure :: read_input          => read_input_diag_spec_qedcyl

end type t_diag_species_qedcyl

contains

!-----------------------------------------------------------------------------------------
! Allocates algorithm specific objects
!-----------------------------------------------------------------------------------------
subroutine allocate_objs_spec_qedcyl( this )

  implicit none

  class( t_species_qedcyl ), intent(inout) :: this

  ! local variables
 
  ! allocate algorithm specific t_diag_species object
  if ( .not. associated( this%diag ) ) then
      allocate( t_diag_species_qedcyl :: this%diag )
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

end subroutine allocate_objs_spec_qedcyl
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Read qed species information from input file
!-----------------------------------------------------------------------------------------
subroutine read_input_spec_qedcyl( this, input_file, def_name, periodic, if_move, grid, &
  dt, read_prof, sim_options )

    class( t_species_qedcyl ),  intent(inout) :: this
    class( t_input_file ), intent(inout) :: input_file

    character(len = *),    intent(in) :: def_name
    logical, dimension(:), intent(in) :: periodic, if_move
    class( t_grid ),               intent(in) :: grid

    real(p_double), intent(in) :: dt
    logical, intent(in) :: read_prof

    type(t_options), intent(in) :: sim_options

    ! call superclass read_input
    call this % t_species_cyl_modes % read_input( input_file, def_name, periodic, if_move, grid, &
                            dt, read_prof, sim_options )

    ! allocate algorithm specific objects
    call this % allocate_objs()

    ! set the omega_p0 and n0 variables
    this % omega_p0 = real( sim_options % omega_p0, p_k_part )
    this % n0       = real( sim_options % n0, p_k_part )

end subroutine read_input_spec_qedcyl
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Initialize QED species specific diagnostics
!-----------------------------------------------------------------------------------------
subroutine read_input_diag_spec_qedcyl( this, input_file, gamma )
  
  use m_species_diagnostics

  class ( t_diag_species_qedcyl ), intent(inout) :: this
  class( t_input_file ), intent(inout) :: input_file
  real(p_double), intent(in) :: gamma

  call this % t_diag_species % read_input ( input_file, gamma )

end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Printout the algorithm used by the pusher
!-----------------------------------------------------------------------------------------
subroutine list_algorithm_spec_qedcyl( this )

  implicit none
  class( t_species_qedcyl ), intent(in) :: this

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

end subroutine list_algorithm_spec_qedcyl
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Cleans up algorithm specific objects
!-----------------------------------------------------------------------------------------
subroutine cleanup_spec_qedcyl( this )

  implicit none

  class(t_species_qedcyl) , intent(inout) :: this

  ! cleanup local data structures
  ! call this % emit % cleanup()

  ! cleanup superclass
  call this % t_species_cyl_modes % cleanup()

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

end subroutine cleanup_spec_qedcyl
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Initializes qed-cyl species
!-----------------------------------------------------------------------------------------
subroutine init_spec_qedcyl( this, sp_id, interpolation, grid_center, grid, g_space, emf, jay, &
  no_co, send_vdf, recv_vdf, bnd_cross, node_cross, send_spec, &
  recv_spec, ndump_fac, restart, restart_handle, sim_options, tstep, tmin, tmax )
!-----------------------------------------------------------------------------------------
!       sets up this data structure from the given information
!-----------------------------------------------------------------------------------------

  implicit none

  class( t_species_qedcyl ), intent(inout) :: this

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
  integer :: n_modes

  n_modes = this % n_cyl_modes

  ! call superclass init
  call this % t_species_cyl_modes % init( sp_id, interpolation, grid_center, grid, g_space, emf, jay, &
      no_co, send_vdf, recv_vdf, bnd_cross, node_cross, send_spec, &
      recv_spec, ndump_fac, restart, restart_handle, sim_options, tstep, tmin, tmax )
  
  ! overwrite pusher option
  call this % set_dudt()

  ! reset radiated energy diagnostic
  this % energy_rad = 0.0_p_double

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
      call setup( this % bremsstrahlung_info % rho_ion_cyl_m, n_modes, this % bremsstrahlung_info % rho_ion )
  endif

  ! setup vdf with ion density for trident coul
  if (this%trident_coul_info%if_trident_coul) then
      ! The number of guard cells required depends on the interpolation level
      gc_num(p_lower,:) = this%interpolation
      gc_num(p_upper,:) = this%interpolation + 1
      dx = 0.
      call this % trident_coul_info % rho_ion % new( p_x_dim, 1, this%my_nx_p(3,1:p_x_dim), gc_num, dx, .true. )
      call setup( this % trident_coul_info % rho_ion_cyl_m, n_modes, this % trident_coul_info % rho_ion )
  endif

end subroutine init_spec_qedcyl
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Return pair energy on local node
!-----------------------------------------------------------------------------------------
function get_energy_pairs_ct_qedcyl( this )

  use m_math, only : pi

  implicit none

  class( t_species_qedcyl ), intent(in) :: this
  real( p_double ) :: get_energy_pairs_ct_qedcyl

  integer :: i
  real( p_double ) :: energy, cell_volume
  
  ! Get energy calculated previously
  energy = this % energy_pairs_ct

  ! Normalize to cell size and charge over mass ratio
  cell_volume = this % dx(1)
  do i = 2, p_x_dim
    cell_volume = cell_volume * this % dx(i)
  enddo
  
  get_energy_pairs_ct_qedcyl = energy * cell_volume * 2.d0 * pi

end function get_energy_pairs_ct_qedcyl
!---------------------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
! redistribute the species particle data through all nodes when node grids change and
! store new grid boundaries
!-----------------------------------------------------------------------------------------
subroutine reshape_spec_qedcyl( this, old_grid, new_grid, msg_patt, no_co, send_msg, recv_msg )
!-----------------------------------------------------------------------------------------

  use m_vdf_comm, only : reshape_copy

  implicit none

  class( t_species_qedcyl ), intent(inout) :: this
  type( t_msg_patt ), intent(in) :: msg_patt
  class( t_grid ), intent(in) :: old_grid, new_grid
  class( t_node_conf ), intent(in) :: no_co
  type(t_vdf_msg), dimension(2), intent(inout) :: send_msg, recv_msg

  call this % t_species_cyl_modes % reshape( old_grid, new_grid, msg_patt, no_co, &
                                             send_msg, recv_msg )

  if ( this%n_emit%f_dim() > 0 ) then
    call reshape_copy( this%n_emit, old_grid, new_grid, no_co, send_msg, recv_msg )
  endif

  if ( this%weight%f_dim() > 0 ) then
    call reshape_copy( this%weight, old_grid, new_grid, no_co, send_msg, recv_msg )
    call reshape_copy( this%weight_chi, old_grid, new_grid, no_co, send_msg, recv_msg )
  endif

  ! TODO: What about these parameters?

  ! real(p_double), dimension(:), pointer :: radspect => null()
  ! real(p_double), dimension(:), pointer :: radspect_oor => null()
  ! real(p_double), dimension(:), pointer :: chispect_photons => null()
  ! real(p_double), dimension(:), pointer :: chispect_leptons => null()
  ! real(p_double), dimension(:,:,:), pointer :: radanglespect => null()
  ! real(p_double), dimension(:,:,:), pointer :: radsphere => null()
  ! real(p_double), dimension(:,:), pointer :: raddetector => null()

end subroutine reshape_spec_qedcyl
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Wrapper for load function
!-----------------------------------------------------------------------------------------
subroutine add_particle_load_spec_qedcyl( this, grid, emf )

  implicit none

  class( t_species_qedcyl ), intent(in) :: this
  class( t_grid ), intent(inout) :: grid
  class( t_emf ), intent(in) :: emf

  call add_particle_load_qed( this, grid, emf, this%omega_p0, this%dlb_fld2_thresh, &
                              this%dlb_part_weight )

end subroutine add_particle_load_spec_qedcyl
!-----------------------------------------------------------------------------------------

end module m_species_define_qedcyl

