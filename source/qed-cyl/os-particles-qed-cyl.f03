!
! The current version does not allow for normal species, only QEDCYL groups
!

#include "os-config.h"
#include "os-preprocess.fpp"

module m_particles_qedcyl

#include "memory/memory.h"

use m_system
use m_parameters
use m_species_define
use m_restart
use m_particles_define,          only: t_particles, p_report_quants
use m_particles_cyl_modes,       only: t_particles_cyl_modes
use m_species_cyl_modes_define,  only: t_species_cyl_modes
use m_emf_cyl_modes,             only: t_emf_cyl_modes
use m_input_file,                only: t_input_file, get_namelist, disp_out
use m_vdf_define,                only: t_vdf, p_max_reports_len, p_max_reports, p_n_report_type, &
                                  p_full, p_savg, p_senv, p_line, p_slice, t_vdf_report
use m_vdf_report,                only: new
use m_species_collisions,        only: read_nml
use m_grid_define,               only: t_grid, t_msg_patt
use m_grid_parallel,             only: new_msg_patt
use m_species_define_qedcyl,        only: t_species_qedcyl
use m_photons_define_qedcyl,        only: t_photons_qedcyl
use m_diagnostic_utilities,      only: p_diag_prec
use m_space,                     only: t_space
use m_emf_define,                only: t_emf
use m_node_conf,                 only: t_node_conf, n_threads
use m_current_define,            only: t_current
use m_bnd,                       only: t_bnd
use m_species,                   only: move_window
use m_vdf_comm,                  only: t_vdf_msg, reshape_nocopy, reshape_copy
use m_particles,                 only: validate
use m_particles_charge,          only: reshape
use m_time_step
use m_betheheitler_qedcyl
use m_diagfile
use m_node_conf
use m_diagnostic_utilities
use m_zpulse,                    only: t_zpulse_list
use m_species_loadbalance,       only: add_density_load

implicit none

private

! Group of 2 species + 1 photon species for QEDCYL code
type :: t_qedcyl_spec_group

  class( t_photons_qedcyl ), pointer :: photons
  class( t_species_qedcyl ), pointer :: electrons ! electrons and positrons point to species
  class( t_species_qedcyl ), pointer :: positrons ! in the main particle set

contains

  procedure :: read_input  => read_input_spec_group_qedcyl
  procedure :: report      => report_spec_group_qedcyl

end type t_qedcyl_spec_group

public :: t_qedcyl_spec_group

! particle type, exteds the cyl_modes type and add arguments/methods from QED
type, extends( t_particles_cyl_modes ) :: t_particles_qedcyl

  integer :: num_qed
  type( t_qedcyl_spec_group ), dimension(:), pointer :: qed_group => NULL()

  ! load balance parameters
  real(p_double) :: dlb_fld2_thresh
  real(p_double) :: dlb_part_weight

contains

  procedure :: allocate_objs => allocate_objs_part_qedcyl
  procedure :: read_input => read_input_part_qedcyl
  procedure :: init => init_part_qedcyl
  procedure :: list_algorithm => list_algorithm_part_qedcyl
  procedure :: move_window => move_window_part_qedcyl
  procedure :: report => report_part_qedcyl
  procedure :: update_boundary => update_boundary_part_qedcyl
  procedure :: report_energy => report_energy_part_qedcyl
  procedure :: advance_deposit => advance_deposit_qedcyl
  procedure :: cleanup => cleanup_part_qedcyl
  procedure :: write_checkpoint => write_checkpoint_part_qedcyl
  procedure :: get_diag_buffer_size => get_diag_buffer_size_qedcyl
  procedure, private :: init_bremsstrahlung => init_bremsstrahlung
  procedure, private :: init_betheheitler => init_betheheitler
  procedure, private :: init_trident_coul => init_trident_coul
  procedure :: add_load => add_load_part_qedcyl
  procedure :: reshape => reshape_part_qedcyl

end type t_particles_qedcyl

interface report_node_load_part_qedcyl
module procedure report_node_load_part_qedcyl
end interface

public :: t_particles_qedcyl, report_node_load_part_qedcyl

contains

!-----------------------------------------------------------------------------------------
! read necessary information from input file
!-----------------------------------------------------------------------------------------
subroutine read_input_spec_group_qedcyl( this, input_file, periodic, if_move, grid, dt, sim_options )
!-----------------------------------------------------------------------------------------
    
    implicit none

    class( t_qedcyl_spec_group ), intent(inout) :: this
    class( t_input_file ), intent(inout) :: input_file
    logical, dimension(:), intent(in) :: periodic, if_move
    class( t_grid ), intent(in) :: grid
    real(p_double), intent(in) :: dt
    type( t_options ), intent(in) :: sim_options

    integer :: ierr
    character(len = p_max_spname_len) :: spname
    logical :: if_damp_classical, if_damp_qed, if_pairprod
    real(p_double) :: qed_g_cutoff, p_emit_cutoff, pp_chi_cutoff

    namelist /nl_qed_group/ if_damp_classical, if_damp_qed, if_pairprod, &
                            qed_g_cutoff, p_emit_cutoff, pp_chi_cutoff

    if_damp_classical = .true.
    if_damp_qed = .true.
    if_pairprod = .true.
    qed_g_cutoff = 10.0_p_double
    p_emit_cutoff = 2.0_p_double
    pp_chi_cutoff = 0.1_p_double

    SCR_ROOT("Reading QED species group configuration...")

    call get_namelist( input_file, "nl_qed_group", ierr )
  
    if ( ierr == 0 ) then
      read (input_file%nml_text, nml = nl_qed_group, iostat = ierr)
      if (ierr /= 0) then
        if ( mpi_node() == 0 ) then
          write(0,*) ""
          write(0,*) "   Error reading QED species group information"
          write(0,*) "   aborting..."
        endif
        stop
      endif
    else
       SCR_ROOT("   - QED species group parameters missing, using default settings")
    endif

    call read_nml_qedcyl_group_diag( this, input_file )

    ! read electrons species
    SCR_ROOT("Reading associated electrons configuration...")
    
    call this % electrons % read_input( input_file , spname, periodic, &
        if_move, grid, dt, .true., sim_options )

    ! force electrons and set emission variables
    this % electrons % rqm = -1.0
    this % electrons % if_damp_classical = if_damp_classical
    this % electrons % if_damp_qed = if_damp_qed
    this % electrons % qed_g_cutoff = qed_g_cutoff
    this % electrons % p_emit_cutoff = p_emit_cutoff

    ! read positrons species
    SCR_ROOT("Reading associated positrons configuration...")
    
    call this % positrons % read_input( input_file , spname, periodic, &
        if_move, grid, dt, .true., sim_options )

    ! force positrons and set emission variables
    this % positrons % rqm = +1.0
    this % positrons % if_damp_classical = if_damp_classical
    this % positrons % if_damp_qed = if_damp_qed
    this % positrons % qed_g_cutoff = qed_g_cutoff
    this % positrons % p_emit_cutoff = p_emit_cutoff

    ! read photons species
    SCR_ROOT("Reading associated photons configuration...")
    
    call this % photons % read_input( input_file, spname,  periodic, if_move, grid, &
                                                dt, .true., sim_options )
    
    ! force photons and set pair production flag
    this % photons % rqm = 0.0_p_k_part
    this % photons % if_pairprod = if_pairprod
    this % photons % pp_chi_cutoff = pp_chi_cutoff

end subroutine read_input_spec_group_qedcyl
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! read necessary information from input file
!-----------------------------------------------------------------------------------------
subroutine read_input_part_qedcyl( this, input_file, periodic, if_move, grid, dt, sim_options )
!-----------------------------------------------------------------------------------------
  
    implicit none

    ! dummy variables

    class( t_particles_qedcyl ), intent(inout) :: this
    class( t_input_file ), intent(inout) :: input_file
    logical, dimension(:), intent(in) :: periodic, if_move
    class( t_grid ),               intent(in) :: grid
    real(p_double), intent(in) :: dt
    type( t_options ), intent(in) :: sim_options

    integer  :: i, j, num_qed, num_species
    logical  :: low_jay_roundoff

    character( len = p_max_reports_len ), dimension( p_max_reports ) :: reports
    integer  :: ndump_fac, ndump_fac_ave, ndump_fac_lineout, ndump_fac_ene
    integer, dimension(p_x_dim) :: n_ave
    integer                     :: prec, n_tavg

    character(len=20) :: interpolation
    logical :: grid_center

    real(p_double) :: dlb_fld2_thresh, dlb_part_weight

    namelist /nl_particles/ num_qed, num_species, low_jay_roundoff, &
                            ndump_fac, ndump_fac_ave, ndump_fac_lineout, ndump_fac_ene, &
                            n_ave, prec, reports, n_tavg, interpolation, grid_center, &
                            dlb_fld2_thresh, dlb_part_weight

    integer, dimension( p_n_report_type ) :: ndump_fac_all
    integer :: ierr

    character(len = p_max_spname_len) :: spname

    class( t_species ), pointer :: species

    ! executable statements

    ! number of standard species
    num_species = 0
    
    ! number of QED groups (electrons + positrons + photons)
    num_qed = 0

    low_jay_roundoff = .false.

    ! Diagnostics
    reports = "-"
    ndump_fac = 0
    ndump_fac_ave = 0
    ndump_fac_lineout = 0
    ndump_fac_ene = 0

    n_ave = -1
    n_tavg = -1
    prec = p_diag_prec

    interpolation = "quadratic"
    grid_center = .false.

    dlb_fld2_thresh = 1.0d6
    dlb_part_weight = 2.0

    ! Get namelist text from input file
    call get_namelist( input_file, "nl_particles", ierr )

    if ( ierr /= 0 ) then
      if (ierr < 0) then
        print *, "Error reading particles parameters"
      else
        print *, "Error: particles parameters missing"
      endif
      print *, "aborting..."
      stop
    endif

    read (input_file%nml_text, nml = nl_particles, iostat = ierr)
    if (ierr /= 0) then
      print *, "Error reading particles parameters"
      print *, "aborting..."
      stop
    endif

    ! process interpolation scheme
    select case ( trim( interpolation ))
      ! case ( "ngp" ) ! not implemented
      !  this%interpolation = p_ngp
    case ( "linear" )
      this%interpolation = p_linear
    case ( "quadratic" )
      this%interpolation = p_quadratic
    case ( "cubic" )
      this%interpolation = p_cubic
    ! case ( "quartic" )
    !   this%interpolation = p_quartic

    case default
      print *, '   Error reading species parameters'
      print *, '   invalid interpolation: "', trim( interpolation ), '"'
      print *, '   Valid values are "linear", "quadratic" and "cubic".'
      print *, '   aborting...'
      stop
    end select

    this%grid_center = grid_center

    if ( this%grid_center ) then
      if ( mpi_node() == 0 ) then
        write(0,*) "   Error reading particles parameters"
        write(0,*) "   grid_center must be set to false when using cyl_modes"
        write(0,*) "   aborting..."
      endif
      stop
    endif

    ! process position type
    select case ( this%interpolation )
    case ( p_linear, p_cubic )
      this%pos_type = p_cell_low
    case ( p_quadratic, p_quartic )
      this%pos_type = p_cell_near
    end select

    this%low_jay_roundoff = low_jay_roundoff


    ! Particle kinetic energy diagnostics
    this%ndump_fac_ene = ndump_fac_ene

    ! global charge diagnostics
    ndump_fac_all(p_full)  = ndump_fac
    ndump_fac_all(p_savg)  = ndump_fac_ave
    ndump_fac_all(p_senv)  = ndump_fac_ave
    ndump_fac_all(p_line)  = ndump_fac_lineout
    ndump_fac_all(p_slice) = ndump_fac_lineout

    ! process normal reports
    call new( this%reports, reports, p_report_quants, &
      ndump_fac_all, n_ave, n_tavg, prec, &
      p_x_dim, ierr )
    if ( ierr /= 0 ) then
      print *, "(*error*) Invalid report"
      print *, "(*error*) aborting..."
      stop
    endif


    ! QED specific species and groups

    if (sim_options%omega_p0 <= 0.0) then
      if ( mpi_node() == 0 ) then
        print *, "(*error*) When using QED algorithm omega_p0 or n0 must be set in the simulation"
        print *, "(*error*) section at the beggining of the input file"
        print *, "(*error*) bailing out..."
      endif
      stop
    endif
    
    ! no cathodes allowed
    this%num_cathode = 0

    ! no ionization allowed
    this%num_neutral = 0

    ! number of QED groups
    this%num_qed = num_qed

    ! total charged particles species
    this%num_species = 2 * num_qed + num_species

    ! add dlb info
    this % dlb_fld2_thresh = dlb_fld2_thresh
    this % dlb_part_weight = dlb_part_weight

    call this % allocate_objs()

    ! QED Species
    if ( num_qed > 0 ) then
        do i=1, num_qed

            ! read input of each qed group
            call this % qed_group(i) % read_input( input_file, periodic, if_move, grid, &
                                        dt, sim_options )

        end do

        ! ----- Read qed collisions info
        ! SCR_ROOT("   - reading associated collisions configuration...")
        ! call this % qed_collide_info % read_input(input_file,this%num_qed)

        species => this % qed_group(num_qed) % positrons % next
    else
      species => this % species
    endif


    if ( num_species > 0 ) then
        do i=1, num_species
            write(spname, '(A,I0)') 'species ',i
            if ( mpi_node() == 0 .and. disp_out(input_file) ) then
                print '(A,I0,A)', " - standard species (",i,") configuration..."
            endif
            call species % read_input( input_file, spname, periodic, if_move, grid, &
                dt, .true., sim_options )
            species => species % next
        enddo  
    endif

    ! Reads Bremsstrahlung bloc in the input deck
    call this % init_bremsstrahlung(input_file)

    ! Reads Bethe Heitler bloc in the input deck
    call this % init_betheheitler(input_file)

    ! Reads trident_coul bloc in the input deck
    call this % init_trident_coul(input_file)

    ! validate species names
    call this % validate_names()

  #ifdef __HAS_COLLISIONS__
  ! Read collision data
    call read_nml( this%coll, input_file, this%species, this%num_species, sim_options )
  #endif

end subroutine read_input_part_qedcyl
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Overrides allocate_objs method and allocates t_species_qed objects instead
!-----------------------------------------------------------------------------------------
subroutine allocate_objs_part_qedcyl( this )

  implicit none

  class( t_particles_qedcyl ), intent(inout) :: this

  integer :: i
  class(t_species_cyl_modes), pointer :: spec, tail

  ! QED groups species
  tail => null()
  if ( this % num_qed > 0 ) then
    
    allocate( this%qed_group( this%num_qed ) )

    ! allocate QED electron and positron species
    do i = 1, this % num_qed
      
      ! allocate electrons and positrons species
      allocate( t_photons_qedcyl :: this % qed_group(i) % photons   )
      allocate( t_species_qedcyl :: this % qed_group(i) % electrons )
      allocate( t_species_qedcyl :: this % qed_group(i) % positrons )

      ! associate photons pointer in electrons and positrons species with qed group photons
      this % qed_group(i) % electrons % photons => this % qed_group(i) % photons
      this % qed_group(i) % positrons % photons => this % qed_group(i) % photons

      ! associate reciprocal species pointer in electrons and positrons species
      this % qed_group(i) % electrons % recip_spec => this % qed_group(i) % positrons
      this % qed_group(i) % positrons % recip_spec => this % qed_group(i) % electrons

      if ( associated(tail) ) tail % next => this % qed_group(i) % electrons
      this % qed_group(i) % electrons % next => this % qed_group(i) % positrons
      tail => this % qed_group(i) % positrons

      ! pass in dlb parameters
      this % qed_group(i) % photons   % dlb_fld2_thresh = this % dlb_fld2_thresh
      this % qed_group(i) % electrons % dlb_fld2_thresh = this % dlb_fld2_thresh
      this % qed_group(i) % positrons % dlb_fld2_thresh = this % dlb_fld2_thresh
      this % qed_group(i) % photons   % dlb_part_weight = this % dlb_part_weight
      this % qed_group(i) % electrons % dlb_part_weight = this % dlb_part_weight
      this % qed_group(i) % positrons % dlb_part_weight = this % dlb_part_weight
    enddo

    ! associate pointers with global particle list
    this % species => this % qed_group(1) % electrons
  endif

  ! Standard species
  do i = 2 * this % num_qed + 1, this % num_species
    ! Allocate standard species
    allocate( t_species_cyl_modes :: spec )
    if ( associated(tail) ) then
      tail % next => spec
    else
      this % species => spec
    endif
    tail => spec
  enddo

  ! currently cathodes and ionization are not supported
  if ( this%num_cathode > 0 .or. this%num_neutral > 0 ) then
    SCR_ROOT('(*warning*) Cathodes and Ionization are currently disabled in QED mode')
  endif


end subroutine allocate_objs_part_qedcyl
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
! read necessary information from input file
!-----------------------------------------------------------------------------------------
subroutine read_nml_qedcyl_group_diag( this, input_file )
!-----------------------------------------------------------------------------------------

    implicit none

    class( t_qedcyl_spec_group ), intent(inout) :: this
    class( t_input_file ), intent(inout) :: input_file

    integer :: ierr, ndump_fac_rad, ndump_fac_pairs, ndump_fac_radspect, radspect_bins,     &
              ndump_fac_radanglespect, radanglespect_bins, ndump_fac_n_emit, ndump_fac_chi, &
              chi_emit_nbins, ndump_fac_chi_emit, ndump_fac_radsphere, radsphere_nbins,     &
              ndump_fac_raddetector, raddetector_nbins, raddetector_dir, ndump_fac_pairs_from_elec

    logical :: if_add_classical_radsphere

    real(p_double) :: radspect_emin, radspect_emax, radanglespect_thresh, chi_emit_min,     &
              chi_emit_max, raddetector_emin, raddetector_emax, raddetector_phimax,         &
              raddetector_phimin, raddetector_thetamax, raddetector_thetamin,               &
              radsphere_emin, radsphere_emax

    character(len = 11) :: detector_direction

    namelist /nl_qed_group_diag/ ndump_fac_rad, ndump_fac_pairs,                        &
                  ndump_fac_radspect, radspect_emin, radspect_emax, radspect_bins,      &
                  ndump_fac_radanglespect, radanglespect_bins, radanglespect_thresh,    &
                  ndump_fac_n_emit, ndump_fac_chi,                                      &
                  ndump_fac_chi_emit, chi_emit_nbins, chi_emit_min, chi_emit_max,       &
                  ndump_fac_radsphere, radsphere_nbins, radsphere_emin, radsphere_emax, if_add_classical_radsphere,      &
                  ndump_fac_raddetector, raddetector_nbins, raddetector_emin, raddetector_emax,                          &
                  raddetector_phimax, raddetector_phimin, raddetector_thetamax, raddetector_thetamin, detector_direction,&
                  ndump_fac_pairs_from_elec

    ndump_fac_rad = 0
    ndump_fac_pairs = 0

    ndump_fac_radspect = 0
    radspect_emin = 1.957e-5
    radspect_emax = 1.957e5
    radspect_bins = 400

    ndump_fac_radanglespect = 0
    radanglespect_bins = 100
    radanglespect_thresh = 1.0

    ndump_fac_chi_emit = 0 
    chi_emit_min = 0.0001 
    chi_emit_max = 100.0 
    chi_emit_nbins = 240 

    ndump_fac_radsphere = 0
    radsphere_nbins = 100
    if_add_classical_radsphere = .false.
    radsphere_emin = 0.0
    radsphere_emax = 1.d15
    
    detector_direction = "x1_positive"
    ndump_fac_raddetector = 0
    raddetector_dir = 1
    raddetector_nbins = 100
    raddetector_emin = 0
    raddetector_emax = 1.d15
    raddetector_phimin = - 0.05
    raddetector_phimax = 0.05
    raddetector_thetamin = - 0.05
    raddetector_thetamax = 0.05

    ndump_fac_n_emit = 0
    ndump_fac_chi = 0

    ndump_fac_pairs_from_elec = 0

    call get_namelist( input_file, "nl_qed_group_diag", ierr )
  
    if ( ierr == 0 ) then
      read (input_file%nml_text, nml = nl_qed_group_diag, iostat = ierr)
      if (ierr /= 0) then
        if ( mpi_node() == 0 ) then
          write(0,*) ""
          write(0,*) "   Error reading QED species group diagnostics information"
          write(0,*) "   aborting..."
        endif
        stop
      endif
    else
       SCR_ROOT("   - QED species group diagnostics parameters missing, using default settings")
    endif

    select case(trim(detector_direction))
      case( "x1_positive" )
        raddetector_dir=1
      case( "x1_negative" )
        raddetector_dir=2
      case( "x2_positive" )
        raddetector_dir=3
      case( "x2_negative" )
        raddetector_dir=4
      case( "x3_positive" )
        raddetector_dir=5
      case( "x3_negative" )
        raddetector_dir=6
      case default
        print *, "Error reading QED diag parameters"
        print *, "Invalid value for the detector_dir parameter"
        print *, "Available detector directions are 'x1_positive', 'x1_negative'," 
        print *,"'x2_positive', 'x2_negative', 'x3_positive', 'x3_negative' "
        print *, "aborting..."
        stop
    end select 

    ! radiated energy diagnostic
    this % electrons % ndump_fac_rad = ndump_fac_rad
    this % positrons % ndump_fac_rad = ndump_fac_rad

    ! todo: validation of energy spectrum diagnostic parameter goes here

    ! radiated energy spectrum diagnostic
    this % electrons % ndump_fac_radspect = ndump_fac_radspect
    this % positrons % ndump_fac_radspect = ndump_fac_radspect

    this % electrons % radspect_emin = radspect_emin
    this % positrons % radspect_emin = radspect_emin
    this % electrons % radspect_emax = radspect_emax
    this % positrons % radspect_emax = radspect_emax
    this % electrons % radspect_bins = radspect_bins
    this % positrons % radspect_bins = radspect_bins

    ! chi parameter spectrum diagnostic
    this % electrons % ndump_fac_chi_emit = ndump_fac_chi_emit
    this % positrons % ndump_fac_chi_emit = ndump_fac_chi_emit

    this % electrons % chi_emit_min = chi_emit_min
    this % positrons % chi_emit_min = chi_emit_min
    this % electrons % chi_emit_max = chi_emit_max
    this % positrons % chi_emit_max = chi_emit_max
    this % electrons % chi_emit_nbins = chi_emit_nbins
    this % positrons % chi_emit_nbins = chi_emit_nbins

    ! spherical radiated energy spectrum diagnostic
    this % electrons % ndump_fac_radanglespect = ndump_fac_radanglespect
    this % positrons % ndump_fac_radanglespect = ndump_fac_radanglespect

    this % electrons % radanglespect_bins = radanglespect_bins
    this % positrons % radanglespect_bins = radanglespect_bins
    this % electrons % radanglespect_thresh = radanglespect_thresh
    this % positrons % radanglespect_thresh = radanglespect_thresh

    ! emitted pairs energy diagnostic
    this % photons % ndump_fac_pairs = ndump_fac_pairs

    ! emission events diagnostic
    this % electrons % ndump_fac_n_emit = ndump_fac_n_emit
    this % positrons % ndump_fac_n_emit = ndump_fac_n_emit
    this % photons % ndump_fac_n_emit = ndump_fac_n_emit

    ! chi parameter grid diagnostic
    this % electrons % ndump_fac_chi = ndump_fac_chi
    this % positrons % ndump_fac_chi = ndump_fac_chi
    this % photons % ndump_fac_chi = ndump_fac_chi

    ! radiation sphere diagnostic
    this % electrons % ndump_fac_radsphere = ndump_fac_radsphere
    this % positrons % ndump_fac_radsphere = ndump_fac_radsphere

    this % electrons % radsphere_nbins = radsphere_nbins
    this % positrons % radsphere_nbins = radsphere_nbins
    this % electrons % radsphere_emin = radsphere_emin
    this % positrons % radsphere_emin = radsphere_emin
    this % electrons % radsphere_emax = radsphere_emax
    this % positrons % radsphere_emax = radsphere_emax
    this % electrons % if_add_classical_radsphere = if_add_classical_radsphere
    this % positrons % if_add_classical_radsphere = if_add_classical_radsphere

    ! radiation detector 
    this % electrons % ndump_fac_raddetector = ndump_fac_raddetector
    this % positrons % ndump_fac_raddetector = ndump_fac_raddetector
    this % electrons % raddetector_nbins = raddetector_nbins
    this % positrons % raddetector_nbins = raddetector_nbins

    this % electrons % raddetector_dir = raddetector_dir
    this % positrons % raddetector_dir = raddetector_dir

    this % electrons % raddetector_emin = raddetector_emin
    this % positrons % raddetector_emin = raddetector_emin
    this % electrons % raddetector_emax = raddetector_emax
    this % positrons % raddetector_emax = raddetector_emax

    this % electrons % raddetector_phimin = raddetector_phimin
    this % positrons % raddetector_phimin = raddetector_phimin
    this % electrons % raddetector_phimax = raddetector_phimax
    this % positrons % raddetector_phimax = raddetector_phimax

    this % electrons % raddetector_thetamin = raddetector_thetamin
    this % positrons % raddetector_thetamin = raddetector_thetamin
    this % electrons % raddetector_thetamax = raddetector_thetamax
    this % positrons % raddetector_thetamax = raddetector_thetamax

    this % electrons % ndump_fac_pairs_from_elec = ndump_fac_pairs_from_elec
    this % positrons % ndump_fac_pairs_from_elec = ndump_fac_pairs_from_elec

end subroutine read_nml_qedcyl_group_diag
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! sets up this data structure from the given information
!-----------------------------------------------------------------------------------------
subroutine init_part_qedcyl( this, g_space, jay, emf, grid, no_co, bnd, zpulse_list, &
                             ndump_fac, restart, restart_handle, t, tstep, tmin, tmax, &
                             sim_options )
!-----------------------------------------------------------------------------------------

  implicit none

  class( t_particles_qedcyl ), intent(inout) :: this

  type( t_space ),     intent(in) :: g_space
  class( t_emf) ,intent(inout) :: emf
  class( t_current ) ,intent(inout) :: jay
  class( t_grid ), intent(in) :: grid
  class( t_node_conf ), intent(in) :: no_co
  class(t_bnd), intent(inout) :: bnd
  class( t_zpulse_list ), intent(in) :: zpulse_list
  integer, intent(in) :: ndump_fac
  logical, intent(in) :: restart
  type( t_restart_handle ), intent(in) :: restart_handle
  real(p_double), intent(in) :: t
  type(t_time_step), intent(in) :: tstep
  real(p_double), intent(in) :: tmin, tmax
  type( t_options ), intent(in) :: sim_options

  integer :: i, size_seed
  integer, dimension(:), pointer :: seed

  ! call superclass init
  ! this will initialize (and read checkpoint data if required) charged particle species
  call this % t_particles_cyl_modes % init( g_space, jay, emf, grid, no_co, bnd, &
                                            zpulse_list, ndump_fac, restart, &
                                            restart_handle, t, tstep, tmin, tmax, &
                                            sim_options )


  ! diagnostics have an additional setup to handle QED specific reports
  ! actually not yet...
  ! call setup_diag_qed( this )

  ! call individual photons init for each qed group
  if ( this%num_qed > 0 ) then
     do i=1, this%num_qed
        call this%qed_group(i)%photons%init( i, &
                    this%interpolation, this%grid_center, grid, g_space, emf, jay, &
                    no_co, bnd%send_vdf, bnd%recv_vdf, bnd%bnd_cross, bnd%node_cross, bnd%send_spec, &
                    bnd%recv_spec, ndump_fac, restart, restart_handle, sim_options, tstep, tmin, tmax )
     enddo

     ! initialize random number generator
     ! (for performance reasons the QED code uses the compiler issued random number gen.)

     call random_seed(SIZE=size_seed)
     call alloc(seed, (/ size_seed /) )
     do i=1, size_seed
        seed(i) = ( no_co%my_aid() + 1)
     enddo
     call random_seed(PUT=seed)

     ! no need to keep this around
     call freemem( seed )
  endif


end subroutine init_part_qedcyl
!-----------------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Printout the algorithm used by each species pusher
!-------------------------------------------------------------------------------
subroutine list_algorithm_part_qedcyl( this )

  implicit none
  class( t_particles_qedcyl ), intent(in) :: this
  integer :: i
  class( t_species ), pointer :: species

  print *, ' '
  print '(A)', 'Particles:'

  print '(A,I0)', '- Number of QED groups: ', this%num_qed

  if ( this%num_qed > 0 ) then

  ! particle shape
  select case ( this%interpolation )
    case (p_linear)
    print '(A)', '- Linear (1st order) interpolation'
    case (p_quadratic)
    print '(A)', '- Quadratic (2nd order) interpolation'
    case (p_cubic)
    print '(A)', '- Cubic (3rd order) interpolation'
    case (p_quartic)
    print '(A)', '- Quartic (4th order) interpolation'
  end select

  endif

  if ( this%grid_center ) then
    print '(A)', '- Using spatially centered grid values for interpolation'
  else
    print '(A)', '- Using staggered Yee grid values for interpolation'
  endif

  ! Current deposition
  if ( this%low_jay_roundoff ) then
    print *, '- Depositing current using low roundoff algorithm'
  endif

  do i = 1, this%num_qed
    print *, ' '
    print '(A,I0)', '- QED group #', i
    call this%qed_group(i)%electrons%list_algorithm()
    call this%qed_group(i)%positrons%list_algorithm()
    call this%qed_group(i)%photons%list_algorithm()
  enddo

  ! List algorithm for standard (non QED) species
  if ( this%num_qed > 0 ) then
    species => this % qed_group(this%num_qed) % positrons % next
  else
    species => this % species
  endif

  if ( associated(species) ) then
    print *, ' '
    print '(A)', '- Standard species'
  endif

  do
    if ( .not. associated(species)) exit
    call species % list_algorithm()
    species => species % next
  enddo


end subroutine list_algorithm_part_qedcyl
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
subroutine move_window_part_qedcyl( this, g_space, grid, jay, no_co, bnd )

  implicit none

  class( t_particles_qedcyl ), intent(inout) :: this

  type( t_space ),     intent(in) :: g_space
  class( t_grid ),  intent(in) :: grid
  class( t_current ), intent(inout)   :: jay
  class( t_node_conf ), intent(in)     :: no_co
  class(t_bnd), intent(inout) :: bnd

  integer :: i

  ! Call superclass method to move charged particle species
  ! No need to call with class t_particles_cyl_modes, as this routine is not overriden
  call this % t_particles_cyl_modes % move_window( g_space, grid, jay, no_co, bnd )

  ! Move photons in QED groups, use of the public routine, not defined in cyl_modes but in the core
  do i=1, this%num_qed
    call move_window( this%qed_group(i)%photons, g_space, grid, jay, no_co, &
        bnd%bnd_cross, bnd%node_cross, bnd%send_spec, bnd%recv_spec )
  enddo

end subroutine move_window_part_qedcyl
!-------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Report charged particle and photon data
!-----------------------------------------------------------------------------------------
subroutine report_part_qedcyl( this, emf, g_space, grid, no_co, tstep, t, tmin, send_msg, recv_msg )

  implicit none

  class( t_particles_qedcyl ), intent(inout) :: this
  class( t_emf ),       intent(inout) :: emf
  type( t_space ),        intent(in) :: g_space
  class( t_grid ),         intent(in) :: grid
  class( t_node_conf ),    intent(in) :: no_co
  type( t_time_step ),    intent(in) :: tstep
  real(p_double),         intent(in) :: t, tmin
  type(t_vdf_msg), dimension(2), intent(inout) :: send_msg, recv_msg

  integer :: i

  ! execute superclass diagnostics
  ! this routine calls species % reports(...) which is overriden in t_species_qedcyl
  call this % t_particles_cyl_modes % report( emf, g_space, grid, no_co, tstep, t, tmin, send_msg, recv_msg )

  ! QED group diagnostics
  do i=1, this%num_qed
    call this % qed_group(i) % report( emf, g_space, grid, no_co, tstep, t, tmin, send_msg, recv_msg )
  enddo

end subroutine report_part_qedcyl
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! reports QED species group specific diagnostics
!-----------------------------------------------------------------------------------------
subroutine report_spec_group_qedcyl( this, emf, g_space, grid, no_co, tstep, t, tmin, send_msg, recv_msg )
!-----------------------------------------------------------------------------------------
    
    implicit none

    ! dummy variables
    class( t_qedcyl_spec_group ), intent(inout) :: this
    class( t_emf ), intent(inout) :: emf
    type( t_space ), intent(in) :: g_space
    class( t_grid ), intent(in) :: grid
    class( t_node_conf ), intent(in) :: no_co
    type( t_time_step ), intent(in) :: tstep
    real(p_double), intent(in) :: t, tmin
    type(t_vdf_msg), dimension(2), intent(inout) :: send_msg, recv_msg
    
    ! photons specific diagnostics
    call this % photons % report( emf, g_space, grid, no_co, tstep, t, tmin, send_msg, recv_msg )

end subroutine report_spec_group_qedcyl
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Update particle data at the boundaries
!-----------------------------------------------------------------------------------------
subroutine update_boundary_part_qedcyl( this, jay, g_space, no_co, dt, bnd )

  implicit none

  class( t_particles_qedcyl ), intent(inout) :: this
  class( t_current ),       intent(inout) :: jay

  type( t_space ),     intent(in) :: g_space
  class( t_node_conf ), intent(in) :: no_co
  real(p_double),     intent(in) :: dt
  class( t_bnd ), intent(inout) :: bnd

  integer :: i

  ! update boundary on superclass data (charged particles)
  call this % t_particles_cyl_modes % update_boundary( jay, g_space, no_co, dt, bnd )

  ! update boundary on photons, which are not part of this % t_particles % species
  do i = 1, this%num_qed
    call this % qed_group(i) % photons % update_boundary(  jay, no_co, dt, &
        bnd%bnd_cross, bnd%node_cross, bnd%send_spec, bnd%recv_spec )
  enddo

end subroutine update_boundary_part_qedcyl
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Report charged particle and photon energy data
!-----------------------------------------------------------------------------------------
subroutine report_energy_part_qedcyl( this, no_co, tstep, t, tmin )

  implicit none

  class( t_particles_qedcyl ), intent(in) :: this
  class( t_node_conf ), intent(in) :: no_co
  type( t_time_step ), intent(in) :: tstep
  real(p_double), intent(in) :: t, tmin

  integer :: i

  ! execute superclass diagnostics
  call this % t_particles_cyl_modes % report_energy( no_co, tstep, t, tmin )

  
  do i=1, this%num_qed
      ! report radiated energy
      call this % qed_group(i) % electrons % report_energy_rad( no_co, tstep, t, tmin )
      call this % qed_group(i) % positrons % report_energy_rad( no_co, tstep, t, tmin )

      ! report pairs created (coulomb Trident)
      call this % qed_group(i) % electrons % report_pairs_ct( no_co, tstep, t, tmin )
      call this % qed_group(i) % positrons % report_pairs_ct( no_co, tstep, t, tmin )

      ! report photon energy
      call this % qed_group(i) % photons % report_energy( no_co, tstep, t, tmin )

      ! report emitted pairs energy
      call this % qed_group(i) % photons % report_pairs( no_co, tstep, t, tmin )
  enddo


end subroutine report_energy_part_qedcyl
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! - Merge particles in QED groups
! - Create pairs and advance photons
! - Advance charged particles and deposit current
!-----------------------------------------------------------------------------------------
subroutine advance_deposit_qedcyl( this, emf, jay, tstep, t, no_co, options )

  implicit none

  class( t_particles_qedcyl ), intent(inout) :: this
  class( t_emf ), intent( inout )  ::  emf
  class( t_current ), intent(inout) :: jay

  real(p_double), intent(in) :: t
  type( t_time_step ), intent(in) :: tstep
  class( t_node_conf ), intent(in) :: no_co
  type( t_options ), intent(in) :: options

  integer :: i

  ! Collides qed particles
  ! if  (this%qed_collide_info%if_qed_collide)   then
  !   call this%collide( n( tstep ), jay % pf(1) % nx_, t, no_co, tstep  )
  ! endif

  ! ! Merge particles in QED groups
  ! do i = 1, this%num_qed
  !   call this%qed_group(i)%electrons%mergedel( n( tstep ), jay % pf(1) % nx_, t, no_co )
  !   call this%qed_group(i)%positrons%mergedel( n( tstep ), jay % pf(1) % nx_, t, no_co )
  !   call this%qed_group(i)%photons  %mergedel( n( tstep ), jay % pf(1) % nx_, t, no_co )
  ! enddo

  select type( emf )
  class is ( t_emf_cyl_modes )

  ! create pairs and advance photons
  do i = 1, this%num_qed
    call this % qed_group(i) % photons % advance( &
            this%qed_group(i)%electrons, &
            this%qed_group(i)%positrons, &
            emf, dt( tstep ), 0, 1 )
  enddo

  end select

  ! call superclass pusher to advance charged particles
  call this % t_particles_cyl_modes % advance_deposit( emf, jay, tstep, t, no_co, options )

end subroutine advance_deposit_qedcyl
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Cleanup particles_qed object
!-----------------------------------------------------------------------------------------
subroutine cleanup_part_qedcyl( this )
!-----------------------------------------------------------------------------------------

  implicit none

  class(t_particles_qedcyl) , intent(inout) :: this

  integer :: i

  ! cleanup superclass
  call this % t_particles_cyl_modes % cleanup()

  ! cleanup local data structures
  if ( this%num_qed > 0 ) then
    ! cleanup photons
    do i=1, this%num_qed
      call this%qed_group(i)%photons%cleanup()
    enddo

    deallocate( this%qed_group )
  endif

end subroutine cleanup_part_qedcyl
!-----------------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
subroutine write_checkpoint_part_qedcyl( this, restart_handle )

  implicit none

  class( t_particles_qedcyl ), intent(in) :: this
  type( t_restart_handle ), intent(inout) :: restart_handle

  integer :: i

  ! Call superclass method to save charged particle checkpoint data
  call this % t_particles_cyl_modes % write_checkpoint( restart_handle )

  ! Save photon data
  do i=1, this%num_qed
    call this%qed_group(i)%photons%write_checkpoint( restart_handle )
  enddo

end subroutine write_checkpoint_part_qedcyl
!-------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Get diagnostic buffer size requirements
!-----------------------------------------------------------------------------------------
subroutine get_diag_buffer_size_qedcyl( this, gnx, diag_buffer_size )

  !use m_species_diagnostics

  implicit none

  class( t_particles_qedcyl ),  intent(in) :: this
  integer, dimension(:), intent(in) :: gnx
  integer, intent(inout) :: diag_buffer_size

  integer :: i

  ! Get dignostic size from superclass
  call this % t_particles_cyl_modes % get_diag_buffer_size( gnx, diag_buffer_size )

  ! Also include photon species diagnostics
  do i=1, this%num_qed
    call this % qed_group(i) % photons % get_diag_buffer_size( gnx, diag_buffer_size )
  enddo


end subroutine get_diag_buffer_size_qedcyl
!-----------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! Sets the particle load.
! For iteration n = 0 the load is set for the whole simulation volume (since this is called
! when no node partitions exist yet), for n > 0 the load is set for the
! local volume/particles only
!---------------------------------------------------------------------------------------------------
subroutine add_load_part_qedcyl( this, grid, emf, n )

  implicit none

  class( t_particles_qedcyl ), intent(inout) :: this
  class( t_grid ), intent(inout) :: grid
  class( t_emf ), intent(in) :: emf
  integer, intent(in) :: n

  class( t_species ), pointer :: species
  integer :: i

  ! loop through all species and add to int_load array
  species => this % species
  if ( n > 0 ) then

    do
      if (.not. associated(species)) exit
      call species % add_particle_load( grid, emf )
      species => species % next
    enddo

    ! Add load for photons
    do i = 1, this % num_qed
      call this % qed_group(i) % photons % add_particle_load( grid, emf )
    enddo

  else

    do
      if (.not. associated(species)) exit
      call add_density_load( species, grid  )
      species => species % next
    enddo

  endif

end subroutine add_load_part_qedcyl
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine reshape_part_qedcyl( this, old_lb, new_lb, no_co, send_msg, recv_msg )
  !-----------------------------------------------------------------------------------------
  ! redistribute particles to new simulation partition
  !-----------------------------------------------------------------------------------------

  implicit none

  class( t_particles_qedcyl ), intent(inout) :: this
  class( t_grid ), intent(in) :: new_lb, old_lb
  class( t_node_conf ), intent(in) :: no_co
  type(t_vdf_msg), dimension(2), intent(inout) :: send_msg, recv_msg

  type( t_msg_patt ) :: msg_patt
  class( t_species ), pointer :: species
  integer :: n, i
  integer, dimension(2,p_max_dim) :: gc_num
  type(t_vdf_report), pointer :: report

  ! reshape species
  if ( this%num_species > 0 ) then
    ! get message pattern (guard cells are not required since no particles
    ! should be in guard cells). This is the same for all species
    ! so we only need to do it once.
    gc_num = 0
    call new_msg_patt( msg_patt, old_lb, new_lb, no_co, gc_num )

    ! loop through all species
    species => this % species
    do
      if (.not. associated(species)) exit
      call species % reshape( old_lb, new_lb, msg_patt, no_co, send_msg, recv_msg )
      species => species % next
    enddo

    ! reshape the photon species
    do i = 1, this % num_qed
      call this % qed_group(i) % photons % reshape( old_lb, new_lb, msg_patt, no_co, &
                                                    send_msg, recv_msg )
    enddo

    ! clear message pattern data
    call msg_patt % cleanup()
  endif

  ! reshape local vdf objects if needed
  if ( this%jay_tmp%x_dim_ > 0 ) call reshape_nocopy( this%jay_tmp, new_lb )

  call reshape( this%charge, old_lb, new_lb, no_co, send_msg, recv_msg )

  ! Reshape tavg_data
  report => this%reports
  do
    if ( .not. associated(report) ) exit
    if ( report%tavg_data%x_dim_ > 0 ) then
      call reshape_copy( report%tavg_data, old_lb, new_lb, no_co, send_msg, recv_msg )
    endif
    report => report%next
  enddo

end subroutine reshape_part_qedcyl
!-----------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! Reports number of particles per node for all nodes (small grid file)
!---------------------------------------------------------------------------------------------------
subroutine report_node_load_part_qedcyl( this, n, ndump, t, no_co )

  implicit none

  class( t_particles ), intent(in) :: this
  integer, intent(in)             :: n, ndump
  real( p_double ), intent(in)    :: t
  class( t_node_conf ), intent(in) :: no_co

  real( p_double ), dimension(:), pointer     :: node_load
  real( p_double ), dimension(:), pointer     :: f1
  real( p_double ), dimension(:,:), pointer   :: f2
  real( p_double ), dimension(:,:,:), pointer :: f3

  class( t_diag_file ), allocatable :: diagFile

  integer, dimension(p_max_dim) :: lnx

  class( t_species ), pointer :: species
  integer :: i, npart, nnodes, ierr
  real( p_double ) :: npart_dbl

  node_load => null()
  f1 => null()
  f2 => null()
  f3 => null()

  ! get total number of particles on node
  npart = 0
  species => this % species
  do
      if (.not. associated(species)) exit
      npart = npart + species%num_par
      species => species % next
  enddo

  ! Add photons
  select type( this )
  class is( t_particles_qedcyl )
    do i = 1, this % num_qed
      npart = npart + this % qed_group(i) % photons % num_par
    end do
  end select  

  npart_dbl = npart

  ! gather data
  nnodes =  no_num( no_co )
  call alloc( node_load, (/nnodes/) )

  if ( nnodes > 1 ) then
      call mpi_gather( npart_dbl, 1, MPI_DOUBLE_PRECISION, &
      node_load, 1, MPI_DOUBLE_PRECISION, &
      0, comm(no_co), ierr )
      if ( ierr /= MPI_SUCCESS ) then
          ERROR('MPI Error')
          call abort_program( p_err_mpi )
      endif
  else
      node_load = npart_dbl
  endif

  ! save data to disk
  if ( root( no_co ) ) then
      
      lnx(1:p_x_dim) = nx( no_co )
      
      ! Open the output file
      call create_diag_file( diagFile )
      diagFile%ftype   = p_diag_grid
      
      diagFile%filepath  = trim(path_mass) // 'LOAD' // p_dir_sep // 'NODE' // p_dir_sep
      diagFile%filename  = trim(get_filename(n/ndump, 'node_load'))
      diagFile%name      = 'Particles per node'
      diagFile%iter%n         = n
      diagFile%iter%t         = t
      diagFile%iter%time_units = '1 / \omega_p'
      
      
      diagFile%grid % ndims = p_x_dim
      diagFile%grid % name  = "load"
      diagFile%grid % label = "particles per node"
      diagFile%grid % units = "particles"
      
      do i = 1, p_x_dim
          diagFile % grid % axis(i) % min = 0
          diagFile % grid% axis(i) % max = lnx(i)
          write( diagFile % grid% axis(i) % name, '(A,I0)' ) 'x', i
          write( diagFile % grid% axis(i) % label, '(A,I0)' ) 'x_', i
          diagFile % grid% axis(i) % units = 'cell'
      enddo
      
      call diagFile % open( p_diag_create )
      
      ! Move node load values to the proper position and write node load values
      select case (p_x_dim)
      case(1)
          call alloc( f1, lnx )
          f1 = -1.0
          do i = 1, nnodes
              f1( no_co%ngp( i, 1 ) ) = node_load(i)
          enddo
          call diagFile % add_dataset( "load", f1 )
          call freemem( f1 )
          
      case(2)
          call alloc( f2, lnx )
          f2 = -1.0
          do i = 1, nnodes
              f2( no_co%ngp( i, 1 ), no_co%ngp( i, 2 ) ) = node_load(i)
          enddo
          call diagFile % add_dataset( "load", f2 )
          call freemem( f2 )
          
      case(3)
          call alloc( f3, lnx )
          f3 = -1.0
          do i = 1, nnodes
              f3( no_co%ngp( i, 1 ), no_co%ngp( i, 2 ), no_co%ngp( i, 3 ) ) = node_load(i)
          enddo
          call diagFile % add_dataset( "load", f3 )
          call freemem( f3 )
      end select
      
      ! close output file
      call diagFile % close( )
      deallocate( diagFile )
  endif

  ! free load array
  call freemem( node_load )

end subroutine report_node_load_part_qedcyl
!---------------------------------------------------------------------------------------------------

subroutine init_bremsstrahlung( this, input_file )

  implicit none

  class( t_particles_qedcyl ), intent(inout) :: this
  class( t_input_file ), intent(inout) :: input_file

  class( t_species ), pointer :: sp1, sp2
  integer :: i, i1, j1, k, num_qed, num_species
  real(p_double), dimension(:,:), allocatable :: couples_br
  logical :: condition, success_linked

  num_qed = this % num_qed
  num_species = this % num_species
  allocate(couples_br(2, num_qed))

  ! Then we activate read the input and eventually activate bremsstrahlung
  if (num_species .gt. 2 * num_qed) then

    ! Read the input file bloc named qed_bremsstrahlung
    SCR_ROOT("   - reading QED bremsstrahlung configuration...")
    call this % qed_group(1) % electrons % bremsstrahlung_info % read_input(input_file, this%num_qed)

    ! Bremsstrahlung photon production is set to false for all electron species
    do i=1, num_qed
        this % qed_group(i) % electrons % bremsstrahlung_info % if_bremsstrahlung = .false.
        this % qed_group(i) % positrons % bremsstrahlung_info % if_bremsstrahlung = .false.
    enddo

    do i=1, num_qed
        couples_br(1,i) = i
        couples_br(2,i) = this % qed_group(1) % electrons % bremsstrahlung_info % i_ion(i)
    enddo

    ! For now, we only allow Bremsstrahlung with the same ion species
    condition = ( (.not. all(couples_br(2,1:num_qed) .eq. couples_br(2,1)) ) .and. all(couples_br(2,1:num_qed) .ge. 1) .and. all(couples_br(2,1:num_qed) .le. num_species-2*num_qed) )
    if (condition) then
        print *, ""
        print *, "   Error reading qed_bremsstrahlung parameters "
        print *, "   All electrons must be linked with same ion species "
        print *, "   aborting..."
        stop
    endif

    ! Copy the info for other qed groups
    if (num_qed .gt. 1) then
      do i=2, num_qed
          this % qed_group(i) % electrons % bremsstrahlung_info % Z_ion = this % qed_group(1) % electrons % bremsstrahlung_info % Z_ion
          this % qed_group(i) % electrons % bremsstrahlung_info % proba_mult = this % qed_group(1) % electrons % bremsstrahlung_info % proba_mult
          this % qed_group(i) % electrons % bremsstrahlung_info % cs_model = this % qed_group(1) % electrons % bremsstrahlung_info % cs_model
          this % qed_group(i) % electrons % bremsstrahlung_info % energy_damp = this % qed_group(1) % electrons % bremsstrahlung_info % energy_damp
          allocate(this % qed_group(i) % electrons % bremsstrahlung_info % i_ion(num_qed))
          this % qed_group(i) % electrons % bremsstrahlung_info % i_ion = this % qed_group(1) % electrons % bremsstrahlung_info % i_ion
      enddo
    endif

    ! We activate Bremsstrahlung by linking the electron and ion species specified in the input
    do k=1, num_qed
      if ( (couples_br(2,k) .ge. 1) .and. (couples_br(2,k) .le. num_species-2*num_qed) ) then
        j1 = int(couples_br(1,k))
        sp1 => this % qed_group(j1) % electrons
        sp2 => this % species; i1 = 1
        do
          if (.not. associated(sp2)) exit
            if (i1 .eq. int(2*num_qed+couples_br(2,k))) then
              select type( sp2 )
                class is( t_species_cyl_modes )
                this % qed_group(j1) % electrons % bremsstrahlung_info % ion_cyl_m => sp2
                this % qed_group(j1) % electrons % bremsstrahlung_info % ion => sp2
                this % qed_group(j1) % electrons % bremsstrahlung_info % if_bremsstrahlung = .true.
                SCR_ROOT("Bremsstrahlung activated between QED group ", k, " and ",trim(sp2%name))
                exit
              end select
            endif
            sp2 => sp2 % next; i1 = i1 + 1
        enddo
      else
        this % qed_group(k) % electrons % bremsstrahlung_info % if_bremsstrahlung = .false.
      endif
    end do

    ! For Bremsstrahlung of positrons, we copy the parameters selected for electrons
    if (num_qed .ge. 1) then
      do i=1, num_qed
          this % qed_group(i) % positrons % bremsstrahlung_info % if_bremsstrahlung = this % qed_group(i) % electrons % bremsstrahlung_info % if_bremsstrahlung
          this % qed_group(i) % positrons % bremsstrahlung_info % Z_ion = this % qed_group(i) % electrons % bremsstrahlung_info % Z_ion
          this % qed_group(i) % positrons % bremsstrahlung_info % proba_mult = this % qed_group(i) % electrons % bremsstrahlung_info % proba_mult
          this % qed_group(i) % positrons % bremsstrahlung_info % cs_model = this % qed_group(i) % electrons % bremsstrahlung_info % cs_model
          this % qed_group(i) % positrons % bremsstrahlung_info % energy_damp = this % qed_group(i) % electrons % bremsstrahlung_info % energy_damp
          allocate(this % qed_group(i) % positrons % bremsstrahlung_info % i_ion(num_qed))
          this % qed_group(i) % positrons % bremsstrahlung_info % i_ion = this % qed_group(i) % electrons % bremsstrahlung_info % i_ion
          this % qed_group(i) % positrons % bremsstrahlung_info % ion => this % qed_group(i) % electrons % bremsstrahlung_info % ion
          this % qed_group(i) % positrons % bremsstrahlung_info % ion_cyl_m => this % qed_group(i) % electrons % bremsstrahlung_info % ion_cyl_m
      enddo
    endif
  else

    do i=1, num_qed
      this % qed_group(i) % electrons % bremsstrahlung_info % if_bremsstrahlung = .false.
      this % qed_group(i) % positrons % bremsstrahlung_info % if_bremsstrahlung = .false.
    enddo
  endif


end subroutine init_bremsstrahlung

subroutine init_betheheitler( this, input_file )

  implicit none

  class( t_particles_qedcyl ), intent(inout) :: this
  class( t_input_file ), intent(inout) :: input_file

  class( t_species ), pointer :: sp1, sp2
  integer :: i, i1, j1, k, num_qed, num_species
  real(p_double), dimension(:,:), allocatable :: couples_bh
  logical :: condition

  num_qed = this % num_qed
  num_species = this % num_species
  allocate(couples_bh(2, num_qed))

  ! Then we read bethe heitler input only if there is one normal species
  if (num_species .gt. 2 * num_qed) then

    ! Read the input file bloc named qed_betheheitler
    SCR_ROOT("   - reading QED bethe heitler configuration...")
    call this % qed_group(1) % photons % betheheitler_info % read_input(input_file, this%num_qed)

    ! Bethe Heitler pair production is set to disactivated by default
    do i=1, num_qed
        this % qed_group(i) % photons % betheheitler_info % if_betheheitler = .false.
    enddo

    do i=1, num_qed
        couples_bh(1,i) = i
        couples_bh(2,i) = this % qed_group(1) % photons % betheheitler_info % i_ion(i)
    enddo

    ! For now, we only allow Bethe Heitler with the same ion species
    condition = ( (.not. all(couples_bh(2,1:num_qed) .eq. couples_bh(2,1)) ) .and. all(couples_bh(2,1:num_qed) .ge. 1) .and. all(couples_bh(2,1:num_qed) .le. num_species-2*num_qed) )
    if (condition) then
        print *, ""
        print *, "   Error reading qed_betheheitler parameters "
        print *, "   All photons must be linked with same ion species "
        print *, "   aborting..."
        stop
    endif

    ! Copy the info for other qed groups
    if (num_qed .gt. 1) then
      do i=2, num_qed
          this % qed_group(i) % photons % betheheitler_info % Z_ion = this % qed_group(1) % photons % betheheitler_info % Z_ion
          this % qed_group(i) % photons % betheheitler_info % proba_mult = this % qed_group(1) % photons % betheheitler_info % proba_mult
          this % qed_group(i) % photons % betheheitler_info % cs_model = this % qed_group(1) % photons % betheheitler_info % cs_model
          this % qed_group(i) % photons % betheheitler_info % decay_photons = this % qed_group(1) % photons % betheheitler_info % decay_photons
          allocate(this % qed_group(i) % photons % betheheitler_info % i_ion(num_qed))
          this % qed_group(i) % photons % betheheitler_info % i_ion = this % qed_group(1) % photons % betheheitler_info % i_ion
      enddo
    endif

    ! Link the photon and ion species
    do k=1, num_qed
      if ( (couples_bh(2,k) .ge. 1) .and. (couples_bh(2,k) .le. num_species-2*num_qed) ) then
        j1 = int(couples_bh(1,k))
        sp1 => this % qed_group(j1) % photons
        sp2 => this % species; i1 = 1
        do
          if (.not. associated(sp2)) exit
            if (i1 .eq. int(2*num_qed+couples_bh(2,k))) then
              select type( sp2 )
                class is( t_species_cyl_modes )
                this % qed_group(j1) % photons % betheheitler_info % ion_cyl_m => sp2
                this % qed_group(j1) % photons % betheheitler_info % ion => sp2
                this % qed_group(j1) % photons % betheheitler_info % if_betheheitler = .true.
                SCR_ROOT("Bethe Heitler activated between QED group ", k, " and ",trim(sp2%name))
                exit
              end select
            endif
            sp2 => sp2 % next; i1 = i1 + 1
        enddo
      else
        this % qed_group(k) % photons % betheheitler_info % if_betheheitler = .false.
      endif
    end do
  else
    do i=1, num_qed
      this % qed_group(i) % photons % betheheitler_info % if_betheheitler = .false.
    enddo
  endif

end subroutine init_betheheitler

subroutine init_trident_coul( this, input_file )

  implicit none

  class( t_particles_qedcyl ), intent(inout) :: this
  class( t_input_file ), intent(inout) :: input_file

  class( t_species ), pointer :: sp1, sp2
  integer :: i, i1, j1, k, num_qed, num_species
  real(p_double), dimension(:,:), allocatable :: couples_ct
  logical :: condition, success_linked

  num_qed = this % num_qed
  num_species = this % num_species
  allocate(couples_ct(2, num_qed))

  ! Then we activate read the input and eventually activate trident coul
  if (num_species .gt. 2 * num_qed) then

    ! Read the input file bloc named qed_trident_coul
    SCR_ROOT("   - reading QED trident_coul configuration...")
    call this % qed_group(1) % electrons % trident_coul_info % read_input(input_file, this%num_qed)

    ! trident_coul pair production is set to false for all electron species
    do i=1, num_qed
        this % qed_group(i) % electrons % trident_coul_info % if_trident_coul = .false.
        this % qed_group(i) % positrons % trident_coul_info % if_trident_coul = .false.
    enddo

    do i=1, num_qed
        couples_ct(1,i) = i
        couples_ct(2,i) = this % qed_group(1) % electrons % trident_coul_info % i_ion(i)
    enddo

    ! For now, we only allow trident_coul with the same ion species
    condition = ( (.not. all(couples_ct(2,1:num_qed) .eq. couples_ct(2,1)) ) .and. all(couples_ct(2,1:num_qed) .ge. 1) .and. all(couples_ct(2,1:num_qed) .le. num_species-2*num_qed) )
    if (condition) then
        print *, ""
        print *, "   Error reading qed_trident_coul parameters "
        print *, "   All electrons must be linked with same ion species "
        print *, "   aborting..."
        stop
    endif

    ! Copy the info for other qed groups
    if (num_qed .gt. 1) then
      do i=2, num_qed
          this % qed_group(i) % electrons % trident_coul_info % if_trident_coul = this % qed_group(1) % electrons % trident_coul_info % if_trident_coul
          this % qed_group(i) % electrons % trident_coul_info % Z_ion = this % qed_group(1) % electrons % trident_coul_info % Z_ion
          this % qed_group(i) % electrons % trident_coul_info % proba_mult = this % qed_group(1) % electrons % trident_coul_info % proba_mult
          this % qed_group(i) % electrons % trident_coul_info % cs_model = this % qed_group(1) % electrons % trident_coul_info % cs_model
          this % qed_group(i) % electrons % trident_coul_info % energy_damp = this % qed_group(1) % electrons % trident_coul_info % energy_damp
          allocate(this % qed_group(i) % electrons % trident_coul_info % i_ion(num_qed))
          this % qed_group(i) % electrons % trident_coul_info % i_ion = this % qed_group(1) % electrons % trident_coul_info % i_ion
      enddo
    endif

    ! We activate trident_coul by linking the electron and ion species specified in the input
    do k=1, num_qed
      if ( (couples_ct(2,k) .ge. 1) .and. (couples_ct(2,k) .le. num_species-2*num_qed) ) then
        j1 = int(couples_ct(1,k))
        sp1 => this % qed_group(j1) % electrons
        sp2 => this % species; i1 = 1
        do
          if (.not. associated(sp2)) exit
            if (i1 .eq. int(2*num_qed+couples_ct(2,k))) then
              select type( sp2 )
                class is( t_species_cyl_modes )
                this % qed_group(j1) % electrons % trident_coul_info % ion_cyl_m => sp2
                this % qed_group(j1) % electrons % trident_coul_info % ion => sp2
                this % qed_group(j1) % electrons % trident_coul_info % if_trident_coul = .true.
                SCR_ROOT("Trident coul activated between QED group ", k, " and ",trim(sp2%name))
                exit
              end select
            endif
            sp2 => sp2 % next; i1 = i1 + 1
        enddo
      else
        this % qed_group(k) % electrons % trident_coul_info % if_trident_coul = .false.
      endif
    end do

    ! For trident_coul of positrons, we copy the parameters selected for electrons
    if (num_qed .ge. 1) then
      do i=1, num_qed
          this % qed_group(i) % positrons % trident_coul_info % if_trident_coul = this % qed_group(i) % electrons % trident_coul_info % if_trident_coul
          this % qed_group(i) % positrons % trident_coul_info % Z_ion = this % qed_group(i) % electrons % trident_coul_info % Z_ion
          this % qed_group(i) % positrons % trident_coul_info % proba_mult = this % qed_group(i) % electrons % trident_coul_info % proba_mult
          this % qed_group(i) % positrons % trident_coul_info % cs_model = this % qed_group(i) % electrons % trident_coul_info % cs_model
          this % qed_group(i) % positrons % trident_coul_info % energy_damp = this % qed_group(i) % electrons % trident_coul_info % energy_damp
          allocate(this % qed_group(i) % positrons % trident_coul_info % i_ion(num_qed))
          this % qed_group(i) % positrons % trident_coul_info % i_ion = this % qed_group(i) % electrons % trident_coul_info % i_ion
          this % qed_group(i) % positrons % trident_coul_info % ion => this % qed_group(i) % electrons % trident_coul_info % ion
          this % qed_group(i) % positrons % trident_coul_info % ion_cyl_m => this % qed_group(i) % electrons % trident_coul_info % ion_cyl_m
      enddo
    endif

  else
    do i=1, num_qed
      this % qed_group(i) % electrons % trident_coul_info % if_trident_coul = .false.
      this % qed_group(i) % positrons % trident_coul_info % if_trident_coul = .false.
    enddo

  endif
end subroutine init_trident_coul

end module m_particles_qedcyl