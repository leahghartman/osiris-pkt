!
! The current version does not allow for normal species, only pkt groups
!

! TODO: match the result of paper as closely as possible (06/12/2026)
! TODO: work on the quasi-3d version of the code (06/13/2026)
! TODO: talk to neil about getting the correctly matched parameters for AWAKE

#include "os-config.h"
#include "os-preprocess.fpp"

module m_particles_pkt

#include "memory/memory.h"

   use m_system
   use m_parameters
   use m_input_file, only: t_input_file, get_namelist, disp_out
   use m_vdf_define, only: t_vdf, p_max_reports_len, p_max_reports, p_n_report_type, &
                           p_full, p_savg, p_senv, p_line, p_slice, t_vdf_report
   use m_vdf_report, only: if_report
   use m_vdf_comm
   use m_particles_define, only: t_particles, pushev, partboundev, diag_part_ev
   use m_particles, only: cleanup

   use m_space, only: t_space
   use m_emf_define, only: t_emf
   use m_restart, only: t_restart_handle, restart_io_write, restart_io_read
   use m_species_define_pkt, only: t_species_pkt
   use m_photons_define_pkt, only: t_photons_pkt
   use m_species_define, only: t_species, p_max_spname_len, p_cell_low, p_cell_near
   use m_particles_define, only: p_report_quants
   use m_particles_charge, only: reshape
   use m_grid_define, only: t_grid, t_msg_patt
   use m_grid_parallel, only: new_msg_patt
   use m_species_collisions, only: read_nml
   use m_vdf_report, only: new
   use m_current_define, only: t_current
   use m_bnd, only: t_bnd
   use m_time_step, only: t_time_step, dt, n, test_if_report, ndump
   use m_species, only: move_window
   use m_diagfile
   use m_node_conf
   use m_diagnostic_utilities
   use m_zpulse, only: t_zpulse_list
   use m_species_loadbalance, only: add_density_load
   use m_logprof, only: begin_event, end_event
   use m_species_current
   use m_dep_density_pkt

   implicit none

   private

! Group of 2 species + 1 photon species for pkt code

   type :: t_pkt_spec_group

      class(t_photons_pkt), pointer :: photons
      class(t_species_pkt), pointer :: electrons ! electrons and positrons point to species

   contains

      procedure :: read_input => read_input_spec_group_pkt
      procedure :: report => report_spec_group_pkt

   end type t_pkt_spec_group

   type, extends(t_particles) :: t_particles_pkt

      ! pkt Species groups
      integer :: num_pkt

      type(t_pkt_spec_group), dimension(:), pointer :: pkt_group => NULL()

      ! load balance parameters
      real(p_double) :: dlb_fld2_thresh
      real(p_double) :: dlb_part_weight

   contains

      procedure :: allocate_objs => allocate_objs_part_pkt
      procedure :: read_input => read_input_part_pkt
      procedure :: list_algorithm => list_algorithm_part_pkt
      procedure :: init => init_part_pkt
      procedure :: cleanup => cleanup_part_pkt

      ! also update boundary on photons
      procedure :: update_boundary => update_boundary_part_pkt

      ! also report on photons
      procedure :: report => report_part_pkt
      procedure :: report_energy => report_energy_part_pkt

      ! include photon diagnostics
      procedure :: get_diag_buffer_size => get_diag_buffer_size_pkt

      procedure :: advance_deposit => advance_deposit_pkt
      procedure :: write_checkpoint => write_checkpoint_part_pkt
      procedure :: move_window => move_window_part_pkt

      procedure :: add_load => add_load_part_pkt
      procedure :: reshape => reshape_part_pkt

   end type t_particles_pkt

   interface report_node_load_part_pkt
      module procedure report_node_load_part_pkt
   end interface

   public :: t_pkt_spec_group, t_particles_pkt, report_node_load_part_pkt

contains

!-----------------------------------------------------------------------------------------
! read necessary information from input file
!-----------------------------------------------------------------------------------------
   subroutine read_input_spec_group_pkt(this, input_file, periodic, if_move, grid, dt, sim_options)
!-----------------------------------------------------------------------------------------

      implicit none

      class(t_pkt_spec_group), intent(inout) :: this
      class(t_input_file), intent(inout) :: input_file
      logical, dimension(:), intent(in) :: periodic, if_move
      class(t_grid), intent(in) :: grid
      real(p_double), intent(in) :: dt
      type(t_options), intent(in) :: sim_options

      integer :: ierr
      character(len=p_max_spname_len) :: spname
      character(len=20) :: pairprod
      logical :: if_damp_classical, if_damp_pkt, if_pairprod
      real(p_double) :: pkt_g_cutoff, p_emit_cutoff, zeta_p, zeta_g, &
                        pairprod_gthr, pairprod_gpair, pp_chi_cutoff

      real(p_k_part) :: omega0

      namelist /nl_pkt_group/ if_damp_classical, if_damp_pkt, if_pairprod, &
         pkt_g_cutoff, p_emit_cutoff, zeta_p, zeta_g, &
         pairprod, pairprod_gthr, pairprod_gpair, &
         pp_chi_cutoff

      if_damp_classical = .false.
      if_damp_pkt = .false.
      if_pairprod = .false.
      pkt_g_cutoff = 10.0_p_double
      p_emit_cutoff = 2.0_p_double
      pp_chi_cutoff = 0.1_p_double
      omega0 = 1.0_p_k_part

      zeta_p = 1.0_p_double
      zeta_g = 1.0_p_double

      pairprod = "pkt"

      pairprod_gthr = 1000._p_double
      pairprod_gpair = 100._p_double

      SCR_ROOT("Reading pkt species group configuration...")

      call get_namelist(input_file, "nl_pkt_group", ierr)

      if (ierr == 0) then
         read (input_file%nml_text, nml=nl_pkt_group, iostat=ierr)
         if (ierr /= 0) then
            if (mpi_node() == 0) then
               write (0, *) ""
               write (0, *) "   Error reading pkt species group information"
               write (0, *) "   aborting..."
            end if
            stop
         end if
      else
         SCR_ROOT("   - pkt species group parameters missing, using default settings")
      end if

      call read_nml_pkt_group_diag(this, input_file)

      ! read electrons species
      SCR_ROOT("Reading associated electrons configuration...")

      call this%electrons%read_input(input_file, spname, periodic, &
                                     if_move, grid, dt, .true., sim_options)

      ! force electrons and set emission variables
      this%electrons%rqm = -1.0
      this%electrons%if_damp_classical = if_damp_classical
      this%electrons%if_damp_pkt = if_damp_pkt
      this%electrons%pkt_g_cutoff = pkt_g_cutoff
      this%electrons%p_emit_cutoff = p_emit_cutoff

      ! read photons species
      SCR_ROOT("Reading associated photons configuration...")

      call this%photons%read_input(input_file, spname, periodic, if_move, grid, &
                                   dt, .true., sim_options)

      ! force photons and set pair production flag
      this%photons%rqm = 0.
      this%photons%if_pairprod = if_pairprod
      this%photons%pp_chi_cutoff = pp_chi_cutoff

      ! pkt rescalling
      this%electrons%pkt_zeta = zeta_p
      this%photons%pkt_zeta = zeta_g

      ! process pair production model
      select case (trim(pairprod))
      case ("pkt")
         this%electrons%pairprod = p_pairprod_qed
         this%photons%pairprod = p_pairprod_qed
      case ("gthr")
         this%electrons%pairprod = p_pairprod_gthr
         this%photons%pairprod = p_pairprod_gthr
      case default
         this%electrons%pairprod = p_pairprod_qed
         this%photons%pairprod = p_pairprod_qed
      end select

      this%electrons%pairprod_gthr = pairprod_gthr

      this%electrons%pairprod_gpair = pairprod_gpair
      this%photons%pairprod_gpair = pairprod_gpair

   end subroutine read_input_spec_group_pkt
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! read necessary information from input file
!-----------------------------------------------------------------------------------------
   subroutine read_nml_pkt_group_diag(this, input_file)
!-----------------------------------------------------------------------------------------

      implicit none

      class(t_pkt_spec_group), intent(inout) :: this
      class(t_input_file), intent(inout) :: input_file

      integer :: ierr, ndump_fac_rad, ndump_fac_pairs, ndump_fac_radspect, radspect_bins, &
                 ndump_fac_radanglespect, radanglespect_bins, &
                 chi_emit_nbins, ndump_fac_radsphere, radsphere_nbins, &
                 ndump_fac_raddetector, raddetector_nbins, raddetector_dir

      logical :: if_add_classical_radsphere

      real(p_double) :: radspect_emin, radspect_emax, radanglespect_thresh, chi_emit_min, &
                        chi_emit_max, raddetector_emin, raddetector_emax, raddetector_phimax, &
                        raddetector_phimin, raddetector_thetamax, raddetector_thetamin, &
                        radsphere_emin, radsphere_emax

      character(len=11) :: detector_direction

      namelist /nl_pkt_group_diag/ ndump_fac_rad, ndump_fac_pairs, &
         ndump_fac_radspect, radspect_emin, radspect_emax, radspect_bins, &
         ndump_fac_radanglespect, radanglespect_bins, radanglespect_thresh, &
         chi_emit_nbins, chi_emit_min, chi_emit_max, &
         ndump_fac_radsphere, radsphere_nbins, radsphere_emin, radsphere_emax, if_add_classical_radsphere, &
         ndump_fac_raddetector, raddetector_nbins, raddetector_emin, raddetector_emax, &
         raddetector_phimax, raddetector_phimin, raddetector_thetamax, raddetector_thetamin, detector_direction

      ndump_fac_rad = 0
      ndump_fac_pairs = 0

      ndump_fac_radspect = 0
      radspect_emin = 1.957e-5
      radspect_emax = 1.957e5
      radspect_bins = 400

      ndump_fac_radanglespect = 0
      radanglespect_bins = 100
      radanglespect_thresh = 1.0

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
      raddetector_phimin = -0.05
      raddetector_phimax = 0.05
      raddetector_thetamin = -0.05
      raddetector_thetamax = 0.05

      call get_namelist(input_file, "nl_pkt_group_diag", ierr)

      if (ierr == 0) then
         read (input_file%nml_text, nml=nl_pkt_group_diag, iostat=ierr)
         if (ierr /= 0) then
            if (mpi_node() == 0) then
               write (0, *) ""
               write (0, *) "   Error reading pkt species group diagnostics information"
               write (0, *) "   aborting..."
            end if
            stop
         end if
      else
         SCR_ROOT("   - pkt species group diagnostics parameters missing, using default settings")
      end if

      select case (trim(detector_direction))
      case ("x1_positive")
         raddetector_dir = 1
      case ("x1_negative")
         raddetector_dir = 2
      case ("x2_positive")
         raddetector_dir = 3
      case ("x2_negative")
         raddetector_dir = 4
      case ("x3_positive")
         raddetector_dir = 5
      case ("x3_negative")
         raddetector_dir = 6
      case default
         print *, "Error reading pkt diag parameters"
         print *, "Invalid value for the detector_dir parameter"
         print *, "Available detector directions are 'x1_positive', 'x1_negative',"
         print *, "'x2_positive', 'x2_negative', 'x3_positive', 'x3_negative' "
         print *, "aborting..."
         stop
      end select

      ! radiated energy diagnostic
      this%electrons%ndump_fac_rad = ndump_fac_rad

      this%electrons%radspect_emin = radspect_emin
      this%electrons%radspect_emax = radspect_emax
      this%electrons%radspect_bins = radspect_bins

      this%electrons%chi_emit_min = chi_emit_min
      this%electrons%chi_emit_max = chi_emit_max
      this%electrons%chi_emit_nbins = chi_emit_nbins

      this%electrons%radanglespect_bins = radanglespect_bins
      this%electrons%radanglespect_thresh = radanglespect_thresh

      this%electrons%radsphere_nbins = radsphere_nbins
      this%electrons%radsphere_emin = radsphere_emin
      this%electrons%radsphere_emax = radsphere_emax
      this%electrons%if_add_classical_radsphere = if_add_classical_radsphere

      ! radiation detector
      this%electrons%raddetector_nbins = raddetector_nbins

      this%electrons%raddetector_dir = raddetector_dir

      this%electrons%raddetector_emin = raddetector_emin
      this%electrons%raddetector_emax = raddetector_emax

      this%electrons%raddetector_phimin = raddetector_phimin
      this%electrons%raddetector_phimax = raddetector_phimax

      this%electrons%raddetector_thetamin = raddetector_thetamin
      this%electrons%raddetector_thetamax = raddetector_thetamax

   end subroutine read_nml_pkt_group_diag
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
! reports pkt species group specific diagnostics
!-----------------------------------------------------------------------------------------
   subroutine report_spec_group_pkt(this, emf, g_space, grid, no_co, tstep, t, tmin, send_msg, recv_msg)
!-----------------------------------------------------------------------------------------

      implicit none

      ! dummy variables
      class(t_pkt_spec_group), intent(inout) :: this
      class(t_emf), intent(inout) :: emf
      type(t_space), intent(in) :: g_space
      class(t_grid), intent(in) :: grid
      class(t_node_conf), intent(in) :: no_co
      type(t_time_step), intent(in) :: tstep
      real(p_double), intent(in) :: t, tmin
      type(t_vdf_msg), dimension(2), intent(inout) :: send_msg, recv_msg

      ! photons specific diagnostics
      call this%photons%report(emf, g_space, grid, no_co, tstep, t, tmin, send_msg, recv_msg)

   end subroutine report_spec_group_pkt
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! read necessary information from input file
!-----------------------------------------------------------------------------------------
   subroutine read_input_part_pkt(this, input_file, periodic, if_move, grid, dt, sim_options)
!-----------------------------------------------------------------------------------------

      implicit none

      ! dummy variables

      class(t_particles_pkt), intent(inout) :: this
      class(t_input_file), intent(inout) :: input_file
      logical, dimension(:), intent(in) :: periodic, if_move
      class(t_grid), intent(in) :: grid
      real(p_double), intent(in) :: dt
      type(t_options), intent(in) :: sim_options

      integer  :: i, j, num_pkt, num_species
      logical  :: low_jay_roundoff

      character(len=p_max_reports_len), dimension(p_max_reports) :: reports
      integer  :: ndump_fac, ndump_fac_ave, ndump_fac_lineout, ndump_fac_ene
      integer, dimension(p_x_dim) :: n_ave
      integer                     :: prec, n_tavg

      character(len=20) :: interpolation
      logical :: grid_center

      real(p_double) :: dlb_fld2_thresh, dlb_part_weight
      type(t_species), pointer :: std_species_head

      namelist /nl_particles/ num_pkt, num_species, low_jay_roundoff, &
         ndump_fac, ndump_fac_ave, ndump_fac_lineout, ndump_fac_ene, &
         n_ave, prec, reports, n_tavg, interpolation, grid_center, &
         dlb_fld2_thresh, dlb_part_weight

      integer, dimension(p_n_report_type) :: ndump_fac_all
      integer :: ierr

      character(len=p_max_spname_len) :: spname

      class(t_species), pointer :: species

      ! executable statements

      ! number of standard species
      num_species = 0

      ! number of pkt groups (electrons + positrons + photons)
      num_pkt = 0

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
      call get_namelist(input_file, "nl_particles", ierr)

      if (ierr /= 0) then
         if (ierr < 0) then
            print *, "Error reading particles parameters"
         else
            print *, "Error: particles parameters missing"
         end if
         print *, "aborting..."
         stop
      end if

      read (input_file%nml_text, nml=nl_particles, iostat=ierr)
      if (ierr /= 0) then
         print *, "Error reading particles parameters"
         print *, "aborting..."
         stop
      end if

      ! process interpolation scheme
      select case (trim(interpolation))
      case ("linear")
         this%interpolation = p_linear
      case ("quadratic")
         this%interpolation = p_quadratic
      case ("cubic")
         this%interpolation = p_cubic
      case ("quartic")
         this%interpolation = p_quartic

      case default
         print *, '   Error reading species parameters'
         print *, '   invalid interpolation: "', trim(interpolation), '"'
         print *, '   Valid values are "linear", "quadratic", "cubic" and "quartic".'
         print *, '   aborting...'
         stop
      end select

      this%grid_center = grid_center

      ! process position type
      select case (this%interpolation)
      case (p_linear, p_cubic)
         this%pos_type = p_cell_low
      case (p_quadratic, p_quartic)
         this%pos_type = p_cell_near
      end select

      this%low_jay_roundoff = low_jay_roundoff

      ! Particle kinetic energy diagnostics
      this%ndump_fac_ene = ndump_fac_ene

      ! global charge diagnostics
      ndump_fac_all(p_full) = ndump_fac
      ndump_fac_all(p_savg) = ndump_fac_ave
      ndump_fac_all(p_senv) = ndump_fac_ave
      ndump_fac_all(p_line) = ndump_fac_lineout
      ndump_fac_all(p_slice) = ndump_fac_lineout

      ! process normal reports
      call new(this%reports, reports, p_report_quants, &
               ndump_fac_all, n_ave, n_tavg, prec, &
               p_x_dim, ierr)
      if (ierr /= 0) then
         print *, "(*error*) Invalid report"
         print *, "(*error*) aborting..."
         stop
      end if

      ! pkt specific species and groups

      if (sim_options%omega_p0 <= 0.0) then
         if (mpi_node() == 0) then
            print *, "(*error*) When using pkt algorithm omega_p0 or n0 must be set in the simulation"
            print *, "(*error*) section at the beggining of the input file"
            print *, "(*error*) bailing out..."
         end if
         stop
      end if

      ! no cathodes allowed
      this%num_cathode = 0

      ! no ionization allowed
      this%num_neutral = 0

      ! number of pkt groups
      this%num_pkt = num_pkt

      ! total charged particles species
      this%num_species = num_pkt + num_species

      ! add dlb info
      this%dlb_fld2_thresh = dlb_fld2_thresh
      this%dlb_part_weight = dlb_part_weight

      call this%allocate_objs()

      ! pkt Species
      if (num_pkt > 0) then
         do i = 1, num_pkt
            ! read input of each pkt group
            call this%pkt_group(i)%read_input(input_file, periodic, if_move, grid, &
                                              dt, sim_options)
         end do
      end if

      ! Advance species pointer past pkt electrons
      species => this%species
      do i = 1, num_pkt
         species => species%next
      end do

      if (num_species > 0) then
         do i = 1, num_species
            write (spname, '(A,I0)') 'species ', i
            if (mpi_node() == 0 .and. disp_out(input_file)) then
               print '(A,I0,A)', " - standard species (", i, ") configuration..."
            end if
            call species%read_input(input_file, spname, periodic, if_move, grid, &
                                    dt, .true., sim_options)
            species => species%next
         end do
      end if

      ! validate species names
      call this%validate_names()

      ! #ifdef __HAS_COLLISIONS__
      ! ! Read collision data
      !   call read_nml( this%coll, input_file, this%species, this%num_species, sim_options )
      ! #endif

   end subroutine read_input_part_pkt
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Overrides allocate_objs method and allocates t_species_pkt objects instead
!-----------------------------------------------------------------------------------------
   subroutine allocate_objs_part_pkt(this)

      implicit none

      class(t_particles_pkt), intent(inout) :: this

      integer :: i
      class(t_species), pointer :: spec, tail

      ! pkt groups species
      tail => null()
      if (this%num_pkt > 0) then
         allocate (this%pkt_group(this%num_pkt))

         ! allocate pkt electron and positron species
         do i = 1, this%num_pkt

            ! allocate electrons and photons species
            allocate (t_photons_pkt :: this%pkt_group(i)%photons)
            allocate (t_species_pkt :: this%pkt_group(i)%electrons)

            ! associate photons pointer in electrons
            this%pkt_group(i)%electrons%photons => this%pkt_group(i)%photons

            if (associated(tail)) tail%next => this%pkt_group(i)%electrons

            ! move tail forward to this electrons
            tail => this%pkt_group(i)%electrons

         end do

         ! associate pointers with global particle list
         this%species => this%pkt_group(1)%electrons
      end if

      ! Standard species
      do i = this%num_pkt + 1, this%num_species
         ! Allocate standard species
         allocate (t_species :: spec)
         if (associated(tail)) then
            tail%next => spec
         else
            this%species => spec
         end if
         tail => spec
      end do

   end subroutine allocate_objs_part_pkt
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! sets up this data structure from the given information
!-----------------------------------------------------------------------------------------
   subroutine init_part_pkt(this, g_space, jay, emf, grid, no_co, bnd, zpulse_list, &
                            ndump_fac, restart, restart_handle, t, tstep, tmin, tmax, &
                            sim_options)
!-----------------------------------------------------------------------------------------

      use m_dep_density_pkt, only: deposit_density_1d

      implicit none

      class(t_particles_pkt), intent(inout) :: this

      type(t_space), intent(in) :: g_space
      class(t_emf), intent(inout) :: emf
      class(t_current), intent(inout) :: jay
      class(t_grid), intent(in) :: grid
      class(t_node_conf), intent(in) :: no_co
      class(t_bnd), intent(inout) :: bnd
      class(t_zpulse_list), intent(in) :: zpulse_list
      integer, intent(in) :: ndump_fac
      logical, intent(in) :: restart
      type(t_restart_handle), intent(in) :: restart_handle
      real(p_double), intent(in) :: t
      type(t_time_step), intent(in) :: tstep
      real(p_double), intent(in) :: tmin, tmax
      type(t_options), intent(in) :: sim_options

      integer :: i, size_seed, num_c
      integer, dimension(:), pointer :: seed

      ! call superclass init
      ! this will initialize (and read checkpoint data if required) charged particle species
      call this%t_particles%init(g_space, jay, emf, grid, no_co, bnd, zpulse_list, &
                                 ndump_fac, restart, restart_handle, t, tstep, tmin, &
                                 tmax, sim_options)

      ! diagnostics have an additional setup to handle pkt specific reports
      ! actually not yet...
      ! call setup_diag_pkt( this )

      ! call individual photons init for each pkt group
      if (this%num_pkt > 0) then
         do i = 1, this%num_pkt
            call this%pkt_group(i)%photons%init(i, &
                                                this%interpolation, this%grid_center, grid, g_space, emf, jay, &
                                                no_co, bnd%send_vdf, bnd%recv_vdf, bnd%bnd_cross, bnd%node_cross, bnd%send_spec, &
                                                bnd%recv_spec, ndump_fac, restart, restart_handle, sim_options, tstep, tmin, tmax)
         end do

         ! initialize random number generator
         ! (for performance reasons the pkt code uses the compiler issued random number gen.)

         call random_seed(SIZE=size_seed)
         call alloc(seed, (/size_seed/))
         do i = 1, size_seed
            seed(i) = (no_co%my_aid() + 1)
         end do
         call random_seed(PUT=seed)

         ! no need to keep this around
         call freemem(seed)
      end if

   end subroutine init_part_pkt
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Cleanup particles_pkt object
!-----------------------------------------------------------------------------------------
   subroutine cleanup_part_pkt(this)
!-----------------------------------------------------------------------------------------

      implicit none

      class(t_particles_pkt), intent(inout) :: this

      integer :: i

      ! cleanup superclass
      call this%t_particles%cleanup()

      ! cleanup local data structures
      if (this%num_pkt > 0) then
         ! cleanup photons
         do i = 1, this%num_pkt
            call this%pkt_group(i)%photons%cleanup()
         end do

         deallocate (this%pkt_group)
      end if

   end subroutine cleanup_part_pkt
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Update particle data at the boundaries
!-----------------------------------------------------------------------------------------
   subroutine update_boundary_part_pkt(this, jay, g_space, no_co, dt, bnd)

      implicit none

      class(t_particles_pkt), intent(inout) :: this
      class(t_current), intent(inout) :: jay

      type(t_space), intent(in) :: g_space
      class(t_node_conf), intent(in) :: no_co
      real(p_double), intent(in) :: dt
      class(t_bnd), intent(inout) :: bnd

      integer :: i

      ! update boundary on superclass data (charged particles)
      call this%t_particles%update_boundary(jay, g_space, no_co, dt, bnd)

      ! update boundary on photons, which are not part of this % t_particles % species
      call begin_event(partboundev)
      do i = 1, this%num_pkt
         call this%pkt_group(i)%photons%update_boundary(jay, no_co, dt, &
                                                        bnd%bnd_cross, bnd%node_cross, bnd%send_spec, bnd%recv_spec)
      end do
      call end_event(partboundev)

   end subroutine update_boundary_part_pkt

!-----------------------------------------------------------------------------------------
! Get diagnostic buffer size requirements
!-----------------------------------------------------------------------------------------
   subroutine get_diag_buffer_size_pkt(this, gnx, diag_buffer_size)

      implicit none

      class(t_particles_pkt), intent(in) :: this
      integer, dimension(:), intent(in) :: gnx
      integer, intent(inout) :: diag_buffer_size

      integer :: i

      ! Get dignostic size from superclass
      call this%t_particles%get_diag_buffer_size(gnx, diag_buffer_size)

      ! Also include photon species diagnostics
      do i = 1, this%num_pkt
         call this%pkt_group(i)%photons%get_diag_buffer_size(gnx, diag_buffer_size)
      end do

   end subroutine get_diag_buffer_size_pkt
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Report charged particle and photon data
!-----------------------------------------------------------------------------------------
   subroutine report_part_pkt(this, emf, g_space, grid, no_co, tstep, t, tmin, send_msg, recv_msg)

      implicit none

      class(t_particles_pkt), intent(inout) :: this
      class(t_emf), intent(inout) :: emf
      type(t_space), intent(in) :: g_space
      class(t_grid), intent(in) :: grid
      class(t_node_conf), intent(in) :: no_co
      type(t_time_step), intent(in) :: tstep
      real(p_double), intent(in) :: t, tmin
      type(t_vdf_msg), dimension(2), intent(inout) :: send_msg, recv_msg

      integer :: i

      ! execute superclass diagnostics
      call this%t_particles%report(emf, g_space, grid, no_co, tstep, t, tmin, send_msg, recv_msg)

      ! pkt group diagnostics
      call begin_event(diag_part_ev)
      do i = 1, this%num_pkt
         call this%pkt_group(i)%report(emf, g_space, grid, no_co, tstep, t, tmin, send_msg, recv_msg)
      end do
      call end_event(diag_part_ev)

   end subroutine report_part_pkt
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Report charged particle and photon energy data
!-----------------------------------------------------------------------------------------
   subroutine report_energy_part_pkt(this, no_co, tstep, t, tmin)

      implicit none

      class(t_particles_pkt), intent(in) :: this
      class(t_node_conf), intent(in) :: no_co
      type(t_time_step), intent(in) :: tstep
      real(p_double), intent(in) :: t, tmin

      integer :: i

      ! execute superclass diagnostics
      call this%t_particles%report_energy(no_co, tstep, t, tmin)

      do i = 1, this%num_pkt
         ! report photon energy
         call this%pkt_group(i)%photons%report_energy(no_co, tstep, t, tmin)

      end do

   end subroutine report_energy_part_pkt
!-----------------------------------------------------------------------------------------

   subroutine advance_deposit_pkt(this, emf, jay, tstep, t, no_co, options)

      implicit none

      class(t_particles_pkt), intent(inout) :: this
      class(t_emf), intent(inout) :: emf
      class(t_current), intent(inout) :: jay
      real(p_double), intent(in)    :: t
      type(t_time_step), intent(in)    :: tstep
      class(t_node_conf), intent(in)    :: no_co
      type(t_options), intent(in)    :: options

      integer :: i, j, np_e, np_ph, interp, nx
      type(t_vdf) :: density, grad_grid
      type(t_vdf), target :: photon_grad_grid
      type(t_vdf), target :: photon_density_grid

      real(p_k_part), pointer :: x_e(:, :), x_ph(:, :), q_e(:)
      real(p_k_part), pointer :: p_e(:, :)
      integer, pointer        :: ix_e(:, :), ix_ph(:, :)

      real(p_k_fld), allocatable :: den_ph(:)
      real(p_k_fld), allocatable :: grad_ph(:, :)

      integer, allocatable :: nx_arr(:)
      integer, allocatable :: gc_arr(:, :)
      real(p_double), allocatable :: dx_arr(:)
      type(t_vdf_msg) :: send_msg(2), recv_msg(2)

      real(p_k_part), allocatable :: p_ph_zero(:, :)
      integer :: ic, ilb, iub

      real(p_k_fld), allocatable :: tmp(:)

      character(len=64) :: filename, filename2
      logical :: if_feedback

      !--------------------------------------------------
      ! Loop over PKT groups
      !--------------------------------------------------
      do i = 1, this%num_pkt

         np_e = this%pkt_group(i)%electrons%num_par
         np_ph = this%pkt_group(i)%photons%num_par
         interp = this%pkt_group(i)%electrons%interpolation
         nx = this%pkt_group(i)%electrons%g_nx(1)
         if_feedback = (this%pkt_group(i)%electrons%push_type == p_std_feedback)

         ! Allocate useful arrays
         allocate (nx_arr(p_x_dim), dx_arr(p_x_dim))
         allocate (gc_arr(2, p_x_dim))
         allocate (p_ph_zero(3, np_ph))
         p_ph_zero = 0.0_p_k_part

         nx_arr(:) = this%pkt_group(i)%electrons%g_nx(1:p_x_dim)
         dx_arr(:) = this%pkt_group(i)%electrons%dx(1:p_x_dim)
         gc_arr(:, :) = interp

         ! Create a new density array by copying rho from the electrons
         call density%new(jay%pf(1), f_dim=1, zero=.true.)
         call photon_density_grid%new(jay%pf(1), f_dim=1, zero=.true.)
         call photon_grad_grid%new(jay%pf(1), f_dim=p_x_dim, zero=.true.)

         select case (this%pkt_group(i)%photons%coordinates)
         case (p_cylindrical_b)
            call grad_grid%new(jay%pf(1), f_dim=2, zero=.true.)
         case default
            call grad_grid%new(jay%pf(1), f_dim=p_x_dim, zero=.true.)
         end select

         ! Allocate particle arrays
         allocate (x_e(p_x_dim, np_e), ix_e(p_x_dim, np_e), q_e(np_e))
         allocate (x_ph(p_x_dim, np_ph), ix_ph(p_x_dim, np_ph))
         allocate (den_ph(np_ph), grad_ph(p_x_dim, np_ph))

         ! Fill electron data
         x_e(:, :) = this%pkt_group(i)%electrons%x(1:p_x_dim, 1:np_e)
         ix_e(:, :) = this%pkt_group(i)%electrons%ix(1:p_x_dim, 1:np_e)
         q_e(:) = this%pkt_group(i)%electrons%q(1:np_e)
         p_e => this%pkt_group(i)%electrons%p(1:3, 1:np_e)

         ! Fill photon data
         x_ph(:, :) = this%pkt_group(i)%photons%x(1:p_x_dim, 1:np_ph)
         ix_ph(:, :) = this%pkt_group(i)%photons%ix(1:p_x_dim, 1:np_ph)

         ! CALL SUBROUTINES BASED ON DIMENSION OF SIMULATION
         select case (this%pkt_group(i)%photons%coordinates)
         case default
            select case (p_x_dim)
            case (1)

               ! ---- 1. Deposit PHOTON density FIRST (electrons need it below) ----
               call deposit_density_1d(photon_density_grid, ix_ph, x_ph, &
                                       this%pkt_group(i)%photons%q(1:np_ph), &
                                       p_ph_zero, np_ph, interp)
               call update_boundary(photon_density_grid, p_vdf_add, no_co, send_msg, recv_msg)
               call cleanup(send_msg)
               call cleanup(recv_msg)
               this%pkt_group(i)%electrons%photon_density_grid => photon_density_grid

               ! --- smoothing: photon_density_grid, feeds electron gamma correction ---
               ilb = lbound(photon_density_grid%f1, 2) + 1
               iub = ubound(photon_density_grid%f1, 2) - 1
               if (allocated(tmp)) deallocate (tmp)
               allocate (tmp(lbound(photon_density_grid%f1, 2):ubound(photon_density_grid%f1, 2)))
               do j = 1, 30
                  tmp(ilb:iub) = 0.25_p_k_fld*photon_density_grid%f1(1, ilb - 1:iub - 1) + &
                                 0.50_p_k_fld*photon_density_grid%f1(1, ilb:iub) + &
                                 0.25_p_k_fld*photon_density_grid%f1(1, ilb + 1:iub + 1)
                  photon_density_grid%f1(1, ilb:iub) = tmp(ilb:iub)
               end do
               deallocate (tmp)

               ! ---- 2. Deposit ELECTRON density, now WITH the <A^2> correction ----
               if (if_feedback) then
                  call deposit_density_1d(density, ix_e, x_e, q_e, p_e, np_e, interp, &
                                          photon_density_grid=photon_density_grid)
               else
                  call deposit_density_1d(density, ix_e, x_e, q_e, p_e, np_e, interp)
               end if
               call update_boundary(density, p_vdf_add, no_co, send_msg, recv_msg)
               call cleanup(send_msg)
               call cleanup(recv_msg)

               ! --- smoothing: density, feeds photon pusher ---
               ilb = lbound(density%f1, 2) + 1
               iub = ubound(density%f1, 2) - 1
               if (allocated(tmp)) deallocate (tmp)
               allocate (tmp(lbound(density%f1, 2):ubound(density%f1, 2)))
               do j = 1, 30
                  tmp(ilb:iub) = 0.25_p_k_fld*density%f1(1, ilb - 1:iub - 1) + &
                                 0.50_p_k_fld*density%f1(1, ilb:iub) + &
                                 0.25_p_k_fld*density%f1(1, ilb + 1:iub + 1)
                  density%f1(1, ilb:iub) = tmp(ilb:iub)
               end do
               deallocate (tmp)

               ! ---- 3. photon_grad_grid (unchanged -- still needs density%dx(1)) ----
               call dep_den_grad_grid_1d(photon_grad_grid, ix_ph, x_ph, &
                                         this%pkt_group(i)%photons%q(1:np_ph), &
                                         p_ph_zero, np_ph, density%dx(1))
               call update_boundary(photon_grad_grid, p_vdf_add, no_co, send_msg, recv_msg)
               call cleanup(send_msg)
               call cleanup(recv_msg)
               this%pkt_group(i)%electrons%photon_gradient_grid => photon_grad_grid

               ! --- smoothing: photon_grad_grid, feeds electron feedback (P0) ---
               ilb = lbound(photon_grad_grid%f1, 2) + 1
               iub = ubound(photon_grad_grid%f1, 2) - 1
               if (allocated(tmp)) deallocate (tmp)
               allocate (tmp(lbound(photon_grad_grid%f1, 2):ubound(photon_grad_grid%f1, 2)))
               do j = 1, 30
                  tmp(ilb:iub) = 0.25_p_k_fld*photon_grad_grid%f1(1, ilb - 1:iub - 1) + &
                                 0.50_p_k_fld*photon_grad_grid%f1(1, ilb:iub) + &
                                 0.25_p_k_fld*photon_grad_grid%f1(1, ilb + 1:iub + 1)
                  photon_grad_grid%f1(1, ilb:iub) = tmp(ilb:iub)
               end do
               deallocate (tmp)

               ! ---- 4. Gather electron density onto photon positions (den_ph) ----
               call gather_density_1d(den_ph, density, ix_ph, x_ph, np_ph, interp)

               ! ---- 5. Electron GRADIENT deposit, also needs the correction ----
               if (if_feedback) then
                  call dep_den_grad_grid_1d(grad_grid, ix_e, x_e, q_e, p_e, np_e, density%dx(1), &
                                            photon_density_grid=photon_density_grid)
               else
                  call dep_den_grad_grid_1d(grad_grid, ix_e, x_e, q_e, p_e, np_e, density%dx(1))
               end if
               call update_boundary(grad_grid, p_vdf_add, no_co, send_msg, recv_msg)
               call cleanup(send_msg)
               call cleanup(recv_msg)

               ! --- smoothing: grad_grid, feeds photon push ---
               ilb = lbound(grad_grid%f1, 2) + 1
               iub = ubound(grad_grid%f1, 2) - 1
               if (allocated(tmp)) deallocate (tmp)
               allocate (tmp(lbound(grad_grid%f1, 2):ubound(grad_grid%f1, 2)))
               do j = 1, 30
                  tmp(ilb:iub) = 0.25_p_k_fld*grad_grid%f1(1, ilb - 1:iub - 1) + &
                                 0.50_p_k_fld*grad_grid%f1(1, ilb:iub) + &
                                 0.25_p_k_fld*grad_grid%f1(1, ilb + 1:iub + 1)
                  grad_grid%f1(1, ilb:iub) = tmp(ilb:iub)
               end do
               deallocate (tmp)

               call gather_density_1d(grad_ph(1, :), grad_grid, ix_ph, x_ph, np_ph, interp)

               ! -- ADVANCE PHOTONS --------------------------------------------
               call this%pkt_group(i)%photons%advance(den_ph, grad_ph, density, &
                                                      grad_grid, dt(tstep), 0, 1)

            case (2)
               call deposit_density_2d(density, ix_e, x_e, q_e, p_e, np_e, interp)
               call gather_density_2d(den_ph, density, ix_ph, x_ph, np_ph, interp)
               call gather_density_grad_2d(grad_ph, density, ix_ph, x_ph, np_ph, interp)

               ! -- ADVANCE PHOTONS --------------------------------------------
               call this%pkt_group(i)%photons%advance(den_ph, grad_ph, density, &
                                                      grad_grid, dt(tstep), 0, 1)

            case (3)
               call deposit_density_3d(density, ix_e, x_e, q_e, np_e, interp)
               call gather_density_3d(den_ph, density, ix_ph, x_ph, np_ph, interp)
               call gather_density_grad_3d(grad_ph, density, ix_ph, x_ph, np_ph, interp)

               ! -- ADVANCE PHOTONS --------------------------------------------
               call this%pkt_group(i)%photons%advance(den_ph, grad_ph, density, &
                                                      grad_grid, dt(tstep), 0, 1)
            case default
               ERROR('Not implemented for x_dim = ', p_x_dim)
               call abort_program(p_err_invalid)
            end select
         case (p_cylindrical_b)
            call deposit_density_2d(density, ix_e, x_e, q_e, p_e, np_e, interp)

            call update_boundary(density, p_vdf_add, no_co, send_msg, recv_msg)
            call cleanup(send_msg)
            call cleanup(recv_msg)

            call gather_density_2d(den_ph, density, ix_ph, x_ph, np_ph, interp)

            ! Deposit both gradient components
            call dep_den_grad_grid_2d(grad_grid, ix_e, x_e, q_e, p_e, np_e, density%dx(1), density%dx(2))

            call update_boundary(grad_grid, p_vdf_add, no_co, send_msg, recv_msg)
            call cleanup(send_msg)
            call cleanup(recv_msg)

            call gather_density_2d(grad_ph(1, :), grad_grid, ix_ph, x_ph, np_ph, interp)
            call gather_density_2d(grad_ph(2, :), grad_grid, ix_ph, x_ph, np_ph, interp, field_comp=2)

            ! -- ADVANCE PHOTONS -----------------------------------------------
            call this%pkt_group(i)%photons%advance(den_ph, grad_ph, density, &
                                                   grad_grid, dt(tstep), 0, 1)

         end select

         ! Cleanup
         call density%cleanup()
         call grad_grid%cleanup()
         deallocate (x_e, ix_e, q_e)
         deallocate (x_ph, ix_ph, den_ph, grad_ph)
         deallocate (nx_arr, dx_arr, gc_arr)
         deallocate (p_ph_zero)

      end do

      ! PUSH ALL REMAINING PARTICLES
      call this%t_particles%advance_deposit(emf, jay, tstep, t, no_co, options)
      call photon_grad_grid%cleanup()

   end subroutine advance_deposit_pkt

!-------------------------------------------------------------------------------
! Printout the algorithm used by each species pusher
!-------------------------------------------------------------------------------
   subroutine list_algorithm_part_pkt(this)

      implicit none
      class(t_particles_pkt), intent(in) :: this
      integer :: i
      class(t_species), pointer :: species

      print *, ' '
      print '(A)', 'Particles:'

      print '(A,I0)', '- Number of pkt groups: ', this%num_pkt

      if (this%num_pkt > 0) then

         ! particle shape
         select case (this%interpolation)
         case (p_linear)
            print '(A)', '- Linear (1st order) interpolation'
         case (p_quadratic)
            print '(A)', '- Quadratic (2nd order) interpolation'
         case (p_cubic)
            print '(A)', '- Cubic (3rd order) interpolation'
         case (p_quartic)
            print '(A)', '- Quartic (4th order) interpolation'
         end select

      end if

      if (this%grid_center) then
         print '(A)', '- Using spatially centered grid values for interpolation'
      else
         print '(A)', '- Using staggered Yee grid values for interpolation'
      end if

      ! Current deposition
      if (this%low_jay_roundoff) then
         print *, '- Depositing current using low roundoff algorithm'
      end if

      do i = 1, this%num_pkt
         print *, ' '
         print '(A,I0)', '- pkt group #', i
         call this%pkt_group(i)%electrons%list_algorithm()
         call this%pkt_group(i)%photons%list_algorithm()
      end do

      ! List algorithm for standard (non pkt) species
      if (this%num_pkt > 0) then
         species => this%pkt_group(this%num_pkt)%electrons%next
      else
         species => this%species
      end if

      if (associated(species)) then
         print *, ' '
         print '(A)', '- Standard species'
      end if

      do
         if (.not. associated(species)) exit
         call species%list_algorithm()
         species => species%next
      end do

   end subroutine list_algorithm_part_pkt
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
   subroutine write_checkpoint_part_pkt(this, restart_handle)

      implicit none

      class(t_particles_pkt), intent(in) :: this
      type(t_restart_handle), intent(inout) :: restart_handle

      integer :: i

      ! Call superclass method to save charged particle checkpoint data
      call this%t_particles%write_checkpoint(restart_handle)

      ! Save photon data
      do i = 1, this%num_pkt
         call this%pkt_group(i)%photons%write_checkpoint(restart_handle)
      end do

   end subroutine write_checkpoint_part_pkt
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
   subroutine move_window_part_pkt(this, g_space, grid, jay, no_co, bnd)

      implicit none

      class(t_particles_pkt), intent(inout) :: this

      type(t_space), intent(in) :: g_space
      class(t_grid), intent(in) :: grid
      class(t_current), intent(inout)   :: jay
      class(t_node_conf), intent(in)     :: no_co
      class(t_bnd), intent(inout) :: bnd

      integer :: i

      ! Call superclass method to move charged particle species
      call this%t_particles%move_window(g_space, grid, jay, no_co, bnd)

      ! Move photons in pkt groups
      do i = 1, this%num_pkt
         call move_window(this%pkt_group(i)%photons, g_space, grid, jay, no_co, &
                          bnd%bnd_cross, bnd%node_cross, bnd%send_spec, bnd%recv_spec)
      end do

   end subroutine move_window_part_pkt
!-------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! Sets the particle load.
! For iteration n = 0 the load is set for the whole simulation volume (since this is called
! when no node partitions exist yet), for n > 0 the load is set for the
! local volume/particles only
!---------------------------------------------------------------------------------------------------
   subroutine add_load_part_pkt(this, grid, emf, n)

      implicit none

      class(t_particles_pkt), intent(inout) :: this
      class(t_grid), intent(inout) :: grid
      class(t_emf), intent(in) :: emf
      integer, intent(in) :: n

      class(t_species), pointer :: species
      integer :: i

      ! loop through all species and add to int_load array
      species => this%species
      if (n > 0) then

         do
            if (.not. associated(species)) exit
            call species%add_particle_load(grid, emf)
            species => species%next
         end do

         ! Add load for photons
         do i = 1, this%num_pkt
            call this%pkt_group(i)%photons%add_particle_load(grid, emf)
         end do

      else

         do
            if (.not. associated(species)) exit
            call add_density_load(species, grid)
            species => species%next
         end do

      end if

   end subroutine add_load_part_pkt
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
   subroutine reshape_part_pkt(this, old_lb, new_lb, no_co, send_msg, recv_msg)
      !-----------------------------------------------------------------------------------------
      ! redistribute particles to new simulation partition
      !-----------------------------------------------------------------------------------------

      implicit none

      class(t_particles_pkt), intent(inout) :: this
      class(t_grid), intent(in) :: new_lb, old_lb
      class(t_node_conf), intent(in) :: no_co
      type(t_vdf_msg), dimension(2), intent(inout) :: send_msg, recv_msg

      type(t_msg_patt) :: msg_patt
      class(t_species), pointer :: species
      integer :: n, i
      integer, dimension(2, p_max_dim) :: gc_num
      type(t_vdf_report), pointer :: report

      ! reshape species
      if (this%num_species > 0) then
         ! get message pattern (guard cells are not required since no particles
         ! should be in guard cells). This is the same for all species
         ! so we only need to do it once.
         gc_num = 0
         call new_msg_patt(msg_patt, old_lb, new_lb, no_co, gc_num)

         ! loop through all species
         species => this%species
         do
            if (.not. associated(species)) exit
            call species%reshape(old_lb, new_lb, msg_patt, no_co, send_msg, recv_msg)
            species => species%next
         end do

         ! reshape the photon species
         do i = 1, this%num_pkt
            call this%pkt_group(i)%photons%reshape(old_lb, new_lb, msg_patt, no_co, &
                                                   send_msg, recv_msg)
         end do

         ! clear message pattern data
         call msg_patt%cleanup()
      end if

      ! reshape local vdf objects if needed
      if (this%jay_tmp%x_dim_ > 0) call reshape_nocopy(this%jay_tmp, new_lb)

      call reshape(this%charge, old_lb, new_lb, no_co, send_msg, recv_msg)

      ! Reshape tavg_data
      report => this%reports
      do
         if (.not. associated(report)) exit
         if (report%tavg_data%x_dim_ > 0) then
            call reshape_copy(report%tavg_data, old_lb, new_lb, no_co, send_msg, recv_msg)
         end if
         report => report%next
      end do

   end subroutine reshape_part_pkt
!-----------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! Reports number of particles per node for all nodes (small grid file)
!---------------------------------------------------------------------------------------------------
   subroutine report_node_load_part_pkt(this, n, ndump, t, no_co)

      implicit none

      class(t_particles), intent(in) :: this
      integer, intent(in)             :: n, ndump
      real(p_double), intent(in)    :: t
      class(t_node_conf), intent(in) :: no_co

      real(p_double), dimension(:), pointer     :: node_load
      real(p_double), dimension(:), pointer     :: f1
      real(p_double), dimension(:, :), pointer   :: f2
      real(p_double), dimension(:, :, :), pointer :: f3

      class(t_diag_file), allocatable :: diagFile

      integer, dimension(p_max_dim) :: lnx

      class(t_species), pointer :: species
      integer :: i, npart, nnodes, ierr
      real(p_double) :: npart_dbl

      node_load => null()
      f1 => null()
      f2 => null()
      f3 => null()

      ! get total number of particles on node
      npart = 0
      species => this%species
      do
         if (.not. associated(species)) exit
         npart = npart + species%num_par
         species => species%next
      end do

      ! Add photons
      select type (this)
      class is (t_particles_pkt)
         do i = 1, this%num_pkt
            npart = npart + this%pkt_group(i)%photons%num_par
         end do
      end select

      npart_dbl = npart

      ! gather data
      nnodes = no_num(no_co)
      call alloc(node_load, (/nnodes/))

      if (nnodes > 1) then
         call mpi_gather(npart_dbl, 1, MPI_DOUBLE_PRECISION, &
                         node_load, 1, MPI_DOUBLE_PRECISION, &
                         0, comm(no_co), ierr)
         if (ierr /= MPI_SUCCESS) then
            ERROR('MPI Error')
            call abort_program(p_err_mpi)
         end if
      else
         node_load = npart_dbl
      end if

      ! save data to disk
      if (root(no_co)) then

         lnx(1:p_x_dim) = nx(no_co)

         ! Open the output file
         call create_diag_file(diagFile)
         diagFile%ftype = p_diag_grid

         diagFile%filepath = trim(path_mass)//'LOAD'//p_dir_sep//'NODE'//p_dir_sep
         diagFile%filename = trim(get_filename(n/ndump, 'node_load'))
         diagFile%name = 'Particles per node'
         diagFile%iter%n = n
         diagFile%iter%t = t
         diagFile%iter%time_units = '1 / \omega_p'

         diagFile%grid%ndims = p_x_dim
         diagFile%grid%name = "load"
         diagFile%grid%label = "particles per node"
         diagFile%grid%units = "particles"

         do i = 1, p_x_dim
            diagFile%grid%axis(i)%min = 0
            diagFile%grid%axis(i)%max = lnx(i)
            write (diagFile%grid%axis(i)%name, '(A,I0)') 'x', i
            write (diagFile%grid%axis(i)%label, '(A,I0)') 'x_', i
            diagFile%grid%axis(i)%units = 'cell'
         end do

         call diagFile%open(p_diag_create)

         ! Move node load values to the proper position and write node load values
         select case (p_x_dim)
         case (1)
            call alloc(f1, lnx)
            f1 = -1.0
            do i = 1, nnodes
               f1(no_co%ngp(i, 1)) = node_load(i)
            end do
            call diagFile%add_dataset("load", f1)
            call freemem(f1)

         case (2)
            call alloc(f2, lnx)
            f2 = -1.0
            do i = 1, nnodes
               f2(no_co%ngp(i, 1), no_co%ngp(i, 2)) = node_load(i)
            end do
            call diagFile%add_dataset("load", f2)
            call freemem(f2)

         case (3)
            call alloc(f3, lnx)
            f3 = -1.0
            do i = 1, nnodes
               f3(no_co%ngp(i, 1), no_co%ngp(i, 2), no_co%ngp(i, 3)) = node_load(i)
            end do
            call diagFile%add_dataset("load", f3)
            call freemem(f3)
         end select

         ! close output file
         call diagFile%close()
         deallocate (diagFile)
      end if

      ! free load array
      call freemem(node_load)

   end subroutine report_node_load_part_pkt
!---------------------------------------------------------------------------------------------------

end module m_particles_pkt
