# 1 "qed-cyl/os-photons-define-qed-cyl.f03"
# 1 "<built-in>" 1
# 1 "<built-in>" 3
# 467 "<built-in>" 3
# 1 "<command line>" 1
# 1 "<built-in>" 2
# 1 "qed-cyl/os-photons-define-qed-cyl.f03" 2
!
! The current version does not allow for normal species, only QEDCYL groups
!

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
# 6 "qed-cyl/os-photons-define-qed-cyl.f03" 2
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
# 7 "qed-cyl/os-photons-define-qed-cyl.f03" 2

module m_photons_define_qedcyl

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
# 11 "qed-cyl/os-photons-define-qed-cyl.f03" 2

use m_system
use m_parameters
use m_species_define
use m_species_define_qed, only: add_particle_load_qed
use m_restart
use m_particles_cyl_modes, only: t_particles_cyl_modes
use m_species_cyl_modes_define, only: t_species_cyl_modes
use m_species_define_qedcyl, only: t_species_qedcyl
use m_input_file, only: t_input_file, get_namelist
use m_vdf_report, only: new
use m_species_collisions, only: read_nml
use m_grid_define, only: t_grid, t_msg_patt
use m_emf_define, only: t_emf
use m_node_conf, only: t_node_conf
use m_current_define, only: t_current
use m_bnd, only: t_bnd
use m_vdf_comm, only : t_vdf_msg
use m_vdf_define
use m_emf_cyl_modes, only: t_emf_cyl_modes
use m_space
use m_time_step
use m_betheheitler_qedcyl, only: t_betheheitler_info_qedcyl
use m_cyl_modes, only: t_cyl_modes, setup

implicit none

private

! string to id restart data
character(len=*), parameter :: p_photons_rst_id = "photons rst data - 0x0000"


!-----------------------------------------------------------------------------------------
! t_photons_qedcyl
! - Photon particle species
!-----------------------------------------------------------------------------------------
type, extends( t_species_cyl_modes ) :: t_photons_qedcyl

  real(p_k_part) :: omega_p0 !reference electron plasma frequency

  ! number of pairs created
  integer :: num_pairs = 0
  real( p_double ) :: weight_pairs = 0.0
  real( p_double ) :: energy_pairs = 0.0

  integer :: num_pairs_bh = 0
  real( p_double ) :: weight_pairs_bh = 0.0
  real( p_double ) :: energy_pairs_bh = 0.0

  ! Additional diagnostics
  integer :: ndump_fac_pairs

  ! pair production flag
  logical :: if_pairprod

  real( p_double ) :: pp_chi_cutoff

  ! emission events diagnostic
  integer :: ndump_fac_n_emit
  type( t_vdf ) :: n_emit

  ! chi parameter diagnostic
  integer :: ndump_fac_chi
  type( t_vdf ) :: weight, weight_chi

  ! Bethe Heitler pair creation
  type( t_betheheitler_info_qedcyl ) :: betheheitler_info

  ! load balance parameters, obtained from particles object
  real(p_double) :: dlb_fld2_thresh
  real(p_double) :: dlb_part_weight

contains

  procedure :: read_input => read_input_photons_qedcyl
  procedure :: init => init_photons_qedcyl
  procedure :: list_algorithm => list_algorithm_photons_qedcyl
  procedure :: report => report_photons_qedcyl
  procedure :: update_boundary => update_boundary_photons_qedcyl
  procedure :: report_energy => report_energy_photons_qedcyl
  procedure :: report_pairs => report_pairs_photons_qedcyl
  procedure :: get_energy => get_energy_photons_qedcyl
  procedure :: get_energy_pairs => get_energy_pairs_qedcyl
  procedure :: advance => advance_photons_qedcyl
  procedure :: inject => inject_photons_qedcyl
  procedure :: validate => validate_photons_qedcyl
  procedure :: read_checkpoint => read_checkpoint_photons_qedcyl
  procedure :: write_checkpoint => write_checkpoint_photons_qedcyl
  procedure :: reshape => reshape_photons_qedcyl
  procedure :: add_particle_load => add_particle_load_photons_qedcyl

end type t_photons_qedcyl

interface
  subroutine advance_photons_qedcyl( this, electrons, positrons, emf, gdt, tid, n_threads )

    import t_photons_qedcyl, t_species_qedcyl, t_emf_cyl_modes, p_double

    class( t_photons_qedcyl ), intent(inout) :: this
    class( t_species_qedcyl ), intent(inout) :: electrons
    class( t_species_qedcyl ), intent(inout) :: positrons
    class( t_emf_cyl_modes ), intent( inout ) :: emf
    real(p_double), intent(in) :: gdt
    integer, intent(in) :: tid ! local thread id
    integer, intent(in) :: n_threads ! total number of threads

  end subroutine
end interface

interface
  subroutine report_photons_qedcyl( this, emf, g_space, grid, no_co, tstep, t, tmin, send_msg, recv_msg )

    import t_photons_qedcyl, t_emf, t_space, t_grid, t_node_conf, t_time_step, p_double, t_vdf_msg

    class( t_photons_qedcyl ), intent(inout) :: this

    class( t_emf ), intent(inout) :: emf
    type( t_space ), intent(in) :: g_space
    class( t_grid ), intent(in) :: grid
    class( t_node_conf ), intent(in) :: no_co
    type( t_time_step ), intent(in) :: tstep
    real(p_double), intent(in) :: t, tmin
    type(t_vdf_msg), dimension(2), intent(inout) :: send_msg, recv_msg

  end subroutine
end interface

interface
  subroutine report_energy_photons_qedcyl( this, no_co, tstep, t, tmin )

    import t_photons_qedcyl, t_node_conf, t_time_step, p_double

    class( t_photons_qedcyl ), intent(in) :: this

    class( t_node_conf ), intent(in) :: no_co
    type( t_time_step ), intent(in) :: tstep
    real(p_double), intent(in) :: t, tmin

  end subroutine
end interface

interface
  subroutine report_pairs_photons_qedcyl( this, no_co, tstep, t, tmin )

    import t_photons_qedcyl, t_node_conf, t_time_step, p_double

    class( t_photons_qedcyl ), intent(in) :: this

    class( t_node_conf ), intent(in) :: no_co
    type( t_time_step ), intent(in) :: tstep
    real(p_double), intent(in) :: t, tmin

  end subroutine
end interface

public :: t_photons_qedcyl

contains

!-----------------------------------------------------------------------------------------
! Read photon species information from input file
!-----------------------------------------------------------------------------------------
subroutine read_input_photons_qedcyl( this, input_file, def_name, periodic, if_move, grid, &
  dt, read_prof, sim_options )

    use m_system
    use m_parameters

    use m_qed
    use m_input_file
    use m_species_define
    use m_grid_define

    use m_profile, only : new_source

    implicit none

    class( t_photons_qedcyl ), intent(inout) :: this
    class( t_input_file ), intent(inout) :: input_file

    character(len = *), intent(in) :: def_name
    logical, dimension(:), intent(in) :: periodic, if_move
    class( t_grid ), intent(in) :: grid

    real(p_double), intent(in) :: dt
    logical, intent(in) :: read_prof

    type(t_options), intent(in) :: sim_options

    integer :: i, j

    ! call superclass read_input
    call this % t_species_cyl_modes % read_input( input_file, def_name, periodic, if_move, grid, &
                            dt, read_prof, sim_options )

    ! reset particle source for now
    !this % source => null()

    ! make sure that boundary conditions are valid
    do i = 1, p_x_dim
        do j = 1, 2
            if (this%bnd_con%type(j,i) == p_bc_thermal) then
                write(err_buf__,*) "Invalid boundary conditions for photon species.";call err__("qed-cyl/os-photons-define-qed-cyl.f03",214)
                call abort_program( p_err_invalid )
            endif
        enddo
    enddo

    ! set the omega_p0 parameter
    this % omega_p0 = real( sim_options % omega_p0, p_k_part )

end subroutine read_input_photons_qedcyl
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Initialize object data
!-----------------------------------------------------------------------------------------
subroutine init_photons_qedcyl( this, sp_id, interpolation, grid_center, grid, g_space, emf, jay, &
  no_co, send_vdf, recv_vdf, bnd_cross, node_cross, send_spec, &
  recv_spec, ndump_fac, restart, restart_handle, sim_options, tstep, tmin, tmax )

    use m_grid_define, only : t_grid
    use m_space, only : t_space
    use m_emf_define, only : t_emf
    use m_current_define, only : t_current
    use m_node_conf

    use m_restart

    use m_species_define
    use m_system, only : t_options
    use m_vdf_comm, only : t_vdf_msg

    implicit none

    ! dummy variables
    class( t_photons_qedcyl ), intent(inout) :: this

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
    integer :: n_modes

    n_modes = this % n_cyl_modes

    ! call superclass initialization
    call this % t_species_cyl_modes % init( sp_id, interpolation, grid_center, grid, g_space, emf, jay, &
        no_co, send_vdf, recv_vdf, bnd_cross, node_cross, send_spec, &
        recv_spec, ndump_fac, restart, restart_handle, sim_options, tstep, tmin, tmax )

    ! initially no pairs have been created
    this%num_pairs = 0
    this%weight_pairs = 0.0
    this%energy_pairs = 0.0

    this%num_pairs_bh = 0
    this%weight_pairs_bh = 0.0
    this%energy_pairs_bh = 0.0

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

    ! setup vdf with ion density for Bethe Heitler
    if (this%betheheitler_info%if_betheheitler) then
        ! The number of guard cells required depends on the interpolation level
        gc_num(p_lower,:) = this%interpolation
        gc_num(p_upper,:) = this%interpolation + 1
        dx = 0.
        call this % betheheitler_info % rho_ion % new( p_x_dim, 1, this%my_nx_p(3,1:p_x_dim), gc_num, dx, .true. )
        call setup( this % betheheitler_info % rho_ion_cyl_m, n_modes, this % betheheitler_info % rho_ion )
    endif

    ! setup local data
    if ( restart ) then
        call this % read_checkpoint( restart_handle )
    endif

end subroutine init_photons_qedcyl
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Printout the algorithm used by the pusher
!-----------------------------------------------------------------------------------------
subroutine list_algorithm_photons_qedcyl( this )

  implicit none
  class( t_photons_qedcyl ), intent(in) :: this

  if (mpi_node()==0) print *,"Species : ", trim(this%name)
  print *, '   - Photon QED pusher'


end subroutine list_algorithm_photons_qedcyl
!-----------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine update_boundary_photons_qedcyl( this, jay, no_co, dt, bnd_cross, node_cross, send_msg, recv_msg )

  use m_node_conf
  use m_current_define

  implicit none

  ! photons
  class( t_photons_qedcyl ), intent(inout) :: this
  class( t_current ), intent(inout) :: jay
  class( t_node_conf ), intent(in) :: no_co
  real(p_double), intent(in) :: dt
  type( t_part_idx ), dimension(2), intent(inout) :: bnd_cross
  type( t_part_idx ), intent(inout) :: node_cross
  type( t_spec_msg ), dimension(2), intent(inout) :: send_msg, recv_msg

  ! Since only p_bc_thermal deposits current and photons cannot have thermal bath boundaries
  ! the jay object will not be altered and could be replaced by a dummy one
  call this % t_species_cyl_modes % update_boundary( jay, no_co, dt, bnd_cross, node_cross, send_msg, recv_msg )

end subroutine update_boundary_photons_qedcyl
!---------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Return photon energy on local node, calculating it if need be
!-----------------------------------------------------------------------------------------
function get_energy_photons_qedcyl( this )

  use m_math, only : pi

  implicit none

  class( t_photons_qedcyl ), intent(in) :: this
  real( p_double ) :: get_energy_photons_qedcyl

  integer :: i
  real( p_double ) :: energy, u2, kin, cell_volume

  ! Check if energy was calculated during the push (time-centered), otherwise
  ! calculate it here
  if ( this % energy(1) < 0 ) then
     ! Energy not available, recalculate it
     energy = 0

     !$omp parallel do private(u2, kin) reduction(+ : energy)
     do i = 1, this%num_par
        u2 = this%p(1,i)**2 + this%p(2,i)**2 + this%p(3,i)**2
        kin = sqrt( u2 )
        energy = energy + this%q(i) * kin
     enddo
     !$omp end parallel do

     ! Store value so it can be reused
     this%energy(1) = energy

     do i = 2, ubound(this % energy, 1)
        this % energy(i) = 0
     enddo


  else

     ! Get energy calculated previously
     energy = this % energy(1)


     ! If doing a multi-threaded run, add contribution from all threads
     do i = 2, ubound(this % energy, 1)
        energy = energy + this % energy(i)
     enddo


  endif

  ! Normalize to cell size and charge over mass ratio
  cell_volume = this % dx(1)
  do i = 2, p_x_dim
    cell_volume = cell_volume * this % dx(i)
  enddo

  cell_volume = cell_volume * 2.d0 * pi

  get_energy_photons_qedcyl = energy * cell_volume

end function get_energy_photons_qedcyl
!---------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Return pair energy on local node
!-----------------------------------------------------------------------------------------
function get_energy_pairs_qedcyl( this, bethe_heitler )

  use m_math, only : pi

  implicit none

  class( t_photons_qedcyl ), intent(in) :: this
  logical, intent(in) :: bethe_heitler
  real( p_double ) :: get_energy_pairs_qedcyl

  integer :: i
  real( p_double ) :: energy, cell_volume

  ! Get energy calculated previously
  energy = this % energy_pairs
  if (bethe_heitler) energy = this % energy_pairs_bh

  ! Normalize to cell size and charge over mass ratio
  cell_volume = this % dx(1)
  do i = 2, p_x_dim
    cell_volume = cell_volume * this % dx(i)
  enddo

  get_energy_pairs_qedcyl = energy * cell_volume * 2.d0 * pi

end function get_energy_pairs_qedcyl
!---------------------------------------------------------------------------------------------------

subroutine inject_photons_qedcyl( this, ig_xbnd_inj, jay, no_co, bnd_cross, node_cross, send_msg, recv_msg )

  use m_current_define

  implicit none

  class( t_photons_qedcyl ), intent(inout) :: this
  integer, dimension(:, :), intent(in) :: ig_xbnd_inj
  class( t_current ), intent(inout) :: jay
  class( t_node_conf ), intent(in) :: no_co
  type( t_part_idx ), dimension(2), intent(inout) :: bnd_cross
  type( t_part_idx ), intent(inout) :: node_cross
  type( t_spec_msg ), dimension(2), intent(inout) :: send_msg, recv_msg

  ! Photons injection is not yet implemented, they are only created from QED processes
  ! This method is called from the move_window routine

  continue

end subroutine inject_photons_qedcyl

!-----------------------------------------------------------------------------------------
! Checks if all particle values are ok.
!-----------------------------------------------------------------------------------------
  subroutine validate_photons_qedcyl( this, msg, over )

    use m_system
    use m_parameters

    implicit none

    ! dummy variables

    class(t_photons_qedcyl), intent(in) :: this
    character( len = * ) , intent(in) :: msg
    logical, intent(in), optional :: over

    ! local variables
    integer :: i_dim, k
    integer, dimension(p_x_dim) :: ilb, iub
    logical :: over_

    ! executable statements

    if ( present( over) ) then
      over_ = over
    else
      over_ = .false.
    endif

    ! get grid boundaries
    ilb = 1
    iub = this%t_species%my_nx_p(3,1:p_x_dim)

    ! allow for 1 cell overflow (ok before update boundary)
    if ( over_ ) then
      ilb = ilb - 1
      iub = iub + 1
    endif


    ! validate positions
    do i_dim = 1, p_x_dim
    do k = 1, this%t_species%num_par
       if ( ( this%x(i_dim, k) < -0.5_p_k_part ) .or. &
        ( this%x(i_dim, k) >= 0.5_p_k_part ) .or. &
        ( this%ix(i_dim, k) < ilb( i_dim ) ) .or. &
        ( this%ix(i_dim, k) > iub( i_dim) ) ) then

            call bad_photon( k, this, ilb, iub, msg // " - Invalid position " )

       endif
    enddo
    enddo

    ! validate momenta
    do k = 1, this%t_species%num_par
      do i_dim = 1, p_p_dim
         if ( isinf( this%p(i_dim, k) ) .or. isnan( this%p(i_dim, k) ) ) then

            call bad_photon( k, this, ilb, iub, msg // " - Invalid momenta " )

         endif
      enddo
    enddo

    ! validate weight
    do k = 1, this%t_species%num_par
     if ( this%q(k) <= 0 .or. isinf( this%q(k) ) .or. isnan( this%q(k) ) ) then

          call bad_photon( k, this, ilb, iub, msg // " - Invalid weight " )

     endif
    enddo

  contains

  subroutine bad_photon( k, self, ilb, iub, msg )

    implicit none

    integer, intent(in) :: k

    class( t_photons_qedcyl ), intent(in) :: self
    integer, dimension(:), intent(in) :: ilb, iub
    character( len = * ), intent(in) :: msg

    write(0,'(A,I0,A,A)') '[', mpi_node(), '] ', trim(msg)

    write(0,'(A,I0,A,I0,A,I0)') "[", mpi_node(), "] Bad photon ", k, " of ", self%t_species%num_par
    write(0,*) "[", mpi_node(), "] p (:)  =", self%p(:, k)
    write(0,*) "[", mpi_node(), "] x (:)  =", self%x(:, k)
    write(0,*) "[", mpi_node(), "] ix (:) =", self%ix(:, k)
    write(0,*) "[", mpi_node(), "] w      =", self%q(k)

    write(0,*) "[", mpi_node(), "] ilb(:) = ", ilb
    write(0,*) "[", mpi_node(), "] iub(:) = ", iub

    write(0,'(A,I0,A,A)') '[', mpi_node(), '] (* error *) Validade photons failed for ', &
                                      trim(self%t_species%name), ' aborting...'

    call abort_program( p_err_invalid )


  end subroutine bad_photon

  end subroutine validate_photons_qedcyl
  !-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
! Read checkpoint data
!-----------------------------------------------------------------------------------------
  subroutine read_checkpoint_photons_qedcyl( this, restart_handle )

    implicit none

    class( t_photons_qedcyl ), intent(inout) :: this
    type( t_restart_handle ), intent(in) :: restart_handle

    character(len=*), parameter :: err_msg = 'error reading restart data for photons object.'
    character(len=len(p_photons_rst_id)) :: rst_id
    integer :: ierr

    call restart_io_read(rst_id, restart_handle, ierr)
    call check_error(ierr,err_msg,p_err_rstrd,"qed-cyl/os-photons-define-qed-cyl.f03",604)

    ! check if restart file is compatible
    if ( rst_id /= p_photons_rst_id) then
      write(err_buf__,*) 'Corrupted restart file, or restart file ';call err__("qed-cyl/os-photons-define-qed-cyl.f03",608)
      write(err_buf__,*) 'from incompatible binary (photons)';call err__("qed-cyl/os-photons-define-qed-cyl.f03",609)
      call abort_program(p_err_rstrd)
    endif

    ! read checkpoint data
    call restart_io_read(this%num_pairs, restart_handle, ierr)
    call check_error(ierr,err_msg,p_err_rstrd,"qed-cyl/os-photons-define-qed-cyl.f03",615)

    call restart_io_read(this%weight_pairs, restart_handle, ierr)
    call check_error(ierr,err_msg,p_err_rstrd,"qed-cyl/os-photons-define-qed-cyl.f03",618)

    call restart_io_read(this%energy_pairs, restart_handle, ierr)
    call check_error(ierr,err_msg,p_err_rstrd,"qed-cyl/os-photons-define-qed-cyl.f03",621)

    call restart_io_read(this%num_pairs_bh, restart_handle, ierr)
    call check_error(ierr,err_msg,p_err_rstrd,"qed-cyl/os-photons-define-qed-cyl.f03",624)

    call restart_io_read(this%weight_pairs_bh, restart_handle, ierr)
    call check_error(ierr,err_msg,p_err_rstrd,"qed-cyl/os-photons-define-qed-cyl.f03",627)

    call restart_io_read(this%energy_pairs_bh, restart_handle, ierr)
    call check_error(ierr,err_msg,p_err_rstrd,"qed-cyl/os-photons-define-qed-cyl.f03",630)

  end subroutine read_checkpoint_photons_qedcyl

!-----------------------------------------------------------------------------------------
  ! Write checkpoint data
  !-----------------------------------------------------------------------------------------
  subroutine write_checkpoint_photons_qedcyl( this, restart_handle )

    use m_species

    implicit none

    class( t_photons_qedcyl ), intent(in) :: this
    type( t_restart_handle ), intent(inout) :: restart_handle

    character(len=*), parameter :: err_msg = 'error writing restart data for photons object.'
    integer :: ierr

    ! write superclass data
    call this % t_species_cyl_modes % restart_write( restart_handle )

    call restart_io_write("p_photons_rst_id", p_photons_rst_id, restart_handle, ierr)
    call check_error(ierr,err_msg,p_err_rstwrt,"qed-cyl/os-photons-define-qed-cyl.f03",653)

    ! write checkpoint data
    call restart_io_write("this%num_pairs", this%num_pairs, restart_handle, ierr)
    call check_error(ierr,err_msg,p_err_rstwrt,"qed-cyl/os-photons-define-qed-cyl.f03",657)

    call restart_io_write("this%weight_pairs", this%weight_pairs, restart_handle, ierr)
    call check_error(ierr,err_msg,p_err_rstwrt,"qed-cyl/os-photons-define-qed-cyl.f03",660)

    call restart_io_write("this%energy_pairs", this%energy_pairs, restart_handle, ierr)
    call check_error(ierr,err_msg,p_err_rstwrt,"qed-cyl/os-photons-define-qed-cyl.f03",663)

    call restart_io_write("this%num_pairs_bh", this%num_pairs_bh, restart_handle, ierr)
    call check_error(ierr,err_msg,p_err_rstwrt,"qed-cyl/os-photons-define-qed-cyl.f03",666)

    call restart_io_write("this%weight_pairs_bh", this%weight_pairs_bh, restart_handle, ierr)
    call check_error(ierr,err_msg,p_err_rstwrt,"qed-cyl/os-photons-define-qed-cyl.f03",669)

    call restart_io_write("this%energy_pairs_bh", this%energy_pairs_bh, restart_handle, ierr)
    call check_error(ierr,err_msg,p_err_rstwrt,"qed-cyl/os-photons-define-qed-cyl.f03",672)

  end subroutine write_checkpoint_photons_qedcyl

  !-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! redistribute the species particle data through all nodes when node grids change and
! store new grid boundaries
!-----------------------------------------------------------------------------------------
subroutine reshape_photons_qedcyl( this, old_grid, new_grid, msg_patt, no_co, send_msg, recv_msg )
!-----------------------------------------------------------------------------------------

  use m_vdf_comm, only : reshape_copy

  implicit none

  class( t_photons_qedcyl ), intent(inout) :: this
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

end subroutine reshape_photons_qedcyl
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Wrapper for load function
!-----------------------------------------------------------------------------------------
subroutine add_particle_load_photons_qedcyl( this, grid, emf )

  implicit none

  class( t_photons_qedcyl ), intent(in) :: this
  class( t_grid ), intent(inout) :: grid
  class( t_emf ), intent(in) :: emf

  call add_particle_load_qed( this, grid, emf, this%omega_p0, this%dlb_fld2_thresh, &
                              this%dlb_part_weight )

end subroutine add_particle_load_photons_qedcyl
!-----------------------------------------------------------------------------------------

end module m_photons_define_qedcyl
