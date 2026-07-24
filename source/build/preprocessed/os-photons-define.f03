# 1 "qed/os-photons-define.f03"
# 1 "<built-in>" 1
# 1 "<built-in>" 3
# 467 "<built-in>" 3
# 1 "<command line>" 1
# 1 "<built-in>" 2
# 1 "qed/os-photons-define.f03" 2
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
# 2 "qed/os-photons-define.f03" 2
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
# 3 "qed/os-photons-define.f03" 2

module m_photons_define

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
# 7 "qed/os-photons-define.f03" 2

!-----------------------------------------------------------------------------------------

use m_system
use m_parameters
use m_species_define, only: t_species, t_part_idx, t_spec_msg
use m_species_define_qed, only: t_species_qed
use m_psource_std
use m_species_udist
use m_mergedel, only: t_merge_info

use m_space
use m_grid_define
use m_node_conf
use m_restart
use m_time_step

use m_emf_define
use m_vdf_define
use m_vdf_comm, only : t_vdf_msg
use m_input_file
use m_betheheitler

implicit none

private

! string to id restart data
character(len=*), parameter :: p_photons_rst_id = "photons rst data - 0x0000"


!-----------------------------------------------------------------------------------------
! t_photons
! - Photon particle species with merging
!-----------------------------------------------------------------------------------------
type, extends( t_species_qed ) :: t_photons

  ! number of pairs created
  integer :: num_pairs = 0
  real( p_double ) :: weight_pairs = 0.0_p_double
  real( p_double ) :: energy_pairs = 0.0_p_double

  integer :: num_pairs_bh = 0
  real( p_double ) :: weight_pairs_bh = 0.0_p_double
  real( p_double ) :: energy_pairs_bh = 0.0_p_double

  ! Additional diagnostics
  integer :: ndump_fac_pairs

  ! pair production flag
  logical :: if_pairprod

  real( p_double ) :: pp_chi_cutoff

  ! Bethe Heitler pair creation
  type( t_betheheitler_info ) :: betheheitler_info

contains
  procedure :: read_input => read_input_photons
  procedure :: list_algorithm => list_algorithm_photons
  procedure :: init => init_photons
  procedure :: advance => advance_photons
  procedure :: update_boundary => update_boundary_photons
  procedure :: inject => inject_photons
  procedure :: validate => validate_photons
  procedure :: report => report_photons
  procedure :: report_energy => report_energy_photons
  procedure :: report_pairs => report_pairs_photons
  procedure :: get_energy => get_energy_photons
  procedure :: get_energy_pairs => get_energy_pairs
  procedure :: mergedel => mergedel_photons

  procedure :: read_checkpoint => read_checkpoint_photons
  procedure :: write_checkpoint => write_checkpoint_photons

end type t_photons

integer, parameter :: p_phot_dens = 1
integer, parameter :: p_phot_ene = 1

! type bound procedure (method) interfaces

interface
  subroutine advance_photons( this, electrons, positrons, emf, gdt, tid, n_threads )

    import t_photons, t_species_qed, t_emf, p_double

    class( t_photons ), intent(inout) :: this
    class( t_species_qed ), intent(inout) :: electrons
    class( t_species_qed ), intent(inout) :: positrons
    class( t_emf ), intent( inout ) :: emf
    real(p_double), intent(in) :: gdt
    integer, intent(in) :: tid ! local thread id
    integer, intent(in) :: n_threads ! total number of threads

  end subroutine
end interface

interface
  subroutine report_photons( this, emf, g_space, grid, no_co, tstep, t, tmin, send_msg, recv_msg )

    import t_photons, t_emf, t_space, t_grid, t_node_conf, t_time_step, p_double, t_vdf_msg

    class( t_photons ), intent(inout) :: this

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
  subroutine report_energy_photons( this, no_co, tstep, t, tmin )

    import t_photons, t_node_conf, t_time_step, p_double

    class( t_photons ), intent(in) :: this

    class( t_node_conf ), intent(in) :: no_co
    type( t_time_step ), intent(in) :: tstep
    real(p_double), intent(in) :: t, tmin

  end subroutine
end interface

interface
  subroutine report_pairs_photons( this, no_co, tstep, t, tmin )

    import t_photons, t_node_conf, t_time_step, p_double

    class( t_photons ), intent(in) :: this

    class( t_node_conf ), intent(in) :: no_co
    type( t_time_step ), intent(in) :: tstep
    real(p_double), intent(in) :: t, tmin

  end subroutine
end interface

public :: t_photons

contains

!-----------------------------------------------------------------------------------------
! Read photon species information from input file
!-----------------------------------------------------------------------------------------
subroutine read_input_photons( this, input_file, def_name, periodic, if_move, grid, &
  dt, read_prof, sim_options )

    use m_system
    use m_parameters

    use m_qed
    use m_input_file
    use m_species_define
    use m_grid_define

    use m_profile, only : new_source

    implicit none

    class( t_photons ), intent(inout) :: this
    class( t_input_file ), intent(inout) :: input_file

    character(len = *), intent(in) :: def_name
    logical, dimension(:), intent(in) :: periodic, if_move
    class( t_grid ), intent(in) :: grid

    real(p_double), intent(in) :: dt
    logical, intent(in) :: read_prof

    type(t_options), intent(in) :: sim_options

    integer :: i, j

    ! call superclass read_input
    call this % t_species % read_input( input_file, def_name, periodic, if_move, grid, &
                            dt, read_prof, sim_options )

    ! @todo: ensure that only basic particle sources are allowed
    ! this % source => null()

    ! make sure that boundary conditions are valid
    do i = 1, p_x_dim
        do j = 1, 2
            if (this%bnd_con%type(j,i) == p_bc_thermal) then
                write(err_buf__,*) "Invalid boundary conditions for photon species.";call err__("qed/os-photons-define.f03",198)
                call abort_program( p_err_invalid )
            endif
        enddo
    enddo

    ! set the omega_p0 and n0 variables
    this % omega_p0 = real( sim_options % omega_p0, p_k_part )
    this % n0 = real( sim_options % n0, p_k_part )

    ! read additional merge information
    call this % merge_info % read_input( input_file )

end subroutine read_input_photons
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Printout the algorithm used by the pusher
!-----------------------------------------------------------------------------------------
subroutine list_algorithm_photons( this )

  implicit none
  class( t_photons ), intent(in) :: this

  if (mpi_node()==0) print *,"Species : ", trim(this%name)
  print *, '   - Photon QED pusher'


end subroutine list_algorithm_photons
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Initialize object data
!-----------------------------------------------------------------------------------------
subroutine init_photons( this, sp_id, interpolation, grid_center, grid, g_space, emf, jay, &
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
    class( t_photons ), intent(inout) :: this

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

    ! call superclass initialization
    call this % t_species % init( sp_id, interpolation, grid_center, grid, g_space, emf, jay, &
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
    endif

    ! setup local data
    if ( restart ) then
        call this % read_checkpoint( restart_handle )
    endif

end subroutine init_photons
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Return photon energy on local node, calculating it if need be
!-----------------------------------------------------------------------------------------
function get_energy_photons( this )

  implicit none

  class( t_photons ), intent(in) :: this
  real( p_double ) :: get_energy_photons

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

  get_energy_photons = energy * cell_volume

end function get_energy_photons
!---------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Return pair energy on local node
!-----------------------------------------------------------------------------------------
function get_energy_pairs( this, bethe_heitler )

  implicit none

  class( t_photons ), intent(in) :: this
  logical, intent(in) :: bethe_heitler
  real( p_double ) :: get_energy_pairs

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

  get_energy_pairs = energy * cell_volume

end function get_energy_pairs
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine update_boundary_photons( this, jay, no_co, dt, bnd_cross, node_cross, send_msg, recv_msg )

  use m_node_conf
  use m_current_define

  implicit none

  ! photons
  class( t_photons ), intent(inout) :: this
  class( t_current ), intent(inout) :: jay
  class( t_node_conf ), intent(in) :: no_co
  real(p_double), intent(in) :: dt
  type( t_part_idx ), dimension(2), intent(inout) :: bnd_cross
  type( t_part_idx ), intent(inout) :: node_cross
  type( t_spec_msg ), dimension(2), intent(inout) :: send_msg, recv_msg

  ! Since only p_bc_thermal deposits current and photons cannot have thermal bath boundaries
  ! the jay object will not be altered and could be replaced by a dummy one
  call this % t_species % update_boundary( jay, no_co, dt, bnd_cross, node_cross, send_msg, recv_msg )

end subroutine update_boundary_photons
!---------------------------------------------------------------------------------------------------

subroutine inject_photons( this, ig_xbnd_inj, jay, no_co, bnd_cross, node_cross, send_msg, recv_msg )

  use m_current_define

  implicit none

  class( t_photons ), intent(inout) :: this
  integer, dimension(:, :), intent(in) :: ig_xbnd_inj
  class( t_current ), intent(inout) :: jay
  class( t_node_conf ), intent(in) :: no_co
  type( t_part_idx ), dimension(2), intent(inout) :: bnd_cross
  type( t_part_idx ), intent(inout) :: node_cross
  type( t_spec_msg ), dimension(2), intent(inout) :: send_msg, recv_msg

  ! Photons injection is not yet implemented, they are only created from QED processes
  ! This method is called from the move_window routine

  continue

end subroutine inject_photons
# 488 "qed/os-photons-define.f03"
!-----------------------------------------------------------------------------------------
! Merge photon particles
!-----------------------------------------------------------------------------------------

  subroutine mergedel_photons( this, n, nxc, t, no_co )

    use m_system
    use m_node_conf
    use m_mergedel

    implicit none


    class( t_photons ), intent(inout) :: this
    integer, intent(in) :: n
    integer, dimension(:), intent(in) :: nxc
    real(p_double), intent(in) :: t
    class( t_node_conf ), intent(in) :: no_co

    ! merge photons
    call mergedel( this, this%merge_info, n, nxc, t, no_co, p_photon )

  end subroutine mergedel_photons
  !-----------------------------------------------------------------------------------------
# 523 "qed/os-photons-define.f03"
!-----------------------------------------------------------------------------------------
! Checks if all particle values are ok.
!-----------------------------------------------------------------------------------------
  subroutine validate_photons( this, msg, over )

    use m_system
    use m_parameters

    implicit none

    ! dummy variables

    class(t_photons), intent(in) :: this
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
    ilb = 1 - this%move_num(1:p_x_dim)
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

            print *, "[",mpi_node(),"] ",'over = ', over_
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

    class( t_photons ), intent(in) :: self
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

  end subroutine validate_photons
  !-----------------------------------------------------------------------------------------
# 646 "qed/os-photons-define.f03"
!-----------------------------------------------------------------------------------------
! Read checkpoint data
!-----------------------------------------------------------------------------------------
  subroutine read_checkpoint_photons( this, restart_handle )

    implicit none

    class( t_photons ), intent(inout) :: this
    type( t_restart_handle ), intent(in) :: restart_handle

    character(len=*), parameter :: err_msg = 'error reading restart data for photons object.'
    character(len=len(p_photons_rst_id)) :: rst_id
    integer :: ierr

    call restart_io_read(rst_id, restart_handle, ierr)
    call check_error(ierr,err_msg,p_err_rstrd,"qed/os-photons-define.f03",661)

    ! check if restart file is compatible
    if ( rst_id /= p_photons_rst_id) then
      write(err_buf__,*) 'Corrupted restart file, or restart file ';call err__("qed/os-photons-define.f03",665)
      write(err_buf__,*) 'from incompatible binary (photons)';call err__("qed/os-photons-define.f03",666)
      call abort_program(p_err_rstrd)
    endif

    ! read checkpoint data
    call restart_io_read(this%num_pairs, restart_handle, ierr)
    call check_error(ierr,err_msg,p_err_rstrd,"qed/os-photons-define.f03",672)

    call restart_io_read(this%weight_pairs, restart_handle, ierr)
    call check_error(ierr,err_msg,p_err_rstrd,"qed/os-photons-define.f03",675)

    call restart_io_read(this%energy_pairs, restart_handle, ierr)
    call check_error(ierr,err_msg,p_err_rstrd,"qed/os-photons-define.f03",678)

    call restart_io_read(this%num_pairs_bh, restart_handle, ierr)
    call check_error(ierr,err_msg,p_err_rstrd,"qed/os-photons-define.f03",681)

    call restart_io_read(this%weight_pairs_bh, restart_handle, ierr)
    call check_error(ierr,err_msg,p_err_rstrd,"qed/os-photons-define.f03",684)

    call restart_io_read(this%energy_pairs_bh, restart_handle, ierr)
    call check_error(ierr,err_msg,p_err_rstrd,"qed/os-photons-define.f03",687)

  end subroutine read_checkpoint_photons

  !-----------------------------------------------------------------------------------------
  ! Write checkpoint data
  !-----------------------------------------------------------------------------------------
  subroutine write_checkpoint_photons( this, restart_handle )

    use m_species

    implicit none

    class( t_photons ), intent(in) :: this
    type( t_restart_handle ), intent(inout) :: restart_handle

    character(len=*), parameter :: err_msg = 'error writing restart data for photons object.'
    integer :: ierr

    ! write superclass data
    call this % t_species % restart_write( restart_handle )

    call restart_io_write("p_photons_rst_id", p_photons_rst_id, restart_handle, ierr)
    call check_error(ierr,err_msg,p_err_rstwrt,"qed/os-photons-define.f03",710)

    ! write checkpoint data
    call restart_io_write("this%num_pairs", this%num_pairs, restart_handle, ierr)
    call check_error(ierr,err_msg,p_err_rstwrt,"qed/os-photons-define.f03",714)

    call restart_io_write("this%weight_pairs", this%weight_pairs, restart_handle, ierr)
    call check_error(ierr,err_msg,p_err_rstwrt,"qed/os-photons-define.f03",717)

    call restart_io_write("this%energy_pairs", this%energy_pairs, restart_handle, ierr)
    call check_error(ierr,err_msg,p_err_rstwrt,"qed/os-photons-define.f03",720)

    call restart_io_write("this%num_pairs_bh", this%num_pairs_bh, restart_handle, ierr)
    call check_error(ierr,err_msg,p_err_rstwrt,"qed/os-photons-define.f03",723)

    call restart_io_write("this%weight_pairs_bh", this%weight_pairs_bh, restart_handle, ierr)
    call check_error(ierr,err_msg,p_err_rstwrt,"qed/os-photons-define.f03",726)

    call restart_io_write("this%energy_pairs_bh", this%energy_pairs_bh, restart_handle, ierr)
    call check_error(ierr,err_msg,p_err_rstwrt,"qed/os-photons-define.f03",729)

  end subroutine write_checkpoint_photons

  !-----------------------------------------------------------------------------------------





end module m_photons_define
