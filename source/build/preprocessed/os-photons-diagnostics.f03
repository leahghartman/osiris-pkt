# 1 "qed/os-photons-diagnostics.f03"
# 1 "<built-in>" 1
# 1 "<built-in>" 3
# 467 "<built-in>" 3
# 1 "<command line>" 1
# 1 "<built-in>" 2
# 1 "qed/os-photons-diagnostics.f03" 2
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
# 2 "qed/os-photons-diagnostics.f03" 2
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
# 3 "qed/os-photons-diagnostics.f03" 2

module m_photons_diagnostics

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
# 7 "qed/os-photons-diagnostics.f03" 2

use m_system
use m_parameters
use m_math

use m_node_conf
use m_time_step
use m_space, only : t_space
use m_grid_define, only : t_grid
use m_vdf_define
use m_diagfile, only : p_time_length

use m_photons_define, only: t_photons

implicit none

private

interface n_emit_hist
    module procedure n_emit_hist
end interface

interface report_n_emit
    module procedure report_n_emit
end interface

interface chi_hist
    module procedure chi_hist
end interface

interface report_chi
    module procedure report_chi
end interface

public :: n_emit_hist, report_n_emit, chi_hist, report_chi

contains

subroutine n_emit_hist(this, ix)

    implicit none

    class(t_photons), intent(inout) :: this
    integer, dimension(p_x_dim), intent(in) :: ix
    integer :: n_emit

    select case ( p_x_dim )
    case ( 1 )
        n_emit = this%n_emit%f1(1, ix(1))
        this%n_emit%f1(1, ix(1)) = n_emit + 1
    case ( 2 )
        n_emit = this%n_emit%f2(1, ix(1), ix(2))
        this%n_emit%f2(1, ix(1), ix(2)) = n_emit + 1
    case ( 3 )
        n_emit = this%n_emit%f3(1, ix(1), ix(2), ix(3))
        this%n_emit%f3(1, ix(1), ix(2), ix(3)) = n_emit + 1
    end select

end subroutine n_emit_hist

subroutine report_n_emit( this, g_space, grid, no_co, tstep, t )

    use m_diagnostic_utilities
    use stringutil, only: idx_string

    implicit none

    class(t_photons), intent(inout) :: this
    type( t_space ), intent(in) :: g_space
    class( t_grid ), intent(in) :: grid
    class( t_node_conf ), intent(in) :: no_co
    type( t_time_step ), intent(in) :: tstep
    real(p_double), intent(in) :: t

    type( t_vdf_report ) :: info
    integer :: ns, nd

    ! radiated energy spectrum diagnostic
    if (test_if_report( tstep, this%ndump_fac_n_emit )) then

        ns = n(tstep)
        nd = ndump(tstep) * this%ndump_fac_n_emit

        ! Prepare the metadata
        info%name = 'n_emit'
        info%ndump = nd

        info%xname = (/'x1', 'x2', 'x3'/)

        info%xlabel = (/'x_1', 'x_2', 'x_3'/)
        info%xunits = (/'c / \omega_p', 'c / \omega_p', 'c / \omega_p'/)

        info%time_units = '1 / \omega_p'
        info%dt = dt(tstep)

        info%label = 'N'
        info%units = ''

        info%n = ns
        info%t = t

        info%fileLabel = ''
        info%path = trim(path_mass) // 'QED' // p_dir_sep // 'n_emit' // &
                     p_dir_sep // trim(this%name) // p_dir_sep
        info%filename = 'n_emit-' // trim(this%name) // '-' // idx_string( ns/nd, p_time_length )

        ! write the file
        call this % n_emit % write( info, 1, g_space, grid, no_co )

        call this % n_emit % zero()
    endif
end subroutine report_n_emit

subroutine chi_hist(this, ix, q, chi)

    implicit none

    class(t_photons), intent(inout) :: this
    integer, dimension(p_x_dim), intent(in) :: ix
    real( p_k_part ), intent(in) :: q, chi
    real( p_k_part ) :: var_01, var_02

    select case ( p_x_dim )

    case ( 1 )

        var_01 = this%weight%f1(1, ix(1))
        this%weight%f1(1, ix(1)) = var_01 + abs( q )
        var_02 = this%weight_chi%f1(1, ix(1))
        this%weight_chi%f1(1, ix(1)) = var_02 + abs( q ) * chi

    case ( 2 )

        var_01 = this%weight%f2(1, ix(1), ix(2))
        this%weight%f2(1, ix(1), ix(2)) = var_01 + abs( q )
        var_02 = this%weight_chi%f2(1, ix(1), ix(2))
        this%weight_chi%f2(1, ix(1), ix(2)) = var_02 + abs( q ) * chi

    case ( 3 )

        var_01 = this%weight%f3(1, ix(1), ix(2), ix(3))
        this%weight%f3(1, ix(1), ix(2), ix(3)) = var_01 + abs( q )
        var_02 = this%weight_chi%f3(1, ix(1), ix(2), ix(3))
        this%weight_chi%f3(1, ix(1), ix(2), ix(3)) = var_02 + abs( q ) * chi

    end select

end subroutine chi_hist

subroutine report_chi( this, g_space, grid, no_co, tstep, t )

    use m_diagnostic_utilities
    use stringutil, only: idx_string
    use m_vdf_math

    implicit none

    class(t_photons), intent(inout) :: this
    type( t_space ), intent(in) :: g_space
    class( t_grid ), intent(in) :: grid
    class( t_node_conf ), intent(in) :: no_co
    type( t_time_step ), intent(in) :: tstep
    real(p_double), intent(in) :: t

    type( t_vdf_report ) :: info
    integer :: ns, nd

    ! radiated energy spectrum diagnostic
    if (test_if_report( tstep, this%ndump_fac_chi )) then

        ns = n(tstep)
        nd = ndump(tstep) * this%ndump_fac_chi

        ! Prepare the metadata
        info%name = 'chi'
        info%ndump = nd

        info%xname = (/'x1', 'x2', 'x3'/)

        info%xlabel = (/'x_1', 'x_2', 'x_3'/)
        info%xunits = (/'c / \omega_p', 'c / \omega_p', 'c / \omega_p'/)

        info%time_units = '1 / \omega_p'
        info%dt = dt(tstep)

        info%label = '\chi'
        info%units = ''

        info%n = ns
        info%t = t

        info%fileLabel = ''
        info%path = trim(path_mass) // 'QED' // p_dir_sep // 'chi_grid' // &
                     p_dir_sep // trim(this%name) // p_dir_sep
        info%filename = 'chi-grid-' // trim(this%name) // '-' // idx_string( ns/nd, p_time_length )

        ! We divided the two vdfs to get the weight-averaged chi parameter in each cell
        call divis( this % weight_chi, this % weight )

        ! write the file
        call this % weight_chi % write( info, 1, g_space, grid, no_co )

    endif

    ! We set it to zero after each time step
    ! whether we wrote it or not
    if ( associated( this % weight % buffer ) ) then
        call this % weight_chi % zero()
        call this % weight % zero()
    endif

end subroutine report_chi

end module m_photons_diagnostics

!------------------------------------------------------------------------------------------
subroutine report_photons( this, emf, g_space, grid, no_co, tstep, t, tmin, send_msg, recv_msg )
!-----------------------------------------------------------------------------------------

    use m_system
    use m_parameters
    use m_photons_define, only: t_photons
    use m_emf_define, only : t_emf
    use m_space, only : t_space
    use m_grid_define, only : t_grid
    use m_node_conf
    use m_time_step
    use m_vdf_comm, only : t_vdf_msg
    use m_species_diagnostics, only: deposit_quant, p_charge
    use m_photons_diagnostics

    implicit none

    ! dummy variables

    class( t_photons ), intent(inout) :: this
    class( t_emf ), intent(inout) :: emf
    type( t_space ), intent(in) :: g_space
    class( t_grid ), intent(in) :: grid
    class( t_node_conf ), intent(in) :: no_co
    type( t_time_step ), intent(in) :: tstep
    real(p_double), intent(in) :: t, tmin
    type(t_vdf_msg), dimension(2), intent(inout) :: send_msg, recv_msg

    ! process superclass diagnostics
    call this % t_species % report( emf, g_space, grid, no_co, tstep, t, tmin, send_msg, recv_msg )

    ! process subclass diagnostics
    call report_n_emit( this, g_space, grid, no_co, tstep, t)
    call report_chi( this, g_space, grid, no_co, tstep, t)

    ! Interpolation of ion density for Bethe Heitler
    ! It is done here because all the variables needed are available
    ! deposit_quant includes the deposit + the update_boundary.
    if (this%betheheitler_info%if_betheheitler ) then
        call this % betheheitler_info % rho_ion % cleanup()
        call deposit_quant( this%betheheitler_info%ion, grid, no_co, this % betheheitler_info % rho_ion, p_charge, send_msg, recv_msg )
    end if

end subroutine report_photons
!-----------------------------------------------------------------------------------------

!------------------------------------------------------------------------------------------
subroutine report_energy_photons( this, no_co, tstep, t, tmin )
!-----------------------------------------------------------------------------------------

    use m_system
    use m_photons_define, only: t_photons
    use m_node_conf
    use m_units
    use m_time_step
    use m_parameters, only: file_id_parene, p_err_diagfile, path_hist
    use stringutil, only: idx_string

    implicit none

    ! dummy variables
    class( t_photons ), intent(in) :: this
    class( t_node_conf ), intent(in) :: no_co
    type( t_time_step ), intent(in) :: tstep
    real(p_double), intent(in) :: t, tmin

    character(80) :: path, full_name
    integer :: ierr
    integer( p_int64 ) :: total_par
    real(p_double) :: energy

    ! Particle energy diagnostic
    if (test_if_report( tstep, this%diag%ndump_fac_ene )) then

        ! Get energy on local node
        energy = this % get_energy( )

        ! sums the results from all nodes
        call reduce(no_co, energy, operation = p_sum)

        ! get the total number of particles in the simulation box
        total_par = this%num_par
        call reduce( no_co, total_par, operation = p_sum)

        if ( root(no_co) ) then

          ! File path
          path = trim(path_hist)

          ! file name including path
          write( full_name, '(A,A,I0.2,A)' ) trim(path), 'pho', this%sp_id, '_ene'


          ! Open file and position at the last record
          if ( t == tmin ) then
              call mkdir( path, ierr )
              call check_error(ierr,"Unable to create directory for particle energy diagnostic",p_err_diagfile,"qed/os-photons-diagnostics.f03",319)

              open (unit=file_id_parene, file=full_name, status = 'REPLACE' , &
              form='formatted')

              ! Write Header
              write(file_id_parene, '(A,A,A)' ) '! Kinetic energy for species "', trim( this%name ), '"'
              write(file_id_parene,'( A6, 2(1X,A15), (1X,A23) )') &
                  'Iter','Time    ','Total Par.','Kin. Energy'

          else
              open (unit=file_id_parene, file=full_name, position = 'append', &
              form='formatted')
          endif

          ! Writes the particle energies to file
          write(file_id_parene, '( I6,1X,g15.8,1X,I14,1(1X,es23.16) )') &
                  n(tstep), t, total_par, energy
          close(file_id_parene)

        endif

    endif

end subroutine report_energy_photons
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Report on number of pairs created
!-----------------------------------------------------------------------------------------
subroutine report_pairs_photons( this, no_co, tstep, t, tmin )
!-----------------------------------------------------------------------------------------

    use m_system
    use m_photons_define, only: t_photons
    use m_node_conf
    use m_units
    use m_time_step
    use m_parameters, only: file_id_parene, p_err_diagfile, path_hist
    use stringutil, only: idx_string

    implicit none

    ! dummy variables
    class( t_photons ), intent(in) :: this
    class( t_node_conf ), intent(in) :: no_co
    type( t_time_step ), intent(in) :: tstep
    real(p_double), intent(in) :: t, tmin

    ! local variables
    character(80) :: path, full_name
    character(2) :: ch_sp_id

    integer :: ierr
    integer( p_int64 ) :: total_pairs, total_pairs_bh
    real( p_double ):: weight_pairs, weight_pairs_bh
    real( p_double ) :: energy_pairs, energy_pairs_bh


    ! Pair creation diagnostics
    if (test_if_report( tstep, this%ndump_fac_pairs )) then

        ! get the total number of particles in the simulation box
        total_pairs = this%num_pairs
        call reduce( no_co, total_pairs, operation = p_sum)

        weight_pairs = this%weight_pairs
        call reduce( no_co, weight_pairs, operation = p_sum)

        energy_pairs = this % get_energy_pairs(.false.)
        call reduce( no_co, energy_pairs, operation = p_sum)

        total_pairs_bh = this%num_pairs_bh
        call reduce( no_co, total_pairs_bh, operation = p_sum)

        weight_pairs_bh = this%weight_pairs_bh
        call reduce( no_co, weight_pairs_bh, operation = p_sum)

        energy_pairs_bh = this % get_energy_pairs(.true.)
        call reduce( no_co, energy_pairs_bh, operation = p_sum)

        if ( no_co%my_aid() == 1 ) then

        ! prepare path and file names
        ch_sp_id = idx_string( this%sp_id, 2 )

        path = trim(path_hist)

        full_name = trim(path) // 'pho' // trim(ch_sp_id) // '_pairs'

        ! Open file and position at the last record
        if ( t == tmin ) then
            call mkdir( path, ierr )

            open (unit=file_id_parene, file=full_name, status = 'REPLACE' , &
            form='formatted')

            ! Write Header
            write(file_id_parene,'( A6, (1X, A15), 6(3X, A20) )') &
                'Iter','Time    ','Total Pairs (BW)', 'Total Weight (BW)', 'Total Energy (BW)','Total Pairs (BH)', 'Total Weight (BH)', 'Total Energy (BH)'

        else
            open (unit=file_id_parene, file=full_name, position = 'append', &
            form='formatted')
        endif

        ! Writes the particle energies to file
        write(file_id_parene, '( I6,1X,g15.8,1X,I14,es28.16,es28.16,I14,es28.16,es28.16 )') n(tstep), t, total_pairs, &
                        weight_pairs, energy_pairs, total_pairs_bh, weight_pairs_bh, energy_pairs_bh
        close(file_id_parene)

        endif

    endif

end subroutine report_pairs_photons
!---------------------------------------------------------------------------------------------------
