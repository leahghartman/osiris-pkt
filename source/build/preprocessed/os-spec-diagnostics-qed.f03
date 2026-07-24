# 1 "qed/os-spec-diagnostics-qed.f03"
# 1 "<built-in>" 1
# 1 "<built-in>" 3
# 467 "<built-in>" 3
# 1 "<command line>" 1
# 1 "<built-in>" 2
# 1 "qed/os-spec-diagnostics-qed.f03" 2
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
# 2 "qed/os-spec-diagnostics-qed.f03" 2
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
# 3 "qed/os-spec-diagnostics-qed.f03" 2

module m_species_diagnostics_qed

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
# 7 "qed/os-spec-diagnostics-qed.f03" 2

use m_system
use m_parameters
use m_math
use m_node_conf
use m_time_step

use m_emf_define, only : t_emf
use m_space, only : t_space
use m_grid_define, only : t_grid
use m_vdf_define
use m_species_define_qed, only: t_species_qed
use m_diagfile, only: p_time_length

implicit none

private

! Physical constants taken from NIST reference data
real(p_double), parameter :: const_e = 1.602176634d-19 ! C
real(p_double), parameter :: const_c = 2.997924580d+08 ! m/s
real(p_double), parameter :: const_m = 9.109383701d-31 ! kg

interface radspect_hist
    module procedure radspect_hist
end interface

interface report_radspect
    module procedure report_radspect
end interface

interface emission_chi_hist
    module procedure emission_chi_hist
end interface

interface report_chi_emit
    module procedure report_chi_emit
end interface

interface radanglespect_hist
    module procedure radanglespect_hist
end interface

interface report_radanglespect
    module procedure report_radanglespect
end interface

interface radsphere_hist
    module procedure radsphere_hist
end interface

interface report_radsphere
    module procedure report_radsphere
end interface

interface detector_hist
    module procedure detector_hist
end interface

interface report_detector
    module procedure report_detector
end interface

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

public :: radspect_hist, report_radspect, radanglespect_hist, &
            report_radanglespect, n_emit_hist, report_n_emit, &
            chi_hist, report_chi, emission_chi_hist, &
            report_chi_emit, radsphere_hist, report_radsphere, &
            detector_hist, report_detector

contains

subroutine radspect_hist(this, p_g, q_g)

    implicit none

    class(t_species_qed), intent(inout) :: this
    real(p_k_part), intent(in) :: p_g, q_g

    real(p_double) :: bin_size
    integer :: index

    bin_size = (log10(this%radspect_emax)-log10(this%radspect_emin))/size(this%radspect)

    if (p_g < this%radspect_emin) then
        this%radspect_oor(1) = this%radspect_oor(1) + 1
    else if (p_g > this%radspect_emax) then
        this%radspect_oor(2) = this%radspect_oor(2) + 1
    else
        index = ceiling((log10(p_g)-log10(this%radspect_emin))/bin_size)
        this%radspect(index) = this%radspect(index) + q_g
    end if

end subroutine radspect_hist

subroutine report_radspect( this, no_co, tstep, t )

    use m_diagnostic_utilities
    use stringutil, only: idx_string

    implicit none

    class(t_species_qed), intent(inout) :: this
    class( t_node_conf ), intent(in) :: no_co
    type( t_time_step ), intent(in) :: tstep
    real(p_double), intent(in) :: t

    real(p_double), dimension(this%radspect_bins) :: radspect_buffer, bin_sizes, &
                                                     x_axis
    real(p_double), dimension(this%radspect_bins + 1) :: x_axis_log, E
    class( t_diag_file ), allocatable :: diagFile
    real(p_double) :: dxvol, fac_weight_to_number, fac_mc2_to_mev, fac_mc2_to_joule, fac
    integer :: i, n_x_dim

    ! radiated energy spectrum diagnostic
    if (test_if_report( tstep, this%ndump_fac_radspect )) then

        ! prepare the logarithmically spaced axis
        call linspace(log10(this%radspect_emin), log10(this%radspect_emax), x_axis_log)
        E = 10**x_axis_log
        ! calculate the energy that each bin represents
        do i = 1, this%radspect_bins
            bin_sizes(i) = E(i+1)-E(i)
            x_axis(i) = ( E(i) + E(i+1) ) / 2
        enddo

        ! We convert the weight into a number
        dxvol = product( this%dx(1 : p_x_dim) )

        ! copy radiation spectrum to buffer array
        radspect_buffer = this%radspect * dxvol

        ! spectrum dE / d \gamma in plasma units
        radspect_buffer = radspect_buffer*x_axis/bin_sizes

        ! convert the weight into a real number of particles
        ! 1.e6 factor to convert n0 in /m3
        fac_weight_to_number = 1.0d6 * this%n0 * ( const_c / this%omega_p0 )**3

        ! convert from mc2 to Joules
        fac_mc2_to_joule = const_m * const_c**2

        ! convert from mc2 to MeV
        fac_mc2_to_mev = 1.0d-6 * const_m * const_c**2 / const_e

        ! total conversion factor from plasma units to J/MeV
        fac = fac_weight_to_number * fac_mc2_to_joule / fac_mc2_to_mev

        ! Apply the conversion factor to the buffer
        radspect_buffer = fac * radspect_buffer

        call reduce_array(no_co, radspect_buffer, operation = p_sum )

        if ( no_co%my_aid() == 1 ) then

            ! create diagFile object
            call create_diag_file( diagFile )

            ! File will be of type grid
            diagFile%ftype = p_diag_grid

            diagFile % filename = 'radspect-' // idx_string( n(tstep)/ndump(tstep), p_time_length )
            diagFile % filepath = trim(path_mass) // 'QED' // p_dir_sep // &
                                  trim("radspect") // p_dir_sep // &
                                  trim(this%name) // p_dir_sep

            diagFile%name = 'radspect'

            ! Iteration info
            diagFile%iter%n = n(tstep)
            diagFile%iter%t = t
            diagFile%iter%time_units = '1 / \omega_p'

            ! Get the dimension of the simulation
            n_x_dim = this % get_n_x_dims()

            ! Initialize grid file metadata
            diagFile%grid%ndims = 1
            diagFile%grid%name = 'radspect'
            diagFile%grid%label = 'dE_\gamma/dE'

            select case ( n_x_dim )
            case( 1 )
                diagFile%grid%units = 'J/MeV/m^2'
            case( 2 )
                diagFile%grid%units = 'J/MeV/m'
            case( 3 )
                diagFile%grid%units = 'J/MeV'
            case default
                diagFile%grid%units = 'J/MeV'
            end select

            diagFile%grid%count(1) = this%radspect_bins

            diagFile%grid%axis(1)%min = log10(this%radspect_emin*fac_mc2_to_mev)
            diagFile%grid%axis(1)%max = log10(this%radspect_emax*fac_mc2_to_mev)

            diagFile%grid%axis(1)%name = 'log_{10}(E)'
            diagFile%grid%axis(1)%label = 'log_{10}(E)'
            diagFile%grid%axis(1)%units = 'MeV'
            diagFile%grid%axis(1)%type = diag_axis_log10

            ! Open the file
            call diagFile % open( p_diag_create )
            call diagFile % add_dataset( 'radspect', radspect_buffer )

            ! Close the file
            call diagFile % close()

            ! Free diagFile object
            deallocate( diagFile )

        endif

    endif
end subroutine report_radspect

subroutine emission_chi_hist(this, eta, chi_g, q_g)

    implicit none

    class(t_species_qed), intent(inout) :: this
    real(p_k_part), intent(in) :: eta, chi_g, q_g

    real(p_k_part) :: bin_size
    integer :: index

    bin_size = (log10(this%chi_emit_max)-log10(this%chi_emit_min))/size(this%chispect_leptons)


    if ((eta < this%chi_emit_max) .and. (eta > this%chi_emit_min)) then
        index = ceiling((log10(eta)-log10(this%chi_emit_min))/bin_size)
        this%chispect_leptons(index) = this%chispect_leptons(index) + q_g
    endif

    if ((chi_g < this%chi_emit_max) .and. (chi_g > this%chi_emit_min)) then
        index = ceiling((log10(chi_g)-log10(this%chi_emit_min))/bin_size)
        this%chispect_photons(index) = this%chispect_photons(index) + q_g
    endif

end subroutine emission_chi_hist

subroutine report_chi_emit( this, no_co, tstep, t )

    use m_diagnostic_utilities
    use stringutil, only: idx_string

    implicit none

    class(t_species_qed), intent(inout) :: this
    class( t_node_conf ), intent(in) :: no_co
    type( t_time_step ), intent(in) :: tstep
    real(p_double), intent(in) :: t

    real(p_double), dimension(this%chi_emit_nbins) :: chi_emit_buffer_leptons, chi_emit_buffer_photons, bin_sizes, &
                                                     x_axis
    real(p_double), dimension(this%chi_emit_nbins + 1) :: x_axis_log, E
    class( t_diag_file ), allocatable :: diagFile, diagFile1
    integer :: i

    ! radiated energy spectrum diagnostic
    if (test_if_report( tstep, this%ndump_fac_chi_emit )) then

        ! prepare the logarithmically spaced axis
        call linspace(log10(this%chi_emit_min), log10(this%chi_emit_max), x_axis_log)
        E = 10**x_axis_log
        ! calculate the energy that each bin represents
        do i = 1, this%chi_emit_nbins
            bin_sizes(i) = E(i+1)-E(i)
            x_axis(i) = ( E(i) + E(i+1) ) / 2
        enddo

        ! copy radiation spectrum to buffer array
        chi_emit_buffer_leptons = this%chispect_leptons
        chi_emit_buffer_photons = this%chispect_photons


        ! this sis normalized to one
        chi_emit_buffer_leptons = chi_emit_buffer_leptons/(bin_sizes)
        chi_emit_buffer_photons = chi_emit_buffer_photons/(bin_sizes)

        call reduce_array(no_co, chi_emit_buffer_leptons, operation = p_sum )
        call reduce_array(no_co, chi_emit_buffer_photons, operation = p_sum )


        if ( no_co%my_aid() == 1 ) then

            ! create diagFile object
            call create_diag_file( diagFile )

            ! File will be of type grid
            diagFile%ftype = p_diag_grid

            diagFile % filename = 'chi-' //trim(this%name) // '-' // idx_string( n(tstep)/ndump(tstep), p_time_length )
            diagFile % filepath = trim(path_mass) // trim("QED") // p_dir_sep // trim("chi_emission") // p_dir_sep // &
                                  trim(this%name) // p_dir_sep

            ! Iteration info
            diagFile%iter%n = n(tstep)
            diagFile%iter%t = t
            diagFile%iter%time_units = '1 / \omega_p'

            ! Initialize grid file metadata
            diagFile%grid%ndims = 1
            diagFile%grid%name = 'chispect'
            diagFile%grid%label = 'dN_e/d\chi'
            diagFile%grid%units = 'arb. u.'

            diagFile%grid%count(1) = this%chi_emit_nbins

            diagFile%grid%axis(1)%min = log10(this%chi_emit_min)
            diagFile%grid%axis(1)%max = log10(this%chi_emit_max)

            diagFile%grid%axis(1)%name = 'log_{10}(\chi)'
            diagFile%grid%axis(1)%label = 'log_{10}(\chi)'
            diagFile%grid%axis(1)%units = ' '
          ! diagFile%grid%axis(1)%type = diag_axis_log10

            ! Open the file
            call diagFile % open( p_diag_create )
            call diagFile % add_dataset( 'chi_leptons', chi_emit_buffer_leptons )

            ! Close the file
            call diagFile % close()

            ! Free diagFile object
            deallocate( diagFile )


            ! create diagFile object
            call create_diag_file( diagFile1 )

            ! File will be of type grid
            diagFile1%ftype = p_diag_grid

            diagFile1 % filename = 'chi-phot-' //trim(this%name) // '-' // idx_string( n(tstep)/ndump(tstep), p_time_length )
            diagFile1 % filepath = trim(path_mass) // trim("QED") // p_dir_sep // trim("chi_emission") // p_dir_sep // &
                                   trim("photons-from-") // trim(this%name) // p_dir_sep

            ! Iteration info
            diagFile1%iter%n = n(tstep)
            diagFile1%iter%t = t
            diagFile1%iter%time_units = '1 / \omega_p'

            ! Initialize grid file metadata
            diagFile1%grid%ndims = 1
            diagFile1%grid%name = 'chispect'
            diagFile1%grid%label = 'dN_\gamma/d\chi'
            diagFile1%grid%units = 'arb. u.'

            diagFile1%grid%count(1) = this%chi_emit_nbins

            diagFile1%grid%axis(1)%min = log10(this%chi_emit_min)
            diagFile1%grid%axis(1)%max = log10(this%chi_emit_max)

            diagFile1%grid%axis(1)%name = 'log_{10}(\chi)'
            diagFile1%grid%axis(1)%label = 'log_{10}(\chi)'
            diagFile1%grid%axis(1)%units = ' '
          ! diagFile1%grid%axis(1)%type = diag_axis_log10

            ! Open the file
            call diagFile1 % open( p_diag_create )
            call diagFile1 % add_dataset( 'chi_photons', chi_emit_buffer_photons )

            ! Close the file
            call diagFile1 % close()

            ! Free diagFile object
            deallocate( diagFile1 )



        endif

    endif
end subroutine report_chi_emit

subroutine linspace(from, to, array)
    use m_parameters
    real(p_double), intent(in) :: from, to
    real(p_double), intent(inout) :: array(:)

    real(p_double) :: range
    integer :: n,i

    n = size(array)
    range = to - from
    if (n == 0) return
    if (n == 1) then
        array(1) = from
        return
    endif
    do i = 1,n
        array(i) = from + range * (i - 1) / (n - 1)
    end do
end subroutine linspace

subroutine radanglespect_hist(this, pph, q_g)
    use m_random

    implicit none

    class(t_species_qed), intent(inout) :: this
    real(p_k_part), dimension(p_p_dim), intent(in) :: pph
    real(p_k_part), intent(in) :: q_g

    real(p_double) :: min, max, ene, phi, theta, bin_size
    integer :: k, index_phi, index_theta

    min = -0.7853981634
    max = 0.7853981634
    bin_size = (max - min) / this%radanglespect_bins

    ene = sqrt( pph(1)**2 + pph(2)**3 + pph(3)**2 )

    if (ene >= this%radanglespect_thresh) then

        ! we look at the sphere from the outside, with x3 up
        if ( abs(pph(1)) >= abs(pph(2)) .and. abs(pph(1)) >= abs(pph(3)) ) then
            if ( pph(1) > 0) then
                k = 1 !we are on positive x1
                phi = atan( pph(2) / pph(1) )
                theta = atan( pph(3) / pph(1) )

                index_phi = ceiling( (phi - min) / bin_size )
                index_theta = ceiling( (theta - min) / bin_size )

                !if indices are within grid
                if ( (index_theta > 0) .and. (index_theta <= this%radanglespect_bins) .and. (index_phi > 0) &
                    .and. (index_phi <= this%radanglespect_bins)) then
                    !add to histogram
                   this%radanglespect( index_phi, index_theta, k) =this%radanglespect(index_phi, index_theta, k) &
                                                                + q_g*ene
                endif
            else
                k = 2 !we are on negative x1
                phi = atan( -pph(2) / (-pph(1)) )
                theta = atan( pph(3) / (-pph(1)) )

                index_phi = ceiling( (phi - min) / bin_size )
                index_theta = ceiling( (theta - min) / bin_size )

               !if indices are within grid
                if ( (index_theta > 0) .and. (index_theta <= this%radanglespect_bins) .and. (index_phi > 0) &
                    .and. (index_phi <= this%radanglespect_bins)) then
                    !add to histogram
                   this%radanglespect( index_phi, index_theta, k) =this%radanglespect(index_phi, index_theta, k) &
                                                                + q_g*ene
                endif
            endif
        endif

        if ( abs(pph(2)) >= abs(pph(1)) .and. abs(pph(2)) >= abs(pph(3)) ) then
            if ( pph(2) > 0 ) then
                k = 3 !we are on positive x2
                phi = ( -pph(1) / pph(2) )
                theta = ( pph(3) / pph(2) )

                index_phi = ceiling( (phi - min) / bin_size )
                index_theta = ceiling( (theta - min) / bin_size )

                !if indices are within grid
                if ( (index_theta > 0) .and. (index_theta <= this%radanglespect_bins) .and. (index_phi > 0) &
                    .and. (index_phi <= this%radanglespect_bins)) then
                    !add to histogram
                   this%radanglespect( index_phi, index_theta, k) =this%radanglespect(index_phi, index_theta, k) &
                                                                + q_g*ene
                endif
            else
                k = 4 !we are on negative x2
                phi = atan( pph(1) / (-pph(2)) )
                theta = atan( pph(3) / (-pph(2)) )

                index_phi = ceiling( (phi - min) / bin_size )
                index_theta = ceiling( (theta - min) / bin_size )

              !if indices are within grid
                if ( (index_theta > 0) .and. (index_theta <= this%radanglespect_bins) .and. (index_phi > 0) &
                    .and. (index_phi <= this%radanglespect_bins)) then
                    !add to histogram
                   this%radanglespect( index_phi, index_theta, k) =this%radanglespect(index_phi, index_theta, k) &
                                                                + q_g*ene
                endif
            endif
        endif

        if ( abs(pph(3)) >= abs(pph(2)) .and. abs(pph(3)) >= abs(pph(1)) ) then
            if ( pph(3) > 0 ) then
                k = 5 !we are on positive x3
                phi = atan( -pph(2) / pph(3) )
                theta = atan( pph(1) / pph(3))
                index_phi = ceiling( (phi - min) / bin_size )
                index_theta = ceiling( (theta - min) / bin_size )

               !if indices are within grid
                if ( (index_theta > 0) .and. (index_theta <= this%radanglespect_bins) .and. (index_phi > 0) &
                    .and. (index_phi <= this%radanglespect_bins)) then
                    !add to histogram
                   this%radanglespect( index_phi, index_theta, k) =this%radanglespect(index_phi, index_theta, k) &
                                                                + q_g*ene
                endif
            else
                k = 6 !we are on negative x3
                phi = atan( pph(2) / (-pph(3)) )
                theta = atan( pph(1) / (-pph(3)) )

                index_phi = ceiling( (phi - min) / bin_size )
                index_theta = ceiling( (theta - min) / bin_size )

                !if indices are within grid
                if ( (index_theta > 0) .and. (index_theta <= this%radanglespect_bins) .and. (index_phi > 0) &
                    .and. (index_phi <= this%radanglespect_bins)) then
                    !add to histogram
                   this%radanglespect( index_phi, index_theta, k) =this%radanglespect(index_phi, index_theta, k) &
                                                                + q_g*ene
                endif
            endif
        endif
    endif
end subroutine radanglespect_hist

subroutine report_radanglespect( this, no_co, tstep, t )

    use m_diagnostic_utilities
    use stringutil, only: idx_string

    implicit none

    class(t_species_qed), intent(inout) :: this
    class( t_node_conf ), intent(in) :: no_co
    type( t_time_step ), intent(in) :: tstep
    real(p_double), intent(in) :: t

    real(p_double), dimension(this%radanglespect_bins, this%radanglespect_bins) :: radanglespect_buffer1, &
                                                                                   &radanglespect_buffer2, &
                                                                                   &radanglespect_buffer3, &
                                                                                   &radanglespect_buffer4, &
                                                                                   &radanglespect_buffer5, &
                                                                                   &radanglespect_buffer6
    class( t_diag_file ), allocatable :: diagFile1,diagFile2, diagFile3, diagFile4, diagFile5, &
                                         diagFile6

    ! radiated energy spectrum diagnostic
    if (test_if_report( tstep, this%ndump_fac_radanglespect )) then

          ! copy radiation spectrum to buffer array
        radanglespect_buffer1 = this%radanglespect(:,:,1)
        radanglespect_buffer2 = this%radanglespect(:,:,2)
        radanglespect_buffer3 = this%radanglespect(:,:,3)
        radanglespect_buffer4 = this%radanglespect(:,:,4)
        radanglespect_buffer5 = this%radanglespect(:,:,5)
        radanglespect_buffer6 = this%radanglespect(:,:,6)

        call reduce_array(no_co, radanglespect_buffer1, operation = p_sum )
        call reduce_array(no_co, radanglespect_buffer2, operation = p_sum )
        call reduce_array(no_co, radanglespect_buffer3, operation = p_sum )
        call reduce_array(no_co, radanglespect_buffer4, operation = p_sum )
        call reduce_array(no_co, radanglespect_buffer5, operation = p_sum )
        call reduce_array(no_co, radanglespect_buffer6, operation = p_sum )

        if ( no_co%my_aid() == 1 ) then
            ! positive x1
            ! create diagFile1 object
            call create_diag_file( diagFile1 )

            ! File will be of type grid
            diagFile1%ftype = p_diag_grid

            diagFile1 % filename = 'radangle_1_1-num-' // idx_string( n(tstep)/ndump(tstep), p_time_length )
            diagFile1 % filepath = trim(path_mass) // 'QED' // p_dir_sep // &
                                   trim("radanglespect") // p_dir_sep // &
                                   trim(this%name) // p_dir_sep

            diagFile1%name = 'radanglespect'

            ! Iteration info
            diagFile1%iter%n = n(tstep)
            diagFile1%iter%t = t
            diagFile1%iter%time_units = '1 / \omega_p'

            ! Initialize grid file metadata
            diagFile1%grid%ndims = 2
            diagFile1%grid%name = 'radanglespect'
            diagFile1%grid%label = 'dN_\gamma/d\Omega'
            diagFile1%grid%units = 'arb. u.'

            diagFile1%grid%axis(1)%min = -0.7853981634 !pi/4
            diagFile1%grid%axis(1)%max = 0.7853981634
            diagFile1%grid%axis(2)%min = -0.7853981634
            diagFile1%grid%axis(2)%max = 0.7853981634

            diagFile1%grid%axis(1)%units = 'rad'
            diagFile1%grid%axis(2)%units = 'rad'

            diagFile1%grid%axis(1)%name = '\theta'
            diagFile1%grid%axis(1)%label = '\theta'

            diagFile1%grid%axis(2)%name = '\phi'
            diagFile1%grid%axis(2)%label = '\phi'

            ! Open the file
            call diagFile1 % open( p_diag_create )
            call diagFile1 % add_dataset( 'radanglespect', radanglespect_buffer1 )

            ! Close the file
            call diagFile1 % close()

            ! Free diagFile1 object
            deallocate( diagFile1 )


            ! negative x1
            ! create diagFile2 object
            call create_diag_file( diagFile2 )

            ! File will be of type grid
            diagFile2%ftype = p_diag_grid

            diagFile2 % filename = 'radangle_1_2-num-' // idx_string( n(tstep)/ndump(tstep), p_time_length )
            diagFile2 % filepath = trim(path_mass) // 'QED' // p_dir_sep // &
                                   trim("radanglespect") // p_dir_sep // &
                                   trim(this%name) // p_dir_sep

            diagFile2%name = 'radanglespect'

            ! Iteration info
            diagFile2%iter%n = n(tstep)
            diagFile2%iter%t = t
            diagFile2%iter%time_units = '1 / \omega_p'

            ! Initialize grid file metadata
            diagFile2%grid%ndims = 2
            diagFile2%grid%name = 'radanglespect'

            diagFile2%grid%count(1) = this%radanglespect_bins
            diagFile2%grid%count(2) = this%radanglespect_bins

            diagFile2%grid%axis(1)%min = -0.7853981634
            diagFile2%grid%axis(1)%max = 0.7853981634
            diagFile2%grid%axis(2)%min = -0.7853981634
            diagFile2%grid%axis(2)%max = 0.7853981634

            diagFile2%grid%axis(1)%units = 'rad'
            diagFile2%grid%axis(2)%units = 'rad'

            diagFile2%grid%axis(1)%name = '\theta'
            diagFile2%grid%axis(1)%label = '\theta'

            diagFile2%grid%axis(2)%name = '\phi'
            diagFile2%grid%axis(2)%label = '\phi'

            ! Open the file
            call diagFile2 % open( p_diag_create )
            call diagFile2 % add_dataset( 'radanglespect', radanglespect_buffer2 )

            ! Close the file
            call diagFile2 % close()

            ! Free diagFile2 object
            deallocate( diagFile2 )

            ! positive x2
            ! create diagFile3 object
            call create_diag_file( diagFile3 )

            ! File will be of type grid
            diagFile3%ftype = p_diag_grid

            diagFile3 % filename = 'radangle_2_1-num-' // idx_string( n(tstep)/ndump(tstep), p_time_length )
            diagFile3 % filepath = trim(path_mass) // 'QED' // p_dir_sep // &
                                   trim("radanglespect") // p_dir_sep // &
                                   trim(this%name) // p_dir_sep

            diagFile3%name = 'radanglespect'

            ! Iteration info
            diagFile3%iter%n = n(tstep)
            diagFile3%iter%t = t
            diagFile3%iter%time_units = '1 / \omega_p'

            ! Initialize grid file metadata
            diagFile3%grid%ndims = 2
            diagFile3%grid%name = 'radanglespect'

            diagFile3%grid%count(1) = this%radanglespect_bins
            diagFile3%grid%count(2) = this%radanglespect_bins

            diagFile3%grid%axis(1)%min = -0.7853981634
            diagFile3%grid%axis(1)%max = 0.7853981634
            diagFile3%grid%axis(2)%min = -0.7853981634
            diagFile3%grid%axis(2)%max = 0.7853981634

            diagFile3%grid%axis(1)%units = 'rad'
            diagFile3%grid%axis(2)%units = 'rad'

            diagFile3%grid%axis(1)%name = '\theta'
            diagFile3%grid%axis(1)%label = '\theta'

            diagFile3%grid%axis(2)%name = '\phi'
            diagFile3%grid%axis(2)%label = '\phi'

            ! Open the file
            call diagFile3 % open( p_diag_create )
            call diagFile3 % add_dataset( 'radanglespect', radanglespect_buffer3 )

            ! Close the file
            call diagFile3 % close()

            ! Free diagFile3 object
            deallocate( diagFile3 )


            ! negative x2
            ! create diagFile4 object
            call create_diag_file( diagFile4 )

            ! File will be of type grid
            diagFile4%ftype = p_diag_grid

            diagFile4 % filename = 'radangle_2_2-num-' // idx_string( n(tstep)/ndump(tstep), p_time_length )
            diagFile4 % filepath = trim(path_mass) // 'QED' // p_dir_sep // &
                                   trim("radanglespect") // p_dir_sep // &
                                   trim(this%name) // p_dir_sep

            diagFile4%name = 'radanglespect'

            ! Iteration info
            diagFile4%iter%n = n(tstep)
            diagFile4%iter%t = t
            diagFile4%iter%time_units = '1 / \omega_p'

            ! Initialize grid file metadata
            diagFile4%grid%ndims = 2
            diagFile4%grid%name = 'radanglespect'

            diagFile4%grid%count(1) = this%radanglespect_bins
            diagFile4%grid%count(2) = this%radanglespect_bins

            diagFile4%grid%axis(1)%min = -0.7853981634
            diagFile4%grid%axis(1)%max = 0.7853981634
            diagFile4%grid%axis(2)%min = -0.7853981634
            diagFile4%grid%axis(2)%max = 0.7853981634

            diagFile4%grid%axis(1)%units = 'rad'
            diagFile4%grid%axis(2)%units = 'rad'

            diagFile4%grid%axis(1)%name = '\theta'
            diagFile4%grid%axis(1)%label = '\theta'

            diagFile4%grid%axis(2)%name = '\phi'
            diagFile4%grid%axis(2)%label = '\phi'

            ! Open the file
            call diagFile4 % open( p_diag_create )
            call diagFile4 % add_dataset( 'radanglespect', radanglespect_buffer4 )

            ! Close the file
            call diagFile4 % close()

            ! Free diagFile4 object
            deallocate( diagFile4 )


            ! positive x3
            ! create diagFile5 object
            call create_diag_file( diagFile5 )

            ! File will be of type grid
            diagFile5%ftype = p_diag_grid

            diagFile5 % filename = 'radangle_3_1-num-' // idx_string( n(tstep)/ndump(tstep), p_time_length )
            diagFile5 % filepath = trim(path_mass) // 'QED' // p_dir_sep // &
                                   trim("radanglespect") // p_dir_sep // &
                                   trim(this%name) // p_dir_sep

            diagFile5%name = 'radanglespect'

            ! Iteration info
            diagFile5%iter%n = n(tstep)
            diagFile5%iter%t = t
            diagFile5%iter%time_units = '1 / \omega_p'

            ! Initialize grid file metadata
            diagFile5%grid%ndims = 2
            diagFile5%grid%name = 'radanglespect'

            diagFile5%grid%count(1) = this%radanglespect_bins
            diagFile5%grid%count(2) = this%radanglespect_bins

            diagFile5%grid%axis(1)%min = -0.7853981634
            diagFile5%grid%axis(1)%max = 0.7853981634
            diagFile5%grid%axis(2)%min = -0.7853981634
            diagFile5%grid%axis(2)%max = 0.7853981634

            diagFile5%grid%axis(1)%units = 'rad'
            diagFile5%grid%axis(2)%units = 'rad'

            diagFile5%grid%axis(1)%name = '\theta'
            diagFile5%grid%axis(1)%label = '\theta'

            diagFile5%grid%axis(2)%name = '\phi'
            diagFile5%grid%axis(2)%label = '\phi'

            ! Open the file
            call diagFile5 % open( p_diag_create )
            call diagFile5 % add_dataset( 'radanglespect', radanglespect_buffer5 )

            ! Close the file
            call diagFile5 % close()

            ! Free diagFile5 object
            deallocate( diagFile5 )


            ! negative x3
            ! create diagFile6 object
            call create_diag_file( diagFile6 )

            ! File will be of type grid
            diagFile6%ftype = p_diag_grid

            diagFile6 % filename = 'radangle_3_2-num-' // idx_string( n(tstep)/ndump(tstep), p_time_length )
            diagFile6 % filepath = trim(path_mass) // 'QED' // p_dir_sep // &
                                   trim("radanglespect") // p_dir_sep // &
                                   trim(this%name) // p_dir_sep

            diagFile6%name = 'radanglespect'

            ! Iteration info
            diagFile6%iter%n = n(tstep)
            diagFile6%iter%t = t
            diagFile6%iter%time_units = '1 / \omega_p'

            ! Initialize grid file metadata
            diagFile6%grid%ndims = 2
            diagFile6%grid%name = 'radanglespect'

            diagFile6%grid%count(1) = this%radanglespect_bins
            diagFile6%grid%count(2) = this%radanglespect_bins

            diagFile6%grid%axis(1)%min = -0.7853981634
            diagFile6%grid%axis(1)%max = 0.7853981634
            diagFile6%grid%axis(2)%min = -0.7853981634
            diagFile6%grid%axis(2)%max = 0.7853981634

            diagFile6%grid%axis(1)%units = 'rad'
            diagFile6%grid%axis(2)%units = 'rad'

            diagFile6%grid%axis(1)%name = '\theta'
            diagFile6%grid%axis(1)%label = '\theta'

            diagFile6%grid%axis(2)%name = '\phi'
            diagFile6%grid%axis(2)%label = '\phi'

            ! Open the file
            call diagFile6 % open( p_diag_create )
            call diagFile6 % add_dataset( 'radanglespect', radanglespect_buffer6 )

            ! Close the file
            call diagFile6 % close()

            ! Free diagFile6 object
            deallocate( diagFile6 )
        endif
    endif
end subroutine report_radanglespect

subroutine radsphere_hist(this, pph, q_g)
    use m_random

    implicit none

    class(t_species_qed), intent(inout) :: this
    real(p_k_part), dimension(p_p_dim), intent(in) :: pph
    real(p_k_part), intent(in) :: q_g

    real(p_double) :: min, max, ene, phi, theta, bin_size
    integer :: k, index_phi, index_theta

    min = -0.7853981634
    max = 0.7853981634
    bin_size = (max - min) / this%radsphere_nbins

    ene = sqrt( pph(1)**2 + pph(2)**2 + pph(3)**2 )

    if ( (ene >= this%radsphere_emin) .and. (ene <= this%radsphere_emax) ) then

        ! we look at the sphere from the outside, with x3 up
        if ( (abs(pph(1)) >= abs(pph(2))) .and. (abs(pph(1)) >= abs(pph(3))) ) then
            if ( pph(1) > 0.0) then
                k = 1 !we are on positive x1
                phi = atan( pph(2) / pph(1) )
                theta = atan( pph(3) / pph(1) )
            else
                k = 2 !we are on negative x1
                phi = atan( -pph(2) / abs(pph(1)) )
                theta = atan( pph(3) / abs(pph(1)) )
            endif
        else
            if ( abs(pph(2)) >= abs(pph(3)) ) then
                if ( pph(2) > 0.0 ) then
                    k = 3 !we are on positive x2
                    phi = atan( -pph(1) / pph(2) )
                    theta = atan( pph(3) / pph(2) )
                else
                    k = 4 !we are on negative x2
                    phi = atan( pph(1) / abs(pph(2)) )
                    theta = atan( pph(3) / abs(pph(2)) )
                endif
            else
                if ( pph(3) > 0.0 ) then
                        k = 5 !we are on positive x3
                        phi = atan( -pph(2) / pph(3) )
                        theta = atan( pph(1) / pph(3))
                    else
                        k = 6 !we are on negative x3
                        phi = atan( pph(2) / abs(pph(3)) )
                        theta = atan( pph(1) / abs(pph(3)) )
                    endif
            endif
        endif

        index_phi = ceiling( (phi - min) / bin_size )
        index_theta = ceiling( (theta - min) / bin_size )

        !if indices are within grid
        if ( ( (index_theta > 0) .and. (index_theta <= this%radsphere_nbins)) .and. ( (index_phi > 0) &
                .and. (index_phi <= this%radsphere_nbins))) then
                !add to histogram
            this%radsphere( index_phi, index_theta, k) =this%radsphere(index_phi, index_theta, k) &
                                                            + q_g*ene
        endif


    endif
end subroutine radsphere_hist

subroutine report_radsphere( this, no_co, tstep, t )

    use m_diagnostic_utilities
    use stringutil, only: idx_string

    implicit none

    class(t_species_qed), intent(inout) :: this
    class( t_node_conf ), intent(in) :: no_co
    type( t_time_step ), intent(in) :: tstep
    real(p_double), intent(in) :: t

    integer :: k, i, j
    character (len=10) :: detector_direction
    character, dimension(6) :: prefix
    integer, dimension(6) :: phi_dir, theta_dir
    real(p_double), dimension(this%radsphere_nbins, this%radsphere_nbins,6) :: radsphere_buff
    real(p_double) :: min, max, phi, theta, bin_size, factor

    class( t_diag_file ), allocatable :: diagFile1


    ! radiated energy spectrum diagnostic
    if (test_if_report( tstep, this%ndump_fac_radsphere )) then

          ! copy radiation spectrum to buffer array
        radsphere_buff = this%radsphere(:,:,:)

        call reduce_array(no_co, radsphere_buff, operation = p_sum )

        min = -0.7853981634
        max = 0.7853981634
        bin_size = (max - min) / this%radsphere_nbins

        ! normalise cells by their size
        do i=1, this%radsphere_nbins
            do j=1, this%radsphere_nbins
                phi = min + bin_size*i
                theta = min + bin_size*j
                ! normalize everything to the pixel size of d\phi d\omega at the main line of sight
                factor = sqrt(1-sin(theta)*sin(theta)*sin(phi)*sin(phi))
                radsphere_buff(i,j,:) = radsphere_buff(i,j,:)*(factor**3)/(cos(theta)*cos(phi))
            enddo
        enddo

        if ( no_co%my_aid() == 1 ) then

           prefix = (/' ', '-', '-', ' ', '-', ' '/)
           phi_dir(1:6) = (/2, 2, 1, 1, 2, 2/)
           theta_dir(1:6) = (/3, 3, 3, 3, 1, 1/)

            do k=1,6

            ! positive x1, negative x1, positive x2, negative x2, positive x3, negative x3
            ! create diagFile1 object
                call create_diag_file( diagFile1 )

                ! File will be of type grid
                diagFile1%ftype = p_diag_grid

                diagFile1 % filename = 'radsphere-'// idx_string( (k+1)/2, 1 ) //'-'// idx_string( mod(k+1,2)+1, 1 )//'-'// &
                                    trim(this%name)//'-' // idx_string( n(tstep)/ndump(tstep), p_time_length )
                diagFile1 % filepath = trim(path_mass) // trim("QED") // p_dir_sep // trim("radsphere") // p_dir_sep // &
                                    trim(this%name) // p_dir_sep

                ! Iteration info
                diagFile1%iter%n = n(tstep)
                diagFile1%iter%t = t
                diagFile1%iter%time_units = '1 / \omega_p'

                ! Initialize grid file metadata

                if ( mod( k, 2 )==1 ) then
                        detector_direction='positive '
                    else
                        detector_direction='negative '
                endif

                diagFile1%grid%ndims = 2
                diagFile1%grid%name = 'radsphere, '//detector_direction//'x_'//idx_string( (k+1)/2, 1 )//' direction'
                diagFile1%grid%label = 'dE_\gamma/d\Omega  along '//detector_direction//'x_'//idx_string( (k+1)/2, 1 )//' direction'
                diagFile1%grid%units = 'arb. u.'

                diagFile1%grid%axis(1)%min = -0.7853981634 !pi/4
                diagFile1%grid%axis(1)%max = 0.7853981634
                diagFile1%grid%axis(2)%min = -0.7853981634
                diagFile1%grid%axis(2)%max = 0.7853981634

                diagFile1%grid%axis(1)%units = 'rad'
                diagFile1%grid%axis(2)%units = 'rad'

                diagFile1%grid%axis(1)%name = '\phi'
                diagFile1%grid%axis(1)%label = prefix(k)//'\theta_'//idx_string(phi_dir(k),1)

                diagFile1%grid%axis(2)%name = '\theta'
                diagFile1%grid%axis(2)%label = '\theta_'//idx_string(theta_dir(k),1)

                ! Open the file
                call diagFile1 % open( p_diag_create )
                call diagFile1 % add_dataset( 'radsphere_'// trim(detector_direction) // '_x' //idx_string( (k+1)/2, 1 ), radsphere_buff(:,:,k) )

                ! Close the file
                call diagFile1 % close()

                ! Free diagFile1 object
                deallocate( diagFile1 )

            enddo



        endif
    endif
end subroutine report_radsphere

subroutine detector_hist(this, pph, q_g)
    use m_random

    implicit none

    class(t_species_qed), intent(inout) :: this
    real(p_k_part), dimension(p_p_dim), intent(in) :: pph
    real(p_k_part), intent(in) :: q_g

    real(p_double) :: min_phi, max_phi, min_theta, max_theta, ene, phi, theta, bin_size_phi, bin_size_theta
    integer :: k, index_phi, index_theta

    ! set the angular range - this should not be larger than (-pi/4, pi/4), but can be smaller and set individually for each direction
    min_phi = this%raddetector_phimin
    max_phi = this%raddetector_phimax
    min_theta = this%raddetector_thetamin
    max_theta = this%raddetector_thetamax
    bin_size_phi = (max_phi - min_phi) / this%raddetector_nbins
    bin_size_theta = (max_theta - min_theta) / this%raddetector_nbins


    ene = sqrt( pph(1)**2 + pph(2)**2 + pph(3)**2 )

    if ( (ene >= this%raddetector_emin) .and. (ene <= this%raddetector_emax) ) then

        ! we look at the sphere from the outside, with x3 up
        if ( (abs(pph(1)) >= abs(pph(2))) .and. (abs(pph(1)) >= abs(pph(3))) ) then
            if ( pph(1) > 0.0) then
                k = 1 !we are on positive x1
                phi = atan( pph(2) / pph(1) )
                theta = atan( pph(3) / pph(1) )
            else
                k = 2 !we are on negative x1
                phi = atan( -pph(2) / abs(pph(1)) )
                theta = atan( pph(3) / abs(pph(1)) )
            endif
        else
            if ( abs(pph(2)) >= abs(pph(3)) ) then
                if ( pph(2) > 0.0 ) then
                    k = 3 !we are on positive x2
                    phi = atan( -pph(1) / pph(2) )
                    theta = atan( pph(3) / pph(2) )
                else
                    k = 4 !we are on negative x2
                    phi = atan( pph(1) / abs(pph(2)) )
                    theta = atan( pph(3) / abs(pph(2)) )
                endif
            else
                if ( pph(3) > 0.0 ) then
                        k = 5 !we are on positive x3
                        phi = atan( -pph(2) / pph(3) )
                        theta = atan( pph(1) / pph(3))
                    else
                        k = 6 !we are on negative x3
                        phi = atan( pph(2) / abs(pph(3)) )
                        theta = atan( pph(1) / abs(pph(3)) )
                    endif
            endif
        endif

        index_phi = ceiling( (phi - min_phi) / bin_size_phi )
        index_theta = ceiling( (theta - min_theta) / bin_size_theta )

        if (this%raddetector_dir == k) then
        !if indices are within grid
            if ( ( (index_theta > 0) .and. (index_theta <= this%raddetector_nbins)) .and. ( (index_phi > 0) &
                    .and. (index_phi <= this%raddetector_nbins))) then
                    !add to histogram
                this%raddetector( index_phi, index_theta) =this%raddetector(index_phi, index_theta) &
                                                                + q_g*ene
            endif
        endif


    endif
end subroutine detector_hist

subroutine report_detector( this, no_co, tstep, t )

    use m_diagnostic_utilities
    use stringutil, only: idx_string

    implicit none

    class(t_species_qed), intent(inout) :: this
    class( t_node_conf ), intent(in) :: no_co
    type( t_time_step ), intent(in) :: tstep
    real(p_double), intent(in) :: t

    integer :: k, i, j
    character (len=10) :: detector_direction
    character, dimension(6) :: prefix
    integer, dimension(6) :: phi_dir, theta_dir
    real(p_double), dimension(this%raddetector_nbins, this%raddetector_nbins) :: raddetector_buff
    real(p_double) :: min, max, phi, theta, bin_size_phi, bin_size_theta, factor

    class( t_diag_file ), allocatable :: diagFile1


    ! radiated energy spectrum diagnostic
    if (test_if_report( tstep, this%ndump_fac_raddetector )) then

          ! copy radiation spectrum to buffer array
        raddetector_buff = this%raddetector(:,:)

        call reduce_array(no_co, raddetector_buff, operation = p_sum )

        bin_size_phi = (this%raddetector_phimax - this%raddetector_phimin) / this%raddetector_nbins
        bin_size_theta = (this%raddetector_thetamax - this%raddetector_thetamin) / this%raddetector_nbins
        ! normalise cells by their size
        do i=1, this%raddetector_nbins
            do j=1, this%raddetector_nbins
                phi = this%raddetector_phimin + bin_size_phi*i
                theta = this%raddetector_thetamin + bin_size_theta*j
                ! normalize everything to the pixel size of d\phi d\omega at the main line of sight
                factor = sqrt(1-sin(theta)*sin(theta)*sin(phi)*sin(phi))
                raddetector_buff(i,j) = raddetector_buff(i,j)*(factor**3)/(cos(theta)*cos(phi))
            enddo
        enddo

        if ( no_co%my_aid() == 1 ) then

           prefix = (/' ', '-', '-', ' ', '-', ' '/)
           phi_dir(1:6) = (/2, 2, 1, 1, 2, 2/)
           theta_dir(1:6) = (/3, 3, 3, 3, 1, 1/)

            k = this%raddetector_dir

            ! positive x1, negative x1, positive x2, negative x2, positive x3, negative x3
            ! create diagFile1 object
                call create_diag_file( diagFile1 )

                ! File will be of type grid
                diagFile1%ftype = p_diag_grid

                diagFile1 % filename = 'rad-detector-'// idx_string( (k+1)/2, 1 ) //'-'// idx_string( mod(k+1,2)+1, 1 )//'-'// &
                                    trim(this%name)//'-' // idx_string( n(tstep)/ndump(tstep), p_time_length )
                diagFile1 % filepath = trim(path_mass) // trim("QED") // p_dir_sep // trim("rad-detector") // p_dir_sep // &
                                    trim(this%name) // p_dir_sep

                ! Iteration info
                diagFile1%iter%n = n(tstep)
                diagFile1%iter%t = t
                diagFile1%iter%time_units = '1 / \omega_p'

                ! Initialize grid file metadata

                if ( mod( k, 2 )==1 ) then
                        detector_direction='positive '
                    else
                        detector_direction='negative '
                endif

                diagFile1%grid%ndims = 2
                diagFile1%grid%name = 'radiation, '//detector_direction//'x_'//idx_string( (k+1)/2, 1 )//' direction'
                diagFile1%grid%label = 'dE_\gamma/d\Omega  along '//detector_direction//'x_'//idx_string( (k+1)/2, 1 )//' direction'
                diagFile1%grid%units = 'arb. u.'

                diagFile1%grid%axis(1)%min = this%raddetector_phimin !pi/4
                diagFile1%grid%axis(1)%max = this%raddetector_phimax
                diagFile1%grid%axis(2)%min = this%raddetector_thetamin
                diagFile1%grid%axis(2)%max = this%raddetector_thetamax

                diagFile1%grid%axis(1)%units = 'rad'
                diagFile1%grid%axis(2)%units = 'rad'

                diagFile1%grid%axis(1)%name = '\phi'
                diagFile1%grid%axis(1)%label = prefix(k)//'\theta_'//idx_string(phi_dir(k),1)

                diagFile1%grid%axis(2)%name = '\theta'
                diagFile1%grid%axis(2)%label = '\theta_'//idx_string(theta_dir(k),1)

                ! Open the file
                call diagFile1 % open( p_diag_create )
                call diagFile1 % add_dataset( 'radiation_'// trim(detector_direction) // '_x' //idx_string( (k+1)/2, 1 ), raddetector_buff(:,:) )

                ! Close the file
                call diagFile1 % close()

                ! Free diagFile1 object
                deallocate( diagFile1 )





        endif
    endif
end subroutine report_detector

subroutine n_emit_hist(this, ix)

    implicit none

    class(t_species_qed), intent(inout) :: this
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

    class(t_species_qed), intent(inout) :: this
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

subroutine chi_hist(this, ix, q, eta)

    implicit none

    class(t_species_qed), intent(inout) :: this
    integer, dimension(p_x_dim), intent(in) :: ix
    real( p_k_part ), intent(in) :: q, eta
    real( p_k_part ) :: weight, weight_chi

    select case ( p_x_dim )

    case ( 1 )

        weight = this%weight%f1(1, ix(1))
        this%weight%f1(1, ix(1)) = weight + abs( q )
        weight_chi = this%weight_chi%f1(1, ix(1))
        this%weight_chi%f1(1, ix(1)) = weight_chi + abs( q ) * eta

    case ( 2 )

        weight = this%weight%f2(1, ix(1), ix(2))
        this%weight%f2(1, ix(1), ix(2)) = weight + abs( q )
        weight_chi = this%weight_chi%f2(1, ix(1), ix(2))
        this%weight_chi%f2(1, ix(1), ix(2)) = weight_chi + abs( q ) * eta

    case ( 3 )

        weight = this%weight%f3(1, ix(1), ix(2), ix(3))
        this%weight%f3(1, ix(1), ix(2), ix(3)) = weight + abs( q )
        weight_chi = this%weight_chi%f3(1, ix(1), ix(2), ix(3))
        this%weight_chi%f3(1, ix(1), ix(2), ix(3)) = weight_chi + abs( q ) * eta

    end select

end subroutine chi_hist

subroutine report_chi( this, g_space, grid, no_co, tstep, t )

    use m_diagnostic_utilities
    use stringutil, only: idx_string
    use m_vdf_math

    implicit none

    class(t_species_qed), intent(inout) :: this
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

end module m_species_diagnostics_qed

!------------------------------------------------------------------------------------------
subroutine report_spec_qed( this, emf, g_space, grid, no_co, tstep, t, tmin, send_msg, recv_msg )
!-----------------------------------------------------------------------------------------

    use m_system
    use m_parameters
    use m_species_define_qed, only: t_species_qed
    use m_emf_define, only : t_emf
    use m_space, only : t_space
    use m_grid_define, only : t_grid
    use m_node_conf
    use m_time_step
    use m_vdf_comm, only : t_vdf_msg
    use m_species_diagnostics
    use m_species_diagnostics_qed

    implicit none

    ! dummy variables

    class( t_species_qed ), intent(inout) :: this
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
    call report_radspect( this, no_co, tstep, t)
    call report_radanglespect( this, no_co, tstep, t)
    call report_n_emit( this, g_space, grid, no_co, tstep, t)
    call report_radsphere( this, no_co, tstep, t)
    call report_detector( this, no_co, tstep, t)
    call report_chi( this, g_space, grid, no_co, tstep, t)
    call report_chi_emit( this, no_co, tstep, t)

    ! Interpolation of ion density for bremsstrahlung
    ! It is done here because all the variables needed are available
    ! deposit_quant includes the deposit + the update_boundary.
    if (this%bremsstrahlung_info%if_bremsstrahlung ) then
        call this % bremsstrahlung_info % rho_ion % cleanup()
        call deposit_quant( this%bremsstrahlung_info%ion, grid, no_co, this % bremsstrahlung_info % rho_ion, p_charge, send_msg, recv_msg )
    end if

    ! Interpolation of ion density for trident_coul
    ! It is done here because all the variables needed are available
    ! deposit_quant includes the deposit + the update_boundary.
    if (this%trident_coul_info%if_trident_coul ) then
        call this % trident_coul_info % rho_ion % cleanup()
        call deposit_quant( this%trident_coul_info%ion, grid, no_co, this % trident_coul_info % rho_ion, p_charge, send_msg, recv_msg )
    end if

end subroutine report_spec_qed
!-----------------------------------------------------------------------------------------

!------------------------------------------------------------------------------------------
subroutine report_energy_rad_spec_qed( this, no_co, tstep, t, tmin )
!-----------------------------------------------------------------------------------------

    use m_system
    use m_species_define_qed, only: t_species_qed
    use m_node_conf
    use m_units
    use m_time_step
    use m_parameters, only: file_id_parene, p_err_diagfile, path_hist
    use stringutil, only: idx_string

    implicit none

    ! dummy variables
    class( t_species_qed ), intent(in) :: this
    class( t_node_conf ), intent(in) :: no_co
    type( t_time_step ), intent(in) :: tstep
    real(p_double), intent(in) :: t, tmin

    character(80) :: path, full_name
    character(2) :: ch_sp_id
    integer :: ierr
    real(p_double), dimension(5) :: energy_rad
    real(p_double) :: dxvol

    ! radiated particle energy diagnostic
    if (test_if_report( tstep, this%ndump_fac_rad )) then

        dxvol = product( this%dx(1 : p_x_dim) )

        ! get normalized radiated energy
        energy_rad = abs(this%energy_rad*this%rqm*dxvol)

        call reduce_array(no_co, energy_rad, operation = p_sum )

        if ( no_co%my_aid() == 1 ) then

            ! prepare path and file names
            ch_sp_id = idx_string( this%sp_id, 2 )

            path = trim(path_hist)

            full_name = trim(path) // 'par' // trim(ch_sp_id) // '_rad'

            ! Open file and position at the last record
            if ( t == tmin ) then
            call mkdir( path, ierr )

            open (unit=file_id_parene, file=full_name, status = 'REPLACE' , &
            form='formatted')

            ! Write Header
            write(file_id_parene,'( A6, (1X, A15), 5(1X, A28) )') &
                    'Iter','Time    ','Rad. Ene. (Class.)', 'Rad. Ene. (QED, non-tracked)', 'Rad. Ene. (QED, tracked)', 'Rad. Ene. (Brem, nn-tracked)', 'Rad. Ene. (Brem, tracked)'

            else
            open (unit=file_id_parene, file=full_name, position = 'append', &
            form='formatted')
            endif

            ! Writes the particle energies to file
            write(file_id_parene, '( I6,1X,g15.8,5(1X,es28.16) )') &
                    n(tstep), t, energy_rad(1), energy_rad(2), energy_rad(3), energy_rad(4), energy_rad(5)
            close(file_id_parene)

        endif
    endif

end subroutine report_energy_rad_spec_qed

!-----------------------------------------------------------------------------------------
! Report on number of pairs created
!-----------------------------------------------------------------------------------------
subroutine report_pairs_ct_spec_qed( this, no_co, tstep, t, tmin )
!-----------------------------------------------------------------------------------------

    use m_system
    use m_species_define_qed, only: t_species_qed
    use m_node_conf
    use m_units
    use m_time_step
    use m_parameters, only: file_id_parene, p_err_diagfile, path_hist
    use stringutil, only: idx_string

    implicit none

    ! dummy variables
    class( t_species_qed ), intent(in) :: this
    class( t_node_conf ), intent(in) :: no_co
    type( t_time_step ), intent(in) :: tstep
    real(p_double), intent(in) :: t, tmin

    ! local variables
    character(80) :: path, full_name
    character(2) :: ch_sp_id

    integer :: ierr
    integer( p_int64 ) :: total_pairs_ct
    real( p_double ):: weight_pairs_ct
    real( p_double ) :: energy_pairs_ct

    ! Pair creation diagnostics
    if (test_if_report( tstep, this%ndump_fac_pairs_from_elec )) then

        total_pairs_ct = this%num_pairs_ct
        call reduce( no_co, total_pairs_ct, operation = p_sum)

        weight_pairs_ct = this%weight_pairs_ct
        call reduce( no_co, weight_pairs_ct, operation = p_sum)

        energy_pairs_ct = this % get_energy_pairs_ct()
        call reduce( no_co, energy_pairs_ct, operation = p_sum)

        if ( no_co%my_aid() == 1 ) then

        ! prepare path and file names
        ch_sp_id = idx_string( this%sp_id, 2 )

        path = trim(path_hist)

        full_name = trim(path) // 'par' // trim(ch_sp_id) // '_pairs'

        ! Open file and position at the last record
        if ( t == tmin ) then
            call mkdir( path, ierr )

            open (unit=file_id_parene, file=full_name, status = 'REPLACE' , &
            form='formatted')

            ! Write Header
            write(file_id_parene,'( A6, (1X, A15), 3(3X, A20) )') &
                'Iter','Time    ','Total Pairs (CT)', 'Total Weight (CT)', 'Total Energy (CT)'

        else
            open (unit=file_id_parene, file=full_name, position = 'append', &
            form='formatted')
        endif

        ! Writes the particle energies to file
        write(file_id_parene, '( I6,1X,g15.8,1X,I14,es28.16,es28.16,I14,es28.16,es28.16 )') n(tstep), t, total_pairs_ct, &
                        weight_pairs_ct, energy_pairs_ct
        close(file_id_parene)

        endif

    endif

end subroutine report_pairs_ct_spec_qed


!-----------------------------------------------------------------------------------------
