#include "os-config.h"
#include "os-preprocess.fpp"

module m_photons

#include "memory/memory.h"

use m_system
use m_parameters
use m_math

use m_species_define, only: p_cell_near
use m_photons_define, only: t_photons
use m_species_define, only : p_ene_recalc
use m_emf_define, only: t_emf
use m_qed
use m_vdf_define
use m_qed_coulomb
use m_species_charge, only: deposit_rho

implicit none

private

public advance_photons_1d, advance_photons_2d, advance_photons_3d, advance_photons_cyl_2d
public create_pairs, ntrim_qed, create_pairs_betheheitler

contains

!---------------------------------------------------------------------------------------------------
subroutine advance_photons_1d( this, gdt, i0, i1)
!---------------------------------------------------------------------------------------------------

    implicit none

    integer, parameter :: rank = 2

    ! dummy variables
    class( t_photons ),    intent(inout) :: this
    real(p_double),   intent(in) :: gdt
    integer,   intent(in) :: i0, i1

    real(p_k_part) :: rdx
    integer :: i, pp, np, ptrcur

    real(p_k_part), dimension(rank,p_cache_size) :: xbuf
    integer, dimension(rank,p_cache_size)        :: dxi
    real(p_k_part), dimension(p_cache_size)      :: rp, rpp
    real(p_k_part) :: dt

    ! executable statements
    rdx = real( 1.0_p_double/this%dx(1), p_k_part )
    dt = real( gdt, p_k_part )

    ! advance position of photons i0 to i1 in chunks of p_cache_size
    do ptrcur = i0, i1, p_cache_size

        ! check if last copy of table and set np
        if( ptrcur + p_cache_size > i1 ) then
            np = i1 - ptrcur + 1
        else
            np = p_cache_size
        endif

        ! this type of loop is actually faster than
        ! using a forall construct
        pp = ptrcur

        do i=1,np
            rp(i) = 1.0_p_k_part / &
            sqrt( this%p(1,pp)**2 + this%p(2,pp)**2 + this%p(3,pp)**2 )
            rpp(i) = dt * rp(i)
        pp = pp + 1
        end do

        ! reference implementation
        pp = ptrcur
        do i=1,np
            xbuf(1,i) = this%x(1,pp) + this%p(1,pp) * rpp(i) * rdx
            dxi(1,i) = ntrim_qed( xbuf(1,i) )

            pp = pp + 1
        end do


        ! copy data from buffer to species data trimming positions
        pp = ptrcur
        do i = 1, np
            this%x(1,pp)  = xbuf(1,i) - dxi(1,i)
            this%ix(1,pp) = this%ix(1,pp) + dxi(1,i)

            pp = pp + 1
        end do

    enddo

end subroutine advance_photons_1d
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine advance_photons_2d( this, gdt, i0, i1)
!---------------------------------------------------------------------------------------------------

    implicit none

    integer, parameter :: rank = 2

    ! dummy variables
    class( t_photons ),    intent(inout) :: this
    real(p_double),   intent(in) :: gdt
    integer,   intent(in) :: i0, i1

    real(p_k_part), dimension(p_x_dim) :: rdx
    integer :: i, pp, np, ptrcur

    real(p_k_part), dimension(rank,p_cache_size) :: xbuf
    integer, dimension(rank,p_cache_size)        :: dxi
    real(p_k_part), dimension(p_cache_size)      :: rp, rpp
    real(p_k_part) :: dt

    ! executable statements
    rdx(1) = real( 1.0_p_double/this%dx(1), p_k_part )
    rdx(2) = real( 1.0_p_double/this%dx(2), p_k_part )
    dt = real( gdt, p_k_part )

    ! advance position of photons i0 to i1 in chunks of p_cache_size
    do ptrcur = i0, i1, p_cache_size

        ! check if last copy of table and set np
        if( ptrcur + p_cache_size > i1 ) then
            np = i1 - ptrcur + 1
        else
            np = p_cache_size
        endif

        ! this type of loop is actually faster than
        ! using a forall construct
        pp = ptrcur

        do i=1,np
            rp(i) = 1.0_p_k_part / &
            sqrt( this%p(1,pp)**2 + this%p(2,pp)**2 + this%p(3,pp)**2 )
            rpp(i) = dt * rp(i)
        pp = pp + 1
        end do

        ! reference implementation
        pp = ptrcur
        do i=1,np
            xbuf(1,i) = this%x(1,pp) + this%p(1,pp) * rpp(i) * rdx(1)
            xbuf(2,i) = this%x(2,pp) + this%p(2,pp) * rpp(i) * rdx(2)
            dxi(1,i) = ntrim_qed( xbuf(1,i) )
            dxi(2,i) = ntrim_qed( xbuf(2,i) )

            pp = pp + 1
        end do


        ! copy data from buffer to species data trimming positions
        pp = ptrcur
        do i = 1, np
            this%x(1,pp)  = xbuf(1,i) - dxi(1,i)
            this%x(2,pp)  = xbuf(2,i) - dxi(2,i)
            this%ix(1,pp) = this%ix(1,pp) + dxi(1,i)
            this%ix(2,pp) = this%ix(2,pp) + dxi(2,i)

            pp = pp + 1
        end do

    enddo

end subroutine advance_photons_2d
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine advance_photons_3d( this, gdt, i0, i1)
!---------------------------------------------------------------------------------------------------

    implicit none

    integer, parameter :: rank = 3

    ! dummy variables

    class( t_photons ),    intent(inout) :: this
    real(p_double),   intent(in) :: gdt
    integer,   intent(in) :: i0, i1

    real(p_k_part), dimension(p_x_dim) :: rdx
    integer :: i, pp, np, ptrcur

    real(p_k_part), dimension(rank,p_cache_size) :: xbuf
    integer, dimension(rank,p_cache_size)        :: dxi
    real(p_k_part), dimension(p_cache_size)      :: rp, rpp
    real(p_k_part) :: dt

    ! executable statements
    rdx(1) = real( 1.0_p_double/this%dx(1), p_k_part )
    rdx(2) = real( 1.0_p_double/this%dx(2), p_k_part )
    rdx(3) = real( 1.0_p_double/this%dx(3), p_k_part )
    dt = real( gdt, p_k_part )

    ! advance position of photons i0 to i1 in chunks of p_cache_size
    do ptrcur = i0, i1, p_cache_size

        ! check if last copy of table and set np
        if( ptrcur + p_cache_size > i1 ) then
            np = i1 - ptrcur + 1
        else
            np = p_cache_size
        endif

        ! this type of loop is actually faster than
        ! using a forall construct
        pp = ptrcur

        do i=1,np
            rp(i) = 1.0_p_k_part / &
            sqrt( this%p(1,pp)**2 + this%p(2,pp)**2 + this%p(3,pp)**2 )
            rpp(i) = dt * rp(i)
            pp = pp + 1
        end do

        ! reference implementation
        pp = ptrcur
        do i=1,np
            xbuf(1,i) = this%x(1,pp) + this%p(1,pp) * rpp(i) * rdx(1)
            xbuf(2,i) = this%x(2,pp) + this%p(2,pp) * rpp(i) * rdx(2)
            xbuf(3,i) = this%x(3,pp) + this%p(3,pp) * rpp(i) * rdx(3)
            dxi(1,i) = ntrim_qed( xbuf(1,i) )
            dxi(2,i) = ntrim_qed( xbuf(2,i) )
            dxi(3,i) = ntrim_qed( xbuf(3,i) )

            pp = pp + 1
        end do


        ! copy data from buffer to species data trimming positions
        pp = ptrcur
        do i = 1, np
            this%x(1,pp)  = xbuf(1,i) - dxi(1,i)
            this%x(2,pp)  = xbuf(2,i) - dxi(2,i)
            this%x(3,pp)  = xbuf(3,i) - dxi(3,i)
            this%ix(1,pp) = this%ix(1,pp) + dxi(1,i)
            this%ix(2,pp) = this%ix(2,pp) + dxi(2,i)
            this%ix(3,pp) = this%ix(3,pp) + dxi(3,i)

            pp = pp + 1
        end do
    enddo

end subroutine advance_photons_3d
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine advance_photons_cyl_2d( this, gdt, i0, i1 )
!---------------------------------------------------------------------------------------------------  

    implicit none

    integer, parameter :: rank = 2

    class( t_photons ), intent(inout) :: this
    real(p_double), intent(in) :: gdt
    integer, intent(in) :: i0, i1

    integer :: i, pp, np, ptrcur
    real(p_double), dimension(p_x_dim)   :: xmin_g
    real(p_k_part), dimension(rank,p_cache_size) :: xbuf
    integer, dimension(rank,p_cache_size)        :: dxi
    real(p_k_part), dimension(p_cache_size)      :: rgamma

    real(p_k_part) :: dt_dx1
    integer :: gix2, shift_ix2
    real(p_double) :: dr, rdr, tmp

    real(p_double) :: x2_new, x3_new, r_old, r_new
    real(p_double) :: r_shift

    ! executable statements

    xmin_g(1) = this%g_box( p_lower, 1 )
    xmin_g(2) = this%g_box( p_lower, 2 )

    shift_ix2 = this%my_nx_p(p_lower, 2) - 2
    dr  = this%dx(p_r_dim)
    rdr = 1.0_p_double/dr

    dt_dx1 = real( gdt / this%dx(1), p_k_part )

    ! Since r = 0 is at the center of cell ix2 = 1 the radial position of particles will be
    ! - Positions defined with regard to the center of the cell (odd interpolation)
    !       r = ( (gix2-1) + x2 ) * dr
    ! - Positions defined with regard to the corner of the cell (even interpolation)
    !       r = ( (gix2-1) + x2 - 0.5 ) * dr

    if ( this%pos_type == p_cell_near ) then
        r_shift = 0.5_p_double
    else
        r_shift = 0.0_p_double
    endif

    ! advance position
    do ptrcur = i0, i1, p_cache_size

        ! check if last copy of table and set np
        if( ptrcur + p_cache_size > i1 ) then
            np = i1 - ptrcur + 1
        else
            np = p_cache_size
        endif

        ! this type of loop is actually faster than
        ! using a forall construct
        pp = ptrcur

        do i=1,np
            rgamma(i) = 1.0_p_k_part / &
            sqrt( ( ( this%p(1,pp)**2 ) + this%p(2,pp)**2 ) + this%p(3,pp)**2 )

            pp = pp + 1
        end do

        ! advance particle position
        pp = ptrcur

        do i=1,np
            xbuf(1,i) = this%x(1,pp) + ( this%p(1,pp) * rgamma(i) ) * dt_dx1

            ! Convert radial "cell" position to "box" position in double precision
            gix2  = this%ix(2,pp)  + shift_ix2
            r_old = ( this%x(2,pp) + (gix2 - r_shift) ) * dr

            ! Particle is pushed in 3D like fashion
            x2_new    = r_old + ( this%p(2,pp) * rgamma(i) ) * gdt
            x3_new    =         ( this%p(3,pp) * rgamma(i) ) * gdt

            r_new     = sqrt( x2_new**2 + x3_new**2 )

            ! Convert new position to "cell" type position

            ! this is a protection against roundoff for cold plasmas
            if ( r_old == r_new ) then
            xbuf(2,i) = this%x(2,pp)
            else
            xbuf(2,i) = real( r_new * rdr - ( gix2 - r_shift ) , p_k_part)
            endif

            ! Correct p_r and p_\theta to conserve angular momentum
            ! there is a potential division by zero here
            tmp          = 1.0_p_double / r_new
            this%p(2,pp) = real( ( this%p(2,pp) * x2_new + this%p(3,pp)*x3_new ) * tmp, p_k_part )
            this%p(3,pp) = real( this%p(3,pp) * ( r_old * tmp ), p_k_part )

            dxi(1,i) = ntrim_qed( xbuf(1,i) )
            dxi(2,i) = ntrim_qed( xbuf(2,i) )

            pp = pp + 1
        end do

        ! copy data from buffer to species data trimming positions
        pp = ptrcur
        do i = 1, np
        this%x(1,pp)  = xbuf(1,i) - dxi(1,i)
        this%x(2,pp)  = xbuf(2,i) - dxi(2,i)
        this%ix(1,pp) = this%ix(1,pp) + dxi(1,i)
        this%ix(2,pp) = this%ix(2,pp) + dxi(2,i)

        pp = pp + 1
        end do

    enddo


end subroutine advance_photons_cyl_2d
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
function ntrim_qed(x)
!---------------------------------------------------------------------------------------------------
! Returns the integer shift (-1, 0 or +1) so that the coordinate remains in the [-0.5, 0.5[
! range. This is the fastest implementation (twice as fast as a sequence of ifs) because
! the two if structures compile as conditional moves and can be processed independently.
! This has no precision problem and is only 12% slower than the previous "int(x+1.5)-1"
! routine that would break for x = nearest( 0.5, -1.0 )
!---------------------------------------------------------------------------------------------------
    implicit none

    real(p_k_part), intent(in) :: x
    integer :: ntrim_qed, a, b

    if ( x < -.5 ) then
    a = -1
    else
    a = 0
    endif

    if ( x >= .5 ) then
    b = +1
    else
    b = 0
    endif

    ntrim_qed = a+b

end function ntrim_qed
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
subroutine create_pairs( this, electrons, positrons, emf, dt, i0, i1)
!---------------------------------------------------------------------------------------------------

    use m_species_define_qed, only : t_species_qed
    use m_species_comm, only : t_part_idx, remove_particles
    use m_photons_diagnostics
    !use m_cross_section

    implicit none

    integer, parameter :: rank = 2

    ! dummy variables

    class( t_photons ),    intent(inout) :: this
    class( t_species_qed ),    intent(inout) :: electrons
    class( t_species_qed ),    intent(inout) :: positrons
    class( t_emf ), intent( in )  ::  emf
    real(p_double),   intent(in) :: dt
    integer,   intent(in) :: i0, i1

    real(p_k_part), dimension(p_p_dim,p_cache_size) :: bp, ep
    integer :: i, pp, np, ptrcur

    ! attributes for new particles
    real(p_k_part)   , dimension(p_p_dim) :: p_new

    ! index of photons to remove
    type(t_part_idx) :: local_idx

    !  QED local variables
    real(p_double) :: chi_g, chi_g_s
    real(p_double) :: coef_QED, omega_p, varand, Wpair, eta, interp, norm_schw, &
    enephot, pdotE2, enephotE_plus_pcrossB2, gam_e, gam_p

    ! dummy variables for chi parameter diagnostic
    real(p_k_part) :: q_g
    integer, dimension(p_x_dim) :: ix_g

    ! Physical constants taken from NIST reference data
    real(p_double), parameter :: hbar  = 1.054571726d-34     ! J s
    real(p_double), parameter :: me    = 9.10938291d-31      ! kg
    real(p_double), parameter :: cl    = 299792458.0         ! m/s
    real(p_double), parameter :: alpha = 7.2973525698d-3     ! (no units)

    real(p_double), parameter :: PI    = 3.141592653589793d0

    ! executable statements
    omega_p = electrons%omega_p0
    norm_schw = hbar * omega_p/(me * cl**2)
    coef_QED = alpha/(sqrt(3.d0)*PI*norm_schw)

    ! advance position of photons i0 to i1 in chunks of p_cache_size
    do ptrcur = i0, i1, p_cache_size

        ! check if last copy of table and set np
        if( ptrcur + p_cache_size > i1 ) then
            np = i1 - ptrcur + 1
        else
            np = p_cache_size
        endif

        call this % get_emf( emf, bp, ep, np, ptrcur )

        ! this type of loop is actually faster than
        ! using a forall construct
        pp = ptrcur
        do i=1, np
            ! calculate the quantum parameter chi for all photons
            pdotE2=((this%p(1,pp)*ep(1,i))+(this%p(2,pp)*ep(2,i))+(this%p(3,pp)*ep(3,i)))**2

            enephot = sqrt(this%p(1,pp)**2+this%p(2,pp)**2+this%p(3,pp)**2)

            enephotE_plus_pcrossB2=(-enephot*ep(1,i)-bp(3,i)*this%p(2,pp)+bp(2,i)*this%p(3,pp))**2 &
            + (-enephot*ep(2,i)+bp(3,i)*this%p(1,pp)-bp(1,i)*this%p(3,pp))**2 &
            + (-enephot*ep(3,i)-bp(2,i)*this%p(1,pp)+bp(1,i)*this%p(2,pp))**2

            chi_g = sqrt(abs(enephotE_plus_pcrossB2-pdotE2))*norm_schw
            chi_g_s = chi_g * this%qed_zeta

            ! condition to radiate significant number of pairs
            if ( (chi_g_s > this%pp_chi_cutoff) .and. (chi_g_s < 5.0d4) .and.(enephot > 2.0d0)  ) then

                interp = interp_spec_pairQED(chi_g_s)
                Wpair = coef_QED*interp/enephot

                call random_number( varand )

                if ( varand < Wpair * dt ) then
                    ! check to be sure that the pair created have at least each an energy > mc^2
                    gam_e = 0.0d0
                    gam_p = 0.0d0

                    eta = specpair_fast(chi_g_s)
                    gam_e = enephot*eta/chi_g_s
                    gam_p = enephot - gam_e

                    q_g = this%q(pp)
                    ix_g = this%ix(1:p_x_dim,pp)

                    ! Diagnostics
                    if (this%ndump_fac_chi > 0) call chi_hist( this, ix_g, real(q_g, p_k_part), real(chi_g_s, p_k_part) )

                    if ( (gam_e.gt.1.0d0) .and. (gam_p.gt.1.0d0) ) then

                        this%num_pairs = this%num_pairs + 1
                        this%weight_pairs = this%weight_pairs + q_g
                        this%energy_pairs = this%energy_pairs + q_g*(gam_e + gam_p)

                        ! create electron
                        p_new  = real( sqrt((gam_e)**2-1)*this%p(:,pp)/enephot, p_k_part )
                        call electrons%create_particle( this%ix(:,pp), this%x(:,pp), p_new, -1.0_p_k_part*q_g )

                        ! create positron
                        p_new = real( sqrt((gam_p)**2-1)*this%p(:,pp)/enephot, p_k_part )
                        call positrons%create_particle( this%ix(:,pp), this%x(:,pp), p_new, +1.0_p_k_part*q_g )

                        ! register emission
                        if (this%ndump_fac_n_emit > 0) call n_emit_hist( this, this%ix(1:p_x_dim,pp) )

                    endif

                    ! mark photon for deletion
                    call add_particle_index( local_idx, pp )

                endif
            endif
            pp = pp + 1
        enddo
    enddo

    ! delete photons
    call remove_particles( this%t_species , local_idx )

    ! cleanup the local_idx buffer for deleted photons
    if ( local_idx%buf_size > 0 ) call freemem( local_idx%idx )
    local_idx%start = 1
    local_idx%nidx = 0
    local_idx%buf_size = 0

    contains

    !---------------------------------------------------------------------------------------------------
    subroutine add_particle_index( list, i )

        use m_species_define

        implicit none

        type(t_part_idx), intent(inout) :: list
        integer, intent(in)              :: i

        integer, dimension(:), pointer :: tmp_idx => null()

        ! grow buffer if necessary
        if ( list%nidx + 1 > list%buf_size ) then
            list%buf_size = list%buf_size + p_spec_buf_block
            call alloc(  tmp_idx, (/ list%buf_size /))
            if ( list%nidx > 0 ) call memcpy( tmp_idx, list%idx, list%nidx )
            call freemem( list%idx )
            list%idx => tmp_idx
        endif

        ! add particle index to send index list
        list%nidx = list%nidx + 1
        list%idx( list%nidx ) = i

    end subroutine add_particle_index
    !---------------------------------------------------------------------------------------------------

end subroutine create_pairs
!---------------------------------------------------------------------------------------------------

  subroutine create_pairs_betheheitler( this, electrons, positrons, dt, i0, i1)

    use m_species_define_qed, only : t_species_qed
    use m_species_define, only : t_species
    use m_species_comm, only : t_part_idx, remove_particles
    use m_photons_diagnostics

    implicit none

    ! dummy variables

    class( t_photons ),    intent(inout) :: this
    class( t_species_qed ),    intent(inout) :: electrons
    class( t_species_qed ),    intent(inout) :: positrons
    real(p_double),   intent(in) :: dt
    integer,   intent(in) :: i0, i1

    integer :: i, pp, np, ptrcur

    ! attributes for new particles
    real(p_k_part)   , dimension(p_p_dim) :: p_new

    ! index of photons to remove
    type(t_part_idx) :: local_idx
    real(p_double) :: proba, varand, enephot, gam_e, gam_p

    ! Bethe Heitler variables
    real(p_double) :: factor_mult, omega_p, n_p, ni 
    real(p_k_part) :: q_pair
    
    integer :: Z_ion, cs_model
    logical :: decay_photons

    ! We introduce shorter names for all Bethe Heitler variables of interest
    factor_mult = this % betheheitler_info % proba_mult
    omega_p = this%omega_p0
    n_p = this%n0
    Z_ion = int(this % betheheitler_info % Z_ion)
    cs_model = int(this % betheheitler_info % cs_model)
    decay_photons = this % betheheitler_info % decay_photons

    ! advance position of photons i0 to i1 in chunks of p_cache_size
    do ptrcur = i0, i1, p_cache_size

        ! check if last copy of table and set np
        if( ptrcur + p_cache_size > i1 ) then
            np = i1 - ptrcur + 1
        else
            np = p_cache_size
        endif

        ! this type of loop is actually faster than
        ! using a forall construct
        pp = ptrcur
        do i=1, np

                ! Photon energy
                enephot = sqrt(this%p(1,pp)**2+this%p(2,pp)**2+this%p(3,pp)**2)

                ! The photon energy must be at least 2 mc2 for pair production
                if (enephot .gt. 2.0d0) then 

                    ! Ion density at the position of the photon
                    select case ( p_x_dim )
                    case ( 1 )
                        ni = this % betheheitler_info % rho_ion % f1(1, this%ix(1, pp))
                    case ( 2 )
                        ni = this % betheheitler_info % rho_ion % f2(1, this%ix(1, pp), this%ix(2, pp))
                    case ( 3 )
                        ni = this % betheheitler_info % rho_ion % f3(1, this%ix(1, pp), this%ix(2, pp), this%ix(3, pp))
                    end select

                    call bh_proba(Z_ion, enephot, ni, dt, omega_p, n_p, factor_mult, proba)
                    call random_number( varand )
                    
                    if ( varand < proba ) then
                        ! check to be sure that the pair created have at least each an energy > mc^2
                        gam_e = 0.0d0
                        gam_p = 0.0d0

                        ! Sample positron energy and deduce electron energy
                        call bh_spec(Z_ion, enephot, gam_p)
                        gam_e = enephot - gam_p
                        q_pair = this%q(pp) / factor_mult

                        ! Diagnostics
                        if ( (gam_e.gt.1.0d0) .and. (gam_p.gt.1.0d0) ) then

                            this%num_pairs_bh = this%num_pairs_bh + 1
                            this%weight_pairs_bh = this%weight_pairs_bh + q_pair
                            this%energy_pairs_bh = this%energy_pairs_bh + q_pair*(gam_e + gam_p)

                            ! create electron
                            p_new  = real( sqrt((gam_e)**2-1)*this%p(:,pp)/enephot, p_k_part )
                            call electrons%create_particle_single_cell_p( this%ix(:,pp), this%x(:,pp), p_new, -1.0_p_k_part*q_pair )

                            ! create positron
                            p_new = real( sqrt((gam_p)**2-1)*this%p(:,pp)/enephot, p_k_part )
                            call positrons%create_particle_single_cell_p( this%ix(:,pp), this%x(:,pp), p_new, 1.0_p_k_part*q_pair )

                        endif

                        ! Decay of the photon (if we increased proba, we don't decay the photon, but just decrease its weight) 
                        if (decay_photons) then
                            if (factor_mult .eq. 1.0_p_double) then
                                call add_particle_index( local_idx, pp )
                            else
                                this%q(pp) = this%q(pp) * (1.0_p_k_part - (gam_e + gam_p - 2.)/(factor_mult*enephot))
                            endif
                        endif

                    endif
                endif
                pp = pp + 1
        enddo
    enddo

    ! delete photons
    if (decay_photons) call remove_particles( this%t_species , local_idx )

    ! cleanup the local_idx buffer for deleted photons
    if ( local_idx%buf_size > 0 ) call freemem( local_idx%idx )
    local_idx%start = 1
    local_idx%nidx = 0
    local_idx%buf_size = 0

    contains

    !---------------------------------------------------------------------------------------------------
    subroutine add_particle_index( list, i )

        use m_species_define

        implicit none

        type(t_part_idx), intent(inout) :: list
        integer, intent(in)              :: i

        integer, dimension(:), pointer :: tmp_idx => null()

        ! grow buffer if necessary
        if ( list%nidx + 1 > list%buf_size ) then
            list%buf_size = list%buf_size + p_spec_buf_block
            call alloc(  tmp_idx, (/ list%buf_size /))
            if ( list%nidx > 0 ) call memcpy( tmp_idx, list%idx, list%nidx )
            call freemem( list%idx )
            list%idx => tmp_idx
        endif

        ! add particle index to send index list
        list%nidx = list%nidx + 1
        list%idx( list%nidx ) = i

    end subroutine add_particle_index
    !---------------------------------------------------------------------------------------------------


  end subroutine create_pairs_betheheitler

end module m_photons

!---------------------------------------------------------------------------------------------------
! advance photons
!---------------------------------------------------------------------------------------------------
subroutine advance_photons( this, electrons, positrons, emf, gdt, tid, n_threads)
    
    use m_system
    use m_parameters
    use m_photons_define, only : t_photons
    use m_species_define, only : p_ene_recalc
    use m_species_define_qed, only : t_species_qed
    use m_emf_define, only: t_emf
    use m_photons
    use m_species_define, only : p_ene_recalc
    use m_species, only : sort

    implicit none
  
    class( t_photons ), intent(inout) :: this
    class( t_species_qed ), intent(inout) :: electrons
    class( t_species_qed ), intent(inout) :: positrons
    class( t_emf ), intent( inout )  ::  emf
    real(p_double), intent(in) :: gdt
    integer, intent(in) :: tid    ! local thread id
    integer, intent(in) :: n_threads  ! total number of threads
  
    integer :: chunk, ip0, ip1

    ! sanity check
    ! Code does not yet support OpenMP parallelism (mainly because of pair production)
    if ( n_threads > 1 ) then
      write(0,*) 'QED photon pusher is not yet supported with multiple threads per node'
      call abort_program()
    endif

    ! Time centered energy diagnostic is not yet implemented, force recalculation later
    ! when get_energy() is called
    ! (the push should not change photon energy)
    this%energy = p_ene_recalc
  
    ! range of photons for each thread
    chunk = ( this%num_par + n_threads - 1 ) / n_threads
    ip0    = tid * chunk + 1
    ip1    = min( (tid+1) * chunk, this%num_par )

    ! Create pairs and absorb photons
    if (this%if_pairprod) call create_pairs(this, electrons, positrons, emf, gdt, ip0, ip1)

    ! Warning : we want to avoid looping on photons that already decayed into pairs
    ip1    = min( (tid+1) * chunk, this%num_par )

    ! call bethe heitler pair creation
    if (this % betheheitler_info % if_betheheitler) then
        call create_pairs_betheheitler(this, electrons, positrons, gdt, ip0, ip1)
    endif

    ! Push photons. Boundary crossings will be checked at update_boundary
    select case ( this%coordinates )
        case default
            select case ( p_x_dim )
                case (1)
                    call advance_photons_1d( this, gdt, ip0, ip1 )

                case (2)
                    call advance_photons_2d( this, gdt, ip0, ip1 )
            
                case (3)
                    call advance_photons_3d( this, gdt, ip0, ip1 )
            
                case default
                    ERROR('Not implemented for x_dim = ',p_x_dim)
                    call abort_program(p_err_invalid)
        
            end select
        case ( p_cylindrical_b )
            call advance_photons_cyl_2d( this, gdt, ip0, ip1 )
    end select
  
end subroutine advance_photons
!---------------------------------------------------------------------------------------------------