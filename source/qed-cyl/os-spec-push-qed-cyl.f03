! m_species_push_qedcyl module
!
! Handles the species advance/deposit for the quasi-3D + QED algorithm

#include "os-config.h"
#include "os-preprocess.fpp"

module m_species_push_qedcyl

use m_system
use m_parameters
use m_math
use m_species_cyl_modes_define, only: t_species_cyl_modes
use m_species_define_qedcyl,       only: t_species_qedcyl
use m_species_define,           only: t_species, p_cell_low, p_cell_near, incr_linear, incr_regressive
use m_emf_define,               only: t_emf
use m_vdf_define,               only: t_vdf
use m_time_step,                only: t_time_step, dt, n
use m_species_current_cyl_modes,only: get_coef, getjr_cyl_m_s1, getjr_cyl_m_s2, getjr_cyl_m_s3
use m_cyl_modes,                only: t_cyl_modes
use m_qed
use m_qed_coulomb

implicit none

private

interface dudt_boris_qedcyl
  module procedure dudt_boris_qedcyl
end interface

interface dudt_exact_qedcyl
  module procedure dudt_exact_qedcyl
end interface

public :: dudt_boris_qedcyl, dudt_exact_qedcyl

integer, parameter :: p_iter_max = 10

contains

!-----------------------------------------------------------------------------------------
subroutine dudt_boris_qedcyl( this, emf, dt, i0, i1, energy, t )
!-----------------------------------------------------------------------------------------
! Advance velocities and calculate time-centered total energy
! process particles in the range [i0,i1]
! This routine is meant to be used for testing simd code and is written to use the same numerics
!-----------------------------------------------------------------------------------------
    
    implicit none

    class( t_species ), intent(inout) :: this
    class( t_emf ), intent(in) ::  emf
    real(p_double), intent(in) :: dt
    integer, intent(in) :: i0, i1
    real(p_double), intent(inout) :: energy
    real(p_double), intent(in) :: t

    ! dummy variables
    real(p_k_part) :: tem
    integer :: i, ptrcur, np, pp
    real(p_k_part), dimension(p_p_dim,p_cache_size) :: bp, ep, utemp
    real(p_k_part), dimension(p_cache_size) ::  gam_tem, otsq
    real(p_k_part) :: gamma, u2
    real(p_double) :: loc_ene

    ! qed specific variables
    real(p_double), dimension(p_cache_size) :: gl, k_damp
    real(p_k_part), dimension(p_p_dim,p_cache_size) :: bpm, epm, um

  ! executable statements

    select type( this )
    class is ( t_species_qedcyl )

  ! get factor that includes timestep and charge-to-mass ratio
  ! note that the charge-to-mass ratio is used since the
  ! momentum p is not the total momentum of a particles but
  ! the momentum per unit restmass (which is the electron mass)

    tem = real( 0.5_p_double * dt / this%rqm, p_k_part )

    !loop through all particles
    do ptrcur = i0, i1, p_cache_size

    ! check if last copy of table and set np
    if( ptrcur + p_cache_size > i1 ) then
        np = i1 - ptrcur + 1
    else
        np = p_cache_size
    endif

    ! modify bp & ep to include timestep and charge-to-mass ratio
    ! and perform half the electric field acceleration.
    ! Result is stored in UTEMP.

    call this % get_emf( emf, bp, ep, np, ptrcur, t=t )

    call emit_photons_qedcyl( this, ep, bp, dt, gl, k_damp, ptrcur, np )

    if (this%if_damp_classical) then
        ! copy the momentum at t = n-1/2
        pp = ptrcur
        do i=1,np
            um(:,i) = this%p(:,pp)
            pp = pp+1
        enddo

        ! copy the electric and magnetic fields
        do i=1, np
        epm(:,i) = ep(:,i)
        bpm(:,i) = bp(:,i)
        end do
    endif

    do i=1, np
        ep(1,i) = ep(1,i) * tem
        ep(2,i) = ep(2,i) * tem
        ep(3,i) = ep(3,i) * tem
    enddo

    loc_ene = 0

    pp = ptrcur
    do i=1,np
        utemp(1,i) = this%p(1,pp) + ep(1,i)
        utemp(2,i) = this%p(2,pp) + ep(2,i)
        utemp(3,i) = this%p(3,pp) + ep(3,i)

        ! Get time centered gamma
        u2 =  (utemp(1,i)**2 + utemp(2,i)**2) + utemp(3,i)**2
        gamma = sqrt( u2 + 1 )

        ! accumulate time centered energy
        ! this is done in double precicion always
        loc_ene = loc_ene + this%q(pp) * u2 / (gamma + 1.0_p_double)

        gam_tem(i)= tem / gamma

        pp = pp + 1
    enddo

    ! accumulate global energy
    energy = energy + loc_ene

    do i=1,np
        bp(1,i) = bp(1,i)*gam_tem(i)
        bp(2,i) = bp(2,i)*gam_tem(i)
        bp(3,i) = bp(3,i)*gam_tem(i)
    enddo

    pp = ptrcur
    do i=1,np
        this%p(1,pp) = utemp(1,i) + utemp(2,i) * bp(3,i)
        this%p(2,pp) = utemp(2,i) + utemp(3,i) * bp(1,i)
        this%p(3,pp) = utemp(3,i) + utemp(1,i) * bp(2,i)
        pp = pp + 1
    enddo

    pp = ptrcur
    do i=1,np
        this%p(1,pp) = this%p(1,pp) - utemp(3,i) * bp(2,i)
        this%p(2,pp) = this%p(2,pp) - utemp(1,i) * bp(3,i)
        this%p(3,pp) = this%p(3,pp) - utemp(2,i) * bp(1,i)
        pp = pp + 1
    enddo

    do i=1,np
        otsq(i) = 2.0_p_k_part / ( ((1.0_p_k_part + bp(1,i)**2) + bp(2,i)**2) + bp(3,i)**2)
    enddo

    do i=1,np
        bp(1,i) = bp(1,i) * otsq(i)
        bp(2,i) = bp(2,i) * otsq(i)
        bp(3,i) = bp(3,i) * otsq(i)
    enddo

    pp = ptrcur
    do i=1,np
        utemp(1,i) = utemp(1,i) + this%p(2,pp) * bp(3,i)
        utemp(2,i) = utemp(2,i) + this%p(3,pp) * bp(1,i)
        utemp(3,i) = utemp(3,i) + this%p(1,pp) * bp(2,i)
        pp = pp + 1
    enddo

    pp = ptrcur
    do i=1,np
        utemp(1,i) = utemp(1,i) - this%p(3,pp) * bp(2,i)
        utemp(2,i) = utemp(2,i) - this%p(1,pp) * bp(3,i)
        utemp(3,i) = utemp(3,i) - this%p(2,pp) * bp(1,i)
        pp = pp + 1
    enddo

    ! Perform second half of electric field acceleration.
    pp = ptrcur
    do i=1,np
        this%p(1,pp) = utemp(1,i) + ep(1,i)
        this%p(2,pp) = utemp(2,i) + ep(2,i)
        this%p(3,pp) = utemp(3,i) + ep(3,i)
        pp = pp + 1
    enddo

!    ! Specific part to photon emission
    if (this%if_damp_qed) call damp_qed_qedcyl( this, k_damp, ptrcur, np )
    if (this%if_damp_classical) call damp_classical_qedcyl( this, epm, bpm, dt, gl, um, ptrcur, np )

    ! bremsstrahlung photon emission
    if (this % bremsstrahlung_info % if_bremsstrahlung) then
        call bremsstrahlung_emission(this, dt, k_damp, ptrcur, i0, i1, np)
    endif

    ! Trident Coulomb pair creation
    if (this % trident_coul_info % if_trident_coul) then
        call trident_coul_pair_creation( this, dt, ptrcur, i0, i1, np)
    endif

  enddo

    class default
        ! this must never happen
        call abort_program( p_err_invalid )
    end select

end subroutine dudt_boris_qedcyl
! ----------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! Emits photons
!---------------------------------------------------------------------------------------------------
subroutine emit_photons_qedcyl( this, ep, bp, dt, gl, k_damp, ptrcur, np )
!---------------------------------------------------------------------------------------------------
    
    use m_random
    use m_species_diagnostics_qedcyl

    implicit none

    class( t_species_qedcyl ),    intent(inout) :: this
    real(p_k_part), dimension(p_p_dim,p_cache_size), intent(in) :: ep, bp
    real(p_double), intent(in) :: dt
    real(p_double), dimension(p_cache_size), intent(inout) :: gl, k_damp
    integer, intent(in) :: ptrcur, np

    ! dummy variables
    integer :: i, pp, n_emit
    real(p_double) :: eta, chi_g, coef_QED, norm_schw, p2, pdotE2, &
                      gamE_plus_pcrossB2, varand, Wrad, interp
    
    real(p_double), dimension(2) :: energy_rad

    real(p_k_part) :: p_frac, p_g, q_g
    real(p_k_part), dimension(p_p_dim) :: p_g_vec
    integer, dimension(p_x_dim) :: ix_g
    
    real( p_double ), parameter :: M_PI  = 3.14159265358979323846264338327950288d0 ! pi

    ! Values from the NIST reference on Constants, Units and Uncertainty
    ! http://physics.nist.gov/cuu/Constants/index.html
    real( p_double ), parameter :: hbar  = 1.054571800d-34 ! Planck constant over 2 pi
    real( p_double ), parameter :: cl    = 299792458.0d0   ! speed of light in vacuum [m/s]
    real( p_double ), parameter :: me    = 9.10938356d-31  ! electron rest mass [kg]
    real( p_double ), parameter :: alpha = 7.2973525664d-3 ! fine-structure constant

    ! executable statements

    ! calculate the quantum coefficient
    norm_schw = hbar*this%omega_p0/(me*cl**2.d0)
    coef_QED = sqrt(3.d0)*alpha/(2.d0*M_PI*norm_schw)

    ! reset energy radiation diagnostic
    energy_rad = 0.0_p_double

    ! reset damping factor
    k_damp = 1.0_p_double

    pp = ptrcur
    do i = 1, np    
        ! compute lorentz factor
        p2 = this%p(1,pp)**2+this%p(2,pp)**2+this%p(3,pp)**2
        gl(i) = sqrt(1.0d0+p2)

        ! compute eta parameter terms
        pdotE2 = (this%p(1,pp)*ep(1,i)+this%p(2,pp)*ep(2,i)+this%p(3,pp)*ep(3,i))**2
        gamE_plus_pcrossB2 = (-gl(i)*ep(1,i)-bp(3,i)*this%p(2,pp)+bp(2,i)*this%p(3,pp))**2 &
                           + (-gl(i)*ep(2,i)+bp(3,i)*this%p(1,pp)-bp(1,i)*this%p(3,pp))**2 &
                           + (-gl(i)*ep(3,i)-bp(2,i)*this%p(1,pp)+bp(1,i)*this%p(2,pp))**2 


        ! compute the quantum parameter eta
        eta = sqrt(abs(gamE_plus_pcrossB2-pdotE2))*norm_schw

        if ( (eta > 1.0d-4) .and. (eta < 5.0d4) .and. (gl(i) > this%qed_g_cutoff) ) then

            if (eta > 300.d0) then
                Wrad = (1.46d0*alpha/norm_schw)*(eta**(2.0d0/3.0d0))/gl(i)
            else
                interp = interp_spec_QED(eta)
                Wrad = interp*coef_QED*eta/gl(i)
            endif

            call random_number( varand )

            if ( varand < Wrad * dt ) then
                ! compute photon chi
                chi_g = find_photon_chi(eta)

                ! compute momentum fraction lost to emitted photon
                p_frac = 2.0d0*chi_g/eta
                k_damp(i) = real(1.0 - p_frac, p_k_part)

                ! compute photon momentum
                p_g = p_frac * gl(i)

                ! compute photon momentum vector
                p_g_vec = p_frac*this%p(:,pp)

                ! compute photon weight
                q_g = abs(this%q(pp))
                ix_g = this%ix(1:p_x_dim,pp)

                if (this%ndump_fac_radspect > 0) call radspect_hist_qedcyl( this, p_g, q_g )
                if (this%ndump_fac_radanglespect > 0) call radanglespect_hist_qedcyl( this, p_g_vec, q_g )
                if (this%ndump_fac_radsphere > 0) call radsphere_hist_qedcyl( this, p_g_vec, q_g )
                if (this%ndump_fac_raddetector>0) call detector_hist_qedcyl( this, p_g_vec, q_g )
                if (this%ndump_fac_chi > 0) call chi_hist_qedcyl( this, ix_g, real(q_g, p_k_part), real(eta, p_k_part) )
                if (this%ndump_fac_chi_emit > 0) call emission_chi_hist_qedcyl( this, real(eta, p_k_part), real(chi_g, p_k_part), real(q_g, p_k_part) )

                ! condition for photons that can emit a pair
                if (p_g > this%p_emit_cutoff) then

                  call this%photons%create_particle_single_cell_p_cyl( &
                          this%ix(1:p_x_dim,pp), &
                          this%x(:,pp), &
                          p_g_vec, q_g &
                        )
                  
                  ! register emission
                  if (this%ndump_fac_n_emit > 0) call n_emit_hist_qedcyl( this, ix_g )
                  
                  ! add photon energy to QED tracked bin
                  energy_rad(2) = energy_rad(2) +  p_g * q_g
                else
                  ! add photon energy to QED non-tracked bin
                  energy_rad(1) = energy_rad(1) +  p_g * q_g
              endif

            endif
        endif

        pp = pp + 1
    end do

    ! accumulate radiated energy
    this % energy_rad(2) = this % energy_rad(2) + energy_rad(1)
    this % energy_rad(3) = this % energy_rad(3) + energy_rad(2)
    
end subroutine emit_photons_qedcyl  
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! Damps particle momentum due to qed emission
!---------------------------------------------------------------------------------------------------
subroutine damp_qed_qedcyl( this, k_damp, ptrcur, np )
!---------------------------------------------------------------------------------------------------

    implicit none

    class( t_species_qedcyl ),    intent(inout) :: this
    real(p_double), dimension(p_cache_size), intent(in) :: k_damp
    integer, intent(in) :: ptrcur, np

    ! dummy variables
    integer :: i, pp
    
    ! QED radiation damping
    pp = ptrcur
    do i=1,np
        this%p(:,pp) = this%p(:,pp)*k_damp(i)
        pp = pp + 1
    end do

end subroutine damp_qed_qedcyl 
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! Damps particle momentum due to classical radiation reaction
!---------------------------------------------------------------------------------------------------
subroutine damp_classical_qedcyl( this, ep, bp, dt, gl, um, ptrcur, np )
!---------------------------------------------------------------------------------------------------
    
    use m_random
    use m_species_diagnostics_qedcyl
    
    implicit none

    class( t_species_qedcyl ),    intent(inout) :: this
    real(p_double), intent(in) :: dt
    real(p_k_part), dimension(p_p_dim,p_cache_size), intent(in) :: ep, bp, um
    real(p_double), dimension(p_cache_size), intent(in) :: gl
    integer, intent(in) :: ptrcur, np

    ! dummy variables
    integer :: i, pp
    real(p_double) :: norm_schw, p2, gf, energy_rad
    real(p_k_part), dimension(p_p_dim,p_cache_size) :: un, frad, frad2, frad3, &
        frad2_aux, frad3_aux
    real(p_k_part), dimension(p_cache_size) :: uE, uB, B2, gamn, inv_gamn, gi
    real(p_k_part), dimension(p_p_dim) :: p_g_vec
    real(p_k_part) :: q_g    

    ! Values from the NIST reference on Constants, Units and Uncertainty
    ! http://physics.nist.gov/cuu/Constants/index.html
    real( p_double ), parameter :: hbar  = 1.054571800d-34 ! Planck constant over 2 pi
    real( p_double ), parameter :: cl    = 299792458.0d0   ! speed of light in vacuum [m/s]
    real( p_double ), parameter :: me    = 9.10938356d-31  ! electron rest mass [kg]
    real( p_double ), parameter :: alpha = 7.2973525664d-3 ! fine-structure constant

    ! executable statements

    ! calculate the quantum coefficient
    norm_schw = hbar*this%omega_p0/(me*cl**2.d0)
    
    if (this%ndump_fac_rad > 0) then
      pp = ptrcur
      do i=1,np
        if (gl(i) < this%qed_g_cutoff) then
          gi(i) = sqrt(1.0_p_k_part + this%p(1,pp)**2 + this%p(2,pp)**2 + this%p(3,pp)**2)
        endif
        pp = pp + 1
      end do
    endif

    ! get time centered momentum and Lorentz factor
    pp = ptrcur
    do i=1,np
      if (gl(i) < this%qed_g_cutoff) then
        un(:,i) = 0.5_p_k_part*(um(:,i) + this%p(:,pp))
        gamn(i) = sqrt(1.0_p_k_part + un(1,i)**2 + un(2,i)**2 + un(3,i)**2)
        inv_gamn(i)= 1.0_p_k_part / gamn(i)
      endif
      pp = pp + 1
    end do

    ! compute auxiliary dot products
    do i=1,np
      if (gl(i) < this%qed_g_cutoff) then
        uE(i) = dot_product(un(:,i), ep(:,i))
        uB(i) = dot_product(un(:,i), bp(:,i))
        B2(i) = dot_product(bp(:,i), bp(:,i))         
      endif
    end do

    do i=1,np
      if (gl(i) < this%qed_g_cutoff) then
        frad2_aux(:,i) = bp(:,i)*uB(i)-un(:,i)*B2(i)+ep(:,i)*uE(i)
        frad3_aux(:,i) = (gamn(i)*ep(1,i)+un(2,i)*bp(3,i)-un(3,i)*bp(2,i))**2
      endif
    end do

    do i=1,np
      if (gl(i) < this%qed_g_cutoff) then
        frad3_aux(:,i) = frad3_aux(:,i) + (gamn(i)*ep(2,i)+un(3,i)*bp(1,i)-un(1,i)*bp(3,i))**2
      endif
    end do

    do i=1,np
      if (gl(i) < this%qed_g_cutoff) then
        frad3_aux(:,i) = frad3_aux(:,i) + (gamn(i)*ep(3,i)+un(1,i)*bp(2,i)-un(2,i)*bp(1,i))**2
      endif
    end do

    do i=1,np
      if (gl(i) < this%qed_g_cutoff) then
        frad3_aux(:,i) = frad3_aux(:,i) - uE(i)**2
      endif
    end do

    do i=1,np
      if (gl(i) < this%qed_g_cutoff) then
        frad3_aux(:,i) = frad3_aux(:,i) * inv_gamn(i)
      endif
    end do

    ! calculate full radiation force from L&L
    pp = ptrcur
    do i=1,np
      if (gl(i) < this%qed_g_cutoff) then
        ! we neglect the first term coming from space time gradients of E and B
        ! second term
        frad2(1,i) = -ep(2,i)*bp(3,i)+ep(3,i)*bp(2,i)- &
                      inv_gamn(i)*frad2_aux(1,i)

        frad2(2,i) = -ep(3,i)*bp(1,i)+ep(1,i)*bp(3,i)- &
                      inv_gamn(i)*frad2_aux(2,i)

        frad2(3,i) = -ep(1,i)*bp(2,i)+ep(2,i)*bp(1,i)- &
                      inv_gamn(i)*frad2_aux(3,i)

        ! third term
        frad3(:,i) = un(:,i)*frad3_aux(:,i)

        frad(:,i) = (0.6666666666666666_p_k_part*alpha*norm_schw)*(frad2(:,i)+frad3(:,i))

        ! calculate classical recoil from L&L
        this%p(:,pp) = this%p(:,pp) - dt*frad(:,i)
      endif
      pp = pp + 1
    end do

    if (this%ndump_fac_rad > 0) then

      energy_rad = 0.0_p_double

      pp = ptrcur
      do i=1,np
        if (gl(i) < this%qed_g_cutoff) then
          ! compute final lorentz factor
          p2 = this%p(1,pp)**2+this%p(2,pp)**2+this%p(3,pp)**2
          gf = sqrt(1.0d0+p2)
          
          energy_rad = energy_rad + (gi(i) - gf) * this%q(pp)
          if (this%if_add_classical_radsphere) then
              if (this%ndump_fac_radsphere > 0) then
                p_g_vec = (this%p(1:3,pp))*(gi(i) - gf) / p2
                q_g=abs(this%q(pp))
                call radsphere_hist_qedcyl( this, p_g_vec, q_g )
              endif
          endif
        endif
        pp = pp + 1
      end do
      
      ! accumulate radiated energy
      this%energy_rad(1) = this%energy_rad(1) + energy_rad

    endif
  
end subroutine damp_classical_qedcyl 
!---------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine dudt_exact_qedcyl( this, emf, dt, i0, i1, energy, time )
!-----------------------------------------------------------------------------------------
! Advance velocities and calculate time-centered total energy
! process particles in the range [i0,i1]
! This pusher follows the ref. https://arxiv.org/abs/2007.07556
!-----------------------------------------------------------------------------------------

  implicit none

  class( t_species ), intent(inout) :: this
  class( t_emf ), intent( in ) :: emf
  real(p_double), intent(in) :: dt
  integer, intent(in) :: i0, i1
  real(p_double), intent(inout) :: energy
  real(p_double), intent(in) :: time

  integer :: i, j, ptrcur, np, pp
  real(p_k_part), dimension(p_p_dim,p_cache_size) :: bp, ep
  real(p_k_part), dimension(4) :: u0, u, fu0, f2u0, f3u0, uk0, fuk0, uk, uo0, fuo0, uo
  real(p_k_part) :: t, cn_o, snc_o, cnc_o, ch_k, shc_k, chc_k, tmp, k2o2, c1, c2
  real(p_k_part) :: inv1, inv2, e2, b2, fac_rr, qm
  real(p_double) :: loc_ene
  real(p_k_part), dimension(p_cache_size) :: k, o, k2, o2
  integer, dimension(p_cache_size) :: ix1, ix2
  integer :: n1, n2

  real(p_double), save :: rel_tol = 1.0d-12
  real(p_double), save :: weak_fld_tol = 1.0d-12

  ! qed specific variables
  real(p_double), dimension(p_cache_size) :: gl, k_damp
  real(p_double), dimension(p_cache_size) :: g_elec
  real(p_double) :: g_elec_tmp
  real(p_double), dimension(2) :: energy_rad

  select type( this )
  class is ( t_species_qedcyl )

  ! executable statements
  qm = 1.0_p_k_part / this%rqm
  fac_rr = this%k_rr * this%q_real * dt * 0.5_p_k_part * qm

  !loop through all particles
  do ptrcur = i0, i1, p_cache_size

    ! check if last copy of table and set np
    if( ptrcur + p_cache_size > i1 ) then
      np = i1 - ptrcur + 1
    else
      np = p_cache_size
    endif

    call this % get_emf( emf, bp, ep, np, ptrcur, t=time )
    call emit_photons_qedcyl( this, ep, bp, dt, gl, k_damp, ptrcur, np )

    loc_ene = 0.0_p_double
    energy_rad = 0.0_p_double

    pp = ptrcur
    do i = 1, np
      ep(1,i) = ep(1,i) * qm
      ep(2,i) = ep(2,i) * qm
      ep(3,i) = ep(3,i) * qm
      bp(1,i) = bp(1,i) * qm
      bp(2,i) = bp(2,i) * qm
      bp(3,i) = bp(3,i) * qm
      pp = pp + 1
    enddo

    ! -------------------------------------------------------------------------
    ! half-dt push of radiation reaction
    ! -------------------------------------------------------------------------
    if ( this%rad_react ) then
      pp = ptrcur
      do i = 1, np

        u0(1:3) = this%p(1:3,pp)
        u0(4) = sqrt( 1.0_p_k_part + u0(1)*u0(1) + u0(2)*u0(2) + u0(3)*u0(3) )

        fu0  = get_fu( ep(:,i), bp(:,i), u0 )
        f2u0 = get_fu( ep(:,i), bp(:,i), fu0 )
        tmp  = dot4( u0, f2u0 )

        u0(4) = fac_rr / u0(4)

        ! Energy of the emitted photon
        g_elec(i) = sqrt( 1.0_p_k_part + u0(1)*u0(1) + u0(2)*u0(2) + u0(3)*u0(3) )

        ! If the photon energy is lower than the threshold for its emission
        ! we keep the damping of the exact pusher
        if(g_elec(i) < this%qed_g_cutoff) then

          this%p(1,pp) = u0(1) + u0(4) * ( f2u0(1) - tmp * u0(1) )
          this%p(2,pp) = u0(2) + u0(4) * ( f2u0(2) - tmp * u0(2) )
          this%p(3,pp) = u0(3) + u0(4) * ( f2u0(3) - tmp * u0(3) )

        endif

        pp = pp + 1
      enddo
    endif

    ! -------------------------------------------------------------------------
    ! categorize particles
    ! -------------------------------------------------------------------------
    n1 = 0; ix1 = 0 
    n2 = 0; ix2 = 0
    do i = 1, np

      e2 = ep(1,i) * ep(1,i) + ep(2,i) * ep(2,i) + ep(3,i) * ep(3,i)
      b2 = bp(1,i) * bp(1,i) + bp(2,i) * bp(2,i) + bp(3,i) * bp(3,i)      

      if ( e2 + b2 < weak_fld_tol ) then

        k2(i) = 0.0_p_k_part; k(i) = 0.0_p_k_part
        o2(i) = 0.0_p_k_part; o(i) = 0.0_p_k_part

      else

        ! Lorentz invariants
        inv1 = e2 - b2
        inv2 = ep(1,i) * bp(1,i) + ep(2,i) * bp(2,i) + ep(3,i) * bp(3,i)

        ! calculate eigenvalues
        tmp = sqrt( inv1 * inv1 + 4.0_p_k_part * inv2 * inv2 )
        k2(i) = 0.5_p_k_part * ( inv1 + tmp ); k(i) = sqrt(k2(i))
        o2(i) = 0.5_p_k_part * (-inv1 + tmp ); o(i) = sqrt(o2(i))

        if ( k2(i) < rel_tol .and. o2(i) < rel_tol ) then
          ! light-like
          n1 = n1 + 1; ix1(n1) = i
        else
          ! general case
          n2 = n2 + 1; ix2(n2) = i
        endif

      endif
    enddo

    ! -------------------------------------------------------------------------
    ! advance particles in light-like field (category 1)
    ! -------------------------------------------------------------------------
    do j = 1, n1

      i = ix1(j); pp = ptrcur + i - 1

      u0(1:3) = this%p(:,pp)
      u0(4) = sqrt( 1.0_p_k_part + u0(1)*u0(1) + u0(2)*u0(2) + u0(3)*u0(3) )

      fu0  = get_fu( ep(:,i), bp(:,i), u0 )
      f2u0 = get_fu( ep(:,i), bp(:,i), fu0 )
      f3u0 = get_fu( ep(:,i), bp(:,i), f2u0 )

      ! solve for the proper time step according to dt
      call get_t_light( dt, u0(4), fu0(4), f2u0(4), f3u0(4), c1, c2, t, this%iter_tol )
      u = u0 + t * fu0 + c1 * f2u0 + c2 * f3u0

      this%p(:,pp) = u(1:3)

    enddo

    ! -------------------------------------------------------------------------
    ! advance particles in arbitrary field (category 2)
    ! -------------------------------------------------------------------------
    do j = 1, n2

      i = ix2(j); pp = ptrcur + i - 1

      u0(1:3) = this%p(:,pp)
      u0(4) = sqrt( 1.0_p_k_part + u0(1)*u0(1) + u0(2)*u0(2) + u0(3)*u0(3) )

      ! subspace decomposition
      k2o2 = k2(i) + o2(i)
      fu0  = get_fu( ep(:,i), bp(:,i), u0 )
      f2u0 = get_fu( ep(:,i), bp(:,i), fu0 )
      uk0  = ( f2u0 + o2(i) * u0 ) / k2o2
      uo0  = u0 - uk0

      fuk0 = get_fu( ep(:,i), bp(:,i), uk0 )
      fuo0 = fu0 - fuk0

      ! solve for the proper time step according to dt
      call get_t( dt, k(i), uk0(4), fuk0(4), o(i), uo0(4), fuo0(4), &
        ch_k, shc_k, chc_k, cn_o, snc_o, cnc_o, t, this%iter_tol )

      uk = uk0 * ch_k + fuk0 * shc_k * t
      uo = uo0 * cn_o + fuo0 * snc_o * t
      u = uk + uo

      this%p(1:3,pp) = u(1:3)

    enddo

    ! -------------------------------------------------------------------------
    ! half-dt push of radiation reaction
    ! -------------------------------------------------------------------------
    if ( this%rad_react ) then
      pp = ptrcur
      do i = 1, np
        u0(1:3) = this%p(1:3,pp)
        u0(4) = sqrt( 1.0_p_k_part + u0(1) * u0(1) + u0(2) * u0(2) + u0(3) * u0(3) )

        fu0  = get_fu( ep(:,i), bp(:,i), u0 )
        f2u0 = get_fu( ep(:,i), bp(:,i), fu0 )
        tmp  = dot4( u0, f2u0 )

        u0(4) = fac_rr / u0(4)

        ! If the photon energy is lower than the threshold for its emission
        ! we keep the damping of the exact pusher
        if(g_elec(i) < this%qed_g_cutoff) then
          this%p(1,pp) = this%p(1,pp) + u0(4) * ( f2u0(1) - tmp * u0(1) )
          this%p(2,pp) = this%p(2,pp) + u0(4) * ( f2u0(2) - tmp * u0(2) )
          this%p(3,pp) = this%p(3,pp) + u0(4) * ( f2u0(3) - tmp * u0(3) )

          ! Diagnostic of energy radiated classically
          if(this%ndump_fac_rad>0) then
            g_elec_tmp = sqrt( 1.0_p_k_part + this%p(1,pp) * this%p(1,pp) + this%p(2,pp) * this%p(2,pp) + this%p(3,pp) * this%p(3,pp) )
            energy_rad(1) = energy_rad(1) + (g_elec(i) - g_elec_tmp) * this%q(pp)
          endif

        endif

        pp = pp + 1
      enddo
    endif

    if(( this%rad_react ) .and. ( this%if_damp_qed ) ) call damp_qed_qedcyl( this, k_damp, ptrcur, np )

    pp = ptrcur
    do i = 1, np
      tmp = sqrt(1.0 + this%p(1,pp)*this%p(1,pp) + this%p(2,pp)*this%p(2,pp) + this%p(3,pp)*this%p(3,pp))
      loc_ene = loc_ene + this%q(pp) * (tmp - 1.0)
      pp = pp + 1
    enddo

    energy = energy + loc_ene

    if( ( this%ndump_fac_rad>0 ) .and. ( this%if_damp_classical ) ) this%energy_rad(1) = this%energy_rad(1) + energy_rad(1)

  enddo

  class default
      ! this must never happen
      call abort_program( p_err_invalid )
  end select

end subroutine dudt_exact_qedcyl
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine get_t_light( dt, u0, fu0, f2u0, f3u0, c1, c2, t, tol )
!-----------------------------------------------------------------------------------------
! Using Newton-Raphson method to calculate proper time step for the light-like case
!-----------------------------------------------------------------------------------------

  implicit none

  real(p_double), intent(in) :: dt, tol
  real(p_k_part), intent(in) :: u0, fu0, f2u0, f3u0
  real(p_k_part), intent(inout) :: t, c1, c2

  real(p_k_part) :: t1, t2, t3, fun, dfun
  integer :: iter

  t1 = dt / u0
  iter = 0

  do

    t2 = t1 * t1
    t3 = t2 * t1
    c1 = 0.5_p_double * t2
    c2 = 0.166666666666667_p_double * t3

    fun  = t1 * u0 + c1 * fu0 + c2 * f2u0 + 0.041666666666667_p_double * t3 * t1 * f3u0 - dt
    dfun = u0 + t * fu0 + c1 * f2u0 + c2 * f3u0

    t = t1 - fun / dfun
    iter = iter + 1

    if ( abs(t - t1) / t1 < tol .or. iter > p_iter_max ) then
      exit
    else
      t1 = t
    endif

  enddo

end subroutine get_t_light
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine get_t( dt, k, uk0, fuk0, o, uo0, fuo0, ch_k, shc_k, chc_k, cn_o, snc_o, cnc_o, t, tol )
!-----------------------------------------------------------------------------------------
! Using Newton-Raphson method to calculate proper time step for general case
!-----------------------------------------------------------------------------------------

  implicit none

  real(p_double), intent(in) :: dt, tol
  real(p_k_part), intent(in) :: k, uk0, fuk0
  real(p_k_part), intent(in) :: o, uo0, fuo0
  real(p_k_part), intent(inout) :: ch_k, shc_k, chc_k, cn_o, snc_o, cnc_o, t

  real(p_k_part) :: fun, dfun, t0
  integer :: iter

  t0 = dt / ( uk0 + uo0 )
  iter = 0

  do

    call  trig_func( o * t0, cn_o, snc_o, cnc_o )
    call hyper_func( k * t0, ch_k, shc_k, chc_k )

    fun = ( uk0 * shc_k + fuk0 * chc_k * t0 ) * t0 + &
          ( uo0 * snc_o - fuo0 * cnc_o * t0 ) * t0 - dt

    dfun = uk0 * ch_k + fuk0 * shc_k * t0 + &
           uo0 * cn_o + fuo0 * snc_o * t0

    t = t0 - fun / dfun
    iter = iter + 1

    if ( abs(t - t0) / t0 < tol .or. iter > p_iter_max ) then
      exit
    else
      t0 = t
    endif

  enddo

end subroutine get_t
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine trig_func( x, cn, snc, cnc, sncc )
!-----------------------------------------------------------------------------------------
! Fast evaluation of the trigonometric functions cos(x), sinc(x), (cos(x)-1)/x^2 and 
! (sinc(x)-1)/x^2. When x is smaller than a threshold value, they are evaluated using
! Taylor expansion, otherwise using tan(x/2) for evaluation.
!-----------------------------------------------------------------------------------------
  implicit none
  real(p_k_part), intent(in) :: x
  real(p_k_part), intent(out) :: cn, snc, cnc
  real(p_k_part), intent(out), optional :: sncc

  real(p_double) :: xd, x2, x4, tmp1, tmp2, tmp3, cn_d, snc_d
  ! threshold making Taylor expansion accurate to machine precision
  real(p_double), save :: eps = 0.007367081148666_p_double

  xd = real(x,p_double)
  x2 = xd * xd

  if ( abs(xd) < eps ) then
    x4 = x2 * x2
    cn  = real( 1.0_p_double - 0.500000000000000_p_double * x2 + 0.041666666666667_p_double * x4, p_k_part )
    snc = real( 1.0_p_double - 0.166666666666667_p_double * x2 + 0.008333333333333_p_double * x4, p_k_part )
    cnc =real( -0.5_p_double + 0.041666666666667_p_double * x2 - 0.001388888888889_p_double * x4, p_k_part )
    if ( present(sncc) ) then
      sncc = real( -0.166666666666667_p_double + 0.008333333333333_p_double * x2 - 1.984126984126984d-4 * x4, p_k_part )
    endif
  else
    tmp1 = tan( 0.5_p_double * xd )
    tmp2 = tmp1 * tmp1
    tmp3 = 1.0_p_double / (1.0_p_double + tmp2)
    cn_d = ( 1.0_p_double - tmp2 ) * tmp3
    cn = real( cn_d, p_k_part )
    snc_d = 2.0_p_double * tmp1 * tmp3 / xd
    snc = real( snc_d, p_k_part )
    cnc = real( ( cn_d - 1.0_p_double ) / x2, p_k_part )
    if ( present(sncc) ) then
      sncc = real( ( snc_d - 1.0_p_double ) / x2, p_k_part )
    endif
  endif

end subroutine trig_func

!-----------------------------------------------------------------------------------------
subroutine hyper_func( x, ch, shc, chc, shcc )
!-----------------------------------------------------------------------------------------
! Fast evaluation of the hyperbolic functions cosh(x), sinhc(x), (cosh(x)-1)/x^2 and 
! (sinhc(x)-1)/x^2. When x is smaller than a threshold value, they are evaluated using
! Taylor expansion, otherwise using tanh(x/2) for evaluation.
!-----------------------------------------------------------------------------------------
  implicit none
  real(p_k_part), intent(in) :: x
  real(p_k_part), intent(out) :: ch, shc, chc
  real(p_k_part), intent(out), optional :: shcc

  real(p_double) :: xd, x2, x4, tmp1, tmp2, tmp3, ch_d, shc_d
  ! threshold making Taylor expansion accurate to machine precision
  real(p_double), save :: eps = 0.007367081148666_p_double

  xd = real(x,p_double)
  x2 = xd * xd

  if ( abs(xd) < eps ) then
    x4 = x2 * x2
    ch  = real( 1.0_p_double + 0.500000000000000_p_double * x2 + 0.041666666666667_p_double * x4, p_k_part )
    shc = real( 1.0_p_double + 0.166666666666667_p_double * x2 + 0.008333333333333_p_double * x4, p_k_part )
    chc = real( 0.5_p_double + 0.041666666666667_p_double * x2 + 0.001388888888889_p_double * x4, p_k_part )
    if ( present(shcc) ) then
      shcc = real( 0.166666666666667_p_double + 0.008333333333333_p_double * x2 + 1.984126984126984d-4 * x4, p_k_part )
    endif
  else
    tmp1 = tanh( 0.5_p_double * xd )
    tmp2 = tmp1 * tmp1
    tmp3 = 1.0_p_double / (1.0_p_double - tmp2)
    ch_d = ( 1.0_p_double + tmp2 ) * tmp3
    ch = real( ch_d, p_k_part )
    shc_d = 2.0_p_double * tmp1 * tmp3 / xd
    shc = real( shc_d, p_k_part )
    chc = real( ( ch_d - 1.0_p_double ) / x2, p_k_part )
    if ( present(shcc) ) then
      shcc = real( ( shc_d - 1.0_p_double ) / x2, p_k_part )
    endif
  endif

end subroutine hyper_func

!-----------------------------------------------------------------------------------------
function dot4( u, v )

  implicit none
  real(p_k_part), intent(in), dimension(4) :: u, v
  real(p_k_part) :: dot4

  dot4 = u(4) * v(4) - u(1) * v(1) - u(2) * v(2) - u(3) * v(3)

end function dot4

!-----------------------------------------------------------------------------------------
function get_fu( ep, bp, u )
!-----------------------------------------------------------------------------------------
! Calculate the product of field tensor F and a four-vector u
!-----------------------------------------------------------------------------------------
  implicit none
  real(p_k_part), intent(in), dimension(3) :: ep, bp
  real(p_k_part), intent(in), dimension(4) :: u
  real(p_k_part), dimension(4) :: get_fu

  get_fu(4) =                ep(1) * u(1) + ep(2) * u(2) + ep(3) * u(3)
  get_fu(1) = ep(1) * u(4) +                bp(3) * u(2) - bp(2) * u(3)
  get_fu(2) = ep(2) * u(4) - bp(3) * u(1) +                bp(1) * u(3)
  get_fu(3) = ep(3) * u(4) + bp(2) * u(1) - bp(1) * u(2)
end function get_fu

subroutine bremsstrahlung_emission(this, dt, k_damp, ptrcur, i0, i1, np)
  use m_parameters

  implicit none

  ! input parameters
  class( t_species_qedcyl ),    intent(inout) :: this
  real(p_double), intent(in) :: dt
  real(p_double), dimension(p_cache_size), intent(in) ::  k_damp

  integer, intent(in) :: ptrcur, np, i0, i1

  ! dummy variables
  integer :: i, pp
  real(p_double) :: proba, varand, ge, gam_k, g_p, erad
  real(p_k_part) :: w_p
  real(p_k_part), dimension(p_p_dim) :: g_p_vec

  ! Bremsstrahlung variables
  real(p_double) :: factor_mult, omega_p, n_p, ni, momentum_damped, normp

  ! Diagnostics
  real(p_k_part), dimension(2) :: energy_rad

  integer :: Z_ion, cs_model
  logical :: energy_damp

  ! Bremsstrahlung variables of interest
  factor_mult = this % bremsstrahlung_info % proba_mult
  omega_p = this%omega_p0
  n_p = this%n0
  Z_ion = int(this % bremsstrahlung_info % Z_ion)
  cs_model = int(this % bremsstrahlung_info % cs_model)
  energy_damp = this % bremsstrahlung_info % energy_damp
  energy_rad = 0.0_p_k_part

  ! loop similar to os-photons
  pp = ptrcur
  do i=1, np

      ge = sqrt(1.0d0 + this%p(1,pp)**2+this%p(2,pp)**2+this%p(3,pp)**2)
      normp = sqrt(this%p(1,pp)**2+this%p(2,pp)**2+this%p(3,pp)**2)

      if ((ge-1.) .gt. 1.0d-3) then

        ! Ion density at the position of the photon
        ni = this % bremsstrahlung_info % rho_ion % f2(1, this%ix(1, pp), this%ix(2, pp))

        call br_proba(Z_ion, ge, ni, dt, omega_p, n_p, factor_mult, proba)

        call random_number( varand )

        if (varand < proba) then

            call br_spec(Z_ion, ge, gam_k)

            ! photon momentum
            g_p_vec = gam_k*this%p(:,pp)/sqrt(this%p(1,pp)**2+this%p(2,pp)**2+this%p(3,pp)**2)

            ! photon weight
            w_p = abs(this%q(pp)) / factor_mult

            ! Add the photon
            if (gam_k > this%p_emit_cutoff) then
              call this%photons%create_particle_single_cell_p_cyl( &
               this%ix(1:p_x_dim,pp), this%x(:,pp), g_p_vec, w_p )
              energy_rad(2) = energy_rad(2) + gam_k * w_p
            else
              energy_rad(1) = energy_rad(1) + gam_k * w_p
            endif

            ! Damp the momentum of the electron, via energy balance
            if (energy_damp) then
              momentum_damped = sqrt((ge - gam_k / factor_mult)**2 - 1.0_p_k_part)
              this%p(:,pp) = momentum_damped * this%p(:,pp) / normp
            end if

        endif

      endif
      pp = pp + 1
  enddo

  ! accumulate radiated energy
  this % energy_rad(4) = this % energy_rad(4) + energy_rad(1)
  this % energy_rad(5) = this % energy_rad(5) + energy_rad(2)

end subroutine bremsstrahlung_emission

subroutine trident_coul_pair_creation( this, dt, ptrcur, i0, i1, np)

      use m_parameters

      implicit none

      ! input parameters
      class( t_species_qedcyl ),    intent(inout) :: this
      real(p_double), intent(in) :: dt
      integer, intent(in) :: ptrcur, i0, i1, np

      ! dummy variables
      integer :: i, pp
      real(p_double) :: proba, varand, gam_1, gam_pa, gam_p, gam_e

      ! Trident_coul variables
      real(p_double) :: factor_mult, omega_p, n_p, ni, normp, momentum_damped
      real(p_k_part) :: q_pair
      integer :: Z_ion, cs_model
      logical :: energy_damp

      ! attributes for new particles
      real(p_k_part)   , dimension(p_p_dim) :: p_new

      ! Trident_coul variables of interest
      factor_mult = this % trident_coul_info % proba_mult
      omega_p = this%omega_p0
      n_p = this%n0
      Z_ion = int(this % trident_coul_info % Z_ion)
      cs_model = int(this % trident_coul_info % cs_model)
      energy_damp = this % trident_coul_info % energy_damp

      pp = ptrcur
      do i=1, np

         gam_1 = sqrt(1.0d0 + this%p(1,pp)**2+this%p(2,pp)**2+this%p(3,pp)**2)
         normp = sqrt(this%p(1,pp)**2+this%p(2,pp)**2+this%p(3,pp)**2)

         ! gam_1 > 6 because the table starts for this value
         if (gam_1 .gt. 6.0d0) then

            ! Ion density at the position of the electron/positron
            ni = this % trident_coul_info % rho_ion % f2(1, this%ix(1, pp), this%ix(2, pp))

            ! Sample probability of an event
            call ct_proba(Z_ion, gam_1, ni, dt, omega_p, n_p, factor_mult, proba)
            call random_number( varand )
            if (varand < proba) then

               ! Sample the total energy of the pair
               call ct_spec_pair(Z_ion, gam_1, gam_pa)

               ! Sample the energy of the positron
               call ct_spec_posi(Z_ion, gam_1, gam_pa, gam_p)

               ! For test purpose only
               if (cs_model .eq. 0.) gam_p = 0.5_p_double * gam_pa

               ! Deduce the energy of the 2ndary electron
               gam_e = gam_pa - gam_p

               ! Do not put absolute value
               q_pair = this%q(pp) / factor_mult

               ! Diagnostics
               if ( (gam_e.gt.1.0d0) .and. (gam_p.gt.1.0d0) ) then

                  ! Only for test purposes
                  if (cs_model .eq. 0.) then
                     if (this % rqm .eq. -1.0_p_double) then
                        p_new  = real( sqrt((gam_p)**2-1)*this%p(:,pp)/normp, p_k_part )
                        call this%recip_spec%create_particle_single_cell_p_cyl( this%ix(:,pp), this%x(:,pp), p_new, -q_pair )
                     endif
                  else
                     ! create a new particle of same species
                     p_new  = real( sqrt((gam_e)**2-1)*this%p(:,pp)/normp, p_k_part )
                     call this%create_particle_single_cell_p_cyl( this%ix(:,pp), this%x(:,pp), p_new, q_pair )

                     ! create a new particle of the complementary species
                     p_new = real( sqrt((gam_p)**2-1)*this%p(:,pp)/normp, p_k_part )
                     call this%recip_spec%create_particle_single_cell_p_cyl( this%ix(:,pp), this%x(:,pp), p_new, -1.0_p_k_part*q_pair )
                  endif

                  this%num_pairs_ct = this%num_pairs_ct + 1
                  this%weight_pairs_ct = this%weight_pairs_ct + abs(q_pair)
                  this%energy_pairs_ct = this%energy_pairs_ct + abs(q_pair)*(gam_e + gam_p)

                  ! Damp the momentum of the electron, via energy balance
                  if (energy_damp) then
                     momentum_damped = sqrt((gam_1 - gam_pa / factor_mult)**2 - 1.0_p_k_part)
                     this%p(:,pp) = momentum_damped * this%p(:,pp) / normp
                  end if

               endif

            endif

         endif
         pp = pp + 1
      enddo

   end subroutine trident_coul_pair_creation

end module m_species_push_qedcyl

!-----------------------------------------------------------------------------------------
! Set the correct dudt function for qed
!-----------------------------------------------------------------------------------------
subroutine set_dudt_spec_qedcyl( this  )

  use m_parameters
  use m_species_define_qedcyl, only : t_species_qedcyl
  use m_species_push_qedcyl, only : dudt_boris_qedcyl, dudt_exact_qedcyl

    implicit none

    class( t_species_qedcyl ), intent(inout) :: this

  ! Set pointer to dudt function
  select case ( this % push_type )
  case(p_std)
    this % dudt => dudt_boris_qedcyl
  case(p_exact)
    this % dudt => dudt_exact_qedcyl
  case default
    this % push_type = p_std
  end select

end subroutine set_dudt_spec_qedcyl
!-----------------------------------------------------------------------------------------
