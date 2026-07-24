!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     cathode class
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module m_cathode_cyl_modes

#include "os-config.h"
#include "os-preprocess.fpp"
#include "memory/memory.h"

use m_parameters
use m_fparser
use m_restart
use m_node_conf
use m_vdf_define
use m_vdf_math
use m_math, only : pi
use m_cyl_modes
use m_species_define, only : t_species, p_cell_near, p_uniform
use m_species_cyl_modes_define
use m_species_current
use m_species_current_cyl_modes
use m_species_udist
use m_cathode

implicit none

private

! Should be the same as in os-cathode.f03
! temporal density profile types
integer, parameter :: p_polynomial = 1 ! Uses t_rise, t_fall and t_flat
integer, parameter :: p_func       = 2 ! Math function

! transverse density profile types
! integer, parameter :: p_uniform = 0  ! uniform, already defined in m_species_define
integer, parameter :: p_gaussian  = 2  ! gaussian
integer, parameter :: p_channel   = 3  ! parabolic channel
integer, parameter :: p_step      = 10 ! Constant within region then zero outside
integer, parameter :: p_thruster  = 11 ! Custom

type, extends(t_cathode) :: t_cathode_cyl_modes

  contains

  ! note that the parent procedures cannot be overwritten
  procedure :: inject_cyl              => inject_cathode_cyl_modes
  procedure :: inject_cyl_lowroundoff  => inject_cathode_cyl_modes_lowroundoff

  procedure :: init                    => init_cathode_cyl_modes

end type t_cathode_cyl_modes

public :: t_cathode_cyl_modes

contains

!---------------------------------------------------

!-----------------------------------------------------------------------------------------
! Initialize object
!-----------------------------------------------------------------------------------------
subroutine init_cathode_cyl_modes(this, cathode_id, restart, restart_handle, t, dt, &
                                  coordinates)
!---------------------------------------------------
!       sets up this data structure from the given information
!---------------------------------------------------

  implicit none

  ! dummy variables
  class(t_cathode_cyl_modes), intent(inout) :: this
  integer, intent(in) :: cathode_id
  logical, intent(in) :: restart
  type(t_restart_handle), intent(in) :: restart_handle
  real(p_double), intent(in) :: t, dt
  integer, intent(in) :: coordinates

  ! local variables
  integer, dimension(p_x_dim) :: lnum_par_x ! number of particles per cell

  ! Variables required for particle initialization inside simulation box
  real(p_double)  :: wall_disp
  integer         :: iwall_disp
  real(p_double), dimension(p_p_dim) :: u ! fluid generalized momenta
  real(p_k_part)  :: t_past

  integer :: i
  logical :: inject_node

  ! executable statements

  ! get particle spacing
  lnum_par_x(1:p_x_dim) = this%species%num_par_x(1:p_x_dim)
  call alloc(this%ppos_cell, (/ 1, lnum_par_x(p_r_dim) /))
  do i = 1, lnum_par_x(p_r_dim)
    this%ppos_cell(1, i) = (2 * i - 1 - lnum_par_x(p_r_dim)) * 0.5_p_k_part / lnum_par_x(p_r_dim)
  enddo


  if (restart) then

    call this%restart_read(restart_handle)

  else

    ! set the internal cathode id
    this%cathode_id = cathode_id

    ! check if associated species is ok
    if (.not. associated(this%species)) then
      ERROR('species is not associated in setup_cathode')
      call abort_program(p_err_invalid)
    endif

    ! get the species id for restart information
    this%sp_id = this%species%sp_id

    select case(this%wall)
    case (p_lower)
      this%inppos = 0
      this%nppos  = 0.5_p_double - 0.5_p_double / lnum_par_x(this%dir)
    case (p_upper)
      this%inppos = this%species%g_nx(this%dir) + 1
      this%nppos  = 0.5_p_double / lnum_par_x(this%dir) - 0.5_p_double
    end select

    ! get fluid velocity from particle momenta for initial particle injection
    u = this%species%udist%ufl
    this%rgamma = 1.0_p_double/sqrt(1.0_p_double + u(1)**2 + u(2)**2 + u(3)**2)
    this%vel = u * this%rgamma / real(this%species%dx(this%dir), p_double)

    ! correct particle injection position with t_start (that can be positive or negative)
    ! Past logic assumed that tmin of the simulation was 0 (most often is).
    ! We have to check the difference between t_start and t of the simulation
    ! Also account for the zero time advance by subtracting dt
    t_past = real(this%t_start - t - dt, p_k_part)

    if (t_past /= 0.0_p_double) then

      wall_disp   = this%vel(this%dir) * t_past
      iwall_disp  = floor(wall_disp)
      wall_disp   = wall_disp - iwall_disp ! wall_disp will be in ]-1,1[

      this%inppos = this%inppos - iwall_disp
      this%nppos  = this%nppos  - wall_disp

      call ntrim_pos(this%nppos, this%inppos)

    endif

    ! check if initial particle injection
    if (t_past < 0.0_p_double) then

      ! adjust multi-node injection
      select case(this%wall)
      case (p_lower)
        inject_node = (this%inppos >= this%species%my_nx_p(p_lower, this%dir))

      case (p_upper)
        inject_node = (this%inppos <= this%species%my_nx_p(p_upper, this%dir))

      case default
        inject_node = .false.
      end select

      ! inject only if necessary
      if (inject_node) then
        ! correct for zero time advance and inject particles
        this%nppos = this%nppos - this%vel(this%dir) * real(dt, p_double)

        call ntrim_pos(this%nppos, this%inppos)

        ! no need to deposit current
        call create_particles_cathode_cyl_modes(this, t, dt, coordinates)
      endif

    endif

    this%t_total = this%t_start + this%t_rise + this%t_flat + this%t_fall

  endif

  if (this%deposit_current) call this % species % bnd_con % init_tmp_buf_current()

end subroutine init_cathode_cyl_modes
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Deposit current from injected cathode particles
!-----------------------------------------------------------------------------------------
subroutine deposit_current_cathode_cyl_modes(this, ip0, ip1, jay_cyl_m, gdt)

  implicit none

  integer, parameter :: p_cm_x = 3, p_cm_y = 4

  class(t_cathode_cyl_modes), intent(inout) :: this  ! cathode object
  integer, intent(in) :: ip0, ip1  ! first and last particle injected
  type(t_cyl_modes), intent(inout) :: jay_cyl_m  ! electrical current object
  real(p_double), intent(in) :: gdt         ! time step

  ! local variables
  integer :: i, k, pp, ip, npar
  integer :: gix2, shift_ix2, mode
  real(p_k_part), dimension(p_x_dim) :: dt_dx
  real(p_double) :: dr, rdr, r_shift
  real(p_double) :: sin_th, cos_th, x_old, y_old, r_old, x_half, y_half, r_half
  class(t_species_cyl_modes), pointer :: spec_cyl

  ! Get a pointer of the correct kind to the species object
  select type(obj => this%species)
  class is (t_species_cyl_modes)
    spec_cyl => obj
  class default
    ERROR('inject_cyl_modes was called with a species of the wrong class')
    call abort_program(p_err_invalid)
  end select

  shift_ix2 = spec_cyl%my_nx_p(p_lower, 2) - 2
  dr = spec_cyl%dx(2)
  rdr = 1.0_p_double / dr

  ! Since r = 0 is at the center of cell ix2 = 1 the radial position of particles will be
  ! - Positions defined with regard to the center of the cell (odd interpolation)
  !       r = ((gix2-1) + x2) * dr
  ! - Positions defined with regard to the corner of the cell (even interpolation)
  !       r = ((gix2-1) + x2 - 0.5) * dr
  if (spec_cyl%pos_type == p_cell_near) then
    r_shift = 0.5_p_double
  else
    r_shift = 0.0_p_double
  endif

  do i = 1, p_x_dim
    dt_dx(i) = gdt / spec_cyl%dx(i)
  enddo

  do ip = ip0, ip1, p_cache_size

    ! check if last copy of table and set np
    if(ip + p_cache_size > ip1) then
      npar = ip1 - ip + 1
    else
      npar = p_cache_size
    endif

    pp = ip
    ! Get corrected particle momenta, corrected charge, and 1 / gamma
    do k = 1, npar
      tmp_xnew(1:p_x_dim + 2, k) = spec_cyl%x(1:p_x_dim + 2, pp)
      tmp_p(1:p_p_dim, k) = spec_cyl%p(1:p_p_dim, pp)
      tmp_q(k) = spec_cyl%q(pp)
      tmp_rg(k) = 1.0_p_k_part / sqrt(1.0_p_k_part + tmp_p(1, k)**2 + tmp_p(2, k)**2 + tmp_p(3, k)**2)

      pp = pp + 1
    enddo

    ! Get previous position and deposit current
    pp = ip
    do k = 1, npar
      tmp_xold(1, k)  = tmp_xnew(1, k) - tmp_p(1, k) * tmp_rg(k) * dt_dx(1)

      ! Convert radial "cell" position to "box" position in double precision
      gix2 = spec_cyl%ix(2, pp) + shift_ix2

      ! Cartesian push in transverse plane
      ! Note: these coordinates are of "box" type
      x_old = tmp_xnew(p_cm_x, k) - tmp_p(2, k) * tmp_rg(k) * gdt
      y_old = tmp_xnew(p_cm_y, k) - tmp_p(3, k) * tmp_rg(k) * gdt

      ! Store cartesian coordinates (used for cyl_modes current deposition)
      tmp_xold(p_cm_x, k) = x_old
      tmp_xold(p_cm_y, k) = y_old
      r_old = sqrt(x_old**2 + y_old**2)

      ! this is a protection against roundoff for cold plasmas
      if (tmp_xnew(p_cm_x, k) == x_old .and. tmp_xnew(p_cm_y, k) == y_old) then
        tmp_xold(2, k) = tmp_xnew(2, k)
      else
        tmp_xold(2, k) = real(r_old * rdr - (gix2 - r_shift), p_k_part)
      endif

      ! Get number of cells moved in each (z,r) direction
      tmp_dxi(1, k) = ntrim(tmp_xold(1, k)) ! returns +1 , 0, or -1
      tmp_dxi(2, k) = ntrim(tmp_xold(2, k))

      ! Half cartesian push in transverse plane (to be centered with momentum)
      x_half = tmp_xnew(p_cm_x, k) - 0.5_p_k_part * tmp_p(2, k) * tmp_rg(k) * gdt
      y_half = tmp_xnew(p_cm_y, k) - 0.5_p_k_part * tmp_p(3, k) * tmp_rg(k) * gdt
      r_half = sqrt(x_half**2 + y_half**2)

      ! convert momentum into cylindrical coordinates before depositing current
      cos_th = x_half / r_half
      sin_th = y_half / r_half

      ! note that tmp_p here is in cylindrical coordinates
      tmp_p(1,i) =  spec_cyl%p(1, pp)
      tmp_p(2,i) =  spec_cyl%p(2, pp) * cos_th + spec_cyl%p(3,pp) * sin_th
      tmp_p(3,i) = -spec_cyl%p(2, pp) * sin_th + spec_cyl%p(3,pp) * cos_th

      pp = pp + 1
    enddo

    ! deposit current
    ! handle mode 0 first - for this part the original 2D deposit should work fine
    call deposit_current_2d(spec_cyl, jay_cyl_m%pf_re(0), tmp_dxi, tmp_xold(1:2, :), &
            spec_cyl%ix(:, ip:), tmp_xnew(1:2, :), tmp_q, tmp_rg, tmp_p, npar, gdt)

    ! handle high order modes - should use special cylindrical mode deposit to remain charge conserving
    do mode = 1, spec_cyl%n_cyl_modes

      ! still need to port cylindrical mode current deposition functions
      select case (spec_cyl%interpolation)
      case(p_linear)

        call getjr_cyl_m_s1(jay_cyl_m%pf_re(mode), jay_cyl_m%pf_im(mode), tmp_dxi, &
                            tmp_xold, spec_cyl%ix(:, ip:), tmp_xnew, &
                            tmp_q, tmp_rg, tmp_p, &
                            npar, gdt, shift_ix2, mode)

      case(p_quadratic)

        call getjr_cyl_m_s2(jay_cyl_m%pf_re(mode), jay_cyl_m%pf_im(mode), tmp_dxi, &
                            tmp_xold, spec_cyl%ix(:, ip:), tmp_xnew, &
                            tmp_q, tmp_rg, tmp_p, &
                            npar, gdt, shift_ix2, mode)

      case(p_cubic)

        call getjr_cyl_m_s3(jay_cyl_m%pf_re(mode), jay_cyl_m%pf_im(mode), tmp_dxi, &
                            tmp_xold, spec_cyl%ix(:, ip:), tmp_xnew, &
                            tmp_q, tmp_rg, tmp_p, &
                            npar, gdt, shift_ix2, mode)

      case default
        ERROR('Not implemented yet')
        call abort_program(p_err_notimplemented)
      end select

    enddo ! mode

  enddo

end subroutine deposit_current_cathode_cyl_modes
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Inject particles from cathode
!-----------------------------------------------------------------------------------------
subroutine inject_cathode_cyl_modes(this, jay_cyl_m, t, gdt, no_co, coords)

  implicit none

  class(t_cathode_cyl_modes), intent(inout) :: this ! cathode object
  type(t_cyl_modes), intent(inout) :: jay_cyl_m ! electrical current object
  real(p_double), intent(in) :: t    ! simulation time
  real(p_double), intent(in) :: gdt  ! time step
  class(t_node_conf), intent(in) :: no_co  ! node configuration
  integer, intent(in) :: coords ! coordinate system being used

  ! local variables
  ! number of particles before injection
  integer :: np

  ! parallel run variables
  integer :: inj_node

  ! check if correct node for launching particles
  select case (this%wall)
  case (1)
    inj_node = 1
  case (2)
    inj_node = nx(no_co, this%dir)
  case default
    inj_node = -1
  end select

  if (my_ngp(no_co, this%dir) == inj_node) then

    ! check if injection is finished
    if (t >= this%t_start .and. (t <= this%t_total .or. this%temp_type == p_func)) then

      ! get number of particles currently in the buffer
      np = this%species%num_par

      call create_particles_cathode_cyl_modes(this, t, gdt, coords)

      ! deposit current corresponding to the injection
      if (this%species%num_par > np .and. this%deposit_current) then
        call deposit_current_cathode_cyl_modes(this, np + 1, this%species%num_par, jay_cyl_m, gdt)
      endif

    endif ! (t <= this%t_total)

  endif ! injecting node


  ! call validate(this%species, "after cathode injection")

end subroutine inject_cathode_cyl_modes
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Inject particles from cathode, using low roundoff current deposition
!-----------------------------------------------------------------------------------------
subroutine inject_cathode_cyl_modes_lowroundoff(this, jay_cyl_m, jay_tmp_cyl_m, t, gdt, &
                                                no_co, coords)

  implicit none

  class(t_cathode_cyl_modes), intent(inout) :: this        ! cathode object
  type(t_cyl_modes), intent(inout) :: jay_cyl_m, jay_tmp_cyl_m  ! electrical current object
  real(p_double), intent(in)    :: t            ! simulation time
  real(p_double), intent(in)    :: gdt          ! time step
  class(t_node_conf), intent(in)    :: no_co       ! node configuration
  integer, intent(in) :: coords ! coordinate system being used

  ! local variables
  ! number of particles before injection
  integer :: np

  ! parallel run variables
  integer :: inj_node

  integer, dimension(2, 3) :: range
  integer :: mode
  class(t_species_cyl_modes), pointer :: spec_cyl

  ! Get a pointer of the correct kind to the species object
  select type(obj => this%species)
  class is (t_species_cyl_modes)
    spec_cyl => obj
  class default
    ERROR('inject_cyl_modes was called with a species of the wrong class')
    call abort_program(p_err_invalid)
  end select

  ! check if correct node for launching particles
  select case (this%wall)
  case (1)
    inj_node = 1
  case (2)
    inj_node = nx(no_co, this%dir)
  case default
    inj_node = -1
  end select

  if (my_ngp(no_co, this%dir) == inj_node) then

    ! check if injection is finished
    if (t >= this%t_start .and. (t <= this%t_total .or. this%temp_type == p_func)) then

      ! get number of particles currently in the buffer
      np = spec_cyl%num_par

      call create_particles_cathode_cyl_modes(this, t, gdt, coords)

      ! deposit current corresponding to the injection
      if ((spec_cyl%num_par > np) .and. (this%deposit_current)) then

        ! optimized implementation - only zero / add relevant cells
        select case (this%wall)
        case (1)
          range(1, this%dir) = 1 - spec_cyl%interpolation
          range(2, this%dir) = 1 + spec_cyl%interpolation
        case (2)
          range(1, this%dir) = jay_cyl_m%pf_re(0)%nx_(this%dir) + 1 - spec_cyl%interpolation
          range(2, this%dir) = jay_cyl_m%pf_re(0)%nx_(this%dir) + 1 + spec_cyl%interpolation
        end select
        range(1, p_r_dim) = jay_cyl_m%pf_re(0)%lbound(p_r_dim)
        range(2, p_r_dim) = jay_cyl_m%pf_re(0)%ubound(p_r_dim)

        ! zero the temporary grid
        jay_tmp_cyl_m%pf_re(0) = 0.0_p_k_fld
        do mode = 1, spec_cyl%n_cyl_modes
          jay_tmp_cyl_m%pf_re(mode) = 0.0_p_k_fld
          jay_tmp_cyl_m%pf_im(mode) = 0.0_p_k_fld
        enddo

        ! Deposit current on temporary grid
        call deposit_current_cathode_cyl_modes(this, np + 1, spec_cyl%num_par, jay_tmp_cyl_m, gdt)

        ! add the values to electric current grid in the relevant range
        call add(jay_cyl_m%pf_re(0), jay_tmp_cyl_m%pf_re(0), range)
        do mode = 1, spec_cyl%n_cyl_modes
          call add(jay_cyl_m%pf_re(mode), jay_tmp_cyl_m%pf_re(mode), range)
          call add(jay_cyl_m%pf_im(mode), jay_tmp_cyl_m%pf_im(mode), range)
        enddo

      endif

    endif ! (t <= this%t_total)

  endif ! injecting node

end subroutine inject_cathode_cyl_modes_lowroundoff
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Advance cathode particle position, put cathode particles that are inside simulation space
! in the species buffer and set their momenta
!-----------------------------------------------------------------------------------------
subroutine create_particles_cathode_cyl_modes(this, t, dt, coords)

  implicit none

  ! dummy variables
  class(t_cathode_cyl_modes), intent(inout) :: this  ! cathode object
  real(p_double), intent(in) :: t   ! simulation time
  real(p_double), intent(in) :: dt  ! time step
  integer, intent(in) :: coords  ! coordinate system being used

  ! local variables
  real(p_k_part) :: dpx        ! particle spacing along the injection direction
  real(p_k_part) :: lden_min   ! den_min for the species to inject
  integer        :: np         ! number of particles to inject along the injection direction

  integer :: i, npar, tr_parser
  real(p_k_part)                       :: qnorm     ! normalization factor for injected charge
  real(p_k_part), dimension(p_x_dim + 2) :: xnewpart  ! position,
  integer, dimension(p_x_dim)        :: ixnewpart ! momenta and
  real(p_k_part)                       :: qnewpart  ! charge of the new particle

  real(p_k_part)                         :: den
  real(p_k_fparse), dimension(p_max_dim) :: eval_var

  integer :: ptrcur

  real(p_k_part), dimension(p_x_dim)  :: gpos_min, dx
  real(p_k_part) :: inj_time, inj_pos, ppos_r, theta
  integer :: j, k, i2, ppcell_r, pptheta

  ! charge normalization coefficients for each particle
  real(p_k_part) :: alpha

  class(t_species_cyl_modes), pointer :: spec_cyl

  ! Get a pointer of the correct kind to the species object
  select type(obj => this%species)
  class is (t_species_cyl_modes)
    spec_cyl => obj
  class default
    ERROR('inject_cyl_modes was called with a species of the wrong class')
    call abort_program(p_err_invalid)
  end select

  lden_min = this%den_min
  dx(1:p_x_dim) = real(spec_cyl%dx(1:p_x_dim), p_k_part)
  gpos_min(1:p_x_dim) = spec_cyl%g_box(p_lower, 1:p_x_dim)

  ! advance the position of the next particle to inject
  this%nppos = this%nppos + this%vel(this%dir) * real(dt, p_k_part)
  call ntrim_pos(this%nppos, this%inppos)

  ! get number of particles in the buffer
  np = spec_cyl%num_par

  ! get number of particles per cell
  ppcell_r = spec_cyl%num_par_x(p_r_dim)
  pptheta = spec_cyl%num_par_theta

  ! get normalization factor for injected charge
  qnorm = sign(1.0_p_k_part / (spec_cyl%num_par_x(1) * spec_cyl%num_par_x(2)), spec_cyl%rqm)
  qnorm = qnorm / pptheta

  ! Get particle spacing along injection direction
  ! This is now done in cell normalized units
  dpx = 1.0 / spec_cyl%num_par_x(this%dir)

  do

    ! Get longitudinal beam density
    select case (this%wall)

    case (p_lower)

      if (this%inppos < spec_cyl%my_nx_p(p_lower, this%dir)) exit
      inj_time = real(t, p_k_part) - this%t_start - (this%inppos - this%nppos - 1) / this%vel(this%dir)

    case (p_upper)

      if (this%inppos > spec_cyl%my_nx_p(p_upper, this%dir)) exit
      inj_time = real(t, p_k_part) - this%t_start - &
                 (this%inppos - spec_cyl%my_nx_p(p_upper, this%dir) + this%nppos) / this%vel(this%dir)
    end select

    if (this%inppos >= spec_cyl%my_nx_p(p_lower, this%dir) .and. &
        this%inppos <= spec_cyl%my_nx_p(p_upper, this%dir)) then

      ! the total_xmoved term corrects for motion of injecting wall (moving window)
      ! after injection started
      inj_time = inj_time + real(spec_cyl%total_xmoved(this%dir), p_k_part)

      ! multiply by global density
      den = this%density * t_env(inj_time, this%temp_type, this%t_rise, this%t_flat, this%t_fall)

      ! inject particles
      xnewpart(this%dir) = this%nppos
      ixnewpart(this%dir) = this%inppos - spec_cyl%my_nx_p(p_lower, this%dir) + 1

      do i2 = spec_cyl%my_nx_p(p_lower, p_r_dim), spec_cyl%my_nx_p(p_upper, p_r_dim)

        ixnewpart(p_r_dim) = i2 - spec_cyl%my_nx_p(p_lower, p_r_dim) + 1
        do i = 1, ppcell_r

          ! get the radial position
          ppos_r = real(gpos_min(2) + (this%ppos_cell(1, i) + (i2 - 1)) * dx(p_r_dim), p_k_part)
          xnewpart(p_r_dim) = this%ppos_cell(1, i)

          ! get the density
          qnewpart = den * tr_profile(this, real(ppos_r - this%center(1), p_k_part))

          ! get the charge correction factor
          alpha = spec_cyl % norm_charge_axis(i2, ppos_r)
          do j = 1, pptheta

            theta = (j - 1) * (2.0_p_k_part * pi / pptheta)
            xnewpart(3) = ppos_r * cos(theta)
            xnewpart(4) = ppos_r * sin(theta)

            ! add particles
            if(abs(qnewpart) > lden_min .and. ppos_r > 0.0_p_k_part) then
              call spec_cyl%create_particle_single_cell_cyl(&
                ixnewpart, xnewpart, qnorm * qnewpart * ppos_r * alpha)
            endif

          enddo
        enddo

      enddo
    endif

    ! process next set of particles
    select case (this%wall)
    case (p_lower)
      this%nppos  = this%nppos - dpx

      ! Check injection halt for when using the moving injector plane
      if(this%mov_inj) then
        inj_pos = spec_cyl%g_box(p_lower, this%dir) + (this%inppos - 1) * dx(this%dir)
        if (inj_pos < int(-this%t_start - t)) exit
      endif

    case (p_upper)
      this%nppos = this%nppos + dpx

      ! Check injection halt for when using the moving injector plane
      if(this%mov_inj) then
        inj_pos = spec_cyl%g_box(p_lower, this%dir) + (this%inppos - 1) * dx(this%dir)
        if(inj_pos > int(spec_cyl%g_box(p_upper, this%dir) + this%t_start + t)) exit
      endif

    end select

    call ntrim_pos(this%nppos, this%inppos)

  enddo

  ! Set momentum of injected particles
  if (this%species%num_par > np) then

    ptrcur = np + 1
    npar = this%species%num_par - np

    ! set momentum of injected particles
    call set_momentum(this%species, ptrcur, ptrcur + npar - 1)

    ! process transversal momentum parser
    eval_var(1) = real(t, p_k_fparse)

    ! test if transverse parser present
    tr_parser = 0
    if(this%tr_u_is_parser(1)) tr_parser = 1
    if(this%tr_u_is_parser(2)) tr_parser = tr_parser + 2

    select case (tr_parser)

    case (1) ! transverse direction 1

      do k = ptrcur, ptrcur + npar - 1

        call spec_cyl%get_position(k, xnewpart)
        eval_var(2) = real(xnewpart(3), p_k_fparse)
        eval_var(3) = real(xnewpart(4), p_k_fparse)

        spec_cyl%p(2, k) = spec_cyl%p(2, k) + real(eval(this%tr_u_parser(1), eval_var), p_k_part)
      enddo

    case (2) ! transverse direction 2

      do k = ptrcur, ptrcur + npar - 1

        call spec_cyl%get_position(k, xnewpart)
        eval_var(2) = real(xnewpart(3), p_k_fparse)
        eval_var(3) = real(xnewpart(4), p_k_fparse)

        spec_cyl%p(3, k) = spec_cyl%p(3, k) + real(eval(this%tr_u_parser(1), eval_var), p_k_part)
      enddo

    case (3) ! transversal directions 1 & 2

      do k = ptrcur, ptrcur + npar - 1

        call spec_cyl%get_position(k, xnewpart)
        eval_var(2) = real(xnewpart(3), p_k_fparse)
        eval_var(3) = real(xnewpart(4), p_k_fparse)

        spec_cyl%p(2, k) = spec_cyl%p(2, k) + real(eval(this%tr_u_parser(1), eval_var), p_k_part)
        spec_cyl%p(3, k) = spec_cyl%p(3, k) + real(eval(this%tr_u_parser(1), eval_var), p_k_part)
      enddo

    case default

    end select
  endif

  ! ------------------------------------------

  contains

  ! function declaration of time envelope

  function t_env( t, temp_type, rise, flat, fall )
    implicit none
    real(p_k_part) :: t_env
    real(p_k_part), intent(in) :: t, rise, flat, fall
    integer, intent(in) :: temp_type

    select case( temp_type )
    case( p_polynomial )

      if (t < 0.0d0 .or. t > rise+flat+fall) then
        t_env = 0.0
      else if (t < rise) then
        t_env = poly_env(t/rise)
      else if ( t < rise+flat ) then
        t_env = 1.0
      else
        t_env = poly_env((rise+flat+fall-t)/fall)
      endif

    case( p_func )

      t_env = real( eval( this%temp_func, real( (/t/), p_k_fparse ) ), p_k_part )

    case default

      t_env = 0.0 ! Should never happen

    end select

  end function t_env


  ! transverse profile functions

  function tr_profile(cathode, r)
    implicit none
    real(p_k_part) :: tr_profile
    class(t_cathode_cyl_modes), intent(in) :: cathode
    real(p_k_part), intent(in) :: r

    select case (cathode%prof_type)
    case(p_uniform)
      tr_profile = 1.0_p_k_part
    case(p_step)
      if (abs(r) <= cathode%gauss_width/2.0) then
        tr_profile = 1.0_p_k_part
      else
        tr_profile = 0.0_p_k_part
      endif
    case(p_gaussian)
      if (abs(r) <= cathode%gauss_width/2.0) then
        !                  x = r
        tr_profile = exp(-2.0*(abs(r)/cathode%gauss_w0)**2)
      else
        tr_profile = 0.0_p_k_part
      endif
    case(p_thruster)
      if ((abs(r) <= cathode%gauss_width/2.0_p_k_part) .and. &
          (abs(r) > cathode%gauss_w0/2.0_p_k_part)) then
        tr_profile = 1.0_p_k_part
      else
        tr_profile = 0.0_p_k_part
      endif
    case (p_channel)

      if (abs(r) <= this%channel_size/2) then ! inside the channel
        tr_profile = this%channel_bottom + &
                     this%channel_depth * &
                     (abs(r)/this%channel_r0)**2
      else                                 ! outside the channel

        if (this%channel_wall <= 0.0) then ! finite channel
          tr_profile = this%channel_bottom + &
                       this%channel_depth * &
                       ((this%channel_size/2)/this%channel_r0)**2

        else                              ! leaky channel

          if (abs(r) < (this%channel_size/2 + this%channel_wall)) then
            tr_profile = (-abs(r) + this%channel_size/2 + this%channel_wall) / &
                         this%channel_wall * (this%channel_bottom + &
                         this%channel_depth * &
                         ((this%channel_size/2)/this%channel_r0)**2)
          else
            tr_profile = 0.0d0
          endif

        endif

      endif

    case(p_polynomial)
      if ( abs(r) <= cathode%gauss_width/2.0 ) then
        tr_profile = 1.0_p_k_part
      else if ( abs(r) < cathode%gauss_width/2.0 + cathode%gauss_w0 ) then
        tr_profile = poly_env(real((cathode%gauss_width/2.0 + cathode%gauss_w0 - abs(r)) &
                          / cathode%gauss_w0, p_k_part))
      else
        tr_profile = 0.0_p_k_part
      endif

    case default

      tr_profile = 1.0_p_k_part

    end select

  end function tr_profile


  ! gaussian like polynomial for time envelope

  function poly_env(x)
    implicit none
    real(p_k_part) :: poly_env
    real(p_k_part), intent(in) :: x

    poly_env = 10 * x**3 - 15 * x**4 + 6 * x**5
  end function poly_env

end subroutine create_particles_cathode_cyl_modes
!-----------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! Returns the integer shift (-1, 0 or +1) so that the coordinate remains in the [-0.5, 0.5[
! range. This is the fastest implementation (twice as fast as a sequence of ifs) because
! the two if structures compile as conditional moves and can be processed independently.
!---------------------------------------------------------------------------------------------------
function ntrim(x)
!---------------------------------------------------------------------------------------------------
  implicit none

  real(p_k_part), intent(in) :: x
  integer :: ntrim, a, b

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

  ntrim = a+b

end function ntrim
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine ntrim_pos( x, ix )

  implicit none

  real(p_k_part), intent(inout) :: x
  integer, intent(inout) :: ix

  integer :: dx

  dx = ntrim(x)
  x  = x - dx
  ix = ix + dx

end subroutine ntrim_pos
!---------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
end module m_cathode_cyl_modes
