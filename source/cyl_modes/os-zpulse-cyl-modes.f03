!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     new cyl modes pulse class
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "os-config.h"
#include "os-preprocess.fpp"

module m_zpulse_cyl_modes

#include "memory/memory.h"

  use m_node_conf,     only: t_node_conf
  use m_restart,       only: t_restart_handle, restart_io_write, restart_io_read
  use m_parameters
  use m_math
  use m_input_file,    only: t_input_file, get_namelist, p_max_nml_section_name
  use m_zpulse,        only: t_zpulse_list
  use m_zpulse_std,    only: p_zpulse_bnormal
  use m_zpulse_std,    only: p_max_chirp_order, p_polynomial, p_sin2, p_gaussian, p_func
  use m_zpulse_std,    only: p_const
  use m_zpulse_std,    only: p_hermite_gaussian, p_hermite_gaussian_astigmatic
  use m_zpulse_std,    only: p_laguerre_gaussian, p_plane, p_gaussian_asym, t_zpulse, p_forward, p_plane
  use m_zpulse_mov_wall, only : t_zpulse_mov_wall 
  use m_vdf_define,    only: t_vdf
  use m_vdf_math,      only: add
  use m_emf_define,    only: t_emf_bound, t_emf
  use m_emf_cyl_modes, only: t_emf_cyl_modes
  use m_space,         only: t_space, x_bnd, get_x_bnd
  use m_fparser,       only: p_max_expr_len, p_k_fparse, eval
  use m_grid_define,   only: t_grid
  use m_grid_cyl_modes,only: t_grid_cyl_modes
  use m_space,         only: xmin, xmax 
  use m_zpulse_mov_wall_cyl_modes, only : t_zpulse_mov_wall_cyl_modes
  use m_zpulse_flyfoc_cyl_modes, only : t_zpulse_flyfoc_cyl_modes
  use m_zpulse_planewaves_cyl_modes, only : t_zpulse_planewaves_cyl_modes

  implicit none

  private

  type, extends( t_zpulse_list ) :: t_zpulse_list_cyl_modes
    ! no new pulse specific parameters
  contains
    procedure :: read_input => read_input_zpulse_list_em_cyl_modes
  end type t_zpulse_list_cyl_modes


  type, extends( t_zpulse ) :: t_zpulse_cyl_modes
    ! no new pulse specific parameters
  contains
    procedure :: read_input => read_input_zpulse_cyl_modes
    procedure :: launch => launch_zpulse_cyl_modes
    procedure :: check_dimensionality => check_dimensionality_zpulse_cyl_modes
  end type t_zpulse_cyl_modes



  public :: t_zpulse_list_cyl_modes

  contains

subroutine read_input_zpulse_list_em_cyl_modes( this, input_file, g_space, bnd_con, periodic, grid, sim_options )

  implicit none

  class( t_zpulse_list_cyl_modes ), intent(inout) :: this
  class( t_input_file ), intent(inout) :: input_file
  type( t_space ), intent(in) :: g_space
  class (t_emf_bound), intent(in) :: bnd_con
  logical, dimension(:), intent(in) :: periodic
  class( t_grid ), intent(in) :: grid
  type( t_options ), intent(in) :: sim_options

  character(len = p_max_nml_section_name) :: section_name
  class( t_zpulse ), pointer :: tail, item

  tail => null()


  do
    section_name = trim(input_file % get_section_name() )
    select case (section_name )
    
    case('nl_zpulse')      
      allocate( t_zpulse_cyl_modes :: item )
    
    case('nl_zpulse_mov_wall' )
      allocate( t_zpulse_mov_wall_cyl_modes :: item )

    case('nl_zpulse_flyfoc' )
      allocate( t_zpulse_flyfoc_cyl_modes :: item )

    case('nl_zpulse_planewaves' )
      allocate( t_zpulse_planewaves_cyl_modes :: item )

    case default 
    
      if ( index( section_name, "zpulse") /= 0 ) then
        if ( mpi_node() == 0 ) then

          write(0,*)  ""
          write(0,*)  "   Error in zpulse parameters"
          write(0,*)  "   '" // trim(section_name) // "' is not yet supported in Quasi-3D simulations"
          write(0,*)  "   aborting..."
        endif
        stop
      endif

      ! otherwise not a zpulse section; exit. 
      exit 
    end select 

    select type( grid )
      class is( t_grid_cyl_modes )
        if ( grid%n_cyl_modes < 1 ) then
          if ( mpi_node() == 0 ) then
            write(0,*)  ""
            write(0,*)  "   Error in zpulse parameters"
            write(0,*)  "   'All zpulses require n_cyl_modes >= 1"
            write(0,*)  "   aborting..."
          endif
          stop
      endif
    end select

    ! Read input for zpulse
    call item % read_input( input_file, g_space, bnd_con, periodic, grid, sim_options )
    call item % check_dimensionality() 
    item % next => null()

    if (associated(tail)) then
      tail % next => item
    else
      this % list => item
    endif

    tail => item

  enddo

end subroutine read_input_zpulse_list_em_cyl_modes


!-----------------------------------------------------------------------------------------
! High order cylindrical modes start
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
function per_envelope_cylm( this, x1, x2 )

!-----------------------------------------------------------------------------------------
!       returns the value of the perpendicular
!       envelope at the requested position
!       this is identical to the 2D version in every way,
!       except for how the amplitude is calculated (which is analogous to 3D)
!-----------------------------------------------------------------------------------------

  implicit none

  class( t_zpulse_cyl_modes ), intent(in) :: this

  integer, parameter :: rank = 2

  ! x1 is the longitudinal coordinate, x2 the perpendicular one

  real(p_double), intent(in) :: x1, x2

  real(p_double) :: per_envelope_cylm
  real(p_double) :: z_center, r_center
  real(p_double) :: k, z, rWl2, Wl2, z0, rho, rho2

  real(p_double) :: gouy_shift, rWu, curv

  ! asymetric laser pulses variables
  real(p_double) :: w0_asym, focus_asym, shp

  integer :: i

  ! Boost variables
  real(p_double) :: uvel, beta, rg1pb

  ! get pulse center for chirp and phase calculations
  z_center = this % lon_center()
  r_center = this % per_center(1)

  ! get wavenumber
  k  = this%omega0

  ! add chirp
  if ( this%propagation == p_forward ) then
    do i = 1, this%chirp_order
      k = k + this%chirp_coefs(i) * ((x1 - z_center)**(i))
    enddo
    do i = 1, this%per_chirp_order
      k = k + this%per_chirp_coefs(i) * ((x2 - r_center)**(i))
    enddo
  else
    do i = 1, this%chirp_order
      k = k + this%chirp_coefs(i) * ((z_center - x1)**(i))
    enddo
    do i = 1, this%per_chirp_order
      k = k + this%per_chirp_coefs(i) * ((r_center - x2)**(i))
    enddo
  endif

  rg1pb = 1.0

  if (this%if_boost) then
    uvel = sqrt(this%gamma**2-1.0_p_double)
    beta = uvel / this%gamma
    !transformations the same like for w
    if (this%propagation == p_forward) then
      ! equivalent to gamma * (1 - beta)
      rg1pb = this%gamma - sqrt(this%gamma**2-1.0_p_double)
    else
      ! equivalent to gamma * (1 + beta)
      rg1pb = this%gamma*(1.0_p_double + beta)
    endif
  endif

  ! get envelope
  select case ( this%per_type )

  case (p_plane)
    per_envelope_cylm = cos( k*(x1-z_center)*rg1pb + this%phase0)

  case (p_hermite_gaussian)

    ! get rayleigh range
    z0 = k * this%per_w0(1)**2 * 0.5_p_double

    ! Boost Rayleigh compression
    if (this%if_boost) z0 = z0/this%gamma

    ! calculate radius
    rho = x2 - r_center
    rho2 = rho**2

    ! calculate gaussian beam parameters
    z  = x1 - this%per_focus(1)

    ! In some situations (i.e. chirps) z0 may be 0 so z = 0 must be treated
    ! as a special case
    if ( z /= 0 ) then
      rWl2 = z0**2 / (z0**2 + z**2)
      curv = 0.5_p_double*rho2*z/(z**2 + z0**2)
      gouy_shift = ( this%per_tem_mode(1) + 1 ) * atan2( z, z0 )
    else
      rWl2 = 1
      curv = 0
      gouy_shift = 0
    endif
    rWu = sqrt(2 * rWl2) / this%per_w0(1)

    if (.not. this%if_boost) then
      ! note that this is different in 2D and 3D ( use 3D algorithm )
      per_envelope_cylm = sqrt(rWl2) * &
                        hermite( this%per_tem_mode(1), rho * rWu ) * &
                        exp( - rho2 * rWl2/this%per_w0(1)**2 ) * &
                        cos( k*(x1-z_center) + k * curv - gouy_shift + this%phase0 )
    else
      per_envelope_cylm = sqrt(rWl2) * &
                          hermite( this%per_tem_mode(1), rho * rWu ) * &
                          exp( - rho2 * rWl2/this%per_w0(1)**2 ) * &
                          cos( k*(x1-z_center)*rg1pb + k * curv/this%gamma - &
                            gouy_shift + this%phase0 )

    endif

  case (p_gaussian_asym)

    ! THIS DOES NOT WORK WITH CHIRPED PULSES
    ! (there is a possibility of a division by 0 in that case)

    ! calculate radius
    rho = x2 - r_center
    rho2 = rho**2

    ! get w0 and focsus as a smooth transition between the 2 regimes
    shp = smooth_heaviside( rho/this%per_asym_trans(1) )
    w0_asym = this%per_w0_asym(1, 1) + shp*( this%per_w0_asym(2, 1) - this%per_w0_asym(1, 1))
    focus_asym = this%per_focus_asym(1,1) + &
                  shp*(this%per_focus_asym(2,1) - this%per_focus_asym(1,1))

    ! get rayleigh range
    z0 = k * w0_asym**2 * 0.5_p_double

    ! Boost Rayleigh compression
    if (this%if_boost) z0 = z0/this%gamma

    ! calculate gaussian beam parameters
    z  = x1 - focus_asym

    Wl2 = 1.0_p_double + (z/z0)**2

    gouy_shift = atan2( z, z0 )

    if (.not. this%if_boost) then
      ! note that this is different in 2D and 3D
      ! TODO: Is the amplitude correct here?
      per_envelope_cylm = sqrt( 1.0_p_double /  sqrt(Wl2) ) * &
                          exp( - rho2/(w0_asym**2 * Wl2) ) * &
                          cos( k*(x1-z_center) + &
                          0.5_p_double*k*rho2*z/(z**2 + z0**2) - &
                          gouy_shift + this%phase0 )
    else
      per_envelope_cylm = sqrt( 1.0_p_double /  sqrt(Wl2) ) * &
                          exp( - rho2/(w0_asym**2 * Wl2) ) * &
                          cos( k*(x1-z_center)*rg1pb + &
                          0.5_p_double*k*rho2*z/(z**2 + z0**2)/this%gamma - &
                          gouy_shift + this%phase0 )

    endif

  case default

    ! This should never happen
    per_envelope_cylm = 0

  end select

end function per_envelope_cylm
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine create_zpulse_cyl_m( this, b_re, e_re, b_im, e_im, g_space, nx_p_min, gir_pos )
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------

  implicit none

  integer, parameter :: rank = 2

  class( t_zpulse_cyl_modes ), intent(inout) :: this
  type( t_vdf ), intent(inout) :: b_re, e_re
  type( t_vdf ), intent(inout) :: b_im, e_im
  type( t_space ), intent(in) :: g_space
  integer, dimension(:), intent(in) :: nx_p_min
  integer, intent(in) :: gir_pos

  real(p_double), dimension(2,rank) :: g_x_range
  real(p_double), dimension(rank) :: ldx

  integer :: i1, i2, gshift_i2

  real(p_double) :: zmin, z, z_2, dz_2
  real(p_double) :: rmin, r, r_2, dr_2
  real(p_double) :: lenv, lenv_2

  real(p_double) :: cos_pol, sin_pol, prop_sign

  real(p_double) :: amp

  ! Boost variables
  real(p_k_fld) :: gam_one_beta, f_t, o_t, n_e, o_e

  gshift_i2 = gir_pos - 2

  call get_x_bnd( g_space, g_x_range )
  ldx(1) = b_re%dx_(1)
  ldx(2) = b_re%dx_(2)

  zmin = g_x_range(p_lower,1) + real(nx_p_min(1)-1, p_double)*ldx(1)
  rmin = g_x_range(p_lower,2) + real(nx_p_min(2)-1, p_double)*ldx(2)

  cos_pol = cos( this%pol )
  sin_pol = sin( this%pol )
  amp = this%omega0 * this%a0

  dz_2 = ldx(1)/2.0_p_double
  dr_2 = ldx(2)/2.0_p_double

  select case( this%interpolation )
  case( p_linear, p_cubic )
    ! Nothing to change here
  case( p_quadratic, p_quartic )
    zmin = zmin + dz_2
    ! rmin = rmin + dr_2
  end select

  if ( this%propagation == p_forward ) then
    prop_sign = 1.0_p_double
  else
    prop_sign = -1.0_p_double
  endif

  ! loop through all grid cells and initialize
  ! no tilted pulses, just normal ones

  ! normal envelope
  do i1 = b_re%lbound(2), b_re%ubound(2)
    z = zmin + real(i1-1, p_double) * ldx(1)
    z_2 = z + dz_2

    lenv   = amp * this % lon_envelope( z   )
    lenv_2 = amp * this % lon_envelope( z_2 )

    do i2 = b_re%lbound(3), b_re%ubound(3)

      r = rmin + real(i2-1, p_double) * ldx(2)
      r_2 = r + dr_2

      ! b_re%f2(1,i1,i2) =  0.0_p_k_fld
      b_re%f2(2,i1,i2) = -sin_pol*lenv_2 * per_envelope_cylm( this, z_2, r   ) * prop_sign
      b_re%f2(3,i1,i2) =  cos_pol*lenv_2 * per_envelope_cylm( this, z_2, r_2 ) * prop_sign
      ! b_im%f2(1,i1,i2) =  0.0_p_k_fld
      b_im%f2(2,i1,i2) =  cos_pol*lenv_2 * per_envelope_cylm( this, z_2, r   ) * prop_sign
      b_im%f2(3,i1,i2) =  sin_pol*lenv_2 * per_envelope_cylm( this, z_2, r_2 ) * prop_sign

      ! e_re%f2(1,i1,i2) =  0.0_p_k_fld
      e_re%f2(2,i1,i2) =  cos_pol*lenv   * per_envelope_cylm( this, z  , r_2 )
      e_re%f2(3,i1,i2) =  sin_pol*lenv   * per_envelope_cylm( this, z  , r   )
      ! e_re%f2(1,i1,i2) =  0.0_p_k_fld
      e_im%f2(2,i1,i2) =  sin_pol*lenv   * per_envelope_cylm( this, z  , r_2 )
      e_im%f2(3,i1,i2) = -cos_pol*lenv   * per_envelope_cylm( this, z  , r   )

    enddo ! i2
  enddo   ! i1

  ! Boost fields
  if (this%if_boost) then

    SCR_ROOT('Boosting zpulse initialization...')

    ! equivalent to gamma * (1 - beta)
    if (this%propagation == p_forward) then
      gam_one_beta = real( this%gamma - sqrt(this%gamma**2-1.0_p_double), p_k_fld )
    else
      gam_one_beta = real( this%gamma + sqrt(this%gamma**2-1.0_p_double), p_k_fld )
    endif
    ! also boost 1 guard cell
    do i1=0, b_re%nx_(1)+1
      do i2=0, b_re%nx_(2)+1

        e_re%f2(2,i1,i2) = gam_one_beta * e_re%f2(2,i1,i2)
        e_re%f2(3,i1,i2) = gam_one_beta * e_re%f2(3,i1,i2)
        e_im%f2(2,i1,i2) = gam_one_beta * e_im%f2(2,i1,i2)
        e_im%f2(3,i1,i2) = gam_one_beta * e_im%f2(3,i1,i2)
        b_re%f2(2,i1,i2) = gam_one_beta * b_re%f2(2,i1,i2)
        b_re%f2(3,i1,i2) = gam_one_beta * b_re%f2(3,i1,i2)
        b_im%f2(2,i1,i2) = gam_one_beta * b_im%f2(2,i1,i2)
        b_im%f2(3,i1,i2) = gam_one_beta * b_im%f2(3,i1,i2)

      enddo
    enddo

  endif

  ! treat axis -- i2 == 1
  if (gshift_i2 < 0) then

    f_t = real(4.0_p_k_fld/3.0_p_k_fld, p_k_fld)
    o_t = real(1.0_p_k_fld/3.0_p_k_fld, p_k_fld)
    n_e = real(9.0_p_k_fld/8.0_p_k_fld, p_k_fld)
    o_e = real(1.0_p_k_fld/8.0_p_k_fld, p_k_fld)

    do i1 = 0, b_re%nx_(1)+1

      ! See comments in os-emf-solver-cyl-modes-yee.f03 for more explanation
      e_re%f2( 2, i1, 1 ) = f_t * e_re%f2( 2, i1, 2 ) - o_t * e_re%f2( 2, i1, 3 )
      e_re%f2( 2, i1, 0 ) = e_re%f2( 2, i1, 2 )
      e_im%f2( 2, i1, 1 ) = f_t * e_im%f2( 2, i1, 2 ) - o_t * e_im%f2( 2, i1, 3 )
      e_im%f2( 2, i1, 0 ) = e_im%f2( 2, i1, 2 )

      e_re%f2( 3, i1, 1 ) = 2.0_p_k_fld * e_im%f2( 2, i1, 1 ) - e_re%f2( 3, i1, 2 )
      e_re%f2( 3, i1, 0 ) = 2.0_p_k_fld * e_im%f2( 2, i1, 1 ) - e_re%f2( 3, i1, 3 )
      e_im%f2( 3, i1, 1 ) =-2.0_p_k_fld * e_re%f2( 2, i1, 1 ) - e_im%f2( 3, i1, 2 )
      e_im%f2( 3, i1, 0 ) =-2.0_p_k_fld * e_re%f2( 2, i1, 1 ) - e_im%f2( 3, i1, 3 )

      b_re%f2( 2, i1, 1 ) = b_re%f2( 2, i1, 2 )
      b_re%f2( 2, i1, 0 ) = b_re%f2( 2, i1, 3 )
      b_im%f2( 2, i1, 1 ) = b_im%f2( 2, i1, 2 )
      b_im%f2( 2, i1, 0 ) = b_im%f2( 2, i1, 3 )

      b_re%f2( 3, i1, 1 ) = n_e * b_im%f2( 2, i1, 2 ) - o_e * b_im%f2( 2, i1, 3 )
      b_re%f2( 3, i1, 0 ) = b_re%f2( 3, i1, 2 )
      b_im%f2( 3, i1, 1 ) =-n_e * b_re%f2( 2, i1, 2 ) + o_e * b_re%f2( 2, i1, 3 )
      b_im%f2( 3, i1, 0 ) = b_im%f2( 3, i1, 2 )

    enddo ! i1

  endif ! gshift_is < 0

end subroutine create_zpulse_cyl_m

!-----------------------------------------------------------------------------------------
! The cylindrical mode laser needs a special divergence correction algorithm, since
! the divergence is now a function of r
!-----------------------------------------------------------------------------------------
subroutine div_corr_zpulse_cyl_m( this, g_space, nx_p_min, g_nx, &
                                  b_re, e_re, b_im, e_im, gir_pos )

  implicit none

  integer, parameter :: rank = 2

  class( t_zpulse_cyl_modes ), intent(inout) :: this
  type( t_space ), intent(in) :: g_space
  integer, dimension(:), intent(in)  :: nx_p_min, g_nx
  type( t_vdf ), intent(inout) :: b_re, e_re, b_im, e_im

  integer, intent(in) :: gir_pos

  real(p_double), dimension(2,rank) :: g_x_range
  integer ::  gshift_i2, i2_0
  real(p_double) :: dx, dr, rp, rm, invr, invr2, rp2, rm2

  integer, dimension(rank) :: lnx
  integer, dimension(2, rank) :: lgc_num

  ! Boost variables
  real(p_double) :: gam_one_beta

  real(p_double) :: e1_re, b1_re, e1_im, b1_im    ! e1, b1 field correction
  real(p_double) :: e2p_re, e2m_re, b2p_re, b2m_re, e2p_im, e2m_im, b2p_im, b2m_im, e3_re, e3_im, b3_re, b3_im
  real(p_double) :: dx1

  real(p_double) :: dx2
  real(p_double) :: zmin, z, z_2, dz_2
  real(p_double) :: rmin, rcent, rcent2

  real(p_double) :: lenv, lenv_2, amp, cos_pol, sin_pol, prop_sign

  integer :: i1, i2

  integer :: idx1 = 1, idx2 = 2

  gshift_i2 = gir_pos - 2

  if ( gshift_i2 < 0 ) then
    i2_0 = 2
  else
    i2_0 = 1
  endif

  lnx(1:rank) = b_re%nx_(1:rank)
  lgc_num(:,1:rank) = b_re%gc_num()

  ! divergence correction is not required in plane waves
  ! TODO: Do we still need the i2 == 1 correction at the bottom even if false?
  if ((this%per_type /= p_plane) .and. (.not. this%no_div_corr)) then

    call get_x_bnd( g_space, g_x_range )

    cos_pol = cos( this%pol )
    sin_pol = sin( this%pol )
    amp = this%omega0 * this%a0
    if ( this%propagation == p_forward ) then
      prop_sign = 1.0_p_double
    else
      prop_sign = -1.0_p_double
    endif

    if (this%if_boost) then
      if (this%propagation == p_forward) then
        ! equivalent to gamma * (1 - beta)
        gam_one_beta = real( this%gamma - sqrt(this%gamma**2-1.0_p_double), p_k_fld )
      else
        ! equivalent to gamma * (1 + beta)
        gam_one_beta = real( this%gamma + sqrt(this%gamma**2-1.0_p_double), p_k_fld )

      endif
    else
      gam_one_beta = 0
    endif

    dx = real( e_re%dx_(1), p_double )
    dr = real( e_re%dx_(2), p_double )

    dx1 = dx
    dx2 = dr
    dz_2 = dx1 / 2.0_p_double

    zmin = g_x_range(p_lower,1) + real(nx_p_min(idx1)-1, p_double)*dx1
    rmin = g_x_range(p_lower,2) + real(nx_p_min(idx2)-1, p_double)*dx2

    select case( this%interpolation )
    case( p_linear, p_cubic )
      ! Nothing to change here
    case( p_quadratic, p_quartic )
      zmin = zmin + dz_2
      ! rmin = rmin + dx2 * 0.5_p_double
    end select

    ! must treat the axis separately, just like in the field solver

    do i2 = i2_0, lnx(idx2) !1, lnx(idx2)

      ! this is for the e-fields
      rp      = rmin + (real(i2,p_double) - 0.5_p_double)*dx2 ! e2
      rcent   = rmin + (real(i2,p_double) - 1.0_p_double)*dx2 ! e3, e1
      rm      = rmin + (real(i2,p_double) - 1.5_p_double)*dx2
      invr    = 1.0_p_double/rcent ! e3, e1

      ! this is for the b-fields
      rp2     = rmin + (real(i2,p_double)               )*dx2
      rcent2  = rmin + (real(i2,p_double) - 0.5_p_double)*dx2 ! b1, b3
      rm2     = rmin + (real(i2,p_double) - 1.0_p_double)*dx2 ! b2
      invr2   = 1.0_p_double/rcent2 ! b1, b3

      e1_re = 0.0_p_double
      b1_re = 0.0_p_double
      e1_im = 0.0_p_double
      b1_im = 0.0_p_double

      ! calculate correction from space on other nodes
      ! coordinates are referenced to local grid coordinates
      ! since we want to avoid unneccessary MPI calls, we'll integrate the out-of-node
      ! values using the known laser profile

      do i1 = g_nx(1) - nx_p_min(1) + 1, lnx(1) + 1, -1

        ! correct e1_re
        z = zmin + real(i1, p_double) * dx1
        lenv = amp * this % lon_envelope( z )

        if (lenv > 0.0_p_double) then

          e2p_re  = lenv * cos_pol * per_envelope_cylm( this, z, rp )
          e2m_re  = lenv * cos_pol * per_envelope_cylm( this, z, rm )
          e2p_im  = lenv * sin_pol * per_envelope_cylm( this, z, rp )
          e2m_im  = lenv * sin_pol * per_envelope_cylm( this, z, rm )

          e3_re =   lenv * sin_pol * per_envelope_cylm( this, z, rcent )
          e3_im = - lenv * cos_pol * per_envelope_cylm( this, z, rcent )

          e1_re = e1_re + invr*dx1*(rp*e2p_re - rm*e2m_re)/dx2 + invr*e3_im*dx1
          e1_im = e1_im + invr*dx1*(rp*e2p_im - rm*e2m_im)/dx2 - invr*e3_re*dx1

        else

          e1_re = 0.0_p_double
          e1_im = 0.0_p_double

        endif

        ! correct b1_re
        z_2 = z - dz_2
        lenv_2 = amp * this % lon_envelope( z_2 )

        if (lenv_2 > 0.0_p_double) then

          b2p_re = -lenv_2 * sin_pol * per_envelope_cylm( this, z_2, rp2 ) * prop_sign
          b2m_re = -lenv_2 * sin_pol * per_envelope_cylm( this, z_2, rm2 ) * prop_sign
          b2p_im =  lenv_2 * cos_pol * per_envelope_cylm( this, z_2, rp2 ) * prop_sign
          b2m_im =  lenv_2 * cos_pol * per_envelope_cylm( this, z_2, rm2 ) * prop_sign

          b3_re = lenv_2 * cos_pol * per_envelope_cylm( this, z_2, rcent2 ) * prop_sign
          b3_im = lenv_2 * sin_pol * per_envelope_cylm( this, z_2, rcent2 ) * prop_sign

          b1_re = b1_re + invr2*dx1*(rp2*b2p_re - rm2*b2m_re)/dx2 + invr2*b3_im*dx1
          b1_im = b1_im + invr2*dx1*(rp2*b2p_im - rm2*b2m_im)/dx2 - invr2*b3_re*dx1

        else

          b1_re = 0.0_p_double
          b1_im = 0.0_p_double

        endif

      enddo ! i1

      ! Boost fields
      if (this%if_boost) then
        e1_re = e1_re*gam_one_beta
        b1_re = b1_re*gam_one_beta
        e1_im = e1_im*gam_one_beta
        b1_im = b1_im*gam_one_beta
      endif

      e_re%f2(1, lnx(1)+1, i2) = real( e1_re, p_k_fld )
      b_re%f2(1, lnx(1)+1, i2) = real( b1_re, p_k_fld )
      e_im%f2(1, lnx(1)+1, i2) = real( e1_im, p_k_fld )
      b_im%f2(1, lnx(1)+1, i2) = real( b1_im, p_k_fld )

      ! local node correction
      do i1 = lnx(1), 1-lgc_num(1,1), -1

        z = zmin + real(i1, p_double) * dx1
        z_2 = z - dz_2

        if ( this % lon_envelope( z ) > 0.0_p_double ) then
          ! correct e1
          e2p_re = e_re%f2(2,i1+1,i2)
          e2m_re = e_re%f2(2,i1+1,i2-1)
          e2p_im = e_im%f2(2,i1+1,i2)
          e2m_im = e_im%f2(2,i1+1,i2-1)

          e3_re = e_re%f2(3,i1+1,i2)
          e3_im = e_im%f2(3,i1+1,i2)

          e1_re = e1_re + invr*(dx1/dx2)*(rp*e2p_re - rm*e2m_re) + invr*e3_im*dx1
          e1_im = e1_im + invr*(dx1/dx2)*(rp*e2p_im - rm*e2m_im) - invr*e3_re*dx1

          e_re%f2(1, i1, i2) = real( e1_re, p_k_fld )
          e_im%f2(1, i1, i2) = real( e1_im, p_k_fld )
        else
          e1_re = 0.0_p_double
          e1_im = 0.0_p_double
        endif

        if ( this % lon_envelope( z_2 ) > 0.0_p_double ) then
          ! correct b1
          b2p_re = b_re%f2(2,i1,i2+1)
          b2m_re = b_re%f2(2,i1,i2)
          b2p_im = b_im%f2(2,i1,i2+1)
          b2m_im = b_im%f2(2,i1,i2)

          b3_re = b_re%f2(3,i1,i2)
          b3_im = b_im%f2(3,i1,i2)

          b1_re = b1_re + invr2*(dx1/dx2)*(rp2*b2p_re - rm2*b2m_re)  + invr2*b3_im*dx1
          b1_im = b1_im + invr2*(dx1/dx2)*(rp2*b2p_im - rm2*b2m_im)  - invr2*b3_re*dx1

          b_re%f2(1, i1, i2) = real( b1_re, p_k_fld )
          b_im%f2(1, i1, i2) = real( b1_im, p_k_fld )
        else
          b1_re = 0.0_p_double
          b1_im = 0.0_p_double
        endif

      enddo ! i1

    enddo ! i2

  endif

  ! treat axis -- i2 == 1
  if (gshift_i2 < 0) then

    ! local node correction
    do i1 = lnx(1), 1-lgc_num(1,1), -1

      ! See comments in os-emf-solver-cyl-modes-yee.f03 for more explanation
      e_re%f2( 1, i1, 1 ) = - e_re%f2( 1, i1, 2 )
      e_re%f2( 1, i1, 0 ) = - e_re%f2( 1, i1, 3 )
      e_im%f2( 1, i1, 1 ) = - e_im%f2( 1, i1, 2 )
      e_im%f2( 1, i1, 0 ) = - e_im%f2( 1, i1, 3 )

      b_re%f2( 1, i1, 1 ) = 0.0_p_k_fld
      b_re%f2( 1, i1, 0 ) = - b_re%f2( 1, i1, 2 )
      b_im%f2( 1, i1, 1 ) = 0.0_p_k_fld
      b_im%f2( 1, i1, 0 ) = - b_im%f2( 1, i1, 2 )

    enddo ! i1

  endif ! gshift_is < 0

end subroutine div_corr_zpulse_cyl_m

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine read_input_zpulse_cyl_modes( this, input_file, g_space, bnd_con, periodic, grid, sim_options )

  implicit none

  class( t_zpulse_cyl_modes ), intent(inout) :: this
  class( t_input_file ), intent(inout) :: input_file
  type( t_space ), intent(in) :: g_space
  class (t_emf_bound), intent(in) :: bnd_con
  logical, dimension(:), intent(in) :: periodic
  class( t_grid ), intent(in) :: grid
  type( t_options ), intent(in) :: sim_options

  ! Call superclass method
  call this % t_zpulse % read_input( input_file, g_space, bnd_con, periodic, grid, sim_options )

  ! Only allow propagation in dimension 1
  if ( this%direction > 1 ) then
    if ( mpi_node() == 0 ) then
      write(0,*)  ""
      write(0,*)  "   Error reading zpulse parameters"
      write(0,*)  "   propagation only supported in direction 1"
      write(0,*)  "   when running with cylindrical modes"
      write(0,*)  "   aborting..."
    endif
    stop
  endif

  ! Only allow normal b_type
  if ( this%b_type /= p_zpulse_bnormal ) then
    if ( mpi_node() == 0 ) then
      write(0,*)  ""
      write(0,*)  "   Error reading zpulse parameters"
      write(0,*)  "   btype must be 'normal' when running with cyl_modes"
      write(0,*)  "   aborting..."
    endif
    stop
  endif

  ! Don't allow tilted pulses for now
  if ( this%if_lon_tilt ) then
    if ( mpi_node() == 0 ) then
      write(0,*)  ""
      write(0,*)  "   Error reading zpulse parameters"
      write(0,*)  "   tilt not implemented when running with cyl_modes"
      write(0,*)  "   aborting..."
    endif
    stop
  endif

  ! if the default per_center was used in regular zpulse, then 
  ! reset the default per_center to be 0 
  if( this%per_center(1) == 0.5_p_double * ( xmin( g_space, 2 ) + xmax( g_space, 2 ) ) ) then 
   ! set per_center to be on axis (otherwise the pulse will be invalid)
    this % per_center(1:2) = 0.0_p_double
  else
    if( this%per_center(1) .ne. 0 ) then 
      SCR_ROOT( "   ERROR: per_center(1) not equal to 0 was requested.")
      SCR_ROOT( "   When running in cyl_modes, per_center(1) must be 0,")
      SCR_ROOT( "   otherwise the pulse will have an infinite number of non-negligible ")
      SCR_ROOT( "   angular harmonics (instead of just m=1 when per_center = 0 for ")
      SCR_ROOT( "   cylindrically symmetric pulses)")
      stop
    endif
  endif


end subroutine read_input_zpulse_cyl_modes
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
! override the regular zpulse subroutine where dimensionality checks 
! are implemented. if a certain pulse type was only valid in 2D or 1D 
! and not in 3D, then logic to prevent that pulse from running would have
! to be implemented here. 
!-----------------------------------------------------------------------------------------
subroutine check_dimensionality_zpulse_cyl_modes( this )
  
  class( t_zpulse_cyl_modes ), intent(inout) :: this

  ! do nothing 

end subroutine check_dimensionality_zpulse_cyl_modes
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! create pulse and add it to e and b fields 
!-----------------------------------------------------------------------------------------
subroutine launch_zpulse_cyl_modes( this, emf, g_space, nx_p_min, g_nx, t, dt, no_co )

  implicit none

  class( t_zpulse_cyl_modes ), intent(inout) :: this
  class( t_emf ), intent(inout) :: emf
  type( t_space ), intent(in) :: g_space
  integer, intent(in), dimension(:) :: nx_p_min, g_nx
  real(p_double), intent(in) :: t
  real(p_double), intent(in) :: dt
  class( t_node_conf ), intent(in) :: no_co

  type( t_vdf ) :: e_pulse_re, b_pulse_re, e_pulse_im, b_pulse_im

  select type( emf )
  class is( t_emf_cyl_modes )

  if ( this%if_launch .and. ( t >= this%launch_time ) ) then

    ! create new vdfs to hold the pulse field
    call b_pulse_re % new( emf%b, zero = .true. )
    call e_pulse_re % new( emf%e, zero = .true. )
    call b_pulse_im % new( emf%b, zero = .true. )
    call e_pulse_im % new( emf%e, zero = .true. )

    ! We assume this%direction == 1

    ! create pulse field in b_pulse and e_pulse
    call create_zpulse_cyl_m( this, b_pulse_re, e_pulse_re, b_pulse_im, e_pulse_im, &
                              g_space, nx_p_min, emf%gix_pos(2) )

    ! correct divergence
    call div_corr_zpulse_cyl_m( this, g_space, nx_p_min, g_nx, b_pulse_re, e_pulse_re, &
                                b_pulse_im, e_pulse_im, emf%gix_pos(2) )

    ! add laser pulse to to real part of mode 1
    call add( emf%b_cyl_m%pf_re(1), b_pulse_re )
    call add( emf%e_cyl_m%pf_re(1), e_pulse_re )

    !add laser to imaginary part of mode 1
    call add( emf%b_cyl_m%pf_im(1), b_pulse_im )
    call add( emf%e_cyl_m%pf_im(1), e_pulse_im )

    ! take care of circular polarization
    if (this%pol_type /= 0) then

      ! clear pulse fields
      call b_pulse_re % zero()
      call e_pulse_re % zero()
      call b_pulse_im % zero()
      call e_pulse_im % zero()

      ! Get phase and polarization of second pulse
      this%phase0 = this%phase0 + real( sign( pi_2, real( this%pol_type, p_double )), p_double )
      this%pol    = this%pol    + real( pi_2, p_double )

      ! We assume this%direction == 1

      ! create pulse field in b_pulse and e_pulse
      call create_zpulse_cyl_m( this, b_pulse_re, e_pulse_re, b_pulse_im, e_pulse_im, &
                                g_space, nx_p_min, emf%gix_pos(2) )

      ! correct divergence
      call div_corr_zpulse_cyl_m( this, g_space, nx_p_min, g_nx, b_pulse_re, e_pulse_re, &
                                  b_pulse_im, e_pulse_im, emf%gix_pos(2) )

      ! add real, phase-shifted part
      call add( emf%b_cyl_m%pf_re(1), b_pulse_re )
      call add( emf%e_cyl_m%pf_re(1), e_pulse_re )

      ! add imaginary, phase-shifted part
      call add( emf%b_cyl_m%pf_im(1), b_pulse_im )
      call add( emf%e_cyl_m%pf_im(1), e_pulse_im )

    endif

    ! free b_pulse and e_pulse vdfs
    call b_pulse_re % cleanup()
    call e_pulse_re % cleanup()
    call b_pulse_im % cleanup()
    call e_pulse_im % cleanup()

    ! pulse has been launched, turn off if_launch
    this%if_launch = .false.

  endif ! this%if_launch .and. ( t >= this%launch_time )

  end select

end subroutine launch_zpulse_cyl_modes
!-----------------------------------------------------------------------------------------



end module m_zpulse_cyl_modes
