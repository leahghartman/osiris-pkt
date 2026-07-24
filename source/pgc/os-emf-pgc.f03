#include "os-config.h"
#include "os-preprocess.fpp"

module m_emf_pgc

#include "memory/memory.h"

use m_system
use m_parameters
use m_emf_pgc_define, only : t_emf_pgc, p_env_nm, p_env_n, p_env_np
use m_emf_pgc_define, only : p_emf_pgc_rst_id
use m_tridiag, only : t_tridiag, t_base_tridiag_solver
use m_tridiag, only : p_tridiag_lower, p_tridiag_middle, p_tridiag_upper
use m_logprof, only : begin_event, end_event

contains

!-------------------------------------------------------------------------------
! Tridiagonal solver routines:
!   - discritization coefficients (2d-cart, 2d-cyl and 2d)
!   - decomposition of coefficients
!   - tdma solver
!   - using reflective boundaries
!-------------------------------------------------------------------------------
subroutine coeff_2d_cart( tridiag_sys , omega , dt , dxi, if_periodic )

  implicit none

  ! input/output parameters
  class(t_base_tridiag_solver), intent(inout) :: tridiag_sys
  real(p_k_fld),                intent(in)    :: omega
  real(p_double),               intent(in)    :: dt , dxi
  logical,                      intent(in)    :: if_periodic

  ! auxiliary variables
  complex(p_double) :: alpha
  real(p_double)    :: dparam
  complex(p_double) :: a, b

  alpha = cmplx( 0.0_p_double, 1.0_p_double, kind=p_double ) * omega / dt
  dparam = 1.0_p_double / (dxi * dxi)

  tridiag_sys%coeff = 0.0_p_double

  ! inner grid values
  tridiag_sys%coeff( p_tridiag_lower , : ) = &
    cmplx(-0.5_p_double * dparam, kind=p_k_fld)
  tridiag_sys%coeff( p_tridiag_middle, : ) = &
    cmplx(alpha + dparam, kind=p_k_fld)
  tridiag_sys%coeff( p_tridiag_upper , : ) = &
    cmplx(-0.5_p_double * dparam, kind=p_k_fld)

  if( if_periodic ) then
    tridiag_sys%gam = &
      cmplx(-1.0_p_double, kind=p_k_fld) * &
      tridiag_sys%coeff( p_tridiag_middle, tridiag_sys%lb )

    a = tridiag_sys%coeff( p_tridiag_upper , tridiag_sys%ub )
    b = tridiag_sys%coeff( p_tridiag_lower , tridiag_sys%lb )

    tridiag_sys%coeff( p_tridiag_middle, tridiag_sys%lb ) = &
      tridiag_sys%coeff( p_tridiag_middle, tridiag_sys%lb ) - tridiag_sys%gam
    tridiag_sys%coeff( p_tridiag_middle, tridiag_sys%ub ) = &
      tridiag_sys%coeff( p_tridiag_middle, tridiag_sys%ub ) - &
      cmplx(a * b / tridiag_sys%gam, kind=p_k_fld)
  else
    ! lower boundary
    tridiag_sys%coeff( p_tridiag_lower, tridiag_sys%lb ) = &
      cmplx(0.0_p_double, kind=p_k_fld)
    tridiag_sys%coeff( p_tridiag_upper, tridiag_sys%lb ) = &
      cmplx(-1.0_p_double * dparam, kind=p_k_fld)
    ! upper boundary
    tridiag_sys%coeff( p_tridiag_lower, tridiag_sys%ub ) = &
      cmplx(-1.0_p_double * dparam, kind=p_k_fld)
    tridiag_sys%coeff( p_tridiag_upper, tridiag_sys%ub ) = &
      cmplx(0.0_p_double, kind=p_k_fld)
  endif

end subroutine coeff_2d_cart
!-------------------------------------------------------------------------------
subroutine coeff_2d_cyl( tridiag_sys , omega , dt , dxi,  if_periodic )

  implicit none

  ! input/output parameters
  class(t_base_tridiag_solver), intent(inout) :: tridiag_sys
  real(p_k_fld),                intent(in)    :: omega
  real(p_double),               intent(in)    :: dt , dxi
  logical,                      intent(in)    :: if_periodic

  ! auxiliary variables
  complex(p_double) :: alpha
  real(p_double) :: dparam
  integer :: j

  alpha = &
    cmplx( 0.0_p_double, 1.0_p_double, kind=p_double ) * &
    real(omega, p_double) / dt
  dparam = 1.0_p_double / (dxi * dxi)

  tridiag_sys%coeff = 0.0_p_double

  ! inner grid
  do j = tridiag_sys%lb+1, tridiag_sys%ub-1
    tridiag_sys%coeff(p_tridiag_lower, j) = &
      cmplx(-0.5_p_double * (1.0_p_double - 0.5_p_double / (j - 1.5)) * dparam, &
        kind=p_k_fld)
    tridiag_sys%coeff(p_tridiag_middle, j) = &
      cmplx(alpha + 1.0_p_double * dparam, &
        kind=p_k_fld)
    tridiag_sys%coeff(p_tridiag_upper, j) = &
      cmplx(-0.5_p_double * (1.0_p_double + 0.5_p_double / (j - 1.5)) * dparam, &
        kind=p_k_fld)
  enddo

  ! lower boundaries
  tridiag_sys%coeff(p_tridiag_upper, tridiag_sys%lb) = &
    cmplx(-1.0_p_double * dparam, kind=p_k_fld)
  tridiag_sys%coeff(p_tridiag_middle, tridiag_sys%lb) = &
    cmplx(alpha + dparam, kind=p_k_fld)

  ! upper boundaries
  tridiag_sys%coeff(p_tridiag_lower, tridiag_sys%ub) = &
    cmplx(-1.0_p_double * dparam, kind=p_k_fld)
  tridiag_sys%coeff(p_tridiag_middle, tridiag_sys%ub) = &
    cmplx(alpha + dparam, kind=p_k_fld)

end subroutine coeff_2d_cyl
!-------------------------------------------------------------------------------
subroutine coeff_3d_cart( tridiag_sys , omega , dt , dxi,  if_periodic )

  implicit none

  ! input/output parameters
  class(t_base_tridiag_solver), intent(inout) :: tridiag_sys
  real(p_k_fld),                intent(in)    :: omega
  real(p_double),               intent(in)    :: dt , dxi
  logical,                      intent(in)    :: if_periodic

  ! auxiliary variables
  complex(p_double) :: alpha
  real(p_double)    :: dparam
  complex(p_double) :: a, b

  alpha = cmplx( 0.0_p_double, 2.0_p_double, kind=p_double ) * omega / dt
  dparam = 1.0_p_double / (dxi**2.0_p_double)

  tridiag_sys%coeff( p_tridiag_lower , : ) = &
    cmplx(-0.5_p_double * dparam, kind=p_k_fld)
  tridiag_sys%coeff( p_tridiag_middle, : ) = &
    cmplx(alpha + dparam, kind=p_k_fld)
  tridiag_sys%coeff( p_tridiag_upper , : ) = &
    cmplx(-0.5_p_double * dparam, kind=p_k_fld)

  if( if_periodic ) then
    tridiag_sys%gam = cmplx(-1.0_p_double * (alpha + dparam), kind=p_k_fld)
    a = tridiag_sys%coeff( p_tridiag_upper , tridiag_sys%ub )
    b = tridiag_sys%coeff( p_tridiag_lower , tridiag_sys%lb )

    tridiag_sys%coeff( p_tridiag_middle, tridiag_sys%lb ) = &
      tridiag_sys%coeff( p_tridiag_middle, tridiag_sys%lb ) - tridiag_sys%gam
    tridiag_sys%coeff( p_tridiag_middle, tridiag_sys%ub ) = &
      tridiag_sys%coeff( p_tridiag_middle, tridiag_sys%ub ) - &
      cmplx(a * b / tridiag_sys%gam, kind=p_k_fld)
  else
    ! lower boundary
    tridiag_sys%coeff( p_tridiag_lower, tridiag_sys%lb ) = &
      cmplx(0.0_p_double, kind=p_k_fld)
    tridiag_sys%coeff( p_tridiag_upper, tridiag_sys%lb ) = &
      cmplx(-1.0_p_double * dparam, kind=p_k_fld)
    ! upper boundary
    tridiag_sys%coeff( p_tridiag_lower, tridiag_sys%ub ) = &
      cmplx(-1.0_p_double * dparam, kind=p_k_fld)
    tridiag_sys%coeff( p_tridiag_upper, tridiag_sys%ub ) = &
      cmplx(0.0_p_double, kind=p_k_fld)
  endif

end subroutine coeff_3d_cart
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! Setup for PGC solver:
!   - initialize laser envelope if not restarting
!   - init required quantaties like array for right-hand side of the envelope
!     equation (source)
!   - init discritization matrix solver with different coefficient functions
!-------------------------------------------------------------------------------
subroutine setup_pgc( this, g_xmin, gc_min, n_x_pmin, dx, dt, grid, no_co, &
                      restart )

  use m_pgc_laser,   only : launch_laser
  use m_grid_define, only : t_grid
  use m_node_conf,   only : t_node_conf, periodic

  implicit none

  class( t_emf_pgc ), intent(inout)                    :: this
  real(p_k_fld),      intent(in),   dimension(p_x_dim) :: g_xmin
  integer,            intent(in),   dimension(1:, 1:)  :: gc_min
  integer,            intent(in),   dimension(:)       :: n_x_pmin
  real(p_double),     intent(in),   dimension(:)       :: dx
  real(p_double),     intent(in)                       :: dt
  class( t_grid ),    intent(in)                       :: grid
  class( t_node_conf ),  intent(in)                    :: no_co
  logical,            intent(in)                       :: restart

  ! auxiliary variables
  logical :: if_periodic(p_x_dim)
  real( p_k_fld ), dimension(3, 5) :: smoothing_kernel

  ! define smoothing kernel
  smoothing_kernel( :, 1 ) = real((/  1.0d0, 2.0d0, 1.0d0 /), p_k_fld)
  smoothing_kernel( :, 2 ) = real((/  1.0d0, 2.0d0, 1.0d0 /), p_k_fld)
  smoothing_kernel( :, 3 ) = real((/  1.0d0, 2.0d0, 1.0d0 /), p_k_fld)
  smoothing_kernel( :, 4 ) = real((/  1.0d0, 2.0d0, 1.0d0 /), p_k_fld)
  smoothing_kernel( :, 5 ) = real((/ -5.0d0, 15.0d0,-5.0d0 /), p_k_fld)

  ! launch only laser if not restarting
  if ( .not.(restart) ) then
    call launch_laser( this, g_xmin, n_x_pmin, dx, dt )
  endif

  ! source term for right-hand-side of envelope equation
  call this%source%new(this%env, f_dim=1)

  if_periodic = periodic(no_co)

  select case ( p_x_dim )
  case(1)
    print * , "Error initializing PGC"
    print * , "Ponderomotive guiding center solver is not implemented in 1D."
    print * , "aborting..."
    stop
  case(2)
    ! check for cylindrical or cartesian coordinates
    select case ( this%t_emf%coordinates )
    case (p_cartesian)
      call this%disc_mat(2)%setup( &
        coeff_2d_cart, this%solver_type, no_co, grid, gc_min, 2 &
      )
    case (p_cylindrical_b)
      ! periodic systems are not allowed for cylindrical coordinates
      if( if_periodic(2) ) then
        print *, "Periodic boundaries in cylindrical geometry is not supported"
        print *, "aborting..."
        stop
      endif

      ! for cylindrical grid index 1 represents -0.5 * dr
      ! we require a shift of lower index by one as the envelope will be
      ! calculated from 0.5*dr ... Nr * dr + 0.5 * dr
      call this%disc_mat(2)%setup( &
        coeff_2d_cyl, this%solver_type, no_co, grid, gc_min, 2, l_shift=1 &
      )
    end select

    call this%disc_mat(2)%solver%coeff_func( &
      this%omega, dt, dx(2), if_periodic(2) &
    )
    call this%disc_mat(2)%solver%decompose()
    if( if_periodic(2) ) then
      call this%disc_mat(2)%periodic_preCalculation()
    endif
    ! pass properties from emf to tridiagonal systems
    this%disc_mat(2)%error_cleaning = this%error_cleaning
    call this%disc_mat(2)%setup_smoothing( &
      smoothing_kernel, this%pre_calc_smoothing, this%post_calc_smoothing &
    )
  case(3)
    ! coefficients and decomposition x2-direction
    call this%disc_mat(2)%setup( &
      coeff_3d_cart, this%solver_type, no_co, grid, gc_min, 2 &
    )
    call this%disc_mat(2)%solver%coeff_func( &
      this%omega, dt, dx(2), if_periodic(2) &
    )
    call this%disc_mat(2)%solver%decompose()
    if( if_periodic(2) ) then
      call this%disc_mat(2)%periodic_preCalculation()
    endif
    ! coefficients and decomposition x3-direction
    call this%disc_mat(3)%setup( &
      coeff_3d_cart, this%solver_type, no_co, grid, gc_min, 3 &
    )
    call this%disc_mat(3)%solver%coeff_func( &
      this%omega, dt, dx(3), if_periodic(3) &
    )
    call this%disc_mat(3)%solver%decompose()
    if( if_periodic(3) ) then
      call this%disc_mat(3)%periodic_preCalculation()
    endif

    ! pass properties from emf to tridiagonal systems
    this%disc_mat(2)%error_cleaning = this%error_cleaning
    this%disc_mat(3)%error_cleaning = this%error_cleaning
    call this%disc_mat(2)%setup_smoothing( &
      smoothing_kernel, this%pre_calc_smoothing, this%post_calc_smoothing &
    )
    call this%disc_mat(3)%setup_smoothing( &
      smoothing_kernel, this%pre_calc_smoothing, this%post_calc_smoothing &
    )
  end select

end subroutine setup_pgc
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! Update particle related quantaties
!   - updates the squared envelope value
!   - updates the ponderomotive force for the required time steps
!   - perform calcualtions over the physical domain -> requires communication
!     afterwards
!-------------------------------------------------------------------------------
subroutine update_pgc_part( this , dxi )

  implicit none

  ! input parameters
  class( t_emf_pgc ), intent(inout)               :: this
  real( p_double ),   intent(in),   dimension(1:) :: dxi
  ! auxiliary variables
  integer          :: i1, i2, i3
  integer          :: nx_min, nx_max, ny_min, ny_max, nz_min, nz_max
  real( p_k_fld ) :: dx(1:size(dxi))
  ! pointers - especially required for openMP parallelization
  integer, dimension(p_x_dim) :: lb

  complex( p_k_fld ), pointer :: a_2d_nm(:,:), a_2d_n(:,:), a_2d_np(:,:)
  complex( p_k_fld ), pointer :: a_3d_nm(:,:,:), a_3d_n(:,:,:), a_3d_np(:,:,:)

  real( p_k_fld ), pointer :: sqrt_2d_nm(:,:), sqrt_2d_n(:,:), sqrt_2d_np(:,:)
  real( p_k_fld ), pointer :: sqrt_3d_nm(:,:,:), sqrt_3d_n(:,:,:), sqrt_3d_np(:,:,:)

  real( p_k_fld ), pointer :: env_2d_mod(:,:), env_3d_mod(:,:,:)

  real( p_k_fld ), pointer :: f_2d_n(:,:,:), f_2d_np(:,:,:)
  real( p_k_fld ), pointer :: f_3d_n(:,:,:,:), f_3d_np(:,:,:,:)

  ! determine physical domain and dxi for each direction
  nx_min = 1
  nx_max = this%env%nx_(1)
  ny_min = 1
  ny_max = this%env%nx_(2)

  dx(1) = 0.5_p_k_fld / real( dxi(1), p_k_fld )
  dx(2) = 0.5_p_k_fld / real( dxi(2), p_k_fld )

  ! for 3d case
  if( p_x_dim == 3 ) then
    nz_min = 1
    nz_max = this%env%nx_(3)

    dx(3) = 0.5_p_k_fld / real( dxi(3), p_k_fld )
  endif

  ! pointer association for openMP
  select case(p_x_dim)
  case(2)
    lb(1) = this%env%lbound(dim=2)
    lb(2) = this%env%lbound(dim=3)

    a_2d_nm(lb(1):, lb(2):) => this%env%z2(p_env_nm, :, :)
    a_2d_n(lb(1):, lb(2):)  => this%env%z2(p_env_n, :, :)
    a_2d_np(lb(1):, lb(2):) => this%env%z2(p_env_np, :, :)

    sqrt_2d_nm(lb(1):, lb(2):) => this%sqrt_nm%f2(1, :, :)
    sqrt_2d_n(lb(1):, lb(2):)  => this%sqrt_n%f2(1, :, :)
    sqrt_2d_np(lb(1):, lb(2):) => this%sqrt_np%f2(1, :, :)

    env_2d_mod(lb(1):, lb(2):) => this%env_mod%f2(1, :, :)

    f_2d_n(1:, lb(1):, lb(2):) => this%f_n%f2
    f_2d_np(1:, lb(1):, lb(2):) => this%f_np%f2
  case(3)
    lb(1) = this%env%lbound(dim=2)
    lb(2) = this%env%lbound(dim=3)
    lb(3) = this%env%lbound(dim=4)

    a_3d_nm(lb(1):, lb(2):, lb(3):) => this%env%z3(p_env_nm, :, :, :)
    a_3d_n(lb(1):, lb(2):, lb(3):)  => this%env%z3(p_env_n, :, :, :)
    a_3d_np(lb(1):, lb(2):, lb(3):) => this%env%z3(p_env_np, :, :, :)

    sqrt_3d_nm(lb(1):, lb(2):, lb(3):) => this%sqrt_nm%f3(1, :, :, :)
    sqrt_3d_n(lb(1):, lb(2):, lb(3):)  => this%sqrt_n%f3(1, :, :, :)
    sqrt_3d_np(lb(1):, lb(2):, lb(3):) => this%sqrt_np%f3(1, :, :, :)

    env_3d_mod(lb(1):, lb(2):, lb(3):) => this%env_mod%f3(1, :, :, :)
    f_3d_n(1:, lb(1):, lb(2):, lb(3):) => this%f_n%f3
    f_3d_np(1:, lb(1):, lb(2):, lb(3):) => this%f_np%f3
  end select


  select case(p_x_dim)
  case(2)
    ! envelope squared and |env|
    !$omp parallel do &
    !$omp&  private(i1) &
    !$omp&  shared(nx_min, nx_max, ny_min, ny_max) &
    !$omp&  shared(sqrt_2d_nm, sqrt_2d_n, sqrt_2d_np, env_2d_mod) &
    !$omp&  shared(a_2d_nm, a_2d_n, a_2d_np)
    do i2 = ny_min-1, ny_max+1
      do i1 = nx_min-1, nx_max+1
        sqrt_2d_nm(i1, i2) = abs(a_2d_nm(i1, i2)) ** 2.0_p_k_fld
        sqrt_2d_n(i1, i2)  = abs(a_2d_n(i1, i2)) ** 2.0_p_k_fld
        sqrt_2d_np(i1, i2) = abs(a_2d_np(i1, i2)) ** 2.0_p_k_fld
        env_2d_mod(i1, i2) = abs(a_2d_n(i1, i2))
      enddo
    enddo
    !$omp end parallel do

    ! calculate ponderomotive force
    !$omp parallel do &
    !$omp&  private(i1) &
    !$omp&  shared(nx_min, nx_max, ny_min, ny_max) &
    !$omp&  shared(f_2d_n, f_2d_np, sqrt_2d_n, sqrt_2d_np)
    do i2 = ny_min, ny_max
      do i1 = nx_min, nx_max
        f_2d_n(1, i1, i2) = dx(1) * (sqrt_2d_n(i1+1, i2) - sqrt_2d_n(i1-1, i2))
        f_2d_n(2, i1, i2) = dx(2) * (sqrt_2d_n(i1, i2+1) - sqrt_2d_n(i1, i2-1))

        f_2d_np(1, i1, i2) = dx(1) * (sqrt_2d_np(i1+1, i2) - sqrt_2d_np(i1-1, i2))
        f_2d_np(2, i1, i2) = dx(2) * (sqrt_2d_np(i1, i2+1) - sqrt_2d_np(i1, i2-1))
      enddo
    enddo
    !$omp end parallel do
  case(3)
    ! envelope squared and |env|
    !$omp parallel do &
    !$omp&  private(i1, i2) &
    !$omp&  shared(nx_min, nx_max, ny_min, ny_max, nz_min, nz_max) &
    !$omp&  shared(sqrt_3d_nm, sqrt_3d_n, sqrt_3d_np, env_3d_mod) &
    !$omp&  shared(a_3d_nm, a_3d_n, a_3d_np)
    do i3 = nz_min-1, nz_max+1
      do i2 = ny_min-1, ny_max+1
        do i1 = nx_min-1, nx_max+1
          sqrt_3d_nm(i1, i2, i3) = abs(a_3d_nm(i1, i2, i3)) ** 2.0_p_k_fld
          sqrt_3d_n(i1, i2, i3)  = abs(a_3d_n(i1, i2, i3)) ** 2.0_p_k_fld
          sqrt_3d_np(i1, i2, i3) = abs(a_3d_np(i1, i2, i3)) ** 2.0_p_k_fld
          env_3d_mod(i1, i2, i3) = abs(a_3d_n(i1, i2, i3))
        enddo
      enddo
    enddo
    !$omp end parallel do

    ! calculate ponderomotive force
    !$omp parallel do &
    !$omp&  private(i1, i2) &
    !$omp&  shared(nx_min, nx_max, ny_min, ny_max, nz_min, nz_max) &
    !$omp&  shared(f_3d_n, f_3d_np, sqrt_3d_n, sqrt_3d_np)
    do i3 = nz_min, nz_max
      do i2 = ny_min, ny_max
        do i1 = nx_min, nx_max
          f_3d_n(1, i1, i2, i3) = &
            dx(1) * (sqrt_3d_n(i1+1, i2, i3) - sqrt_3d_n(i1-1, i2, i3))
          f_3d_n(2, i1, i2, i3) = &
            dx(2) * (sqrt_3d_n(i1, i2+1, i3) - sqrt_3d_n(i1, i2-1, i3))
          f_3d_n(3, i1, i2, i3) = &
            dx(3) * (sqrt_3d_n(i1, i2, i3+1) - sqrt_3d_n(i1, i2, i3-1))

          f_3d_np(1, i1, i2, i3) = &
            dx(1) * (sqrt_3d_np(i1+1, i2, i3) - sqrt_3d_np(i1-1, i2, i3))
          f_3d_np(2, i1, i2, i3) = &
            dx(2) * (sqrt_3d_np(i1, i2+1, i3) - sqrt_3d_np(i1, i2-1, i3))
          f_3d_np(3, i1, i2, i3) = &
            dx(3) * (sqrt_3d_np(i1, i2, i3+1) - sqrt_3d_np(i1, i2, i3-1))
        enddo
      enddo
    enddo
    !$omp end parallel do
  end select

  ! nullify pointers used throughout the calculations
  select case(p_x_dim)
  case(2)
    nullify(a_2d_nm, a_2d_n, a_2d_np, sqrt_2d_nm, sqrt_2d_n, sqrt_2d_np)
    nullify(env_2d_mod, f_2d_n, f_2d_np)
  case(3)
    nullify(a_3d_nm, a_3d_n, a_3d_np, sqrt_3d_nm, sqrt_3d_n, sqrt_3d_np)
    nullify(env_3d_mod, f_3d_n, f_3d_np)
  end select

end subroutine update_pgc_part
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! Initialize longitudinal boundaries
!   - use a mask in longitudinal direction for envelope and ponderomotive force
!   - check how big the mask on a local node has to be to be applied
!-------------------------------------------------------------------------------
subroutine init_lon_bounds(this, no_co, grid)

  use m_emf_pgc_define, only : p_mask_length, p_mask_values
  use m_node_conf,      only : t_node_conf, p_upper, p_lower!, my_aid
  use m_grid_define,    only : t_grid
  use m_grid,           only : nx_p, g_nx

  implicit none

  class( t_emf_pgc ),   intent( inout ) :: this
  class( t_node_conf ), intent( in )    :: no_co
  class( t_grid ),      intent( in )    :: grid

  integer, parameter :: bnd_dir = 1       ! direction of boundaries
  ! paramerter to access grid info ( lower, upper grid index and num of cells)
  integer, parameter :: l_ind = 1, u_ind = 2, num = 3

  ! auxiliary
  integer, dimension(1:3, 1:p_x_dim) :: local_grid
  integer :: mask_left, mask_right
  integer :: grid_length

  integer :: i1

  ! local_grid = nx_p(grid, no_co, no_co%my_aid())
  local_grid = nx_p(grid, no_co, no_co%my_aid())
  ! - grid_length <-> length of global grid in x1
  grid_length = g_nx(grid, bnd_dir)

  ! determines the number of mask to apply
  mask_left = p_mask_length - local_grid(l_ind, bnd_dir) + 1
  mask_right = local_grid(u_ind, bnd_dir) - (grid_length - p_mask_length)

  ! truncate masks for each side
  if( mask_left > local_grid(num, bnd_dir) ) then
    mask_left = local_grid(num, bnd_dir)
  endif
  if( mask_right > local_grid(num, bnd_dir) ) then
    mask_right = local_grid(num, bnd_dir)
  endif

  ! apply left boundary mask
  if( mask_left > 0 ) then
    this%lb_apply = .true.

    ! include guard cells in filter mask only on the outer nodes
    if( no_co%on_edge( bnd_dir, p_lower) ) then
      this%bnd_lb_low = this%env%lbound(2)
      this%bnd_lb_upp = mask_left
    else
      this%bnd_lb_low = 1
      this%bnd_lb_upp = mask_left
    endif

    call alloc(this%bd_val_lb, (/ this%bnd_lb_low /), (/ this%bnd_lb_upp /))

    ! fill filter mask - guard cells are zero
    if( no_co%on_edge( bnd_dir, p_lower) ) then
      do i1 = this%bnd_lb_low, this%bnd_lb_upp
        if( i1 < 1 ) then
          this%bd_val_lb(i1) = 0.0
        else
          this%bd_val_lb(i1) = p_mask_values(i1)
        endif
      enddo
    else
      do i1 = this%bnd_lb_low, this%bnd_lb_upp
        this%bd_val_lb(i1) = p_mask_values( i1 + (p_mask_length - mask_left) )
      enddo
    endif
  else
    this%lb_apply = .false.
  endif

  ! apply right boundary mask
  if( mask_right > 0 ) then
    this%ub_apply = .true.

    ! include guard cells in filter mask only on the outer nodes
    if( no_co%on_edge( bnd_dir, p_upper) ) then
      this%bnd_ub_low = this%env%nx_(1) - mask_right + 1
      this%bnd_ub_upp = this%env%ubound(2)
    else
      this%bnd_ub_low = this%env%nx_(1) - mask_right + 1
      this%bnd_ub_upp = this%env%nx_(1)
    endif

    call alloc(this%bd_val_ub, (/ this%bnd_ub_low /), (/ this%bnd_ub_upp /))

    ! fill filter mask - guard cells are zero
    if( no_co%on_edge( bnd_dir, p_upper) ) then
      do i1 = this%bnd_ub_low, this%bnd_ub_upp
        if( i1 > this%env%nx_(1) ) then
          this%bd_val_ub(i1) = 0.0
        else
          this%bd_val_ub(i1) = p_mask_values((this%env%nx_(1) - i1) + 1)
        endif
      enddo
    else
      do i1 = this%bnd_ub_low, this%bnd_ub_upp
        this%bd_val_ub(i1) = p_mask_values(p_mask_length - (this%env%nx_(1) - i1))
      enddo
    endif
  else
    this%ub_apply = .false.
  endif

end subroutine init_lon_bounds
!-------------------------------------------------------------------------------
! Checks if envelope overlaps with longitudinal mask
!   - envelope should have a value lower than `p_env_bnd_overlap` otherwise
!     simulation stops
!-------------------------------------------------------------------------------
subroutine check_mask_env_overlap( this, low, upp )

  use m_emf_pgc_define, only : p_env_bnd_overlap

  implicit none

  ! subroutine parameters
  class(t_emf_pgc) , intent(inout) :: this
  integer, intent(in) :: low, upp

  ! auxiliary variables
  integer :: i1, i2, i3
  integer :: ny_min, ny_max, nz_min, nz_max
  logical :: unfulfilled_condition = .false.
  real( p_k_fld ), pointer :: env_mod_2d(:,:,:), env_mod_3d(:,:,:,:)

  select case( p_x_dim )
  case( 2 )

    env_mod_2d => this%env_mod%f2

    ny_min = 1
    ny_max = this%env_mod%nx_(2)

    do i2 = ny_min, ny_max
      do i1 = low, upp
        if( env_mod_2d(1, i1, i2) > p_env_bnd_overlap ) then
          unfulfilled_condition = .true.
          exit
        endif
      enddo
    enddo

  case( 3 )

    env_mod_3d => this%env_mod%f3

    ny_min = 1
    ny_max = this%env_mod%nx_(2)
    nz_min = 1
    nz_max = this%env_mod%nx_(3)

    do i3 = nz_min, nz_max
      do i2 = ny_min, ny_max
        do i1 = low, upp
          if( env_mod_3d(1, i1, i2, i3) > p_env_bnd_overlap ) then
            unfulfilled_condition = .true.
            exit
          endif
        enddo
      enddo
    enddo

  end select

  if( unfulfilled_condition ) then
    print *, "(* error::pgc *) wrong laser initialization"
    print *, "(* error::pgc *) laser envelope to close to x1 boundaries"
    print *, "aborting..."
    stop
  endif

end subroutine check_mask_env_overlap
!-------------------------------------------------------------------------------
! Applies mask values in longitudinal direction - low-level subroutine
!   - goes through mask array and applies it onto a vdf
!   - works for complex valued vdf and real valued vdf
!-------------------------------------------------------------------------------
subroutine apply_lon_mask( vdf_arr, mask, low, upp )

  use m_vdf_define,  only : t_vdf
  use m_cmplx_vdf,   only : t_cmplx_vdf

  implicit none

  class( t_vdf ),   intent(inout) :: vdf_arr
  integer,          intent(in)    :: low, upp
  real( p_double ), intent(in),   dimension(low:upp)  :: mask

  integer :: i1, i2, i3, f
  integer :: ny_min, ny_max, nz_min, nz_max, f_min, f_max
  real( p_k_fld ),    pointer :: vdf_2d(:,:,:), vdf_3d(:,:,:,:)
  complex( p_k_fld ), pointer :: zvdf_2d(:,:,:), zvdf_3d(:,:,:,:)

  f_min = 1
  f_max = vdf_arr%f_dim()

  select type (vdf_arr)
  ! for complex valued vdfs
  class is ( t_cmplx_vdf )
    select case( p_x_dim )
    case( 2 )
      zvdf_2d => vdf_arr%z2

      ny_min = lbound(zvdf_2d, dim=3)
      ny_max = ubound(zvdf_2d, dim=3)

      do i2 = ny_min, ny_max
        do i1 = low, upp
          do f = f_min, f_max
            zvdf_2d(f, i1, i2) = &
              zvdf_2d(f, i1, i2) * cmplx(mask(i1), kind=p_k_fld)
          enddo
        enddo
      enddo
    case( 3 )
      zvdf_3d => vdf_arr%z3

      ny_min = lbound(zvdf_3d, dim=3)
      ny_max = ubound(zvdf_3d, dim=3)
      nz_min = lbound(zvdf_3d, dim=4)
      nz_max = ubound(zvdf_3d, dim=4)

      do i3 = nz_min, nz_max
        do i2 = ny_min, ny_max
          do i1 = low, upp
            do f = f_min, f_max
              zvdf_3d(f, i1, i2, i3) = &
                zvdf_3d(f, i1, i2, i3) * cmplx(mask(i1), kind=p_k_fld)
            enddo
          enddo
        enddo
      enddo
    end select
  ! for regular vdfs
  class default
    select case( p_x_dim )
    case( 2 )
      vdf_2d => vdf_arr%f2

      ny_min = lbound(vdf_2d, dim=3)
      ny_max = ubound(vdf_2d, dim=3)

      do i2 = ny_min, ny_max
        do i1 = low, upp
          do f = f_min, f_max
            vdf_2d(f, i1, i2) = &
              vdf_2d(f, i1, i2) * real(mask(i1), p_k_fld)
          enddo
        enddo
      enddo
    case( 3 )
      vdf_3d => vdf_arr%f3

      ny_min = lbound(vdf_3d, dim=3)
      ny_max = ubound(vdf_3d, dim=3)
      nz_min = lbound(vdf_3d, dim=4)
      nz_max = ubound(vdf_3d, dim=4)

      do i3 = nz_min, nz_max
        do i2 = ny_min, ny_max
          do i1 = low, upp
            do f = f_min, f_max
              vdf_3d(f, i1, i2, i3) = &
                vdf_3d(f, i1, i2, i3) * real(mask(i1), p_k_fld)
            enddo
          enddo
        enddo
      enddo
    end select
  end select

end subroutine apply_lon_mask
!-------------------------------------------------------------------------------
! Update axial envelope values
!   - for cylindrical geometries the physical grid is given only for 2 - Nr in
!     transversal direction
!   - update is done through mirroring values into guard cells from index 2 on
!-------------------------------------------------------------------------------
subroutine update_axial_bnd( this )

  use m_tridiag, only : root

  implicit none

  class(t_emf_pgc), intent(inout) :: this

  integer :: f, i, j
  integer :: mirror_point

  ! mirroring index (env[mirror_point - 1] == env[mirror_point])
  mirror_point = this%disc_mat(2)%lb

  ! update only required on root node - closest to axial axis
  if( root(this%disc_mat(2)%nconf) ) then
    do j = this%env%lbound(dim=3), mirror_point - 1
      do i = this%env%lbound(dim=2), this%env%ubound(dim=2)
        do f = this%env%lbound(dim=1), this%env%ubound(dim=1)
          this%env%z2(f, i, j) = this%env%z2(f, i, 2 * mirror_point - j - 1)
        enddo
      enddo
    enddo
  endif

end subroutine update_axial_bnd
!-------------------------------------------------------------------------------
! Applies mask values in longitudinal direction for the envelope
!   - checks if boundaries have to be applied on node and applies them
!-------------------------------------------------------------------------------
subroutine apply_lon_bnds_env( this )

  implicit none

  class(t_emf_pgc) , intent(inout) :: this

  if( this%lb_apply ) then
    call apply_lon_mask( this%env, this%bd_val_lb, &
                         this%bnd_lb_low, this%bnd_lb_upp )
  endif

  if( this%ub_apply ) then
    call apply_lon_mask( this%env, this%bd_val_ub,  &
                         this%bnd_ub_low, this%bnd_ub_upp )
  endif

end subroutine apply_lon_bnds_env
!-------------------------------------------------------------------------------
! Applies mask values in longitudinal direction for the ponderomotive force
!   - checks if boundaries have to be applied on node and applies them
!-------------------------------------------------------------------------------
subroutine apply_lon_bnds_fp( this )

  implicit none

  class(t_emf_pgc) , intent(inout) :: this

  if( this%lb_apply ) then
    call apply_lon_mask( this%f_np, this%bd_val_lb, &
                         this%bnd_lb_low, this%bnd_lb_upp )
    call apply_lon_mask( this%f_n, this%bd_val_lb,  &
                         this%bnd_lb_low, this%bnd_lb_upp )
  endif

  if( this%ub_apply ) then
    call apply_lon_mask( this%f_np, this%bd_val_ub, &
                         this%bnd_ub_low, this%bnd_ub_upp )
    call apply_lon_mask( this%f_n, this%bd_val_ub,  &
                         this%bnd_ub_low, this%bnd_ub_upp )
  endif

end subroutine apply_lon_bnds_fp
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! Normalize charge for cylindrical geometry
!-------------------------------------------------------------------------------
subroutine norm_charge_cyl_pgc( chi, gir_pos, gdr )

  use m_vdf_define, only : t_vdf

  implicit none

  class( t_vdf ),   intent(inout) :: chi
  integer,          intent(in)    :: gir_pos ! local grid position on global grid
  real( p_double ), intent(in)    :: gdr

  integer         :: ir
  real( p_k_fld ) :: r, dr

  ! normalize for 'ring' particles
  ! Note that the lower radial spatial boundary of a cylindrical geometry
  ! simulation is always -dr/2, where dr is the radial cell size. Also note
  ! that charge beyond axial boundary is reversed since r is considered
  ! to be < 0.

  dr = real( gdr, p_k_fld )
  do ir = lbound( chi%f2, 3), ubound( chi%f2, 3)
     r = ( (ir+ gir_pos - 2) - 0.5_p_k_fld ) * dr
     chi%f2(1,:,ir) = chi%f2(1,:,ir) / r
  enddo

  ! Fold axial guard cells back into simulation space
  if ( gir_pos == 1 ) then
    do ir = 0, (1 - lbound( chi%f2, 3))
      chi%f2(1,:,ir+2) = chi%f2(1,:,ir+2) + chi%f2(1,:,1-ir)
      chi%f2(1,:,1-ir) = chi%f2(1,:,ir+2)
    enddo
  endif

end subroutine norm_charge_cyl_pgc
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! Advancing interface for envelope in 2d
!   - calls after construction of source term tridiagonal solver
!   - passes solution to corresponding time fields in envelope
!-------------------------------------------------------------------------------
subroutine advance_env_2d( this )

  implicit none

  ! input/output parameters
  class( t_emf_pgc ), intent(inout)  ::  this

  ! auxiliary variables
  integer :: i, j, nx_min, ny_min, nx_max, ny_max
  complex( p_k_fld ), pointer :: a_nm(:, :), a_n(:, :), a_np(:, :), source(:, :)
  integer :: imin, jmin, jmax

  nx_min = 1
  nx_max = this%env%nx_(1)

  ny_min = this%disc_mat(2)%ext_min
  ny_max = this%disc_mat(2)%ext_max

  call this%disc_mat(2)%solve(this%source%z2(1, nx_min:nx_max, ny_min:ny_max))

  imin = this%source%lbound(dim=2)

  jmin = this%source%lbound(dim=3)
  jmax = this%source%ubound(dim=3)

  a_nm(imin:, jmin:)   => this%env%z2(p_env_nm, :, :)
  a_n(imin:, jmin:)    => this%env%z2(p_env_n, :, :)
  a_np(imin:, jmin:)   => this%env%z2(p_env_np, :, :)
  source(imin:, jmin:) => this%source%z2( 1, :, :)

  ! calculate envelope at corresponding time step
  !$omp parallel do &
  !$omp&  private(i) &
  !$omp&  shared(jmin, jmax, a_nm, a_n, a_np, source)
  do j = jmin, jmax
    do i = 1, nx_max
      a_nm(i, j) = a_n(i, j)
      a_n(i, j)  = a_np(i, j)
      a_np(i, j) = source(i, j)
    enddo
  enddo
  !$omp end parallel do

end subroutine advance_env_2d
!-------------------------------------------------------------------------------
! Construct the right-hand-side - for 2d envelope propagation
!-------------------------------------------------------------------------------
subroutine build_source_2d( this , dx , dt , coordinates , if_periodic )

  use m_tridiag, only : root, last

  implicit none

  class( t_emf_pgc ), intent(inout)  ::  this
  real( p_double ), dimension(:), intent(in) :: dx
  real( p_double ), intent(in)  ::  dt
  integer, intent(in) :: coordinates
  logical, intent(in) :: if_periodic

  ! auxiliary parameters
  complex(p_double) :: alpha, beta
  integer :: i, j, grid_offset
  integer, dimension(2) :: lb, grid_lb, grid_ub
  real(p_double) :: dxi, ddy
  ! pointer helper
  real( p_k_fld ),    pointer :: chi(:,:)
  complex( p_k_fld ), pointer :: src(:,:)
  complex( p_k_fld ), pointer :: a_nm(:,:), a_n(:,:), a_np(:,:)

  ! coefficient of envelope equation
  alpha = cmplx(0.0_p_double, 1.0_p_double, kind=p_double) * this%omega / dt
  beta = cmplx(0.0_p_double, 1.0_p_double, kind=p_double) / this%omega

  ! finite difference coefficient
  dxi = 0.5_p_k_fld / dx(1)             ! for CDS - use 2 * \delta\xi
  ddy = 1.0_p_k_fld / (dx(2) * dx(2))

  ! required for pointer arithmetics
  lb(1) = this%source%lbound(dim=2)
  lb(2) = this%source%lbound(dim=3)

  ! physical grid
  ! - transversal grid is defined by tridiagonal system
  ! - longitudinal grid is given by grid length `nx` of vdf
  grid_lb(1) = 1
  grid_lb(2) = this%disc_mat(2)%ext_min

  grid_ub(1) = this%source%nx_(1)
  grid_ub(2) = this%disc_mat(2)%ext_max


  ! pointer for better reading and for openMP
  chi(lb(1):, lb(2):) => this%chi(1)%f2(1, :, :)
  src(lb(1):, lb(2):) => this%source%z2(1, :, :)
  a_nm(lb(1):, lb(2):) => this%env%z2(p_env_nm, :, :)
  a_n(lb(1):, lb(2):) => this%env%z2(p_env_n, :, :)
  a_np(lb(1):, lb(2):) => this%env%z2(p_env_np, :, :)


  select case ( coordinates )
  case(p_cartesian)
    ! inner grid
    !$omp parallel do &
    !$omp&  private(i) &
    !$omp&  shared(grid_lb, grid_ub, a_np, a_n, a_nm, src, chi) &
    !$omp&  shared(alpha, beta, dxi, ddy)
    do j = grid_lb(2), grid_ub(2)
      do i = grid_lb(1), grid_ub(1)
        src(i, j) = &
          chi(i, j) * a_np(i, j) - &
          beta * dxi * ( chi(i+1, j) * a_np(i+1, j) - chi(i-1, j) * a_np(i-1, j) ) + &
          alpha * a_n(i, j) + &
          0.5_p_k_fld * ddy * ( &
            a_n(i, j-1) - &
            2.0_p_k_fld * a_n(i, j) + &
            a_n(i, j+1) &
          ) + &
          beta * dxi * ( &
            ddy * ( a_np(i+1, j+1) - 2.0_p_k_fld * a_np(i+1, j) + a_np(i+1, j-1) ) - &
            ddy * ( a_np(i-1, j+1) - 2.0_p_k_fld * a_np(i-1, j) + a_np(i-1, j-1) ) &
          )
      enddo
    enddo
    !$omp end parallel do

    if( .not. if_periodic ) then
      ! lower boundary
      if( root(this%disc_mat(2)%nconf) ) then
        !$omp parallel do &
        !$omp&  shared(grid_lb, grid_ub, a_np, a_n, a_nm, src, chi) &
        !$omp&  shared(alpha, beta, dxi, ddy)
        do i = grid_lb(1), grid_ub(1)
          src(i, grid_lb(2)) = &
            chi(i, grid_lb(2)) * a_np(i, grid_lb(2)) - &
            beta * dxi * ( &
              chi(i+1, grid_lb(2)) * a_np(i+1, grid_lb(2)) - &
              chi(i-1, grid_lb(2)) * a_np(i-1, grid_lb(2)) &
            ) + &
            alpha * a_n(i, grid_lb(2)) + &
            ddy * ( a_n(i, grid_lb(2)+1) - a_n(i, grid_lb(2)) ) + &
            beta * dxi * ( &
              2.0_p_k_fld * ddy * ( &
                a_np(i+1, grid_lb(2)+1) - a_np(i+1, grid_lb(2)) &
              ) - &
              2.0_p_k_fld * ddy * ( &
                a_np(i-1, grid_lb(2)+1) - a_np(i-1, grid_lb(2)) &
              ) &
            )
          enddo
          !$omp end parallel do
      endif
      ! upper boundary
      if( last(this%disc_mat(2)%nconf) ) then
        !$omp parallel do &
        !$omp&  shared(grid_lb, grid_ub, a_np, a_n, a_nm, src, chi) &
        !$omp&  shared(alpha, beta, dxi, ddy)
        do i = grid_lb(1), grid_ub(1)
          src(i, grid_ub(2)) = &
            chi(i, grid_ub(2)) * a_np(i, grid_ub(2)) - &
            beta * dxi * ( &
              chi(i+1, grid_ub(2)) * a_np(i+1, grid_ub(2)) - &
              chi(i-1, grid_ub(2)) * a_np(i-1, grid_ub(2)) &
            ) + &
            alpha * a_n(i, grid_ub(2)) + &
            ddy * ( &
              a_n(i, grid_ub(2)-1) - a_n(i, grid_ub(2)) &
            ) + &
            beta * dxi * ( &
              2.0_p_k_fld * ddy * ( &
                a_np(i+1, grid_ub(2)-1) - a_np(i+1, grid_ub(2)) &
              ) - &
              2.0_p_k_fld * ddy * ( &
                a_np(i-1, grid_ub(2)-1) - a_np(i-1, grid_ub(2)) &
              ) &
            )
        enddo
        !$omp end parallel do
      endif
    endif

  case(p_cylindrical_b)
    ! grid offset for laplace
    if( .not. root(this%disc_mat(2)%nconf)) then
      grid_offset = this%disc_mat(2)%id_min - 1
    else
      grid_offset = 0
    endif

    ! inner grid
    !$omp parallel do &
    !$omp&  private(i) &
    !$omp&  shared(grid_lb, grid_ub, a_np, a_n, a_nm, src, chi) &
    !$omp&  shared(alpha, beta, dxi, ddy, grid_offset)
    do j = grid_lb(2), grid_ub(2)
      do i = grid_lb(1), grid_ub(1)
        src(i, j) = &
          chi(i, j) * a_np(i, j) - &
          beta * dxi * ( chi(i+1, j) * a_np(i+1, j) - chi(i-1, j) * a_np(i-1, j) ) + &
          alpha * a_n(i, j) + &
          0.5_p_k_fld * ddy * ( &
            (1.0_p_k_fld - 0.5_p_k_fld / ((j + grid_offset) - 1.5)) * a_n(i, j-1) - &
            2.0_p_k_fld * a_n(i, j) + &
            (1.0_p_k_fld + 0.5_p_k_fld / ((j + grid_offset) - 1.5)) * a_n(i, j+1) &
          ) + &
          beta * dxi * ddy * ( &
            ( &
              (1.0_p_k_fld - 0.5_p_k_fld / ((j + grid_offset) - 1.5)) * a_np(i+1, j-1) - &
              2.0_p_k_fld * a_np(i+1, j) + &
              (1.0_p_k_fld + 0.5_p_k_fld / ((j + grid_offset) - 1.5)) * a_np(i+1, j+1) &
            ) - &
            ( &
              (1.0_p_k_fld - 0.5_p_k_fld / ((j + grid_offset) - 1.5)) * a_np(i-1, j-1) - &
              2.0_p_k_fld * a_np(i-1, j) + &
              (1.0_p_k_fld + 0.5_p_k_fld / ((j + grid_offset) - 1.5)) * a_np(i-1, j+1) &
            ) &
         )
      enddo
    enddo
    !$omp end parallel do

    ! lower boundary
    if( root(this%disc_mat(2)%nconf) ) then
      !$omp parallel do &
      !$omp&  shared(grid_lb, grid_ub, a_np, a_n, a_nm, src, chi) &
      !$omp&  shared(alpha, beta, dxi, ddy)
      do i = grid_lb(1), grid_ub(1)
        src(i, grid_lb(2)) = &
          chi(i, grid_lb(2)) * a_np(i, grid_lb(2)) - &
          beta * dxi * ( &
            chi(i+1, grid_lb(2)) * a_np(i+1, grid_lb(2)) - &
            chi(i-1, grid_lb(2)) * a_np(i-1, grid_lb(2)) &
          ) + &
          alpha * a_n(i, 2) + &
          ddy * (a_n(i, 3) - a_n(i, 2)) + &
          beta * dxi * ddy * 2.0_p_k_fld * ( &
            (a_np(i+1, 2) - a_np(i+1, 1)) - &
            (a_np(i-1, 2) - a_np(i-1, 1))   &
          )
      enddo
      !$omp end parallel do
    endif
    ! upper boundary
    if( last(this%disc_mat(2)%nconf) ) then
      !$omp parallel do &
      !$omp&  shared(grid_lb, grid_ub, a_np, a_n, a_nm, src, chi) &
      !$omp&  shared(alpha, beta, dxi, ddy)
      do i = grid_lb(1), grid_ub(1)
        src(i, grid_ub(2)) = &
          chi(i, grid_ub(2)) * a_np(i, grid_ub(2)) - &
          beta * dxi * ( &
            chi(i+1, grid_ub(2)) * a_np(i+1, grid_ub(2)) - &
            chi(i-1, grid_ub(2)) * a_np(i-1, grid_ub(2))   &
          ) + &
          alpha * a_n(i, grid_ub(2)) + &
          ddy * (a_n(i, grid_ub(2)-1) - a_n(i, grid_ub(2))) + &
          beta * dxi * ddy * 2.0_p_k_fld * ( &
            (a_np(i+1, grid_ub(2)-1) - a_np(i+1, grid_ub(2))) - &
            (a_np(i-1, grid_ub(2)-1) - a_np(i-1, grid_ub(2)))   &
          )
      enddo
      !$omp end parallel do
    endif
  end select

  nullify( chi, src, a_n, a_np )

end subroutine build_source_2d
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! Advancing interface for envelope in 3d in x2 direction
!   - advances the envelope in x2 direction
!   - ADI method allows to use same methods as for 2d advancing by doing two
!     individual steps
!   - passes solution to corresponding time fields in envelope
!-------------------------------------------------------------------------------
subroutine advance_env_3d_x2( this )

  implicit none

  ! input/output parameters
  class( t_emf_pgc ), intent(inout)  ::  this

  ! auxiliary variables
  integer, dimension(3) :: lb, ub
  integer :: i, j, k, nx_min, ny_min, nx_max, ny_max, nz_min, nz_max
  complex( p_k_fld ), pointer :: src(:,:,:)
  complex( p_k_fld ), pointer :: a_nm(:,:,:), a_n(:,:,:), a_np(:,:,:)

  ! extract informations by accessing the tridiagonal system informations
  nx_min = 1
  nx_max = this%env%nx_(1)

  ny_min = this%disc_mat(2)%ext_min
  ny_max = this%disc_mat(2)%ext_max

  nz_min = this%disc_mat(3)%ext_min
  nz_max = this%disc_mat(3)%ext_max

  lb(1) = this%source%lbound(dim=2)
  lb(2) = this%source%lbound(dim=3)
  lb(3) = this%source%lbound(dim=4)

  ub(1) = this%source%ubound(dim=2)
  ub(2) = this%source%ubound(dim=3)
  ub(3) = this%source%ubound(dim=4)

  a_nm(lb(1):, lb(2):, lb(3):) => this%env%z3(p_env_nm, :, :, :)
  a_n(lb(1):, lb(2):, lb(3):)  => this%env%z3(p_env_n, :, :, :)
  a_np(lb(1):, lb(2):, lb(3):) => this%env%z3(p_env_np, :, :, :)
  src(lb(1):, lb(2):, lb(3):)  => this%source%z3(1, :, :, :)

  do k = nz_min, nz_max
    call this%disc_mat(2)%solve(src(nx_min:nx_max, ny_min:ny_max, k))
  enddo

  ! envelope update - for further calculations
  !$omp parallel do &
  !$omp&  private(i, j) &
  !$omp&  shared(lb, ub, a_nm, a_n, a_np, src)
  do k = lb(3), ub(3)
    do j = lb(2), ub(2)
      do i = lb(1), ub(1)
        a_nm(i, j, k) = a_np(i, j, k)
        a_n(i, j, k)  = a_np(i, j, k)
        a_np(i, j, k) = src(i, j, k)
      enddo
    enddo
  enddo
  !$omp end parallel do

end subroutine advance_env_3d_x2
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! Construction of right hand side for 3d advancing of envelope in x2 direction
!-------------------------------------------------------------------------------
subroutine build_source_3d_x2( this , dx , dt , if_periodic )

  use m_tridiag, only : root, last

  implicit none

  class( t_emf_pgc ), intent(inout), target  ::  this
  real( p_double ), dimension(:), intent(in) :: dx
  real( p_double ), intent(in)  ::  dt
  logical, intent(in) :: if_periodic

  complex(p_double) :: alpha, beta
  integer :: i, j, k
  integer, dimension(3) :: lb, grid_lb, grid_ub
  real(p_double) :: dxi, ddy, ddz
  ! pointer helper
  real( p_k_fld ),    pointer :: chi(:,:,:)
  complex( p_k_fld ), pointer :: src(:,:,:), a_n(:,:,:), a_np(:,:,:)

  ! finite difference coefficients
  alpha = cmplx(0.0_p_double, 2.0_p_double, kind=p_double) * this%omega / dt
  beta = cmplx(0.0_p_double, 1.0_p_double, kind=p_double) / this%omega

  ! finite difference coefficient
  dxi = 0.5_p_k_fld / dx(1)             ! for CDS - use 2 * \delta\xi
  ddy = 1.0_p_k_fld / (dx(2) * dx(2))
  ddz = 1.0_p_k_fld / (dx(3) * dx(3))

  ! required for pointer arithmetics
  lb(1) = this%source%lbound(dim=2)
  lb(2) = this%source%lbound(dim=3)
  lb(3) = this%source%lbound(dim=4)

  ! physical grid
  ! - in x1 direction, use local grid and corresponding boundaries
  ! - in x2 direction, use tridiagonal system grid
  ! - in x3 direction, use local grid and corresponding boundaries
  grid_lb(1) = 1
  grid_lb(2) = this%disc_mat(2)%ext_min
  grid_lb(3) = 1

  grid_ub(1) = this%source%nx_(1)
  grid_ub(2) = this%disc_mat(2)%ext_max
  grid_ub(3) = this%source%nx_(3)

  ! pointer for better reading and for openMP
  chi(lb(1):, lb(2):, lb(3):)  => this%chi(1)%f3(1, :, :, :)
  src(lb(1):, lb(2):, lb(3):)  => this%source%z3(1, :, :, :)
  a_n(lb(1):, lb(2):, lb(3):)  => this%env%z3(p_env_n, :, :, :)
  a_np(lb(1):, lb(2):, lb(3):) => this%env%z3(p_env_np, :, :, :)

  ! inner grid
  ! - boundaries for guard cells in transversal direction are set to zero
  !$omp parallel do &
  !$omp&  private(i, j) &
  !$omp&  shared(grid_lb, grid_ub, a_np, a_n, src, chi) &
  !$omp&  shared(alpha, beta, dxi, ddy, ddz)
  do k = grid_lb(3), grid_ub(3)
    do j = grid_lb(2), grid_ub(2)
      do i = grid_lb(1), grid_ub(1)
        src(i, j, k) = &
          chi(i, j, k) * a_np(i, j, k) - &
          beta * dxi * ( chi(i+1, j, k) * a_np(i+1, j, k) - chi(i-1, j, k) * a_np(i-1, j, k) ) + &
          alpha * a_n(i, j, k) + &
          0.5_p_k_fld * ddy * ( a_n(i, j+1, k) - 2.0_p_k_fld * a_n(i, j, k) + a_n(i, j-1, k) ) + &
          ddz * ( a_np(i, j, k+1) - 2.0_p_k_fld * a_np(i, j, k) + a_np(i, j, k-1) ) - &
          beta * dxi * ( &
            ddy * ( a_np(i+1, j+1, k) - 2.0_p_k_fld * a_np(i+1, j, k) + a_np(i+1, j-1, k) ) - &
            ddy * ( a_np(i-1, j+1, k) - 2.0_p_k_fld * a_np(i-1, j, k) + a_np(i-1, j-1, k) ) + &
            ddz * ( a_np(i+1, j, k+1) - 2.0_p_k_fld * a_np(i+1, j, k) + a_np(i+1, j, k-1) ) - &
            ddz * ( a_np(i-1, j, k+1) - 2.0_p_k_fld * a_np(i-1, j, k) + a_np(i-1, j, k-1) ) &
          )
      enddo
    enddo
  enddo
  !$omp end parallel do

  if( .not. if_periodic ) then
    ! lower boundary
    if( root(this%disc_mat(2)%nconf) ) then
      !$omp parallel do &
      !$omp&  private(i) &
      !$omp&  shared(grid_lb, grid_ub, a_np, a_n, src, chi) &
      !$omp&  shared(alpha, beta, dxi, ddy, ddz)
      do k = grid_lb(3), grid_ub(3)
        do i = grid_lb(1), grid_ub(1)
          src(i, grid_lb(2), k) = &
            chi(i, grid_lb(2), k) * a_np(i, grid_lb(2), k) - &
            beta * dxi * ( &
              chi(i+1, grid_lb(2), k) * a_np(i+1, grid_lb(2), k) - &
              chi(i-1, grid_lb(2), k) * a_np(i-1, grid_lb(2), k) &
            ) + &
            alpha * a_n(i, grid_lb(2), k) + &
            ddy * ( a_n(i, grid_lb(2)+1, k) - a_n(i, grid_lb(2), k) ) + &
            ddz * ( &
              a_np(i, grid_lb(2), k+1) - &
              2.0_p_k_fld * a_np(i, grid_lb(2), k) + &
              a_np(i, grid_lb(2), k-1) &
            ) - &
            beta * dxi * ( &
              2.0_p_k_fld * ddy * ( &
                a_np(i+1, grid_lb(2)+1, k) - a_np(i+1, grid_lb(2), k) &
              ) - &
              2.0_p_k_fld * ddy * ( &
                a_np(i-1, grid_lb(2)+1, k) - a_np(i-1, grid_lb(2), k) &
              ) + &
              ddz * ( &
                a_np(i+1, grid_lb(2), k+1) - &
                2.0_p_k_fld * a_np(i+1, grid_lb(2), k) + &
                a_np(i+1, grid_lb(2), k-1) &
              ) - &
              ddz * ( &
                a_np(i-1, grid_lb(2), k+1) - &
                2.0_p_k_fld * a_np(i-1, grid_lb(2), k) + &
                a_np(i-1, grid_lb(2), k-1) &
              ) &
          )
        enddo
      enddo
      !$omp end parallel do
    endif

    ! upper boundary
    if( last(this%disc_mat(2)%nconf) ) then
      !$omp parallel do &
      !$omp&  private(i) &
      !$omp&  shared(grid_lb, grid_ub, a_np, a_n, src, chi) &
      !$omp&  shared(alpha, beta, dxi, ddy, ddz)
      do k = grid_lb(3), grid_ub(3)
        do i = grid_lb(1), grid_ub(1)
          src(i, grid_ub(2), k) = &
            chi(i, grid_ub(2), k) * a_np(i, grid_ub(2), k) - &
            beta * dxi * ( &
              chi(i+1, grid_ub(2), k) * a_np(i+1, grid_ub(2), k) - &
              chi(i-1, grid_ub(2), k) * a_np(i-1, grid_ub(2), k) &
            ) + &
            alpha * a_n(i, grid_ub(2), k) + &
            ddy * ( a_n(i, grid_ub(2)-1, k) - a_n(i, grid_ub(2), k) ) + &
            ddz * ( &
              a_np(i, grid_ub(2), k+1) - &
              2.0_p_k_fld * a_np(i, grid_ub(2), k) + &
              a_np(i, grid_ub(2), k-1) &
            ) - &
            beta * dxi * ( &
              2.0_p_k_fld * ddy * ( &
                a_np(i+1, grid_ub(2)-1, k) - a_np(i+1, grid_ub(2), k) &
              ) - &
              2.0_p_k_fld * ddy * ( &
                a_np(i-1, grid_ub(2)-1, k) - a_np(i-1, grid_ub(2), k) &
              ) + &
              ddz * ( &
                a_np(i+1, grid_ub(2), k+1) - &
                2.0_p_k_fld * a_np(i+1, grid_ub(2), k) + &
                a_np(i+1, grid_ub(2), k-1) &
              ) - &
              ddz * ( &
                a_np(i-1, grid_ub(2), k+1) - &
                2.0_p_k_fld * a_np(i-1, grid_ub(2), k) + &
                a_np(i-1, grid_ub(2), k-1) &
              ) &
          )
        enddo
      enddo
      !$omp end parallel do
    endif
  endif

  nullify( chi, src, a_n, a_np )

end subroutine build_source_3d_x2
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! Advancing interface for envelope in 3d in x3 direction
!   - advances the envelope in x3 direction
!   - ADI method allows to use same methods as for 2d advancing by doing two
!     individual steps
!   - passes solution to corresponding time fields in envelope
!-------------------------------------------------------------------------------
subroutine advance_env_3d_x3( this )

  implicit none

  ! input/output parameters
  class( t_emf_pgc ), intent(inout)  ::  this

  ! auxiliary variables
  integer, dimension(3) :: lb, ub
  integer :: i, j, k, nx_min, ny_min, nx_max, ny_max, nz_min, nz_max
  complex( p_k_fld ), pointer :: src(:,:,:)
  complex( p_k_fld ), pointer :: a_nm(:,:,:), a_n(:,:,:), a_np(:,:,:)

  ! extract informations by accessing the tridiagonal system informations
  nx_min = 1
  nx_max = this%env%nx_(1)

  ny_min = this%disc_mat(2)%ext_min
  ny_max = this%disc_mat(2)%ext_max

  nz_min = this%disc_mat(3)%ext_min
  nz_max = this%disc_mat(3)%ext_max

  lb(1) = this%source%lbound(dim=2)
  lb(2) = this%source%lbound(dim=3)
  lb(3) = this%source%lbound(dim=4)

  ub(1) = this%source%ubound(dim=2)
  ub(2) = this%source%ubound(dim=3)
  ub(3) = this%source%ubound(dim=4)

  a_nm(lb(1):, lb(2):, lb(3):) => this%env%z3(p_env_nm, :, :, :)
  a_n(lb(1):, lb(2):, lb(3):)  => this%env%z3(p_env_n, :, :, :)
  a_np(lb(1):, lb(2):, lb(3):) => this%env%z3(p_env_np, :, :, :)
  src(lb(1):, lb(2):, lb(3):)  => this%source%z3(1, :, :, :)

  do j = ny_min, ny_max
    call this%disc_mat(3)%solve(src(nx_min:nx_max, j, nz_min:nz_max))
  enddo

  ! update envelope to match the advancing
  !$omp parallel do &
  !$omp&  private(i, j) &
  !$omp&  shared(lb, ub, a_nm, a_n, a_np, src)
  do k = lb(3), ub(3)
    do j = lb(2), ub(2)
      do i = lb(1), ub(1)
        a_np(i, j, k) = src(i, j, k)
        a_n(i,j,k) = 0.5_p_k_fld * (a_nm(i, j, k) + a_np(i, j, k))
      enddo
    enddo
  enddo
  !$omp end parallel do

end subroutine advance_env_3d_x3
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! Construction of right hand side for 3d advancing of envelope in x2 direction
!-------------------------------------------------------------------------------
subroutine build_source_3d_x3( this , dx , dt , if_periodic )

  use m_tridiag, only : root, last

  implicit none

  class( t_emf_pgc ), intent(inout), target  ::  this
  real( p_double ), dimension(:), intent(in) :: dx
  real( p_double ), intent(in)  ::  dt
  logical, intent(in) :: if_periodic

  complex(p_double) :: alpha, beta
  integer :: i, j, k
  integer, dimension(3) :: lb, grid_lb, grid_ub
  real(p_double) :: dxi, ddy, ddz
  ! pointer helper
  real( p_k_fld ),    pointer :: chi(:,:,:)
  complex( p_k_fld ), pointer :: src(:,:,:), a_n(:,:,:), a_np(:,:,:)

  ! finite difference coefficients
  alpha = cmplx(0.0_p_double, 2.0_p_double, kind=p_double) * this%omega / dt
  beta = cmplx(0.0_p_double, 1.0_p_double, kind=p_double) / this%omega

  ! finite difference coefficient
  dxi = 0.5_p_k_fld / dx(1)             ! for CDS - use 2 * \delta\xi
  ddy = 1.0_p_k_fld / (dx(2) * dx(2))
  ddz = 1.0_p_k_fld / (dx(3) * dx(3))

  ! required for pointer arithmetics
  lb(1) = this%source%lbound(dim=2)
  lb(2) = this%source%lbound(dim=3)
  lb(3) = this%source%lbound(dim=4)

  ! physical grid
  ! - in x1 direction, use local grid and corresponding boundaries
  ! - in x2 direction, use tridiagonal system grid
  ! - in x3 direction, use local grid and corresponding boundaries
  grid_lb(1) = 1
  grid_lb(2) = this%disc_mat(2)%ext_min
  grid_lb(3) = 1

  grid_ub(1) = this%source%nx_(1)
  grid_ub(2) = this%disc_mat(2)%ext_max
  grid_ub(3) = this%source%nx_(3)

  ! pointer for better reading and for openMP
  chi(lb(1):, lb(2):, lb(3):)  => this%chi(1)%f3(1, :, :, :)
  src(lb(1):, lb(2):, lb(3):)  => this%source%z3(1, :, :, :)
  a_n(lb(1):, lb(2):, lb(3):)  => this%env%z3(p_env_n, :, :, :)
  a_np(lb(1):, lb(2):, lb(3):) => this%env%z3(p_env_np, :, :, :)

  ! inner grid
  ! - boundaries for guard cells in transversal direction are set to zero
  !$omp parallel do &
  !$omp&  private(i, j) &
  !$omp&  shared(grid_lb, grid_ub, a_np, a_n, src, chi) &
  !$omp&  shared(alpha, beta, dxi, ddy, ddz)
  do k = grid_lb(3), grid_ub(3)
    do j = grid_lb(2), grid_ub(2)
      do i = grid_lb(1), grid_ub(1)
        src(i, j, k) = &
          chi(i, j, k) * a_np(i, j, k) - &
          beta * dxi * ( chi(i+1, j, k) * a_np(i+1, j, k) - chi(i-1, j, k) * a_np(i-1, j, k) ) + &
          alpha * a_n(i, j, k) + &
          0.5_p_k_fld * ddz * ( a_n(i, j, k+1) - 2.0_p_k_fld * a_n(i, j, k) + a_n(i, j, k-1) ) + &
          ddy * ( a_np(i, j+1, k) - 2.0_p_k_fld * a_np(i, j, k) + a_np(i, j-1, k) ) - &
          beta * dxi * ( &
            ddy * ( a_np(i+1, j+1, k) - 2.0_p_k_fld * a_np(i+1, j, k) + a_np(i+1, j-1, k) ) - &
            ddy * ( a_np(i-1, j+1, k) - 2.0_p_k_fld * a_np(i-1, j, k) + a_np(i-1, j-1, k) ) + &
            ddz * ( a_np(i+1, j, k+1) - 2.0_p_k_fld * a_np(i+1, j, k) + a_np(i+1, j, k-1) ) - &
            ddz * ( a_np(i-1, j, k+1) - 2.0_p_k_fld * a_np(i-1, j, k) + a_np(i-1, j, k-1) ) &
          )
      enddo
    enddo
  enddo
  !$omp end parallel do

  if( .not. if_periodic ) then
    ! lower boundary
    if( root(this%disc_mat(3)%nconf) ) then
      !$omp parallel do &
      !$omp&  private(i) &
      !$omp&  shared(grid_lb, grid_ub, a_np, a_n, src, chi) &
      !$omp&  shared(alpha, beta, dxi, ddy, ddz)
      do j = grid_lb(2), grid_ub(2)
        do i = grid_lb(1), grid_ub(1)
          src(i, j, grid_lb(3)) = &
            chi(i, j, grid_lb(3)) * a_np(i, j, grid_lb(3)) - &
            beta * dxi * ( &
              chi(i+1, j, grid_lb(3)) * a_np(i+1, j, grid_lb(3)) - &
              chi(i-1, j, grid_lb(3)) * a_np(i-1, j, grid_lb(3)) &
            ) + &
            alpha * a_n(i, j, grid_lb(3)) + &
            ddz * ( a_n(i, j, grid_lb(3)+1) - a_n(i, j, grid_lb(3)) ) + &
            ddy * ( &
              a_np(i, j+1, grid_lb(3)) - &
              2.0_p_k_fld * a_np(i, j, grid_lb(3)) + &
              a_np(i, j-1, grid_lb(3)) &
            ) - &
            beta * dxi * ( &
              ddy * ( &
                a_np(i+1, j+1, grid_lb(3)) - &
                2.0_p_k_fld * a_np(i+1, j, grid_lb(3)) + &
                a_np(i+1, j-1, grid_lb(3)) &
              ) - &
              ddy * ( &
                a_np(i-1, j+1, grid_lb(3)) - &
                2.0_p_k_fld * a_np(i-1, j, grid_lb(3)) + &
                a_np(i-1, j-1, grid_lb(3)) &
              ) + &
              2.0_p_k_fld * ddz * ( &
                a_np(i+1, j, grid_lb(3)+1) - a_np(i+1, j, grid_lb(3)) &
              ) - &
              2.0_p_k_fld * ddz * ( &
                a_np(i-1, j, grid_lb(3)+1) - a_np(i-1, j, grid_lb(3)) &
            ) &
          )
        enddo
      enddo
      !$omp end parallel do

    endif

    ! upper boundary
    if( last(this%disc_mat(3)%nconf) ) then
      !$omp parallel do &
      !$omp&  private(i) &
      !$omp&  shared(grid_lb, grid_ub, a_np, a_n, src, chi) &
      !$omp&  shared(alpha, beta, dxi, ddy, ddz)
      do j = grid_lb(2), grid_ub(2)
        do i = grid_lb(1), grid_ub(1)
          src(i, j, grid_ub(3)) = &
            chi(i, j, grid_ub(3)) * a_np(i, j, grid_ub(3)) - &
            beta * dxi * ( &
              chi(i+1, j, grid_ub(3)) * a_np(i+1, j, grid_ub(3)) - &
              chi(i-1, j, grid_ub(3)) * a_np(i-1, j, grid_ub(3)) &
            ) + &
            alpha * a_n(i, j, grid_ub(3)) + &
            2.0_p_k_fld * ddz * ( &
              a_n(i, j, grid_ub(3)-1) - a_n(i, j, grid_ub(3)) &
            ) + &
            ddy * ( &
              a_np(i, j+1, grid_ub(3)) - &
              2.0_p_k_fld * a_np(i, j, grid_ub(3)) + &
              a_np(i, j-1, grid_ub(3)) &
            ) - &
            beta * dxi * ( &
              ddy * ( &
                a_np(i+1, j+1, grid_ub(3)) - &
                2.0_p_k_fld * a_np(i+1, j, grid_ub(3)) + &
                a_np(i+1, j-1, grid_ub(3)) &
              ) - &
              ddy * ( &
                a_np(i-1, j+1, grid_ub(3)) - &
                2.0_p_k_fld * a_np(i-1, j, grid_ub(3)) + &
                a_np(i-1, j-1, grid_ub(3)) &
              ) + &
              2.0_p_k_fld * ddz * ( &
                a_np(i+1, j, grid_ub(3)-1) - a_np(i+1, j, grid_ub(3)) &
              ) - &
              2.0_p_k_fld * ddz * ( &
                a_np(i-1, j, grid_ub(3)-1) - a_np(i-1, j, grid_ub(3)) &
              ) &
            )
        enddo
      enddo
      !$omp end parallel do
    endif
  endif

  nullify( chi, src, a_n, a_np )

end subroutine build_source_3d_x3
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!       read object information from a restart file
!-------------------------------------------------------------------------------
subroutine restart_read_pgc( this, restart_handle )

  use m_restart

  implicit none

  class( t_emf_pgc ), intent(inout) :: this
  type( t_restart_handle ), intent(in) :: restart_handle

  character(len=len(p_emf_pgc_rst_id)) :: rst_id
  integer :: ierr

  restart_io_rd( rst_id, restart_handle, ierr )
  if ( ierr/=0 ) then
  ERROR('error reading restart data for emf object.')
  call abort_program(p_err_rstrd)
  endif

  ! check if restart file is compatible
  if ( rst_id /= p_emf_pgc_rst_id) then
    ERROR('Corrupted restart file, or restart file')
    ERROR('from incompatible binary (emf)')
    call abort_program(p_err_rstrd)
  endif

  call this%env%read_checkpoint( restart_handle )

  ! Write restart data for chi
  call this%chi(1)%read_checkpoint( restart_handle )

end subroutine restart_read_pgc
!-------------------------------------------------------------------------------
end module m_emf_pgc
!-------------------------------------------------------------------------------
!===============================================================================
! t_emf_pgc type bound procedures
!===============================================================================
!-------------------------------------------------------------------------------
! Read input file data for PGC algorithm
!   - this routine also calls the superclass (t_emf) read_nml routine to read
!     that section first
!-------------------------------------------------------------------------------
subroutine read_input_pgc( this, input_file, periodic, if_move, grid, dx, dt, gamma )
  use m_emf_pgc_define
  use m_emf_pgc
  use m_tridiag, only : p_lapack_alg, p_lu_alg, p_thomas_alg
  use m_input_file, only : t_input_file, get_namelist, disp_out
  use m_grid_define, only : t_grid
  use m_vdf_smooth, only : read_nml

  implicit none

  class( t_emf_pgc ), intent( inout )  :: this
  class( t_input_file ), intent(inout)  :: input_file
  logical, dimension(:), intent(in)    :: periodic, if_move
  class( t_grid ), intent(in)           :: grid
  real(p_double), dimension(:), intent(in) :: dx
  real(p_double), intent(in) :: dt
  real(p_double), intent(in) :: gamma 

  ! namelist parameters
  real(p_k_fld)               :: a0, omega
  logical                     :: free_stream
  character(len=16)           :: per_type, lon_type
  real(p_k_fld)               :: lon_center, lon_duration, lon_range
  real(p_k_fld)               :: lon_rise, lon_fall, lon_flat, lon_start
  real(p_k_fld)               :: w0, per_focus
  real(p_k_fld), dimension(2) :: per_center, w0_asym, per_focus_asym
  character(len=16)           :: solver_type
  logical                     :: error_cleaning
  character(len=16)           :: source_smoothing

  logical :: chi_smooth, fp_smooth, env_smooth
  integer :: fp_smooth_niter, fp_smooth_nmax

  integer :: ierr

  namelist /nl_pgc/ a0, omega, free_stream, per_type, lon_type, &
                    lon_center, lon_duration, lon_range, &
                    lon_rise, lon_fall, lon_flat, lon_start, &
                    w0, per_focus, per_center, &
                    w0_asym, per_focus_asym, &
                    chi_smooth, env_smooth, &
                    fp_smooth, fp_smooth_niter, fp_smooth_nmax, &
                    solver_type, error_cleaning, source_smoothing

  ! -------------- Allocate any sub-objects needed --------------------
  call this%allocate_objs()

  ! first read the superclass (emf) section
  call this%t_emf%read_input(input_file, periodic, if_move, grid, dx, dt, gamma)

  ! message for processing PGC section
  if (disp_out(input_file)) then
    SCR_ROOT("Reading pgc configuration...")
  endif

  !> default values <
  ! general pgc parameters
  a0             = 0.0_p_k_fld
  omega          = 10.0_p_k_fld
  free_stream    = .false.
  per_type       = "gaussian"
  lon_type       = "sin2"

  ! longitudinal parameters
  lon_center     = 0.0_p_k_fld
  lon_start      = 0.0_p_k_fld
  lon_duration   = 2.0_p_k_fld
  lon_range      = 4.0_p_k_fld
  lon_rise       = 1.0_p_k_fld
  lon_fall       = 1.0_p_k_fld
  lon_flat       = 0.0_p_k_fld

  ! parameters symmetric gaussian beam
  per_center     = 0.0_p_k_fld
  w0             = 1.0_p_k_fld
  per_focus      = 0.0_p_k_fld
  w0_asym        = 1.0_p_k_fld
  per_focus_asym = 0.0_p_k_fld

  ! chi smoothing default parameters
  chi_smooth = .false.

  ! fp smoothing default parameters
  fp_smooth       = .false.
  fp_smooth_niter = 1
  fp_smooth_nmax  = -1

  ! chi smoothing default parameters
  env_smooth = .false.

  ! solver type for the advancing
  solver_type = "thomas"
  ! error cleaning and source term smoothing for tridiagonal system
  error_cleaning   = .false.
  source_smoothing = "none"

  ! read values from file
  call get_namelist( input_file, "nl_pgc", ierr )

  if (ierr == 0) then
    read (input_file%nml_text, nml = nl_pgc, iostat = ierr)
    if (ierr /= 0) then
      print *, "Error reading pgc parameters"
      print *, "aborting..."
      stop
    endif
  else
    SCR_ROOT(" - no pgc parameters specified")
    stop
  endif

  ! store namelist parameters
  this%a0          = a0
  this%omega       = omega
  this%free_stream = free_stream

  select case ( trim(lon_type) )
  case ("gaussian")
    this%lon_type     = p_pgc_lon_gaussian

    this%lon_center   = lon_center
    this%lon_duration = lon_duration
    this%lon_range    = lon_range

  case ("polynomial")
    this%lon_type    = p_pgc_lon_polynomial

    this%lon_start   = lon_start
    this%lon_rise    = lon_rise
    this%lon_fall    = lon_fall
    this%lon_flat    = lon_flat

  case ("sin2")
    this%lon_type    = p_pgc_lon_sin2

    this%lon_rise    = lon_rise
    this%lon_fall    = lon_fall
    this%lon_flat    = lon_flat
    this%lon_center  = lon_center

  case default
    SCR_ROOT("")
    SCR_ROOT("Error in pgc parameters :: lon_type")
    SCR_ROOT("Was not able to find a valid longitudinal profile")
    SCR_ROOT("aborting ...")
    stop
  end select

  select case ( trim(per_type) )
  case ("plane")
    this%per_type    = p_pgc_per_plane

  case ("gaussian")
    this%per_type    = p_pgc_per_gaussian

    this%w0          = w0
    this%per_focus   = per_focus
    this%per_center  = per_center

  case ("asymmetric")
    this%per_type       = p_pgc_per_gaussian_asym

    this%w0_asym        = w0_asym
    this%per_focus_asym = per_focus_asym

  case default
    SCR_ROOT("")
    SCR_ROOT("Error in pgc parameters :: per_type")
    SCR_ROOT("Was not able to find a valid perpendicular profile")
    SCR_ROOT("aborting ...")
    SCR_ROOT("")
    stop
  end select

  ! checks for smothing flags
  if( chi_smooth ) then
    this%if_smooth_chi = .true.
  else
    this%if_smooth_chi = .false.
  endif

  if( fp_smooth ) then
    this%if_smooth_fp = .true.

    if (fp_smooth_niter <= 0) then
      SCR_ROOT("")
      SCR_ROOT("Error in pgc parameters :: fp_smooth_niter")
      SCR_ROOT("'fp_smooth_niter' must be defined and > 0")
      SCR_ROOT("aborting ...")
      SCR_ROOT("")
      stop
    else
      this%fp_smooth_niter = fp_smooth_niter
      this%fp_smooth_nmax = fp_smooth_nmax
    endif
  else
    this%if_smooth_fp = .false.
  endif

  if( env_smooth ) then
    this%if_smooth_env = .true.
  else
    this%if_smooth_env = .false.
  endif

  call read_nml( this%chi_smooth, input_file, nml_name="chi_smooth" )
  call read_nml( this%fp_smooth, input_file, nml_name="fp_smooth" )
  call read_nml( this%env_smooth, input_file, nml_name="env_smooth" )

  select case ( trim(solver_type) )
  case("lapack")
    if (disp_out(input_file)) then
      SCR_ROOT(" - using LAPACK for tridiagonal solver")
    endif
    this%solver_type = p_lapack_alg
  case("LU")
    if (disp_out(input_file)) then
      SCR_ROOT(" - using custom LU decomposition for tridiagonal solver")
    endif
    this%solver_type = p_lu_alg
  case("thomas")
    if (disp_out(input_file)) then
      SCR_ROOT(" - using Thomas algorithm for tridiagonal solver")
    endif
    this%solver_type = p_thomas_alg
  case default
    SCR_ROOT("")
    SCR_ROOT("Error in pgc parameters :: solver_type")
    SCR_ROOT("solver_type should have the value")
    SCR_ROOT("'thomas', 'LU' and 'lapack'")
    SCR_ROOT("'lapack' parameter requires compilation with LAPACK")
    SCR_ROOT("aborting ...")
    SCR_ROOT("")
    stop
  end select

  this%error_cleaning   = error_cleaning
  select case ( trim(source_smoothing) )
  case("both")
    this%pre_calc_smoothing = .true.
    this%post_calc_smoothing = .true.
    if (disp_out(input_file)) then
      SCR_ROOT(" - smooths right-hand before and after advancing")
    endif
  case("before")
    this%pre_calc_smoothing = .true.
    this%post_calc_smoothing = .false.
    if (disp_out(input_file)) then
      SCR_ROOT(" - smooths right-hand side before advancing")
    endif
  case("after")
    this%pre_calc_smoothing = .false.
    this%post_calc_smoothing = .true.
    if (disp_out(input_file)) then
      SCR_ROOT(" - smooths right-hand side after advancing")
    endif
  case("none")
    this%pre_calc_smoothing = .false.
    this%post_calc_smoothing = .false.
    if (disp_out(input_file)) then
      SCR_ROOT(" - no smoothing of right-hand side for advancing")
    endif
  case default
    SCR_ROOT("")
    SCR_ROOT("Error in pgc parameters :: source_smoothing")
    SCR_ROOT("source smoothing should have the value both/before/after/none")
    SCR_ROOT("aborting ...")
    SCR_ROOT("")
    stop
  end select

end subroutine read_input_pgc
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! Initialize emf_pgc object
!   - this routine also calls the superclass (t_emf) setup routine first to
!     setup the electric magnetic fields
!   - setup smooothing for PGC solver
!   - setup envelope solver quantaties
!   -
!-------------------------------------------------------------------------------
subroutine init_pgc( this, part_grid_center, part_interpolation, &
                     g_space, grid, gc_min, dx, tstep, tmin, tmax, &
                     no_co, send_msg, recv_msg, restart,  restart_handle, sim_options )

  use m_emf_pgc_define
  use m_emf_pgc
  use m_space,       only : t_space, xmin
  use m_grid_define, only : t_grid
  use m_node_conf,   only : t_node_conf, n_threads
  use m_restart,     only : t_restart_handle
  use m_vdf_define,  only : t_vdf
  use m_vdf_memory,  only : alloc
  use m_vdf_smooth,  only : if_smooth, setup
  use m_vdf_comm,    only : update_boundary, t_vdf_msg
  use m_time_step,   only : t_time_step, dt
#ifdef _OPENMP
  use omp_lib
#endif

  implicit none

  ! input parameters
  class( t_emf_pgc ), intent( inout ), target :: this
  logical, intent(in)                         :: part_grid_center
  integer, intent(in)                         :: part_interpolation
  type( t_space ), intent(in)                 :: g_space
  class( t_grid ), intent(in)                 :: grid
  integer, dimension(:,:), intent(in)         :: gc_min
  real( p_double ), dimension(:), intent(in)  :: dx
  type( t_time_step ), intent(in) :: tstep
  real(p_double), intent(in) :: tmin, tmax
  class( t_node_conf ), intent(in), target    :: no_co
  type(t_vdf_msg), dimension(2), intent(inout) :: send_msg, recv_msg
  logical, intent(in)                         :: restart
  type( t_restart_handle ), intent(in)        :: restart_handle
  type( t_options ), intent(in)               :: sim_options
  ! auxiliary variables
  real(p_k_fld), dimension(p_max_dim) :: g_xmin
#ifdef _OPENMP
  integer :: nt, i
#endif

  ! setup superclass data - this will also setup the diagnostics
  call this % t_emf % init( part_grid_center, part_interpolation, &
                            g_space, grid, gc_min, dx, tstep, tmin, tmax, &
                            no_co, send_msg, recv_msg, &
                            restart, restart_handle, sim_options )

  ! Store a pointer to no_co so we can use it in the advance method
  this % no_co => no_co

  ! setup smoothing
  call setup(this%chi_smooth, dx, sim_options%gamma)
  call setup(this%fp_smooth, dx, sim_options%gamma)
  call setup(this%env_smooth, dx, sim_options%gamma)

  !executable statements
  g_xmin(1:p_x_dim) = real(xmin(g_space), p_k_fld)

  ! allocate squared envelope and ponderomotive force for particle advancing
  call this%sqrt_nm%new(this%e , f_dim = 1)
  call this%sqrt_n%new(this%e , f_dim = 1)
  call this%sqrt_np%new(this%e , f_dim = 1)
  call this%env_mod%new(this%e, f_dim = 1)
  call this%f_n%new(this%e)
  call this%f_np%new(this%e)

  ! chi - allocate for each thread for chi deposition
#ifdef _OPENMP
  nt = n_threads( no_co )
  call alloc( this%chi, (/ nt /) )

  !$omp parallel private (i)
  i = omp_get_thread_num() + 1
  !$omp critical
  call this%chi(i)%new( this%e, f_dim = 1 )
  !$omp end critical
  !$omp end parallel
#else
  call alloc( this%chi, (/ 1 /) )
  call this%chi(1)%new( this%e, f_dim = 1 )
#endif

  ! allocate discritization matrix ( 2 -> x2, 3 -> x3)
  allocate(this%disc_mat(2:p_x_dim ))

  if ( restart ) then
    call restart_read_pgc( this, restart_handle )
  else
    ! Initialize `env` arrays - check for `p_env_nm, p_env_n, p_env_np`
    ! for corresponding time step
    call this%env%new( this%e , f_dim = 3 )
    call this%env%zero()
    !zero ponderomotive force case e is not zero at t = 0
    call this%f_n%zero()
    call this%f_np%zero()
    call this%reset_chi()
  endif

  call setup_pgc( this, g_xmin, gc_min, grid%my_nx( p_lower, : ), dx, dt(tstep), &
                  grid, no_co, restart)

  ! update envelope values after initialization
  call update_boundary(this%env, p_vdf_replace, no_co, send_msg, recv_msg)

  ! calculate particle advancing quantaties and communicate them
  call update_pgc_part( this, dx )
  call update_boundary(this%sqrt_nm, p_vdf_replace, no_co, send_msg, recv_msg)
  call update_boundary(this%sqrt_n, p_vdf_replace, no_co, send_msg, recv_msg)
  call update_boundary(this%sqrt_np, p_vdf_replace, no_co, send_msg, recv_msg)
  call update_boundary(this%f_np, p_vdf_replace, no_co, send_msg, recv_msg)
  call update_boundary(this%f_n, p_vdf_replace, no_co, send_msg, recv_msg)

  ! initialize longitudinal boundaries
  call init_lon_bounds( this, no_co, grid )

  ! check if envelope is initialized outside the boundary mask
  if( this%lb_apply ) then
    call check_mask_env_overlap( this, this%bnd_lb_low, this%bnd_lb_upp )
  endif

  if( this%ub_apply ) then
    call check_mask_env_overlap( this, this%bnd_ub_low, this%bnd_ub_upp )
  endif

  if( if_smooth(this%chi_smooth) ) then
    if( .not. this%if_smooth_chi ) then
      print *, "Error in pgc initialization"
      print *, "specified a chi smoothing section but no chi smoothing"
      print *, "aborting..."
      stop
    endif
  else
    if( this%if_smooth_chi ) then
      print *, "Error in pgc initialization"
      print *, "specified chi smoothing but no smoothing section found"
      print *, "aborting..."
      stop
    endif
  endif

  ! check for ponderomotive force smoothing sections
  ! and if ponderomotive force should be smoothed
  if( if_smooth(this%fp_smooth) ) then
    if( .not. this%if_smooth_fp ) then
      print *, "Error in pgc initialization"
      print *, "specified a fp smoothing section but no fp smoothing"
      print *, "aborting..."
      stop
    endif
  else
    if( this%if_smooth_fp ) then
      print *, "Error in pgc initialization"
      print *, "specified fp smoothing but no smoothing section found"
      print *, "aborting..."
      stop
    endif
  endif

  ! check for chi smoothing sections and if chi should be smoothed
  if( if_smooth(this%env_smooth) ) then
    if( .not. this%if_smooth_env ) then
      print *, "Error in pgc initialization"
      print *, "specified a envelope smoothing section but no envelope smoothing"
      print *, "aborting..."
      stop
    endif
  else
    if( this%if_smooth_env ) then
      print *, "Error in pgc initialization"
      print *, "specified envelope smoothing but no smoothing section found"
      print *, "aborting..."
      stop
    endif
  endif

end subroutine init_pgc
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! Advance EMF fields and PGC laser
!-------------------------------------------------------------------------------
subroutine advance_pgc( this, jay, dt, send_msg, recv_msg )

  use m_emf_pgc_define
  use m_emf_pgc
  use m_emf_pgc_diag,   only : pgc_advance_env
  use m_node_conf,      only : t_node_conf, periodic
  use m_vdf_define,     only : t_vdf
  use m_vdf_comm,       only : update_boundary, t_vdf_msg
  use m_current_define, only : t_current
  use m_vdf_smooth,     only : smooth

  implicit none

  class( t_emf_pgc ), intent( inout )  ::  this
  class( t_current ), intent(inout) :: jay
  real(p_double), intent(in) :: dt

!  class( t_node_conf ), intent(in) :: no_co
!  integer , intent(in) :: coordinates
  type(t_vdf_msg), dimension(2), intent(inout) :: send_msg, recv_msg

  !local variables
  real( p_double ), dimension(p_x_dim) :: dx_pgc
  integer, dimension(p_x_dim) :: nx
  logical :: if_periodic(p_x_dim)

  ! Advance the EMF fields
  call this % t_emf % advance(jay, dt, send_msg, recv_msg)

  if_periodic = periodic( this % no_co )

  ! advance the pgc laser
  call begin_event( pgc_advance_env )
  if ( this%free_stream ) then
    SCR_ROOT( "Using free streaming PGC laser.")

  else
    ! ensure chi is fully updated
    call update_boundary(this%chi(1), p_vdf_add, this % no_co, send_msg, recv_msg)

    ! perform smoothing of chi
    if ( this%if_smooth_chi ) then
      call smooth(this%chi(1), this%chi_smooth)
    endif

    dx_pgc = this%e%dx_(1:p_x_dim)
    nx     = this%e%nx_(1:p_x_dim)

    ! normalize Chi if running in cylindrical coordinates
    if( this % coordinates == p_cylindrical_b ) then
      call norm_charge_cyl_pgc(this%chi(1) , 1 , dx_pgc(2))
    endif

    ! select dimension for the solver
    select case(p_x_dim)
    case (2)
      call build_source_2d(this , dx_pgc , dt , this % coordinates, if_periodic(2))
      call advance_env_2d(this)

      if( this % coordinates == p_cylindrical_b ) then
        call update_axial_bnd(this)
      endif
    case (3)
      ! Leapfrog for advancing of envelope (x2 -> x3)
      call build_source_3d_x2(this , dx_pgc , dt, if_periodic(2))
      call advance_env_3d_x2(this)

      ! update envelope values
      call update_boundary(this%env, p_vdf_replace, this % no_co, send_msg, recv_msg)
      call apply_lon_bnds_env(this)

      ! ADI - 2nd in x3-direction
      call build_source_3d_x3(this , dx_pgc , dt, if_periodic(3))
      call advance_env_3d_x3(this)
    end select

    call apply_lon_bnds_env(this)
    call update_boundary(this%env, p_vdf_replace, this % no_co, send_msg, recv_msg)

    ! perform smoothing of envelope
    if (this%if_smooth_env) then
      call smooth(this%env , this%env_smooth)
    endif

    ! calculate requred quantaties for particle pushing
    ! apply longitudinal boundaries on ponderomotive force
    call update_pgc_part(this, dx_pgc)
    call apply_lon_bnds_fp(this)
    ! update quataties for particle advancing
    call update_boundary(this%sqrt_nm, p_vdf_replace, this % no_co, send_msg, recv_msg)
    call update_boundary(this%sqrt_n, p_vdf_replace, this % no_co, send_msg, recv_msg)
    call update_boundary(this%sqrt_np, p_vdf_replace, this % no_co, send_msg, recv_msg)
    call update_boundary(this%f_np, p_vdf_replace, this % no_co, send_msg, recv_msg)
    call update_boundary(this%f_n, p_vdf_replace, this % no_co, send_msg, recv_msg)

  endif

  call end_event(pgc_advance_env)
  ! TODO: remove: this%n_current_emf = this%n_current_emf + 1
  
end subroutine advance_pgc
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!       report on electro-magnetic field - diagnostic
!-------------------------------------------------------------------------------
subroutine report_diag_pgc( this, g_space, grid, no_co, tstep, t, send_msg, recv_msg )


  use m_emf_pgc_diag
  use m_emf_pgc
  use m_emf_pgc_define
  use m_time_step,    only : t_time_step
  use m_space,        only : t_space
  use m_grid_define,  only : t_grid
  use m_node_conf,    only : t_node_conf
  use m_vdf_define,   only : t_vdf_report
  use m_vdf_comm,     only : t_vdf_msg

  use m_vdf_report

  implicit none

  class( t_emf_pgc ), intent(inout) :: this
  type( t_space ), intent(in)       :: g_space
  class( t_grid ), intent(in)       :: grid
  class( t_node_conf ), intent(in)   :: no_co
  type( t_time_step ), intent(in)   :: tstep
  real( p_double ), intent(in)      :: t
  type(t_vdf_msg), dimension(2), intent(inout) :: send_msg, recv_msg

  type( t_vdf_report ), pointer :: rep
  integer :: quant

  ! run superclass diagnostics
  call this % t_emf % report( g_space, grid, no_co, tstep, t, send_msg, recv_msg )

  ! run PGC specific diagnostics
  rep => this%diag%reports
  do
    if ( .not. associated( rep ) ) exit

    if ( if_report( rep, tstep ) ) then

     quant = rep%quant - this % diag % quant_offset( )

     select case ( quant )

      case ( p_a_mod)
      ! PGC laser envelope
      call report_vdf( rep, this%env_mod, 1, g_space, grid, no_co, tstep, t )

      case ( p_fp1 )
      ! PGC ponderomotive force fp1
      call report_vdf( rep, this%f_n, 1, g_space, grid, no_co, tstep, t )

      case ( p_fp2 )
      ! PGC ponderomotive force fp2
      call report_vdf( rep, this%f_n, 2, g_space, grid, no_co, tstep, t )

      case ( p_fp3 )
      ! PGC ponderomotive force fp3
      call report_vdf( rep, this%f_n, 3, g_space, grid, no_co, tstep, t )

      case ( p_chi )
      ! PGC plasma \chi (susceptibility)
      call report_vdf( rep, this%chi(1), 1, g_space, grid, no_co, tstep, t )

      case default
      ! must be a superclass diagnostic, ignore
      continue
    end select

    endif

    rep => rep%next
  enddo

end subroutine report_diag_pgc
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! Cleanup emf_pgc object
!-------------------------------------------------------------------------------
subroutine cleanup_pgc( this )

  use m_emf_pgc_define
  use m_vdf_memory, only : freemem
  use m_vdf_smooth, only : cleanup

  implicit none

  class(t_emf_pgc) , intent(inout) :: this

  integer :: i

  ! cleanup superclass
  call this%t_emf%cleanup()

  ! cleanup local data structures
  call this%env%cleanup()
  call this%sqrt_nm%cleanup()
  call this%sqrt_n%cleanup()
  call this%sqrt_np%cleanup()
  call this%env_mod%cleanup()

  call this%f_np%cleanup()
  call this%f_n%cleanup()

  do i = lbound(this%chi, dim=1), ubound(this%chi, dim=1)
    call this%chi(i)%cleanup()
  enddo
  call freemem(this%chi)

  do i = lbound(this%disc_mat, dim=1), ubound(this%disc_mat, dim=1)
    call this%disc_mat(i)%cleanup()
  enddo
  deallocate(this%disc_mat)

  call this%source%cleanup()

  if(this%lb_apply) call freemem(this%bd_val_lb)
  if(this%ub_apply) call freemem(this%bd_val_ub)

  call cleanup( this%chi_smooth )
  call cleanup( this%fp_smooth )
  call cleanup( this%env_smooth )

end subroutine cleanup_pgc
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! Updates on particle fields - used for applying smoothing on ponderomotive
! force for the particles
!-------------------------------------------------------------------------------
subroutine update_particle_fld_pgc( this, n, t, dt, space, nx_p_min )

  use m_emf_pgc_define
  use m_emf_pgc

  use m_vdf_smooth, only : smooth
  use m_space,      only : t_space

  implicit none

  class( t_emf_pgc ), intent( inout )  ::  this
  integer,     intent(in) :: n
  real( p_double ), intent(in) :: t, dt
  type( t_space ),     intent(in) :: space
  integer, dimension(:), intent(in) :: nx_p_min

  call this%t_emf%update_particle_fld( n, t, dt, space, nx_p_min )

  ! do local smooth of fields if necessary
  if ( this%if_smooth_fp ) then
     ! Stop smoothing after given iteration
     if (this%fp_smooth_nmax > 0 .and. n >= this%fp_smooth_nmax) then
       this%if_smooth_fp = .false.
     else if ( mod( n, this%fp_smooth_niter ) == 0) then

       call smooth( this%f_np, this%fp_smooth )
       call smooth( this%f_n , this%fp_smooth )

     endif
  endif

end subroutine update_particle_fld_pgc
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!       write object information into a restart file
!-------------------------------------------------------------------------------
subroutine write_checkpoint_pgc( this, restart_handle )

  use m_emf_pgc_define
  use m_emf_pgc

  use m_vdf_smooth, only : smooth
  use m_restart,    only : t_restart_handle, restart_io_write

  implicit none

  class( t_emf_pgc ), intent(in) :: this
  type( t_restart_handle ), intent(inout) :: restart_handle

  character(len=*), parameter :: &
      err_msg = 'error writing restart data for emf_pgc object.'
  integer :: ierr

  !print *, '[pgc] in write_checkpoint_pgc'

  ! write superclass checkpoint data first
  call this%t_emf%write_checkpoint( restart_handle )

  ! write checkpoint id
  restart_io_wr( p_emf_pgc_rst_id, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  ! Write restart data for laser vector potential
  call this%env%write_checkpoint( restart_handle )

  ! Write restart data for chi
  call this%chi(1)%write_checkpoint( restart_handle )

end subroutine write_checkpoint_pgc
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!       Allocate any objects contained within the EMF (t_emf) object
!-------------------------------------------------------------------------------
subroutine allocate_objs_pgc( this )

    use m_emf_pgc_define, only : t_emf_pgc
    use m_emf_pgc_diag,   only : t_diag_emf_pgc

    implicit none

    class( t_emf_pgc ), intent(inout) :: this

    ! Allocate t_diag_emf_pgc object instead of t_diag_emf
    if ( .not. associated( this%diag ) ) then
      allocate( t_diag_emf_pgc :: this%diag )
    endif

    ! Allocate superclass objects
    call this%t_emf%allocate_objs()

end subroutine
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! Sets the chi buffer to zero
!-------------------------------------------------------------------------------
subroutine reset_chi_pgc( this )

  use m_emf_pgc_define, only : t_emf_pgc

  implicit none

  class( t_emf_pgc ), intent( inout )  ::  this

  integer :: i

  !$omp parallel do
  do i = lbound(this%chi, dim=1), ubound(this%chi, dim=1)
    call this%chi(i)%zero()
  enddo
  !$omp end parallel do

end subroutine reset_chi_pgc
!-------------------------------------------------------------------------------
