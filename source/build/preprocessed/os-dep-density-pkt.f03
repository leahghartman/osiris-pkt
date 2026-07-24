# 1 "pkt/os-dep-density-pkt.f03"
# 1 "<built-in>" 1
# 1 "<built-in>" 3
# 467 "<built-in>" 3
# 1 "<command line>" 1
# 1 "<built-in>" 2
# 1 "pkt/os-dep-density-pkt.f03" 2
!===============================================================================
! File: os-dep-density-pkt.f03
!
! Module: m_dep_density_pkt
!
! Purpose:
! Implements particle-to-grid deposition and grid-to-particle interpolation
! of charged particle number density for 1D, 2D, and 3D Cartesian grids.
!
! Author: leahghartman
! Created: 12.04.2025
!===============================================================================



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
# 17 "pkt/os-dep-density-pkt.f03" 2
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
# 18 "pkt/os-dep-density-pkt.f03" 2

module m_dep_density_pkt

   use m_parameters
   use m_vdf_define, only: t_vdf

   implicit none

   private

   !> 1D structure to hold particle data
   type :: t_vp_pkt_1D
      !> Particle exact position
      real(p_k_part) :: x0
      !> Particle charge (weight)
      real(p_k_part) :: q
      !> Particle grid index
      integer :: i
      ! Particle momentum
      real(p_k_part), dimension(3) :: p
   end type t_vp_pkt_1D

   !> 2D structure to hold particle data
   type :: t_vp_pkt_2D
      !> Particle exact position
      real(p_k_part) :: x0, y0
      !> Particle charge (weight)
      real(p_k_part) :: q
      !> Particle grid index
      integer :: i, j
      !> Particle momentum
      real(p_k_part), dimension(3) :: p
   end type t_vp_pkt_2D

   !> 3D structure to hold particle data
   type :: t_vp_pkt_3D
      !> Particle exact position
      real(p_k_part) :: x0, y0, z0
      !> Particle charge (weight)
      real(p_k_part) :: q
      !> Particle grid index
      integer :: i, j, k
   end type t_vp_pkt_3D

   !-------------------------------------------
   ! Public interfaces: density deposition
   !-------------------------------------------
   interface deposit_density_1d
      module procedure deposit_density_1d
   end interface

   interface deposit_density_2d
      module procedure deposit_density_2d
   end interface

   interface deposit_density_3d
      module procedure deposit_density_3d
   end interface

   !-------------------------------------------
   ! Public interfaces: density interpolation
   !-------------------------------------------
   interface gather_density_1d
      module procedure gather_density_1d
   end interface

   interface gather_density_2d
      module procedure gather_density_2d
   end interface

   interface gather_density_3d
      module procedure gather_density_3d
   end interface

   !----------------------------------------------------
   ! Public interfaces: density gradient interpolation
   !----------------------------------------------------
   interface gather_density_grad_1d
      module procedure gather_density_grad_1d
   end interface

   interface gather_density_grad_2d
      module procedure gather_density_grad_2d
   end interface

   interface gather_density_grad_3d
      module procedure gather_density_grad_3d
   end interface

   public :: deposit_density_1d, deposit_density_2d, deposit_density_3d
   public :: gather_density_1d, gather_density_2d, gather_density_3d
   public :: gather_density_grad_1d, gather_density_grad_2d, gather_density_grad_3d
   public :: dep_den_grad_grid_1d, dep_den_grad_grid_2d

contains
   !=============================================================================
   ! B-spline shape functions and derivatives
   !
   ! This section defines 1D B-spline shape functions S(x) and their derivatives
   ! dS/dx for spline orders 1 through 4 (linear to quartic).
   !
   ! Conventions:
   ! - x is the particle position relative to the cell center
   ! and is assumed to lie in [-0.5, 0.5]
   ! - Output arrays are indexed by grid offset relative to the
   ! particle cell index
   ! - Returned derivatives are with respect to x (not divided by dx)
   !
   ! Usage:
   ! These routines are used internally by deposition and gather kernels.
   ! They are not intended to be called directly by high-level code.
   !
   !=============================================================================

   !------------------------------
   ! Linear (1st-order) B-spline
   !------------------------------
   subroutine spline_s1(x, s)

      implicit none

      real(p_k_fld), intent(in) :: x
      real(p_k_fld), dimension(0:1), intent(inout) :: s

      s(0) = 0.5_p_k_fld - x
      s(1) = 0.5_p_k_fld + x

   end subroutine spline_s1

   subroutine dspline_s1(x, ds)

      implicit none

      real(p_k_fld), intent(in) :: x
      real(p_k_fld), dimension(0:1), intent(inout) :: ds

      ds(0) = -1.0_p_k_fld
      ds(1) = 1.0_p_k_fld

   end subroutine dspline_s1

   !---------------------------------
   ! Quadratic (2nd-order) B-spline
   !---------------------------------
   subroutine spline_s2(x, s)

      implicit none

      real(p_k_fld), intent(in) :: x
      real(p_k_fld), dimension(-1:1), intent(inout) :: s

      real(p_k_fld) :: t0, t1

      t0 = 0.5_p_k_fld - x
      t1 = 0.5_p_k_fld + x

      s(-1) = 0.5_p_k_fld*t0**2
      s(0) = 0.5_p_k_fld + t0*t1
      s(1) = 0.5_p_k_fld*t1**2

   end subroutine spline_s2

   subroutine dspline_s2(x, ds)

      implicit none

      real(p_k_fld), intent(in) :: x
      real(p_k_fld), dimension(-1:1), intent(inout) :: ds

      real(p_k_fld) :: t0, t1

      t0 = 0.5_p_k_fld - x
      t1 = 0.5_p_k_fld + x

      ds(-1) = -t0
      ds(0) = -2.0_p_k_fld*x
      ds(1) = t1

   end subroutine dspline_s2

   !-----------------------------
   ! Cubic (3rd-order) B-spline
   !-----------------------------
   subroutine spline_s3(x, s)

      implicit none

      real(p_k_fld), intent(in) :: x
      real(p_k_fld), dimension(-1:2), intent(inout) :: s

      real(p_k_fld) :: t0, t1, t2, t3
      real(p_k_fld), parameter :: c1_6 = 1.0_p_k_fld/6.0_p_k_fld
      real(p_k_fld), parameter :: c2_3 = 2.0_p_k_fld/3.0_p_k_fld

      t0 = 0.5_p_k_fld - x
      t1 = 0.5_p_k_fld + x

      t2 = t0*t0
      t3 = t1*t1

      t0 = t0*t2
      t1 = t1*t3

      s(-1) = c1_6*t0
      s(0) = (c2_3 - t3) + 0.5_p_k_fld*t1
      s(1) = (c2_3 - t2) + 0.5_p_k_fld*t0
      s(2) = c1_6*t1

   end subroutine spline_s3

   subroutine dspline_s3(x, ds)

      implicit none

      real(p_k_fld), intent(in) :: x
      real(p_k_fld), dimension(-1:2), intent(inout) :: ds

      real(p_k_fld) :: t0, t1

      t0 = 0.5_p_k_fld - x
      t1 = 0.5_p_k_fld + x

      ds(-1) = -0.5_p_k_fld*t0*t0
      ds(0) = -2.0_p_k_fld*t1 + 1.5_p_k_fld*t1*t1
      ds(1) = 2.0_p_k_fld*t0 - 1.5_p_k_fld*t0*t0
      ds(2) = 0.5_p_k_fld*t1*t1

   end subroutine dspline_s3

   !-------------------------------
   ! Quartic (4th-order) B-spline
   !-------------------------------
   subroutine spline_s4(x, s)

      implicit none

      real(p_k_fld), intent(in) :: x
      real(p_k_fld), dimension(-2:2), intent(inout) :: s

      real(p_k_fld) :: t0, t1, t2

      real(p_k_fld), parameter :: c1_6 = 1.0_p_k_fld/6.0_p_k_fld
      real(p_k_fld), parameter :: c1_24 = 1.0_p_k_fld/24.0_p_k_fld
      real(p_k_fld), parameter :: c11_24 = 11.0_p_k_fld/24.0_p_k_fld

      t0 = 0.5_p_k_fld - x
      t1 = 0.5_p_k_fld + x
      t2 = t0*t1

      s(-2) = c1_24*((t0**2)*(t0**2))
      s(-1) = c1_6*(0.25_p_k_fld + t0*(1.0_p_k_fld + t0*(1.5_p_k_fld + t2)))
      s(0) = c11_24 + t2*(0.5_p_k_fld + 0.25_p_k_fld*t2)
      s(1) = c1_6*(0.25_p_k_fld + t1*(1.0_p_k_fld + t1*(1.5_p_k_fld + t2)))
      s(2) = c1_24*((t1**2)*(t1**2))

   end subroutine spline_s4

   subroutine dspline_s4(x, ds)

      implicit none

      real(p_k_fld), intent(in) :: x
      real(p_k_fld), dimension(-2:2), intent(inout) :: ds

      real(p_k_fld) :: t0, t1, t2
      real(p_k_fld), parameter :: c1_6 = 1.0_p_k_fld/6.0_p_k_fld

      t0 = 0.5_p_k_fld - x
      t1 = 0.5_p_k_fld + x
      t2 = t0*t1

      ds(-2) = -c1_6*t0**3
      ds(-1) = c1_6*(-1.0_p_k_fld - 3.0_p_k_fld*t0 + t0**2*(t0 - 3.0_p_k_fld*t1))
      ds(0) = -x*(1.0_p_k_fld + t2)
      ds(1) = c1_6*(1.0_p_k_fld + 3.0_p_k_fld*t1 + t1**2*(3.0_p_k_fld*t0 - t1))
      ds(2) = c1_6*t1**3

   end subroutine dspline_s4

!-----------------------------------------------------------------------
!> Deposits particle charge onto a 1D density grid.
!!
!! Maps particle charges to the grid using spline interpolation of the
!! specified order. The resulting density field is defined in grid space
!! and is consistent with the corresponding gather routines.
!!
!! This routine performs deposition only; no normalization by dx is
!! applied here. Any physical-space interpretation (e.g. gradients)
!! must be handled in subsequent gather operations.
!!
!! @param[inout] density Density field on grid
!! @param[inout] ix Particle cell indices
!! @param[inout] x Particle fractional positions in cell
!! @param[inout] q Particle charges
!! @param[in] np Number of particles
!! @param[in] interpolation Interpolation order (linear, quadratic, ...)
!!
!! @note Density is accumulated in grid units.
!! @warning Assumes uniform Cartesian grid.
!!
!-----------------------------------------------------------------------
   subroutine deposit_density_1d(density, ix, x, q, p, np, interpolation, photon_density_grid)

      implicit none

      type(t_vdf), intent(inout) :: density
      integer, dimension(:, :), intent(inout) :: ix
      real(p_k_part), dimension(:, :), intent(inout) :: x
      real(p_k_part), dimension(:), intent(inout) :: q
      real(p_k_part), dimension(:, :), intent(inout) :: p
      integer, intent(inout) :: np
      integer, intent(in) :: interpolation
      type(t_vdf), intent(in), optional :: photon_density_grid

      select case (interpolation)
      case (p_linear)
         call dep_density_1d_s1(density, ix, x, q, p, np, photon_density_grid)

      case (p_quadratic)
         call dep_density_1d_s2(density, ix, x, q, p, np, photon_density_grid)

      case (p_cubic)
         call dep_density_1d_s3(density, ix, x, q, p, np, photon_density_grid)

      case (p_quartic)
         call dep_density_1d_s4(density, ix, x, q, p, np, photon_density_grid)

      case default
         write(err_buf__,*) 'Not implemented';call err__("pkt/os-dep-density-pkt.f03",346)
         call abort_program(p_err_notimplemented)
      end select

   end subroutine deposit_density_1d
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!> Handles the 2D deposition of particle number density onto a grid using a chosen interpolation order
!!
!! @param density Particle number density VDF
!! @param ix Particle grid indices
!! @param x Particle simulation positions
!! @param q Particle charges (weights)
!! @param np Number of particles to deposit
!! @param interpolation Type of interpolation used for the species
   subroutine deposit_density_2d(density, ix, x, q, p, np, interpolation)

      implicit none

      type(t_vdf), intent(inout) :: density
      integer, dimension(:, :), intent(inout) :: ix
      real(p_k_part), dimension(:, :), intent(inout) :: x
      real(p_k_part), dimension(:), intent(inout) :: q
      real(p_k_part), dimension(:, :), intent(inout) :: p
      integer, intent(inout) :: np
      integer, intent(in) :: interpolation

      select case (interpolation)
      case (p_linear)
         call dep_density_2d_s1(density, ix, x, q, p, np)

      case (p_quadratic)
         call dep_density_2d_s2(density, ix, x, q, p, np)

      case (p_cubic)
         call dep_density_2d_s3(density, ix, x, q, p, np)

      case (p_quartic)
         call dep_density_2d_s4(density, ix, x, q, p, np)

      case default
         write(err_buf__,*) 'Not implemented';call err__("pkt/os-dep-density-pkt.f03",388)
         call abort_program(p_err_notimplemented)
      end select

   end subroutine deposit_density_2d
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!> Handles the 3D deposition of particle number density onto a grid using a chosen interpolation order
!!
!! @param density Particle number density VDF
!! @param ix Particle grid indices
!! @param x Particle simulation positions
!! @param q Particle charges (weights)
!! @param np Number of particles to deposit
!! @param interpolation Type of interpolation used for the species
   subroutine deposit_density_3d(density, ix, x, q, np, interpolation)

      implicit none

      type(t_vdf), intent(inout) :: density
      integer, dimension(:, :), intent(inout) :: ix
      real(p_k_part), dimension(:, :), intent(inout) :: x
      real(p_k_part), dimension(:), intent(inout) :: q
      integer, intent(inout) :: np
      integer, intent(in) :: interpolation

      select case (interpolation)
      case (p_linear)
         call dep_density_3d_s1(density, ix, x, q, np)

      case (p_quadratic)
         call dep_density_3d_s2(density, ix, x, q, np)

      case (p_cubic)
         call dep_density_3d_s3(density, ix, x, q, np)

      case (p_quartic)
         call dep_density_3d_s4(density, ix, x, q, np)

      case default
         write(err_buf__,*) 'Not implemented';call err__("pkt/os-dep-density-pkt.f03",429)
         call abort_program(p_err_notimplemented)
      end select

   end subroutine deposit_density_3d
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!> Handles the 1D interpolation of particle density back to particle positions using a chosen interpolation order
!!
!! @param den_part Container for the number density at particle positions
!! @param den_grid Particle number density VDF on the grid
!! @param ix Particle grid indices
!! @param x Particle simulation positions
!! @param q Particle charges (weights)
!! @param np Number of particles to deposit
!! @param interpolation Type of interpolation used for the species
   subroutine gather_density_1d(den_part, den_grid, ix, x, np, interpolation)

      implicit none

      real(p_k_fld), dimension(:), intent(inout) :: den_part
      type(t_vdf), intent(inout) :: den_grid
      integer, dimension(:, :), intent(inout) :: ix
      real(p_k_part), dimension(:, :), intent(inout) :: x
      integer, intent(in) :: np
      integer, intent(in) :: interpolation

      select case (interpolation)
      case (p_linear)
         call gat_density_1d_s1(den_part, den_grid, ix, x, np)

      case (p_quadratic)
         call gat_density_1d_s2(den_part, den_grid, ix, x, np)

      case (p_cubic)
         call gat_density_1d_s3(den_part, den_grid, ix, x, np)

      case (p_quartic)
         call gat_density_1d_s4(den_part, den_grid, ix, x, np)

      case default
         write(err_buf__,*) 'Not implemented';call err__("pkt/os-dep-density-pkt.f03",471)
         call abort_program(p_err_notimplemented)
      end select

   end subroutine gather_density_1d
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!> Handles the 2D interpolation of particle density back to particle positions using a chosen interpolation order
!!
!! @param den_part Container for the number density at particle positions
!! @param den_grid Particle number density VDF on the grid
!! @param ix Particle grid indices
!! @param x Particle simulation positions
!! @param q Particle charges (weights)
!! @param np Number of particles to deposit
!! @param interpolation Type of interpolation used for the species
   subroutine gather_density_2d(den_part, den_grid, ix, x, np, interpolation, field_comp)

      implicit none

      real(p_k_fld), dimension(:), intent(inout) :: den_part
      type(t_vdf), intent(inout) :: den_grid
      integer, dimension(:, :), intent(inout) :: ix
      real(p_k_part), dimension(:, :), intent(inout) :: x
      integer, intent(in) :: np
      integer, intent(in) :: interpolation
      integer, intent(in), optional :: field_comp

      integer :: fc
      fc = 1
      if (present(field_comp)) fc = field_comp

      select case (interpolation)
      case (p_linear)
         call gat_density_2d_s1(den_part, den_grid, ix, x, np, fc)

      case (p_quadratic)
         call gat_density_2d_s2(den_part, den_grid, ix, x, np, fc)

      case (p_cubic)
         call gat_density_2d_s3(den_part, den_grid, ix, x, np, fc)

      case (p_quartic)
         call gat_density_2d_s4(den_part, den_grid, ix, x, np, fc)

      case default
         write(err_buf__,*) 'Not implemented';call err__("pkt/os-dep-density-pkt.f03",518)
         call abort_program(p_err_notimplemented)
      end select

   end subroutine gather_density_2d

!-------------------------------------------------------------------------------
!> Handles the 2D interpolation of particle density back to particle positions using a chosen interpolation order
!!
!! @param den_part Container for the number density at particle positions
!! @param den_grid Particle number density VDF on the grid
!! @param ix Particle grid indices
!! @param x Particle simulation positions
!! @param q Particle charges (weights)
!! @param np Number of particles to deposit
!! @param interpolation Type of interpolation used for the species
   subroutine gather_density_3d(den_part, den_grid, ix, x, np, interpolation)

      implicit none

      real(p_k_fld), dimension(:), intent(inout) :: den_part
      type(t_vdf), intent(inout) :: den_grid
      integer, dimension(:, :), intent(inout) :: ix
      real(p_k_part), dimension(:, :), intent(inout) :: x
      integer, intent(in) :: np
      integer, intent(in) :: interpolation

      select case (interpolation)
      case (p_linear)
         call gat_density_3d_s1(den_part, den_grid, ix, x, np)

      case (p_quadratic)
         call gat_density_3d_s2(den_part, den_grid, ix, x, np)

      case (p_cubic)
         call gat_density_3d_s3(den_part, den_grid, ix, x, np)

      case (p_quartic)
         call gat_density_3d_s4(den_part, den_grid, ix, x, np)

      case default
         write(err_buf__,*) 'Not implemented';call err__("pkt/os-dep-density-pkt.f03",559)
         call abort_program(p_err_notimplemented)
      end select

   end subroutine gather_density_3d

!-------------------------------------------------------------------------------
!> Handles the 1D interpolation of particle density back to particle positions using a chosen interpolation order.
!> This subroutine will return the gradient of the particle density at the particle positions.
!!
!! @param grad_part Container for the gradient of the number density at the particle positions
!! @param den_grid Particle number density VDF on the grid
!! @param ix Particle grid indices
!! @param x Particle simulation positions
!! @param q Particle charges (weights)
!! @param np Number of particles to deposit
!! @param interpolation Type of interpolation used for the species
   subroutine gather_density_grad_1d(grad_part, den_grid, ix, x, np, dx, interpolation)

      implicit none

      real(p_k_fld), dimension(:, :), intent(inout) :: grad_part
      type(t_vdf), intent(inout) :: den_grid
      integer, dimension(:, :), intent(inout) :: ix
      real(p_k_part), dimension(:, :), intent(inout) :: x
      integer, intent(in) :: np
      real(p_double), intent(in) :: dx
      integer, intent(in) :: interpolation

      select case (interpolation)
      case (p_linear)
         call gat_density_grad_1d_s1(grad_part, den_grid, ix, x, dx, np)

      case (p_quadratic)
         call gat_density_grad_1d_s2(grad_part, den_grid, ix, x, dx, np)

      case (p_cubic)
         call gat_density_grad_1d_s3(grad_part, den_grid, ix, x, dx, np)

      case (p_quartic)
         call gat_density_grad_1d_s4(grad_part, den_grid, ix, x, dx, np)

      case default
         write(err_buf__,*) 'Not implemented';call err__("pkt/os-dep-density-pkt.f03",602)
         call abort_program(p_err_notimplemented)
      end select

   end subroutine gather_density_grad_1d
!-------------------------------------------------------------------------------

   subroutine gather_density_grad_2d(grad_part, den_grid, ix, x, np, interpolation)

      implicit none

      real(p_k_fld), dimension(:, :), intent(inout) :: grad_part
      type(t_vdf), intent(inout) :: den_grid
      integer, dimension(:, :), intent(inout) :: ix
      real(p_k_part), dimension(:, :), intent(inout) :: x
      integer, intent(in) :: np
      integer, intent(in) :: interpolation

      select case (interpolation)
      case (p_linear)
         call gat_density_grad_2d_s1(grad_part, den_grid, ix, x, np)

      case (p_quadratic)
         call gat_density_grad_2d_s2(grad_part, den_grid, ix, x, np)

      case (p_cubic)
         call gat_density_grad_2d_s3(grad_part, den_grid, ix, x, np)

      case (p_quartic)
         call gat_density_grad_2d_s4(grad_part, den_grid, ix, x, np)

      case default
         write(err_buf__,*) 'Not implemented';call err__("pkt/os-dep-density-pkt.f03",634)
         call abort_program(p_err_notimplemented)
      end select

   end subroutine gather_density_grad_2d

   subroutine gather_density_grad_3d(grad_part, den_grid, ix, x, np, interpolation)

      implicit none

      real(p_k_fld), dimension(:, :), intent(inout) :: grad_part
      type(t_vdf), intent(inout) :: den_grid
      integer, dimension(:, :), intent(inout) :: ix
      real(p_k_part), dimension(:, :), intent(inout) :: x
      integer, intent(in) :: np
      integer, intent(in) :: interpolation

      select case (interpolation)
      case (p_linear)
         call gat_density_grad_3d_s1(grad_part, den_grid, ix, x, np)

      case (p_quadratic)
         call gat_density_grad_3d_s2(grad_part, den_grid, ix, x, np)

      case (p_cubic)
         call gat_density_grad_3d_s3(grad_part, den_grid, ix, x, np)

      case (p_quartic)
         call gat_density_grad_3d_s4(grad_part, den_grid, ix, x, np)

      case default
         write(err_buf__,*) 'Not implemented';call err__("pkt/os-dep-density-pkt.f03",665)
         call abort_program(p_err_notimplemented)
      end select

   end subroutine gather_density_grad_3d

!-------------------------------------------------------------------------------
!> Gathers all the particle data together to localize memory. Removed all of the particle splitting etc. because no particles are moving across a border.
!!
!! @param ix Particle grid indices
!! @param x Particle simulation positions
!! @param q Particle charges (weights)
!! @param v Particle speeds
!! @param p Particle momenta
!! @param np Number of particles to deposit
!! @param vpbuf1D Particle buffer to localize information
!! @param nptotal Total number of particles in vpbuf1D (can differ based on vectorization with ghost particles added)
   subroutine split_1d(ix, x, q, p, np, vpbuf1D, nptotal)

      implicit none

      integer, dimension(:, :), intent(in) :: ix
      real(p_k_part), dimension(:, :), intent(in) :: x
      real(p_k_part), dimension(:), intent(in) :: q
      real(p_k_part), dimension(:, :), intent(in) :: p

      integer, intent(in) :: np ! number of particles to process

      type(t_vp_pkt_1D), dimension(:), intent(out) :: vpbuf1D
      integer, intent(out) :: nptotal

      integer :: i

      nptotal = np

      do i = 1, nptotal
         vpbuf1D(i)%x0 = x(1, i)
         vpbuf1D(i)%q = abs(q(i))
         vpbuf1D(i)%i = ix(1, i)

         ! Copy momentum
         vpbuf1D(i)%p(1) = p(1, i)
         vpbuf1D(i)%p(2) = p(2, i)
         vpbuf1D(i)%p(3) = p(3, i)
      end do

   end subroutine split_1d
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!> Gathers all the particle data together to localize memory. Removed all of the particle splitting etc. because no particles are moving across a border.
!!
!! @param ix Particle grid indices
!! @param x Particle simulation positions
!! @param q Particle charges (weights)
!! @param v Particle speeds
!! @param p Particle momenta
!! @param np Number of particles to deposit
!! @param vpbuf1D Particle buffer to localize information
!! @param nptotal Total number of particles in vpbuf1D (can differ based on vectorization with ghost particles added)
   subroutine split_2d(ix, x, q, np, vpbuf2D, nptotal, p)

      implicit none

      integer, dimension(:, :), intent(in) :: ix
      real(p_k_part), dimension(:, :), intent(in) :: x
      real(p_k_part), dimension(:), intent(in) :: q
      real(p_k_part), dimension(:, :), intent(in), optional :: p

      integer, intent(in) :: np ! number of particles to process

      type(t_vp_pkt_2D), dimension(:), intent(out) :: vpbuf2D
      integer, intent(out) :: nptotal

      integer :: i

      nptotal = np

      do i = 1, nptotal
         vpbuf2D(i)%x0 = x(1, i)
         vpbuf2D(i)%y0 = x(2, i)
         vpbuf2D(i)%q = abs(q(i))
         vpbuf2D(i)%i = ix(1, i)
         vpbuf2D(i)%j = ix(2, i)
         if (present(p)) then
            vpbuf2D(i)%p(1) = p(1, i)
            vpbuf2D(i)%p(2) = p(2, i)
            vpbuf2D(i)%p(3) = p(3, i)
         end if
      end do

   end subroutine split_2d
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!> Gathers all the particle data together to localize memory. Removed all of the particle splitting etc. because no particles are moving across a border.
!!
!! @param ix Particle grid indices
!! @param x Particle simulation positions
!! @param q Particle charges (weights)
!! @param v Particle speeds
!! @param p Particle momenta
!! @param np Number of particles to deposit
!! @param vpbuf1D Particle buffer to localize information
!! @param nptotal Total number of particles in vpbuf1D (can differ based on vectorization with ghost particles added)
   subroutine split_3d(ix, x, q, np, vpbuf3D, nptotal)

      implicit none

      integer, dimension(:, :), intent(in) :: ix
      real(p_k_part), dimension(:, :), intent(in) :: x
      real(p_k_part), dimension(:), intent(in) :: q

      integer, intent(in) :: np ! number of particles to process

      type(t_vp_pkt_3D), dimension(:), intent(out) :: vpbuf3D
      integer, intent(out) :: nptotal

      integer :: i

      nptotal = np

      do i = 1, nptotal
         vpbuf3D(i)%x0 = x(1, i)
         vpbuf3D(i)%y0 = x(2, i)
         vpbuf3D(i)%z0 = x(3, i)

         vpbuf3D(i)%q = abs(q(i))

         vpbuf3D(i)%i = ix(1, i)
         vpbuf3D(i)%j = ix(2, i)
         vpbuf3D(i)%k = ix(2, i)

      end do

   end subroutine split_3d
!-------------------------------------------------------------------------------

   ! subroutine dep_den_grad_grid_1d(grad_grid, ix, x, q, p, np, dx, photon_density_grid)
   !
   ! implicit none
   !
   ! type(t_vdf), intent(inout) :: grad_grid
   ! integer, dimension(:, :), intent(inout) :: ix
   ! real(p_k_part), dimension(:, :), intent(inout) :: x
   ! real(p_k_part), dimension(:), intent(inout) :: q
   ! real(p_k_part), dimension(:, :), intent(inout) :: p
   ! integer, intent(inout) :: np
   ! real(p_double), intent(in) :: dx
   ! type(t_vdf), intent(in), optional :: photon_density_grid
   !
   ! ! Local variables
   ! type(t_vp_pkt_1D), allocatable :: vpbuf1D(:)
   ! integer :: ip, i, k1
   ! real(p_k_part) :: x0, qi, qi_eff, u2, gamma
   ! real(p_k_fld), dimension(-1:1) :: w1p
   !
   ! ! Allocate particle buffer
   ! allocate (vpbuf1D(np))
   !
   ! ! Gather particle info
   ! call split_1d(ix, x, q, p, np, vpbuf1D, np)
   !
   ! ! Initialize gradient grid
   ! grad_grid%f1 = 0.0_p_k_fld
   !
   ! ! Loop over particles
   ! do ip = 1, np
   !
   ! i = vpbuf1D(ip)%i ! Particle grid index
   ! x0 = vpbuf1D(ip)%x0 ! Particle position
   ! qi = vpbuf1D(ip)%q ! Particle charge
   !
   ! ! Relativistic correction
   ! u2 = vpbuf1D(ip)%p(1)**2 + vpbuf1D(ip)%p(2)**2 + vpbuf1D(ip)%p(3)**2
   ! gamma = sqrt(1.0_p_k_part + u2)
   ! qi_eff = qi/gamma
   !
   ! ! Derivative weights at particle position
   ! call dspline_s2(x0, w1p)
   !
   ! ! Deposit derivative onto grid
   ! do k1 = -1, 1
   ! grad_grid%f1(1, i + k1) = grad_grid%f1(1, i + k1) - qi_eff*w1p(k1)/dx
   ! end do
   !
   ! end do
   !
   ! end subroutine dep_den_grad_grid_1d

   subroutine dep_den_grad_grid_1d(grad_grid, ix, x, q, p, np, dx, photon_density_grid)

      implicit none

      type(t_vdf), intent(inout) :: grad_grid
      integer, dimension(:, :), intent(inout) :: ix
      real(p_k_part), dimension(:, :), intent(inout) :: x
      real(p_k_part), dimension(:), intent(inout) :: q
      real(p_k_part), dimension(:, :), intent(inout) :: p
      integer, intent(inout) :: np
      real(p_double), intent(in) :: dx
      type(t_vdf), intent(in), optional :: photon_density_grid

      ! Local variables
      type(t_vp_pkt_1D), allocatable :: vpbuf1D(:)
      integer :: ip, i, k1
      real(p_k_part) :: x0, qi, qi_eff, u2, gamma, A2_local
      real(p_k_fld), dimension(-1:1) :: w1p
      real(p_k_fld), dimension(-1:1) :: w1

      ! Allocate particle buffer
      allocate (vpbuf1D(np))

      ! Gather particle info
      call split_1d(ix, x, q, p, np, vpbuf1D, np)

      ! Initialize gradient grid
      grad_grid%f1 = 0.0_p_k_fld

      ! Loop over particles
      do ip = 1, np

         i = vpbuf1D(ip)%i ! Particle grid index
         x0 = vpbuf1D(ip)%x0 ! Particle position
         qi = vpbuf1D(ip)%q ! Particle charge

         ! Relativistic correction
         u2 = vpbuf1D(ip)%p(1)**2 + vpbuf1D(ip)%p(2)**2 + vpbuf1D(ip)%p(3)**2

         ! Derivative weights at particle position (for depositing the gradient)
         call dspline_s2(x0, w1p)

         A2_local = 0.0_p_k_part
         if (present(photon_density_grid)) then
            ! Need the VALUE weights (not derivative) to interpolate
            ! the local photon density at this particle's position
            call spline_s2(x0, w1)
            do k1 = -1, 1
               A2_local = A2_local + w1(k1)*photon_density_grid%f1(1, i + k1)
            end do
         end if

         gamma = sqrt(1.0_p_k_part + u2 + 1.0_p_k_part*A2_local)
         qi_eff = qi/gamma

         ! Deposit derivative onto grid
         do k1 = -1, 1
            grad_grid%f1(1, i + k1) = grad_grid%f1(1, i + k1) - qi_eff*w1p(k1)/dx
         end do

      end do

   end subroutine dep_den_grad_grid_1d

   subroutine dep_den_grad_grid_2d(grad_grid, ix, x, q, p, np, dx1, dx2)

      implicit none

      type(t_vdf), intent(inout) :: grad_grid
      integer, dimension(:, :), intent(inout) :: ix
      real(p_k_part), dimension(:, :), intent(inout) :: x
      real(p_k_part), dimension(:), intent(inout) :: q
      real(p_k_part), dimension(:, :), intent(inout) :: p
      integer, intent(inout) :: np
      real(p_double), intent(in) :: dx1, dx2

      type(t_vp_pkt_2D), allocatable :: vpbuf2D(:)
      integer :: ip, i, j, k1, k2
      real(p_k_part) :: x0, y0, qi, qi_eff, u2, gamma
      real(p_k_fld), dimension(-1:1) :: w1, w1p, w2, w2p

      allocate (vpbuf2D(np))
      call split_2d(ix, x, q, np, vpbuf2D, np, p)

      ! grad_grid%f2 = 0.0_p_k_fld

      do ip = 1, np
         i = vpbuf2D(ip)%i
         j = vpbuf2D(ip)%j
         x0 = vpbuf2D(ip)%x0
         y0 = vpbuf2D(ip)%y0
         qi = vpbuf2D(ip)%q

         u2 = vpbuf2D(ip)%p(1)**2 + vpbuf2D(ip)%p(2)**2 + vpbuf2D(ip)%p(3)**2
         gamma = sqrt(1.0_p_k_part + u2)
         qi_eff = qi/gamma

         call spline_s2(x0, w1)
         call dspline_s2(x0, w1p)
         call spline_s2(y0, w2)
         call dspline_s2(y0, w2p)

         ! NOTE: this is hard-coded for quadratic splines at the moment
         do k1 = -1, 1
            do k2 = -1, 1
               ! Axial gradient component: W'(x1) * W(r)
               grad_grid%f2(1, i + k1, j + k2) = grad_grid%f2(1, i + k1, j + k2) &
                                                 - qi_eff*w1p(k1)*w2(k2)/dx1

               ! Radial gradient component: W(x1) * W'(r)
               grad_grid%f2(2, i + k1, j + k2) = grad_grid%f2(2, i + k1, j + k2) &
                                                 - qi_eff*w1(k1)*w2p(k2)/dx2
            end do
         end do
      end do

      deallocate (vpbuf2D)

   end subroutine dep_den_grad_grid_2d

!-------------------------------------------------------------------------------
! Generate specific template functions for linear, quadratic, cubic and quartic
! interpolation levels.
!-------------------------------------------------------------------------------






!**************************** Linear interpolation ****************************



! Lower point

! Upper point


# 1 "pkt/os-dep-density-pkt.f03" 1
!===============================================================================
! File: os-dep-density-pkt.f03
!
! Module: m_dep_density_pkt
!
! Purpose:
! Implements particle-to-grid deposition and grid-to-particle interpolation
! of charged particle number density for 1D, 2D, and 3D Cartesian grids.
!
! Author: leahghartman
! Created: 12.04.2025
!===============================================================================
# 1037 "pkt/os-dep-density-pkt.f03"
!*******************************************************************************
!
! Template function definitions for rate deposition
!
!*******************************************************************************
# 1058 "pkt/os-dep-density-pkt.f03"
!-------------------------------------------------------------------------------
!> Handles 1D particle-to-grid rate and momentum interpolation and proper momentum removal from current species' particles
!!
!! @param rate Collisional ionization rate VDF
!! @param ix Particle grid indices
!! @param x Particle simulation positions
!! @param q Particle charges (weights)
!! @param np Number of particles to deposit
subroutine dep_density_1d_s1(density, ix, x, q, p, np, photon_density_grid)

   implicit none

   type(t_vdf), intent(inout) :: density
   integer, dimension(:, :), intent(inout) :: ix
   real(p_k_part), dimension(:, :), intent(inout) :: x
   real(p_k_part), dimension(:), intent(inout) :: q
   real(p_k_part), dimension(:, :), intent(inout) :: p
   integer, intent(inout) :: np

   type(t_vdf), intent(in), optional :: photon_density_grid

   ! TODO: Add p_cache_size back in here for the dimension. Was throwing errors, so just changed for now
   type(t_vp_pkt_1D), allocatable :: vpbuf1D(:)
   ! type( t_vp_pkt_1D ), dimension(p_cache_size) :: vpbuf1D

   integer :: ip
   integer :: i, k1
   real(p_k_part) :: x0, qi, qi_eff
   real(p_k_part) :: u2, gamma, A2_local
   real(p_k_fld), dimension(0:1) :: w1

   ! Allocate our particle information container with size equal to the number of particles
   allocate (vpbuf1D(np))

   ! Call split_1d to gather all of this information into our vpbuf1D object
   call split_1d(ix, x, q, p, np, vpbuf1D, np)

   do ip = 1, np

      x0 = vpbuf1D(ip)%x0
      i = vpbuf1D(ip)%i
      qi = vpbuf1D(ip)%q

      u2 = vpbuf1D(ip)%p(1)**2 + vpbuf1D(ip)%p(2)**2 + vpbuf1D(ip)%p(3)**2

      call spline_s1(x0, w1)

      A2_local = 0.0_p_k_part
      if (present(photon_density_grid)) then
         do k1 = 0, 1
            A2_local = A2_local + w1(k1)*photon_density_grid%f1(1, i + k1)
         end do
      end if

      gamma = sqrt(1.0_p_k_part + u2 + 1.0_p_k_part*A2_local)
      qi_eff = qi/gamma

      do k1 = 0, 1
         density%f1(1, i + k1) = density%f1(1, i + k1) + qi_eff*w1(k1)
      end do

   end do

end subroutine dep_density_1d_s1
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!> Handles 1D particle-to-grid rate and momentum interpolation and proper momentum removal from current species' particles
!!
!! @param rate Collisional ionization rate VDF
!! @param ix Particle grid indices
!! @param x Particle simulation positions
!! @param q Particle charges (weights)
!! @param np Number of particles to deposit
subroutine dep_density_2d_s1(density, ix, x, q, p, np)

   implicit none

   type(t_vdf), intent(inout) :: density
   integer, dimension(:, :), intent(inout) :: ix
   real(p_k_part), dimension(:, :), intent(inout) :: x
   real(p_k_part), dimension(:), intent(inout) :: q
   real(p_k_part), dimension(:, :), intent(in) :: p
   integer, intent(inout) :: np

   ! TODO: Add p_cache_size back in here for the dimension. Was throwing errors, so just changed for now
   type(t_vp_pkt_2D), allocatable :: vpbuf2D(:)

   integer :: ip
   integer :: i, j, k1, k2
   real(p_k_part) :: x0, y0, qi, qi_eff
   real(p_k_fld), dimension(0:1) :: w1, w2
   real(p_k_part) :: u2, gamma
   real(p_k_part) :: gamma_sum, q_sum, gamma_avg

   ! Allocate our particle information container with size equal to the number of particles
   allocate (vpbuf2D(np))

   ! Call split_2d to gather all of this information into our vpbuf1D object
   call split_2d(ix, x, q, np, vpbuf2D, np, p)

   gamma_sum = 0.0
   q_sum = 0.0
   do ip = 1, np

      x0 = vpbuf2D(ip)%x0 ! 2D particle positions
      y0 = vpbuf2D(ip)%y0
      i = vpbuf2D(ip)%i ! 2D particle grid indices
      j = vpbuf2D(ip)%j
      qi = vpbuf2D(ip)%q ! Particle charge (weight)

      ! --- Calculate the relativistic correction ---
      u2 = vpbuf2D(ip)%p(1)**2 + &
           vpbuf2D(ip)%p(2)**2 + &
           vpbuf2D(ip)%p(3)**2

      gamma = sqrt(1.0_p_k_part + u2) ! NOTE: Keep the calculation of gamma the same
      qi_eff = qi/gamma
      ! ---------------------------------------------

      ! Compute the spline weights given position in two dimensions
      call spline_s1(x0, w1)
      call spline_s1(y0, w2)

      ! Deposit particle density onto grid
      do k1 = 0, 1
         do k2 = 0, 1
            density%f2(1, i + k1, j + k2) = density%f2(1, i + k1, j + k2) + qi_eff*w1(k1)*w2(k2)
         end do
      end do

   end do

end subroutine dep_density_2d_s1
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!> Handles 3D particle-to-grid rate and momentum interpolation and proper momentum removal from current species' particles
!!
!! @param rate Collisional ionization rate VDF
!! @param ix Particle grid indices
!! @param x Particle simulation positions
!! @param q Particle charges (weights)
!! @param np Number of particles to deposit
subroutine dep_density_3d_s1(density, ix, x, q, np)

   implicit none

   type(t_vdf), intent(inout) :: density
   integer, dimension(:, :), intent(inout) :: ix
   real(p_k_part), dimension(:, :), intent(inout) :: x
   real(p_k_part), dimension(:), intent(inout) :: q
   integer, intent(inout) :: np

   ! TODO: Add p_cache_size back in here for the dimension. Was throwing errors, so just changed for now
   type(t_vp_pkt_3D), allocatable :: vpbuf3D(:)

   integer :: ip
   integer :: i, j, k, k1, k2, k3
   real(p_k_part) :: x0, y0, z0, qi
   real(p_k_fld), dimension(0:1) :: w1, w2, w3

   ! Allocate our particle information container with size equal to the number of particles
   allocate (vpbuf3D(np))

   ! Call split_3d to gather all of this information into our vpbuf1D object
   call split_3d(ix, x, q, np, vpbuf3D, np)

   do ip = 1, np

      x0 = vpbuf3D(ip)%x0 ! 3D particle positions
      y0 = vpbuf3D(ip)%y0
      z0 = vpbuf3D(ip)%z0

      i = vpbuf3D(ip)%i ! 3D particle grid indices
      j = vpbuf3D(ip)%j
      k = vpbuf3D(ip)%k

      qi = vpbuf3D(ip)%q ! Particle charge (weight)

      ! Compute the spline weights given position in two dimensions
      call spline_s1(x0, w1)
      call spline_s1(y0, w2)
      call spline_s1(z0, w3)

      ! Deposit particle density onto grid
      do k1 = 0, 1
         do k2 = 0, 1
            do k3 = 0, 1
               density%f3(1, i + k1, j + k2, k + k3) = &
                  density%f3(1, i + k1, j + k2, k + k3) &
                  + qi*w1(k1)*w2(k2)*w3(k3)
            end do
         end do
      end do

   end do

end subroutine dep_density_3d_s1
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!> Handles 1D particle-to-grid rate and momentum interpolation and proper momentum removal from current species' particles
!!
!! @param rate Collisional ionization rate VDF
!! @param ix Particle grid indices
!! @param x Particle simulation positions
!! @param q Particle charges (weights)
!! @param np Number of particles to deposit
subroutine gat_density_1d_s1(den_part, den_grid, ix, x, np)

   implicit none

   real(p_k_fld), dimension(:), intent(inout) :: den_part
   type(t_vdf), intent(inout) :: den_grid
   integer, dimension(:, :), intent(inout) :: ix
   real(p_k_part), dimension(:, :), intent(inout) :: x
   integer, intent(in) :: np

   type(t_vp_pkt_1D), dimension(p_cache_size) :: vpbuf1D

   integer :: i, ip, k1
   real(p_k_part) :: x0
   real(p_k_fld), dimension(0:1) :: w1

   ! Loop over particles
   do ip = 1, np

      x0 = x(1, ip)
      i = ix(1, ip)

      ! Compute spline weights
      call spline_s1(x0, w1)

      ! Initialize the output
      den_part(ip) = 0.0_p_k_fld

      ! Interpolate from grid to particle
      do k1 = 0, 1
         den_part(ip) = den_part(ip) + w1(k1)*den_grid%f1(1, i + k1)
      end do

   end do

end subroutine gat_density_1d_s1
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!> Handles 2D particle-to-grid rate and momentum interpolation and proper momentum removal from current species' particles
!!
!! @param rate Collisional ionization rate VDF
!! @param ix Particle grid indices
!! @param x Particle simulation positions
!! @param q Particle charges (weights)
!! @param np Number of particles to deposit
subroutine gat_density_2d_s1(den_part, den_grid, ix, x, np, fc)

   implicit none

   real(p_k_fld), dimension(:), intent(inout) :: den_part
   type(t_vdf), intent(inout) :: den_grid
   integer, dimension(:, :), intent(inout) :: ix
   real(p_k_part), dimension(:, :), intent(inout) :: x
   integer, intent(in) :: np
   integer, intent(in) :: fc

   integer :: ip
   integer :: i, j, k1, k2
   real(p_k_part) :: x0, y0
   real(p_k_fld), dimension(0:1) :: w1, w2

   ! Loop over particles
   do ip = 1, np

      x0 = x(1, ip)
      y0 = x(2, ip)

      i = ix(1, ip)
      j = ix(2, ip)

      ! Compute spline weights
      call spline_s1(x0, w1)
      call spline_s1(y0, w2)

      ! Initialize the output
      den_part(ip) = 0.0_p_k_fld

      ! Interpolate from grid to particle
      do k1 = 0, 1
         do k2 = 0, 1
            den_part(ip) = den_part(ip) &
                           + den_grid%f2(fc, i + k1, j + k2)*w1(k1)*w2(k2)
         end do
      end do

   end do

end subroutine gat_density_2d_s1
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!> Handles 3D particle-to-grid rate and momentum interpolation and proper momentum removal from current species' particles
!!
!! @param rate Collisional ionization rate VDF
!! @param ix Particle grid indices
!! @param x Particle simulation positions
!! @param q Particle charges (weights)
!! @param np Number of particles to deposit
subroutine gat_density_3d_s1(den_part, den_grid, ix, x, np)

   implicit none

   real(p_k_fld), dimension(:), intent(inout) :: den_part
   type(t_vdf), intent(inout) :: den_grid
   integer, dimension(:, :), intent(inout) :: ix
   real(p_k_part), dimension(:, :), intent(inout) :: x
   integer, intent(in) :: np

   integer :: ip
   integer :: i, j, k, k1, k2, k3
   real(p_k_part) :: x0, y0, z0
   real(p_k_fld), dimension(0:1) :: w1, w2, w3

   ! Loop over particles
   do ip = 1, np

      x0 = x(1, ip)
      y0 = x(2, ip)
      z0 = x(3, ip)

      i = ix(1, ip)
      j = ix(2, ip)
      k = ix(3, ip)

      ! Compute spline weights
      call spline_s1(x0, w1)
      call spline_s1(y0, w2)
      call spline_s1(z0, w3)

      ! Initialize the output
      den_part(ip) = 0.0_p_k_fld

      ! Interpolate from grid to particle
      do k1 = 0, 1
         do k2 = 0, 1
            do k3 = 0, 1
               den_part(ip) = den_part(ip) &
                              + den_grid%f3(1, i + k1, j + k2, k + k3)*w1(k1)*w2(k2)*w3(k3)
            end do
         end do
      end do

   end do

end subroutine gat_density_3d_s1
!-------------------------------------------------------------------------------

subroutine gat_density_grad_1d_s1(grad_part, den_grid, ix, x, dx, np)

   implicit none

   real(p_k_fld), dimension(:, :), intent(out) :: grad_part
   type(t_vdf), intent(inout) :: den_grid
   integer, dimension(:, :), intent(in) :: ix
   real(p_k_part), dimension(:, :), intent(in) :: x
   real(p_double), intent(in) :: dx
   integer, intent(in) :: np

   integer :: i
   integer :: ip, k1
   real(p_k_part) :: x0
   real(p_k_fld), dimension(0:1) :: w1p

   do ip = 1, np

      i = ix(1, ip)
      x0 = x(1, ip)

      call dspline_s1(x0, w1p)

      grad_part(1, ip) = 0.0_p_k_fld

      do k1 = 0, 1
         grad_part(1, ip) = grad_part(1, ip) + w1p(k1)*den_grid%f1(1, i + k1)/dx
      end do
   end do

end subroutine gat_density_grad_1d_s1

subroutine gat_density_grad_2d_s1(grad_part, den_grid, ix, x, np)

   implicit none

   real(p_k_fld), dimension(:, :), intent(out) :: grad_part
   type(t_vdf), intent(inout) :: den_grid
   integer, dimension(:, :), intent(in) :: ix
   real(p_k_part), dimension(:, :), intent(in) :: x
   integer, intent(in) :: np

   integer :: i, j
   integer :: ip, k1, k2
   real(p_k_part) :: x0, y0
   real(p_k_fld), dimension(0:1) :: w1, w2
   real(p_k_fld), dimension(0:1) :: w1p, w2p

   do ip = 1, np

      x0 = x(1, ip)
      y0 = x(2, ip)

      i = ix(1, ip)
      j = ix(2, ip)

      call spline_s1(x0, w1)
      call spline_s1(y0, w2)
      call dspline_s1(x0, w1p)
      call dspline_s1(y0, w2p)

      grad_part(:, ip) = 0.0_p_k_fld

      do k1 = 0, 1
         do k2 = 0, 1

            grad_part(1, ip) = grad_part(1, ip) &
                               + w1(k1)*w2p(k2)*den_grid%f2(1, i + k1, j + k2)

            grad_part(2, ip) = grad_part(2, ip) &
                               + w1p(k1)*w2(k2)*den_grid%f2(1, i + k1, j + k2)

         end do
      end do

   end do

end subroutine gat_density_grad_2d_s1

subroutine gat_density_grad_3d_s1(grad_part, den_grid, ix, x, np)

   implicit none

   real(p_k_fld), dimension(:, :), intent(out) :: grad_part
   type(t_vdf), intent(inout) :: den_grid
   integer, dimension(:, :), intent(in) :: ix
   real(p_k_part), dimension(:, :), intent(in) :: x
   integer, intent(in) :: np

   integer :: i, j, k
   integer :: ip, k1, k2, k3
   real(p_k_part) :: x0, y0, z0
   real(p_k_fld), dimension(0:1) :: w1, w2, w3
   real(p_k_fld), dimension(0:1) :: w1p, w2p, w3p

   do ip = 1, np

      x0 = x(1, ip)
      y0 = x(2, ip)
      z0 = x(3, ip)

      i = ix(1, ip)
      j = ix(2, ip)
      k = ix(3, ip)

      call spline_s1(x0, w1)
      call spline_s1(y0, w2)
      call spline_s1(z0, w3)
      call dspline_s1(x0, w1p)
      call dspline_s1(y0, w2p)
      call dspline_s1(z0, w3p)

      grad_part(:, ip) = 0.0_p_k_fld

      do k1 = 0, 1
         do k2 = 0, 1
            do k3 = 0, 1

               grad_part(1, ip) = grad_part(1, ip) &
                                  + w1p(k1)*w2(k2)*w3(k3) &
                                  *den_grid%f3(1, i + k1, j + k2, k + k3)

               grad_part(2, ip) = grad_part(2, ip) &
                                  + w1(k1)*w2p(k2)*w3(k3) &
                                  *den_grid%f3(1, i + k1, j + k2, k + k3)

               grad_part(3, ip) = grad_part(3, ip) &
                                  + w1(k1)*w2(k2)*w3p(k3) &
                                  *den_grid%f3(1, i + k1, j + k2, k + k3)

            end do
         end do
      end do

   end do

end subroutine gat_density_grad_3d_s1

!-------------------------------------------------------------------------------

! Clear template definitions
# 995 "pkt/os-dep-density-pkt.f03" 2

!*************************** Quadratic interpolation ***************************



! Lower point

! Upper point


# 1 "pkt/os-dep-density-pkt.f03" 1
!===============================================================================
! File: os-dep-density-pkt.f03
!
! Module: m_dep_density_pkt
!
! Purpose:
! Implements particle-to-grid deposition and grid-to-particle interpolation
! of charged particle number density for 1D, 2D, and 3D Cartesian grids.
!
! Author: leahghartman
! Created: 12.04.2025
!===============================================================================
# 1037 "pkt/os-dep-density-pkt.f03"
!*******************************************************************************
!
! Template function definitions for rate deposition
!
!*******************************************************************************
# 1058 "pkt/os-dep-density-pkt.f03"
!-------------------------------------------------------------------------------
!> Handles 1D particle-to-grid rate and momentum interpolation and proper momentum removal from current species' particles
!!
!! @param rate Collisional ionization rate VDF
!! @param ix Particle grid indices
!! @param x Particle simulation positions
!! @param q Particle charges (weights)
!! @param np Number of particles to deposit
subroutine dep_density_1d_s2(density, ix, x, q, p, np, photon_density_grid)

   implicit none

   type(t_vdf), intent(inout) :: density
   integer, dimension(:, :), intent(inout) :: ix
   real(p_k_part), dimension(:, :), intent(inout) :: x
   real(p_k_part), dimension(:), intent(inout) :: q
   real(p_k_part), dimension(:, :), intent(inout) :: p
   integer, intent(inout) :: np

   type(t_vdf), intent(in), optional :: photon_density_grid

   ! TODO: Add p_cache_size back in here for the dimension. Was throwing errors, so just changed for now
   type(t_vp_pkt_1D), allocatable :: vpbuf1D(:)
   ! type( t_vp_pkt_1D ), dimension(p_cache_size) :: vpbuf1D

   integer :: ip
   integer :: i, k1
   real(p_k_part) :: x0, qi, qi_eff
   real(p_k_part) :: u2, gamma, A2_local
   real(p_k_fld), dimension(-1:1) :: w1

   ! Allocate our particle information container with size equal to the number of particles
   allocate (vpbuf1D(np))

   ! Call split_1d to gather all of this information into our vpbuf1D object
   call split_1d(ix, x, q, p, np, vpbuf1D, np)

   do ip = 1, np

      x0 = vpbuf1D(ip)%x0
      i = vpbuf1D(ip)%i
      qi = vpbuf1D(ip)%q

      u2 = vpbuf1D(ip)%p(1)**2 + vpbuf1D(ip)%p(2)**2 + vpbuf1D(ip)%p(3)**2

      call spline_s2(x0, w1)

      A2_local = 0.0_p_k_part
      if (present(photon_density_grid)) then
         do k1 = -1, 1
            A2_local = A2_local + w1(k1)*photon_density_grid%f1(1, i + k1)
         end do
      end if

      gamma = sqrt(1.0_p_k_part + u2 + 1.0_p_k_part*A2_local)
      qi_eff = qi/gamma

      do k1 = -1, 1
         density%f1(1, i + k1) = density%f1(1, i + k1) + qi_eff*w1(k1)
      end do

   end do

end subroutine dep_density_1d_s2
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!> Handles 1D particle-to-grid rate and momentum interpolation and proper momentum removal from current species' particles
!!
!! @param rate Collisional ionization rate VDF
!! @param ix Particle grid indices
!! @param x Particle simulation positions
!! @param q Particle charges (weights)
!! @param np Number of particles to deposit
subroutine dep_density_2d_s2(density, ix, x, q, p, np)

   implicit none

   type(t_vdf), intent(inout) :: density
   integer, dimension(:, :), intent(inout) :: ix
   real(p_k_part), dimension(:, :), intent(inout) :: x
   real(p_k_part), dimension(:), intent(inout) :: q
   real(p_k_part), dimension(:, :), intent(in) :: p
   integer, intent(inout) :: np

   ! TODO: Add p_cache_size back in here for the dimension. Was throwing errors, so just changed for now
   type(t_vp_pkt_2D), allocatable :: vpbuf2D(:)

   integer :: ip
   integer :: i, j, k1, k2
   real(p_k_part) :: x0, y0, qi, qi_eff
   real(p_k_fld), dimension(-1:1) :: w1, w2
   real(p_k_part) :: u2, gamma
   real(p_k_part) :: gamma_sum, q_sum, gamma_avg

   ! Allocate our particle information container with size equal to the number of particles
   allocate (vpbuf2D(np))

   ! Call split_2d to gather all of this information into our vpbuf1D object
   call split_2d(ix, x, q, np, vpbuf2D, np, p)

   gamma_sum = 0.0
   q_sum = 0.0
   do ip = 1, np

      x0 = vpbuf2D(ip)%x0 ! 2D particle positions
      y0 = vpbuf2D(ip)%y0
      i = vpbuf2D(ip)%i ! 2D particle grid indices
      j = vpbuf2D(ip)%j
      qi = vpbuf2D(ip)%q ! Particle charge (weight)

      ! --- Calculate the relativistic correction ---
      u2 = vpbuf2D(ip)%p(1)**2 + &
           vpbuf2D(ip)%p(2)**2 + &
           vpbuf2D(ip)%p(3)**2

      gamma = sqrt(1.0_p_k_part + u2) ! NOTE: Keep the calculation of gamma the same
      qi_eff = qi/gamma
      ! ---------------------------------------------

      ! Compute the spline weights given position in two dimensions
      call spline_s2(x0, w1)
      call spline_s2(y0, w2)

      ! Deposit particle density onto grid
      do k1 = -1, 1
         do k2 = -1, 1
            density%f2(1, i + k1, j + k2) = density%f2(1, i + k1, j + k2) + qi_eff*w1(k1)*w2(k2)
         end do
      end do

   end do

end subroutine dep_density_2d_s2
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!> Handles 3D particle-to-grid rate and momentum interpolation and proper momentum removal from current species' particles
!!
!! @param rate Collisional ionization rate VDF
!! @param ix Particle grid indices
!! @param x Particle simulation positions
!! @param q Particle charges (weights)
!! @param np Number of particles to deposit
subroutine dep_density_3d_s2(density, ix, x, q, np)

   implicit none

   type(t_vdf), intent(inout) :: density
   integer, dimension(:, :), intent(inout) :: ix
   real(p_k_part), dimension(:, :), intent(inout) :: x
   real(p_k_part), dimension(:), intent(inout) :: q
   integer, intent(inout) :: np

   ! TODO: Add p_cache_size back in here for the dimension. Was throwing errors, so just changed for now
   type(t_vp_pkt_3D), allocatable :: vpbuf3D(:)

   integer :: ip
   integer :: i, j, k, k1, k2, k3
   real(p_k_part) :: x0, y0, z0, qi
   real(p_k_fld), dimension(-1:1) :: w1, w2, w3

   ! Allocate our particle information container with size equal to the number of particles
   allocate (vpbuf3D(np))

   ! Call split_3d to gather all of this information into our vpbuf1D object
   call split_3d(ix, x, q, np, vpbuf3D, np)

   do ip = 1, np

      x0 = vpbuf3D(ip)%x0 ! 3D particle positions
      y0 = vpbuf3D(ip)%y0
      z0 = vpbuf3D(ip)%z0

      i = vpbuf3D(ip)%i ! 3D particle grid indices
      j = vpbuf3D(ip)%j
      k = vpbuf3D(ip)%k

      qi = vpbuf3D(ip)%q ! Particle charge (weight)

      ! Compute the spline weights given position in two dimensions
      call spline_s2(x0, w1)
      call spline_s2(y0, w2)
      call spline_s2(z0, w3)

      ! Deposit particle density onto grid
      do k1 = -1, 1
         do k2 = -1, 1
            do k3 = -1, 1
               density%f3(1, i + k1, j + k2, k + k3) = &
                  density%f3(1, i + k1, j + k2, k + k3) &
                  + qi*w1(k1)*w2(k2)*w3(k3)
            end do
         end do
      end do

   end do

end subroutine dep_density_3d_s2
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!> Handles 1D particle-to-grid rate and momentum interpolation and proper momentum removal from current species' particles
!!
!! @param rate Collisional ionization rate VDF
!! @param ix Particle grid indices
!! @param x Particle simulation positions
!! @param q Particle charges (weights)
!! @param np Number of particles to deposit
subroutine gat_density_1d_s2(den_part, den_grid, ix, x, np)

   implicit none

   real(p_k_fld), dimension(:), intent(inout) :: den_part
   type(t_vdf), intent(inout) :: den_grid
   integer, dimension(:, :), intent(inout) :: ix
   real(p_k_part), dimension(:, :), intent(inout) :: x
   integer, intent(in) :: np

   type(t_vp_pkt_1D), dimension(p_cache_size) :: vpbuf1D

   integer :: i, ip, k1
   real(p_k_part) :: x0
   real(p_k_fld), dimension(-1:1) :: w1

   ! Loop over particles
   do ip = 1, np

      x0 = x(1, ip)
      i = ix(1, ip)

      ! Compute spline weights
      call spline_s2(x0, w1)

      ! Initialize the output
      den_part(ip) = 0.0_p_k_fld

      ! Interpolate from grid to particle
      do k1 = -1, 1
         den_part(ip) = den_part(ip) + w1(k1)*den_grid%f1(1, i + k1)
      end do

   end do

end subroutine gat_density_1d_s2
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!> Handles 2D particle-to-grid rate and momentum interpolation and proper momentum removal from current species' particles
!!
!! @param rate Collisional ionization rate VDF
!! @param ix Particle grid indices
!! @param x Particle simulation positions
!! @param q Particle charges (weights)
!! @param np Number of particles to deposit
subroutine gat_density_2d_s2(den_part, den_grid, ix, x, np, fc)

   implicit none

   real(p_k_fld), dimension(:), intent(inout) :: den_part
   type(t_vdf), intent(inout) :: den_grid
   integer, dimension(:, :), intent(inout) :: ix
   real(p_k_part), dimension(:, :), intent(inout) :: x
   integer, intent(in) :: np
   integer, intent(in) :: fc

   integer :: ip
   integer :: i, j, k1, k2
   real(p_k_part) :: x0, y0
   real(p_k_fld), dimension(-1:1) :: w1, w2

   ! Loop over particles
   do ip = 1, np

      x0 = x(1, ip)
      y0 = x(2, ip)

      i = ix(1, ip)
      j = ix(2, ip)

      ! Compute spline weights
      call spline_s2(x0, w1)
      call spline_s2(y0, w2)

      ! Initialize the output
      den_part(ip) = 0.0_p_k_fld

      ! Interpolate from grid to particle
      do k1 = -1, 1
         do k2 = -1, 1
            den_part(ip) = den_part(ip) &
                           + den_grid%f2(fc, i + k1, j + k2)*w1(k1)*w2(k2)
         end do
      end do

   end do

end subroutine gat_density_2d_s2
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!> Handles 3D particle-to-grid rate and momentum interpolation and proper momentum removal from current species' particles
!!
!! @param rate Collisional ionization rate VDF
!! @param ix Particle grid indices
!! @param x Particle simulation positions
!! @param q Particle charges (weights)
!! @param np Number of particles to deposit
subroutine gat_density_3d_s2(den_part, den_grid, ix, x, np)

   implicit none

   real(p_k_fld), dimension(:), intent(inout) :: den_part
   type(t_vdf), intent(inout) :: den_grid
   integer, dimension(:, :), intent(inout) :: ix
   real(p_k_part), dimension(:, :), intent(inout) :: x
   integer, intent(in) :: np

   integer :: ip
   integer :: i, j, k, k1, k2, k3
   real(p_k_part) :: x0, y0, z0
   real(p_k_fld), dimension(-1:1) :: w1, w2, w3

   ! Loop over particles
   do ip = 1, np

      x0 = x(1, ip)
      y0 = x(2, ip)
      z0 = x(3, ip)

      i = ix(1, ip)
      j = ix(2, ip)
      k = ix(3, ip)

      ! Compute spline weights
      call spline_s2(x0, w1)
      call spline_s2(y0, w2)
      call spline_s2(z0, w3)

      ! Initialize the output
      den_part(ip) = 0.0_p_k_fld

      ! Interpolate from grid to particle
      do k1 = -1, 1
         do k2 = -1, 1
            do k3 = -1, 1
               den_part(ip) = den_part(ip) &
                              + den_grid%f3(1, i + k1, j + k2, k + k3)*w1(k1)*w2(k2)*w3(k3)
            end do
         end do
      end do

   end do

end subroutine gat_density_3d_s2
!-------------------------------------------------------------------------------

subroutine gat_density_grad_1d_s2(grad_part, den_grid, ix, x, dx, np)

   implicit none

   real(p_k_fld), dimension(:, :), intent(out) :: grad_part
   type(t_vdf), intent(inout) :: den_grid
   integer, dimension(:, :), intent(in) :: ix
   real(p_k_part), dimension(:, :), intent(in) :: x
   real(p_double), intent(in) :: dx
   integer, intent(in) :: np

   integer :: i
   integer :: ip, k1
   real(p_k_part) :: x0
   real(p_k_fld), dimension(-1:1) :: w1p

   do ip = 1, np

      i = ix(1, ip)
      x0 = x(1, ip)

      call dspline_s2(x0, w1p)

      grad_part(1, ip) = 0.0_p_k_fld

      do k1 = -1, 1
         grad_part(1, ip) = grad_part(1, ip) + w1p(k1)*den_grid%f1(1, i + k1)/dx
      end do
   end do

end subroutine gat_density_grad_1d_s2

subroutine gat_density_grad_2d_s2(grad_part, den_grid, ix, x, np)

   implicit none

   real(p_k_fld), dimension(:, :), intent(out) :: grad_part
   type(t_vdf), intent(inout) :: den_grid
   integer, dimension(:, :), intent(in) :: ix
   real(p_k_part), dimension(:, :), intent(in) :: x
   integer, intent(in) :: np

   integer :: i, j
   integer :: ip, k1, k2
   real(p_k_part) :: x0, y0
   real(p_k_fld), dimension(-1:1) :: w1, w2
   real(p_k_fld), dimension(-1:1) :: w1p, w2p

   do ip = 1, np

      x0 = x(1, ip)
      y0 = x(2, ip)

      i = ix(1, ip)
      j = ix(2, ip)

      call spline_s2(x0, w1)
      call spline_s2(y0, w2)
      call dspline_s2(x0, w1p)
      call dspline_s2(y0, w2p)

      grad_part(:, ip) = 0.0_p_k_fld

      do k1 = -1, 1
         do k2 = -1, 1

            grad_part(1, ip) = grad_part(1, ip) &
                               + w1(k1)*w2p(k2)*den_grid%f2(1, i + k1, j + k2)

            grad_part(2, ip) = grad_part(2, ip) &
                               + w1p(k1)*w2(k2)*den_grid%f2(1, i + k1, j + k2)

         end do
      end do

   end do

end subroutine gat_density_grad_2d_s2

subroutine gat_density_grad_3d_s2(grad_part, den_grid, ix, x, np)

   implicit none

   real(p_k_fld), dimension(:, :), intent(out) :: grad_part
   type(t_vdf), intent(inout) :: den_grid
   integer, dimension(:, :), intent(in) :: ix
   real(p_k_part), dimension(:, :), intent(in) :: x
   integer, intent(in) :: np

   integer :: i, j, k
   integer :: ip, k1, k2, k3
   real(p_k_part) :: x0, y0, z0
   real(p_k_fld), dimension(-1:1) :: w1, w2, w3
   real(p_k_fld), dimension(-1:1) :: w1p, w2p, w3p

   do ip = 1, np

      x0 = x(1, ip)
      y0 = x(2, ip)
      z0 = x(3, ip)

      i = ix(1, ip)
      j = ix(2, ip)
      k = ix(3, ip)

      call spline_s2(x0, w1)
      call spline_s2(y0, w2)
      call spline_s2(z0, w3)
      call dspline_s2(x0, w1p)
      call dspline_s2(y0, w2p)
      call dspline_s2(z0, w3p)

      grad_part(:, ip) = 0.0_p_k_fld

      do k1 = -1, 1
         do k2 = -1, 1
            do k3 = -1, 1

               grad_part(1, ip) = grad_part(1, ip) &
                                  + w1p(k1)*w2(k2)*w3(k3) &
                                  *den_grid%f3(1, i + k1, j + k2, k + k3)

               grad_part(2, ip) = grad_part(2, ip) &
                                  + w1(k1)*w2p(k2)*w3(k3) &
                                  *den_grid%f3(1, i + k1, j + k2, k + k3)

               grad_part(3, ip) = grad_part(3, ip) &
                                  + w1(k1)*w2(k2)*w3p(k3) &
                                  *den_grid%f3(1, i + k1, j + k2, k + k3)

            end do
         end do
      end do

   end do

end subroutine gat_density_grad_3d_s2

!-------------------------------------------------------------------------------

! Clear template definitions
# 1006 "pkt/os-dep-density-pkt.f03" 2

!***************************** Cubic interpolation *****************************



! Lower point

! Upper point


# 1 "pkt/os-dep-density-pkt.f03" 1
!===============================================================================
! File: os-dep-density-pkt.f03
!
! Module: m_dep_density_pkt
!
! Purpose:
! Implements particle-to-grid deposition and grid-to-particle interpolation
! of charged particle number density for 1D, 2D, and 3D Cartesian grids.
!
! Author: leahghartman
! Created: 12.04.2025
!===============================================================================
# 1037 "pkt/os-dep-density-pkt.f03"
!*******************************************************************************
!
! Template function definitions for rate deposition
!
!*******************************************************************************
# 1058 "pkt/os-dep-density-pkt.f03"
!-------------------------------------------------------------------------------
!> Handles 1D particle-to-grid rate and momentum interpolation and proper momentum removal from current species' particles
!!
!! @param rate Collisional ionization rate VDF
!! @param ix Particle grid indices
!! @param x Particle simulation positions
!! @param q Particle charges (weights)
!! @param np Number of particles to deposit
subroutine dep_density_1d_s3(density, ix, x, q, p, np, photon_density_grid)

   implicit none

   type(t_vdf), intent(inout) :: density
   integer, dimension(:, :), intent(inout) :: ix
   real(p_k_part), dimension(:, :), intent(inout) :: x
   real(p_k_part), dimension(:), intent(inout) :: q
   real(p_k_part), dimension(:, :), intent(inout) :: p
   integer, intent(inout) :: np

   type(t_vdf), intent(in), optional :: photon_density_grid

   ! TODO: Add p_cache_size back in here for the dimension. Was throwing errors, so just changed for now
   type(t_vp_pkt_1D), allocatable :: vpbuf1D(:)
   ! type( t_vp_pkt_1D ), dimension(p_cache_size) :: vpbuf1D

   integer :: ip
   integer :: i, k1
   real(p_k_part) :: x0, qi, qi_eff
   real(p_k_part) :: u2, gamma, A2_local
   real(p_k_fld), dimension(-1:2) :: w1

   ! Allocate our particle information container with size equal to the number of particles
   allocate (vpbuf1D(np))

   ! Call split_1d to gather all of this information into our vpbuf1D object
   call split_1d(ix, x, q, p, np, vpbuf1D, np)

   do ip = 1, np

      x0 = vpbuf1D(ip)%x0
      i = vpbuf1D(ip)%i
      qi = vpbuf1D(ip)%q

      u2 = vpbuf1D(ip)%p(1)**2 + vpbuf1D(ip)%p(2)**2 + vpbuf1D(ip)%p(3)**2

      call spline_s3(x0, w1)

      A2_local = 0.0_p_k_part
      if (present(photon_density_grid)) then
         do k1 = -1, 2
            A2_local = A2_local + w1(k1)*photon_density_grid%f1(1, i + k1)
         end do
      end if

      gamma = sqrt(1.0_p_k_part + u2 + 1.0_p_k_part*A2_local)
      qi_eff = qi/gamma

      do k1 = -1, 2
         density%f1(1, i + k1) = density%f1(1, i + k1) + qi_eff*w1(k1)
      end do

   end do

end subroutine dep_density_1d_s3
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!> Handles 1D particle-to-grid rate and momentum interpolation and proper momentum removal from current species' particles
!!
!! @param rate Collisional ionization rate VDF
!! @param ix Particle grid indices
!! @param x Particle simulation positions
!! @param q Particle charges (weights)
!! @param np Number of particles to deposit
subroutine dep_density_2d_s3(density, ix, x, q, p, np)

   implicit none

   type(t_vdf), intent(inout) :: density
   integer, dimension(:, :), intent(inout) :: ix
   real(p_k_part), dimension(:, :), intent(inout) :: x
   real(p_k_part), dimension(:), intent(inout) :: q
   real(p_k_part), dimension(:, :), intent(in) :: p
   integer, intent(inout) :: np

   ! TODO: Add p_cache_size back in here for the dimension. Was throwing errors, so just changed for now
   type(t_vp_pkt_2D), allocatable :: vpbuf2D(:)

   integer :: ip
   integer :: i, j, k1, k2
   real(p_k_part) :: x0, y0, qi, qi_eff
   real(p_k_fld), dimension(-1:2) :: w1, w2
   real(p_k_part) :: u2, gamma
   real(p_k_part) :: gamma_sum, q_sum, gamma_avg

   ! Allocate our particle information container with size equal to the number of particles
   allocate (vpbuf2D(np))

   ! Call split_2d to gather all of this information into our vpbuf1D object
   call split_2d(ix, x, q, np, vpbuf2D, np, p)

   gamma_sum = 0.0
   q_sum = 0.0
   do ip = 1, np

      x0 = vpbuf2D(ip)%x0 ! 2D particle positions
      y0 = vpbuf2D(ip)%y0
      i = vpbuf2D(ip)%i ! 2D particle grid indices
      j = vpbuf2D(ip)%j
      qi = vpbuf2D(ip)%q ! Particle charge (weight)

      ! --- Calculate the relativistic correction ---
      u2 = vpbuf2D(ip)%p(1)**2 + &
           vpbuf2D(ip)%p(2)**2 + &
           vpbuf2D(ip)%p(3)**2

      gamma = sqrt(1.0_p_k_part + u2) ! NOTE: Keep the calculation of gamma the same
      qi_eff = qi/gamma
      ! ---------------------------------------------

      ! Compute the spline weights given position in two dimensions
      call spline_s3(x0, w1)
      call spline_s3(y0, w2)

      ! Deposit particle density onto grid
      do k1 = -1, 2
         do k2 = -1, 2
            density%f2(1, i + k1, j + k2) = density%f2(1, i + k1, j + k2) + qi_eff*w1(k1)*w2(k2)
         end do
      end do

   end do

end subroutine dep_density_2d_s3
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!> Handles 3D particle-to-grid rate and momentum interpolation and proper momentum removal from current species' particles
!!
!! @param rate Collisional ionization rate VDF
!! @param ix Particle grid indices
!! @param x Particle simulation positions
!! @param q Particle charges (weights)
!! @param np Number of particles to deposit
subroutine dep_density_3d_s3(density, ix, x, q, np)

   implicit none

   type(t_vdf), intent(inout) :: density
   integer, dimension(:, :), intent(inout) :: ix
   real(p_k_part), dimension(:, :), intent(inout) :: x
   real(p_k_part), dimension(:), intent(inout) :: q
   integer, intent(inout) :: np

   ! TODO: Add p_cache_size back in here for the dimension. Was throwing errors, so just changed for now
   type(t_vp_pkt_3D), allocatable :: vpbuf3D(:)

   integer :: ip
   integer :: i, j, k, k1, k2, k3
   real(p_k_part) :: x0, y0, z0, qi
   real(p_k_fld), dimension(-1:2) :: w1, w2, w3

   ! Allocate our particle information container with size equal to the number of particles
   allocate (vpbuf3D(np))

   ! Call split_3d to gather all of this information into our vpbuf1D object
   call split_3d(ix, x, q, np, vpbuf3D, np)

   do ip = 1, np

      x0 = vpbuf3D(ip)%x0 ! 3D particle positions
      y0 = vpbuf3D(ip)%y0
      z0 = vpbuf3D(ip)%z0

      i = vpbuf3D(ip)%i ! 3D particle grid indices
      j = vpbuf3D(ip)%j
      k = vpbuf3D(ip)%k

      qi = vpbuf3D(ip)%q ! Particle charge (weight)

      ! Compute the spline weights given position in two dimensions
      call spline_s3(x0, w1)
      call spline_s3(y0, w2)
      call spline_s3(z0, w3)

      ! Deposit particle density onto grid
      do k1 = -1, 2
         do k2 = -1, 2
            do k3 = -1, 2
               density%f3(1, i + k1, j + k2, k + k3) = &
                  density%f3(1, i + k1, j + k2, k + k3) &
                  + qi*w1(k1)*w2(k2)*w3(k3)
            end do
         end do
      end do

   end do

end subroutine dep_density_3d_s3
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!> Handles 1D particle-to-grid rate and momentum interpolation and proper momentum removal from current species' particles
!!
!! @param rate Collisional ionization rate VDF
!! @param ix Particle grid indices
!! @param x Particle simulation positions
!! @param q Particle charges (weights)
!! @param np Number of particles to deposit
subroutine gat_density_1d_s3(den_part, den_grid, ix, x, np)

   implicit none

   real(p_k_fld), dimension(:), intent(inout) :: den_part
   type(t_vdf), intent(inout) :: den_grid
   integer, dimension(:, :), intent(inout) :: ix
   real(p_k_part), dimension(:, :), intent(inout) :: x
   integer, intent(in) :: np

   type(t_vp_pkt_1D), dimension(p_cache_size) :: vpbuf1D

   integer :: i, ip, k1
   real(p_k_part) :: x0
   real(p_k_fld), dimension(-1:2) :: w1

   ! Loop over particles
   do ip = 1, np

      x0 = x(1, ip)
      i = ix(1, ip)

      ! Compute spline weights
      call spline_s3(x0, w1)

      ! Initialize the output
      den_part(ip) = 0.0_p_k_fld

      ! Interpolate from grid to particle
      do k1 = -1, 2
         den_part(ip) = den_part(ip) + w1(k1)*den_grid%f1(1, i + k1)
      end do

   end do

end subroutine gat_density_1d_s3
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!> Handles 2D particle-to-grid rate and momentum interpolation and proper momentum removal from current species' particles
!!
!! @param rate Collisional ionization rate VDF
!! @param ix Particle grid indices
!! @param x Particle simulation positions
!! @param q Particle charges (weights)
!! @param np Number of particles to deposit
subroutine gat_density_2d_s3(den_part, den_grid, ix, x, np, fc)

   implicit none

   real(p_k_fld), dimension(:), intent(inout) :: den_part
   type(t_vdf), intent(inout) :: den_grid
   integer, dimension(:, :), intent(inout) :: ix
   real(p_k_part), dimension(:, :), intent(inout) :: x
   integer, intent(in) :: np
   integer, intent(in) :: fc

   integer :: ip
   integer :: i, j, k1, k2
   real(p_k_part) :: x0, y0
   real(p_k_fld), dimension(-1:2) :: w1, w2

   ! Loop over particles
   do ip = 1, np

      x0 = x(1, ip)
      y0 = x(2, ip)

      i = ix(1, ip)
      j = ix(2, ip)

      ! Compute spline weights
      call spline_s3(x0, w1)
      call spline_s3(y0, w2)

      ! Initialize the output
      den_part(ip) = 0.0_p_k_fld

      ! Interpolate from grid to particle
      do k1 = -1, 2
         do k2 = -1, 2
            den_part(ip) = den_part(ip) &
                           + den_grid%f2(fc, i + k1, j + k2)*w1(k1)*w2(k2)
         end do
      end do

   end do

end subroutine gat_density_2d_s3
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!> Handles 3D particle-to-grid rate and momentum interpolation and proper momentum removal from current species' particles
!!
!! @param rate Collisional ionization rate VDF
!! @param ix Particle grid indices
!! @param x Particle simulation positions
!! @param q Particle charges (weights)
!! @param np Number of particles to deposit
subroutine gat_density_3d_s3(den_part, den_grid, ix, x, np)

   implicit none

   real(p_k_fld), dimension(:), intent(inout) :: den_part
   type(t_vdf), intent(inout) :: den_grid
   integer, dimension(:, :), intent(inout) :: ix
   real(p_k_part), dimension(:, :), intent(inout) :: x
   integer, intent(in) :: np

   integer :: ip
   integer :: i, j, k, k1, k2, k3
   real(p_k_part) :: x0, y0, z0
   real(p_k_fld), dimension(-1:2) :: w1, w2, w3

   ! Loop over particles
   do ip = 1, np

      x0 = x(1, ip)
      y0 = x(2, ip)
      z0 = x(3, ip)

      i = ix(1, ip)
      j = ix(2, ip)
      k = ix(3, ip)

      ! Compute spline weights
      call spline_s3(x0, w1)
      call spline_s3(y0, w2)
      call spline_s3(z0, w3)

      ! Initialize the output
      den_part(ip) = 0.0_p_k_fld

      ! Interpolate from grid to particle
      do k1 = -1, 2
         do k2 = -1, 2
            do k3 = -1, 2
               den_part(ip) = den_part(ip) &
                              + den_grid%f3(1, i + k1, j + k2, k + k3)*w1(k1)*w2(k2)*w3(k3)
            end do
         end do
      end do

   end do

end subroutine gat_density_3d_s3
!-------------------------------------------------------------------------------

subroutine gat_density_grad_1d_s3(grad_part, den_grid, ix, x, dx, np)

   implicit none

   real(p_k_fld), dimension(:, :), intent(out) :: grad_part
   type(t_vdf), intent(inout) :: den_grid
   integer, dimension(:, :), intent(in) :: ix
   real(p_k_part), dimension(:, :), intent(in) :: x
   real(p_double), intent(in) :: dx
   integer, intent(in) :: np

   integer :: i
   integer :: ip, k1
   real(p_k_part) :: x0
   real(p_k_fld), dimension(-1:2) :: w1p

   do ip = 1, np

      i = ix(1, ip)
      x0 = x(1, ip)

      call dspline_s3(x0, w1p)

      grad_part(1, ip) = 0.0_p_k_fld

      do k1 = -1, 2
         grad_part(1, ip) = grad_part(1, ip) + w1p(k1)*den_grid%f1(1, i + k1)/dx
      end do
   end do

end subroutine gat_density_grad_1d_s3

subroutine gat_density_grad_2d_s3(grad_part, den_grid, ix, x, np)

   implicit none

   real(p_k_fld), dimension(:, :), intent(out) :: grad_part
   type(t_vdf), intent(inout) :: den_grid
   integer, dimension(:, :), intent(in) :: ix
   real(p_k_part), dimension(:, :), intent(in) :: x
   integer, intent(in) :: np

   integer :: i, j
   integer :: ip, k1, k2
   real(p_k_part) :: x0, y0
   real(p_k_fld), dimension(-1:2) :: w1, w2
   real(p_k_fld), dimension(-1:2) :: w1p, w2p

   do ip = 1, np

      x0 = x(1, ip)
      y0 = x(2, ip)

      i = ix(1, ip)
      j = ix(2, ip)

      call spline_s3(x0, w1)
      call spline_s3(y0, w2)
      call dspline_s3(x0, w1p)
      call dspline_s3(y0, w2p)

      grad_part(:, ip) = 0.0_p_k_fld

      do k1 = -1, 2
         do k2 = -1, 2

            grad_part(1, ip) = grad_part(1, ip) &
                               + w1(k1)*w2p(k2)*den_grid%f2(1, i + k1, j + k2)

            grad_part(2, ip) = grad_part(2, ip) &
                               + w1p(k1)*w2(k2)*den_grid%f2(1, i + k1, j + k2)

         end do
      end do

   end do

end subroutine gat_density_grad_2d_s3

subroutine gat_density_grad_3d_s3(grad_part, den_grid, ix, x, np)

   implicit none

   real(p_k_fld), dimension(:, :), intent(out) :: grad_part
   type(t_vdf), intent(inout) :: den_grid
   integer, dimension(:, :), intent(in) :: ix
   real(p_k_part), dimension(:, :), intent(in) :: x
   integer, intent(in) :: np

   integer :: i, j, k
   integer :: ip, k1, k2, k3
   real(p_k_part) :: x0, y0, z0
   real(p_k_fld), dimension(-1:2) :: w1, w2, w3
   real(p_k_fld), dimension(-1:2) :: w1p, w2p, w3p

   do ip = 1, np

      x0 = x(1, ip)
      y0 = x(2, ip)
      z0 = x(3, ip)

      i = ix(1, ip)
      j = ix(2, ip)
      k = ix(3, ip)

      call spline_s3(x0, w1)
      call spline_s3(y0, w2)
      call spline_s3(z0, w3)
      call dspline_s3(x0, w1p)
      call dspline_s3(y0, w2p)
      call dspline_s3(z0, w3p)

      grad_part(:, ip) = 0.0_p_k_fld

      do k1 = -1, 2
         do k2 = -1, 2
            do k3 = -1, 2

               grad_part(1, ip) = grad_part(1, ip) &
                                  + w1p(k1)*w2(k2)*w3(k3) &
                                  *den_grid%f3(1, i + k1, j + k2, k + k3)

               grad_part(2, ip) = grad_part(2, ip) &
                                  + w1(k1)*w2p(k2)*w3(k3) &
                                  *den_grid%f3(1, i + k1, j + k2, k + k3)

               grad_part(3, ip) = grad_part(3, ip) &
                                  + w1(k1)*w2(k2)*w3p(k3) &
                                  *den_grid%f3(1, i + k1, j + k2, k + k3)

            end do
         end do
      end do

   end do

end subroutine gat_density_grad_3d_s3

!-------------------------------------------------------------------------------

! Clear template definitions
# 1017 "pkt/os-dep-density-pkt.f03" 2

!**************************** Quartic interpolation ****************************



! Lower point

! Upper point


# 1 "pkt/os-dep-density-pkt.f03" 1
!===============================================================================
! File: os-dep-density-pkt.f03
!
! Module: m_dep_density_pkt
!
! Purpose:
! Implements particle-to-grid deposition and grid-to-particle interpolation
! of charged particle number density for 1D, 2D, and 3D Cartesian grids.
!
! Author: leahghartman
! Created: 12.04.2025
!===============================================================================
# 1037 "pkt/os-dep-density-pkt.f03"
!*******************************************************************************
!
! Template function definitions for rate deposition
!
!*******************************************************************************
# 1058 "pkt/os-dep-density-pkt.f03"
!-------------------------------------------------------------------------------
!> Handles 1D particle-to-grid rate and momentum interpolation and proper momentum removal from current species' particles
!!
!! @param rate Collisional ionization rate VDF
!! @param ix Particle grid indices
!! @param x Particle simulation positions
!! @param q Particle charges (weights)
!! @param np Number of particles to deposit
subroutine dep_density_1d_s4(density, ix, x, q, p, np, photon_density_grid)

   implicit none

   type(t_vdf), intent(inout) :: density
   integer, dimension(:, :), intent(inout) :: ix
   real(p_k_part), dimension(:, :), intent(inout) :: x
   real(p_k_part), dimension(:), intent(inout) :: q
   real(p_k_part), dimension(:, :), intent(inout) :: p
   integer, intent(inout) :: np

   type(t_vdf), intent(in), optional :: photon_density_grid

   ! TODO: Add p_cache_size back in here for the dimension. Was throwing errors, so just changed for now
   type(t_vp_pkt_1D), allocatable :: vpbuf1D(:)
   ! type( t_vp_pkt_1D ), dimension(p_cache_size) :: vpbuf1D

   integer :: ip
   integer :: i, k1
   real(p_k_part) :: x0, qi, qi_eff
   real(p_k_part) :: u2, gamma, A2_local
   real(p_k_fld), dimension(-2:2) :: w1

   ! Allocate our particle information container with size equal to the number of particles
   allocate (vpbuf1D(np))

   ! Call split_1d to gather all of this information into our vpbuf1D object
   call split_1d(ix, x, q, p, np, vpbuf1D, np)

   do ip = 1, np

      x0 = vpbuf1D(ip)%x0
      i = vpbuf1D(ip)%i
      qi = vpbuf1D(ip)%q

      u2 = vpbuf1D(ip)%p(1)**2 + vpbuf1D(ip)%p(2)**2 + vpbuf1D(ip)%p(3)**2

      call spline_s4(x0, w1)

      A2_local = 0.0_p_k_part
      if (present(photon_density_grid)) then
         do k1 = -2, 2
            A2_local = A2_local + w1(k1)*photon_density_grid%f1(1, i + k1)
         end do
      end if

      gamma = sqrt(1.0_p_k_part + u2 + 1.0_p_k_part*A2_local)
      qi_eff = qi/gamma

      do k1 = -2, 2
         density%f1(1, i + k1) = density%f1(1, i + k1) + qi_eff*w1(k1)
      end do

   end do

end subroutine dep_density_1d_s4
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!> Handles 1D particle-to-grid rate and momentum interpolation and proper momentum removal from current species' particles
!!
!! @param rate Collisional ionization rate VDF
!! @param ix Particle grid indices
!! @param x Particle simulation positions
!! @param q Particle charges (weights)
!! @param np Number of particles to deposit
subroutine dep_density_2d_s4(density, ix, x, q, p, np)

   implicit none

   type(t_vdf), intent(inout) :: density
   integer, dimension(:, :), intent(inout) :: ix
   real(p_k_part), dimension(:, :), intent(inout) :: x
   real(p_k_part), dimension(:), intent(inout) :: q
   real(p_k_part), dimension(:, :), intent(in) :: p
   integer, intent(inout) :: np

   ! TODO: Add p_cache_size back in here for the dimension. Was throwing errors, so just changed for now
   type(t_vp_pkt_2D), allocatable :: vpbuf2D(:)

   integer :: ip
   integer :: i, j, k1, k2
   real(p_k_part) :: x0, y0, qi, qi_eff
   real(p_k_fld), dimension(-2:2) :: w1, w2
   real(p_k_part) :: u2, gamma
   real(p_k_part) :: gamma_sum, q_sum, gamma_avg

   ! Allocate our particle information container with size equal to the number of particles
   allocate (vpbuf2D(np))

   ! Call split_2d to gather all of this information into our vpbuf1D object
   call split_2d(ix, x, q, np, vpbuf2D, np, p)

   gamma_sum = 0.0
   q_sum = 0.0
   do ip = 1, np

      x0 = vpbuf2D(ip)%x0 ! 2D particle positions
      y0 = vpbuf2D(ip)%y0
      i = vpbuf2D(ip)%i ! 2D particle grid indices
      j = vpbuf2D(ip)%j
      qi = vpbuf2D(ip)%q ! Particle charge (weight)

      ! --- Calculate the relativistic correction ---
      u2 = vpbuf2D(ip)%p(1)**2 + &
           vpbuf2D(ip)%p(2)**2 + &
           vpbuf2D(ip)%p(3)**2

      gamma = sqrt(1.0_p_k_part + u2) ! NOTE: Keep the calculation of gamma the same
      qi_eff = qi/gamma
      ! ---------------------------------------------

      ! Compute the spline weights given position in two dimensions
      call spline_s4(x0, w1)
      call spline_s4(y0, w2)

      ! Deposit particle density onto grid
      do k1 = -2, 2
         do k2 = -2, 2
            density%f2(1, i + k1, j + k2) = density%f2(1, i + k1, j + k2) + qi_eff*w1(k1)*w2(k2)
         end do
      end do

   end do

end subroutine dep_density_2d_s4
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!> Handles 3D particle-to-grid rate and momentum interpolation and proper momentum removal from current species' particles
!!
!! @param rate Collisional ionization rate VDF
!! @param ix Particle grid indices
!! @param x Particle simulation positions
!! @param q Particle charges (weights)
!! @param np Number of particles to deposit
subroutine dep_density_3d_s4(density, ix, x, q, np)

   implicit none

   type(t_vdf), intent(inout) :: density
   integer, dimension(:, :), intent(inout) :: ix
   real(p_k_part), dimension(:, :), intent(inout) :: x
   real(p_k_part), dimension(:), intent(inout) :: q
   integer, intent(inout) :: np

   ! TODO: Add p_cache_size back in here for the dimension. Was throwing errors, so just changed for now
   type(t_vp_pkt_3D), allocatable :: vpbuf3D(:)

   integer :: ip
   integer :: i, j, k, k1, k2, k3
   real(p_k_part) :: x0, y0, z0, qi
   real(p_k_fld), dimension(-2:2) :: w1, w2, w3

   ! Allocate our particle information container with size equal to the number of particles
   allocate (vpbuf3D(np))

   ! Call split_3d to gather all of this information into our vpbuf1D object
   call split_3d(ix, x, q, np, vpbuf3D, np)

   do ip = 1, np

      x0 = vpbuf3D(ip)%x0 ! 3D particle positions
      y0 = vpbuf3D(ip)%y0
      z0 = vpbuf3D(ip)%z0

      i = vpbuf3D(ip)%i ! 3D particle grid indices
      j = vpbuf3D(ip)%j
      k = vpbuf3D(ip)%k

      qi = vpbuf3D(ip)%q ! Particle charge (weight)

      ! Compute the spline weights given position in two dimensions
      call spline_s4(x0, w1)
      call spline_s4(y0, w2)
      call spline_s4(z0, w3)

      ! Deposit particle density onto grid
      do k1 = -2, 2
         do k2 = -2, 2
            do k3 = -2, 2
               density%f3(1, i + k1, j + k2, k + k3) = &
                  density%f3(1, i + k1, j + k2, k + k3) &
                  + qi*w1(k1)*w2(k2)*w3(k3)
            end do
         end do
      end do

   end do

end subroutine dep_density_3d_s4
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!> Handles 1D particle-to-grid rate and momentum interpolation and proper momentum removal from current species' particles
!!
!! @param rate Collisional ionization rate VDF
!! @param ix Particle grid indices
!! @param x Particle simulation positions
!! @param q Particle charges (weights)
!! @param np Number of particles to deposit
subroutine gat_density_1d_s4(den_part, den_grid, ix, x, np)

   implicit none

   real(p_k_fld), dimension(:), intent(inout) :: den_part
   type(t_vdf), intent(inout) :: den_grid
   integer, dimension(:, :), intent(inout) :: ix
   real(p_k_part), dimension(:, :), intent(inout) :: x
   integer, intent(in) :: np

   type(t_vp_pkt_1D), dimension(p_cache_size) :: vpbuf1D

   integer :: i, ip, k1
   real(p_k_part) :: x0
   real(p_k_fld), dimension(-2:2) :: w1

   ! Loop over particles
   do ip = 1, np

      x0 = x(1, ip)
      i = ix(1, ip)

      ! Compute spline weights
      call spline_s4(x0, w1)

      ! Initialize the output
      den_part(ip) = 0.0_p_k_fld

      ! Interpolate from grid to particle
      do k1 = -2, 2
         den_part(ip) = den_part(ip) + w1(k1)*den_grid%f1(1, i + k1)
      end do

   end do

end subroutine gat_density_1d_s4
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!> Handles 2D particle-to-grid rate and momentum interpolation and proper momentum removal from current species' particles
!!
!! @param rate Collisional ionization rate VDF
!! @param ix Particle grid indices
!! @param x Particle simulation positions
!! @param q Particle charges (weights)
!! @param np Number of particles to deposit
subroutine gat_density_2d_s4(den_part, den_grid, ix, x, np, fc)

   implicit none

   real(p_k_fld), dimension(:), intent(inout) :: den_part
   type(t_vdf), intent(inout) :: den_grid
   integer, dimension(:, :), intent(inout) :: ix
   real(p_k_part), dimension(:, :), intent(inout) :: x
   integer, intent(in) :: np
   integer, intent(in) :: fc

   integer :: ip
   integer :: i, j, k1, k2
   real(p_k_part) :: x0, y0
   real(p_k_fld), dimension(-2:2) :: w1, w2

   ! Loop over particles
   do ip = 1, np

      x0 = x(1, ip)
      y0 = x(2, ip)

      i = ix(1, ip)
      j = ix(2, ip)

      ! Compute spline weights
      call spline_s4(x0, w1)
      call spline_s4(y0, w2)

      ! Initialize the output
      den_part(ip) = 0.0_p_k_fld

      ! Interpolate from grid to particle
      do k1 = -2, 2
         do k2 = -2, 2
            den_part(ip) = den_part(ip) &
                           + den_grid%f2(fc, i + k1, j + k2)*w1(k1)*w2(k2)
         end do
      end do

   end do

end subroutine gat_density_2d_s4
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!> Handles 3D particle-to-grid rate and momentum interpolation and proper momentum removal from current species' particles
!!
!! @param rate Collisional ionization rate VDF
!! @param ix Particle grid indices
!! @param x Particle simulation positions
!! @param q Particle charges (weights)
!! @param np Number of particles to deposit
subroutine gat_density_3d_s4(den_part, den_grid, ix, x, np)

   implicit none

   real(p_k_fld), dimension(:), intent(inout) :: den_part
   type(t_vdf), intent(inout) :: den_grid
   integer, dimension(:, :), intent(inout) :: ix
   real(p_k_part), dimension(:, :), intent(inout) :: x
   integer, intent(in) :: np

   integer :: ip
   integer :: i, j, k, k1, k2, k3
   real(p_k_part) :: x0, y0, z0
   real(p_k_fld), dimension(-2:2) :: w1, w2, w3

   ! Loop over particles
   do ip = 1, np

      x0 = x(1, ip)
      y0 = x(2, ip)
      z0 = x(3, ip)

      i = ix(1, ip)
      j = ix(2, ip)
      k = ix(3, ip)

      ! Compute spline weights
      call spline_s4(x0, w1)
      call spline_s4(y0, w2)
      call spline_s4(z0, w3)

      ! Initialize the output
      den_part(ip) = 0.0_p_k_fld

      ! Interpolate from grid to particle
      do k1 = -2, 2
         do k2 = -2, 2
            do k3 = -2, 2
               den_part(ip) = den_part(ip) &
                              + den_grid%f3(1, i + k1, j + k2, k + k3)*w1(k1)*w2(k2)*w3(k3)
            end do
         end do
      end do

   end do

end subroutine gat_density_3d_s4
!-------------------------------------------------------------------------------

subroutine gat_density_grad_1d_s4(grad_part, den_grid, ix, x, dx, np)

   implicit none

   real(p_k_fld), dimension(:, :), intent(out) :: grad_part
   type(t_vdf), intent(inout) :: den_grid
   integer, dimension(:, :), intent(in) :: ix
   real(p_k_part), dimension(:, :), intent(in) :: x
   real(p_double), intent(in) :: dx
   integer, intent(in) :: np

   integer :: i
   integer :: ip, k1
   real(p_k_part) :: x0
   real(p_k_fld), dimension(-2:2) :: w1p

   do ip = 1, np

      i = ix(1, ip)
      x0 = x(1, ip)

      call dspline_s4(x0, w1p)

      grad_part(1, ip) = 0.0_p_k_fld

      do k1 = -2, 2
         grad_part(1, ip) = grad_part(1, ip) + w1p(k1)*den_grid%f1(1, i + k1)/dx
      end do
   end do

end subroutine gat_density_grad_1d_s4

subroutine gat_density_grad_2d_s4(grad_part, den_grid, ix, x, np)

   implicit none

   real(p_k_fld), dimension(:, :), intent(out) :: grad_part
   type(t_vdf), intent(inout) :: den_grid
   integer, dimension(:, :), intent(in) :: ix
   real(p_k_part), dimension(:, :), intent(in) :: x
   integer, intent(in) :: np

   integer :: i, j
   integer :: ip, k1, k2
   real(p_k_part) :: x0, y0
   real(p_k_fld), dimension(-2:2) :: w1, w2
   real(p_k_fld), dimension(-2:2) :: w1p, w2p

   do ip = 1, np

      x0 = x(1, ip)
      y0 = x(2, ip)

      i = ix(1, ip)
      j = ix(2, ip)

      call spline_s4(x0, w1)
      call spline_s4(y0, w2)
      call dspline_s4(x0, w1p)
      call dspline_s4(y0, w2p)

      grad_part(:, ip) = 0.0_p_k_fld

      do k1 = -2, 2
         do k2 = -2, 2

            grad_part(1, ip) = grad_part(1, ip) &
                               + w1(k1)*w2p(k2)*den_grid%f2(1, i + k1, j + k2)

            grad_part(2, ip) = grad_part(2, ip) &
                               + w1p(k1)*w2(k2)*den_grid%f2(1, i + k1, j + k2)

         end do
      end do

   end do

end subroutine gat_density_grad_2d_s4

subroutine gat_density_grad_3d_s4(grad_part, den_grid, ix, x, np)

   implicit none

   real(p_k_fld), dimension(:, :), intent(out) :: grad_part
   type(t_vdf), intent(inout) :: den_grid
   integer, dimension(:, :), intent(in) :: ix
   real(p_k_part), dimension(:, :), intent(in) :: x
   integer, intent(in) :: np

   integer :: i, j, k
   integer :: ip, k1, k2, k3
   real(p_k_part) :: x0, y0, z0
   real(p_k_fld), dimension(-2:2) :: w1, w2, w3
   real(p_k_fld), dimension(-2:2) :: w1p, w2p, w3p

   do ip = 1, np

      x0 = x(1, ip)
      y0 = x(2, ip)
      z0 = x(3, ip)

      i = ix(1, ip)
      j = ix(2, ip)
      k = ix(3, ip)

      call spline_s4(x0, w1)
      call spline_s4(y0, w2)
      call spline_s4(z0, w3)
      call dspline_s4(x0, w1p)
      call dspline_s4(y0, w2p)
      call dspline_s4(z0, w3p)

      grad_part(:, ip) = 0.0_p_k_fld

      do k1 = -2, 2
         do k2 = -2, 2
            do k3 = -2, 2

               grad_part(1, ip) = grad_part(1, ip) &
                                  + w1p(k1)*w2(k2)*w3(k3) &
                                  *den_grid%f3(1, i + k1, j + k2, k + k3)

               grad_part(2, ip) = grad_part(2, ip) &
                                  + w1(k1)*w2p(k2)*w3(k3) &
                                  *den_grid%f3(1, i + k1, j + k2, k + k3)

               grad_part(3, ip) = grad_part(3, ip) &
                                  + w1(k1)*w2(k2)*w3p(k3) &
                                  *den_grid%f3(1, i + k1, j + k2, k + k3)

            end do
         end do
      end do

   end do

end subroutine gat_density_grad_3d_s4

!-------------------------------------------------------------------------------

! Clear template definitions
# 1028 "pkt/os-dep-density-pkt.f03" 2

!*******************************************************************************

end module m_dep_density_pkt

!#endif
