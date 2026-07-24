#include "os-config.h"
#include "os-preprocess.fpp"

module m_photons_pkt

#include "memory/memory.h"

   use m_system
   use m_parameters
   use m_math

   use m_species_define, only: p_cell_near
   use m_photons_define_pkt, only: t_photons_pkt
   use m_species_define, only: p_ene_recalc
   use m_emf_define, only: t_emf
   use m_pkt
   use m_vdf_define
   use m_species_charge, only: deposit_rho
   use m_dep_density_pkt

   implicit none

   private

   public advance_photons_1d, advance_photons_2d, advance_photons_3d, advance_photons_cyl_2d
   public ntrim_pkt

contains

   subroutine advance_photons_1d(this, den_ph, grad_ph, den_grid, grad_grid, gdt, i0, i1)

      implicit none
      integer, parameter :: rank = 1

      class(t_photons_pkt), intent(inout) :: this
      real(p_k_fld), intent(inout)           :: den_ph(:)
      real(p_k_fld), intent(inout)           :: grad_ph(:, :)
      type(t_vdf), intent(inout)             :: den_grid, grad_grid
      real(p_double), intent(in)           :: gdt
      integer, intent(in)                  :: i0, i1

      real(p_k_part) :: rdx, dt
      integer :: ptrcur, np, pp, i
      real(p_k_part), dimension(rank, p_cache_size) :: xbuf
      integer, dimension(rank, p_cache_size) :: dxi
      real(p_k_part), dimension(p_cache_size)      :: invsq, dt_invsq
      real(p_k_part) :: k2, omega2, omega02
      real(p_k_part) :: arg
      real(p_k_part), parameter :: eps = 1.0e-14_p_k_part
      real(p_k_fld), dimension(p_cache_size)           :: den_ph_pred
      real(p_k_fld), dimension(p_x_dim, p_cache_size)  :: grad_ph_pred
      real(p_k_part), dimension(rank, p_cache_size)     :: x_pred
      integer, dimension(rank, p_cache_size)     :: ix_pred
      real(p_k_part), dimension(p_cache_size)           :: p1_pred
      real(p_k_part), dimension(p_cache_size)           :: invsq_pred

      rdx = real(1.0_p_double/this%dx(1), p_k_part)
      dt = real(gdt, p_k_part)

      do ptrcur = i0, i1, p_cache_size
         if (ptrcur + p_cache_size > i1) then
            np = i1 - ptrcur + 1
         else
            np = p_cache_size
         end if

         ! Compute 1/omega for the current positions and momenta
         pp = ptrcur
         do i = 1, np
            k2 = this%p(1, pp)**2 + this%p(2, pp)**2 + this%p(3, pp)**2
            omega2 = den_ph(pp)
            invsq(i) = 1.0_p_k_part/sqrt(k2 + omega2)
            dt_invsq(i) = dt*invsq(i)
            pp = pp + 1
         end do

         ! FIRST: use forward Euler as a predictor for the position and momentum
         pp = ptrcur
         do i = 1, np
            p1_pred(i) = this%p(1, pp) - 0.5_p_k_part*dt_invsq(i)*grad_ph(1, pp)
            x_pred(1, i) = this%x(1, pp) + this%p(1, pp)*dt_invsq(i)*rdx
            ix_pred(1, i) = this%ix(1, pp) + ntrim_pkt(x_pred(1, i))
            x_pred(1, i) = x_pred(1, i) - ntrim_pkt(x_pred(1, i))
            pp = pp + 1
         end do

         ! Next, use the predicted positions and momenta to...
         ! TODO: Right now this call is only using quadratic interpolation
         call gather_density_1d(den_ph_pred, den_grid, ix_pred, x_pred, np, 2)
         call gather_density_1d(grad_ph_pred(1, :), grad_grid, ix_pred, x_pred, np, 2)

         ! Finally, use the trapezoidal method to update both the momenta and positions
         pp = ptrcur
         do i = 1, np
            k2 = p1_pred(i)**2 + this%p(2, pp)**2 + this%p(3, pp)**2
            omega2 = den_ph_pred(i)
            invsq_pred(i) = 1.0_p_k_part/sqrt(k2 + omega2)

            xbuf(1, i) = this%x(1, pp) + 0.5_p_k_part*(this%p(1, pp)*invsq(i) + &
                                                       p1_pred(i)*invsq_pred(i))*dt*rdx
            dxi(1, i) = ntrim_pkt(xbuf(1, i))

            this%p(1, pp) = this%p(1, pp) - 0.5_p_k_part*dt &
                            *0.5_p_k_part*(invsq(i)*grad_ph(1, pp) + &
                                                    invsq_pred(i)*grad_ph_pred(1, i))

            pp = pp + 1
         end do

         ! Do the final position update with the weird stuff OSIRIS does
         pp = ptrcur
         do i = 1, np
            this%x(1, pp) = xbuf(1, i) - dxi(1, i)
            this%ix(1, pp) = this%ix(1, pp) + dxi(1, i)
            pp = pp + 1
         end do

      end do

   end subroutine advance_photons_1d

!---------------------------------------------------------------------------------------------------
   subroutine advance_photons_2d(this, gdt, i0, i1)
!---------------------------------------------------------------------------------------------------

      implicit none

      integer, parameter :: rank = 2

      ! dummy variables
      class(t_photons_pkt), intent(inout) :: this
      real(p_double), intent(in) :: gdt
      integer, intent(in) :: i0, i1

      real(p_k_part), dimension(rank) :: rdx
      integer :: i, pp, np, ptrcur

      real(p_k_part), dimension(rank, p_cache_size) :: xbuf
      integer, dimension(rank, p_cache_size)        :: dxi
      real(p_k_part), dimension(p_cache_size)      :: rp, rpp
      real(p_k_part) :: dt

      ! executable statements
      rdx(1) = real(1.0_p_double/this%dx(1), p_k_part)
      rdx(2) = real(1.0_p_double/this%dx(2), p_k_part)
      dt = real(gdt, p_k_part)

      ! advance position of photons i0 to i1 in chunks of p_cache_size
      do ptrcur = i0, i1, p_cache_size

         ! check if last copy of table and set np
         if (ptrcur + p_cache_size > i1) then
            np = i1 - ptrcur + 1
         else
            np = p_cache_size
         end if

         ! this type of loop is actually faster than
         ! using a forall construct
         pp = ptrcur

         do i = 1, np
            rp(i) = 1.0_p_k_part/ &
                    sqrt(this%p(1, pp)**2 + this%p(2, pp)**2 + this%p(3, pp)**2)
            rpp(i) = dt*rp(i)
            pp = pp + 1
         end do

         ! reference implementation
         pp = ptrcur
         do i = 1, np
            xbuf(1, i) = this%x(1, pp) + this%p(1, pp)*rpp(i)*rdx(1)
            xbuf(2, i) = this%x(2, pp) + this%p(2, pp)*rpp(i)*rdx(2)
            dxi(1, i) = ntrim_pkt(xbuf(1, i))
            dxi(2, i) = ntrim_pkt(xbuf(2, i))

            pp = pp + 1
         end do

         ! copy data from buffer to species data trimming positions
         pp = ptrcur
         do i = 1, np
            this%x(1, pp) = xbuf(1, i) - dxi(1, i)
            this%x(2, pp) = xbuf(2, i) - dxi(2, i)
            this%ix(1, pp) = this%ix(1, pp) + dxi(1, i)
            this%ix(2, pp) = this%ix(2, pp) + dxi(2, i)

            pp = pp + 1
         end do

      end do

   end subroutine advance_photons_2d
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
   subroutine advance_photons_3d(this, gdt, i0, i1)
!---------------------------------------------------------------------------------------------------

      implicit none

      integer, parameter :: rank = 3

      ! dummy variables

      class(t_photons_pkt), intent(inout) :: this
      real(p_double), intent(in) :: gdt
      integer, intent(in) :: i0, i1

      real(p_k_part), dimension(rank) :: rdx
      integer :: i, pp, np, ptrcur

      real(p_k_part), dimension(rank, p_cache_size) :: xbuf
      integer, dimension(rank, p_cache_size)        :: dxi
      real(p_k_part), dimension(p_cache_size)      :: rp, rpp
      real(p_k_part) :: dt

      ! executable statements
      rdx(1) = real(1.0_p_double/this%dx(1), p_k_part)
      rdx(2) = real(1.0_p_double/this%dx(2), p_k_part)
      rdx(3) = real(1.0_p_double/this%dx(3), p_k_part)
      dt = real(gdt, p_k_part)

      ! advance position of photons i0 to i1 in chunks of p_cache_size
      do ptrcur = i0, i1, p_cache_size

         ! check if last copy of table and set np
         if (ptrcur + p_cache_size > i1) then
            np = i1 - ptrcur + 1
         else
            np = p_cache_size
         end if

         ! this type of loop is actually faster than
         ! using a forall construct
         pp = ptrcur

         do i = 1, np
            rp(i) = 1.0_p_k_part/ &
                    sqrt(this%p(1, pp)**2 + this%p(2, pp)**2 + this%p(3, pp)**2)
            rpp(i) = dt*rp(i)
            pp = pp + 1
         end do

         ! reference implementation
         pp = ptrcur
         do i = 1, np
            xbuf(1, i) = this%x(1, pp) + this%p(1, pp)*rpp(i)*rdx(1)
            xbuf(2, i) = this%x(2, pp) + this%p(2, pp)*rpp(i)*rdx(2)
            xbuf(3, i) = this%x(3, pp) + this%p(3, pp)*rpp(i)*rdx(3)
            dxi(1, i) = ntrim_pkt(xbuf(1, i))
            dxi(2, i) = ntrim_pkt(xbuf(2, i))
            dxi(3, i) = ntrim_pkt(xbuf(3, i))

            pp = pp + 1
         end do

         ! copy data from buffer to species data trimming positions
         pp = ptrcur
         do i = 1, np
            this%x(1, pp) = xbuf(1, i) - dxi(1, i)
            this%x(2, pp) = xbuf(2, i) - dxi(2, i)
            this%x(3, pp) = xbuf(3, i) - dxi(3, i)
            this%ix(1, pp) = this%ix(1, pp) + dxi(1, i)
            this%ix(2, pp) = this%ix(2, pp) + dxi(2, i)
            this%ix(3, pp) = this%ix(3, pp) + dxi(3, i)

            pp = pp + 1
         end do
      end do

   end subroutine advance_photons_3d
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
   subroutine advance_photons_cyl_2d(this, den_ph, grad_ph, den_grid, grad_grid, gdt, i0, i1)
!---------------------------------------------------------------------------------------------------

      implicit none

      integer, parameter :: rank = 2

      class(t_photons_pkt), intent(inout) :: this
      real(p_k_fld), intent(inout)           :: den_ph(:)
      real(p_k_fld), intent(inout)           :: grad_ph(:, :)
      type(t_vdf), intent(inout)             :: den_grid, grad_grid
      real(p_double), intent(in) :: gdt
      integer, intent(in) :: i0, i1

      integer :: i, pp, np, ptrcur
      real(p_double), dimension(p_x_dim)            :: xmin_g
      real(p_k_part), dimension(rank, p_cache_size) :: xbuf
      integer, dimension(rank, p_cache_size)        :: dxi

      ! 1/omega-related arrays
      real(p_k_part), dimension(p_cache_size) :: invsq, dt_invsq, invsq_pred
      
      ! Predictor arrays
      real(p_k_fld), dimension(p_cache_size)            :: den_ph_pred
      real(p_k_fld), dimension(p_x_dim, p_cache_size)   :: grad_ph_pred
      real(p_k_part), dimension(rank, p_cache_size)     :: x_pred
      integer, dimension(rank, p_cache_size)            :: ix_pred
      real(p_k_part), dimension(p_cache_size)           :: p1_pred, p2_pred

      ! Other random arrays that are needed for PKT-specific updates
      real(p_k_part) :: k2, omega2, dt, rdx
      integer :: gix2, shift_ix2
      real(p_double) :: dr, rdr, tmp
      real(p_double) :: x2_new, x3_new, r_old, r_new
      real(p_double) :: x2_pred, x3_pred, r_pred
      real(p_double) :: r_shift

      ! -- EXECUTABLE STATEMENTS -------------------------
      xmin_g(1) = this%g_box(p_lower, 1)
      xmin_g(2) = this%g_box(p_lower, 2)

      shift_ix2 = this%my_nx_p(p_lower, 2) - 2
      dr = this%dx(p_r_dim)
      rdr = 1.0_p_double/dr

      rdx = real(gdt/this%dx(1), p_k_part)

      ! Since r = 0 is at the center of cell ix2 = 1 the radial position of particles will be
      ! - Positions defined with regard to the center of the cell (odd interpolation)
      !       r = ( (gix2-1) + x2 ) * dr
      ! - Positions defined with regard to the corner of the cell (even interpolation)
      !       r = ( (gix2-1) + x2 - 0.5 ) * dr

      if (this%pos_type == p_cell_near) then
         r_shift = 0.5_p_double
      else
         r_shift = 0.0_p_double
      end if

      ! -- ADVANCE PARTICLES -----------------------------
      do ptrcur = i0, i1, p_cache_size

         ! Check if there's a last copy of table and set up
         if (ptrcur + p_cache_size > i1) then
            np = i1 - ptrcur + 1
         else
            np = p_cache_size
         end if

         ! Compute 1 / omega for current positions and momenta
         pp = ptrcur
         do i = 1, np
            k2 = this%p(1, pp)**2 + this%p(2, pp)**2 + this%p(3, pp)**2
            omega2 = den_ph(pp)
            invsq(i) = 1.0_p_k_part/sqrt(k2 + omega2)
            dt_invsq(i) = dt*invsq(i)
            pp = pp + 1
         end do

         ! FIRST: use forward Euler as a predictor step
         pp = ptrcur
         do i = 1, np
            ! Predicted momenta --  kicked by both gradient components
            p1_pred(i) = this%p(1, pp) - 0.5_p_k_part*dt_invsq(i)*grad_ph(1, pp)
            p2_pred(i) = this%p(2, pp) - 0.5_p_k_part*dt_invsq(i)*grad_ph(2, pp)

            ! Predicted axial position -- same as 1D case
            x_pred(1, i) = this%x(1, pp) + this%p(1, pp)*dt_invsq(i)*rdx
            ix_pred(1, i) = this%ix(1, pp) + ntrim_pkt(x_pred(1, i))
            x_pred(1, i) = x_pred(1, i) - ntrim_pkt(x_pred(1, i))

            ! Predicted radial position -- use the same conversion as the 
            ! original cylindrical subroutine
            gix2  = this%ix(2, pp) + shift_ix2
            r_old = ( this%x(2, pp) + (gix2 - r_shift) ) * dr

            x2_pred = r_old + ( this%p(2, pp) * dt_invsq(i) )
            x3_pred =         ( this%p(3, pp) * dt_invsq(i) )
            r_pred  = sqrt(x2_pred**2 + x3_pred**2)

            x_pred(2, i) = real(r_pred*rdr - (gix2 - r_shift), p_k_part)
            ix_pred(2, i) = this%ix(2, pp) + ntrim_pkt(x_pred(2, i))
            x_pred(2, i) = x_pred(2, i) - ntrim_pkt(x_pred(2, i))

            pp = pp + 1
         end do

         ! SECOND: Use the predicted positions and momenta to calculate the
         ! updated values for the gradient(s) and densities at the photon 
         ! positions
         ! TODO: need to fix this such that it works for any interpolation order
         call gather_density_2d(den_ph_pred, den_grid, ix_pred, x_pred, np, 2)
         call gather_density_2d(grad_ph_pred(1,:), grad_grid, ix_pred, x_pred, np, 2)
         call gather_density_2d(grad_ph_pred(2,:), grad_grid, ix_pred, x_pred, np, 2, field_comp=2)
         
         ! THIRD: Use the trapezoidal method as the corrector method to update
         ! the momenta of all of the photons.
         pp = ptrcur
         do i = 1, np
            ! Recalculate 1/omega with the updated momenta and density at
            ! photon positions
            k2 = p1_pred(i)**2 + p2_pred(i)**2 + this%p(3, pp)**2
            omega2 = den_ph_pred(i)
            invsq_pred(i) = 1.0_p_k_part/sqrt(k2 + omega2)

            ! Calculate the updated axial momentum (should be the SAME as 1D case)
            this%p(1, pp) = this%p(1, pp) - 0.5_p_k_part*dt &
                            *0.5_p_k_part*(invsq(i)*grad_ph(1, pp) + &
                                           invsq_pred(i)*grad_ph_pred(1, i))

            ! Calculate the updated radial momentum
            this%p(2, pp) = this%p(2, pp) - 0.5_p_k_part*dt &
                            *0.5_p_k_part*(invsq(i)*grad_ph(2, pp) + &
                                           invsq_pred(i)*grad_ph_pred(2, i))

            ! Update the axial position -- this should be identical to 1D
            xbuf(1, i) = this%x(1, pp) + 0.5_p_k_part*(this%p(1, pp)*invsq(i) + &
                                                       p1_pred(i)*invsq_pred(i))*dt*rdx
            dxi(1, i) = ntrim_pkt(xbuf(1, i))

            ! For the radial position, do like the original cylindrical subroutine
            ! and how it's done above in the predictor step
            gix2 = this%ix(2, pp) + shift_ix2
            r_old = (this%x(2, pp) + (gix2 - r_shift))*dr
            
            ! Push in a 3d-like fashion using the trapezoidal average of the 
            ! radial velocity
            x2_new = r_old + 0.5_p_k_part*(this%p(2, pp)*invsq(i) + &
                                           p2_pred(i)*invsq_pred(i)) * gdt
            x3_new = this%p(3, pp) * 0.5_p_k_part * (invsq(i) + invsq_pred(i)) * gdt
            r_new = sqrt(x2_new**2 + x3_new**2)

            ! Use the original cylindrical code and convert back to a cell-
            ! fractional representation.
            if ( r_old == r_new ) then
                xbuf(2, i) = this%x(2, pp)
            else
                xbuf(2, i) = real(r_new * rdr - (gix2 - r_shift), p_k_part)
            endif

            ! Angular momentum conservation -- same as the original cylindrical
            ! subroutine
            tmp          = 1.0_p_double / r_new
            this%p(2,pp) = real((this%p(2,pp) * x2_new + this%p(3,pp)*x3_new ) * tmp, p_k_part)
            this%p(3,pp) = real(this%p(3,pp) * ( r_old * tmp ), p_k_part)

            dxi(2,i) = ntrim_pkt( xbuf(2,i) )

            pp = pp + 1
         end do

         ! copy data from buffer to species data trimming positions
         pp = ptrcur
         do i = 1, np
            this%x(1, pp) = xbuf(1, i) - dxi(1, i)
            this%x(2, pp) = xbuf(2, i) - dxi(2, i)
            this%ix(1, pp) = this%ix(1, pp) + dxi(1, i)
            this%ix(2, pp) = this%ix(2, pp) + dxi(2, i)
            pp = pp + 1
         end do
      end do

   end subroutine advance_photons_cyl_2d
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
   function ntrim_pkt(x)
!---------------------------------------------------------------------------------------------------
! Returns the integer shift (-1, 0 or +1) so that the coordinate remains in the [-0.5, 0.5[
! range. This is the fastest implementation (twice as fast as a sequence of ifs) because
! the two if structures compile as conditional moves and can be processed independently.
! This has no precision problem and is only 12% slower than the previous "int(x+1.5)-1"
! routine that would break for x = nearest( 0.5, -1.0 )
!---------------------------------------------------------------------------------------------------
      implicit none

      real(p_k_part), intent(in) :: x
      integer :: ntrim_pkt, a, b

      if (x < -.5) then
         a = -1
      else
         a = 0
      end if

      if (x >= .5) then
         b = +1
      else
         b = 0
      end if

      ntrim_pkt = a + b

   end function ntrim_pkt
!---------------------------------------------------------------------------------------------------

end module m_photons_pkt

!---------------------------------------------------------------------------------------------------
! advance photons
!---------------------------------------------------------------------------------------------------
subroutine advance_photons_pkt(this, den_ph, grad_ph, den_grid, grad_grid, gdt, tid, n_threads)

   use m_system
   use m_parameters
   use m_photons_define_pkt, only: t_photons_pkt
   use m_species_define, only: p_ene_recalc
   use m_species_define_pkt, only: t_species_pkt
   use m_emf_define, only: t_emf
   use m_photons_pkt
   use m_species_define, only: p_ene_recalc
   use m_species, only: sort
   use m_vdf_define

   implicit none

   class(t_photons_pkt), intent(inout) :: this
   real(p_k_fld), intent(inout) :: den_ph(:)
   real(p_k_fld), intent(inout) :: grad_ph(:, :)
   type(t_vdf), intent(inout)   :: den_grid, grad_grid
   real(p_double), intent(in)   :: gdt
   integer, intent(in)    :: tid        ! local thread id
   integer, intent(in)    :: n_threads  ! total number of threads
   integer :: chunk, ip0, ip1

   ! sanity check
   ! Code does not yet support OpenMP parallelism (mainly because of pair production)
   if (n_threads > 1) then
      write (0, *) 'pkt photon pusher is not yet supported with multiple threads per node'
      call abort_program()
   end if

   ! range of photons for each thread
   chunk = (this%num_par + n_threads - 1)/n_threads
   ip0 = tid*chunk + 1
   ip1 = min((tid + 1)*chunk, this%num_par)

   ! Push photons. Boundary crossings will be checked at update_boundary
   select case (this%coordinates)
   case default
      select case (p_x_dim)
      case (1)
         call advance_photons_1d(this, den_ph, grad_ph, den_grid, grad_grid, gdt, ip0, ip1)

      case (2)
         call advance_photons_2d(this, gdt, ip0, ip1)

      case (3)
         call advance_photons_3d(this, gdt, ip0, ip1)

      case default
         ERROR('Not implemented for x_dim = ', p_x_dim)
         call abort_program(p_err_invalid)
      end select
   case (p_cylindrical_b)
      call advance_photons_cyl_2d(this, den_ph, grad_ph, den_grid, grad_grid, gdt, ip0, ip1)
   end select

end subroutine advance_photons_pkt
!---------------------------------------------------------------------------------------------------
