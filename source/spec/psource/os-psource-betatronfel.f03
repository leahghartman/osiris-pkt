#include "os-config.h"
#include "os-preprocess.fpp"

module m_psource_btfel

#include "memory/memory.h"

use m_species_define
use m_parameters
use m_psource_std
use m_current_define
use m_node_conf

implicit none

private

type, extends( t_psource_std ) :: t_psource_betatronfel

  real(p_k_part)                        :: gamma_0
  real(p_k_part)                        :: K_number
  real(p_k_part)                        :: DE_E
  real(p_k_part)                        :: DK_K
  real(p_k_part), dimension(p_p_dim)    :: rotation_center
  logical, dimension(p_p_dim)           :: init_random
  logical                               :: in_boosted_frame

contains

  procedure :: inject => inject_betatronfel
  procedure :: read_input => read_input_betatronfel

end type t_psource_betatronfel

public :: t_psource_betatronfel

contains

!---------------------------------------------------------------------------------------------------
function inject_betatronfel( this, species, ig_xbnd_inj, jay, no_co, bnd_cross, node_cross, &
                             send_msg, recv_msg ) result(num_inj)
!---------------------------------------------------------------------------------------------------
! X.D. 26 July
! This routine injects particles for the Betatron FEL psource
! usefull (for the study of ICL)
!---------------------------------------------------------------------------------------------------
  use m_random
  use m_species_tag

  implicit none

  class(t_psource_betatronfel), intent(inout) :: this
  class(t_species), intent(inout), target :: species
  integer, dimension(:,:), intent(in) :: ig_xbnd_inj
  class( t_current ), intent(inout)   :: jay
  class( t_node_conf ), intent(in)     :: no_co
  type( t_part_idx ), dimension(2), intent(inout) :: bnd_cross
  type( t_part_idx ), intent(inout) :: node_cross
  type( t_spec_msg ), dimension(2), intent(inout) :: send_msg, recv_msg

  integer :: num_inj

  integer, parameter :: p_max_x_dim = 3

  ! cell size (p_x_dim)
  real(p_double), dimension(p_max_x_dim) :: dx
  real(p_double), dimension(2, p_max_x_dim) :: g_xbnd

  ! increment in position for the particles
  real(p_double), dimension(p_max_x_dim) :: dxpart
  ! number of particles per cell for each direction
  integer, dimension(p_max_x_dim) :: npx

  integer :: i, ierr, i1, i2, i3, ipart

  ! volume that each particle occupies
  real(p_k_part) :: pvol

  ! number of particles per cell
  integer :: ppcell

  ! particle positions (global / inside cell )
  real(p_k_part), dimension(:,:), allocatable :: ppos, ppos_cell

  ! particle charges
  real(p_k_part), dimension(:), allocatable :: pcharge

  real(p_double) :: dt

  ! betatron parameters
  real(p_k_part)  :: g0, K, B_beam, g_beam, r0, DE_E, DK_K, dk_k_res, de_e_res, rnd
  real(p_k_part)  :: pos_x1, pos_x2, pos_x3, pos_xr, p1_lab, p2_lab, p3_lab, pr_lab
  real(p_k_part)  :: center_x2, center_x3, theta_r, theta_t, ssin, gamma, d_theta_r
  real(p_k_part)  :: pos_x2_min, pos_x2_max, pos_x3_min, pos_x3_max
  logical, dimension(p_max_x_dim) :: init_random
  integer, dimension(p_max_x_dim) :: pos_ix 
  logical :: in_boosted_frame
  integer :: offset

  integer :: max_up_bound

  real(p_k_part), parameter :: pi = 3.14159265358979323846264_p_k_part
  ! executable statements

  ! THIS IS A HACK PUT IN BY ADAM
  ! TODO: make this real
  !
  ! Calculate the integer portion of each particle to pass to the random number generator
  !   so that consistent random numbers are generated no matter what spatial distibution is used.
  !   - 'pos_ix' is a place holder for now so that the code will compile
  !   - until we fix things, this code, if run, will prevent random numbers from being generated consistenly
  !       across all spatial decompositions. TODO: We should alert to warn the user of this fact (but for now
  !       this code isn't even called at all so no rush)
  !call set_rng_constancy_is_spoiled()
  pos_ix = 1
  
  dt = species%dt

  num_inj = 0

  ! check if num_par_x is > 0 for all directions
  ! if not return silently
  do i=1, p_x_dim
    if (species%num_par_x(i) <= 0) return
  enddo

  do i = 1, p_x_dim
    ! get cell size
    dx(i) = species%dx(i)
    g_xbnd( :, i ) = species%g_box( :, i )
    ! get distance between particles
    dxpart(i) = 1.0_p_double/species%num_par_x(i)
    npx(i) = species%num_par_x(i)
    init_random(i) = this%init_random(i)
  enddo

  ! betatron parameters
  g0 = this%gamma_0
  K = this%K_number
  DE_E = this%DE_E
  DK_K = this%DK_K
  if ( p_x_dim > 1 ) center_x2 = this%rotation_center(1)
  if ( p_x_dim == 3 ) center_x3 = this%rotation_center(2)
  in_boosted_frame = this%in_boosted_frame

  B_beam = 1.0_p_k_part-(2.0_p_k_part+K**2)/(4.0_p_k_part*g0**2)
  g_beam = 1.0_p_k_part/sqrt(1.0_p_k_part-B_beam**2)
  r0 = (4.0_p_k_part*g0*K)/(2.0_p_k_part+K**2)
  if (in_boosted_frame) then
     r0 = r0 / (g_beam*(1+B_beam))
  endif

  ipart = species%num_par + 1

  select case ( p_x_dim )

     case (1) ! -----------------------------------------------------------

       print *, "1D geometry cannot be used with vdist_type = 'betatronfel'."
       print *, "Aborting..."
       stop

     case (2) ! -----------------------------------------------------------

       ! limit min et max for x2 and theta, depending on the CPU boundary
       pos_x2_min = real( g_xbnd(p_lower, 2) + (species%my_nx_p(p_lower,2)-1)*dx(2) , p_k_part )
       pos_x2_max = real( g_xbnd(p_lower, 2) + species%my_nx_p(p_upper,2)*dx(2) , p_k_part )

       ! find normalization factor
       pvol = sign( real((2*r0/dx(2))/(npx(1)*npx(2)), p_k_part), species%rqm )

       ! find ppcell
       ppcell=0
       theta_r = real( -pi , p_k_part )
       do i2 = 1, npx(2)
          pos_x2 = center_x2 + r0 * cos( theta_r )
          if ( (pos_x2.ge.pos_x2_min).and.(pos_x2<pos_x2_max) ) then
             ppcell = ppcell + 1
          endif
          p2_lab = K * sin( theta_r )
          gamma = sqrt( g0**2 + p2_lab**2 )
          if (in_boosted_frame) then
             gamma = gamma*g_beam - sqrt(g0**2-1.0_p_k_part)*g_beam*B_beam
          endif
          d_theta_r = real( 2.0*pi*K / ( npx(2)*r0*gamma ) , p_k_part )
          theta_r = theta_r + d_theta_r
       enddo
       ppcell = ppcell * npx(1)

       if ( ppcell == 0 ) return

       ! initialize temp buffers
       allocate( ppos(p_x_dim, ppcell), &
                 ppos_cell(p_x_dim, ppcell), &
                 pcharge(ppcell), stat = ierr)
       if ( ierr/=0 ) then
          ERROR('Allocation Failed!')
          call abort_program( p_err_alloc )
       endif

       i=0
       theta_r = real( -pi , p_k_part )
       do i2 = 1, npx(2)
          pos_x2 = center_x2 + r0 * cos( theta_r )
          if ( (pos_x2.ge.pos_x2_min).and.(pos_x2<pos_x2_max) ) then
             do i1 = 1, npx(1)
                i = i + 1
                ppos_cell(1,i) = 0.5_p_k_part
                ppos_cell(2,i) = theta_r
             enddo
          endif
          p2_lab = K * sin( theta_r )
          gamma = sqrt( g0**2 + p2_lab**2 )
          if (in_boosted_frame) then
             gamma = gamma*g_beam - sqrt(g0**2-1.0_p_k_part)*g_beam*B_beam
          endif
          d_theta_r = real( 2.0*pi*K / ( npx(2)*r0*gamma ) , p_k_part )
          theta_r = theta_r + d_theta_r
       enddo

       do i1 = ig_xbnd_inj(p_lower,1), ig_xbnd_inj(p_upper,1)
          ppos(1,:) = real( g_xbnd(p_lower, 1) + (ppos_cell(1,:) + (i1-1))*dx(1), p_k_part )
          ppos(2,:) = ppos_cell(2,:)
          call this%get_den_value( ppos, ppcell, pcharge )
          do i = 1, ppcell
             if (pcharge(i) >= this%den_min ) then
                num_inj = num_inj + 1
             endif
          enddo
       enddo

       if ( num_inj > 0 ) then

          ! check if the buffer size is sufficient and grow it if necessary
          if ( num_inj > species%num_par_max - species%num_par ) then
             call species % grow_buffer( species%num_par_max + num_inj + p_spec_buf_block )
          endif


          ! loop through all the injection cells and
          ! inject particles, normalizing charge
          ! cartesian coordinates

           ! offset needed to calculate the integer ix when ix<0
           offset = int( max( (species%g_box(p_upper,1)-species%g_box(p_lower,1))/dx(1) , (species%g_box(p_upper,2)-species%g_box(p_lower,2))/dx(2) ) ) + 1

           do i1 = ig_xbnd_inj(p_lower,1), ig_xbnd_inj(p_upper,1)
              ppos(1,:) = real( g_xbnd(p_lower, 1) + (ppos_cell(1,:) + (i1-1))*dx(1), p_k_part )
              ppos(2,:) = ppos_cell(2,:)

              call this%get_den_value( ppos, ppcell, pcharge )

              do i = 1, ppcell
                 if (pcharge(i) >= this%den_min ) then

                    ! definition of x1 with spread
                    pos_x1 = ppos(1,i)
                    if (init_random(1)) then
                       call rng%harvest_real3( rnd, pos_ix )
                       pos_x1 = pos_x1 - rnd * real(dx(1), p_k_part)
                    endif

                    ! definition of p1 with spread in the lab frame
                    de_e_res = DE_E*real( rng%genrand_gaussian( pos_ix ), p_k_part )
                    p1_lab = sqrt(g0**2-1.0_p_k_part) * ( 1.0_p_k_part + de_e_res )

                    ! definition of x2 and p2 with spread
                    theta_r = ppos(2,i)
                    if (init_random(2)) then
                       p2_lab = K * sin( theta_r )
                       gamma = sqrt( g0**2 + p2_lab**2 )
                       if (in_boosted_frame) then
                          gamma = gamma*g_beam - sqrt(g0**2-1.0_p_k_part)*g_beam*B_beam
                       endif
                       d_theta_r = real( 2.0*pi*K / ( npx(2)*r0*gamma ) , p_k_part )
                       call rng%harvest_real3( rnd, pos_ix )
                       theta_r = theta_r + ( rnd - 0.5_p_k_part ) * d_theta_r
                    endif
                    dk_k_res = DK_K*real( rng%genrand_gaussian( pos_ix ), p_k_part )
                    pos_x2 = center_x2 + r0 * ( 1.0_p_k_part + dk_k_res ) * cos( theta_r )
                    ssin = sin( theta_r )
                    p2_lab = K * ssin * sqrt( 1 + (K*ssin/(2.0_p_k_part*g0))**2 ) * ( 1.0_p_k_part + dk_k_res )

                    ! definition of p1 in the boosted frame if necessary + initialisation of p1 and p2 at the time step n=-1/2
                    if (in_boosted_frame) then
                       p1_lab = g_beam*( p1_lab - B_beam*sqrt(1.0_p_k_part+p1_lab**2+p2_lab**2) )
                    endif
                    ssin = real( sin( theta_r + dt/2.0_p_k_part * K/(r0*sqrt(1+p1_lab**2+p2_lab**2)) ), p_k_part )
                    p2_lab = K * ssin * sqrt( 1 + (K*ssin/(2.0_p_k_part*g0))**2 ) * ( 1.0_p_k_part + dk_k_res )
                    p1_lab = sqrt(g0**2-1.0_p_k_part) * ( 1.0_p_k_part + de_e_res )
                    if (in_boosted_frame) then
                       p1_lab = g_beam*( p1_lab - B_beam*sqrt(1.0_p_k_part+p1_lab**2+p2_lab**2) )
                    endif

                    ! add particle
                    species%q(ipart)  = pcharge(i) * pvol
                    species%ix(1, ipart) = int( (pos_x1 - species%g_box(p_lower,1)) / dx(1) + 2.5_p_k_part - species%my_nx_p(p_lower,1) + offset ) - offset
                    species% x(1, ipart) = (pos_x1 - real(species%g_box(p_lower,1),p_k_part)) / real(dx(1),p_k_part) + 2.0_p_k_part - species%my_nx_p(p_lower,1) - species%ix(1,ipart)
                    species%ix(2, ipart) = int( (pos_x2 - species%g_box(p_lower,2)) / dx(2) + 2.5_p_k_part - species%my_nx_p(p_lower,2) + offset ) - offset
                    species% x(2, ipart) = (pos_x2 - real(species%g_box(p_lower,2), p_k_part)) / real(dx(2),p_k_part) + 2.0_p_k_part - species%my_nx_p(p_lower,2) - species%ix(2,ipart)
                    species%p(1, ipart) = p1_lab
                    species%p(2, ipart) = p2_lab
                    if ( p_p_dim == 3 ) then
                       species%p(3, ipart) = 0.0_p_k_part
                    endif

                    ipart = ipart + 1

                 endif
              enddo
           enddo

       endif


     case (3) ! -----------------------------------------------------------

       ! limit min et max for x2 and x3, depending on the CPU boundary
       pos_x2_min = real( g_xbnd(p_lower, 2) + (species%my_nx_p(p_lower,2)-1)*dx(2) , p_k_part )
       pos_x2_max = real( g_xbnd(p_lower, 2) + species%my_nx_p(p_upper,2)*dx(2) , p_k_part )
       pos_x3_min = real( g_xbnd(p_lower, 3) + (species%my_nx_p(p_lower,3)-1)*dx(3) , p_k_part )
       pos_x3_max = real( g_xbnd(p_lower, 3) + species%my_nx_p(p_upper,3)*dx(3) , p_k_part )

       ! find normalization factor
       pvol = sign( real(( pi*r0**2/(dx(2)*dx(3)) ) / ( npx(1)*npx(2)*npx(3) ), p_k_part), species%rqm )

       ! find ppcell
       ppcell=0
       theta_r = real( -pi , p_k_part )      ! theta_r, the angle of the part. in the xrpr plane
       do i2 = 1, npx(2)
          pos_xr = r0 * cos( theta_r )
          do i3 = 1, npx(3)
             theta_t = real( -pi/2.0 + i3*pi/npx(3) , p_k_part )       ! theta_t, the angle of the part. in the x2x3 plane
             pos_x2 = pos_xr * cos( theta_t ) + center_x2
             pos_x3 = pos_xr * sin( theta_t ) + center_x3
             if ( (pos_x2.ge.pos_x2_min).and.(pos_x2<pos_x2_max).and.(pos_x3.ge.pos_x3_min).and.(pos_x3<pos_x3_max) ) then
                ppcell = ppcell + 1
             endif
          enddo
          pr_lab = K * sin( theta_r )
          gamma = sqrt( g0**2 + pr_lab**2 )
          if (in_boosted_frame) then
             gamma = gamma*g_beam - sqrt(g0**2-1.0_p_k_part)*g_beam*B_beam
          endif
          d_theta_r = real( 2.0*pi*K / ( npx(2)*r0*gamma ) , p_k_part )
          theta_r = theta_r + d_theta_r
       enddo
       ppcell = ppcell * npx(1)

       if ( ppcell == 0 ) return

       ! initialize temp buffers
       allocate( ppos(p_x_dim, ppcell), &
                 ppos_cell(p_x_dim, ppcell), &
                 pcharge(ppcell), stat = ierr)
       if ( ierr/=0 ) then
          ERROR('Allocation Failed!')
          call abort_program( p_err_alloc )
       endif

       i=0
       theta_r = real( -pi , p_k_part )      ! theta_r, the angle of the part. in the xrpr plane
       do i2 = 1, npx(2)
          pos_xr = r0 * cos( theta_r )
          do i3 = 1, npx(3)
             theta_t = real( -pi/2.0 + i3*pi/npx(3) , p_k_part )       ! theta_t, the angle of the part. in the x2x3 plane
             pos_x2 = pos_xr * cos( theta_t ) + center_x2
             pos_x3 = pos_xr * sin( theta_t ) + center_x3
             if ( (pos_x2.ge.pos_x2_min).and.(pos_x2<pos_x2_max).and.(pos_x3.ge.pos_x3_min).and.(pos_x3<pos_x3_max) ) then
                do i1 = 1, npx(1)
                   i = i + 1
                   ppos_cell(1,i) = 0.5_p_k_part
                   ppos_cell(2,i) = theta_r
                   ppos_cell(3,i) = theta_t
                enddo
             endif
          enddo
          pr_lab = K * sin( theta_r )
          gamma = sqrt( g0**2 + pr_lab**2 )
          if (in_boosted_frame) then
             gamma = gamma*g_beam - sqrt(g0**2-1.0_p_k_part)*g_beam*B_beam
          endif
          d_theta_r = real( 2.0*pi*K / ( npx(2)*r0*gamma ) , p_k_part )
          theta_r = theta_r + d_theta_r
       enddo

       ! loop through all cells and count number of particles to inject
       do i1 = ig_xbnd_inj(p_lower,1), ig_xbnd_inj(p_upper,1)
          ppos(1,:) = real( g_xbnd(p_lower, 1) + (ppos_cell(1,:) + (i1-1))*dx(1), p_k_part )
          ppos(2,:) = ppos_cell(2,:)
          ppos(3,:) = ppos_cell(3,:)
          call this%get_den_value( ppos, ppcell, pcharge )
          do i = 1, ppcell
             if (pcharge(i) >= this%den_min ) then
                num_inj = num_inj + 1
             endif
          enddo
       enddo

       if ( num_inj > 0 ) then

          ! check if the buffer size is sufficient and grow it if necessary
          if ( num_inj > species%num_par_max - species%num_par ) then
             call species % grow_buffer( species%num_par_max + num_inj + p_spec_buf_block )
          endif


          ! loop through all the injection cells and
          ! inject particles, normalizing charge
          ! cartesian coordinates

           ! offset needed to calculate the integer ix when ix<0
           offset = int( max( (species%g_box(p_upper,1)-species%g_box(p_lower,1))/dx(1) , (species%g_box(p_upper,2)-species%g_box(p_lower,2))/dx(2) , (species%g_box(p_upper,3)-species%g_box(p_lower,3))/dx(3) ) ) + 1

           do i1 = ig_xbnd_inj(p_lower,1), ig_xbnd_inj(p_upper,1)
              ppos(1,:) = real( g_xbnd(p_lower, 1) + (ppos_cell(1,:) + (i1-1))*dx(1), p_k_part )
              ppos(2,:) = ppos_cell(2,:)
              ppos(3,:) = ppos_cell(3,:)

              call this%get_den_value( ppos, ppcell, pcharge )

              do i = 1, ppcell
                 if (pcharge(i) >= this%den_min ) then

                    ! definition of x1 with spread
                    pos_x1 = ppos(1,i)
                    if (init_random(1)) then
                       call rng%harvest_real3( rnd, pos_ix )
                       pos_x1 = pos_x1 - real( rnd * dx(1), p_k_part )
                    endif

                    ! definition of p1 with spread in the lab frame
                    de_e_res = DE_E*real( rng%genrand_gaussian( pos_ix ), p_k_part )
                    p1_lab = sqrt(g0**2-1.0_p_k_part) * ( 1.0_p_k_part + de_e_res )

                    ! definition of xr and pr with spread
                    theta_r = ppos(2,i)
                    if (init_random(2)) then
                       pr_lab = K * sin( theta_r )
                       gamma = sqrt( g0**2 + pr_lab**2 )
                       if (in_boosted_frame) then
                          gamma = gamma*g_beam - sqrt(g0**2-1.0_p_k_part)*g_beam*B_beam
                       endif
                       d_theta_r = real( 2.0*pi*K / ( npx(2)*r0*gamma ) , p_k_part )
                       call rng%harvest_real3( rnd, pos_ix )
                       theta_r = theta_r + rnd * d_theta_r
                    endif
                    dk_k_res = DK_K*real( rng%genrand_gaussian( pos_ix ), p_k_part )
                    pos_xr = r0 * ( 1.0_p_k_part + dk_k_res ) * cos( theta_r )
                    ssin = sin( theta_r )
                    pr_lab = K * ssin * sqrt( 1 + (K*ssin/(2.0_p_k_part*g0))**2 ) * ( 1.0_p_k_part + dk_k_res )

                    ! definition of p1 in the boosted frame if necessary + definition of p1 and pr at the time step n=-1/2
                    if (in_boosted_frame) then
                       p1_lab = g_beam*( p1_lab - B_beam*sqrt(1.0_p_k_part+p1_lab**2+pr_lab**2) )
                    endif
                    ssin = real(sin( theta_r + dt/2.0_p_k_part * K/(r0*sqrt(1+p1_lab**2+pr_lab**2)) ), p_k_part)
                    pr_lab = K * ssin * sqrt( 1 + (K*ssin/(2.0_p_k_part*g0))**2 ) * ( 1.0_p_k_part + dk_k_res )
                    p1_lab = sqrt(g0**2-1.0_p_k_part) * ( 1.0_p_k_part + de_e_res )
                    if (in_boosted_frame) then
                       p1_lab = g_beam*( p1_lab - B_beam*sqrt(1.0_p_k_part+p1_lab**2+pr_lab**2) )
                    endif

                    ! definition of x2, x3 ,p2 ,p3
                    theta_t = ppos(3,i)
                    if (init_random(3)) then
                       call rng%harvest_real3( rnd, pos_ix )
                       rnd = ( rnd - 0.5_p_k_part ) * pi/npx(3)
                       theta_t = theta_t + rnd
                    endif
                    pos_x2 = center_x2 + pos_xr * cos( theta_t )
                    pos_x3 = center_x3 + pos_xr * sin( theta_t )
                    p2_lab = pr_lab * cos( theta_t )
                    p3_lab = pr_lab * sin( theta_t )

                    ! add particle
                    species%q(ipart)  = pcharge(i) * pvol
                    species%ix(1, ipart) = int( (pos_x1 - species%g_box(p_lower,1)) / dx(1) + 2.5_p_k_part - species%my_nx_p(p_lower,1) + offset ) - offset
                    species% x(1, ipart) = (pos_x1 - species%g_box(p_lower,1)) / dx(1) + 2.0_p_k_part - species%my_nx_p(p_lower,1) - species%ix(1,ipart)
                    species%ix(2, ipart) = int( (pos_x2 - species%g_box(p_lower,2)) / dx(2) + 2.5_p_k_part - species%my_nx_p(p_lower,2) + offset ) - offset
                    species% x(2, ipart) = (pos_x2 - species%g_box(p_lower,2)) / dx(2) + 2.0_p_k_part - species%my_nx_p(p_lower,2) - species%ix(2,ipart)
                    species%ix(3, ipart) = int( (pos_x3 - species%g_box(p_lower,3)) / dx(3) + 2.5_p_k_part - species%my_nx_p(p_lower,3) + offset ) - offset
                    species% x(3, ipart) = (pos_x3 - species%g_box(p_lower,3)) / dx(3) + 2.0_p_k_part - species%my_nx_p(p_lower,3) - species%ix(3,ipart)
                    species%p(1, ipart) = p1_lab
                    species%p(2, ipart) = p2_lab
                    species%p(3, ipart) = p3_lab

                    ipart = ipart + 1

                 endif
              enddo
           enddo

       endif

  end select

  deallocate( ppos, ppos_cell, pcharge, stat = ierr )
  if ( ierr /= 0 ) then
    ERROR('Deallocation failed')
    call abort_program( p_err_dealloc )
  endif

  if ( num_inj > 0 ) then

    ! if required set tags of particles
    if (species%add_tag) then
      call set_tags( species, species%num_par+1, species%num_par + num_inj )
    endif

    ! increase num_created counter
    species%num_created = species%num_created + num_inj
    species%num_par = species%num_par + num_inj

    ! Some particles may have been injected outside of local boundaries
    max_up_bound = nx( no_co, 2 )
    if (p_x_dim==3) max_up_bound = max( max_up_bound , nx( no_co, 3 ) )
    do i=1, max_up_bound
      call species % update_boundary( jay, no_co, dt, bnd_cross, node_cross, send_msg, recv_msg )
    enddo

    ! Reset the number of injected particles so no further processing occurs
    num_inj = 0

  endif


end function inject_betatronfel
!-------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! Read information from input file
!---------------------------------------------------------------------------------------------------
subroutine read_input_betatronfel( this, input_file, coordinates )

  use m_input_file

  implicit none

  class( t_psource_betatronfel ), intent(inout) :: this
  class( t_input_file ), intent(inout) :: input_file
  integer, intent(in) :: coordinates

  real(p_k_part)                        :: gamma_0
  real(p_k_part)                        :: K_number
  real(p_k_part)                        :: DE_E
  real(p_k_part)                        :: DK_K
  real(p_k_part), dimension(p_p_dim)    :: rotation_center
  logical, dimension(p_p_dim)           :: init_random
  logical                               :: in_boosted_frame

  integer :: ierr

  namelist /nl_betatron_fel/ gamma_0, K_number, DE_E, DK_K, rotation_center, &
						     init_random, in_boosted_frame

  ! Read betatron FEL specific parameters
  gamma_0 = 50.0              ! X.D. 29 July
  K_number = 1.0              ! X.D. 29 July
  DE_E = 0.0                  ! X.D. 29 July
  DK_K = 0.0                  ! X.D. 29 July
  rotation_center = 0         ! X.D. 29 July
  init_random = .true.        ! X.D. 29 July
  in_boosted_frame = .true.   ! X.D. 29 July

  ! Get namelist text from input file
  call get_namelist( input_file, "nl_betatron_fel", ierr )

  if ( ierr == 0 ) then
    read (input_file%nml_text, nml = nl_betatron_fel, iostat = ierr)
    if (ierr /= 0) then
      if ( mpi_node() == 0 ) then
        write(0,*) ""
        write(0,*) "   Error reading betatron_fel parameters"
        write(0,*) "   aborting..."
      endif
      stop
    endif
  else
    if ( mpi_node() == 0 ) then
      write(0,*) ""
      write(0,*) "   Error, betatron_fel parameters missing"
      write(0,*) "   aborting..."
    endif
    stop
  endif

  this%gamma_0 = gamma_0                    ! X.D. 29 July
  this%K_number = K_number                  ! X.D. 29 July
  this%DE_E = DE_E                          ! X.D. 29 July
  this%DK_K = DK_K                          ! X.D. 29 July
  this%rotation_center = rotation_center    ! X.D. 29 July
  this%init_random = init_random            ! X.D. 29 July
  this%in_boosted_frame = in_boosted_frame  ! X.D. 29 July

  ! Read standard density psource parameters
  call this % t_psource_std % read_input( input_file, coordinates )

end subroutine read_input_betatronfel
!---------------------------------------------------------------------------------------------------

end module m_psource_btfel
