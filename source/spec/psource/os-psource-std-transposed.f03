#include "os-config.h"
#include "os-preprocess.fpp"

module m_psource_std_transposed

#include "memory/memory.h"

use m_species_define
use m_parameters
use m_current_define
use m_node_conf
use m_psource_std

implicit none

private

type, extends(t_psource_std) :: t_psource_std_transposed

contains

  procedure :: inject => inject_std_transposed

end type t_psource_std_transposed

public :: t_psource_std_transposed

contains

!-------------------------------------------------------------------------------
! transposed version of the standard injection, this is used with algorithm cuda
!-------------------------------------------------------------------------------
function inject_std_transposed( this, species, ig_xbnd_inj, jay, no_co, bnd_cross, node_cross, &
     send_msg, recv_msg ) result(num_inj)

  use m_species_udist
#ifdef __HAS_SPIN__
  use m_species_sdist
#endif

  implicit none

  class(t_psource_std_transposed), intent(inout) :: this
  class(t_species), intent(inout), target :: species
  integer, dimension(:, :), intent(in) :: ig_xbnd_inj
  class( t_current ), intent(inout)   :: jay
  class( t_node_conf ), intent(in)     :: no_co
  type( t_part_idx ), dimension(2), intent(inout) :: bnd_cross
  type( t_part_idx ), intent(inout) :: node_cross
  type( t_spec_msg ), dimension(2), intent(inout) :: send_msg, recv_msg

  integer :: num_inj

  integer, parameter :: p_max_x_dim = 3

  ! cell size (p_x_dim)
  real(p_double), dimension(p_max_x_dim) :: dx
  real(p_double), dimension(p_max_x_dim) :: g_xmin

  ! half distance between particles
  real(p_k_part), dimension(p_max_x_dim) :: dxp_2

  integer :: i, i1_grid, i2_grid, i3_grid, i1_cell, i2_cell, i3_cell, ipart

  real(p_k_part)  :: x1, x2, x3 

  ! volume that each particle occupies
  real(p_k_part) :: pvol

  ! number of particles per cell
  integer :: ppcell

  ! total number of cells to inject
  integer :: ncells_total

  ! particle positions (global / inside cell )
  real(p_k_part), dimension(:,:), pointer :: ppos, pidx, ppos_cell

  ! particle charges
  real(p_k_part), dimension(:), pointer :: pcharge

  ! Set total number of injected particles to 0
  num_inj = 0

  ! check if num_par_x is > 0 for all directions
  ! if not return silently
  do i=1, p_x_dim
     if (species%num_par_x(i) <= 0) return
  enddo

  do i = 1, p_x_dim
     ! get cell size
     dx(i) = species%dx(i)

     ! get global mininum (this is shifted by +0.5 cells from global simulation values)
     g_xmin( i ) = species%g_box( p_lower , i )

     ! get half distance between particles
     dxp_2(i) = 0.5_p_k_part/species%num_par_x(i)
  enddo

  ! find total number particles per cell
  ppcell = species%num_par_x(1)
  do i = 2, p_x_dim
     ppcell = ppcell * species%num_par_x(i)
  enddo

  ! find total number of cells to inject 
  ncells_total = ig_xbnd_inj(p_upper,1) - ig_xbnd_inj(p_lower,1) + 1
  do i = 2, p_x_dim
     ncells_total = ncells_total * ( ig_xbnd_inj(p_upper,i) - ig_xbnd_inj(p_lower,i) + 1 )
  enddo

  ! initialize temp buffers
  call alloc( ppos, (/p_x_dim, ncells_total * ppcell/) )
  call alloc( pidx, (/p_x_dim, ncells_total * ppcell/) )
  call alloc( ppos_cell, (/p_x_dim, ncells_total * ppcell/) )
  call alloc( pcharge, (/ncells_total * ppcell/))

  ! find normalization factor
  pvol = sign( 1.0_p_k_part/ppcell, species%rqm )

  ! Get position of particles inside the cell
  ! The position will always be in the range [-0.5, +0.5 [ regardless of interpolation type

  i = 0
  select case ( p_x_dim )
  case (1)

     do i1_cell = 1, species%num_par_x(1)
        x1 = (2*i1_cell - 1 - species%num_par_x(1)) * dxp_2(1)

        do i1_grid = ig_xbnd_inj(p_lower,1), ig_xbnd_inj(p_upper,1) 
           
           i = i + 1
           
           pidx(1,i) = i1_grid

           ppos_cell(1,i) = x1
                 
           ppos(1,i) = real( g_xmin(1) + (i1_grid-1)*dx(1) + x1, p_k_part )

        enddo
     enddo

  case (2)

     do i1_cell = 1, species%num_par_x(1)
        x1 = (2*i1_cell - 1 - species%num_par_x(1)) * dxp_2(1)

        do i2_cell = 1, species%num_par_x(2)
           x2 = (2*i2_cell - 1 - species%num_par_x(2)) * dxp_2(2)

           do i1_grid = ig_xbnd_inj(p_lower,1), ig_xbnd_inj(p_upper,1) 

              do i2_grid = ig_xbnd_inj(p_lower,2), ig_xbnd_inj(p_upper,2) 

                 i = i + 1

                 pidx(1,i) = i1_grid
                 pidx(2,i) = i2_grid

                 ppos_cell(1,i) = x1
                 ppos_cell(2,i) = x2 

                 ppos(1,i) = real( g_xmin(1) + (i1_grid-1)*dx(1) + x1, p_k_part )
                 ppos(2,i) = real( g_xmin(2) + (i2_grid-1)*dx(2) + x2, p_k_part )

              enddo
           enddo
        enddo
     enddo

  case (3)

     do i1_cell = 1, species%num_par_x(1)
        x1 = (2*i1_cell - 1 - species%num_par_x(1)) * dxp_2(1)

        do i2_cell = 1, species%num_par_x(2)
           x2 = (2*i2_cell - 1 - species%num_par_x(2)) * dxp_2(2)

           do i3_cell = 1, species%num_par_x(3)
              x3 = (2*i3_cell - 1 - species%num_par_x(3)) * dxp_2(3)

              do i1_grid = ig_xbnd_inj(p_lower,1), ig_xbnd_inj(p_upper,1) 

                 do i2_grid = ig_xbnd_inj(p_lower,2), ig_xbnd_inj(p_upper,2) 

                    do i3_grid = ig_xbnd_inj(p_lower,3), ig_xbnd_inj(p_upper,3) 

                       i = i + 1

                       pidx(1,i) = i1_grid
                       pidx(2,i) = i2_grid
                       pidx(3,i) = i3_grid

                       ppos_cell(1,i) = x1
                       ppos_cell(2,i) = x2 
                       ppos_cell(3,i) = x3 

                       ppos(1,i) = real( g_xmin(1) + (i1_grid-1)*dx(1) + x1, p_k_part )
                       ppos(2,i) = real( g_xmin(2) + (i2_grid-1)*dx(2) + x2, p_k_part )
                       ppos(3,i) = real( g_xmin(3) + (i3_grid-1)*dx(3) + x3, p_k_part )

                    enddo
                 enddo
              enddo
           enddo
        enddo
     enddo

  end select

  ! Inject particles
  ipart = species%num_par + 1

  select case ( p_x_dim )

  case (1) ! -----------------------------------------------------------

     call this%get_den_value( ppos, ncells_total * ppcell, pcharge )

     do i = 1, ncells_total * ppcell
        if (pcharge(i) > this % den_min ) then
           num_inj = num_inj + 1
        endif
     enddo

     if ( num_inj > 0 ) then

        ! check if the buffer size is sufficient and grow it if necessary
        if ( num_inj > species%num_par_max - species%num_par ) then
           call species%grow_buffer( species%num_par_max + num_inj + p_spec_buf_block )
        endif

        ! loop through all the injection cells and inject particles
        do i=1, ncells_total * ppcell

           if (pcharge(i) > this % den_min) then

              ! add particle
              species%q(ipart)  = pcharge(i) * pvol
              species%x(1, ipart) = ppos_cell(1,i)
              species%ix(1, ipart) = pidx(1,i) - species%my_nx_p( p_lower, 1 ) + 1

              ipart = ipart + 1

           endif

        enddo

     endif   
  
  case (2) ! -----------------------------------------------------------


     ! loop through all cells and count number of particles to inject
     select case ( species%coordinates )

     case default  ! cartesian coordinates

        call this%get_den_value( ppos, ncells_total * ppcell, pcharge )

        do i = 1, ncells_total * ppcell 
           if (pcharge(i) > this % den_min ) then
              num_inj = num_inj + 1
           endif
        enddo

     case ( p_cylindrical_b ) ! cyl. coord. -> B1 on axis
        ! don't inject on axis

        call this%get_den_value( ppos, ncells_total * ppcell, pcharge )

        do i = 1, ncells_total * ppcell
           if ((pcharge(i) > this % den_min) .and. &
                (ppos(p_r_dim,i) > 0.0_p_k_part)) then

              num_inj = num_inj + 1
           endif
        enddo

     end select

     if ( num_inj > 0 ) then

        ! check if the buffer size is sufficient and grow it if necessary
        if ( num_inj > species%num_par_max - species%num_par ) then
           call species%grow_buffer( species%num_par_max + num_inj + p_spec_buf_block )
        endif


        ! loop through all the injection cells and
        ! inject particles, normalizing charge
        select case ( species%coordinates )
        case default  ! cartesian coordinates

           call this%get_den_value( ppos, ncells_total * ppcell, pcharge )

           do i=1, ncells_total * ppcell
              if (pcharge(i) > this % den_min) then
                 ! add particle
                 species%x(1, ipart) = ppos_cell(1,i)
                 species%x(2, ipart) = ppos_cell(2,i)
                 species%ix(1, ipart) = pidx(1,i) - species%my_nx_p( p_lower, 1 ) + 1
                 species%ix(2, ipart) = pidx(2,i) - species%my_nx_p( p_lower, 2 ) + 1
                 species%q(ipart)  = pcharge(i) * pvol
                 ipart = ipart + 1
              endif

           enddo

        case ( p_cylindrical_b ) ! cyl. coord. -> B1 on axis
           ! don't inject on axis

           call this%get_den_value( ppos, ncells_total * ppcell, pcharge )

           do i=1, ncells_total * ppcell
              if ((pcharge(i) > this % den_min) .and. &
                   (ppos(p_r_dim,i) > 0.0_p_k_part)) then
                 ! add particle
                 species%x(1, ipart) = ppos_cell(1,i)
                 species%x(2, ipart) = ppos_cell(2,i)
                 species%ix(1, ipart) = pidx(1,i) - species%my_nx_p( p_lower, 1 ) + 1
                 species%ix(2, ipart) = pidx(2,i) - species%my_nx_p( p_lower, 2 ) + 1
                 species%q(ipart)  = pcharge(i) * pvol * ppos( p_r_dim, i )
                 ipart = ipart + 1
              endif
           enddo

        end select

     endif


  case (3) ! -----------------------------------------------------------

     call this%get_den_value( ppos, ncells_total * ppcell, pcharge )

     do i = 1, ncells_total * ppcell
        if (pcharge(i) > this % den_min ) then
           num_inj = num_inj + 1
        endif
     enddo

     if ( num_inj > 0 ) then

        ! check if the buffer size is sufficient and grow it if necessary
        if ( num_inj > species%num_par_max - species%num_par ) then
           call species%grow_buffer( species%num_par_max + num_inj + p_spec_buf_block )
        endif

        ! loop through all the injection cells and inject particles
        do i=1, ncells_total * ppcell

           if (pcharge(i) > this % den_min) then

              ! add particle
              species%q(ipart)  = pcharge(i) * pvol

              species%x(1, ipart) = ppos_cell(1,i)
              species%x(2, ipart) = ppos_cell(2,i)
              species%x(3, ipart) = ppos_cell(3,i)

              species%ix(1, ipart) = pidx(1,i) - species%my_nx_p( p_lower, 1 ) + 1
              species%ix(2, ipart) = pidx(2,i) - species%my_nx_p( p_lower, 2 ) + 1
              species%ix(3, ipart) = pidx(3,i) - species%my_nx_p( p_lower, 3 ) + 1

              ipart = ipart + 1

           endif

        enddo

     endif

  end select

  ! print *, "num_inj", num_inj
  
  ! free temporary memory
  call freemem( ppos )
  call freemem( ppos_cell )
  call freemem( pcharge )
  call freemem( pidx )

  ! Set momentum of injected particles
  if ( num_inj > 0 ) then
     ! Initialize particle momentum
     call set_momentum( species, species%num_par+1, species%num_par + num_inj )
#ifdef __HAS_SPIN__
     call set_spin( species, species%num_par+1, species%num_par + num_inj )
#endif
  endif

end function inject_std_transposed
!-------------------------------------------------------------------------------

end module m_psource_std_transposed
