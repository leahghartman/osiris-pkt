#include "os-preprocess.fpp"
#include "os-config.h"

module m_mergedel

#include "memory/memory.h"

use m_parameters
use m_species_define
use m_system
use m_node_conf
use m_species

use m_math

implicit none

private


type( t_part_idx ), save, private :: delidx

!-------------------------------------------------------------------------------
! Merging class
!-------------------------------------------------------------------------------

type :: t_merge_info

  integer :: n_check_max_sim_particles
  integer :: max_sim_particles
  integer, dimension(p_x_dim) :: nx_merge_cells
  integer, dimension(p_p_dim) :: np_per_dir

  real( p_k_part ) :: max_weight

contains
  procedure :: read_input => read_input_merge_info

end type

type :: t_merge
  integer, dimension(p_x_dim) :: nx_merge_cells                     ! number of pic cells within one merging cell
  integer, dimension(p_p_dim) :: np_per_dir                         ! number of cells in momentum space for each direction
  integer, dimension(:), pointer :: mcell_i => null()               ! first index of this merge cell in species
  integer, dimension(:), pointer:: mcell_numpar => null()           ! total number of particles in this merge cell
  real(p_k_part), dimension(:,:), pointer :: pmin => null()         ! min p in each direction for each merge cell
  real(p_k_part), dimension(:,:), pointer :: pmax => null()         ! max p in each direction for each merge cell
  integer, dimension(:), pointer :: part_pspace => null()           ! index of momentum cell each particle belongs to
  real( p_k_part ), dimension(:,:), pointer :: ptot => null()       ! total momentum within a momentum cell
  real( p_k_part ), dimension(:), pointer :: enetot => null()       ! total energy within a momentum cell
  real( p_k_part ), dimension(:), pointer :: weighttot => null()    ! total particle weight within a momentum cell
  integer, dimension(:), pointer :: numpar_per_pc => null()         ! number of particles per momemtum cell
  real( p_k_part ), dimension(:,:), pointer :: pfirst => null()     ! first final momentum within a momentum cell
  real( p_k_part ), dimension(:), pointer :: enefirst => null()     ! first final total energy within a momentum cell
  real( p_k_part ), dimension(:), pointer :: weightfirst => null()  ! first final particle weight within a momentum cell
  real( p_k_part ), dimension(:,:), pointer :: psecond => null()    ! first final momentum within a momentum cell
  real( p_k_part ), dimension(:), pointer :: enesecond => null()    ! first final total energy within a momentum cell
  real( p_k_part ), dimension(:), pointer :: weightsecond => null() ! first final particle weight within a momentum cell
  integer, dimension(:), pointer :: ifirst  => null()
  integer, dimension(:), pointer :: isecond  => null()

contains
  procedure :: init => init_merge
  procedure :: cleanup => cleanup_merge

end type t_merge

interface mergedel
  module procedure merge_particles
end interface

integer, parameter :: p_charge = 1
integer, parameter :: p_photon = 2

public :: t_merge_info, mergedel, p_charge, p_photon, sort_merge_species

contains

subroutine read_input_merge_info( this, input_file )

  use m_input_file, only : t_input_file, get_namelist

  implicit none

  class( t_merge_info ), intent(inout) :: this
  class( t_input_file ), intent(inout) :: input_file

  integer :: n_check_max_sim_particles
  integer :: max_sim_particles
  integer, dimension(p_x_dim) :: nx_merge_cells
  integer, dimension(p_p_dim) :: np_per_dir
  real( p_k_part ) :: max_weight

  integer :: i, ierr

  namelist /nl_merge_info/ &
  n_check_max_sim_particles,  max_sim_particles, nx_merge_cells, np_per_dir, max_weight

  ! default values
  max_sim_particles = -1
  n_check_max_sim_particles = 0
  nx_merge_cells = -1
  np_per_dir     = -1

  max_weight = 1.e4


  ! Get namelist text from input file
  call get_namelist( input_file, "nl_merge_info", ierr )

  if ( ierr == 0 ) then
    read (input_file%nml_text, nml = nl_merge_info, iostat = ierr)
    if (ierr /= 0) then
      print *, ""
      print *, "   Error reading merge_info parameters "
      print *, "   aborting..."
      stop
    endif
  else
    SCR_ROOT("   - no particle merging specified")
  endif

  ! particle deletion
  this%max_sim_particles = max_sim_particles
  if ( max_sim_particles > 0 ) then
    if ( n_check_max_sim_particles < 1 ) then
      print *, "   Error reading merge_info parameters"
      print *, "   When setting max_sim_particles > 0, n_check_max_sim_particles must be >= 1"
      print *, "   aborting..."
      stop
    endif
    this%n_check_max_sim_particles = n_check_max_sim_particles
  else
    this%n_check_max_sim_particles = 0
  endif

  if ( this%n_check_max_sim_particles > 0 ) then
    do i = 1, p_x_dim
      if ( nx_merge_cells(i) <= 0 ) then
        print *, "   Error reading merge_info parameters"
        print *, "   Size of the merging cell must be >= 1 in all directions:"
        print "(A,I0,A,I0)", "    -> nx_merge_cells(",i,") = ", nx_merge_cells(i)
        print *, "   aborting..."
        stop
      endif
      this%nx_merge_cells(i) = nx_merge_cells(i)
    enddo

    do i = 1, p_p_dim
      if ( np_per_dir(i) <= 0 ) then
        print *, "   Error reading merge_info parameters"
        print *, "   Size of the momentum cell must be >= 1 in all directions:"
        print "(A,I0,A,I0)", "    -> np_per_dir(",i,") = ", np_per_dir(i)
        print *, "   aborting..."
        stop
      endif
      this%np_per_dir(i) = np_per_dir(i)
    enddo
  endif

  if ( max_weight <= 0) then
    print *, "   Error reading merge_info parameters"
    print *, "   When specified, max_weight must be > 0"
    print *, "   aborting..."
    stop
  endif
  this % max_weight = max_weight

end subroutine


subroutine init_merge( this, num_merge_cells, num_par_max_merge_cell, npc_total )
  implicit none

  class( t_merge ), intent(inout) :: this
  integer, intent(in) :: num_merge_cells
  integer, intent(in) :: num_par_max_merge_cell
  integer, intent(in) :: npc_total

  call alloc( this%pmin, (/ p_p_dim, num_merge_cells /))
  call alloc( this%pmax, (/ p_p_dim, num_merge_cells /))
  call alloc( this%mcell_i, (/ num_merge_cells /) )
  call alloc( this%mcell_numpar, (/ num_merge_cells /) )

  ! this is the only one that varies from call to call
  call alloc( this%part_pspace, (/ num_par_max_merge_cell /) )

  ! I switched the indexes here to keep all momentum components close in memory
  call alloc( this%ptot, (/ p_p_dim, npc_total /) )
  call alloc( this%pfirst, (/ p_p_dim, npc_total /) )
  call alloc( this%psecond, (/ p_p_dim, npc_total /) )

  call alloc( this%enetot, (/ npc_total /)  )
  call alloc( this%weighttot, (/ npc_total /)  )
  call alloc( this%numpar_per_pc, (/ npc_total /) )
  call alloc( this%enefirst, (/ npc_total /)  )
  call alloc( this%weightfirst, (/ npc_total /)  )
  call alloc( this%enesecond, (/ npc_total /)  )
  call alloc( this%weightsecond, (/ npc_total /)  )
  call alloc( this%ifirst, (/ npc_total /)  )
  call alloc( this%isecond, (/ npc_total /)  )

end subroutine init_merge

subroutine cleanup_merge( this )
  implicit none

  class( t_merge ), intent(inout) :: this

  call freemem( this%pmin )
  call freemem( this%pmax )
  call freemem( this%mcell_i )
  call freemem( this%mcell_numpar )
  call freemem( this%part_pspace )
  call freemem( this%ptot )
  call freemem( this%enetot )
  call freemem( this%weighttot )
  call freemem( this%numpar_per_pc )
  call freemem( this%pfirst )
  call freemem( this%enefirst )
  call freemem( this%weightfirst )
  call freemem( this%psecond )
  call freemem( this%enesecond )
  call freemem( this%weightsecond )
  call freemem( this%ifirst )
  call freemem( this%isecond )

end subroutine cleanup_merge


!-----------------------------------------------------------------------------------------
! Merge/delete excess particles that are close in phase space (Marija & Thomas)
! This works for both charges and photons
!-----------------------------------------------------------------------------------------
subroutine merge_particles( species, merge_info, n, nxc, t, no_co, part_type )

  use m_species_comm

  implicit none

  class( t_species ), intent(inout) :: species
  class( t_merge_info), intent(in) :: merge_info
  integer, intent(in) :: n
  integer,        dimension(:), intent(in) :: nxc
  real(p_double),               intent(in) :: t
  class( t_node_conf ), intent(in) :: no_co
  integer, intent(in) :: part_type

  integer :: i

  ! variables for collision-like sorting part
  integer :: i_dim, i_cell, i_s_temp
  integer :: num_merge_cells
  integer :: num_par_max_merge_cell
  integer :: aid0
  integer, dimension(:), pointer :: sap_merge_aid ! array for indexes
  integer :: npc_total ! total number of momentum cells

  type( t_merge ) :: merge


  do i_dim=1, p_x_dim
    merge%nx_merge_cells(i_dim) = merge_info%nx_merge_cells(i_dim)
  enddo

  npc_total = 1
  do i_dim=1, p_p_dim
    merge%np_per_dir(i_dim) = merge_info%np_per_dir(i_dim)
    npc_total=npc_total * merge%np_per_dir(i_dim)
  enddo


  if ( merge_info%n_check_max_sim_particles > 0 ) then
    if ( ( n > 0 ) .and. ( mod( n, merge_info%n_check_max_sim_particles ) == 0 ) ) then

      ! If larger than threshold then remove particles
      if ( (species%num_par > merge_info%max_sim_particles) .and. (species % num_par > 2) ) then
        
        !------------------------------------------------------------------------
        ! Beginning of space sorting routines
        !-----------------------------------------------------------------------
        ! calculate the number of merge cells and their size
        num_merge_cells = 1

        do i_dim=1, p_x_dim
          num_merge_cells = num_merge_cells * (nxc(i_dim) / merge%nx_merge_cells(i_dim))
        enddo

        ! allocate array for pointers to last particle in collision cell coll_aid
        ! polulate coll_aid by sorting species and callculate num_par_max_merge_cell

        call alloc( sap_merge_aid, (/num_merge_cells/) )

        call sort_merge_species( species, t, merge%nx_merge_cells, sap_merge_aid )

        ! maximum number of particles per collision cell
        num_par_max_merge_cell=0
        aid0=0
        do i_cell=1,num_merge_cells
          i_s_temp = sap_merge_aid(i_cell) - aid0
          aid0 = sap_merge_aid(i_cell)
          if( i_s_temp > num_par_max_merge_cell ) num_par_max_merge_cell = i_s_temp
        enddo ! merge cells

        !--------------------------------------------------------------
        ! end of space-sorting
        !-------------------------------------------------------------

        !------------------------------------------------------------------
        ! momentum space handling
        !----------------------------------------------------------------

        ! Allocate the merge info
        call merge % init( num_merge_cells, num_par_max_merge_cell, npc_total )

        ! Calculate the number of particles within the specific merge cell
        merge%mcell_i(1) = 1
        merge%mcell_numpar(1) = sap_merge_aid( 1 )

        do i = 2, num_merge_cells
          merge%mcell_numpar(i) = sap_merge_aid( i ) - sap_merge_aid( i-1 )
          merge%mcell_i(i) = sap_merge_aid( i-1 ) + 1
        enddo

        ! Merge particles
        ! Indexes of particles to be removed are stored in the delidx global variable
        call momentum_space( merge, species, num_merge_cells, npc_total, merge_info % max_weight, part_type )

        ! remove merged particles
        call remove_particles( species, delidx )

        ! free temporary memory
        call merge % cleanup()
        call freemem( sap_merge_aid )

      endif

    endif
  endif


end subroutine merge_particles
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
! Subroutine that handles the momentum space within the merge cell
!---------------------------------------------------------------------------------------------------
subroutine momentum_space(merge, species, num_merge_cells, npc_total, max_weight, part_type )

  use m_species_comm

  implicit none

  ! dummy variables
  type( t_merge ), intent(inout):: merge
  class( t_species ), intent(inout) :: species
  integer, intent(in) :: num_merge_cells 
  integer, intent(in) :: npc_total ! total number of momentum cells per merging cell
  real (p_k_part), intent(in) :: max_weight ! max weight of a particle participating in merging

  integer, intent(in) :: part_type


  ! local variables
  integer :: i, j, k
  integer, dimension(p_p_dim) :: ibox
  real ( p_k_part ), dimension(p_p_dim) :: pstep, rpstep, bin
  real ( p_k_part ) :: enetemp
  real( p_single ) :: randval(2)
  integer :: randint

  real ( p_k_part ) :: rqm, rm

  ! If merging photons the energy calculation is different
  select case ( part_type )
  case ( p_charge )
    rqm = species%rqm
    rm  = 1.0_p_k_part

  case ( p_photon )
    rqm = 1.0_p_k_part
    rm  = 0.0_p_k_part

  case default
    rqm = 0
    rm  = 0
    write(0,*) '(*error*) In os-mergedel.f03:momentum_space, invalid part_type.'
    call abort_program( p_err_invalid )
  end select

  ! allocate or adjust the buffer for the array of indexes to be deleted
  if ( delidx%buf_size < species%num_par ) then
    call freemem( delidx%idx )
    call alloc( delidx%idx, (/ species%num_par + 16384 /) )
    delidx%buf_size = size( delidx%idx )
  endif
  delidx%nidx = 0

  ! find the pmin and pmax of the particles in every merging cell
  do i=1,num_merge_cells

    merge%pmax(1,i) = species%p( 1 , merge%mcell_i(i) )
    merge%pmin(1,i) = species%p( 1 , merge%mcell_i(i) )
    merge%pmax(2,i) = species%p( 2 , merge%mcell_i(i) )
    merge%pmin(2,i) = species%p( 2 , merge%mcell_i(i) )
    merge%pmax(3,i) = species%p( 3 , merge%mcell_i(i) )
    merge%pmin(3,i) = species%p( 3 , merge%mcell_i(i) )
    do j=2,merge%mcell_numpar(i)

      if (species%p( 1 , j+merge%mcell_i(i)-1 ) > merge%pmax(1,i) ) then
        merge%pmax(1,i) = species%p( 1 , j+merge%mcell_i(i)-1 )
      endif
      if (species%p( 1 , j+merge%mcell_i(i)-1 ) < merge%pmin(1,i) ) then
        merge%pmin(1,i) = species%p( 1 , j+merge%mcell_i(i)-1 )
      endif

      if (species%p( 2 , j+merge%mcell_i(i)-1 ) > merge%pmax(2,i) ) then
        merge%pmax(2,i) = species%p( 2 , j+merge%mcell_i(i)-1 )
      endif
      if (species%p( 2 , j+merge%mcell_i(i)-1 ) < merge%pmin(2,i) ) then
        merge%pmin(2,i) = species%p( 2 , j+merge%mcell_i(i)-1 )
      endif

      if (species%p( 3 , j+merge%mcell_i(i)-1 ) > merge%pmax(3,i) ) then
        merge%pmax(3,i) = species%p( 3 , j+merge%mcell_i(i)-1 )
      endif
      if (species%p( 3 , j+merge%mcell_i(i)-1 ) < merge%pmin(3,i) ) then
        merge%pmin(3,i) = species%p( 3 , j+merge%mcell_i(i)-1 )
      endif
    enddo
  enddo



  do i=1,num_merge_cells

    ! binning the momemtum space
    do j = 1, p_p_dim
      pstep(j) = merge%pmax(j,i)-merge%pmin(j,i)
      if ( pstep(j) /= 0 ) then
        rpstep(j) = merge%np_per_dir(j) / pstep(j)
        bin(j) =  pstep(j) / merge%np_per_dir(j)
      else
        rpstep(j) = 0
        bin(j) = 0
      endif
    enddo

    do j = 1, npc_total
      merge%enetot(j) = 0
      merge%ptot(1,j) = 0
      merge%ptot(2,j) = 0
      merge%ptot(3,j) = 0
      merge%weighttot(j) = 0
      merge%numpar_per_pc(j) = 0
    enddo

    do j=1,merge%mcell_numpar(i)

      if (abs(species%q( j+merge%mcell_i(i)-1 )) < max_weight) then
        ! Calculate which momemtum cell the particle belongs to
        ibox(1) = int( (species%p( 1 , j+merge%mcell_i(i)-1 ) - merge%pmin(1,i) ) * rpstep(1) )+1
        if ( ibox(1) > merge%np_per_dir(1) ) ibox(1) = merge%np_per_dir(1)

        ibox(2)= int( (species%p( 2 , j+merge%mcell_i(i)-1 ) - merge%pmin(2,i) ) * rpstep(2) )+1
        if ( ibox(2) > merge%np_per_dir(2) ) ibox(2) = merge%np_per_dir(2)

        ibox(3)= int( (species%p( 3 , j+merge%mcell_i(i)-1 ) - merge%pmin(3,i) ) * rpstep(3) )+1
        if ( ibox(3) > merge%np_per_dir(3) ) ibox(3) = merge%np_per_dir(3)

      else
        ibox(1)=0
        ibox(2)=0
        ibox(3)=0
      endif
      ! the momentum cell index is stored in merge%part_pspace(j)
      ! for j-th particle in the i-th merging cell
      !  we check if species particle is supposed to merge
      if (abs(species%q( j+merge%mcell_i(i)-1 )) < max_weight) then
        merge%part_pspace(j) = ibox(1) + (ibox(2)-1)*merge%np_per_dir(1) + &
                              (ibox(3)-1) * merge%np_per_dir(1) * merge%np_per_dir(2)

        ! Add particle quantities (momemtum,energy, weight) to the corresponding momentum cell total
        merge%ptot( 1, merge%part_pspace(j)) = merge%ptot( 1, merge%part_pspace(j) ) + &
              rqm * species%q(j+merge%mcell_i(i)-1)*species%p( 1 , j+merge%mcell_i(i)-1 )
        merge%ptot( 2, merge%part_pspace(j)) = merge%ptot( 2, merge%part_pspace(j) ) + &
              rqm * species%q(j+merge%mcell_i(i)-1)*species%p( 2 , j+merge%mcell_i(i)-1 )
        merge%ptot( 3, merge%part_pspace(j)) = merge%ptot( 3, merge%part_pspace(j)) + &
              rqm * species%q(j+merge%mcell_i(i)-1)*species%p( 3 , j+merge%mcell_i(i)-1 )

        enetemp = sqrt( rm + species%p( 1 , j+merge%mcell_i(i)-1 )**2 + &
                             species%p( 2 , j+merge%mcell_i(i)-1 )**2 + &
                             species%p( 3 , j+merge%mcell_i(i)-1 )**2 )

        merge%enetot( merge%part_pspace(j) ) = merge%enetot( merge%part_pspace(j) ) + &
              rqm * species%q(j+merge%mcell_i(i)-1)*enetemp

        merge%weighttot( merge%part_pspace(j) ) = merge%weighttot( merge%part_pspace(j) ) + &
              rqm * species%q( j+merge%mcell_i(i)-1 )
        !increase the value on counter
        merge%numpar_per_pc( merge%part_pspace(j) ) = merge%numpar_per_pc( merge%part_pspace(j) ) + 1

      else
        merge%part_pspace(j)=0
      endif
    enddo

    ! loop to identify particles to be deleted within species merge cell and add them to delete list
    do k=1, npc_total

      if (merge%numpar_per_pc(k) > 2) then

        call part_mergeinto2( bin, merge%ptot(:,k), merge%enetot(k), merge%weighttot(k), &
          merge%pfirst(:,k), merge%psecond(:,k), merge%weightfirst(k), rm) 

        merge%weightfirst(k)=merge%weightfirst(k) * rqm
        merge%weightsecond(k)=merge%weightfirst(k)

        merge%enefirst(k) = sqrt(rm + merge%pfirst(1,k)**2 + merge%pfirst(2,k)**2 + &
          merge%pfirst(3,k)**2)

        merge%enesecond(k) = sqrt(rm + merge%psecond(1,k)**2 + merge%psecond(2,k)**2 + &
          merge%psecond(3,k)**2)

        ! the indexes of the two new particles correspond random former indexes of
        ! the array of the particles that were merged

        call random_number( randval )
        randint = int(randval(1)*merge%numpar_per_pc(k)) + 1
        merge%ifirst(k) = randint

        randint = int(randval(2)*merge%numpar_per_pc(k)) + 1

        if (randint == merge%ifirst(k)) then
          if (randint > 1) then
            merge%isecond(k) = randint - 1
          else
            merge%isecond(k) = 2
          endif
        else
          merge%isecond(k) = randint
        endif


      endif
    enddo

    ! looping over all particles of the merge cell
    do j=1,merge%mcell_numpar(i)
      if (merge%part_pspace(j)>0) then           !check if species is a valid particle for merging
        ! if the momentum cell of species particle has more than 2 particles
        if ( ( merge%numpar_per_pc( merge%part_pspace(j) ) > 2 ) .and.( merge%weightfirst(merge%part_pspace(j)) /= 0.0 )) then
          ! if species particle is one of the two particles to keep it gets new momentum value and weight
          if ( merge%ifirst( merge%part_pspace(j) ) == 1) then
            species%p( 1 , j+merge%mcell_i(i)-1 ) = merge%pfirst(1,merge%part_pspace(j))
            species%p( 2 , j+merge%mcell_i(i)-1 ) = merge%pfirst(2,merge%part_pspace(j))
            species%p( 3 , j+merge%mcell_i(i)-1 ) = merge%pfirst(3,merge%part_pspace(j))
            !                    species%ene( 1 , j+merge%mcell_i(i)-1 ) = merge%enefirst(merge%part_pspace(j))
            species%q( j+merge%mcell_i(i)-1 ) = merge%weightfirst(merge%part_pspace(j))
          elseif ( merge%isecond( merge%part_pspace(j) ) == 1) then
            species%p( 1 , j+merge%mcell_i(i)-1 ) = merge%psecond(1,merge%part_pspace(j))
            species%p( 2 , j+merge%mcell_i(i)-1 ) = merge%psecond(2,merge%part_pspace(j))
            species%p( 3 , j+merge%mcell_i(i)-1 ) = merge%psecond(3,merge%part_pspace(j))
            !                    species%ene( 1 , j+merge%mcell_i(i)-1 ) = merge%enesecond(merge%part_pspace(j))
            species%q( j+merge%mcell_i(i)-1 ) = merge%weightsecond(merge%part_pspace(j))
          else
            !if species particle is not one of the two, it is added to the deletion list
            delidx%nidx = delidx%nidx + 1   ! increasing the total number of particles to be deleted
            delidx%idx( delidx%nidx ) = j+merge%mcell_i(i)-1   !passing the absolute index of the particle to delete
          endif

          merge%ifirst( merge%part_pspace(j) ) = merge%ifirst( merge%part_pspace(j) ) - 1
          merge%isecond( merge%part_pspace(j) ) = merge%isecond( merge%part_pspace(j) ) - 1
        endif
      endif
    enddo
  enddo

end subroutine momentum_space

!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------
! Collision cell sorting routines
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine sort_merge_species( species, t, nx_merge_cells, coll_aid )
!---------------------------------------------------------------------------------------------------
! Sort the species along collision cells
!---------------------------------------------------------------------------------------------------

  implicit none

  class( t_species ), intent(inout) :: species

  real(p_double), intent(in) :: t
  integer, dimension(:), intent(in) :: nx_merge_cells
  integer, dimension(:), pointer          :: coll_aid

  ! local variables

  integer, dimension(:), pointer :: idx

  ! sort particles at certain timesteps if required

  ! need to sort for collisions even if numpart==0
  if (  (t > species%push_start_time)) then

    ! call begin_event(sort_merge_ev)
    call alloc( idx, (/species%num_par/) )

    ! call begin_event( sort_coll_genidx_ev )

    select case (p_x_dim)
    case (1)
      call generate_sort_merge_idx_1d(species, idx, nx_merge_cells, coll_aid )
    case (2)
      call generate_sort_merge_idx_2d(species, idx, nx_merge_cells, coll_aid )
    case (3)
      call generate_sort_merge_idx_3d(species, idx, nx_merge_cells, coll_aid )
    end select

    ! call end_event( sort_coll_genidx_ev )

    ! call begin_event( sort_coll_rearrange_ev )

    call rearrange( species, idx )

    !  call end_event( sort_coll_rearrange_ev )

    call freemem( idx )

    ! call end_event(sort_merge_ev)
  endif



end subroutine sort_merge_species
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine generate_sort_merge_idx_1d( species, ip, nx_collision_cells, coll_aid )
  !---------------------------------------------------------------------------------------------------
  !       generate sort indexes for a 1d run
  !---------------------------------------------------------------------------------------------------
  
    implicit none
  
    class( t_species ), intent(in) :: species
    integer, dimension(:), intent(out) :: ip
  
    integer, dimension(p_x_dim), intent(in) :: nx_collision_cells
    integer, dimension(:), pointer          :: coll_aid
  
  
    ! local variables
  
    integer, pointer, dimension(:) :: npic
    integer :: i
    integer :: index
  
    integer :: index_xpic
    integer :: index_xcoll
    integer :: index_dxpic
    integer :: n_coll_cells_x, n_coll_cells
    integer :: s_coll_cell_x, s_coll_cell
  
    integer :: isum, ist
    integer :: n_grid, n_grid_x
  
    n_grid_x = species%my_nx_p(3, 1)
    n_grid = n_grid_x
  
    s_coll_cell_x = nx_collision_cells(1)
    s_coll_cell = s_coll_cell_x
    n_coll_cells_x = n_grid_x / s_coll_cell_x
    n_coll_cells = n_coll_cells_x
  
  
    call alloc(npic, (/ n_grid /) )
  
    npic = 0
  
    ! This sorts by collision cell first and interpolation cell second
    do i=1,species%num_par
      index_xpic = species%ix(1,i)-1
      index_xcoll = index_xpic / s_coll_cell_x
      index_dxpic = index_xpic - index_xcoll * s_coll_cell_x
      index = s_coll_cell * index_xcoll + index_dxpic + 1
  
      npic(index) = npic(index) + 1
      ip(i)=index
    end do
  
    isum=0
    do i=1,n_grid
      ist = npic(i)
      npic(i) = isum
      isum = isum + ist
    end do
  
    ! isum must be the same as the total number of particles
    ! ASSERT(isum == species%num_par)
  
    do i=1,species%num_par
      index=ip(i)
      npic(index) = npic(index) + 1
      ip(i) = npic(index)
    end do
  
    !return array with the last indexes per collision cell
    do i=1,n_coll_cells
      coll_aid(i) = npic(i*s_coll_cell)
    enddo
  
    call freemem(npic)
  
  end subroutine generate_sort_merge_idx_1d
  !---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine generate_sort_merge_idx_2d( species, ip, nx_collision_cells, coll_aid )
!---------------------------------------------------------------------------------------------------
!       generate sort indexes for a 2d run
!---------------------------------------------------------------------------------------------------

  implicit none

  class( t_species ), intent(in) :: species
  integer, dimension(:), intent(out) :: ip

  integer, dimension(p_x_dim), intent(in) :: nx_collision_cells
  integer, dimension(:), pointer          :: coll_aid


  ! local variables

  integer, pointer, dimension(:) :: npic
  integer :: i
  integer :: index

  integer :: index_xpic, index_ypic
  integer :: index_xcoll, index_ycoll
  integer :: index_dxpic, index_dypic
  integer :: n_coll_cells_x, n_coll_cells_y, n_coll_cells
  integer :: s_coll_cell_x, s_coll_cell_y, s_coll_cell

  integer :: isum,ist
  integer :: n_grid, n_grid_x, n_grid_y

  n_grid_x = species%my_nx_p(3, 1)
  n_grid_y = species%my_nx_p(3, 2)
  n_grid = n_grid_x * n_grid_y

  s_coll_cell_x = nx_collision_cells(1)
  s_coll_cell_y = nx_collision_cells(2)
  s_coll_cell = s_coll_cell_x * s_coll_cell_y
  n_coll_cells_x = n_grid_x / s_coll_cell_x
  n_coll_cells_y = n_grid_y / s_coll_cell_y
  n_coll_cells = n_coll_cells_x * n_coll_cells_y


  call alloc(npic, (/ n_grid /) )

  npic = 0

  ! This sorts by collision cell first and interpolation cell second
  do i=1,species%num_par
    index_xpic = species%ix(1,i)-1
    index_ypic = species%ix(2,i)-1
    index_xcoll = index_xpic / s_coll_cell_x
    index_ycoll = index_ypic / s_coll_cell_y
    index_dxpic = index_xpic - index_xcoll * s_coll_cell_x
    index_dypic = index_ypic - index_ycoll * s_coll_cell_y
    !index = species%ix(i1,i) + n_grid_x * (species%ix(i2,i)-1)       ! sort by y then x
    index = s_coll_cell * ( index_xcoll + n_coll_cells_x * index_ycoll ) + &
    index_dxpic + s_coll_cell_x * index_dypic + 1

    npic(index) = npic(index) + 1
    ip(i)=index
  end do

  isum=0
  do i=1,n_grid
    ist = npic(i)
    npic(i) = isum
    isum = isum + ist
  end do

  ! isum must be the same as the total number of particles
  ! ASSERT(isum == species%num_par)

  do i=1,species%num_par
    index=ip(i)
    npic(index) = npic(index) + 1
    ip(i) = npic(index)
  end do

  !return array with the last indexes per collision cell
  do i=1,n_coll_cells
    coll_aid(i) = npic(i*s_coll_cell)
  enddo

  call freemem(npic)

end subroutine generate_sort_merge_idx_2d
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine generate_sort_merge_idx_3d( species, ip, nx_collision_cells, coll_aid )
!---------------------------------------------------------------------------------------------------
!       generate sort indexes for a 3d run
!---------------------------------------------------------------------------------------------------

  implicit none

  ! dummy variables

  class( t_species ), intent(in) :: species
  integer, dimension(:), intent(out) :: ip

  integer, dimension(p_x_dim), intent(in) :: nx_collision_cells
  integer, dimension(:), pointer          :: coll_aid

  ! local variables

  integer, pointer, dimension(:) :: npic
  integer :: i
  integer :: index
  integer :: isum,ist

  integer :: n_grid_x,    n_grid_y,    n_grid_z,   n_grid
  integer :: index_xpic,  index_ypic,  index_zpic
  integer :: index_xcoll, index_ycoll, index_zcoll
  integer :: index_dxpic, index_dypic, index_dzpic
  integer :: n_coll_cells_x, n_coll_cells_y, n_coll_cells_z, n_coll_cells
  integer :: s_coll_cell_x, s_coll_cell_y, s_coll_cell_z, s_coll_cell

  n_grid_x = species%my_nx_p(3, 1)
  n_grid_y = species%my_nx_p(3, 2)
  n_grid_z = species%my_nx_p(3, 3)
  n_grid = n_grid_x * n_grid_y * n_grid_z

  s_coll_cell_x = nx_collision_cells(1)
  s_coll_cell_y = nx_collision_cells(2)
  s_coll_cell_z = nx_collision_cells(3)

  s_coll_cell    =  s_coll_cell_x  * s_coll_cell_y  * s_coll_cell_z

  n_coll_cells_x = n_grid_x / s_coll_cell_x
  n_coll_cells_y = n_grid_y / s_coll_cell_y
  n_coll_cells_z = n_grid_z / s_coll_cell_z
  n_coll_cells = n_coll_cells_x * n_coll_cells_y * n_coll_cells_z


  call alloc(npic, (/ n_grid /) )

  npic = 0

  ! This sorts by collision cell first and interpolation cell second
  do i=1,species%num_par
    index_xpic = species%ix(1,i)-1
    index_ypic = species%ix(2,i)-1
    index_zpic = species%ix(3,i)-1

    index_xcoll = index_xpic / s_coll_cell_x
    index_ycoll = index_ypic / s_coll_cell_y
    index_zcoll = index_zpic / s_coll_cell_z

    index_dxpic = index_xpic - index_xcoll * s_coll_cell_x
    index_dypic = index_ypic - index_ycoll * s_coll_cell_y
    index_dzpic = index_zpic - index_zcoll * s_coll_cell_z

    index = s_coll_cell * ( index_xcoll + n_coll_cells_x * &
      ( index_ycoll + n_coll_cells_y * index_zcoll ) ) + &
    index_dxpic + s_coll_cell_x * ( index_dypic + s_coll_cell_y  * index_dzpic ) + 1

    npic(index) = npic(index) + 1
    ip(i)=index
  end do

  isum=0
  do i=1,n_grid
    ist = npic(i)
    npic(i) = isum
    isum = isum + ist
  end do

  do i=1,species%num_par
    index=ip(i)
    npic(index) = npic(index) + 1
    ip(i) = npic(index)
  end do

  !return array with the last indexes per collision cell
  do i=1,n_coll_cells
    coll_aid(i) = npic(i*s_coll_cell)
  enddo

  call freemem(npic)

end subroutine generate_sort_merge_idx_3d
!---------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Produce 2 macro-particles from the sum of the momentum, energy and weight of N particles.
! The procedure conserves exactly by default the sumed quantities
! The plane in which the two new vectors lie is random
!-----------------------------------------------------------------------------------------
subroutine part_mergeinto2( bin, psin, enesin, weights, p1, p2, weightnew,rm)

  implicit none

  real( p_k_part ), dimension(:), intent(in) :: bin
  real( p_k_part ), dimension(:), intent(in) :: psin
  real( p_k_part ), intent(in) :: enesin, weights
  real( p_k_part ), dimension(:), intent(inout) :: p1, p2
  real( p_k_part ), intent(inout) :: weightnew
  real( p_k_part ),intent(in) :: rm


  real( p_k_part ) :: psnorm, psnorm2, randval0, randval, randval2, randval3, r1, r2, absrand
  real( p_k_part ) :: costheta1, sintheta1, costheta2, sintheta2
  real( p_k_part ) :: vpcxperp, vpcyperp, vpczperp, vpcx, vpcy, vpcz, vpcperpnorm, projection
  real( p_k_part ) :: ex, ey, ez
  real( p_k_part ) :: rweights, pxs, pys, pzs, enes

  integer :: randint

  ! new weight of the macro photons
  weightnew = 0.5_p_k_part*weights
  rweights = 1.0_p_k_part / weights
  ! calculate the average momenta and energy within the cell (dividing by the total weight)
  pxs = psin(1) * rweights
  pys = psin(2) * rweights
  pzs = psin(3) * rweights
  enes = enesin * rweights

  ! norm of the average momentum p
  psnorm = sqrt(pxs**2 + pys**2 + pzs**2)
  psnorm2 = pxs**2 + pys**2 + pzs**2

  if ( ( enes**2 >= (psnorm**2 + rm) ) .and. ( psnorm > 0 ) ) then ! problem: not for photons FGD

    !the new momemtum norm is the same for both new particles
    ! the energy of the new particles is equal to the average energy of the momemtum cell

    r1 = sqrt(enes**2-rm)
    r2 = r1

    costheta1 = psnorm/r1
    sintheta1 = sqrt(1 - costheta1**2)
    costheta2 = psnorm/r2
    sintheta2 = sqrt(1 - costheta2**2)

    ! selection of the plane perpendicular to the momemtum vector parallel with the chosen pcell vector

    call random_number( randval0 )
    randint = int(randval0*4.0d0) + 1

    select case ( randint )
      case( 1 )
        vpcx = bin(1)
        vpcy = bin(2)
        vpcz = bin(3)
      case( 2 )
        vpcx = bin(1)
        vpcy = bin(2)
        vpcz = -bin(3)
      case( 3 )
        vpcx = bin(1)
        vpcy = -bin(2)
        vpcz = bin(3)
      case( 4 )
        vpcx = -bin(1)
        vpcy = bin(2)
        vpcz = bin(3)
      case default
        ! This is only used to avoid compiler warnings stating that vpc? may be used uninitialized
        vpcx = 0
        vpcy = 0
        vpcz = 0
    end select

    projection = ( vpcx * pxs + vpcy * pys + vpcz * pzs ) / psnorm

    vpcxperp = vpcx - projection * pxs / psnorm
    vpcyperp = vpcy - projection * pys / psnorm
    vpczperp = vpcz - projection * pzs / psnorm

    vpcperpnorm = sqrt( vpcxperp**2 + vpcyperp**2 + vpczperp**2 )

    if ( vpcperpnorm > 1.0d-12) then

      !vector perpendicular to the total p
      ex = vpcxperp / vpcperpnorm
      ey = vpcyperp / vpcperpnorm
      ez = vpczperp / vpcperpnorm

      ! calculating components of the momenta of the two particles
      p1(1) =  r1 * sintheta1 * ex + pxs
      p1(2) =  r1 * sintheta1 * ey + pys
      p1(3) =  r1 * sintheta1 * ez + pzs

      p2(1) =  -r2 * sintheta1 * ex + pxs
      p2(2) =  -r2 * sintheta1 * ey + pys
      p2(3) =  -r2 * sintheta1 * ez + pzs

    else
      weightnew = 0.0d0
      p1(1) = 0.0d0
      p1(2) = 0.0d0
      p1(3) = 0.0d0

      p2(1) = 0.0d0
      p2(2) = 0.0d0
      p2(3) = 0.0d0

    endif

  elseif ( ( enes > psnorm ) .and. ( psnorm == 0.0d0 ) ) then
    call random_number( randval)
    call random_number( randval2)
    call random_number( randval3)
    absrand = sqrt(randval**2 + randval2**2 + randval3**2)

    p1(1) = enes*randval / absrand
    p1(2) = enes*randval2 / absrand
    p1(3) = enes*randval3 / absrand

    p2(1) = -p1(1)
    p2(2) = -p1(2)
    p2(3) = -p1(3)

  else
    weightnew = 0.0d0

    p1(1) =  0.0d0
    p1(2) =  0.0d0
    p1(3) =  0.0d0

    p2(1) =  0.0d0
    p2(2) =  0.0d0
    p2(3) =  0.0d0

  endif


end subroutine part_mergeinto2

end module m_mergedel
