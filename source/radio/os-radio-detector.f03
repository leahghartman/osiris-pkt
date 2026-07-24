!----------------------------------------------------------------------------------------!
!--------------------------------- DETECTOR BASE CLASS ----------------------------------!
!                                                                                        !
! Here we define a detector objecte with the required methods and parameters.            !
!                                                                                        !
! This is an abstract class, we then use polymorphism to define the derived classes      !
! (spherical detector, cartesian detector, etc...)                                       !
!                                                                                        !
! Every function that you may use in the derived class must be defined here as a         !
! deferred procedure with a matching definition in the interface block                   !
!                                                                                        !
!----------------------------------------------------------------------------------------!

module Detector_cls

#include "memory/memory.h"
use m_system
use m_parameters
use m_input_file, only : t_input_file
use m_time_step, only : n, ndump, t_time_step, test_if_report
use m_diagnostic_utilities, only : t_diag_file, p_diag_create, p_diag_grid, &
                                   create_diag_file
use m_diagfile, only : p_time_length, diag_axis_linear
use stringutil, only : idx_string
use m_node_conf, only : t_node_conf, comm, root, n_threads
use m_restart, only : t_restart_handle

implicit none

private

type, abstract, public :: Detector

  !number of spacetime cell
  integer :: n1, n2, n3, n4, nf(4)
  integer :: ncmp
  !Detector caracteriation
  real(p_k_part) :: dmin(4), d2, d3, d4, dmax(4)
  real(p_k_part), pointer :: dmin_(:) => null(), dmax_(:) => null()
  real(p_double) :: d1, sim_dt, d1m1
  integer, pointer :: dims(:) => null()
  integer(p_int64), pointer :: iter_bnds(:,:) => null()
  integer :: ndims
  integer :: place
  !stuff for diagnostics
  integer :: ndump_fac
  logical :: cmpts(6)
  real(p_single), pointer :: Emf(:,:,:,:) => null()
  real(p_double), pointer :: tpos(:) => null()
  integer :: dim_pos(2)
  !variables for output stuff
  character(len=21), dimension(4) :: bas_ax_nms, bas_ax_lnms, bas_ax_uns
  !name of members
  character(len=21), dimension(6) :: bas_fld_nms, bas_cmp_nam, bas_cmp_lnam
  character(len=21), dimension(1) :: bas_fld_uns
  character(len=256) :: basePath = ''
  class( Detector ), pointer :: next => null()

contains

  procedure(destr),deferred :: destruct
  procedure(filp),deferred :: FillPos
  procedure(read_inp),deferred :: read_input
  procedure(initt),deferred :: init
  procedure(rstr_wr), deferred :: restart_write
  procedure(rstr_rd), deferred :: restart_read
  procedure(fill_emf_det), deferred :: fill_emf
  procedure :: report => report_det

end type Detector

abstract interface
  subroutine destr(this)
    import Detector
    class(Detector), intent (INOUT) :: this
  end subroutine
  subroutine filp(this)
    import Detector
    class(Detector), INTENT(INOUT):: this
  end subroutine
  subroutine read_inp(this,input_file,dt)
    import Detector, t_input_file, p_double
    class(Detector), INTENT(INOUT):: this
    class( t_input_file ), intent(inout) :: input_file
    real(p_double),intent(in)::dt
  end subroutine
  subroutine initt(this,name,no_co,restart,restart_handle)
    import Detector, t_node_conf, t_restart_handle
    class( Detector ), intent(inout) :: this
    character(len = 64) ,intent(in):: name
    class( t_node_conf ), intent(in) :: no_co
    logical, intent(in) :: restart
    type(t_restart_handle), intent(in) :: restart_handle
  end subroutine
  subroutine rstr_wr( this, restart_handle )
    import Detector, t_restart_handle
    class(Detector), intent(in) :: this
    type(t_restart_handle), intent(inout) :: restart_handle
  end subroutine
  subroutine rstr_rd( this, restart_handle )
    import Detector, t_restart_handle
    class(Detector), intent(inout) :: this
    type(t_restart_handle), intent(in) :: restart_handle
  end subroutine
  subroutine fill_emf_det(this,x,beta,time,xprev,betaprev,charge,np,tid)
    import Detector, p_k_part, p_cache_size, p_p_dim, p_double
    class(Detector), intent(inout):: this
    real(p_k_part), dimension(3,p_cache_size), intent(in) :: x, xprev
    real(p_k_part), dimension(p_p_dim,p_cache_size), intent(in) :: beta, betaprev
    real(p_k_part), dimension(p_cache_size), intent(in) :: charge
    real(p_double), intent(in) :: time
  integer, intent(in) :: np,tid
  end subroutine
end interface

contains

!-----------------------------------------------------------------------------------------
subroutine outp_det_file( det, t, timestep, dir, write_buffer )

  implicit none

  class(Detector), intent(in) :: det
  integer, intent(in) :: dir
  real(p_double), intent(in) :: t
  type( t_time_step ), intent(in) :: timestep
  real(p_single), dimension(:,:,:), pointer, intent(in) :: write_buffer

  integer :: i, ierr
  class( t_diag_file ), allocatable :: diagFile

  ! Sanity check
  if ( det % ndims > 3 ) then
    if ( mpi_node() == 0 ) then
      write(0,"(A)") "(*error*) Unable to write output for radio detector"
      write(0,"(A,I0,A)") "(*error*) Invalid dimensions (", det % ndims, " > 3 )"
    endif
    call abort_program( p_err_invalid )
  endif

  ! Create diagFile object of type grid
  call create_diag_file( diagFile )
  diagFile%ftype = p_diag_grid

  i=dir
  diagFile % filename = trim(det%bas_cmp_nam(dir)) // '-' // &
                        idx_string( n(timestep)/ndump(timestep), p_time_length )
  diagFile % filepath = trim(det%basePath)//trim(det%bas_cmp_nam(i))//p_dir_sep
  diagFile % name     = 'Radiated '//det%bas_cmp_lnam(dir)

  ! Iteration info
  diagFile%iter%n          = n(timestep)
  diagFile%iter%t          = t
  diagFile%iter%time_units = "1 / \omega_p"

  ! Initialize grid dataset parameters
  diagFile%grid%ndims = det%ndims
  diagFile%grid%name  = det%bas_cmp_nam(dir)
  diagFile%grid%label = det%bas_cmp_lnam(dir)
  diagFile%grid%units = det%bas_fld_uns(1)

  do i = 1, det%ndims
    diagFile%grid%count(i) = det%dims(i)
    diagFile%grid%axis(i)%type = diag_axis_linear
    diagFile%grid%axis(i)%min = det%dmin_(i)
    diagFile%grid%axis(i)%max = det%dmax_(i)
    diagFile%grid%axis(i)%name  = det%bas_ax_nms(i)
    diagFile%grid%axis(i)%label = det%bas_ax_lnms(i)
    diagFile%grid%axis(i)%units = det%bas_ax_uns(i)
  enddo

  ! Create directory
  call mkdir( diagFile%filepath, ierr )

  ! Open the file
  call diagFile % open( p_diag_create )
  select case( det%ndims )
  case(1)
    call diagFile % add_dataset( det%bas_cmp_nam(dir), write_buffer(1,1,:))
  case(2)
    call diagFile % add_dataset( det%bas_cmp_nam(dir), write_buffer(:,1,:))
  case(3)
    call diagFile % add_dataset( det%bas_cmp_nam(dir), write_buffer(:,:,:))
  end select

  ! Close the file
  call diagFile % close( )

  deallocate( diagFile )

end subroutine outp_det_file
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine report_det(this,t,tstep,no_co)

  implicit none

  class( Detector ), intent(inout) :: this
  real(p_double), intent(in) :: t
  type( t_time_step ), intent(in) :: tstep
  class( t_node_conf ), intent(in) :: no_co

  integer :: i, j, k, l, dp1, dp2
  integer :: dir, count, ierr, mpi_type
  real(p_single), dimension(:,:,:), pointer :: write_buffer

  write_buffer => null()

  if ( test_if_report(tstep,this%ndump_fac) ) then

    dp1 = this%dim_pos(1)
    dp2 = this%dim_pos(2)

    call alloc( write_buffer, (/ this%nf(dp1), this%nf(dp2), this%nf(1)+1 /) )

    ! Get size of message buffer
    count = size( write_buffer )
    i=0
    do dir=1,this%ncmp
      i = i+1
      do while (.not. this%cmpts(i))
        i = i+1
      end do

      ! Fill write_buffer array
      do l=1,this%nf(1)+1
        do k=1,this%nf(dp2)
          do j=1,this%nf(dp1)
            write_buffer(j,k,l) = this%Emf(dir,j,k,l)
          enddo
        enddo
      enddo

      ! If running in parallel reduce data to node 0
      if ( no_co % no_num > 1 ) then
        mpi_type = mpi_real_type( p_single )
        if ( root(no_co) ) then
          call MPI_REDUCE( MPI_IN_PLACE, write_buffer, count, mpi_type, MPI_SUM, 0, &
                           comm(no_co), ierr)
        else
          call MPI_REDUCE( write_buffer, 0, count, mpi_type, MPI_SUM, 0, &
                           comm(no_co), ierr )
        endif
      endif

      if ( root(no_co) ) then
        do l=1,this%nf(1)+1
          do k=1,this%nf(dp2)
            do j=1,this%nf(dp1)
              write_buffer(j,k,l) = real( write_buffer(j,k,l) * this%d1m1, p_single )
            enddo
          enddo
        enddo
        call outp_det_file(this,t,tstep,i,write_buffer)
      endif

    enddo

    call freemem( write_buffer )

  endif

end subroutine report_det
!-----------------------------------------------------------------------------------------

end module Detector_cls
