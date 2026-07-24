! if __TEMPLATE__ is defined then read the template definition at the end of the file
#ifndef __TEMPLATE__

module class_CartDet

#include "os-config.h"
#include "os-preprocess.fpp"
#include "memory/memory.h"

use m_system
use m_parameters
use detector_cls
use m_input_file, only : t_input_file, get_namelist
use m_restart
use stringutil, only : replace_blanks
use m_node_conf, only : t_node_conf, n_threads

implicit none

private

! string to id restart data
character(len=*), parameter :: p_cart_det_rst_id = "cart_det rst data - 0x0000"

type, public, extends(Detector) :: cart_Detector
  real(p_k_part), pointer :: xpos(:) => null(), ypos(:) => null(), zpos(:) => null()
  ! Varies with field component
  procedure(fill_emf_loc_cart_det), pointer :: fill_emf_loc => null()
contains
  procedure :: destruct => destruct_cart_det
  procedure :: FillPos => FillPos_cart_det
  procedure :: read_input => read_input_cart_det
  procedure :: init => init_cart_det
  procedure :: restart_write => restart_write_cart_det
  procedure :: restart_read => restart_read_cart_det
  procedure :: fill_emf => fill_emf_cart_det
  procedure :: comp_select => comp_select_cart_det
end type cart_Detector

abstract interface
  subroutine fill_emf_loc_cart_det(this,x,beta,time,xprev,betaprev,charge,np,tid)
    import cart_Detector, p_k_part, p_cache_size, p_p_dim, p_double
    class(cart_Detector), intent(inout):: this
    real(p_k_part), dimension(3,p_cache_size), intent(in) :: x, xprev
    real(p_k_part), dimension(p_p_dim,p_cache_size), intent(in) :: beta, betaprev
    real(p_k_part), dimension(p_cache_size), intent(in) :: charge
    real(p_double), intent(in) :: time
    integer, intent(in) :: np, tid
  end subroutine
end interface

contains

!-----------------------------------------------------------------------------------------
subroutine read_input_cart_det(this,input_file,dt)

  implicit none

  real(p_double), intent(in) :: dt
  class( cart_Detector ), intent(inout) :: this
  class( t_input_file ), intent(inout) :: input_file
  integer :: ndump_fac_det
  integer :: nf(4), ierr
  real(p_k_part) :: dmin(4), dmax(4)
  logical :: comps(6)
#ifndef __HAS_RAD_BFLD__
  logical :: has_b
  integer :: i
#endif
  namelist /nl_cartesian/ nf, dmin, dmax, comps, ndump_fac_det

  ndump_fac_det = 0
  nf = 0
  dmin = 0
  dmax = 0
  comps = .false.

  ! Get namelist text from input file
  call get_namelist( input_file, "nl_cartesian", ierr )
  if (ierr /= 0) then
    if(mpi_node()==0) then
      if (ierr < 0) then
        write(0,*) "Error reading detector parameters"
      else
        ! This should never happen because of the detector_list
        write(0,*) "Error: detector parameters missing"
      endif
      write(0,*) "aborting..."
    endif
    stop
  endif

  read (input_file%nml_text, nml = nl_cartesian, iostat = ierr)
  if (ierr /= 0) then
    if(mpi_node()==0) then
      write(0,*) "Error reading detector parameters"
      write(0,*) "aborting..."
    endif
    stop
  endif

  this%nf=nf
  this%dmin=dmin
  this%dmax=dmax
  this%cmpts=comps
  this%ndump_fac=ndump_fac_det
  this%sim_dt=dt

  ! Report any issues
#ifndef __HAS_RAD_BFLD__
  ! If any B component is required
  has_b = .false.
  do i = 4, 6
    if (this%cmpts(i)) then
      has_b = .true.
      this%cmpts(i)=.False.
    endif
  enddo
  if(has_b.and.mpi_node()==0 )then
    write(0,*) "Warning: __HAS_RAD_BFLD__ not defined"
    write(0,*) "RaDiO was not compiled with B-Field support"
    write(0,*) "Selected B-field components will be ignored"
  endif
#endif

  this%ncmp=count(this%cmpts)

  this%ndims=count(this%nf(2:4)>1)
  if (this%ndims>2) then
    if(mpi_node()==0) then
      write(0,*) "Error: detector grid has more than 2 spatial dimensions"
      write(0,*) "aborting..."
    endif
    stop
  endif

  if (this%nf(1)<2) then
    if(mpi_node()==0) then
      write(0,*) "Error: detector time grid has more less than 2 cells"
      write(0,*) "aborting..."
    endif
    stop
  endif
  this%ndims=this%ndims+1

end subroutine read_input_cart_det
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine init_cart_det(this,name,no_co,restart,restart_handle)

  class( cart_Detector ), intent(inout) :: this
  character(len = 64), intent(in):: name
  class( t_node_conf ), intent(in) :: no_co
  logical, intent(in) :: restart
  type( t_restart_handle ), intent(in) :: restart_handle

  integer :: i, j, dp1, dp2, nt
  integer(p_int64) :: ncells, chunk, i0, i1
  integer :: sizes(4)
  real(p_k_part) :: lbounds(4), ubounds(4)
  character(len=256) :: fileLabel = '',aux2
  character(len=20), dimension(4) :: bas_ax_nms, bas_ax_lnms, bas_ax_uns ! name of members

  call alloc(this%dims,(/ this%ndims /))
  call alloc(this%dmin_,(/ this%ndims /))
  call alloc(this%dmax_,(/ this%ndims /))

  this%bas_cmp_nam=(/"E1    ","E2    ","E3    ","B1    ","B2    ","B3    "/)
  this%bas_cmp_lnam=(/"E_1    ","E_2    ","E_3    ","B_1    ","B_2    ","B_3    "/)
  call this % comp_select()

  this%bas_fld_uns = (/ "m_e c \omega_p e^{-1}" /)
  write(aux2,"(A4,I2.2)") "_det",this%place
  fileLabel = replace_blanks(trim(name))//aux2
  this%basePath = trim(path_mass) // trim("RADIAT")  // p_dir_sep // &
                  trim(fileLabel) // p_dir_sep

  sizes(4)=this%nf(1)+1
  sizes(:3)=this%nf(2:4)
  lbounds(4)=this%dmin(1) - (this%dmax(1)-this%dmin(1))/(this%nf(1)-1)
  lbounds(:3)=this%dmin(2:4)
  ubounds(4)=this%dmax(1)
  ubounds(:3)=this%dmax(2:4)

  !set axis names
  bas_ax_nms  =(/"x1","x2","x3","t "/)
  bas_ax_lnms =(/"x_2","x_2","x_3","t  "/)
  bas_ax_uns  =(/"c/\omega_p","c/\omega_p","c/\omega_p","1/\omega_p" /)

  j=1
  this%dim_pos=(/2,3/)
  do i=1,this%ndims
    do while(sizes(j)<=1)
      j=j+1
    end do
    this%bas_ax_nms(i)=bas_ax_nms(j)
    this%bas_ax_lnms(i)=bas_ax_lnms(j)
    this%bas_ax_uns(i)=bas_ax_uns(j)

    this%dims(i)=sizes(j)
    this%dmin_(i)=lbounds(j)
    this%dmax_(i)=ubounds(j)
    if(i<this%ndims) then
      this%dim_pos(2)=j+2
      this%dim_pos(i)=j+1
    endif
    j=j+1
  end do
  ! Correct for case where we only have spatial cells in the last dimension
  if (this%dim_pos(2)==5) this%dim_pos(2) = 3

  !1,2,3,4=t,x,y,z
  this%n2  =this%nf(2) !x
  this%n3  =this%nf(3) !y
  this%n4  =this%nf(4) !z
  this%n1  =this%nf(1) !t

  this%d1  =(this%dmax(1)-this%dmin(1))
  this%d2  =(this%dmax(2)-this%dmin(2))
  this%d3  =(this%dmax(3)-this%dmin(3))
  this%d4  =(this%dmax(4)-this%dmin(4))

  if(this%nf(1)>1) this%d1 = (this%dmax(1)-this%dmin(1))/(this%nf(1)-1)
  if(this%nf(2)>1) this%d2 = (this%dmax(2)-this%dmin(2))/(this%nf(2)-1)
  if(this%nf(3)>1) this%d3 = (this%dmax(3)-this%dmin(3))/(this%nf(3)-1)
  if(this%nf(4)>1) this%d4 = (this%dmax(4)-this%dmin(4))/(this%nf(4)-1)
  this%d1m1=1/this%d1

  call alloc(this%xpos,(/this%nf(2)/))
  call alloc(this%ypos,(/this%nf(3)/))
  call alloc(this%zpos,(/this%nf(4)/))
  call alloc(this%tpos,(/this%nf(1)+1/))

  if (this%ncmp>=1) then
    dp1 = this%dim_pos(1)
    dp2 = this%dim_pos(2)

    call alloc( this%Emf, (/ this%ncmp, this%nf(dp1), this%nf(dp2), this%nf(1)+1 /) )
    this%Emf = 0.0_p_single

    if (restart) then
      call this % restart_read( restart_handle )
    endif

  endif

  call this%FillPos()

  ncells = this%nf(this%dim_pos(1))*this%nf(this%dim_pos(2))
  nt = n_threads(no_co)

  if (ncells < nt) then
    if(mpi_node()==0) then
      write(0,*) "Error: in the current detector implementation, when using multiple"
      write(0,*) "threads, the total number of spatial detector cells must be greater"
      write(0,*) "than or equal to the number of threads."
      write(0,*) "aborting..."
    endif
    stop
  endif

  call alloc( this%iter_bnds, (/ 2, nt /) )
  do i = 0, nt-1
    chunk = ( ncells + nt - 1 ) / nt
    i0    = i * chunk + 1
    i1    = min( (i+1) * chunk, ncells )
    this%iter_bnds(p_lower,i+1) = i0
    this%iter_bnds(p_upper,i+1) = i1
  enddo

end subroutine init_cart_det
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine FillPos_cart_det(this)

  class(cart_Detector), INTENT(INOUT):: this
  integer :: k
  !efficient array
  !do i
  !  do j
  !    x(j,i)

  !$omp parallel
  !$omp do private(k)
  do k=0,this%n2-1
    this%xpos(k+1)=this%dmin(2)+k*this%d2
  enddo
  !$omp enddo
  !$omp do private(k)
  do k=0,this%n3-1
    this%ypos(k+1)=this%dmin(3)+k*this%d3
  enddo
  !$omp enddo

  !$omp do private(k)
  do k=0,this%n4-1
    this%zpos(k+1)=this%dmin(4)+k*this%d4
  enddo
  !$omp enddo

  !$omp do private(k)
  do k=-1,this%n1-1
    this%tpos(k+2)=this%dmin(1)+k*this%d1
  enddo
  !$omp enddo
  !$omp end parallel
end subroutine FillPos_cart_det
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine destruct_cart_det(this)
  class(cart_Detector), INTENT(INOUT):: this

  call freemem(this%xpos)
  call freemem(this%ypos)
  call freemem(this%zpos)
  call freemem(this%tpos)
  call freemem(this%Emf)
  call freemem(this%dims)
  call freemem(this%dmin_)
  call freemem(this%dmax_)
  call freemem(this%iter_bnds)

end subroutine destruct_cart_det
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine restart_write_cart_det( this, restart_handle )

  implicit none

  class(cart_Detector), intent(in) :: this
  type( t_restart_handle ), intent(inout) :: restart_handle

  character(len=*), parameter :: err_msg='error writing restart data for cart_det object.'
  integer :: ierr

  ! Write data
  restart_io_wr( p_cart_det_rst_id, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

  restart_io_wr( this%Emf, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstwrt )

end subroutine restart_write_cart_det
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine restart_read_cart_det( this, restart_handle )

  implicit none

  class(cart_Detector), intent(inout) :: this
  type( t_restart_handle ), intent(in) :: restart_handle

  character(len=*), parameter :: err_msg='error reading restart data for cart_det object.'
  character(len=len(p_cart_det_rst_id)) :: rst_id
  integer :: ierr

  ! Read data
  restart_io_rd( rst_id, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd )

  ! check if restart file is compatible
  if ( rst_id /= p_cart_det_rst_id ) then
    ERROR('Corrupted restart file, or restart file ')
    ERROR('from incompatible binary (cart_det)')
    call abort_program(p_err_rstrd)
  endif

  restart_io_rd( this%Emf, restart_handle, ierr )
  CHECK_ERROR( ierr, err_msg, p_err_rstrd )

end subroutine restart_read_cart_det
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine fill_emf_cart_det(this,x,beta,time,xprev,betaprev,charge,np,tid)
  implicit none
  class(cart_Detector), intent(inout) :: this
  real(p_k_part), dimension(3,p_cache_size), intent(in) :: x, xprev
  real(p_k_part), dimension(p_p_dim,p_cache_size), intent(in) :: beta, betaprev
  real(p_k_part), dimension(p_cache_size), intent(in) :: charge
  real(p_double), intent(in) :: time
  integer, intent(in) :: np, tid

  ! To allow for using approx as true or false
  call this%fill_emf_loc(x,beta,time,xprev,betaprev,charge,np,tid)

end subroutine fill_emf_cart_det
!-----------------------------------------------------------------------------------------

! comp_select defined here
#include "os-cart-sel.f03"

! Now define component functions
#define __TEMPLATE__

#define FNAME(a) a ## _E1
#define __NC__ 1
#define __IE1__ 1
#include __FILE__

#define FNAME(a) a ## _E2
#define __NC__ 1
#define __IE2__ 1
#include __FILE__

#define FNAME(a) a ## _E3
#define __NC__ 1
#define __IE3__ 1
#include __FILE__

#define FNAME(a) a ## _E1E2
#define __NC__ 2
#define __IE1__ 1
#define __IE2__ 2
#include __FILE__

#define FNAME(a) a ## _E1E3
#define __NC__ 2
#define __IE1__ 1
#define __IE3__ 2
#include __FILE__

#define FNAME(a) a ## _E2E3
#define __NC__ 2
#define __IE2__ 1
#define __IE3__ 2
#include __FILE__

#define FNAME(a) a ## _E1E2E3
#define __NC__ 3
#define __IE1__ 1
#define __IE2__ 2
#define __IE3__ 3
#include __FILE__

#ifdef __HAS_RAD_BFLD__

#define FNAME(a) a ## _B1
#define __NC__ 1
#define __IB1__ 1
#include __FILE__

#define FNAME(a) a ## _B2
#define __NC__ 1
#define __IB2__ 1
#include __FILE__

#define FNAME(a) a ## _B3
#define __NC__ 1
#define __IB3__ 1
#include __FILE__

#define FNAME(a) a ## _E1B1
#define __NC__ 2
#define __IE1__ 1
#define __IB1__ 2
#include __FILE__

#define FNAME(a) a ## _E1B2
#define __NC__ 2
#define __IE1__ 1
#define __IB2__ 2
#include __FILE__

#define FNAME(a) a ## _E1B3
#define __NC__ 2
#define __IE1__ 1
#define __IB3__ 2
#include __FILE__

#define FNAME(a) a ## _E2B1
#define __NC__ 2
#define __IE2__ 1
#define __IB1__ 2
#include __FILE__

#define FNAME(a) a ## _E2B2
#define __NC__ 2
#define __IE2__ 1
#define __IB2__ 2
#include __FILE__

#define FNAME(a) a ## _E2B3
#define __NC__ 2
#define __IE2__ 1
#define __IB3__ 2
#include __FILE__

#define FNAME(a) a ## _E3B1
#define __NC__ 2
#define __IE3__ 1
#define __IB1__ 2
#include __FILE__

#define FNAME(a) a ## _E3B2
#define __NC__ 2
#define __IE3__ 1
#define __IB2__ 2
#include __FILE__

#define FNAME(a) a ## _E3B3
#define __NC__ 2
#define __IE3__ 1
#define __IB3__ 2
#include __FILE__

#define FNAME(a) a ## _B1B2
#define __NC__ 2
#define __IB1__ 1
#define __IB2__ 2
#include __FILE__

#define FNAME(a) a ## _B1B3
#define __NC__ 2
#define __IB1__ 1
#define __IB3__ 2
#include __FILE__

#define FNAME(a) a ## _B2B3
#define __NC__ 2
#define __IB2__ 1
#define __IB3__ 2
#include __FILE__

#define FNAME(a) a ## _E1E2B1
#define __NC__ 3
#define __IE1__ 1
#define __IE2__ 2
#define __IB1__ 3
#include __FILE__

#define FNAME(a) a ## _E1E2B2
#define __NC__ 3
#define __IE1__ 1
#define __IE2__ 2
#define __IB2__ 3
#include __FILE__

#define FNAME(a) a ## _E1E2B3
#define __NC__ 3
#define __IE1__ 1
#define __IE2__ 2
#define __IB3__ 3
#include __FILE__

#define FNAME(a) a ## _E1E3B1
#define __NC__ 3
#define __IE1__ 1
#define __IE3__ 2
#define __IB1__ 3
#include __FILE__

#define FNAME(a) a ## _E1E3B2
#define __NC__ 3
#define __IE1__ 1
#define __IE3__ 2
#define __IB2__ 3
#include __FILE__

#define FNAME(a) a ## _E1E3B3
#define __NC__ 3
#define __IE1__ 1
#define __IE3__ 2
#define __IB3__ 3
#include __FILE__

#define FNAME(a) a ## _E1B1B2
#define __NC__ 3
#define __IE1__ 1
#define __IB1__ 2
#define __IB2__ 3
#include __FILE__

#define FNAME(a) a ## _E1B1B3
#define __NC__ 3
#define __IE1__ 1
#define __IB1__ 2
#define __IB3__ 3
#include __FILE__

#define FNAME(a) a ## _E1B2B3
#define __NC__ 3
#define __IE1__ 1
#define __IB2__ 2
#define __IB3__ 3
#include __FILE__

#define FNAME(a) a ## _E2E3B1
#define __NC__ 3
#define __IE2__ 1
#define __IE3__ 2
#define __IB1__ 3
#include __FILE__

#define FNAME(a) a ## _E2E3B2
#define __NC__ 3
#define __IE2__ 1
#define __IE3__ 2
#define __IB2__ 3
#include __FILE__

#define FNAME(a) a ## _E2E3B3
#define __NC__ 3
#define __IE2__ 1
#define __IE3__ 2
#define __IB3__ 3
#include __FILE__

#define FNAME(a) a ## _E2B1B2
#define __NC__ 3
#define __IE2__ 1
#define __IB1__ 2
#define __IB2__ 3
#include __FILE__

#define FNAME(a) a ## _E2B1B3
#define __NC__ 3
#define __IE2__ 1
#define __IB1__ 2
#define __IB3__ 3
#include __FILE__

#define FNAME(a) a ## _E2B2B3
#define __NC__ 3
#define __IE2__ 1
#define __IB2__ 2
#define __IB3__ 3
#include __FILE__

#define FNAME(a) a ## _E3B1B2
#define __NC__ 3
#define __IE3__ 1
#define __IB1__ 2
#define __IB2__ 3
#include __FILE__

#define FNAME(a) a ## _E3B1B3
#define __NC__ 3
#define __IE3__ 1
#define __IB1__ 2
#define __IB3__ 3
#include __FILE__

#define FNAME(a) a ## _E3B2B3
#define __NC__ 3
#define __IE3__ 1
#define __IB2__ 2
#define __IB3__ 3
#include __FILE__

#define FNAME(a) a ## _B1B2B3
#define __NC__ 3
#define __IB1__ 1
#define __IB2__ 2
#define __IB3__ 3
#include __FILE__

#define FNAME(a) a ## _E1E2E3B1
#define __NC__ 4
#define __IE1__ 1
#define __IE2__ 2
#define __IE3__ 3
#define __IB1__ 4
#include __FILE__

#define FNAME(a) a ## _E1E2E3B2
#define __NC__ 4
#define __IE1__ 1
#define __IE2__ 2
#define __IE3__ 3
#define __IB2__ 4
#include __FILE__

#define FNAME(a) a ## _E1E2E3B3
#define __NC__ 4
#define __IE1__ 1
#define __IE2__ 2
#define __IE3__ 3
#define __IB3__ 4
#include __FILE__

#define FNAME(a) a ## _E1E2B1B2
#define __NC__ 4
#define __IE1__ 1
#define __IE2__ 2
#define __IB1__ 3
#define __IB2__ 4
#include __FILE__

#define FNAME(a) a ## _E1E2B1B3
#define __NC__ 4
#define __IE1__ 1
#define __IE2__ 2
#define __IB1__ 3
#define __IB3__ 4
#include __FILE__

#define FNAME(a) a ## _E1E2B2B3
#define __NC__ 4
#define __IE1__ 1
#define __IE2__ 2
#define __IB2__ 3
#define __IB3__ 4
#include __FILE__

#define FNAME(a) a ## _E1E3B1B2
#define __NC__ 4
#define __IE1__ 1
#define __IE3__ 2
#define __IB1__ 3
#define __IB2__ 4
#include __FILE__

#define FNAME(a) a ## _E1E3B1B3
#define __NC__ 4
#define __IE1__ 1
#define __IE3__ 2
#define __IB1__ 3
#define __IB3__ 4
#include __FILE__

#define FNAME(a) a ## _E1E3B2B3
#define __NC__ 4
#define __IE1__ 1
#define __IE3__ 2
#define __IB2__ 3
#define __IB3__ 4
#include __FILE__

#define FNAME(a) a ## _E1B1B2B3
#define __NC__ 4
#define __IE1__ 1
#define __IB1__ 2
#define __IB2__ 3
#define __IB3__ 4
#include __FILE__

#define FNAME(a) a ## _E2E3B1B2
#define __NC__ 4
#define __IE2__ 1
#define __IE3__ 2
#define __IB1__ 3
#define __IB2__ 4
#include __FILE__

#define FNAME(a) a ## _E2E3B1B3
#define __NC__ 4
#define __IE2__ 1
#define __IE3__ 2
#define __IB1__ 3
#define __IB3__ 4
#include __FILE__

#define FNAME(a) a ## _E2E3B2B3
#define __NC__ 4
#define __IE2__ 1
#define __IE3__ 2
#define __IB2__ 3
#define __IB3__ 4
#include __FILE__

#define FNAME(a) a ## _E2B1B2B3
#define __NC__ 4
#define __IE2__ 1
#define __IB1__ 2
#define __IB2__ 3
#define __IB3__ 4
#include __FILE__

#define FNAME(a) a ## _E3B1B2B3
#define __NC__ 4
#define __IE3__ 1
#define __IB1__ 2
#define __IB2__ 3
#define __IB3__ 4
#include __FILE__

#define FNAME(a) a ## _E1E2E3B1B2
#define __NC__ 5
#define __IE1__ 1
#define __IE2__ 2
#define __IE3__ 3
#define __IB1__ 4
#define __IB2__ 5
#include __FILE__

#define FNAME(a) a ## _E1E2E3B1B3
#define __NC__ 5
#define __IE1__ 1
#define __IE2__ 2
#define __IE3__ 3
#define __IB1__ 4
#define __IB3__ 5
#include __FILE__

#define FNAME(a) a ## _E1E2E3B2B3
#define __NC__ 5
#define __IE1__ 1
#define __IE2__ 2
#define __IE3__ 3
#define __IB2__ 4
#define __IB3__ 5
#include __FILE__

#define FNAME(a) a ## _E1E2B1B2B3
#define __NC__ 5
#define __IE1__ 1
#define __IE2__ 2
#define __IB1__ 3
#define __IB2__ 4
#define __IB3__ 5
#include __FILE__

#define FNAME(a) a ## _E1E3B1B2B3
#define __NC__ 5
#define __IE1__ 1
#define __IE3__ 2
#define __IB1__ 3
#define __IB2__ 4
#define __IB3__ 5
#include __FILE__

#define FNAME(a) a ## _E2E3B1B2B3
#define __NC__ 5
#define __IE2__ 1
#define __IE3__ 2
#define __IB1__ 3
#define __IB2__ 4
#define __IB3__ 5
#include __FILE__

#define FNAME(a) a ## _E1E2E3B1B2B3
#define __NC__ 6
#define __IE1__ 1
#define __IE2__ 2
#define __IE3__ 3
#define __IB1__ 4
#define __IB2__ 5
#define __IB3__ 6
#include __FILE__

#endif

end module class_CartDet


#else

!-----------------------------------------------------------------------------------------
subroutine FNAME(calc_comp)(this,x,beta,time,xprev,betaprev,charge,np,tid)
  implicit none
  class(cart_Detector), intent(inout) :: this
  real(p_k_part), dimension(3,p_cache_size), intent(in) :: x, xprev
  real(p_k_part), dimension(p_p_dim,p_cache_size), intent(in) :: beta, betaprev
  real(p_k_part), dimension(p_cache_size), intent(in) :: charge
  real(p_double), intent(in) :: time
  integer, intent(in) :: np, tid

  ! n is the direction vector and R is the distance
  ! tdep time of deposition
  real(p_k_part) :: Emf(__NC__,p_cache_size), n(3,p_cache_size), xyzpos(3)
  real(p_double) :: tdep(p_cache_size), tdepprev(p_cache_size), tprev, rsim_dt, temp
  real(p_double) :: tdep_c(p_cache_size), tdepprev_c(p_cache_size)
  real(p_k_part) :: betaprev_c(p_p_dim,p_cache_size), betadot_c(p_p_dim,p_cache_size)
  real(p_k_part) :: R(p_cache_size), rR, betadot(p_p_dim,p_cache_size)
  integer :: tit, it, titprev, c, i, l, npp
  integer :: iters(3), iter1, iter2, nf1
  integer(p_int64) :: iit, nf1_long
  real(p_k_part) :: div_charge(p_cache_size)
#if defined (__IE1__) || defined (__IE2__) || defined (__IE3__)
  real(p_k_part), dimension(p_cache_size) :: nmb1, nmb2, nmb3
#endif
#ifdef __IB1__
  real(p_k_part), dimension(p_cache_size) :: onemb1n1, onemn12
#endif
#ifdef __IB2__
  real(p_k_part), dimension(p_cache_size) :: onemb2n2, onemn22
#endif
#ifdef __IB3__
  real(p_k_part), dimension(p_cache_size) :: onemb3n3, onemn32
#endif

  rsim_dt = 1.0_p_double / this%sim_dt
  do i = 1, np
    betadot(1,i) = (beta(1,i)-betaprev(1,i)) * rsim_dt
    betadot(2,i) = (beta(2,i)-betaprev(2,i)) * rsim_dt
    betadot(3,i) = (beta(3,i)-betaprev(3,i)) * rsim_dt
  enddo

  tprev=time-this%sim_dt
  iters=1
  nf1 = this%nf(this%dim_pos(1))
  nf1_long = int( nf1, p_int64 )
  ! Set initial values of iteration variables
  iit = this%iter_bnds(p_lower,tid+1)
  iter1 = int( mod( iit-1, nf1_long ) ) + 1
  iter2 = int( (iit-1) / nf1_long ) + 1

  do iit = this%iter_bnds(p_lower,tid+1), this%iter_bnds(p_upper,tid+1)

    iters(this%dim_pos(1)-1)=iter1
    iters(this%dim_pos(2)-1)=iter2

    xyzpos(1) = this%xpos(iters(1))
    xyzpos(2) = this%ypos(iters(2))
    xyzpos(3) = this%zpos(iters(3))

    !Compute the non-normalized direction vector x-pos
    do i = 1, np
      n(1,i) = x(1,i) - xyzpos(1)
      n(2,i) = x(2,i) - xyzpos(2)
      n(3,i) = x(3,i) - xyzpos(3)
    enddo

    !Compute the Distance
    do i = 1, np
      R(i) = sqrt( n(1,i)**2 + n(2,i)**2 + n(3,i)**2 )
    enddo
    do i = 1, np
      tdep(i) = time + R(i)
    enddo

    !Compute the non-normalized direction vector x-pos for prev
    do i = 1, np
      n(1,i) = xprev(1,i) - xyzpos(1)
      n(2,i) = xprev(2,i) - xyzpos(2)
      n(3,i) = xprev(3,i) - xyzpos(3)
    enddo

    !Compute the Distance
    do i = 1, np
      R(i) = sqrt( n(1,i)**2 + n(2,i)**2 + n(3,i)**2 )
    enddo
    do i = 1, np
      tdepprev(i) = tprev + R(i)
    enddo

    l = 1
    do i = 1, np
      if ((tdepprev(i) >= this%tpos(1)) .and. (tdep(i) <= this%tpos(this%n1+1))) then

        rR = 1.0_p_k_part / R(i)
        n(1,l) = -n(1,i) * rR
        n(2,l) = -n(2,i) * rR
        n(3,l) = -n(3,i) * rR
        div_charge(l) = charge(i) * rR ! will be modified again later
        betaprev_c(1,l) = betaprev(1,i)
        betaprev_c(2,l) = betaprev(2,i)
        betaprev_c(3,l) = betaprev(3,i)
        betadot_c(1,l) = betadot(1,i)
        betadot_c(2,l) = betadot(2,i)
        betadot_c(3,l) = betadot(3,i)
        tdep_c(l) = tdep(i)
        tdepprev_c(l) = tdepprev(i)

        l = l + 1
      endif
    enddo
    npp = l - 1

    do i = 1, npp
#if defined (__IE1__) || defined (__IE2__) || defined (__IE3__)
      nmb1(i)=n(1,i)-betaprev_c(1,i)
      nmb2(i)=n(2,i)-betaprev_c(2,i)
      nmb3(i)=n(3,i)-betaprev_c(3,i)
#endif
#ifdef __IB1__
      onemb1n1(i)=1-betaprev_c(1,i)*n(1,i)
      onemn12(i)=1-n(1,i)**2
#endif
#ifdef __IB2__
      onemb2n2(i)=1-betaprev_c(2,i)*n(2,i)
      onemn22(i)=1-n(2,i)**2
#endif
#ifdef __IB3__
      onemb3n3(i)=1-betaprev_c(3,i)*n(3,i)
      onemn32(i)=1-n(3,i)**2
#endif
      div_charge(i)=div_charge(i)/((1-betaprev_c(1,i)*n(1,i)-betaprev_c(2,i)*n(2,i)&
                                     -betaprev_c(3,i)*n(3,i))**3)
    enddo

    do i = 1, npp
#ifdef __IE1__
      Emf(__IE1__,i)=div_charge(i)*(n(2,i)*(betadot_c(2,i)*nmb1(i)-betadot_c(1,i)*nmb2(i)) &
            +n(3,i)*(betadot_c(3,i)*nmb1(i)-betadot_c(1,i)*nmb3(i)))
#endif
#ifdef __IE2__
      Emf(__IE2__,i)=div_charge(i)*(nmb2(i)*(betadot_c(1,i)*n(1,i)+betadot_c(3,i)*n(3,i)) &
             + betadot_c(2,i)*(-nmb1(i)*n(1,i)-nmb3(i)*n(3,i)))
#endif
#ifdef __IE3__
      Emf(__IE3__,i)=div_charge(i)*(betadot_c(3,i)*(-nmb1(i)*n(1,i)-nmb2(i)*n(2,i)) &
             + nmb3(i)*(betadot_c(1,i)*n(1,i)+betadot_c(2,i)*n(2,i)))
#endif
#ifdef __IB1__
      Emf(__IB1__,i)=div_charge(i)*(betadot_c(3,i)*(betaprev_c(2,i)*onemn12(i)-n(2,i)*onemb1n1(i)) &
            +betadot_c(1,i)*n(1,i)*(betaprev_c(2,i)*n(3,i)-betaprev_c(3,i)*n(2,i)) &
            +betadot_c(2,i)*(n(3,i)*onemb1n1(i)-betaprev_c(3,i)*onemn12(i)))
#endif
#ifdef __IB2__
      Emf(__IB2__,i)=div_charge(i)*(betadot_c(3,i)*(n(1,i)*onemb2n2(i)-betaprev_c(1,i)*onemn22(i)) &
            +betadot_c(2,i)*n(2,i)*(betaprev_c(3,i)*n(1,i)-betaprev_c(1,i)*n(3,i)) &
            +betadot_c(1,i)*(betaprev_c(3,i)*onemn22(i)-n(3,i)*onemb2n2(i)))
#endif
#ifdef __IB3__
      Emf(__IB3__,i)=div_charge(i)*(betadot_c(2,i)*(betaprev_c(1,i)*onemn32(i)-n(1,i)*onemb3n3(i)) &
            +betadot_c(3,i)*n(3,i)*(betaprev_c(1,i)*n(2,i)-betaprev_c(2,i)*n(1,i)) &
            +betadot_c(1,i)*(n(2,i)*onemb3n3(i)-betaprev_c(2,i)*onemn32(i)))
#endif
    enddo

    do i = 1, npp
      tit=int((tdep_c(i)-this%tpos(1))*this%d1m1)+1
      titprev=int((tdepprev_c(i)-this%tpos(1))*this%d1m1)+1
      !interpolation shenanigans. refer to report for better understanding
      do it=titprev,tit-1
        temp = (this%tpos(it+1)-tdepprev_c(i))
        do c = 1, __NC__
          this%Emf(c,iter1,iter2,it) = this%Emf(c,iter1,iter2,it) &
                                     + real(Emf(c,i)*temp,p_single)
        enddo
        tdepprev_c(i)=this%tpos(it+1)
      enddo
      temp = (tdep_c(i)-tdepprev_c(i))
      do c = 1, __NC__
        this%Emf(c,iter1,iter2,tit) = this%Emf(c,iter1,iter2,tit) &
                                    + real(Emf(c,i)*temp,p_single)
      enddo
    enddo

    iter1 = iter1 + 1
    if ( iter1 > nf1 ) then
      iter1 = iter1 - nf1
      iter2 = iter2 + 1
    endif

  enddo

end subroutine FNAME(calc_comp)
!-----------------------------------------------------------------------------------------

#undef __NC__
#undef __IE1__
#undef __IE2__
#undef __IE3__
#undef __IB1__
#undef __IB2__
#undef __IB3__
#undef FNAME

#endif
