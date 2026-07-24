# 1 "radio/os-radio-cartdet.f03"
# 1 "<built-in>" 1
# 1 "<built-in>" 3
# 467 "<built-in>" 3
# 1 "<command line>" 1
# 1 "<built-in>" 2
# 1 "radio/os-radio-cartdet.f03" 2
! if __TEMPLATE__ is defined then read the template definition at the end of the file


module class_CartDet

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
# 7 "radio/os-radio-cartdet.f03" 2
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
# 8 "radio/os-radio-cartdet.f03" 2
# 1 "./memory/memory.h" 1
! Include file for the memory module
! Files must include this file using #include "memory/memory.h" rather than just using the module
! the #include must be placed where the module use statement would usually be i.e.
!
! module module2
!
! use module1
! #include "memory/memory.h"
!
# 19 "./memory/memory.h"
use memory
# 9 "radio/os-radio-cartdet.f03" 2

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

  logical :: has_b
  integer :: i

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

  call alloc(this%dims,(/ this%ndims /),"radio/os-radio-cartdet.f03",165)
  call alloc(this%dmin_,(/ this%ndims /),"radio/os-radio-cartdet.f03",166)
  call alloc(this%dmax_,(/ this%ndims /),"radio/os-radio-cartdet.f03",167)

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
  bas_ax_nms =(/"x1","x2","x3","t "/)
  bas_ax_lnms =(/"x_2","x_2","x_3","t  "/)
  bas_ax_uns =(/"c/\omega_p","c/\omega_p","c/\omega_p","1/\omega_p" /)

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
  this%n2 =this%nf(2) !x
  this%n3 =this%nf(3) !y
  this%n4 =this%nf(4) !z
  this%n1 =this%nf(1) !t

  this%d1 =(this%dmax(1)-this%dmin(1))
  this%d2 =(this%dmax(2)-this%dmin(2))
  this%d3 =(this%dmax(3)-this%dmin(3))
  this%d4 =(this%dmax(4)-this%dmin(4))

  if(this%nf(1)>1) this%d1 = (this%dmax(1)-this%dmin(1))/(this%nf(1)-1)
  if(this%nf(2)>1) this%d2 = (this%dmax(2)-this%dmin(2))/(this%nf(2)-1)
  if(this%nf(3)>1) this%d3 = (this%dmax(3)-this%dmin(3))/(this%nf(3)-1)
  if(this%nf(4)>1) this%d4 = (this%dmax(4)-this%dmin(4))/(this%nf(4)-1)
  this%d1m1=1/this%d1

  call alloc(this%xpos,(/this%nf(2)/),"radio/os-radio-cartdet.f03",230)
  call alloc(this%ypos,(/this%nf(3)/),"radio/os-radio-cartdet.f03",231)
  call alloc(this%zpos,(/this%nf(4)/),"radio/os-radio-cartdet.f03",232)
  call alloc(this%tpos,(/this%nf(1)+1/),"radio/os-radio-cartdet.f03",233)

  if (this%ncmp>=1) then
    dp1 = this%dim_pos(1)
    dp2 = this%dim_pos(2)

    call alloc(this%Emf, (/ this%ncmp, this%nf(dp1), this%nf(dp2), this%nf(1)+1 /),"radio/os-radio-cartdet.f03",239)
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

  call alloc(this%iter_bnds, (/ 2, nt /),"radio/os-radio-cartdet.f03",263)
  do i = 0, nt-1
    chunk = ( ncells + nt - 1 ) / nt
    i0 = i * chunk + 1
    i1 = min( (i+1) * chunk, ncells )
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
  ! do j
  ! x(j,i)

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

  call freemem(this%xpos,"radio/os-radio-cartdet.f03",316)
  call freemem(this%ypos,"radio/os-radio-cartdet.f03",317)
  call freemem(this%zpos,"radio/os-radio-cartdet.f03",318)
  call freemem(this%tpos,"radio/os-radio-cartdet.f03",319)
  call freemem(this%Emf,"radio/os-radio-cartdet.f03",320)
  call freemem(this%dims,"radio/os-radio-cartdet.f03",321)
  call freemem(this%dmin_,"radio/os-radio-cartdet.f03",322)
  call freemem(this%dmax_,"radio/os-radio-cartdet.f03",323)
  call freemem(this%iter_bnds,"radio/os-radio-cartdet.f03",324)

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
  call restart_io_write("p_cart_det_rst_id", p_cart_det_rst_id, restart_handle, ierr)
  call check_error(ierr,err_msg,p_err_rstwrt,"radio/os-radio-cartdet.f03",342)

  call restart_io_write("this%Emf", this%Emf, restart_handle, ierr)
  call check_error(ierr,err_msg,p_err_rstwrt,"radio/os-radio-cartdet.f03",345)

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
  call restart_io_read(rst_id, restart_handle, ierr)
  call check_error(ierr,err_msg,p_err_rstrd,"radio/os-radio-cartdet.f03",364)

  ! check if restart file is compatible
  if ( rst_id /= p_cart_det_rst_id ) then
    write(err_buf__,*) 'Corrupted restart file, or restart file ';call err__("radio/os-radio-cartdet.f03",368)
    write(err_buf__,*) 'from incompatible binary (cart_det)';call err__("radio/os-radio-cartdet.f03",369)
    call abort_program(p_err_rstrd)
  endif

  call restart_io_read(this%Emf, restart_handle, ierr)
  call check_error(ierr,err_msg,p_err_rstrd,"radio/os-radio-cartdet.f03",374)

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
# 1 "radio/os-cart-sel.f03" 1
!-----------------------------------------------------------------------------------------
subroutine calc_dflt(this,x,beta,time,xprev,betaprev,charge,np,tid)
  implicit none
  class(cart_Detector), intent(inout) :: this
  real(p_k_part), dimension(3,p_cache_size), intent(in) :: x, xprev
  real(p_k_part), dimension(p_p_dim,p_cache_size), intent(in) :: beta, betaprev
  real(p_k_part), dimension(p_cache_size), intent(in) :: charge
  real(p_double), intent(in) :: time
  integer, intent(in) :: np, tid

  ! Do nothing

end subroutine calc_dflt
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine comp_select_cart_det( this )

  implicit none

  class(cart_Detector), intent(inout) :: this

  integer :: sele

  sele = 0
  if (this%cmpts(6)) sele = sele + 32
  if (this%cmpts(5)) sele = sele + 16
  if (this%cmpts(4)) sele = sele + 8
  if (this%cmpts(3)) sele = sele + 4
  if (this%cmpts(2)) sele = sele + 2
  if (this%cmpts(1)) sele = sele + 1

  select case(sele)
  case(1)
    this % fill_emf_loc => calc_comp_E1
  case(2)
    this % fill_emf_loc => calc_comp_E2
  case(3)
    this % fill_emf_loc => calc_comp_E1E2
  case(4)
    this % fill_emf_loc => calc_comp_E3
  case(5)
    this % fill_emf_loc => calc_comp_E1E3
  case(6)
    this % fill_emf_loc => calc_comp_E2E3
  case(7)
    this % fill_emf_loc => calc_comp_E1E2E3
# 162 "radio/os-cart-sel.f03"
  case default
    this % fill_emf_loc => calc_dflt
  end select

end subroutine comp_select_cart_det
!-----------------------------------------------------------------------------------------
# 397 "radio/os-radio-cartdet.f03" 2

! Now define component functions





# 1 "radio/os-radio-cartdet.f03" 1
! if is defined then read the template definition at the end of the file
# 854 "radio/os-radio-cartdet.f03"
!-----------------------------------------------------------------------------------------
subroutine calc_comp_E1(this,x,beta,time,xprev,betaprev,charge,np,tid)
  implicit none
  class(cart_Detector), intent(inout) :: this
  real(p_k_part), dimension(3,p_cache_size), intent(in) :: x, xprev
  real(p_k_part), dimension(p_p_dim,p_cache_size), intent(in) :: beta, betaprev
  real(p_k_part), dimension(p_cache_size), intent(in) :: charge
  real(p_double), intent(in) :: time
  integer, intent(in) :: np, tid

  ! n is the direction vector and R is the distance
  ! tdep time of deposition
  real(p_k_part) :: Emf(1,p_cache_size), n(3,p_cache_size), xyzpos(3)
  real(p_double) :: tdep(p_cache_size), tdepprev(p_cache_size), tprev, rsim_dt, temp
  real(p_double) :: tdep_c(p_cache_size), tdepprev_c(p_cache_size)
  real(p_k_part) :: betaprev_c(p_p_dim,p_cache_size), betadot_c(p_p_dim,p_cache_size)
  real(p_k_part) :: R(p_cache_size), rR, betadot(p_p_dim,p_cache_size)
  integer :: tit, it, titprev, c, i, l, npp
  integer :: iters(3), iter1, iter2, nf1
  integer(p_int64) :: iit, nf1_long
  real(p_k_part) :: div_charge(p_cache_size)

  real(p_k_part), dimension(p_cache_size) :: nmb1, nmb2, nmb3
# 888 "radio/os-radio-cartdet.f03"
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

      nmb1(i)=n(1,i)-betaprev_c(1,i)
      nmb2(i)=n(2,i)-betaprev_c(2,i)
      nmb3(i)=n(3,i)-betaprev_c(3,i)
# 984 "radio/os-radio-cartdet.f03"
      div_charge(i)=div_charge(i)/((1-betaprev_c(1,i)*n(1,i)-betaprev_c(2,i)*n(2,i)&
                                     -betaprev_c(3,i)*n(3,i))**3)
    enddo

    do i = 1, npp

      Emf(1,i)=div_charge(i)*(n(2,i)*(betadot_c(2,i)*nmb1(i)-betadot_c(1,i)*nmb2(i)) &
            +n(3,i)*(betadot_c(3,i)*nmb1(i)-betadot_c(1,i)*nmb3(i)))
# 1016 "radio/os-radio-cartdet.f03"
    enddo

    do i = 1, npp
      tit=int((tdep_c(i)-this%tpos(1))*this%d1m1)+1
      titprev=int((tdepprev_c(i)-this%tpos(1))*this%d1m1)+1
      !interpolation shenanigans. refer to report for better understanding
      do it=titprev,tit-1
        temp = (this%tpos(it+1)-tdepprev_c(i))
        do c = 1, 1
          this%Emf(c,iter1,iter2,it) = this%Emf(c,iter1,iter2,it) &
                                     + real(Emf(c,i)*temp,p_single)
        enddo
        tdepprev_c(i)=this%tpos(it+1)
      enddo
      temp = (tdep_c(i)-tdepprev_c(i))
      do c = 1, 1
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

end subroutine calc_comp_E1
!-----------------------------------------------------------------------------------------
# 405 "radio/os-radio-cartdet.f03" 2




# 1 "radio/os-radio-cartdet.f03" 1
! if is defined then read the template definition at the end of the file
# 854 "radio/os-radio-cartdet.f03"
!-----------------------------------------------------------------------------------------
subroutine calc_comp_E2(this,x,beta,time,xprev,betaprev,charge,np,tid)
  implicit none
  class(cart_Detector), intent(inout) :: this
  real(p_k_part), dimension(3,p_cache_size), intent(in) :: x, xprev
  real(p_k_part), dimension(p_p_dim,p_cache_size), intent(in) :: beta, betaprev
  real(p_k_part), dimension(p_cache_size), intent(in) :: charge
  real(p_double), intent(in) :: time
  integer, intent(in) :: np, tid

  ! n is the direction vector and R is the distance
  ! tdep time of deposition
  real(p_k_part) :: Emf(1,p_cache_size), n(3,p_cache_size), xyzpos(3)
  real(p_double) :: tdep(p_cache_size), tdepprev(p_cache_size), tprev, rsim_dt, temp
  real(p_double) :: tdep_c(p_cache_size), tdepprev_c(p_cache_size)
  real(p_k_part) :: betaprev_c(p_p_dim,p_cache_size), betadot_c(p_p_dim,p_cache_size)
  real(p_k_part) :: R(p_cache_size), rR, betadot(p_p_dim,p_cache_size)
  integer :: tit, it, titprev, c, i, l, npp
  integer :: iters(3), iter1, iter2, nf1
  integer(p_int64) :: iit, nf1_long
  real(p_k_part) :: div_charge(p_cache_size)

  real(p_k_part), dimension(p_cache_size) :: nmb1, nmb2, nmb3
# 888 "radio/os-radio-cartdet.f03"
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

      nmb1(i)=n(1,i)-betaprev_c(1,i)
      nmb2(i)=n(2,i)-betaprev_c(2,i)
      nmb3(i)=n(3,i)-betaprev_c(3,i)
# 984 "radio/os-radio-cartdet.f03"
      div_charge(i)=div_charge(i)/((1-betaprev_c(1,i)*n(1,i)-betaprev_c(2,i)*n(2,i)&
                                     -betaprev_c(3,i)*n(3,i))**3)
    enddo

    do i = 1, npp





      Emf(1,i)=div_charge(i)*(nmb2(i)*(betadot_c(1,i)*n(1,i)+betadot_c(3,i)*n(3,i)) &
             + betadot_c(2,i)*(-nmb1(i)*n(1,i)-nmb3(i)*n(3,i)))
# 1016 "radio/os-radio-cartdet.f03"
    enddo

    do i = 1, npp
      tit=int((tdep_c(i)-this%tpos(1))*this%d1m1)+1
      titprev=int((tdepprev_c(i)-this%tpos(1))*this%d1m1)+1
      !interpolation shenanigans. refer to report for better understanding
      do it=titprev,tit-1
        temp = (this%tpos(it+1)-tdepprev_c(i))
        do c = 1, 1
          this%Emf(c,iter1,iter2,it) = this%Emf(c,iter1,iter2,it) &
                                     + real(Emf(c,i)*temp,p_single)
        enddo
        tdepprev_c(i)=this%tpos(it+1)
      enddo
      temp = (tdep_c(i)-tdepprev_c(i))
      do c = 1, 1
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

end subroutine calc_comp_E2
!-----------------------------------------------------------------------------------------
# 410 "radio/os-radio-cartdet.f03" 2




# 1 "radio/os-radio-cartdet.f03" 1
! if is defined then read the template definition at the end of the file
# 854 "radio/os-radio-cartdet.f03"
!-----------------------------------------------------------------------------------------
subroutine calc_comp_E3(this,x,beta,time,xprev,betaprev,charge,np,tid)
  implicit none
  class(cart_Detector), intent(inout) :: this
  real(p_k_part), dimension(3,p_cache_size), intent(in) :: x, xprev
  real(p_k_part), dimension(p_p_dim,p_cache_size), intent(in) :: beta, betaprev
  real(p_k_part), dimension(p_cache_size), intent(in) :: charge
  real(p_double), intent(in) :: time
  integer, intent(in) :: np, tid

  ! n is the direction vector and R is the distance
  ! tdep time of deposition
  real(p_k_part) :: Emf(1,p_cache_size), n(3,p_cache_size), xyzpos(3)
  real(p_double) :: tdep(p_cache_size), tdepprev(p_cache_size), tprev, rsim_dt, temp
  real(p_double) :: tdep_c(p_cache_size), tdepprev_c(p_cache_size)
  real(p_k_part) :: betaprev_c(p_p_dim,p_cache_size), betadot_c(p_p_dim,p_cache_size)
  real(p_k_part) :: R(p_cache_size), rR, betadot(p_p_dim,p_cache_size)
  integer :: tit, it, titprev, c, i, l, npp
  integer :: iters(3), iter1, iter2, nf1
  integer(p_int64) :: iit, nf1_long
  real(p_k_part) :: div_charge(p_cache_size)

  real(p_k_part), dimension(p_cache_size) :: nmb1, nmb2, nmb3
# 888 "radio/os-radio-cartdet.f03"
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

      nmb1(i)=n(1,i)-betaprev_c(1,i)
      nmb2(i)=n(2,i)-betaprev_c(2,i)
      nmb3(i)=n(3,i)-betaprev_c(3,i)
# 984 "radio/os-radio-cartdet.f03"
      div_charge(i)=div_charge(i)/((1-betaprev_c(1,i)*n(1,i)-betaprev_c(2,i)*n(2,i)&
                                     -betaprev_c(3,i)*n(3,i))**3)
    enddo

    do i = 1, npp
# 998 "radio/os-radio-cartdet.f03"
      Emf(1,i)=div_charge(i)*(betadot_c(3,i)*(-nmb1(i)*n(1,i)-nmb2(i)*n(2,i)) &
             + nmb3(i)*(betadot_c(1,i)*n(1,i)+betadot_c(2,i)*n(2,i)))
# 1016 "radio/os-radio-cartdet.f03"
    enddo

    do i = 1, npp
      tit=int((tdep_c(i)-this%tpos(1))*this%d1m1)+1
      titprev=int((tdepprev_c(i)-this%tpos(1))*this%d1m1)+1
      !interpolation shenanigans. refer to report for better understanding
      do it=titprev,tit-1
        temp = (this%tpos(it+1)-tdepprev_c(i))
        do c = 1, 1
          this%Emf(c,iter1,iter2,it) = this%Emf(c,iter1,iter2,it) &
                                     + real(Emf(c,i)*temp,p_single)
        enddo
        tdepprev_c(i)=this%tpos(it+1)
      enddo
      temp = (tdep_c(i)-tdepprev_c(i))
      do c = 1, 1
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

end subroutine calc_comp_E3
!-----------------------------------------------------------------------------------------
# 415 "radio/os-radio-cartdet.f03" 2





# 1 "radio/os-radio-cartdet.f03" 1
! if is defined then read the template definition at the end of the file
# 854 "radio/os-radio-cartdet.f03"
!-----------------------------------------------------------------------------------------
subroutine calc_comp_E1E2(this,x,beta,time,xprev,betaprev,charge,np,tid)
  implicit none
  class(cart_Detector), intent(inout) :: this
  real(p_k_part), dimension(3,p_cache_size), intent(in) :: x, xprev
  real(p_k_part), dimension(p_p_dim,p_cache_size), intent(in) :: beta, betaprev
  real(p_k_part), dimension(p_cache_size), intent(in) :: charge
  real(p_double), intent(in) :: time
  integer, intent(in) :: np, tid

  ! n is the direction vector and R is the distance
  ! tdep time of deposition
  real(p_k_part) :: Emf(2,p_cache_size), n(3,p_cache_size), xyzpos(3)
  real(p_double) :: tdep(p_cache_size), tdepprev(p_cache_size), tprev, rsim_dt, temp
  real(p_double) :: tdep_c(p_cache_size), tdepprev_c(p_cache_size)
  real(p_k_part) :: betaprev_c(p_p_dim,p_cache_size), betadot_c(p_p_dim,p_cache_size)
  real(p_k_part) :: R(p_cache_size), rR, betadot(p_p_dim,p_cache_size)
  integer :: tit, it, titprev, c, i, l, npp
  integer :: iters(3), iter1, iter2, nf1
  integer(p_int64) :: iit, nf1_long
  real(p_k_part) :: div_charge(p_cache_size)

  real(p_k_part), dimension(p_cache_size) :: nmb1, nmb2, nmb3
# 888 "radio/os-radio-cartdet.f03"
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

      nmb1(i)=n(1,i)-betaprev_c(1,i)
      nmb2(i)=n(2,i)-betaprev_c(2,i)
      nmb3(i)=n(3,i)-betaprev_c(3,i)
# 984 "radio/os-radio-cartdet.f03"
      div_charge(i)=div_charge(i)/((1-betaprev_c(1,i)*n(1,i)-betaprev_c(2,i)*n(2,i)&
                                     -betaprev_c(3,i)*n(3,i))**3)
    enddo

    do i = 1, npp

      Emf(1,i)=div_charge(i)*(n(2,i)*(betadot_c(2,i)*nmb1(i)-betadot_c(1,i)*nmb2(i)) &
            +n(3,i)*(betadot_c(3,i)*nmb1(i)-betadot_c(1,i)*nmb3(i)))


      Emf(2,i)=div_charge(i)*(nmb2(i)*(betadot_c(1,i)*n(1,i)+betadot_c(3,i)*n(3,i)) &
             + betadot_c(2,i)*(-nmb1(i)*n(1,i)-nmb3(i)*n(3,i)))
# 1016 "radio/os-radio-cartdet.f03"
    enddo

    do i = 1, npp
      tit=int((tdep_c(i)-this%tpos(1))*this%d1m1)+1
      titprev=int((tdepprev_c(i)-this%tpos(1))*this%d1m1)+1
      !interpolation shenanigans. refer to report for better understanding
      do it=titprev,tit-1
        temp = (this%tpos(it+1)-tdepprev_c(i))
        do c = 1, 2
          this%Emf(c,iter1,iter2,it) = this%Emf(c,iter1,iter2,it) &
                                     + real(Emf(c,i)*temp,p_single)
        enddo
        tdepprev_c(i)=this%tpos(it+1)
      enddo
      temp = (tdep_c(i)-tdepprev_c(i))
      do c = 1, 2
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

end subroutine calc_comp_E1E2
!-----------------------------------------------------------------------------------------
# 421 "radio/os-radio-cartdet.f03" 2





# 1 "radio/os-radio-cartdet.f03" 1
! if is defined then read the template definition at the end of the file
# 854 "radio/os-radio-cartdet.f03"
!-----------------------------------------------------------------------------------------
subroutine calc_comp_E1E3(this,x,beta,time,xprev,betaprev,charge,np,tid)
  implicit none
  class(cart_Detector), intent(inout) :: this
  real(p_k_part), dimension(3,p_cache_size), intent(in) :: x, xprev
  real(p_k_part), dimension(p_p_dim,p_cache_size), intent(in) :: beta, betaprev
  real(p_k_part), dimension(p_cache_size), intent(in) :: charge
  real(p_double), intent(in) :: time
  integer, intent(in) :: np, tid

  ! n is the direction vector and R is the distance
  ! tdep time of deposition
  real(p_k_part) :: Emf(2,p_cache_size), n(3,p_cache_size), xyzpos(3)
  real(p_double) :: tdep(p_cache_size), tdepprev(p_cache_size), tprev, rsim_dt, temp
  real(p_double) :: tdep_c(p_cache_size), tdepprev_c(p_cache_size)
  real(p_k_part) :: betaprev_c(p_p_dim,p_cache_size), betadot_c(p_p_dim,p_cache_size)
  real(p_k_part) :: R(p_cache_size), rR, betadot(p_p_dim,p_cache_size)
  integer :: tit, it, titprev, c, i, l, npp
  integer :: iters(3), iter1, iter2, nf1
  integer(p_int64) :: iit, nf1_long
  real(p_k_part) :: div_charge(p_cache_size)

  real(p_k_part), dimension(p_cache_size) :: nmb1, nmb2, nmb3
# 888 "radio/os-radio-cartdet.f03"
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

      nmb1(i)=n(1,i)-betaprev_c(1,i)
      nmb2(i)=n(2,i)-betaprev_c(2,i)
      nmb3(i)=n(3,i)-betaprev_c(3,i)
# 984 "radio/os-radio-cartdet.f03"
      div_charge(i)=div_charge(i)/((1-betaprev_c(1,i)*n(1,i)-betaprev_c(2,i)*n(2,i)&
                                     -betaprev_c(3,i)*n(3,i))**3)
    enddo

    do i = 1, npp

      Emf(1,i)=div_charge(i)*(n(2,i)*(betadot_c(2,i)*nmb1(i)-betadot_c(1,i)*nmb2(i)) &
            +n(3,i)*(betadot_c(3,i)*nmb1(i)-betadot_c(1,i)*nmb3(i)))






      Emf(2,i)=div_charge(i)*(betadot_c(3,i)*(-nmb1(i)*n(1,i)-nmb2(i)*n(2,i)) &
             + nmb3(i)*(betadot_c(1,i)*n(1,i)+betadot_c(2,i)*n(2,i)))
# 1016 "radio/os-radio-cartdet.f03"
    enddo

    do i = 1, npp
      tit=int((tdep_c(i)-this%tpos(1))*this%d1m1)+1
      titprev=int((tdepprev_c(i)-this%tpos(1))*this%d1m1)+1
      !interpolation shenanigans. refer to report for better understanding
      do it=titprev,tit-1
        temp = (this%tpos(it+1)-tdepprev_c(i))
        do c = 1, 2
          this%Emf(c,iter1,iter2,it) = this%Emf(c,iter1,iter2,it) &
                                     + real(Emf(c,i)*temp,p_single)
        enddo
        tdepprev_c(i)=this%tpos(it+1)
      enddo
      temp = (tdep_c(i)-tdepprev_c(i))
      do c = 1, 2
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

end subroutine calc_comp_E1E3
!-----------------------------------------------------------------------------------------
# 427 "radio/os-radio-cartdet.f03" 2





# 1 "radio/os-radio-cartdet.f03" 1
! if is defined then read the template definition at the end of the file
# 854 "radio/os-radio-cartdet.f03"
!-----------------------------------------------------------------------------------------
subroutine calc_comp_E2E3(this,x,beta,time,xprev,betaprev,charge,np,tid)
  implicit none
  class(cart_Detector), intent(inout) :: this
  real(p_k_part), dimension(3,p_cache_size), intent(in) :: x, xprev
  real(p_k_part), dimension(p_p_dim,p_cache_size), intent(in) :: beta, betaprev
  real(p_k_part), dimension(p_cache_size), intent(in) :: charge
  real(p_double), intent(in) :: time
  integer, intent(in) :: np, tid

  ! n is the direction vector and R is the distance
  ! tdep time of deposition
  real(p_k_part) :: Emf(2,p_cache_size), n(3,p_cache_size), xyzpos(3)
  real(p_double) :: tdep(p_cache_size), tdepprev(p_cache_size), tprev, rsim_dt, temp
  real(p_double) :: tdep_c(p_cache_size), tdepprev_c(p_cache_size)
  real(p_k_part) :: betaprev_c(p_p_dim,p_cache_size), betadot_c(p_p_dim,p_cache_size)
  real(p_k_part) :: R(p_cache_size), rR, betadot(p_p_dim,p_cache_size)
  integer :: tit, it, titprev, c, i, l, npp
  integer :: iters(3), iter1, iter2, nf1
  integer(p_int64) :: iit, nf1_long
  real(p_k_part) :: div_charge(p_cache_size)

  real(p_k_part), dimension(p_cache_size) :: nmb1, nmb2, nmb3
# 888 "radio/os-radio-cartdet.f03"
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

      nmb1(i)=n(1,i)-betaprev_c(1,i)
      nmb2(i)=n(2,i)-betaprev_c(2,i)
      nmb3(i)=n(3,i)-betaprev_c(3,i)
# 984 "radio/os-radio-cartdet.f03"
      div_charge(i)=div_charge(i)/((1-betaprev_c(1,i)*n(1,i)-betaprev_c(2,i)*n(2,i)&
                                     -betaprev_c(3,i)*n(3,i))**3)
    enddo

    do i = 1, npp





      Emf(1,i)=div_charge(i)*(nmb2(i)*(betadot_c(1,i)*n(1,i)+betadot_c(3,i)*n(3,i)) &
             + betadot_c(2,i)*(-nmb1(i)*n(1,i)-nmb3(i)*n(3,i)))


      Emf(2,i)=div_charge(i)*(betadot_c(3,i)*(-nmb1(i)*n(1,i)-nmb2(i)*n(2,i)) &
             + nmb3(i)*(betadot_c(1,i)*n(1,i)+betadot_c(2,i)*n(2,i)))
# 1016 "radio/os-radio-cartdet.f03"
    enddo

    do i = 1, npp
      tit=int((tdep_c(i)-this%tpos(1))*this%d1m1)+1
      titprev=int((tdepprev_c(i)-this%tpos(1))*this%d1m1)+1
      !interpolation shenanigans. refer to report for better understanding
      do it=titprev,tit-1
        temp = (this%tpos(it+1)-tdepprev_c(i))
        do c = 1, 2
          this%Emf(c,iter1,iter2,it) = this%Emf(c,iter1,iter2,it) &
                                     + real(Emf(c,i)*temp,p_single)
        enddo
        tdepprev_c(i)=this%tpos(it+1)
      enddo
      temp = (tdep_c(i)-tdepprev_c(i))
      do c = 1, 2
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

end subroutine calc_comp_E2E3
!-----------------------------------------------------------------------------------------
# 433 "radio/os-radio-cartdet.f03" 2






# 1 "radio/os-radio-cartdet.f03" 1
! if is defined then read the template definition at the end of the file
# 854 "radio/os-radio-cartdet.f03"
!-----------------------------------------------------------------------------------------
subroutine calc_comp_E1E2E3(this,x,beta,time,xprev,betaprev,charge,np,tid)
  implicit none
  class(cart_Detector), intent(inout) :: this
  real(p_k_part), dimension(3,p_cache_size), intent(in) :: x, xprev
  real(p_k_part), dimension(p_p_dim,p_cache_size), intent(in) :: beta, betaprev
  real(p_k_part), dimension(p_cache_size), intent(in) :: charge
  real(p_double), intent(in) :: time
  integer, intent(in) :: np, tid

  ! n is the direction vector and R is the distance
  ! tdep time of deposition
  real(p_k_part) :: Emf(3,p_cache_size), n(3,p_cache_size), xyzpos(3)
  real(p_double) :: tdep(p_cache_size), tdepprev(p_cache_size), tprev, rsim_dt, temp
  real(p_double) :: tdep_c(p_cache_size), tdepprev_c(p_cache_size)
  real(p_k_part) :: betaprev_c(p_p_dim,p_cache_size), betadot_c(p_p_dim,p_cache_size)
  real(p_k_part) :: R(p_cache_size), rR, betadot(p_p_dim,p_cache_size)
  integer :: tit, it, titprev, c, i, l, npp
  integer :: iters(3), iter1, iter2, nf1
  integer(p_int64) :: iit, nf1_long
  real(p_k_part) :: div_charge(p_cache_size)

  real(p_k_part), dimension(p_cache_size) :: nmb1, nmb2, nmb3
# 888 "radio/os-radio-cartdet.f03"
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

      nmb1(i)=n(1,i)-betaprev_c(1,i)
      nmb2(i)=n(2,i)-betaprev_c(2,i)
      nmb3(i)=n(3,i)-betaprev_c(3,i)
# 984 "radio/os-radio-cartdet.f03"
      div_charge(i)=div_charge(i)/((1-betaprev_c(1,i)*n(1,i)-betaprev_c(2,i)*n(2,i)&
                                     -betaprev_c(3,i)*n(3,i))**3)
    enddo

    do i = 1, npp

      Emf(1,i)=div_charge(i)*(n(2,i)*(betadot_c(2,i)*nmb1(i)-betadot_c(1,i)*nmb2(i)) &
            +n(3,i)*(betadot_c(3,i)*nmb1(i)-betadot_c(1,i)*nmb3(i)))


      Emf(2,i)=div_charge(i)*(nmb2(i)*(betadot_c(1,i)*n(1,i)+betadot_c(3,i)*n(3,i)) &
             + betadot_c(2,i)*(-nmb1(i)*n(1,i)-nmb3(i)*n(3,i)))


      Emf(3,i)=div_charge(i)*(betadot_c(3,i)*(-nmb1(i)*n(1,i)-nmb2(i)*n(2,i)) &
             + nmb3(i)*(betadot_c(1,i)*n(1,i)+betadot_c(2,i)*n(2,i)))
# 1016 "radio/os-radio-cartdet.f03"
    enddo

    do i = 1, npp
      tit=int((tdep_c(i)-this%tpos(1))*this%d1m1)+1
      titprev=int((tdepprev_c(i)-this%tpos(1))*this%d1m1)+1
      !interpolation shenanigans. refer to report for better understanding
      do it=titprev,tit-1
        temp = (this%tpos(it+1)-tdepprev_c(i))
        do c = 1, 3
          this%Emf(c,iter1,iter2,it) = this%Emf(c,iter1,iter2,it) &
                                     + real(Emf(c,i)*temp,p_single)
        enddo
        tdepprev_c(i)=this%tpos(it+1)
      enddo
      temp = (tdep_c(i)-tdepprev_c(i))
      do c = 1, 3
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

end subroutine calc_comp_E1E2E3
!-----------------------------------------------------------------------------------------
# 440 "radio/os-radio-cartdet.f03" 2
# 849 "radio/os-radio-cartdet.f03"
end module class_CartDet
