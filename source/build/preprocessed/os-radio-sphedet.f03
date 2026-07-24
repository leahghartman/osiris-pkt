# 1 "radio/os-radio-sphedet.f03"
# 1 "<built-in>" 1
# 1 "<built-in>" 3
# 467 "<built-in>" 3
# 1 "<command line>" 1
# 1 "<built-in>" 2
# 1 "radio/os-radio-sphedet.f03" 2
! if __TEMPLATE__ is defined then read the template definition at the end of the file


module class_SpheDet

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
# 7 "radio/os-radio-sphedet.f03" 2
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
# 8 "radio/os-radio-sphedet.f03" 2
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
# 9 "radio/os-radio-sphedet.f03" 2

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
character(len=*), parameter :: p_sphe_det_rst_id = "sphe_det rst data - 0x0000"

type, public, extends(Detector) :: sphe_Detector
  real(p_k_part), pointer :: e4(:,:,:) => null(), e3(:,:,:) => null(), e2(:,:,:) => null()
  logical :: approx
  ! Varies based on approx and with field component
  procedure(fill_emf_loc_sphe_det), pointer :: fill_emf_loc => null()
contains
  procedure :: destruct => destruct_sphe_det
  procedure :: FillPos => FillPos_sphe_det
  procedure :: read_input => read_input_sphe_det
  procedure :: init => init_sphe_det
  procedure :: restart_write => restart_write_sphe_det
  procedure :: restart_read => restart_read_sphe_det
  procedure :: fill_emf => fill_emf_sphe_det
  procedure :: comp_select => comp_select_sphe_det
end type sphe_Detector

abstract interface
  subroutine fill_emf_loc_sphe_det(this,x,beta,time,xprev,betaprev,charge,np,tid)
    import sphe_Detector, p_k_part, p_cache_size, p_p_dim, p_double
    class(sphe_Detector), intent(inout):: this
    real(p_k_part), dimension(3,p_cache_size), intent(in) :: x, xprev
    real(p_k_part), dimension(p_p_dim,p_cache_size), intent(in) :: beta, betaprev
    real(p_k_part), dimension(p_cache_size), intent(in) :: charge
    real(p_double), intent(in) :: time
    integer, intent(in) :: np, tid
  end subroutine
end interface

contains

!-----------------------------------------------------------------------------------------
subroutine read_input_sphe_det(this,input_file,dt)

  implicit none

  real(p_double), intent(in) :: dt
  class( sphe_Detector ), intent(inout) :: this
  class( t_input_file ), intent(inout) :: input_file
  integer :: ndump_fac_det
  integer :: nf(4), ierr
  real(p_k_part) :: dmin(4), dmax(4)
  logical :: comps(6), apprx, approx

  logical :: has_b
  integer :: i

  namelist /nl_spherical/ apprx, approx, nf, dmin, dmax, comps, ndump_fac_det

  ndump_fac_det = 0
  nf = 0
  dmin = 0
  dmax = 0
  comps = .false.
  apprx = .false. ! kept for legacy reasons
  approx = .false.

  ! Get namelist text from input file
  call get_namelist( input_file, "nl_spherical", ierr )
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

  read (input_file%nml_text, nml = nl_spherical, iostat = ierr)
  if (ierr /= 0) then
    if(mpi_node()==0) then
      write(0,*) "Error reading detect parameters"
      write(0,*) "aborting..."
    endif
    stop
  endif

  this%approx=(apprx .or. approx)
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

  if (this%nf(4)>1 .and. this%approx) then
    if(mpi_node()==0) then
      write(0,*) "Error: only one cell is allowed in the r-direction when"
      write(0,*) "approx is set to true."
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

end subroutine read_input_sphe_det
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine init_sphe_det(this,name,no_co,restart,restart_handle)

  implicit none

  class( sphe_Detector ), intent(inout) :: this
  character(len = 64), intent(in) :: name
  class( t_node_conf ), intent(in) :: no_co
  logical, intent(in) :: restart
  type( t_restart_handle ), intent(in) :: restart_handle

  integer :: i, j, dp1, dp2, nt
  integer(p_int64) :: ncells, chunk, i0, i1
  integer :: sizes(4)
  real(p_k_part) :: lbounds(4), ubounds(4)
  character(len=256) :: fileLabel = '',aux2
  character(len=20), dimension(4) :: bas_ax_nms, bas_ax_lnms, bas_ax_uns ! name of members

  call alloc(this%dims,(/ this%ndims /),"radio/os-radio-sphedet.f03",180)
  call alloc(this%dmin_,(/ this%ndims /),"radio/os-radio-sphedet.f03",181)
  call alloc(this%dmax_,(/ this%ndims /),"radio/os-radio-sphedet.f03",182)

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
  bas_ax_nms =(/"\theta","\phi  ","r     ","t     "/)
  bas_ax_lnms =(/"\theta","\phi  ","r     ","t     "/)
  bas_ax_uns =(/"rad       ","rad       ","rad       ","1/\omega_p" /)

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

  !1,2,3,4=t,theta,phi,r
  this%n2 =this%nf(2) !theta
  this%n3 =this%nf(3) !phi
  this%n4 =this%nf(4) !r
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

  call alloc(this%e2,(/3,this%nf(2),this%nf(3)/),"radio/os-radio-sphedet.f03",245)
  call alloc(this%e3,(/3,this%nf(2),this%nf(3)/),"radio/os-radio-sphedet.f03",246)
  call alloc(this%e4,(/3,this%nf(2),this%nf(3)/),"radio/os-radio-sphedet.f03",247)
  call alloc(this%tpos,(/this%nf(1)+1/),"radio/os-radio-sphedet.f03",248)

  if (this%ncmp>=1) then
    dp1 = this%dim_pos(1)
    dp2 = this%dim_pos(2)

    call alloc(this%Emf, (/ this%ncmp, this%nf(dp1), this%nf(dp2), this%nf(1)+1 /),"radio/os-radio-sphedet.f03",254)
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

  call alloc(this%iter_bnds, (/ 2, nt /),"radio/os-radio-sphedet.f03",278)
  do i = 0, nt-1
    chunk = ( ncells + nt - 1 ) / nt
    i0 = i * chunk + 1
    i1 = min( (i+1) * chunk, ncells )
    this%iter_bnds(p_lower,i+1) = i0
    this%iter_bnds(p_upper,i+1) = i1
  enddo

end subroutine init_sphe_det
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine FillPos_sphe_det(this)

  implicit none

  class(sphe_Detector), INTENT(INOUT):: this
  integer :: j,k
  !efficient array
  !do i
  ! do j
  ! x(j,i)

  ! th->2
  ! phi->3
  ! r->1

  !$omp parallel
  !$omp do private(k,j)
  do j=0,this%nf(3)-1
    do k=0,this%nf(2)-1
      this%e4(:,k+1,j+1)=(/ sin(this%dmin(2)+k*this%d2)*cos(this%dmin(3)+j*this%d3),&
                            sin(this%dmin(2)+k*this%d2)*sin(this%dmin(3)+j*this%d3),&
                            cos(this%dmin(2)+k*this%d2) /)
      this%e2(:,k+1,j+1)=(/ cos(this%dmin(2)+k*this%d2)*cos(this%dmin(3)+j*this%d3),&
                            cos(this%dmin(2)+k*this%d2)*sin(this%dmin(3)+j*this%d3),&
                           -sin(this%dmin(2)+k*this%d2) /)
      this%e3(:,k+1,j+1)=(/-sin(this%dmin(3)+j*this%d3),&
                            cos(this%dmin(3)+j*this%d3), real(0.0,p_k_part) /)
    enddo
  enddo
  !$omp enddo
  !$omp end parallel
  ! t->1
  if(this%approx) then
    this%tpos=(/ (real(this%dmin(1),p_double)+k*this%d1,k=-1,this%nf(1)-1) /)
  else
    this%tpos=(/ (real(this%dmin(4),p_double)+real(this%dmin(1),p_double)+k*this%d1,k=-1,&
                  this%nf(1)-1) /)
  endif

end subroutine FillPos_sphe_det
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine destruct_sphe_det(this)
  implicit none
  class(sphe_Detector), INTENT(INOUT):: this

  call freemem(this%e4,"radio/os-radio-sphedet.f03",338)
  call freemem(this%e3,"radio/os-radio-sphedet.f03",339)
  call freemem(this%e2,"radio/os-radio-sphedet.f03",340)
  call freemem(this%tpos,"radio/os-radio-sphedet.f03",341)
  call freemem(this%Emf,"radio/os-radio-sphedet.f03",342)
  call freemem(this%dims,"radio/os-radio-sphedet.f03",343)
  call freemem(this%dmin_,"radio/os-radio-sphedet.f03",344)
  call freemem(this%dmax_,"radio/os-radio-sphedet.f03",345)
  call freemem(this%iter_bnds,"radio/os-radio-sphedet.f03",346)

end subroutine destruct_sphe_det
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine restart_write_sphe_det( this, restart_handle )

  implicit none

  class(sphe_Detector), intent(in) :: this
  type( t_restart_handle ), intent(inout) :: restart_handle

  character(len=*), parameter :: err_msg='error writing restart data for sphe_det object.'
  integer :: ierr

  ! Write data
  call restart_io_write("p_sphe_det_rst_id", p_sphe_det_rst_id, restart_handle, ierr)
  call check_error(ierr,err_msg,p_err_rstwrt,"radio/os-radio-sphedet.f03",364)

  call restart_io_write("this%Emf", this%Emf, restart_handle, ierr)
  call check_error(ierr,err_msg,p_err_rstwrt,"radio/os-radio-sphedet.f03",367)

end subroutine restart_write_sphe_det
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine restart_read_sphe_det( this, restart_handle )

  implicit none

  class(sphe_Detector), intent(inout) :: this
  type( t_restart_handle ), intent(in) :: restart_handle

  character(len=*), parameter :: err_msg='error reading restart data for sphe_det object.'
  character(len=len(p_sphe_det_rst_id)) :: rst_id
  integer :: ierr

  ! Read data
  call restart_io_read(rst_id, restart_handle, ierr)
  call check_error(ierr,err_msg,p_err_rstrd,"radio/os-radio-sphedet.f03",386)

  ! check if restart file is compatible
  if ( rst_id /= p_sphe_det_rst_id ) then
    write(err_buf__,*) 'Corrupted restart file, or restart file ';call err__("radio/os-radio-sphedet.f03",390)
    write(err_buf__,*) 'from incompatible binary (sphe_det)';call err__("radio/os-radio-sphedet.f03",391)
    call abort_program(p_err_rstrd)
  endif

  call restart_io_read(this%Emf, restart_handle, ierr)
  call check_error(ierr,err_msg,p_err_rstrd,"radio/os-radio-sphedet.f03",396)

end subroutine restart_read_sphe_det
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine fill_emf_sphe_det(this,x,beta,time,xprev,betaprev,charge,np,tid)
  implicit none
  class(sphe_Detector), intent(inout) :: this
  real(p_k_part), dimension(3,p_cache_size), intent(in) :: x, xprev
  real(p_k_part), dimension(p_p_dim,p_cache_size), intent(in) :: beta, betaprev
  real(p_k_part), dimension(p_cache_size), intent(in) :: charge
  real(p_double), intent(in) :: time
  integer, intent(in) :: np, tid

  ! To allow for using approx as true or false and for component selection
  call this%fill_emf_loc(x,beta,time,xprev,betaprev,charge,np,tid)

end subroutine fill_emf_sphe_det
!-----------------------------------------------------------------------------------------

! comp_select defined here
# 1 "radio/os-sphe-sel.f03" 1
!-----------------------------------------------------------------------------------------
subroutine calc_dflt(this,x,beta,time,xprev,betaprev,charge,np,tid)
  implicit none
  class(sphe_Detector), intent(inout) :: this
  real(p_k_part), dimension(3,p_cache_size), intent(in) :: x, xprev
  real(p_k_part), dimension(p_p_dim,p_cache_size), intent(in) :: beta, betaprev
  real(p_k_part), dimension(p_cache_size), intent(in) :: charge
  real(p_double), intent(in) :: time
  integer, intent(in) :: np, tid

  ! Do nothing

end subroutine calc_dflt
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine comp_select_sphe_det( this )

  implicit none

  class(sphe_Detector), intent(inout) :: this

  integer :: sele

  sele = 0
  if (this%cmpts(6)) sele = sele + 32
  if (this%cmpts(5)) sele = sele + 16
  if (this%cmpts(4)) sele = sele + 8
  if (this%cmpts(3)) sele = sele + 4
  if (this%cmpts(2)) sele = sele + 2
  if (this%cmpts(1)) sele = sele + 1

  if (this%approx) then

    select case(sele)
    case(1)
      this % fill_emf_loc => calc_comp_approx_E1
    case(2)
      this % fill_emf_loc => calc_comp_approx_E2
    case(3)
      this % fill_emf_loc => calc_comp_approx_E1E2
    case(4)
      this % fill_emf_loc => calc_comp_approx_E3
    case(5)
      this % fill_emf_loc => calc_comp_approx_E1E3
    case(6)
      this % fill_emf_loc => calc_comp_approx_E2E3
    case(7)
      this % fill_emf_loc => calc_comp_approx_E1E2E3
# 164 "radio/os-sphe-sel.f03"
    case default
      this % fill_emf_loc => calc_dflt
    end select

  else

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
# 299 "radio/os-sphe-sel.f03"
    case default
      this % fill_emf_loc => calc_dflt
    end select

  endif

end subroutine comp_select_sphe_det
!-----------------------------------------------------------------------------------------
# 419 "radio/os-radio-sphedet.f03" 2

! Now define component functions





# 1 "radio/os-radio-sphedet.f03" 1
! if is defined then read the template definition at the end of the file
# 876 "radio/os-radio-sphedet.f03"
!-----------------------------------------------------------------------------------------
subroutine calc_comp_E1(this,x,beta,time,xprev,betaprev,charge,np,tid)
  implicit none
  class(sphe_Detector), intent(inout) :: this
  real(p_k_part), dimension(3,p_cache_size), intent(in) :: x, xprev
  real(p_k_part), dimension(p_p_dim,p_cache_size), intent(in) :: beta, betaprev
  real(p_k_part), dimension(p_cache_size), intent(in) :: charge
  real(p_double), intent(in) :: time
  integer, intent(in) :: np, tid

  ! n is the direction vector and R is the distance
  ! tdep time of deposition
  real(p_k_part) :: Emf(1,p_cache_size), n(3,p_cache_size), e4_radi(3)
  real(p_double) :: tdep(p_cache_size), tdepprev(p_cache_size), tprev, rsim_dt, temp
  real(p_double) :: tdep_c(p_cache_size), tdepprev_c(p_cache_size)
  real(p_k_part) :: betaprev_c(p_p_dim,p_cache_size), betadot_c(p_p_dim,p_cache_size)
  real(p_k_part) :: R(p_cache_size), rR, betadot(p_p_dim,p_cache_size)
  integer :: tit, it, titprev, c, i, l, npp
  integer :: iters(3), iter1, iter2, nf1
  integer(p_int64) :: iit, nf1_long
  real(p_k_part) :: div_charge(p_cache_size), radi

  real(p_k_part), dimension(p_cache_size) :: nmb1, nmb2, nmb3
# 910 "radio/os-radio-sphedet.f03"
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

    !Compute the non-normalized direction vector x-pos
    radi=this%dmin(4)+(iters(3)-1)*this%d4
    e4_radi = this%e4(:,iters(1),iters(2)) * radi
    do i = 1, np
      n(1,i) = x(1,i) - e4_radi(1)
      n(2,i) = x(2,i) - e4_radi(2)
      n(3,i) = x(3,i) - e4_radi(3)
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
      n(1,i) = xprev(1,i) - e4_radi(1)
      n(2,i) = xprev(2,i) - e4_radi(2)
      n(3,i) = xprev(3,i) - e4_radi(3)
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
# 1004 "radio/os-radio-sphedet.f03"
      div_charge(i)=div_charge(i)/((1-betaprev_c(1,i)*n(1,i)-betaprev_c(2,i)*n(2,i)&
                                     -betaprev_c(3,i)*n(3,i))**3)
    enddo

    do i = 1, npp

      Emf(1,i)=div_charge(i)*(n(2,i)*(betadot_c(2,i)*nmb1(i)-betadot_c(1,i)*nmb2(i)) &
            +n(3,i)*(betadot_c(3,i)*nmb1(i)-betadot_c(1,i)*nmb3(i)))
# 1036 "radio/os-radio-sphedet.f03"
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

!-----------------------------------------------------------------------------------------
subroutine calc_comp_approx_E1(this,x,beta,time,xprev,betaprev,charge,np,tid)
  implicit none
  class(sphe_Detector), intent(inout) :: this
  real(p_k_part), dimension(3,p_cache_size), intent(in) :: x, xprev
  real(p_k_part), dimension(p_p_dim,p_cache_size), intent(in) :: beta, betaprev
  real(p_k_part), dimension(p_cache_size), intent(in) :: charge
  real(p_double), intent(in) :: time
  integer, intent(in) :: np, tid

  ! n is the direction vector and R is the distance
  ! tdep time of deposition
  real(p_k_part) :: Emf(1,p_cache_size), n(3,p_cache_size)
  real(p_double) :: tdep(p_cache_size), tdepprev(p_cache_size), tprev, rsim_dt, temp
  real(p_double) :: tdep_c(p_cache_size), tdepprev_c(p_cache_size)
  real(p_k_part) :: betaprev_c(p_p_dim,p_cache_size), betadot_c(p_p_dim,p_cache_size)
  real(p_k_part) :: R, betadot(p_p_dim,p_cache_size)
  integer :: tit, it, titprev, c, i, l, npp
  integer :: iters(3), iter1, iter2, nf1
  integer(p_int64) :: iit, nf1_long
  real(p_k_part) :: div_charge(p_cache_size)

  real(p_k_part), dimension(p_cache_size) :: nmb1, nmb2, nmb3
# 1102 "radio/os-radio-sphedet.f03"
  rsim_dt = 1.0_p_double / this%sim_dt
  do i = 1, np
    betadot(1,i) = (beta(1,i)-betaprev(1,i)) * rsim_dt
    betadot(2,i) = (beta(2,i)-betaprev(2,i)) * rsim_dt
    betadot(3,i) = (beta(3,i)-betaprev(3,i)) * rsim_dt
  enddo

  tprev=time-this%sim_dt
  R=this%dmin(4)
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

    do i = 1, np
      n(1,i) = this%e4(1,iters(1),iters(2))
      n(2,i) = this%e4(2,iters(1),iters(2))
      n(3,i) = this%e4(3,iters(1),iters(2))
    enddo
    do i = 1, np
      tdep(i) = time - n(1,i)*x(1,i) - n(2,i)*x(2,i) - n(3,i)*x(3,i)
      tdepprev(i) = tprev - n(1,i)*xprev(1,i) - n(2,i)*xprev(2,i) - n(3,i)*xprev(3,i)
    enddo

    l = 1
    do i = 1, np
      if ((tdepprev(i) >= this%tpos(1)) .and. (tdep(i) <= this%tpos(this%n1+1))) then

        betaprev_c(1,l) = betaprev(1,i)
        betaprev_c(2,l) = betaprev(2,i)
        betaprev_c(3,l) = betaprev(3,i)
        betadot_c(1,l) = betadot(1,i)
        betadot_c(2,l) = betadot(2,i)
        betadot_c(3,l) = betadot(3,i)
        div_charge(l) = charge(i) ! will be modified again later
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
# 1171 "radio/os-radio-sphedet.f03"
      div_charge(i)=div_charge(i)/(R*(1-betaprev_c(1,i)*n(1,i)-betaprev_c(2,i)*n(2,i)&
                                       -betaprev_c(3,i)*n(3,i))**3)
    enddo

    do i = 1, npp

      Emf(1,i)=div_charge(i)*(n(2,i)*(betadot_c(2,i)*nmb1(i)-betadot_c(1,i)*nmb2(i)) &
            +n(3,i)*(betadot_c(3,i)*nmb1(i)-betadot_c(1,i)*nmb3(i)))
# 1203 "radio/os-radio-sphedet.f03"
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

end subroutine calc_comp_approx_E1
!-----------------------------------------------------------------------------------------
# 427 "radio/os-radio-sphedet.f03" 2




# 1 "radio/os-radio-sphedet.f03" 1
! if is defined then read the template definition at the end of the file
# 876 "radio/os-radio-sphedet.f03"
!-----------------------------------------------------------------------------------------
subroutine calc_comp_E2(this,x,beta,time,xprev,betaprev,charge,np,tid)
  implicit none
  class(sphe_Detector), intent(inout) :: this
  real(p_k_part), dimension(3,p_cache_size), intent(in) :: x, xprev
  real(p_k_part), dimension(p_p_dim,p_cache_size), intent(in) :: beta, betaprev
  real(p_k_part), dimension(p_cache_size), intent(in) :: charge
  real(p_double), intent(in) :: time
  integer, intent(in) :: np, tid

  ! n is the direction vector and R is the distance
  ! tdep time of deposition
  real(p_k_part) :: Emf(1,p_cache_size), n(3,p_cache_size), e4_radi(3)
  real(p_double) :: tdep(p_cache_size), tdepprev(p_cache_size), tprev, rsim_dt, temp
  real(p_double) :: tdep_c(p_cache_size), tdepprev_c(p_cache_size)
  real(p_k_part) :: betaprev_c(p_p_dim,p_cache_size), betadot_c(p_p_dim,p_cache_size)
  real(p_k_part) :: R(p_cache_size), rR, betadot(p_p_dim,p_cache_size)
  integer :: tit, it, titprev, c, i, l, npp
  integer :: iters(3), iter1, iter2, nf1
  integer(p_int64) :: iit, nf1_long
  real(p_k_part) :: div_charge(p_cache_size), radi

  real(p_k_part), dimension(p_cache_size) :: nmb1, nmb2, nmb3
# 910 "radio/os-radio-sphedet.f03"
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

    !Compute the non-normalized direction vector x-pos
    radi=this%dmin(4)+(iters(3)-1)*this%d4
    e4_radi = this%e4(:,iters(1),iters(2)) * radi
    do i = 1, np
      n(1,i) = x(1,i) - e4_radi(1)
      n(2,i) = x(2,i) - e4_radi(2)
      n(3,i) = x(3,i) - e4_radi(3)
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
      n(1,i) = xprev(1,i) - e4_radi(1)
      n(2,i) = xprev(2,i) - e4_radi(2)
      n(3,i) = xprev(3,i) - e4_radi(3)
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
# 1004 "radio/os-radio-sphedet.f03"
      div_charge(i)=div_charge(i)/((1-betaprev_c(1,i)*n(1,i)-betaprev_c(2,i)*n(2,i)&
                                     -betaprev_c(3,i)*n(3,i))**3)
    enddo

    do i = 1, npp





      Emf(1,i)=div_charge(i)*(nmb2(i)*(betadot_c(1,i)*n(1,i)+betadot_c(3,i)*n(3,i)) &
             + betadot_c(2,i)*(-nmb1(i)*n(1,i)-nmb3(i)*n(3,i)))
# 1036 "radio/os-radio-sphedet.f03"
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

!-----------------------------------------------------------------------------------------
subroutine calc_comp_approx_E2(this,x,beta,time,xprev,betaprev,charge,np,tid)
  implicit none
  class(sphe_Detector), intent(inout) :: this
  real(p_k_part), dimension(3,p_cache_size), intent(in) :: x, xprev
  real(p_k_part), dimension(p_p_dim,p_cache_size), intent(in) :: beta, betaprev
  real(p_k_part), dimension(p_cache_size), intent(in) :: charge
  real(p_double), intent(in) :: time
  integer, intent(in) :: np, tid

  ! n is the direction vector and R is the distance
  ! tdep time of deposition
  real(p_k_part) :: Emf(1,p_cache_size), n(3,p_cache_size)
  real(p_double) :: tdep(p_cache_size), tdepprev(p_cache_size), tprev, rsim_dt, temp
  real(p_double) :: tdep_c(p_cache_size), tdepprev_c(p_cache_size)
  real(p_k_part) :: betaprev_c(p_p_dim,p_cache_size), betadot_c(p_p_dim,p_cache_size)
  real(p_k_part) :: R, betadot(p_p_dim,p_cache_size)
  integer :: tit, it, titprev, c, i, l, npp
  integer :: iters(3), iter1, iter2, nf1
  integer(p_int64) :: iit, nf1_long
  real(p_k_part) :: div_charge(p_cache_size)

  real(p_k_part), dimension(p_cache_size) :: nmb1, nmb2, nmb3
# 1102 "radio/os-radio-sphedet.f03"
  rsim_dt = 1.0_p_double / this%sim_dt
  do i = 1, np
    betadot(1,i) = (beta(1,i)-betaprev(1,i)) * rsim_dt
    betadot(2,i) = (beta(2,i)-betaprev(2,i)) * rsim_dt
    betadot(3,i) = (beta(3,i)-betaprev(3,i)) * rsim_dt
  enddo

  tprev=time-this%sim_dt
  R=this%dmin(4)
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

    do i = 1, np
      n(1,i) = this%e4(1,iters(1),iters(2))
      n(2,i) = this%e4(2,iters(1),iters(2))
      n(3,i) = this%e4(3,iters(1),iters(2))
    enddo
    do i = 1, np
      tdep(i) = time - n(1,i)*x(1,i) - n(2,i)*x(2,i) - n(3,i)*x(3,i)
      tdepprev(i) = tprev - n(1,i)*xprev(1,i) - n(2,i)*xprev(2,i) - n(3,i)*xprev(3,i)
    enddo

    l = 1
    do i = 1, np
      if ((tdepprev(i) >= this%tpos(1)) .and. (tdep(i) <= this%tpos(this%n1+1))) then

        betaprev_c(1,l) = betaprev(1,i)
        betaprev_c(2,l) = betaprev(2,i)
        betaprev_c(3,l) = betaprev(3,i)
        betadot_c(1,l) = betadot(1,i)
        betadot_c(2,l) = betadot(2,i)
        betadot_c(3,l) = betadot(3,i)
        div_charge(l) = charge(i) ! will be modified again later
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
# 1171 "radio/os-radio-sphedet.f03"
      div_charge(i)=div_charge(i)/(R*(1-betaprev_c(1,i)*n(1,i)-betaprev_c(2,i)*n(2,i)&
                                       -betaprev_c(3,i)*n(3,i))**3)
    enddo

    do i = 1, npp





      Emf(1,i)=div_charge(i)*(nmb2(i)*(betadot_c(1,i)*n(1,i)+betadot_c(3,i)*n(3,i)) &
             + betadot_c(2,i)*(-nmb1(i)*n(1,i)-nmb3(i)*n(3,i)))
# 1203 "radio/os-radio-sphedet.f03"
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

end subroutine calc_comp_approx_E2
!-----------------------------------------------------------------------------------------
# 432 "radio/os-radio-sphedet.f03" 2




# 1 "radio/os-radio-sphedet.f03" 1
! if is defined then read the template definition at the end of the file
# 876 "radio/os-radio-sphedet.f03"
!-----------------------------------------------------------------------------------------
subroutine calc_comp_E3(this,x,beta,time,xprev,betaprev,charge,np,tid)
  implicit none
  class(sphe_Detector), intent(inout) :: this
  real(p_k_part), dimension(3,p_cache_size), intent(in) :: x, xprev
  real(p_k_part), dimension(p_p_dim,p_cache_size), intent(in) :: beta, betaprev
  real(p_k_part), dimension(p_cache_size), intent(in) :: charge
  real(p_double), intent(in) :: time
  integer, intent(in) :: np, tid

  ! n is the direction vector and R is the distance
  ! tdep time of deposition
  real(p_k_part) :: Emf(1,p_cache_size), n(3,p_cache_size), e4_radi(3)
  real(p_double) :: tdep(p_cache_size), tdepprev(p_cache_size), tprev, rsim_dt, temp
  real(p_double) :: tdep_c(p_cache_size), tdepprev_c(p_cache_size)
  real(p_k_part) :: betaprev_c(p_p_dim,p_cache_size), betadot_c(p_p_dim,p_cache_size)
  real(p_k_part) :: R(p_cache_size), rR, betadot(p_p_dim,p_cache_size)
  integer :: tit, it, titprev, c, i, l, npp
  integer :: iters(3), iter1, iter2, nf1
  integer(p_int64) :: iit, nf1_long
  real(p_k_part) :: div_charge(p_cache_size), radi

  real(p_k_part), dimension(p_cache_size) :: nmb1, nmb2, nmb3
# 910 "radio/os-radio-sphedet.f03"
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

    !Compute the non-normalized direction vector x-pos
    radi=this%dmin(4)+(iters(3)-1)*this%d4
    e4_radi = this%e4(:,iters(1),iters(2)) * radi
    do i = 1, np
      n(1,i) = x(1,i) - e4_radi(1)
      n(2,i) = x(2,i) - e4_radi(2)
      n(3,i) = x(3,i) - e4_radi(3)
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
      n(1,i) = xprev(1,i) - e4_radi(1)
      n(2,i) = xprev(2,i) - e4_radi(2)
      n(3,i) = xprev(3,i) - e4_radi(3)
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
# 1004 "radio/os-radio-sphedet.f03"
      div_charge(i)=div_charge(i)/((1-betaprev_c(1,i)*n(1,i)-betaprev_c(2,i)*n(2,i)&
                                     -betaprev_c(3,i)*n(3,i))**3)
    enddo

    do i = 1, npp
# 1018 "radio/os-radio-sphedet.f03"
      Emf(1,i)=div_charge(i)*(betadot_c(3,i)*(-nmb1(i)*n(1,i)-nmb2(i)*n(2,i)) &
             + nmb3(i)*(betadot_c(1,i)*n(1,i)+betadot_c(2,i)*n(2,i)))
# 1036 "radio/os-radio-sphedet.f03"
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

!-----------------------------------------------------------------------------------------
subroutine calc_comp_approx_E3(this,x,beta,time,xprev,betaprev,charge,np,tid)
  implicit none
  class(sphe_Detector), intent(inout) :: this
  real(p_k_part), dimension(3,p_cache_size), intent(in) :: x, xprev
  real(p_k_part), dimension(p_p_dim,p_cache_size), intent(in) :: beta, betaprev
  real(p_k_part), dimension(p_cache_size), intent(in) :: charge
  real(p_double), intent(in) :: time
  integer, intent(in) :: np, tid

  ! n is the direction vector and R is the distance
  ! tdep time of deposition
  real(p_k_part) :: Emf(1,p_cache_size), n(3,p_cache_size)
  real(p_double) :: tdep(p_cache_size), tdepprev(p_cache_size), tprev, rsim_dt, temp
  real(p_double) :: tdep_c(p_cache_size), tdepprev_c(p_cache_size)
  real(p_k_part) :: betaprev_c(p_p_dim,p_cache_size), betadot_c(p_p_dim,p_cache_size)
  real(p_k_part) :: R, betadot(p_p_dim,p_cache_size)
  integer :: tit, it, titprev, c, i, l, npp
  integer :: iters(3), iter1, iter2, nf1
  integer(p_int64) :: iit, nf1_long
  real(p_k_part) :: div_charge(p_cache_size)

  real(p_k_part), dimension(p_cache_size) :: nmb1, nmb2, nmb3
# 1102 "radio/os-radio-sphedet.f03"
  rsim_dt = 1.0_p_double / this%sim_dt
  do i = 1, np
    betadot(1,i) = (beta(1,i)-betaprev(1,i)) * rsim_dt
    betadot(2,i) = (beta(2,i)-betaprev(2,i)) * rsim_dt
    betadot(3,i) = (beta(3,i)-betaprev(3,i)) * rsim_dt
  enddo

  tprev=time-this%sim_dt
  R=this%dmin(4)
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

    do i = 1, np
      n(1,i) = this%e4(1,iters(1),iters(2))
      n(2,i) = this%e4(2,iters(1),iters(2))
      n(3,i) = this%e4(3,iters(1),iters(2))
    enddo
    do i = 1, np
      tdep(i) = time - n(1,i)*x(1,i) - n(2,i)*x(2,i) - n(3,i)*x(3,i)
      tdepprev(i) = tprev - n(1,i)*xprev(1,i) - n(2,i)*xprev(2,i) - n(3,i)*xprev(3,i)
    enddo

    l = 1
    do i = 1, np
      if ((tdepprev(i) >= this%tpos(1)) .and. (tdep(i) <= this%tpos(this%n1+1))) then

        betaprev_c(1,l) = betaprev(1,i)
        betaprev_c(2,l) = betaprev(2,i)
        betaprev_c(3,l) = betaprev(3,i)
        betadot_c(1,l) = betadot(1,i)
        betadot_c(2,l) = betadot(2,i)
        betadot_c(3,l) = betadot(3,i)
        div_charge(l) = charge(i) ! will be modified again later
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
# 1171 "radio/os-radio-sphedet.f03"
      div_charge(i)=div_charge(i)/(R*(1-betaprev_c(1,i)*n(1,i)-betaprev_c(2,i)*n(2,i)&
                                       -betaprev_c(3,i)*n(3,i))**3)
    enddo

    do i = 1, npp
# 1185 "radio/os-radio-sphedet.f03"
      Emf(1,i)=div_charge(i)*(betadot_c(3,i)*(-nmb1(i)*n(1,i)-nmb2(i)*n(2,i)) &
             + nmb3(i)*(betadot_c(1,i)*n(1,i)+betadot_c(2,i)*n(2,i)))
# 1203 "radio/os-radio-sphedet.f03"
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

end subroutine calc_comp_approx_E3
!-----------------------------------------------------------------------------------------
# 437 "radio/os-radio-sphedet.f03" 2





# 1 "radio/os-radio-sphedet.f03" 1
! if is defined then read the template definition at the end of the file
# 876 "radio/os-radio-sphedet.f03"
!-----------------------------------------------------------------------------------------
subroutine calc_comp_E1E2(this,x,beta,time,xprev,betaprev,charge,np,tid)
  implicit none
  class(sphe_Detector), intent(inout) :: this
  real(p_k_part), dimension(3,p_cache_size), intent(in) :: x, xprev
  real(p_k_part), dimension(p_p_dim,p_cache_size), intent(in) :: beta, betaprev
  real(p_k_part), dimension(p_cache_size), intent(in) :: charge
  real(p_double), intent(in) :: time
  integer, intent(in) :: np, tid

  ! n is the direction vector and R is the distance
  ! tdep time of deposition
  real(p_k_part) :: Emf(2,p_cache_size), n(3,p_cache_size), e4_radi(3)
  real(p_double) :: tdep(p_cache_size), tdepprev(p_cache_size), tprev, rsim_dt, temp
  real(p_double) :: tdep_c(p_cache_size), tdepprev_c(p_cache_size)
  real(p_k_part) :: betaprev_c(p_p_dim,p_cache_size), betadot_c(p_p_dim,p_cache_size)
  real(p_k_part) :: R(p_cache_size), rR, betadot(p_p_dim,p_cache_size)
  integer :: tit, it, titprev, c, i, l, npp
  integer :: iters(3), iter1, iter2, nf1
  integer(p_int64) :: iit, nf1_long
  real(p_k_part) :: div_charge(p_cache_size), radi

  real(p_k_part), dimension(p_cache_size) :: nmb1, nmb2, nmb3
# 910 "radio/os-radio-sphedet.f03"
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

    !Compute the non-normalized direction vector x-pos
    radi=this%dmin(4)+(iters(3)-1)*this%d4
    e4_radi = this%e4(:,iters(1),iters(2)) * radi
    do i = 1, np
      n(1,i) = x(1,i) - e4_radi(1)
      n(2,i) = x(2,i) - e4_radi(2)
      n(3,i) = x(3,i) - e4_radi(3)
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
      n(1,i) = xprev(1,i) - e4_radi(1)
      n(2,i) = xprev(2,i) - e4_radi(2)
      n(3,i) = xprev(3,i) - e4_radi(3)
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
# 1004 "radio/os-radio-sphedet.f03"
      div_charge(i)=div_charge(i)/((1-betaprev_c(1,i)*n(1,i)-betaprev_c(2,i)*n(2,i)&
                                     -betaprev_c(3,i)*n(3,i))**3)
    enddo

    do i = 1, npp

      Emf(1,i)=div_charge(i)*(n(2,i)*(betadot_c(2,i)*nmb1(i)-betadot_c(1,i)*nmb2(i)) &
            +n(3,i)*(betadot_c(3,i)*nmb1(i)-betadot_c(1,i)*nmb3(i)))


      Emf(2,i)=div_charge(i)*(nmb2(i)*(betadot_c(1,i)*n(1,i)+betadot_c(3,i)*n(3,i)) &
             + betadot_c(2,i)*(-nmb1(i)*n(1,i)-nmb3(i)*n(3,i)))
# 1036 "radio/os-radio-sphedet.f03"
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

!-----------------------------------------------------------------------------------------
subroutine calc_comp_approx_E1E2(this,x,beta,time,xprev,betaprev,charge,np,tid)
  implicit none
  class(sphe_Detector), intent(inout) :: this
  real(p_k_part), dimension(3,p_cache_size), intent(in) :: x, xprev
  real(p_k_part), dimension(p_p_dim,p_cache_size), intent(in) :: beta, betaprev
  real(p_k_part), dimension(p_cache_size), intent(in) :: charge
  real(p_double), intent(in) :: time
  integer, intent(in) :: np, tid

  ! n is the direction vector and R is the distance
  ! tdep time of deposition
  real(p_k_part) :: Emf(2,p_cache_size), n(3,p_cache_size)
  real(p_double) :: tdep(p_cache_size), tdepprev(p_cache_size), tprev, rsim_dt, temp
  real(p_double) :: tdep_c(p_cache_size), tdepprev_c(p_cache_size)
  real(p_k_part) :: betaprev_c(p_p_dim,p_cache_size), betadot_c(p_p_dim,p_cache_size)
  real(p_k_part) :: R, betadot(p_p_dim,p_cache_size)
  integer :: tit, it, titprev, c, i, l, npp
  integer :: iters(3), iter1, iter2, nf1
  integer(p_int64) :: iit, nf1_long
  real(p_k_part) :: div_charge(p_cache_size)

  real(p_k_part), dimension(p_cache_size) :: nmb1, nmb2, nmb3
# 1102 "radio/os-radio-sphedet.f03"
  rsim_dt = 1.0_p_double / this%sim_dt
  do i = 1, np
    betadot(1,i) = (beta(1,i)-betaprev(1,i)) * rsim_dt
    betadot(2,i) = (beta(2,i)-betaprev(2,i)) * rsim_dt
    betadot(3,i) = (beta(3,i)-betaprev(3,i)) * rsim_dt
  enddo

  tprev=time-this%sim_dt
  R=this%dmin(4)
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

    do i = 1, np
      n(1,i) = this%e4(1,iters(1),iters(2))
      n(2,i) = this%e4(2,iters(1),iters(2))
      n(3,i) = this%e4(3,iters(1),iters(2))
    enddo
    do i = 1, np
      tdep(i) = time - n(1,i)*x(1,i) - n(2,i)*x(2,i) - n(3,i)*x(3,i)
      tdepprev(i) = tprev - n(1,i)*xprev(1,i) - n(2,i)*xprev(2,i) - n(3,i)*xprev(3,i)
    enddo

    l = 1
    do i = 1, np
      if ((tdepprev(i) >= this%tpos(1)) .and. (tdep(i) <= this%tpos(this%n1+1))) then

        betaprev_c(1,l) = betaprev(1,i)
        betaprev_c(2,l) = betaprev(2,i)
        betaprev_c(3,l) = betaprev(3,i)
        betadot_c(1,l) = betadot(1,i)
        betadot_c(2,l) = betadot(2,i)
        betadot_c(3,l) = betadot(3,i)
        div_charge(l) = charge(i) ! will be modified again later
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
# 1171 "radio/os-radio-sphedet.f03"
      div_charge(i)=div_charge(i)/(R*(1-betaprev_c(1,i)*n(1,i)-betaprev_c(2,i)*n(2,i)&
                                       -betaprev_c(3,i)*n(3,i))**3)
    enddo

    do i = 1, npp

      Emf(1,i)=div_charge(i)*(n(2,i)*(betadot_c(2,i)*nmb1(i)-betadot_c(1,i)*nmb2(i)) &
            +n(3,i)*(betadot_c(3,i)*nmb1(i)-betadot_c(1,i)*nmb3(i)))


      Emf(2,i)=div_charge(i)*(nmb2(i)*(betadot_c(1,i)*n(1,i)+betadot_c(3,i)*n(3,i)) &
             + betadot_c(2,i)*(-nmb1(i)*n(1,i)-nmb3(i)*n(3,i)))
# 1203 "radio/os-radio-sphedet.f03"
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

end subroutine calc_comp_approx_E1E2
!-----------------------------------------------------------------------------------------
# 443 "radio/os-radio-sphedet.f03" 2





# 1 "radio/os-radio-sphedet.f03" 1
! if is defined then read the template definition at the end of the file
# 876 "radio/os-radio-sphedet.f03"
!-----------------------------------------------------------------------------------------
subroutine calc_comp_E1E3(this,x,beta,time,xprev,betaprev,charge,np,tid)
  implicit none
  class(sphe_Detector), intent(inout) :: this
  real(p_k_part), dimension(3,p_cache_size), intent(in) :: x, xprev
  real(p_k_part), dimension(p_p_dim,p_cache_size), intent(in) :: beta, betaprev
  real(p_k_part), dimension(p_cache_size), intent(in) :: charge
  real(p_double), intent(in) :: time
  integer, intent(in) :: np, tid

  ! n is the direction vector and R is the distance
  ! tdep time of deposition
  real(p_k_part) :: Emf(2,p_cache_size), n(3,p_cache_size), e4_radi(3)
  real(p_double) :: tdep(p_cache_size), tdepprev(p_cache_size), tprev, rsim_dt, temp
  real(p_double) :: tdep_c(p_cache_size), tdepprev_c(p_cache_size)
  real(p_k_part) :: betaprev_c(p_p_dim,p_cache_size), betadot_c(p_p_dim,p_cache_size)
  real(p_k_part) :: R(p_cache_size), rR, betadot(p_p_dim,p_cache_size)
  integer :: tit, it, titprev, c, i, l, npp
  integer :: iters(3), iter1, iter2, nf1
  integer(p_int64) :: iit, nf1_long
  real(p_k_part) :: div_charge(p_cache_size), radi

  real(p_k_part), dimension(p_cache_size) :: nmb1, nmb2, nmb3
# 910 "radio/os-radio-sphedet.f03"
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

    !Compute the non-normalized direction vector x-pos
    radi=this%dmin(4)+(iters(3)-1)*this%d4
    e4_radi = this%e4(:,iters(1),iters(2)) * radi
    do i = 1, np
      n(1,i) = x(1,i) - e4_radi(1)
      n(2,i) = x(2,i) - e4_radi(2)
      n(3,i) = x(3,i) - e4_radi(3)
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
      n(1,i) = xprev(1,i) - e4_radi(1)
      n(2,i) = xprev(2,i) - e4_radi(2)
      n(3,i) = xprev(3,i) - e4_radi(3)
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
# 1004 "radio/os-radio-sphedet.f03"
      div_charge(i)=div_charge(i)/((1-betaprev_c(1,i)*n(1,i)-betaprev_c(2,i)*n(2,i)&
                                     -betaprev_c(3,i)*n(3,i))**3)
    enddo

    do i = 1, npp

      Emf(1,i)=div_charge(i)*(n(2,i)*(betadot_c(2,i)*nmb1(i)-betadot_c(1,i)*nmb2(i)) &
            +n(3,i)*(betadot_c(3,i)*nmb1(i)-betadot_c(1,i)*nmb3(i)))






      Emf(2,i)=div_charge(i)*(betadot_c(3,i)*(-nmb1(i)*n(1,i)-nmb2(i)*n(2,i)) &
             + nmb3(i)*(betadot_c(1,i)*n(1,i)+betadot_c(2,i)*n(2,i)))
# 1036 "radio/os-radio-sphedet.f03"
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

!-----------------------------------------------------------------------------------------
subroutine calc_comp_approx_E1E3(this,x,beta,time,xprev,betaprev,charge,np,tid)
  implicit none
  class(sphe_Detector), intent(inout) :: this
  real(p_k_part), dimension(3,p_cache_size), intent(in) :: x, xprev
  real(p_k_part), dimension(p_p_dim,p_cache_size), intent(in) :: beta, betaprev
  real(p_k_part), dimension(p_cache_size), intent(in) :: charge
  real(p_double), intent(in) :: time
  integer, intent(in) :: np, tid

  ! n is the direction vector and R is the distance
  ! tdep time of deposition
  real(p_k_part) :: Emf(2,p_cache_size), n(3,p_cache_size)
  real(p_double) :: tdep(p_cache_size), tdepprev(p_cache_size), tprev, rsim_dt, temp
  real(p_double) :: tdep_c(p_cache_size), tdepprev_c(p_cache_size)
  real(p_k_part) :: betaprev_c(p_p_dim,p_cache_size), betadot_c(p_p_dim,p_cache_size)
  real(p_k_part) :: R, betadot(p_p_dim,p_cache_size)
  integer :: tit, it, titprev, c, i, l, npp
  integer :: iters(3), iter1, iter2, nf1
  integer(p_int64) :: iit, nf1_long
  real(p_k_part) :: div_charge(p_cache_size)

  real(p_k_part), dimension(p_cache_size) :: nmb1, nmb2, nmb3
# 1102 "radio/os-radio-sphedet.f03"
  rsim_dt = 1.0_p_double / this%sim_dt
  do i = 1, np
    betadot(1,i) = (beta(1,i)-betaprev(1,i)) * rsim_dt
    betadot(2,i) = (beta(2,i)-betaprev(2,i)) * rsim_dt
    betadot(3,i) = (beta(3,i)-betaprev(3,i)) * rsim_dt
  enddo

  tprev=time-this%sim_dt
  R=this%dmin(4)
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

    do i = 1, np
      n(1,i) = this%e4(1,iters(1),iters(2))
      n(2,i) = this%e4(2,iters(1),iters(2))
      n(3,i) = this%e4(3,iters(1),iters(2))
    enddo
    do i = 1, np
      tdep(i) = time - n(1,i)*x(1,i) - n(2,i)*x(2,i) - n(3,i)*x(3,i)
      tdepprev(i) = tprev - n(1,i)*xprev(1,i) - n(2,i)*xprev(2,i) - n(3,i)*xprev(3,i)
    enddo

    l = 1
    do i = 1, np
      if ((tdepprev(i) >= this%tpos(1)) .and. (tdep(i) <= this%tpos(this%n1+1))) then

        betaprev_c(1,l) = betaprev(1,i)
        betaprev_c(2,l) = betaprev(2,i)
        betaprev_c(3,l) = betaprev(3,i)
        betadot_c(1,l) = betadot(1,i)
        betadot_c(2,l) = betadot(2,i)
        betadot_c(3,l) = betadot(3,i)
        div_charge(l) = charge(i) ! will be modified again later
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
# 1171 "radio/os-radio-sphedet.f03"
      div_charge(i)=div_charge(i)/(R*(1-betaprev_c(1,i)*n(1,i)-betaprev_c(2,i)*n(2,i)&
                                       -betaprev_c(3,i)*n(3,i))**3)
    enddo

    do i = 1, npp

      Emf(1,i)=div_charge(i)*(n(2,i)*(betadot_c(2,i)*nmb1(i)-betadot_c(1,i)*nmb2(i)) &
            +n(3,i)*(betadot_c(3,i)*nmb1(i)-betadot_c(1,i)*nmb3(i)))






      Emf(2,i)=div_charge(i)*(betadot_c(3,i)*(-nmb1(i)*n(1,i)-nmb2(i)*n(2,i)) &
             + nmb3(i)*(betadot_c(1,i)*n(1,i)+betadot_c(2,i)*n(2,i)))
# 1203 "radio/os-radio-sphedet.f03"
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

end subroutine calc_comp_approx_E1E3
!-----------------------------------------------------------------------------------------
# 449 "radio/os-radio-sphedet.f03" 2





# 1 "radio/os-radio-sphedet.f03" 1
! if is defined then read the template definition at the end of the file
# 876 "radio/os-radio-sphedet.f03"
!-----------------------------------------------------------------------------------------
subroutine calc_comp_E2E3(this,x,beta,time,xprev,betaprev,charge,np,tid)
  implicit none
  class(sphe_Detector), intent(inout) :: this
  real(p_k_part), dimension(3,p_cache_size), intent(in) :: x, xprev
  real(p_k_part), dimension(p_p_dim,p_cache_size), intent(in) :: beta, betaprev
  real(p_k_part), dimension(p_cache_size), intent(in) :: charge
  real(p_double), intent(in) :: time
  integer, intent(in) :: np, tid

  ! n is the direction vector and R is the distance
  ! tdep time of deposition
  real(p_k_part) :: Emf(2,p_cache_size), n(3,p_cache_size), e4_radi(3)
  real(p_double) :: tdep(p_cache_size), tdepprev(p_cache_size), tprev, rsim_dt, temp
  real(p_double) :: tdep_c(p_cache_size), tdepprev_c(p_cache_size)
  real(p_k_part) :: betaprev_c(p_p_dim,p_cache_size), betadot_c(p_p_dim,p_cache_size)
  real(p_k_part) :: R(p_cache_size), rR, betadot(p_p_dim,p_cache_size)
  integer :: tit, it, titprev, c, i, l, npp
  integer :: iters(3), iter1, iter2, nf1
  integer(p_int64) :: iit, nf1_long
  real(p_k_part) :: div_charge(p_cache_size), radi

  real(p_k_part), dimension(p_cache_size) :: nmb1, nmb2, nmb3
# 910 "radio/os-radio-sphedet.f03"
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

    !Compute the non-normalized direction vector x-pos
    radi=this%dmin(4)+(iters(3)-1)*this%d4
    e4_radi = this%e4(:,iters(1),iters(2)) * radi
    do i = 1, np
      n(1,i) = x(1,i) - e4_radi(1)
      n(2,i) = x(2,i) - e4_radi(2)
      n(3,i) = x(3,i) - e4_radi(3)
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
      n(1,i) = xprev(1,i) - e4_radi(1)
      n(2,i) = xprev(2,i) - e4_radi(2)
      n(3,i) = xprev(3,i) - e4_radi(3)
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
# 1004 "radio/os-radio-sphedet.f03"
      div_charge(i)=div_charge(i)/((1-betaprev_c(1,i)*n(1,i)-betaprev_c(2,i)*n(2,i)&
                                     -betaprev_c(3,i)*n(3,i))**3)
    enddo

    do i = 1, npp





      Emf(1,i)=div_charge(i)*(nmb2(i)*(betadot_c(1,i)*n(1,i)+betadot_c(3,i)*n(3,i)) &
             + betadot_c(2,i)*(-nmb1(i)*n(1,i)-nmb3(i)*n(3,i)))


      Emf(2,i)=div_charge(i)*(betadot_c(3,i)*(-nmb1(i)*n(1,i)-nmb2(i)*n(2,i)) &
             + nmb3(i)*(betadot_c(1,i)*n(1,i)+betadot_c(2,i)*n(2,i)))
# 1036 "radio/os-radio-sphedet.f03"
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

!-----------------------------------------------------------------------------------------
subroutine calc_comp_approx_E2E3(this,x,beta,time,xprev,betaprev,charge,np,tid)
  implicit none
  class(sphe_Detector), intent(inout) :: this
  real(p_k_part), dimension(3,p_cache_size), intent(in) :: x, xprev
  real(p_k_part), dimension(p_p_dim,p_cache_size), intent(in) :: beta, betaprev
  real(p_k_part), dimension(p_cache_size), intent(in) :: charge
  real(p_double), intent(in) :: time
  integer, intent(in) :: np, tid

  ! n is the direction vector and R is the distance
  ! tdep time of deposition
  real(p_k_part) :: Emf(2,p_cache_size), n(3,p_cache_size)
  real(p_double) :: tdep(p_cache_size), tdepprev(p_cache_size), tprev, rsim_dt, temp
  real(p_double) :: tdep_c(p_cache_size), tdepprev_c(p_cache_size)
  real(p_k_part) :: betaprev_c(p_p_dim,p_cache_size), betadot_c(p_p_dim,p_cache_size)
  real(p_k_part) :: R, betadot(p_p_dim,p_cache_size)
  integer :: tit, it, titprev, c, i, l, npp
  integer :: iters(3), iter1, iter2, nf1
  integer(p_int64) :: iit, nf1_long
  real(p_k_part) :: div_charge(p_cache_size)

  real(p_k_part), dimension(p_cache_size) :: nmb1, nmb2, nmb3
# 1102 "radio/os-radio-sphedet.f03"
  rsim_dt = 1.0_p_double / this%sim_dt
  do i = 1, np
    betadot(1,i) = (beta(1,i)-betaprev(1,i)) * rsim_dt
    betadot(2,i) = (beta(2,i)-betaprev(2,i)) * rsim_dt
    betadot(3,i) = (beta(3,i)-betaprev(3,i)) * rsim_dt
  enddo

  tprev=time-this%sim_dt
  R=this%dmin(4)
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

    do i = 1, np
      n(1,i) = this%e4(1,iters(1),iters(2))
      n(2,i) = this%e4(2,iters(1),iters(2))
      n(3,i) = this%e4(3,iters(1),iters(2))
    enddo
    do i = 1, np
      tdep(i) = time - n(1,i)*x(1,i) - n(2,i)*x(2,i) - n(3,i)*x(3,i)
      tdepprev(i) = tprev - n(1,i)*xprev(1,i) - n(2,i)*xprev(2,i) - n(3,i)*xprev(3,i)
    enddo

    l = 1
    do i = 1, np
      if ((tdepprev(i) >= this%tpos(1)) .and. (tdep(i) <= this%tpos(this%n1+1))) then

        betaprev_c(1,l) = betaprev(1,i)
        betaprev_c(2,l) = betaprev(2,i)
        betaprev_c(3,l) = betaprev(3,i)
        betadot_c(1,l) = betadot(1,i)
        betadot_c(2,l) = betadot(2,i)
        betadot_c(3,l) = betadot(3,i)
        div_charge(l) = charge(i) ! will be modified again later
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
# 1171 "radio/os-radio-sphedet.f03"
      div_charge(i)=div_charge(i)/(R*(1-betaprev_c(1,i)*n(1,i)-betaprev_c(2,i)*n(2,i)&
                                       -betaprev_c(3,i)*n(3,i))**3)
    enddo

    do i = 1, npp





      Emf(1,i)=div_charge(i)*(nmb2(i)*(betadot_c(1,i)*n(1,i)+betadot_c(3,i)*n(3,i)) &
             + betadot_c(2,i)*(-nmb1(i)*n(1,i)-nmb3(i)*n(3,i)))


      Emf(2,i)=div_charge(i)*(betadot_c(3,i)*(-nmb1(i)*n(1,i)-nmb2(i)*n(2,i)) &
             + nmb3(i)*(betadot_c(1,i)*n(1,i)+betadot_c(2,i)*n(2,i)))
# 1203 "radio/os-radio-sphedet.f03"
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

end subroutine calc_comp_approx_E2E3
!-----------------------------------------------------------------------------------------
# 455 "radio/os-radio-sphedet.f03" 2






# 1 "radio/os-radio-sphedet.f03" 1
! if is defined then read the template definition at the end of the file
# 876 "radio/os-radio-sphedet.f03"
!-----------------------------------------------------------------------------------------
subroutine calc_comp_E1E2E3(this,x,beta,time,xprev,betaprev,charge,np,tid)
  implicit none
  class(sphe_Detector), intent(inout) :: this
  real(p_k_part), dimension(3,p_cache_size), intent(in) :: x, xprev
  real(p_k_part), dimension(p_p_dim,p_cache_size), intent(in) :: beta, betaprev
  real(p_k_part), dimension(p_cache_size), intent(in) :: charge
  real(p_double), intent(in) :: time
  integer, intent(in) :: np, tid

  ! n is the direction vector and R is the distance
  ! tdep time of deposition
  real(p_k_part) :: Emf(3,p_cache_size), n(3,p_cache_size), e4_radi(3)
  real(p_double) :: tdep(p_cache_size), tdepprev(p_cache_size), tprev, rsim_dt, temp
  real(p_double) :: tdep_c(p_cache_size), tdepprev_c(p_cache_size)
  real(p_k_part) :: betaprev_c(p_p_dim,p_cache_size), betadot_c(p_p_dim,p_cache_size)
  real(p_k_part) :: R(p_cache_size), rR, betadot(p_p_dim,p_cache_size)
  integer :: tit, it, titprev, c, i, l, npp
  integer :: iters(3), iter1, iter2, nf1
  integer(p_int64) :: iit, nf1_long
  real(p_k_part) :: div_charge(p_cache_size), radi

  real(p_k_part), dimension(p_cache_size) :: nmb1, nmb2, nmb3
# 910 "radio/os-radio-sphedet.f03"
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

    !Compute the non-normalized direction vector x-pos
    radi=this%dmin(4)+(iters(3)-1)*this%d4
    e4_radi = this%e4(:,iters(1),iters(2)) * radi
    do i = 1, np
      n(1,i) = x(1,i) - e4_radi(1)
      n(2,i) = x(2,i) - e4_radi(2)
      n(3,i) = x(3,i) - e4_radi(3)
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
      n(1,i) = xprev(1,i) - e4_radi(1)
      n(2,i) = xprev(2,i) - e4_radi(2)
      n(3,i) = xprev(3,i) - e4_radi(3)
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
# 1004 "radio/os-radio-sphedet.f03"
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
# 1036 "radio/os-radio-sphedet.f03"
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

!-----------------------------------------------------------------------------------------
subroutine calc_comp_approx_E1E2E3(this,x,beta,time,xprev,betaprev,charge,np,tid)
  implicit none
  class(sphe_Detector), intent(inout) :: this
  real(p_k_part), dimension(3,p_cache_size), intent(in) :: x, xprev
  real(p_k_part), dimension(p_p_dim,p_cache_size), intent(in) :: beta, betaprev
  real(p_k_part), dimension(p_cache_size), intent(in) :: charge
  real(p_double), intent(in) :: time
  integer, intent(in) :: np, tid

  ! n is the direction vector and R is the distance
  ! tdep time of deposition
  real(p_k_part) :: Emf(3,p_cache_size), n(3,p_cache_size)
  real(p_double) :: tdep(p_cache_size), tdepprev(p_cache_size), tprev, rsim_dt, temp
  real(p_double) :: tdep_c(p_cache_size), tdepprev_c(p_cache_size)
  real(p_k_part) :: betaprev_c(p_p_dim,p_cache_size), betadot_c(p_p_dim,p_cache_size)
  real(p_k_part) :: R, betadot(p_p_dim,p_cache_size)
  integer :: tit, it, titprev, c, i, l, npp
  integer :: iters(3), iter1, iter2, nf1
  integer(p_int64) :: iit, nf1_long
  real(p_k_part) :: div_charge(p_cache_size)

  real(p_k_part), dimension(p_cache_size) :: nmb1, nmb2, nmb3
# 1102 "radio/os-radio-sphedet.f03"
  rsim_dt = 1.0_p_double / this%sim_dt
  do i = 1, np
    betadot(1,i) = (beta(1,i)-betaprev(1,i)) * rsim_dt
    betadot(2,i) = (beta(2,i)-betaprev(2,i)) * rsim_dt
    betadot(3,i) = (beta(3,i)-betaprev(3,i)) * rsim_dt
  enddo

  tprev=time-this%sim_dt
  R=this%dmin(4)
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

    do i = 1, np
      n(1,i) = this%e4(1,iters(1),iters(2))
      n(2,i) = this%e4(2,iters(1),iters(2))
      n(3,i) = this%e4(3,iters(1),iters(2))
    enddo
    do i = 1, np
      tdep(i) = time - n(1,i)*x(1,i) - n(2,i)*x(2,i) - n(3,i)*x(3,i)
      tdepprev(i) = tprev - n(1,i)*xprev(1,i) - n(2,i)*xprev(2,i) - n(3,i)*xprev(3,i)
    enddo

    l = 1
    do i = 1, np
      if ((tdepprev(i) >= this%tpos(1)) .and. (tdep(i) <= this%tpos(this%n1+1))) then

        betaprev_c(1,l) = betaprev(1,i)
        betaprev_c(2,l) = betaprev(2,i)
        betaprev_c(3,l) = betaprev(3,i)
        betadot_c(1,l) = betadot(1,i)
        betadot_c(2,l) = betadot(2,i)
        betadot_c(3,l) = betadot(3,i)
        div_charge(l) = charge(i) ! will be modified again later
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
# 1171 "radio/os-radio-sphedet.f03"
      div_charge(i)=div_charge(i)/(R*(1-betaprev_c(1,i)*n(1,i)-betaprev_c(2,i)*n(2,i)&
                                       -betaprev_c(3,i)*n(3,i))**3)
    enddo

    do i = 1, npp

      Emf(1,i)=div_charge(i)*(n(2,i)*(betadot_c(2,i)*nmb1(i)-betadot_c(1,i)*nmb2(i)) &
            +n(3,i)*(betadot_c(3,i)*nmb1(i)-betadot_c(1,i)*nmb3(i)))


      Emf(2,i)=div_charge(i)*(nmb2(i)*(betadot_c(1,i)*n(1,i)+betadot_c(3,i)*n(3,i)) &
             + betadot_c(2,i)*(-nmb1(i)*n(1,i)-nmb3(i)*n(3,i)))


      Emf(3,i)=div_charge(i)*(betadot_c(3,i)*(-nmb1(i)*n(1,i)-nmb2(i)*n(2,i)) &
             + nmb3(i)*(betadot_c(1,i)*n(1,i)+betadot_c(2,i)*n(2,i)))
# 1203 "radio/os-radio-sphedet.f03"
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

end subroutine calc_comp_approx_E1E2E3
!-----------------------------------------------------------------------------------------
# 462 "radio/os-radio-sphedet.f03" 2
# 871 "radio/os-radio-sphedet.f03"
end module class_spheDet
