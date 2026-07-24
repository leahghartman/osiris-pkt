# 1 "vdf/os-vdf.f03"
# 1 "<built-in>" 1
# 1 "<built-in>" 3
# 467 "<built-in>" 3
# 1 "<command line>" 1
# 1 "<built-in>" 2
# 1 "vdf/os-vdf.f03" 2
! vdf module

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
# 4 "vdf/os-vdf.f03" 2
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
# 5 "vdf/os-vdf.f03" 2


module m_vdf

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
# 10 "vdf/os-vdf.f03" 2

  use m_vdf_define

  use m_parameters

  implicit none

  private

  interface link_vdf
    module procedure link_vdf
  end interface

  interface move_window
    module procedure move_window_vdf
  end interface

  interface copy
    module procedure copy_scalar_double_vdf_range
    module procedure copy_scalar_single_vdf_range
  end interface

  ! local methods
  public :: move_window, link_vdf, copy

contains

subroutine link_vdf( this, vdf_source)
  type( t_vdf ), intent( inout ) :: this
  type( t_vdf ), intent( in ) :: vdf_source

  this%x_dim_ = vdf_source%x_dim_
  this%f_dim_ = vdf_source%f_dim_
  this%nx_ = vdf_source%nx_
  this%gc_num_ = vdf_source%gc_num_
  this%dx_ = vdf_source%dx_

  this%f1 => vdf_source%f1
  this%f2 => vdf_source%f2
  this%f3 => vdf_source%f3
  this%buffer => vdf_source%buffer

end subroutine link_vdf


!---------------------------------------------------
subroutine move_window_vdf( this, space )
!---------------------------------------------------
! Shifts the data in the vdf object according to the
! space object. The cells left blank are set to the
! supplied value ( or 0.0 if no value is supplied )
!---------------------------------------------------

  use m_space

  implicit none

  ! dummy variables

  class( t_vdf ), intent( inout ) :: this
  type( t_space ), intent( in ) :: space ! local or global, only nx_move is needed

  ! local variables
  integer, dimension(p_max_dim) :: lnx_move
  integer :: i_dim, i1, i2, i3, lb, ub

  ! executable statements

  lnx_move( 1:this%x_dim_ ) = nx_move( space )

  ! check if move is legal
  ! this is for debug purposes and should be removed from production code
  if (x_dim(space) /= this%x_dim_) then
    write(err_buf__,*) "The dimensions of the space object, ", x_dim(space) ," do not";call err__("vdf/os-vdf.f03",83)
    write(err_buf__,*) "match the dimensions of the vdf object, ", this%x_dim_;call err__("vdf/os-vdf.f03",84)
    call abort_program( p_err_invalid )
  endif

  do i_dim = 1, this%x_dim_
    if(lnx_move(i_dim) < 0) then
      write(err_buf__,*) "Illegal move, only positive moves allowed";call err__("vdf/os-vdf.f03",90)
      call abort_program( p_err_invalid )
    endif
  enddo

  do i_dim = 1, this%x_dim_
    if ((lnx_move(i_dim) > this%gc_num_(1,i_dim)) .or. &
        (lnx_move(i_dim) > this%gc_num_(2,i_dim))) then
       write(err_buf__,*) "Illegal move, too many cells along direction ",i_dim;call err__("vdf/os-vdf.f03",98)
       write(err_buf__,*) "gc_num(:,",i_dim,") = ",this%gc_num_(:,i_dim);call err__("vdf/os-vdf.f03",99)
       write(err_buf__,*) "lnx_move(",i_dim,")  = ",lnx_move(i_dim);call err__("vdf/os-vdf.f03",100)
       call abort_program(p_err_invalid)
    endif
  enddo


  ! shift the data according to the supplied parameters

  select case(this%x_dim_)

  case(1) ! 1D vdf

     if (lnx_move(1) > 0) then ! shift left
       lb = lbound(this%f1,2)
       ub = ubound(this%f1,2) - lnx_move(1)
       do i1 = lb, ub
         this%f1(:,i1) = this%f1(:,i1+lnx_move(1))
       enddo

         ! fill empty cells with 0.0
         this%f1(:,ub+1:) = 0.0_p_k_fld
     endif
  case(2) ! 2D vdf

     if (lnx_move(1) > 0) then ! shift left
       ! new version
       lb = lbound(this%f2,2)
       ub = ubound(this%f2,2) - lnx_move(1)
       do i1 = lb, ub
         this%f2(:,i1,:) = this%f2(:,i1+lnx_move(1),:)
       enddo
       this%f2(:,ub+1:,:) = 0.0_p_k_fld

     endif

     if (lnx_move(2) > 0) then ! shift left
       lb = lbound(this%f2,3)
       ub = ubound(this%f2,3) - lnx_move(2)
       do i2 = lb, ub
         this%f2(:,:,i2) = this%f2(:,:,i2+lnx_move(2))
       enddo

     this%f2(:,:,ub+1: ) = 0.0_p_k_fld
     endif

  case(3) ! 3D vdf

     if (lnx_move(1) > 0) then ! shift left
     lb = lbound(this%f3,2)
     ub = ubound(this%f3,2) - lnx_move(1)

       if ( this%f_dim_ == 3 ) then

          !$omp parallel do private(i2,i1)
          do i3 = lbound(this%f3,4), ubound(this%f3,4)
        do i2 = lbound(this%f3,3), ubound(this%f3,3)
          do i1 = lb, ub
             this%f3(1,i1,i2,i3) = this%f3(1,i1 + lnx_move(1),i2,i3)
             this%f3(2,i1,i2,i3) = this%f3(2,i1 + lnx_move(1),i2,i3)
             this%f3(3,i1,i2,i3) = this%f3(3,i1 + lnx_move(1),i2,i3)
          enddo
          this%f3(1,ub+1,i2,i3) = 0.0_p_k_fld
          this%f3(2,ub+1,i2,i3) = 0.0_p_k_fld
          this%f3(3,ub+1,i2,i3) = 0.0_p_k_fld
        enddo
          enddo
          !$omp end parallel do

       else

          do i3 = lbound(this%f3,4), ubound(this%f3,4)
        do i2 = lbound(this%f3,3), ubound(this%f3,3)
          do i1 = lb, ub
             this%f3(:,i1,i2,i3) = this%f3(:,i1 + lnx_move(1),i2,i3)
          enddo
          this%f3(:,ub+1,i2,i3) = 0.0_p_k_fld
        enddo
          enddo

         endif

     endif

     if (lnx_move(2) > 0) then ! shift left
       lb = lbound(this%f3,3)
       ub = ubound(this%f3,3) - lnx_move(2)
       do i2 = lb, ub
         this%f3(:,:,i2,:) = this%f3(:,:,i2+lnx_move(2),:)
       enddo

       this%f3(:,:,ub+1:,:) = 0.0_p_k_fld
     endif

     if (lnx_move(3) > 0) then ! shift left
       lb = lbound(this%f3,4)
       ub = ubound(this%f3,4) - lnx_move(3)
       do i3 = lb, ub
         this%f3(:,:,:,i3) = this%f3(:,:,:,i3+lnx_move(3))
       enddo

       this%f3(:,:,:,ub+1:) = 0.0_p_k_fld
     endif

  end select

end subroutine move_window_vdf
!---------------------------------------------------

!---------------------------------------------------
subroutine copy_scalar_double_vdf_range( vdf_a, double_b, range )
!---------------------------------------------------
! copies the value of double_b to vdf_a
!---------------------------------------------------

  implicit none

! dummy variables

  type(t_vdf), intent(inout) :: vdf_a
  real(p_double), intent(in) :: double_b
  integer, dimension(:,:), intent(in) :: range

! local variables - none

! executable statements

  select case (vdf_a%x_dim_)

  case(1)
    vdf_a%f1(:,range(1,1):range(2,1)) = real(double_b, p_k_fld)

  case(2)
    vdf_a%f2(:,range(1,1):range(2,1), &
         range(1,2):range(2,2)) = real(double_b, p_k_fld)

  case(3)
    vdf_a%f3(:,range(1,1):range(2,1), &
         range(1,2):range(2,2), &
         range(1,3):range(2,3)) = real(double_b, p_k_fld)
  end select


end subroutine copy_scalar_double_vdf_range
!---------------------------------------------------


!---------------------------------------------------
subroutine copy_scalar_single_vdf_range( vdf_a, single_b, range )
!---------------------------------------------------
! copies the value of double_b to vdf_a
!---------------------------------------------------

  implicit none

! dummy variables

  type(t_vdf), intent(inout) :: vdf_a
  real(p_single), intent(in) :: single_b
  integer, dimension(:,:), intent(in) :: range

! local variables - none

! executable statements

  select case (vdf_a%x_dim_)

  case(1)
    vdf_a%f1(:,range(1,1):range(2,1)) = real(single_b, p_k_fld)

  case(2)
    vdf_a%f2(:,range(1,1):range(2,1), &
         range(1,2):range(2,2)) = real(single_b, p_k_fld)

  case(3)
    vdf_a%f3(:,range(1,1):range(2,1), &
         range(1,2):range(2,2), &
         range(1,3):range(2,3)) = real(single_b, p_k_fld)
  end select


end subroutine copy_scalar_single_vdf_range
!---------------------------------------------------




end module m_vdf
