# 1 "spec/os-spec-piston.f90"
# 1 "<built-in>" 1
# 1 "<built-in>" 3
# 467 "<built-in>" 3
# 1 "<command line>" 1
# 1 "<built-in>" 2
# 1 "spec/os-spec-piston.f90" 2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! piston class
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
# 6 "spec/os-spec-piston.f90" 2
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
# 7 "spec/os-spec-piston.f90" 2

module m_piston

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
# 11 "spec/os-spec-piston.f90" 2

  use m_system

  use m_parameters

  use m_space
  use m_vdf

  use m_utilities
  use m_random

  use m_species_define
  use m_species_current

  use m_fparser



  implicit none

! restrict access to things explicitly declared public
  private

  ! piston parameters
  character(1), dimension(2), parameter :: p_sp_piston_transverse_symbols=(/"A","B"/)
  integer, parameter :: p_sp_piston_uniform = 1
  integer, parameter :: p_sp_piston_fparser = 2
  character(20), parameter :: p_sp_piston_uniform_str = 'uniform'
  character(20), parameter :: p_sp_piston_fparser_str = 'fparser'
  integer, parameter :: p_sp_piston_up = 1
  integer, parameter :: p_sp_piston_down = 2

  ! buffers for particle functions

  !buffers for piston
  ! indexes of prticles that crossed
  integer, pointer, dimension(:) :: piston_sp_idx => null()
  ! number of particles in piston_sp_idx
  integer :: num_idx
  !for current deposit of piston
  real(p_k_part), pointer, dimension(:,:) :: pst_curr_xold => null(), &
                                            pst_curr_xnew => null()
  real(p_k_part), pointer, dimension(:,:) :: pst_curr_udeposit => null()

  interface read_nml_piston
    module procedure read_nml_pst
  end interface

  interface check_piston
    module procedure check_pst
  end interface

  interface update_piston
    module procedure update_pst
  end interface

  interface cross_particles_piston
    module procedure cross_particles_pst
  end interface

! interface process_piston_1d
! module procedure process_pst_1d
! end interface

! interface process_piston_2d
! module procedure process_pst_2d
! end interface

! interface process_piston_3d
! module procedure process_pst_3d
! end interface

  interface init_buffers_piston
    module procedure init_buffers_pst
  end interface

  interface cleanup_buffers_piston
    module procedure cleanup_buffers_pst
  end interface

! declare things that should be public
  public :: t_piston
  public :: read_nml_piston
  public :: check_piston
  public :: cross_particles_piston
  public :: update_piston
! public :: process_piston_1d, process_piston_2d, process_piston_3d
  public :: init_buffers_piston, cleanup_buffers_piston

interface alloc
  module procedure alloc_piston
  module procedure alloc_1d_piston
  module procedure alloc_bound_1d_piston
end interface

interface freemem
  module procedure freemem_piston
  module procedure freemem_1d_piston
end interface

  public :: alloc, freemem

 contains


!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------
subroutine read_nml_pst( this, input_file )

  use m_input_file

  implicit none

  type( t_piston ), intent(inout) :: this
  class( t_input_file ), intent(inout) :: input_file

  integer :: dim
  integer :: updown
  real(p_k_part) :: v
  real(p_k_part) :: start_pos
  real(p_k_part) :: start_time
  real(p_k_part) :: stop_time
  real(p_k_part) :: opacity_factor
  character(20) :: profile_type
  character(p_max_expr_len) :: profile_expression


  namelist /nl_piston/ dim, updown, v, start_pos, start_time, stop_time, &
    opacity_factor, profile_type, profile_expression

  integer :: ierr

! executable statements

! set defaults for variables that can be omitted in
! the input file
  dim = 1
  updown = p_sp_piston_up
  v = 0.0_p_k_part
  start_pos = 0.0_p_k_part
  start_time = 0.0_p_k_part
  stop_time = 0.0_p_k_part
  opacity_factor = 1.0_p_k_part
  profile_type = "uniform"
  profile_expression = ""

  ! Get namelist text from input file
  call get_namelist( input_file, "nl_piston", ierr )
  if ( ierr /= 0 ) then
 print *, "   No piston parameters specified."
 print *, "   aborting..."
 stop
  endif

  read (input_file%nml_text, nml = nl_piston, iostat = ierr)
  if (ierr /= 0) then
 print *, "   Error reading piston parameters."
 print *, "   aborting..."
 stop
  endif

  ! setup of piston values

  if ( v <= 0.0_p_k_part) then
 print *, ""
 print *, "   Error reading piston parameters - v."
 print *, "   v should be greater than zero."
 print *, "   aborting..."
 stop
  endif

  if ( updown == p_sp_piston_down ) v = -v

  this%dim = dim
  this%updown = updown
  this%u = v
  this%start_pos = start_pos
  this%start_time = start_time
  this%stop_time = stop_time
  this%opacity_factor = opacity_factor


  if ( start_time >= stop_time) then
 print *, ""
 print *, "   Error reading piston parameters - start_/stop_ time."
 print *, "   Stop_time must be greater than start_time. "
 print *, "   aborting..."
 stop
  endif

  if ( dim > p_x_dim) then
 print *, ""
 print *, "   Error reading piston parameters - dim."
 print *, "   dim > p_x_dim. "
 print *, "   aborting..."
 stop
  endif

  if ( updown /= p_sp_piston_up .AND. updown /= p_sp_piston_down) then
 print *, ""
 print *, "   Error reading piston parameters - updown."
 print *, "   updown must be either of 1, 2. "
 print *, "   aborting..."
 stop
  endif
  select case ( profile_type )
 case ( p_sp_piston_uniform_str )
   this%profile_type = p_sp_piston_uniform
 case ( p_sp_piston_fparser_str )
   if( p_x_dim == 1) then
  print *, ""
  print *, "   Error reading piston parameters - profile_type"
  print *, "   It does not make sense to use f_parser for a 1D run."
  print *, "   aborting..."
  stop
   endif
   this%profile_type = p_sp_piston_fparser

   ! set up profile
   call setup( this%profile_parser, profile_expression, &
       p_sp_piston_transverse_symbols(1:p_x_dim-1), ierr )
   if (ierr /= 0) then
  write(err_buf__,*) 'Error seting up opacity profile for piston.';call err__("spec/os-spec-piston.f90",233)
  call abort_program(p_err_invalid)
   endif
 case default
   write(err_buf__,*) 'Invalid profile_type_piston specified';call err__("spec/os-spec-piston.f90",237)
   stop
  end select

  ! calculate constannt values

  this%squ = this%u**2 !u^2
  this%sqg = 1.0_p_k_part + this%u**2 !gamma^2
  this%g = sqrt( this%sqg ) ! gamma
  this%v = this%u / this%g ! v

  !this%sqv = this%v**2

end subroutine read_nml_pst
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
subroutine init_buffers_pst( )
!---------------------------------------------------------------------------------------------------
! initializes cache buffers for piston pusher
! routines
!---------------------------------------------------------------------------------------------------
  implicit none

  ! setup buffers for advance_deposit routines

  ! setup buffers for piston
  call alloc(piston_sp_idx, (/ p_cache_size /),"spec/os-spec-piston.f90",265)

  call alloc(pst_curr_xold, (/ p_x_dim, p_cache_size /),"spec/os-spec-piston.f90",267)

  call alloc(pst_curr_xnew, (/ p_x_dim, p_cache_size /),"spec/os-spec-piston.f90",269)

  call alloc(pst_curr_udeposit, (/ p_p_dim, p_cache_size /),"spec/os-spec-piston.f90",271)

end subroutine init_buffers_pst
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
subroutine cleanup_buffers_pst( )
!---------------------------------------------------------------------------------------------------

  implicit none

  ! setup buffers for advance_deposit routines

  ! setup buffers for piston
  call freemem(piston_sp_idx,"spec/os-spec-piston.f90",285)

  call freemem(pst_curr_xold,"spec/os-spec-piston.f90",287)

  call freemem(pst_curr_xnew,"spec/os-spec-piston.f90",289)

  call freemem(pst_curr_udeposit,"spec/os-spec-piston.f90",291)

end subroutine cleanup_buffers_pst
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
        subroutine cross_particles_pst( piston_all, np, &
                       & x, xbuf, ptrcur )
!---------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------

          implicit none

! dummy variables

! class(t_species), intent(inout) :: this
          type(t_piston),dimension(:), pointer :: piston_all
          real(p_k_part), dimension(:,:), intent(in) :: x, xbuf
          integer, intent(in) :: ptrcur
          integer, intent(in) :: np

! local variables

          integer :: i, pst_id, num_pistons, dim
          logical :: piston_crossed


! executable statements

          num_pistons=size(piston_all)

          ! generate list of particles to be processed by piston
          num_idx = 0
          do i = 1, np ! all particles
            do pst_id=1, num_pistons !all piston
              dim = piston_all(pst_id)%dim
              if ( piston_all(pst_id)%inbox ) then
                if( piston_all(pst_id)%updown == p_sp_piston_up ) then
                  piston_crossed = x(dim,ptrcur+i) > piston_all(pst_id)%pos_before .AND. &
      & xbuf(dim,i) <= piston_all(pst_id)%pos_after
                else
                  piston_crossed = x(dim,ptrcur+i) < piston_all(pst_id)%pos_before .AND. &
      & xbuf(dim,i) >= piston_all(pst_id)%pos_after
                endif !updown
                if( piston_crossed ) then
                  ! add particle to list of particles to be processed
                  num_idx = num_idx + 1
                  piston_sp_idx( num_idx ) = i
                  exit !of loop all piston
                endif !crossed
              endif !inbox
            enddo !all piston
          enddo ! all particles

        end subroutine cross_particles_pst
!---------------------------------------------------------------------------------------------------
# 905 "spec/os-spec-piston.f90"
!---------------------------------------------------------------------------------------------------
        subroutine update_pst( this, t, dt, l_space )
!---------------------------------------------------------------------------------------------------
! updates pos_before, pos_after, inbox for piston
!---------------------------------------------------------------------------------------------------

! <use-statements for this subroutine>

          implicit none

! dummy variables

          type(t_piston), intent(inout) :: this
          real(p_double), intent(in) :: t
          real(p_k_part), intent(in) :: dt
          type( t_space ), intent(in) :: l_space

! local variables

          real(p_k_part), dimension(p_x_dim) :: l_xmin, l_xmax
          real(p_k_part) :: t_before, t_after

! executable statements

          l_xmin = real( xmin(l_space), p_k_part )
          l_xmax = real( xmax(l_space), p_k_part )

          t_before = real( t - dt, p_k_part )
          t_after = real( t, p_k_part )

          if ( this%start_time <= t_after .AND. this%stop_time >= t_after ) then ! piston active during this timestep ?

            ! update pos_before and pos_after
            this%pos_before = this%start_pos + this%v * ( t_before - this%start_time )
            this%pos_after = this%start_pos + this%v * ( t_after - this%start_time )

            ! update inbox, adjust boundaries by 3 dt so we include guard cells
            this%inbox = this%pos_after >= l_xmin(this%dim) - 3 * dt .AND. &
     & this%pos_after < l_xmax(this%dim) + 3 * dt

          else
            this%inbox=.false.
          endif ! piston active in this dt

       end subroutine update_pst
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
        subroutine check_pst( this, dt )
!---------------------------------------------------------------------------------------------------
! updates pos_before, pos_after, inbox for piston
!---------------------------------------------------------------------------------------------------

! <use-statements for this subroutine>

          implicit none

! dummy variables

          type(t_piston), dimension(:), intent(in) :: this
          real(p_double), intent(in) :: dt

! local variables

          integer :: pst_id_a, pst_id_b, dima, dimb, updowna, updownb
          real(p_k_part) :: edge_a, edge_b, va, vb, tsa, tsb, intervall, t_cross

! executable statements

          do pst_id_a = 1, size(this)

            dima = this(pst_id_a)%dim
            va = this(pst_id_a)%v
            tsa = this(pst_id_a)%start_time
            updowna = this(pst_id_a)%updown
            edge_a = this(pst_id_a)%start_pos

            do pst_id_b = pst_id_a + 1, size(this)

              dimb = this(pst_id_b)%dim
              vb = this(pst_id_b)%v
              tsb = this(pst_id_b)%start_time
              updownb = this(pst_id_b)%updown
              edge_b = this(pst_id_b)%start_pos

              if( dima == dimb .AND. updowna /= updownb) then
                t_cross = (edge_b - edge_a + va*tsa - vb*tsb) / (va-vb)
                intervall = abs(real(2 * dt, p_k_part)/(va -vb))
                if( overlap( t_cross - intervall, t_cross + intervall, this(pst_id_a)%start_time, this(pst_id_a)%stop_time ) .AND. &
              & overlap( t_cross - intervall, t_cross + intervall, this(pst_id_b)%start_time, this(pst_id_b)%stop_time ) ) then
                  print *, "Two or more piston of the same species and oposite direction"
                  print *, "get too close to each other (<c * dt)"
                  print *, "Abort..."
                  stop
                endif
              endif !oposit rirections

            enddo ! remaining piston

          enddo ! all piston
       end subroutine check_pst
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
        function overlap( al, ah, bl, bh )
!---------------------------------------------------------------------------------------------------
! returns true if [al,ah] and [bl,bh] intersect
!---------------------------------------------------------------------------------------------------

! <use-statements for this subroutine>

          implicit none

! dummy variables

          real(p_k_part), intent(in) :: al,ah,bl,bh
          logical :: overlap

          overlap = ( bl <= ah ) .AND. ( al <= bh )

       end function overlap
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! Generate memory allocation / deallocation routines





# 1 "./memory/mem-template.h" 1
! Generated automatically by memory.py on 2019-03-06 18:36:58
# 24 "./memory/mem-template.h"
subroutine alloc_piston( p, file, line, stat )

  use iso_c_binding

  implicit none

  type( t_piston ), pointer :: p

  character(len=*), intent(in) :: file
  integer, intent(in) :: line
  integer, intent(out), optional :: stat

  integer :: ierr
  integer( p_int64 ) :: addr, bsize







  ! Nullify pointer
  p => null()

  ! Allocate memory



  allocate( p, stat = ierr )


  ! Check allocation
  if ( ierr /= 0 ) then
  write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
         "t_piston" // ', ' // &
           trim(strfileline(file,line))
    if ( present(stat) ) then
      stat = ierr
      return
    else
      call exit(ierr)
    endif
  endif






  ! log allocation
  addr = loc( p )
  bsize = sizeof( p )
  call log_allocation( addr, "t_piston" , &
                        bsize, &
                        file, line )

end subroutine alloc_piston



subroutine alloc_1d_piston( p, n, file, line, stat )

  use iso_c_binding

  implicit none

  type( t_piston ), dimension(:), pointer :: p
  integer, dimension(:), intent(in) :: n

  character(len=*), intent(in) :: file
  integer, intent(in) :: line
  integer, intent(out), optional :: stat

  integer :: i, ierr
  integer( p_int64 ) :: addr, bsize

  ! Error constants
  integer, parameter :: p_err_invalid_dims = -1001
# 111 "./memory/mem-template.h"
  ! Nullify pointer
  p => null()

  ! Verify requested size
  do i = 1, 1
    if ( n(i) <= 0 ) then
      write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
             'invalid dimensions requested, ' //&
             "t_piston" // '(' // &
             trim(strrange(1,n)) // '), ' // &
             trim(strfileline(file,line))

      if ( present(stat) ) then
        stat = p_err_invalid_dims
        return
      else
        call exit(p_err_invalid_dims)
      endif

    endif
  enddo

  ! Allocate memory
# 143 "./memory/mem-template.h"
  allocate( p( n(1) ), stat = ierr )


  ! Check allocation
  if ( ierr /= 0 ) then
  write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
         "t_piston" // '(' // &
         trim(strrange(1,n)) // '), ' // &
           trim(strfileline(file,line))
  if ( present(stat) ) then
    stat = ierr
    return
  else
    call exit(ierr)
  endif
  endif






  ! log allocation
  addr = loc(p)
  bsize = sizeof( p )
  call log_allocation( addr, "t_piston" // '('// trim(strrange(1,n)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_1d_piston
# 270 "./memory/mem-template.h"
subroutine alloc_bound_1d_piston( p, lb, ub, file, line, stat )

  use iso_c_binding

  implicit none

  type( t_piston ), dimension(:), pointer :: p
  integer, dimension(:), intent(in) :: lb, ub

  character(len=*), intent(in) :: file
  integer, intent(in) :: line
  integer, intent(out), optional :: stat

  integer :: i,ierr
  integer( p_int64 ) :: addr, bsize

  ! Error constants
  integer, parameter :: p_err_invalid_bounds = -1002
# 298 "./memory/mem-template.h"
  ! Nullify pointer
  p => null()

  ! Verify requested boundaries
  do i = 1, 1
    if ( ub(i) < lb(i) ) then
      write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
             'invalid boundaries requested, ' // &
             "t_piston" // '(' // &
             trim(strrange(1,lb,ub)) // '), ' // &
             trim(strfileline(file,line))

      if ( present(stat) ) then
        stat = p_err_invalid_bounds
        return
      else
        call exit(p_err_invalid_bounds)
      endif

    endif
  enddo

  ! Allocate memory
# 329 "./memory/mem-template.h"
  allocate( p( lb(1):ub(1) ), stat = ierr )


  ! Check allocation
  if ( ierr /= 0 ) then
    write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
           "t_piston" // '(' // &
           trim(strrange(1,lb,ub)) // '), ' // &
           trim(strfileline(file,line))

    if ( present(stat) ) then
      stat = ierr
      return
    else
      call exit(ierr)
    endif
  endif







  ! log allocation
  addr = loc(p)
  bsize = sizeof(p)
  call log_allocation( addr, "t_piston" // '('// trim(strrange(1,lb,ub)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_bound_1d_piston



subroutine freemem_piston( p, file, line, stat )

  use iso_c_binding

  implicit none

  type( t_piston ), pointer :: p

  character(len=*), intent(in) :: file
  integer, intent(in) :: line
  integer, intent(out), optional :: stat

  integer :: ierr
  integer( p_int64 ) :: addr, bsize

  ! check if p is associated, otherwise return silently
  if ( associated(p) ) then
   ! get memory size to deallocate
   bsize = sizeof(p)
   addr = loc(p)

   ! Decrease global memory counter
   call log_deallocation( addr, "t_piston", &
                        bsize, file, line )

   ! deallocate memory





   deallocate( p, stat = ierr )


   ! Check allocation
   if ( ierr /= 0 ) then
     write( p_stderr, * ) '(*error*) Unable to deallocate memory, ' // &
            "t_piston" // ', ' // &
        trim(strfileline(file,line))
     if ( present(stat) ) then
     stat = ierr
     return
     else
     call exit(ierr)
     endif
   endif

   ! Nullify pointer
   p => null()

  endif

end subroutine freemem_piston



subroutine freemem_1d_piston( p, file, line, stat )

  use iso_c_binding

  implicit none

  type( t_piston ), dimension(:), pointer :: p

  character(len=*), intent(in) :: file
  integer, intent(in) :: line
  integer, intent(out), optional :: stat

  integer :: ierr
  integer( p_int64 ) :: addr, bsize

  ! check if p is associated, otherwise return silently
  if ( associated(p) ) then
   ! get memory size to deallocate
   bsize = sizeof(p)
   addr = loc(p)

   ! Decrease global memory counter
     call log_deallocation( addr, "t_piston" // '(:)', &
                        bsize, file, line )

   ! deallocate memory





   deallocate( p, stat = ierr )


   ! Check allocation
   if ( ierr /= 0 ) then
     write( p_stderr, * ) '(*error*) Unable to deallocate memory, ' // &
            "t_piston" // '(' // &
        trim(strrange(1,lbound(p),ubound(p))) // '), ' // &
        trim(strfileline(file,line))
     if ( present(stat) ) then
     stat = ierr
     return
     else
     call exit(ierr)
     endif
   endif

   ! Nullify pointer
   p => null()

  endif

end subroutine freemem_1d_piston
# 1038 "spec/os-spec-piston.f90" 2

!---------------------------------------------------------------------------------------------------

      end module m_piston
