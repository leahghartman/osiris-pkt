# 1 "spec/os-spec-memory.f90"
# 1 "<built-in>" 1
# 1 "<built-in>" 3
# 467 "<built-in>" 3
# 1 "<command line>" 1
# 1 "<built-in>" 2
# 1 "spec/os-spec-memory.f90" 2
module m_species_memory

use m_species_define, only : t_species, t_part_idx

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
# 6 "spec/os-spec-memory.f90" 2

private

integer, private, parameter :: p_stderr = 0
integer, private, parameter :: p_int64 = selected_int_kind(10)

interface alloc
  module procedure alloc_species
  module procedure alloc_1d_species
  module procedure alloc_bound_1d_species
  module procedure alloc_part_idx
  module procedure alloc_1d_part_idx
  module procedure alloc_bound_1d_part_idx
end interface

interface freemem
  module procedure freemem_species
  module procedure freemem_1d_species
  module procedure freemem_part_idx
  module procedure freemem_1d_part_idx
end interface

public :: alloc, freemem

contains

!---------------------------------------------------------------------------------------------------
! Generate memory allocation / deallocation routines for type t_species





# 1 "./memory/mem-template.h" 1
! Generated automatically by memory.py on 2019-03-06 18:36:58
# 24 "./memory/mem-template.h"
subroutine alloc_species( p, file, line, stat )

  use iso_c_binding

  implicit none

  class( t_species ), pointer :: p

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
         "t_species" // ', ' // &
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
  call log_allocation( addr, "t_species" , &
                        bsize, &
                        file, line )

end subroutine alloc_species



subroutine alloc_1d_species( p, n, file, line, stat )

  use iso_c_binding

  implicit none

  class( t_species ), dimension(:), pointer :: p
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
             "t_species" // '(' // &
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
         "t_species" // '(' // &
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
  call log_allocation( addr, "t_species" // '('// trim(strrange(1,n)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_1d_species
# 270 "./memory/mem-template.h"
subroutine alloc_bound_1d_species( p, lb, ub, file, line, stat )

  use iso_c_binding

  implicit none

  class( t_species ), dimension(:), pointer :: p
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
             "t_species" // '(' // &
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
           "t_species" // '(' // &
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
  call log_allocation( addr, "t_species" // '('// trim(strrange(1,lb,ub)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_bound_1d_species



subroutine freemem_species( p, file, line, stat )

  use iso_c_binding

  implicit none

  class( t_species ), pointer :: p

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
   call log_deallocation( addr, "t_species", &
                        bsize, file, line )

   ! deallocate memory





   deallocate( p, stat = ierr )


   ! Check allocation
   if ( ierr /= 0 ) then
     write( p_stderr, * ) '(*error*) Unable to deallocate memory, ' // &
            "t_species" // ', ' // &
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

end subroutine freemem_species



subroutine freemem_1d_species( p, file, line, stat )

  use iso_c_binding

  implicit none

  class( t_species ), dimension(:), pointer :: p

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
     call log_deallocation( addr, "t_species" // '(:)', &
                        bsize, file, line )

   ! deallocate memory





   deallocate( p, stat = ierr )


   ! Check allocation
   if ( ierr /= 0 ) then
     write( p_stderr, * ) '(*error*) Unable to deallocate memory, ' // &
            "t_species" // '(' // &
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

end subroutine freemem_1d_species
# 40 "spec/os-spec-memory.f90" 2

!---------------------------------------------------------------------------------------------------
! Generate memory allocation / deallocation routines for type t_part_idx





# 1 "./memory/mem-template.h" 1
! Generated automatically by memory.py on 2019-03-06 18:36:58
# 24 "./memory/mem-template.h"
subroutine alloc_part_idx( p, file, line, stat )

  use iso_c_binding

  implicit none

  type( t_part_idx ), pointer :: p

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
         "t_part_idx" // ', ' // &
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
  call log_allocation( addr, "t_part_idx" , &
                        bsize, &
                        file, line )

end subroutine alloc_part_idx



subroutine alloc_1d_part_idx( p, n, file, line, stat )

  use iso_c_binding

  implicit none

  type( t_part_idx ), dimension(:), pointer :: p
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
             "t_part_idx" // '(' // &
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
         "t_part_idx" // '(' // &
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
  call log_allocation( addr, "t_part_idx" // '('// trim(strrange(1,n)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_1d_part_idx
# 270 "./memory/mem-template.h"
subroutine alloc_bound_1d_part_idx( p, lb, ub, file, line, stat )

  use iso_c_binding

  implicit none

  type( t_part_idx ), dimension(:), pointer :: p
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
             "t_part_idx" // '(' // &
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
           "t_part_idx" // '(' // &
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
  call log_allocation( addr, "t_part_idx" // '('// trim(strrange(1,lb,ub)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_bound_1d_part_idx



subroutine freemem_part_idx( p, file, line, stat )

  use iso_c_binding

  implicit none

  type( t_part_idx ), pointer :: p

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
   call log_deallocation( addr, "t_part_idx", &
                        bsize, file, line )

   ! deallocate memory





   deallocate( p, stat = ierr )


   ! Check allocation
   if ( ierr /= 0 ) then
     write( p_stderr, * ) '(*error*) Unable to deallocate memory, ' // &
            "t_part_idx" // ', ' // &
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

end subroutine freemem_part_idx



subroutine freemem_1d_part_idx( p, file, line, stat )

  use iso_c_binding

  implicit none

  type( t_part_idx ), dimension(:), pointer :: p

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
     call log_deallocation( addr, "t_part_idx" // '(:)', &
                        bsize, file, line )

   ! deallocate memory





   deallocate( p, stat = ierr )


   ! Check allocation
   if ( ierr /= 0 ) then
     write( p_stderr, * ) '(*error*) Unable to deallocate memory, ' // &
            "t_part_idx" // '(' // &
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

end subroutine freemem_1d_part_idx
# 49 "spec/os-spec-memory.f90" 2

!---------------------------------------------------------------------------------------------------

end module m_species_memory
