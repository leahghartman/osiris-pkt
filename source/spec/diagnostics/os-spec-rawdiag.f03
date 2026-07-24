#include "os-config.h"
#include "os-preprocess.fpp"

module m_species_rawdiag

#include "memory/memory.h"

use m_species_define
use m_fparser

use stringutil
use m_diagnostic_utilities
use m_node_conf

use m_space
use m_grid_define
use m_grid, only : my_nx_p
use m_time_step, only : t_time_step, n
use m_random
use m_boosted_diag, only : t_boosted_diag, generate_sort_idx

use m_parameters

implicit none

private

interface write_raw
  module procedure write_raw_
end interface

interface extract_boosted_raw
  module procedure extract_boosted_raw_
end interface

interface write_boosted_raw
  module procedure write_boosted_raw_
end interface

public :: write_raw, write_raw_select_data, write_raw_create_file
public :: create_datasets, close_datasets
public :: extract_boosted_raw, write_boosted_raw

contains

!-------------------------------------------------------------------------------
! Get index of selected particles for output
!-------------------------------------------------------------------------------
subroutine write_raw_select_data( spec, t, idx, num_par )
!-------------------------------------------------------------------------------

  implicit none

  class( t_species ), intent(inout) :: spec
  real(p_double),      intent(in) :: t
  integer, dimension(:), intent(out) :: idx
  integer, intent(out) :: num_par

  ! local variables
  real(p_k_part) :: gamma_limit_sqr, gamma_par_sqr, rnd
  real(p_k_part), dimension(4) :: pos ! was 3
  real(p_k_fparse) :: raw_eval
  real(p_k_fparse), dimension(p_p_dim + spec%get_n_x_dims() + 2) :: raw_var
  logical :: has_raw_math_expr
  integer :: l, lbuf, n_x_dim

  n_x_dim = spec%get_n_x_dims()

  ! initialize gamma limit
  gamma_limit_sqr = spec%diag%raw_gamma_limit**2

  ! loop over all particles

  l = 1
  lbuf = 0

  ! Copy all the particles to the temp buffer

  ! time for raw_math evaluation
  raw_var(n_x_dim + p_p_dim + 2) = t

  has_raw_math_expr = (spec%diag%raw_math_expr /= '')

  do
     if (l > spec%num_par)  exit
     gamma_par_sqr = 1.0_p_k_part + spec%p(1,l)**2 + spec%p(2,l)**2 + spec%p(3,l)**2

     if ( gamma_par_sqr >= gamma_limit_sqr ) then

        raw_eval = 1
        if (has_raw_math_expr) then

           ! fill evaluation variables
           call spec % get_position( l, pos )                    ! x1-x3
           raw_var(1:n_x_dim) = real( pos(1:n_x_dim), p_k_fparse )
           raw_var(n_x_dim+1) = spec%p(1,l)                     ! p1
           raw_var(n_x_dim+2) = spec%p(2,l)                     ! p2
           raw_var(n_x_dim+3) = spec%p(3,l)                     ! p3
           raw_var(n_x_dim+p_p_dim + 1) = sqrt(gamma_par_sqr)   ! g

           ! evaluate
           raw_eval = eval( spec%diag%raw_func, raw_var )

        endif

        if (raw_eval > 0) then

              call rng % harvest_real2( rnd, spec%ix(:,l) )
              if ( rnd <= spec%diag%raw_fraction ) then

                 lbuf = lbuf + 1
                 idx(lbuf) = l

              endif

         endif
     endif
     l = l+1
  enddo

  num_par = lbuf

end subroutine write_raw_select_data
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
subroutine write_raw_create_file( spec, total_par, n, t, n_aux, diagFile, comm )
!-------------------------------------------------------------------------------

  implicit none

  class( t_species ), intent(inout) :: spec
  integer,                   intent(in) :: n
  integer, intent(in) :: total_par
  real(p_double),            intent(in) :: t
  integer,                   intent(in) :: n_aux
  integer, intent(in) :: comm

  class( t_diag_file ), allocatable, intent(inout) :: diagFile

  ! local variables
  integer :: i, j, n_x_dim

  n_x_dim = spec%get_n_x_dims()
  ! create file

  ! Prepare path and file name
  call create_diag_file( diagFile )
  diagFile % ftype = p_diag_particles

  diagFile%filepath = trim(path_mass)//'RAW'//p_dir_sep//replace_blanks(trim(spec%name))//p_dir_sep
  diagFile%filename = get_filename( n_aux, '' // 'RAW'  // '-' // replace_blanks(trim(spec%name)) )

  diagFile%name = spec%name

  diagFile%iter%n         = n
  diagFile%iter%t         = t
  diagFile%iter%time_units = '1 / \omega_p'

  diagFile%particles%name = spec%name

  ! Species object doesn't have a label (yet) so just set it to the name
  diagFile%particles%label = spec%name

  diagFile%particles%np = total_par
  diagFile%particles%has_tags = spec%add_tag

  j = 0

  ! positions
  do i = 1, n_x_dim
    j = j+1
    diagFile%particles%quants(j) = 'x'//(char(iachar('0')+i))
    diagFile%particles%qlabels(j) = 'x_'//(char(iachar('0')+i))
    diagFile%particles%qunits(j) = 'c/\omega_p'
    diagFile%particles%offset_t(j) = 0.0_p_double
  enddo

  ! momenta
  do i = 1, p_p_dim
    j = j+1
    diagFile%particles%quants(j) = 'p'//(char(iachar('0')+i))
    diagFile%particles%qlabels(j) = 'p_'//(char(iachar('0')+i))
    diagFile%particles%qunits(j) = 'm_e c'
    diagFile%particles%offset_t(j) = -0.5_p_double
  enddo

  ! charge
  j = j+1
  diagFile%particles%quants(j) = 'q'
  diagFile%particles%qlabels(j) = 'q'
  diagFile%particles%qunits(j) = 'e'
  diagFile%particles%offset_t(j) = 0.0_p_double

  ! energy
  j = j+1
  diagFile%particles%quants(j) = 'ene'
  diagFile%particles%qlabels(j) = 'Ene'
  diagFile%particles%qunits(j) = 'm_e c^2'
  diagFile%particles%offset_t(j) = -0.5_p_double

#ifdef __HAS_SPIN__

  ! spin
  do i = 1, p_s_dim
    j = j+1
    diagFile%particles%quants(j) = 's'//(char(iachar('0')+i))
    diagFile%particles%qlabels(j) = 's_'//(char(iachar('0')+i))
    diagFile%particles%qunits(j) = 'h/4\pi'
    diagFile%particles%offset_t(j) = -0.5_p_double
  enddo

#endif

  if ( diagFile%particles%has_tags ) then
    j = j+1
    diagFile%particles%quants(j) = 'tag'
    diagFile%particles%qlabels(j) = 'Tag'
    diagFile%particles%qunits(j) = ''
    diagFile%particles%offset_t(j) = 0.0_p_double
  endif

  diagFile%particles%nquants = j

  ! create the file
  call diagFile % open( p_diag_create, comm )

end subroutine write_raw_create_file
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
! Create datasets in the file
!-------------------------------------------------------------------------------
subroutine create_datasets( diagFile, total_par, datasets )
!-------------------------------------------------------------------------------

  implicit none

  class( t_diag_file ), intent(inout) :: diagFile
  integer, intent(in) :: total_par

  type( t_diag_dataset ), dimension(:), intent(inout) :: datasets

  integer :: i, j, data_type

  select case (p_diag_prec)
  case ( p_single )
    data_type = p_diag_float32
  case ( p_double )
    data_type = p_diag_float64
  end select

  if ( diagFile%particles%has_tags ) then
    j = diagFile % particles % nquants - 1
  else
    j = diagFile % particles % nquants
  endif

  ! particle quantities
  do i = 1, j
    call diagFile % start_cdset( diagFile % particles % quants(i), 1, [total_par], &
                            data_type, datasets(i) )
  enddo

  ! particle tags
  if ( diagFile%particles%has_tags ) then
    i = diagFile % particles % nquants
    call diagFile % start_cdset( 'tag', 2, [2,total_par], p_diag_int32, &
                            datasets(i) )
  endif

end subroutine create_datasets
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Close datasets in the file
!-------------------------------------------------------------------------------
subroutine close_datasets( diagFile, datasets )
!-------------------------------------------------------------------------------
  implicit none

  class( t_diag_file ), intent(inout) :: diagFile
  type( t_diag_dataset ), dimension(:), intent(inout) :: datasets

  integer :: i

  do i = 1, diagFile % particles % nquants
    call diagFile % end_cdset( datasets(i) )
  enddo

end subroutine close_datasets
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
! write raw particle data
!-------------------------------------------------------------------------------

subroutine write_raw_( spec, no_co, n, t, n_aux )

  use m_units
  !use mpi

  implicit none

  class( t_species ), intent(inout) :: spec
  class( t_node_conf ),               intent(in) :: no_co
  integer,                   intent(in) :: n
  real(p_double), intent(in) :: t
  integer,                   intent(in) :: n_aux

  class( t_diag_file ), allocatable :: diagFile
#ifdef __HAS_SPIN__
  type( t_diag_dataset ), dimension( spec%get_n_x_dims() + p_p_dim + p_s_dim + 3 ) :: datasets
#else
  type( t_diag_dataset ), dimension( spec%get_n_x_dims() + p_p_dim + 3 ) :: datasets
#endif
  type( t_diag_chunk ) :: chunk, chunk_tag
  integer :: start, rank, num_par, total_par, color, comm, ierr
  integer :: ndims, i, j, k
  integer, dimension(:), pointer :: idx
  real( p_diag_prec ), dimension(:), pointer :: write_buffer
  integer, dimension( :, : ), pointer :: write_tag_buffer
  real( p_k_part ) :: u2

  ndims = spec%get_n_x_dims()

  ! get index of particles to output
  if ( spec%num_par > 0 ) then
    call alloc(idx,  [spec%num_par] )
    call write_raw_select_data( spec, t, idx, num_par )
  else
    num_par = 0
  endif


  if ( no_co % comm_size() > 1 ) then
    comm = MPI_COMM_NULL

    call MPI_ALLREDUCE( num_par, total_par, 1, MPI_INTEGER, MPI_SUM, no_co % comm, ierr )

    if ( total_par > 0 ) then
      if ( num_par > 0 ) then
        color = 1
      else
        color = MPI_UNDEFINED
      endif
      call MPI_COMM_SPLIT( no_co % comm, color, 0, comm, ierr )

      if ( num_par > 0 ) then
        call MPI_COMM_RANK( comm, rank, ierr )
        call MPI_EXSCAN( num_par, start, 1, MPI_INTEGER, MPI_SUM, comm, ierr )
        if ( rank == 0 ) start = 0
      endif
    endif
  else
    comm = MPI_COMM_SELF
    total_par = num_par
    start = 0
  endif
  
  ! Write data, if any
  if ( total_par > 0 ) then

    ! Only nodes with particles are involved in this section
    if ( num_par > 0 ) then

      ! Create the file
      call write_raw_create_file( spec, total_par, n, t, n_aux, diagFile, comm )

      ! Create the datasets
      call create_datasets( diagFile, total_par, datasets )

      call alloc( write_buffer, [ num_par ] )
      chunk % start(1) = start
      chunk % count(1) = num_par
      chunk % stride(1) = 1
      chunk % data = c_loc( write_buffer )

      ! Positions
      k = 1
      do j = 1, ndims
        call spec % get_position( j, idx, num_par, write_buffer, spec%diag%raw_ref_box(j))
        call diagFile % write_par_cdset( datasets(k), chunk, start )
        k = k+1
      enddo

      ! Momenta
      do j = 1, p_p_dim
        do i = 1, num_par
          write_buffer(i) = real(spec%p(j,idx(i)), p_diag_prec)
        enddo
        call diagFile % write_par_cdset( datasets(k), chunk, start )
        k = k+1
      enddo

      ! Charge
      do i = 1, num_par
        write_buffer(i) = real(spec%q(idx(i)), p_diag_prec)
      enddo
      call diagFile % write_par_cdset( datasets(k), chunk, start )
      k = k+1

      ! Energy
      do i = 1, num_par
        u2 = spec%p( 1, idx(i) )**2 + spec%p( 2, idx(i) )**2 + spec%p( 3, idx(i) )**2
        write_buffer(i) = real(u2 / ( sqrt(1 + u2) + 1) , p_diag_prec)
      enddo
      call diagFile % write_par_cdset( datasets(k), chunk, start )
      k = k+1

#ifdef __HAS_SPIN__

      ! Spin
      do j = 1, p_s_dim
        do i = 1, num_par
          write_buffer(i) = real(spec%s(j,idx(i)), p_diag_prec)
        enddo
        call diagFile % write_par_cdset( datasets(k), chunk, start )
        k = k+1
      enddo

#endif

      call freemem( write_buffer )

      ! Tags
      if ( spec%add_tag ) then
        call alloc( write_tag_buffer, [2, num_par] )
        chunk_tag % count(1) = 2
        chunk_tag % count(2) = num_par
        chunk_tag % start(1) = 0
        chunk_tag % start(2) = start
        chunk_tag % stride(1) = 1
        chunk_tag % stride(2) = 1
        chunk_tag % data = c_loc(write_tag_buffer)

        do i = 1, num_par
          write_tag_buffer(1,i) = spec%tag(1,idx(i))
          write_tag_buffer(2,i) = spec%tag(2,idx(i))
        enddo

        call diagFile % write_par_cdset( datasets(k), chunk_tag, 2*start )

        call freemem( write_tag_buffer )

      endif

      call close_datasets( diagFile, datasets )

      ! close file
      call diagFile % close( )
      deallocate( diagFile )
    endif
  else
    ! No particles, root node writes an empty file
    if ( no_co % local_rank() == 0 ) then
      call write_raw_create_file( spec, 0, n, t, n_aux, diagFile, MPI_COMM_NULL ) 
      call create_datasets( diagFile, num_par, datasets )
      call close_datasets( diagFile, datasets )
      call diagFile % close( )
      deallocate( diagFile )
    endif
  endif

  ! Free particle indexes
  if ( spec%num_par > 0 ) then
    call freemem( idx )
  endif

  if ( comm /= MPI_COMM_NULL .and. comm /= MPI_COMM_SELF ) then
    call MPI_COMM_FREE( comm, ierr )
  endif

end subroutine write_raw_
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
! extract raw data for the boosted diagnostics
!-------------------------------------------------------------------------------
subroutine extract_boosted_raw_( spec, bd, tstep, t1, grid, g_space )

  use m_units

  implicit none

  class( t_species ), intent(inout) :: spec
  class( t_boosted_diag ), intent(inout) :: bd
  type(t_time_step), intent(in) :: tstep
  real(p_double), intent(in) :: t1
  class( t_grid ), intent(in) :: grid
  type( t_space ),     intent(in) :: g_space

  real(p_k_part), dimension(p_cache_size) :: x1, rgamma
  integer, dimension(p_cache_size) :: idx_sub
  integer, dimension(:), pointer :: idx
  real(p_k_part), dimension(spec%get_n_x_dims()) :: x

  integer :: i1, j, k, raw_idx, l, lbuf, i, ibuf, n_x_dims
  integer :: t2_idx, x1_idx_min, x1_idx_max
  real(p_double) :: xmin2, xmax2, rat_dx1
  real(p_double) :: t1_m, x1_m, t2
  real(p_double) :: gamma, gamma_times_beta, beta

  integer :: ptrcur, np, pp
  real(p_k_part) :: gamma_limit_sqr, gamma_par_sqr, v1, u2, x2, rnd, gamma_p
  real(p_k_part) :: px_p
  logical :: has_raw_math_expr
  real(p_k_fparse) :: raw_eval
  real(p_k_fparse), dimension(p_p_dim + spec%get_n_x_dims() + 2) :: raw_var

  ! Set up some parameters
  gamma = bd%grid_mgr%gamma
  gamma_times_beta = bd%grid_mgr%gamma_times_beta
  beta  = bd%grid_mgr%beta
  n_x_dims = spec%get_n_x_dims()
  rat_dx1 = bd%grid_mgr%dx1 / bd%grid_mgr%dx2
  has_raw_math_expr = (bd%raw_math_expr /= '')

  ! initialize gamma limit
  gamma_limit_sqr = spec%diag%raw_gamma_limit**2

  ! i1 has the property t = (i1 - 1) * dt1 + tmin
  i1 = n( tstep ) + 1

  ! Resize idx buffer if necessary
  if ( .not. associated(bd%idx) ) then
    call alloc( bd%idx, [spec%num_par_max] )
  else if ( spec%num_par_max > size(bd%idx) ) then
    call freemem( bd%idx )
    call alloc( bd%idx, [spec%num_par_max] )
  endif
  idx => bd%idx

  ! Get indices of all particles that satisfy the gamma limit
  lbuf = 0
  if ( gamma_limit_sqr > 1.0_p_k_part ) then

    ! Need to compute gamma^2 in the lab frame
    do l = 1, spec%num_par
      ! gamma in boosted frame
      gamma_p = sqrt(1.0_p_k_part + spec%p(1,l)**2 + spec%p(2,l)**2 + spec%p(3,l)**2)
      ! p1 in lab frame
      px_p = real( gamma * ( spec%p(1,l) - beta * gamma_p ), p_k_part )
      ! gamma^2 in lab frame
      gamma_par_sqr = 1.0_p_k_part + px_p**2 + spec%p(2,l)**2 + spec%p(3,l)**2

      if ( gamma_par_sqr >= gamma_limit_sqr ) then
        lbuf = lbuf + 1
        idx(lbuf) = l
      endif
    enddo

  else

    ! All particles will be included
    do l = 1, spec%num_par
      lbuf = lbuf + 1
      idx(lbuf) = l
    enddo

  endif

  ! Loop through these particles and put them into the correct diagnostics
  do ptrcur = 1, lbuf, p_cache_size

    ! check if last copy of table and set np
    if( ptrcur + p_cache_size > lbuf ) then
      np = lbuf - ptrcur + 1
    else
      np = p_cache_size
    endif

    raw_idx = bd%grid_mgr%starting_raw_idx_at_t1(i1)

    ! check which particles should be tabled in each boosted diagnostic
    do j = 1, bd % grid_mgr % num_t2_points_at_t1_raw(i1)

      x1_idx_min = bd%grid_mgr%boosted_raw_coords(raw_idx+j-1)%x1_idx_min
      x1_idx_max = bd%grid_mgr%boosted_raw_coords(raw_idx+j-1)%x1_idx_max

      ! loop over particles and collect those that have valid x1 indices
      ibuf = 0
      pp = ptrcur
      do i = 1, np
        if ( spec%ix(1,idx(pp)) >= x1_idx_min .and. &
             spec%ix(1,idx(pp)) <= x1_idx_max ) then
          ibuf = ibuf + 1
          idx_sub(ibuf) = idx(pp)
        endif
        pp = pp + 1
      enddo

      if ( ibuf == 0 ) cycle

      t2_idx = bd%grid_mgr%boosted_raw_coords(raw_idx+j-1)%t2_idx
      t2 = bd%grid_mgr%boosted_raw_coords(raw_idx+j-1)%t2
      xmin2 = bd%grid_mgr%boosted_raw_coords(raw_idx+j-1)%xmin2
      xmax2 = bd%grid_mgr%boosted_raw_coords(raw_idx+j-1)%xmax2
      raw_var(n_x_dims + p_p_dim + 2) = real( t2, p_k_fparse )

      ! get the x1 positions of the valid particles
      call spec % get_position( 1, idx_sub(1:ibuf), ibuf, x1(1:ibuf) )

      ! compute the x2 position of each particle at time t2
      do i = 1, ibuf
        l = idx_sub(i)
        rgamma(i) = 1.0_p_k_part / &
          sqrt( ( ( 1.0_p_k_part + spec%p(1,l)**2 ) + spec%p(2,l)**2 ) + spec%p(3,l)**2 )
        v1 = spec%p(1,l) * rgamma(i)

        ! assuming a straight-line trajectory at constant velocity, this is the time in
        ! the simulation frame at which the diagnostic frame time equals t2
        t1_m = ( t2 + gamma_times_beta * x1(i) - gamma_times_beta * v1 * t1 ) / &
               ( gamma - gamma_times_beta * v1 )

        ! compute the x1 position at this time
        x1_m = x1(i) + v1 * ( t1_m - t1 )

        ! compute the diagnostic-frame position
        x2 = real( gamma * x1_m - gamma_times_beta * t1_m, p_k_part )

        ! check if this position is in bounds for the diagnostic
        if ( x2 > xmin2 .and. x2 <= xmax2 ) then

          raw_eval = 1
          if ( has_raw_math_expr ) then

            ! fill evaluation variables
            call spec % get_position( l, x )
            x(1) = x2
            call bd % push_back( x, spec%p(:,l), rgamma(i), t1_m, t1 )

            gamma_p = sqrt(1.0_p_k_part + spec%p(1,l)**2 + spec%p(2,l)**2 + spec%p(3,l)**2)
            px_p = real( gamma * ( spec%p(1,l) - beta * gamma_p ), p_k_part )
            gamma_par_sqr = 1.0_p_k_part + px_p**2 + spec%p(2,l)**2 + spec%p(3,l)**2

            raw_var(1:n_x_dims) = real( x(1:n_x_dims), p_k_fparse )
            raw_var(n_x_dims+1) = real( px_p, p_k_fparse )
            raw_var(n_x_dims+2:n_x_dims+3) = real( spec%p(2:3,l), p_k_fparse )
            raw_var(n_x_dims+p_p_dim+1) = real( sqrt(gamma_par_sqr), p_k_fparse )

            ! evaluate
            raw_eval = eval( bd%raw_func, raw_var )

          endif

          if ( raw_eval > 0 ) then
            call rng % harvest_real2( rnd, spec%ix(:,l) )
            if ( rnd > bd%raw_fraction ) cycle
          else
            cycle
          endif

          ! grow the buffers if necessary
          if ( bd%num_par == bd%num_par_max ) then
            call bd % grow_buffer()
          endif

          if ( .not. has_raw_math_expr ) then

            ! need to fill in some parameters

            ! compute the correct particle position data
            call spec % get_position( l, x )

            x(1) = x2
            call bd % push_back( x, spec%p(:,l), rgamma(i), t1_m, t1 )

            gamma_p = sqrt(1.0_p_k_part + spec%p(1,l)**2 + spec%p(2,l)**2 + spec%p(3,l)**2)
            px_p = real( gamma * ( spec%p(1,l) - beta * gamma_p ), p_k_part )

          endif

          ! store the particle data
          bd%num_par = bd%num_par + 1
          k = bd%num_par
          bd%buffered_part(1:n_x_dims,k) = x

          u2 = px_p**2 + spec%p(2,l)**2 + spec%p(3,l)**2

          bd%buffered_part(n_x_dims+1,k) = px_p
          bd%buffered_part(n_x_dims+2:n_x_dims+3,k) = spec%p(2:3,l)
          ! We can think of dividing the charge by gamma, then multiplying by dx1 in the
          ! lab frame, which is expanded by a factor of gamma.  We then divide by dx2 so
          ! that we can multiply these charges by dx2 in post-processing to get the total
          ! charge.  (Remember, q in OSIRIS is really charge density.)  The two factors of
          ! gamma cancel, and we are left with just multiplying by the ratio of the
          ! spatial resolutions.
          bd%buffered_part(n_x_dims+4,k) = real( spec%q(l) * rat_dx1, p_k_part )
          bd%buffered_part(n_x_dims+5,k) = u2 / ( sqrt(1 + u2) + 1)

          bd%buffered_diag_idx(k) = t2_idx

          if ( spec%add_tag ) then
            bd%buffered_tags(:,k) = spec%tag(:,l)
          endif

        endif
      enddo

    enddo

  enddo

end subroutine extract_boosted_raw_
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! write boosted raw particle data
!-------------------------------------------------------------------------------
subroutine write_boosted_raw_( spec, bd, no_co, grid, n )

  use m_units

  implicit none

  class( t_species ), intent(inout) :: spec
  class( t_boosted_diag ), intent(inout) :: bd
  class( t_node_conf ), intent(in) :: no_co
  class( t_grid ), intent(in) :: grid
  integer, intent(in) :: n

  class( t_diag_file ), allocatable :: diagFile
  type( t_diag_dataset ), dimension( spec%get_n_x_dims() + p_p_dim + 3 ) :: datasets
  type( t_diag_chunk ) :: chunk, chunk_tag
  integer :: start, rank, num_par, total_par, my_par, color, comm, ierr, num_par_saved
  integer :: ndims, i, j, k, nt2, i2, i0, data_type
  integer, dimension(:), pointer :: idx, n_file, idx_sorted
  real( p_diag_prec ), dimension(:), pointer :: write_buffer
  integer, dimension( :, : ), pointer :: write_tag_buffer

  if ( n == 0 ) then
    call bd%create_files_raw( grid )
  endif

  ndims = spec%get_n_x_dims()
  nt2 = bd%grid_mgr%nt2
  num_par = bd%num_par

  if ( no_co % comm_size() == 1 .and. num_par == 0 ) return

  ! get a list of sorted particles and the number in each output file
  if ( num_par > 0 ) then
    call alloc( idx, (/ bd%num_par /) )
    call alloc( idx_sorted, (/ bd%num_par /) )
    call alloc( n_file, (/ nt2 /) )
    call generate_sort_idx( bd, idx, n_file, idx_sorted )
  endif

  ! loop through the diagnostics
  do i2 = 1, nt2
    if ( no_co % comm_size() > 1 ) then
      comm = MPI_COMM_NULL

      if ( num_par > 0 ) then
        if ( i2 == 1 ) then
          my_par = n_file(1)
        else
          my_par = n_file(i2) - n_file(i2-1)
        endif
      else
        my_par = 0
      endif

      call MPI_ALLREDUCE( my_par, total_par, 1, MPI_INTEGER, MPI_SUM, no_co % comm, ierr )

      if ( total_par > 0 ) then
        if ( my_par > 0 ) then
          color = 1
        else
          color = MPI_UNDEFINED
        endif
        call MPI_COMM_SPLIT( no_co % comm, color, 0, comm, ierr )

        if ( my_par > 0 ) then
          call MPI_COMM_RANK( comm, rank, ierr )
          call MPI_EXSCAN( my_par, start, 1, MPI_INTEGER, MPI_SUM, comm, ierr )
          if ( rank == 0 ) start = 0
        endif
      endif
    else
      comm = MPI_COMM_SELF
      if ( i2 == 1 ) then
        my_par = n_file(1)
      else
        my_par = n_file(i2) - n_file(i2-1)
      endif
      total_par = my_par
      start = 0
      rank = 0
    endif

    ! Write data, if any
    if ( total_par > 0 ) then

      num_par_saved = bd%num_par_diag(i2)

      ! Only nodes with particles are involved in this section
      if ( my_par > 0 ) then

        ! Create the file
        call create_diag_file( diagFile )
        call bd%populate_raw_diagFile( diagFile, i2 - 1 )
        call diagFile % open( p_diag_update, comm )

        ! Open the datasets
        select case (p_diag_prec)
        case ( p_single )
          data_type = p_diag_float32
        case ( p_double )
          data_type = p_diag_float64
        end select

        if ( diagFile%particles%has_tags ) then
          j = diagFile % particles % nquants - 1
        else
          j = diagFile % particles % nquants
        endif

        ! particle quantities
        do i = 1, j
          datasets(i)%name = diagFile % particles % quants(i)
          datasets(i)%ndims = 1
          datasets(i)%data_type = data_type
          call diagFile % open_cdset( datasets(i) )
          if ( rank == 0 ) then
            call diagFile % extend_cdset( datasets(i), int([num_par_saved + total_par], p_int64) )
          endif
        enddo

        ! particle tags
        if ( diagFile%particles%has_tags ) then
          i = diagFile % particles % nquants
          datasets(i)%name = 'tag'
          datasets(i)%ndims = 2
          datasets(i)%data_type = p_diag_int32
          call diagFile % open_cdset( datasets(i) )
          if ( rank == 0 ) then
            call diagFile % extend_cdset( datasets(i), int([2, num_par_saved + total_par], p_int64) )
          endif
        endif

        call alloc( write_buffer, [ my_par ] )
        chunk % start(1) = start + num_par_saved
        chunk % count(1) = my_par
        chunk % stride(1) = 1
        chunk % data = c_loc( write_buffer )

        ! Prepare to index particle list
        if ( i2 == 1 ) then
          i0 = 0
        else
          i0 = n_file(i2-1)
        endif

        ! Positions
        k = 1
        do j = 1, ndims
          do i = 1, my_par
            write_buffer(i) = real( bd%buffered_part(j,idx_sorted(i+i0)), &
                                    p_diag_prec )
          enddo
          call diagFile % write_par_cdset( datasets(k), chunk, start )
          k = k+1
        enddo

        ! Momenta
        do j = 1, p_p_dim
          do i = 1, my_par
            write_buffer(i) = real(bd%buffered_part(j+ndims,idx_sorted(i+i0)), &
                                   p_diag_prec)
          enddo
          call diagFile % write_par_cdset( datasets(k), chunk, start )
          k = k+1
        enddo

        ! Charge
        do i = 1, my_par
          write_buffer(i) = real(bd%buffered_part(p_p_dim+ndims+1,idx_sorted(i+i0)), &
                                 p_diag_prec)
        enddo
        call diagFile % write_par_cdset( datasets(k), chunk, start )
        k = k+1

        ! Energy
        do i = 1, my_par
          write_buffer(i) = real(bd%buffered_part(p_p_dim+ndims+2,idx_sorted(i+i0)), &
                                 p_diag_prec)
        enddo
        call diagFile % write_par_cdset( datasets(k), chunk, start )
        k = k+1

        call freemem( write_buffer )

        ! Tags
        if ( spec%add_tag ) then
          call alloc( write_tag_buffer, [2, my_par] )
          chunk_tag % count(1) = 2
          chunk_tag % count(2) = my_par
          chunk_tag % start(1) = 0
          chunk_tag % start(2) = start + num_par_saved
          chunk_tag % stride(1) = 1
          chunk_tag % stride(2) = 1
          chunk_tag % data = c_loc(write_tag_buffer)

          do i = 1, my_par
            write_tag_buffer(1,i) = bd%buffered_tags(1,idx_sorted(i+i0))
            write_tag_buffer(2,i) = bd%buffered_tags(2,idx_sorted(i+i0))
          enddo

          call diagFile % write_par_cdset( datasets(k), chunk_tag, 2*start )

          call freemem( write_tag_buffer )

        endif

        call close_datasets( diagFile, datasets )

        ! close file
        call diagFile % close( )
        deallocate( diagFile )
      endif

      bd%num_par_diag(i2) = num_par_saved + total_par

    endif

  enddo

  ! Free particle indexes
  if ( num_par > 0 ) then
    call freemem( idx )
    call freemem( idx_sorted )
    call freemem( n_file )
  endif

  ! Remove particles from boosted diagnostic
  bd%num_par = 0

  if ( comm /= MPI_COMM_NULL .and. comm /= MPI_COMM_SELF ) then
    call MPI_COMM_FREE( comm, ierr )
  endif

end subroutine write_boosted_raw_
!-------------------------------------------------------------------------------

end module m_species_rawdiag
