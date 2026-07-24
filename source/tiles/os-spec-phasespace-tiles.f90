#include "os-preprocess.fpp"
#include "os-config.h"

module m_species_phasespace_tiles

#include "memory/memory.h"

use m_parameters
use m_species_define, only : t_species, t_spec_arr, t_phasespace, t_phasespace_diagnostics
use m_species_define, only : t_phasespace_params
use m_species_phasespace
use m_space, only : t_space
use m_node_conf, only : t_node_conf, root, comm, no_num, p_sum, reduce_array_size
use m_diagnostic_utilities
use m_time_step, only : t_time_step, ndump, n, test_if_report

private

interface report
  module procedure report_phasespace_tiles
end interface

public :: report

contains

!---------------------------------------------------------------------------------------------------
subroutine report_phasespace_tiles( this, spec_arr, no_co, g_space, tstep, t )
!---------------------------------------------------------------------------------------------------

  implicit none

  class( t_phasespace_diagnostics ), intent(inout) :: this
  type( t_spec_arr ), dimension(:), pointer :: spec_arr
  class( t_node_conf ), intent(in) :: no_co
  type( t_space ), intent(in) :: g_space
  type( t_time_step ), intent(in) :: tstep
  real( p_double ), intent(in) :: t

  type( t_phasespace ), pointer :: phasespace

  if (test_if_report( tstep, this%ndump_fac )) then

    ! autorange phasespaces ( p, gamma, l )
    call acquire_range_tiles( spec_arr, no_co )

    ! Process phasespace list
    phasespace => this%phasespace_list%head
    do
      if (.not. associated(phasespace)) exit
      call write_phasespace_tiles( phasespace, spec_arr, no_co, g_space, tstep, t )

      phasespace => phasespace%next
    enddo

    ! pha_ene_bin, pha_cell_avg, and pha_time_avg lists are not currently implemented

  endif

end subroutine report_phasespace_tiles

!-----------------------------------------------------------------------------------------
 subroutine acquire_range_tiles( spec, no_co )
!-----------------------------------------------------------------------------------------
! acquires the ranges for autorange of phasespaces for tiled runs
!-----------------------------------------------------------------------------------------

  implicit none

  integer, parameter :: max_quant = 7 ! 3 linear momentum, 3 angular momentum, 1 gamma

  type(t_spec_arr), dimension(:), pointer :: spec
  class(t_node_conf), intent(in) :: no_co

  class(t_species), pointer :: species
  integer :: i, l, j, t, num_til, ierr, np, nquant, quant_idx
  real(p_k_part), dimension(2, max_quant) :: range_local, range_global
  real(p_k_part) :: dx, gamma
  integer, parameter :: p_ang_p_npar = 128
  real(p_diag_prec), dimension(p_ang_p_npar) :: ang_p

  ! First do those things which can be done once for all tiles
  species => spec(1)%s
  num_til = size(spec)

  ! check wich ranges to get
  nquant = 0

  ! check momenta
  do i=1, p_p_dim
    if (species%diag%phasespaces%if_p_auto(i)) nquant = nquant + 1
  enddo

  ! check gamma
  if (species%diag%phasespaces%if_gamma_auto) nquant = nquant + 1

  ! check angular momenta
  do i=1, 3
    if (species%diag%phasespaces%if_l_auto(i)) nquant = nquant + 1
  enddo

  ! If any ranges necessary...
  if ( nquant > 0) then

    ! This insures that any value will be less than the set minimum and

    ! more than the set maximum

    range_local( 1, 1:nquant ) = +huge( 1.0_p_k_part )
    range_local( 2, 1:nquant ) = -huge( 1.0_p_k_part )

    quant_idx = 1

    ! find local momenta limits
    do i=1, p_p_dim
      if (species%diag%phasespaces%if_p_auto(i)) then

        ! Loop over all tiles
        do t = 1, num_til
          species => spec(t)%s
          do l = 1, species%num_par
            if ( species%p(i,l) < range_local( 1, quant_idx ) ) then
              range_local( 1, quant_idx ) =  species%p(i,l)
            else if ( species%p(i,l) > range_local( 2, quant_idx ) ) then
              range_local( 2, quant_idx ) =  species%p(i,l)
            endif
          enddo
        enddo
        species => spec(1)%s

        quant_idx = quant_idx + 1
      endif
    enddo

    ! find local gamma limits
    if (species%diag%phasespaces%if_gamma_auto) then

      ! Loop over all tiles
      do t = 1, num_til
        species => spec(t)%s
        do l = 1, species%num_par
          gamma = sqrt( species%p(1,l)**2 + species%p(2,l)**2 + species%p(3,l)**2 + 1 )
          if ( gamma < range_local( 1, quant_idx ) ) then
            range_local( 1, quant_idx ) =  gamma
          else if ( gamma > range_local( 2, quant_idx ) ) then
            range_local( 2, quant_idx ) =  gamma
          endif
        enddo
      enddo
      species => spec(1)%s

      quant_idx = quant_idx + 1
    endif

    ! find local angular momenta limits
    do i=1, 3
      if (species%diag%phasespaces%if_l_auto(i)) then

        ! particles are processed in p_ang_p_npar batches

        ! Loop over all tiles
        do t = 1, num_til
          species => spec(t)%s

          ! process first batch
          np = min ( p_ang_p_npar, species%num_par )
          do j = 1, np
            if ( ang_p(j) < range_local( 1, quant_idx ) ) then
              range_local( 1, quant_idx ) =  ang_p(j)
            else if ( ang_p(j) > range_local( 2, quant_idx ) ) then
              range_local( 2, quant_idx ) =  ang_p(j)
            endif
          enddo

          ! process remaining batches
          do l = p_ang_p_npar+1, species%num_par, p_ang_p_npar
            np = min ( l + p_ang_p_npar - 1, species%num_par )
            call get_l( ang_p, i, species, l, np )
            do j = 1, np-l+1
              if ( ang_p(j) < range_local( 1, quant_idx ) ) then
                range_local( 1, quant_idx ) =  ang_p(j)
              else if ( ang_p(j) > range_local( 2, quant_idx ) ) then
                range_local( 2, quant_idx ) =  ang_p(j)
              endif
            enddo
          enddo

        enddo
        species => spec(1)%s

        quant_idx = quant_idx + 1
      endif

    enddo

    ! If no particles are present on the node the autorange values would now be set to
    ! the values defined in the input file

    ! Flip sign of minimum values to use a single MPI_REDUCE operation
    range_local( 1, : ) = - range_local( 1, : )

    ! find global ranges
    if ( no_num( no_co ) > 1 ) then
      call MPI_ALLREDUCE(range_local, range_global, 2 * nquant, mpi_real_type( p_k_part ), &
        MPI_MAX, comm( no_co ), ierr)
    else
      range_global( :, 1:nquant) = range_local( :, 1:nquant )
    endif

    ! check if at least 1 node had particles in it (as a result the maximum range of quantity
    ! 1 would have been set)
    if ( range_global( 2, 1 ) /= -huge( 1.0_p_k_part ) ) then

      ! Flip sign of minimum values back to get correct values
      range_global( 1, : ) = - range_global( 1, : )

      ! store results. if using autorange values try to place the highest valued particle in
      ! the middle of the cell

      ! Not necessary to loop over tiles here since first tile contains information
      quant_idx = 1

      ! momenta

      do i=1, p_p_dim
        if (species%diag%phasespaces%if_p_auto(i)) then

          dx = (range_global(2, quant_idx) - range_global(1, quant_idx)) / species%diag%phasespaces%np(i)
          range_global( 1, quant_idx ) = range_global( 1, quant_idx ) - dx/2
          range_global( 2, quant_idx ) = range_global( 2, quant_idx ) + dx/2

          if ( species%diag%phasespaces%pmin(i) > range_global( 1, quant_idx ) ) then
            species%diag%phasespaces%pmin(i) = range_global( 1, quant_idx )
          endif

          if ( species%diag%phasespaces%pmax(i) < range_global( 2, quant_idx ) ) then
            species%diag%phasespaces%pmax(i) = range_global( 2, quant_idx )
          endif

          quant_idx = quant_idx + 1
        endif
      enddo

      ! gamma
      if (species%diag%phasespaces%if_gamma_auto) then

        dx = (range_global(2, quant_idx) - range_global(1, quant_idx)) / species%diag%phasespaces%ngamma
        range_global( 1, quant_idx ) = range_global( 1, quant_idx ) - dx/2
        range_global( 2, quant_idx ) = range_global( 2, quant_idx ) + dx/2

        if ( species%diag%phasespaces%gammamin > range_global( 1, quant_idx ) ) then
          species%diag%phasespaces%gammamin = range_global( 1, quant_idx )
        endif

        if ( species%diag%phasespaces%gammamax < range_global( 2, quant_idx ) ) then
          species%diag%phasespaces%gammamax = range_global( 2, quant_idx )
        endif

        quant_idx = quant_idx + 1
      endif

      ! ang. momenta

      do i=1, 3
        if (species%diag%phasespaces%if_l_auto(i)) then
          dx = (range_global(2, quant_idx) - range_global(1, quant_idx)) / species%diag%phasespaces%nl(i)
          range_global( 1, quant_idx ) = range_global( 1, quant_idx ) - dx/2
          range_global( 2, quant_idx ) = range_global( 2, quant_idx ) + dx/2

          if ( species%diag%phasespaces%lmin(i) > range_global( 1, quant_idx ) ) then
            species%diag%phasespaces%lmin(i) = range_global( 1, quant_idx )
          endif

          if ( species%diag%phasespaces%lmax(i) < range_global( 2, quant_idx ) ) then
            species%diag%phasespaces%lmax(i) = range_global( 2, quant_idx )
          endif

          quant_idx = quant_idx + 1
        endif
      enddo

    endif

  endif

 end subroutine acquire_range_tiles
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Normal phasespaces for tiled runs
!-----------------------------------------------------------------------------------------
subroutine write_phasespace_tiles( phasespace, spec_arr, no_co, g_space, tstep, t )
!-----------------------------------------------------------------------------------------

  implicit none

  type(t_phasespace), intent(in) :: phasespace
  type(t_spec_arr), dimension(:), pointer :: spec_arr
  class( t_node_conf ), intent(in) :: no_co
  type( t_space ), intent(in) :: g_space
  type( t_time_step ), intent(in) :: tstep
  real(p_double), intent(in) :: t

  class(t_species), pointer :: species
  type(t_phasespace_params) :: params
  real(p_diag_prec), dimension(:), pointer :: buffer
  integer :: bsize, chunksize, nchunks, data_type
  class( t_diag_file ), allocatable :: diagFile
  type( t_diag_dataset ) :: dset

  species => spec_arr(1)%s

  ! select the right parameters for this phasespace
  call species%diag%phasespaces%get_params( phasespace, g_space, params )

  ! Open the file

  if ( root( no_co ) ) then
    ! create the file
    call init_phasespace_file( params, species, n(tstep)/ndump(tstep), n(tstep), t, diagFile )
    call diagFile % open( p_diag_create )

    ! create dataset
    select case (p_diag_prec)
    case (p_double)
      data_type = p_diag_float64
    case (p_single)
      data_type = p_diag_float32
    end select

    call diagFile % start_cdset( params % name, params % ndims, params % n, data_type, dset)

  endif

  ! get buffer for diagnostics
  call get_diag_buffer( buffer, bsize )

  ! get division of phasespace
  ! division must be on the last coordinate
  if ( species%diag%phasespaces%size( phasespace ) > bsize ) then
    ! multiple chunks
    nchunks = ceiling( float( species%diag%phasespaces%size( phasespace ) ) / bsize )
    if ( nchunks > params%n( params%ndims ) ) then
      if ( root(no_co) ) then
        write(0,'(A,A)') '(*error*) phasespace too large: ', trim(params%name)
        write(0,'(A)')   '(*error*) try smaller phasespace or increasing diag. buffer size'
      endif
      call abort_program( p_err_alloc )
    endif
    chunksize = params%n( params%ndims ) / nchunks
  else
    ! single chunk
    chunksize = params%n( params%ndims )
  endif

  ! process phasespace
  select case(phasespace%ndims)
  case(1)
    call save_phasespace_1D_tiles( phasespace, params, spec_arr, no_co, &
                                    chunkSize, buffer, diagFile, dset )

  case(2)
    call save_phasespace_2D_tiles( phasespace, params, spec_arr, no_co, &
                                    chunkSize, buffer, diagFile, dset )

  case(3)
    call save_phasespace_3D_tiles( phasespace, params, spec_arr, no_co, &
                                    chunkSize, buffer, diagFile, dset )

  end select

  if ( root( no_co ) ) then

    ! close dataset on file
    call diagFile % end_cdset( dset )

    ! close the file
    call diagFile % close( )

    deallocate( diagFile )

  endif

end subroutine write_phasespace_tiles
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine save_phasespace_1D_tiles( phasespace, params, spec_arr, no_co, &
                                      chunkSize, buffer, diagFile, dset, ene_bin )
!-----------------------------------------------------------------------------------------

  implicit none

  integer, parameter :: rank = 1

  type(t_phasespace), intent(in) :: phasespace    ! phasespace to save
  type(t_phasespace_params), intent(in) :: params  ! parsed parameters for phasespace
  type(t_spec_arr), dimension(:), pointer :: spec_arr  ! species object array
  class( t_node_conf ), intent(in) :: no_co      ! parallel node configuration
  integer, intent(in) :: chunkSize          ! max. chunk size
  real(p_diag_prec), dimension(:), pointer :: buffer  ! buffer for diagnostic
  class( t_diag_file ), allocatable, intent(inout) :: diagFile
  type( t_diag_dataset ), intent(inout) :: dset
  real(p_diag_prec), dimension(:), optional, intent(in) :: ene_bin

  class(t_species), pointer :: species
  type( t_diag_chunk ) :: chunk
  integer :: lchunk, uchunk    ! lower/upper boundary for chunk
  integer :: bsize        ! buffer size in use (in array elements)
  integer :: l, lp, np      ! particle index, lower, number of

  ! arrays to hold temp particle data
  real(p_diag_prec), dimension(p_par_buf_size) :: xp1, charge
  integer :: num_til, i

  num_til = size(spec_arr)

  ! loop through all chunks
  lchunk = 1
  do

    if ( lchunk > params%n(rank) ) exit

    uchunk = min( lchunk + chunkSize - 1, params%n(rank) )

    bsize = ( uchunk - lchunk + 1 )

    buffer( 1 : bsize ) = 0

    ! loop through all tiles
    do i = 1, num_til
      species => spec_arr(i)%s

      ! loop through all particles
      if ( species%num_par > 0 ) then

        l = 1
        do

          if ( l > species%num_par ) exit

          lp = min(l + p_par_buf_size - 1, species%num_par)
          np = lp - l + 1

          call species%get_phasespace_axis( xp1, l, lp, phasespace%x_or_p(1), phasespace%xp_dim(1))

          if (present(ene_bin)) then
            call get_phasespace_charge( phasespace%ps_type, charge, species, l, lp, ene_bin )
          else
            call get_phasespace_charge( phasespace%ps_type, charge, species, l, lp )
          endif

          ! deposit on slab
          call deposit_1D( buffer, params%n, &
            xp1, charge, np, &
            params%min, params%max, &
            (/ lchunk /), (/ uchunk /) )

          l = l + p_par_buf_size
        enddo
      endif
    enddo
    species => spec_arr(1)%s

    ! accumulate data on node 0
    call reduce_array_size( no_co, buffer, bsize, p_sum )

    ! write chunk
    if ( root( no_co ) ) then
      ! normalize data
      call normalize_1D( buffer, (/ lchunk /), (/ uchunk /), phasespace, params, &
                          species%dx, species%coordinates )

      ! chunk start coordinates are 0 indexed
      chunk % start( 1 ) = 0
      chunk % count( 1 ) = uchunk - lchunk
      chunk % stride( 1 ) = 1

      chunk % data = c_loc( buffer )

      call diagFile % write_cdset( dset, chunk )

    endif

    ! Process next chunk
    lchunk = lchunk + chunkSize
  enddo

end subroutine save_phasespace_1D_tiles
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine save_phasespace_2D_tiles( phasespace, params, spec_arr, no_co, &
                                      chunkSize, buffer, diagFile, dset, ene_bin )
!-----------------------------------------------------------------------------------------

  implicit none

  integer, parameter :: rank = 2

  type(t_phasespace), intent(in) :: phasespace    ! phasespace to save
  type(t_phasespace_params), intent(in) :: params  ! parsed parameters for phasespace
  type(t_spec_arr), dimension(:), pointer :: spec_arr  ! species object array
  class( t_node_conf ), intent(in) :: no_co      ! parallel node configuration
  integer, intent(in) :: chunkSize          ! max. chunk size
  real(p_diag_prec), dimension(:), pointer :: buffer  ! buffer for diagnostic
  class( t_diag_file ), allocatable, intent(inout) :: diagFile
  type( t_diag_dataset ), intent(inout) :: dset
  real(p_diag_prec), dimension(:), optional, intent(in) :: ene_bin

  class(t_species), pointer :: species
  type( t_diag_chunk ) :: chunk
  integer :: lchunk, uchunk    ! lower/upper boundary for chunk
  integer :: bsize        ! buffer size in use (in array elements)
  integer :: l, lp, np      ! particle index, lower, number of

  ! arrays to hold temp particle data
  real(p_diag_prec), dimension(p_par_buf_size) :: xp1, xp2, charge
  integer :: num_til, i

  num_til = size(spec_arr)

  ! loop through all chunks
  lchunk = 1
  do

    if ( lchunk > params%n(rank) ) exit

    uchunk = min( lchunk + chunkSize - 1, params%n(rank) )

    bsize = params%n(1) * ( uchunk - lchunk + 1 )

    buffer( 1 : bsize ) = 0

    ! loop through all tiles
    do i = 1, num_til
      species => spec_arr(i)%s

      ! loop through all particles
      if ( species%num_par > 0 ) then

        l = 1
        do

          if ( l > species%num_par ) exit

          lp = min(l + p_par_buf_size - 1, species%num_par)
          np = lp - l + 1

          call species%get_phasespace_axis( xp1, l, lp, phasespace%x_or_p(1), phasespace%xp_dim(1))
          call species%get_phasespace_axis( xp2, l, lp, phasespace%x_or_p(2), phasespace%xp_dim(2))

          if (present(ene_bin)) then
            call get_phasespace_charge( phasespace%ps_type, charge, species, l, lp, ene_bin )
          else
            call get_phasespace_charge( phasespace%ps_type, charge, species, l, lp )
          endif

          ! deposit on slab
          call deposit_2D( buffer, params%n, &
            xp1, xp2, charge, np, &
            params%min, params%max, &
            (/ 1, lchunk /), (/ params%n(1), uchunk /) )

          l = l + p_par_buf_size
        enddo
      endif
    enddo
    species => spec_arr(1)%s

    ! accumulate data on node 0
    call reduce_array_size( no_co, buffer, bsize, p_sum )

    ! write chunk
    if ( root( no_co ) ) then
      ! normalize data
      call normalize_2D( buffer, (/ 1, lchunk /), (/ params%n(1), uchunk /), &
                          phasespace, params, &
                          species%dx, species%coordinates )

      ! chunk start coordinates are 0 indexed
      chunk % start( 1:2 ) = [0,lchunk-1]
      chunk % count( 1:2 ) = [params%n(1), uchunk - lchunk + 1]
      chunk % stride( 1:2 ) = 1

      chunk % data = c_loc( buffer )

      call diagFile % write_cdset( dset, chunk )

    endif

    ! Process next chunk
    lchunk = lchunk + chunkSize
  enddo

end subroutine save_phasespace_2D_tiles
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine save_phasespace_3D_tiles( phasespace, params, spec_arr, no_co, &
                                      chunkSize, buffer, diagFile, dset, ene_bin )
!-----------------------------------------------------------------------------------------

  implicit none

  integer, parameter :: rank = 3

  type(t_phasespace), intent(in) :: phasespace    ! phasespace to save
  type(t_phasespace_params), intent(in) :: params  ! parsed parameters for phasespace
  type(t_spec_arr), dimension(:), pointer :: spec_arr  ! species object array
  class( t_node_conf ), intent(in) :: no_co      ! parallel node configuration
  integer, intent(in) :: chunkSize          ! max. chunk size
  real(p_diag_prec), dimension(:), pointer :: buffer  ! buffer for diagnostic
  class( t_diag_file ), allocatable, intent(inout) :: diagFile
  type( t_diag_dataset ), intent(inout) :: dset
  real(p_diag_prec), dimension(:), optional, intent(in) :: ene_bin

  class(t_species), pointer :: species
  type( t_diag_chunk ) :: chunk
  integer :: lchunk, uchunk    ! lower/upper boundary for chunk
  integer :: bsize        ! buffer size in use (in array elements)
  integer :: l, lp, np      ! particle index, lower, number of

  ! arrays to hold temp particle data
  real(p_diag_prec), dimension(p_par_buf_size) :: xp1, xp2, xp3, charge
  integer :: num_til, i

  num_til = size(spec_arr)

  ! loop through all chunks
  lchunk = 1
  do

    if ( lchunk > params%n(rank) ) exit

    uchunk = min( lchunk + chunkSize - 1, params%n(rank) )

    bsize = params%n(1) * params%n(2) * ( uchunk - lchunk + 1 )

    buffer( 1 : bsize ) = 0

    ! loop through all tiles
    do i = 1, num_til
      species => spec_arr(i)%s

      ! loop through all particles
      if ( species%num_par > 0 ) then

        l = 1
        do

          if ( l > species%num_par ) exit

          lp = min(l + p_par_buf_size - 1, species%num_par)
          np = lp - l + 1

          call species%get_phasespace_axis( xp1, l, lp, phasespace%x_or_p(1), phasespace%xp_dim(1))
          call species%get_phasespace_axis( xp2, l, lp, phasespace%x_or_p(2), phasespace%xp_dim(2))
          call species%get_phasespace_axis( xp3, l, lp, phasespace%x_or_p(3), phasespace%xp_dim(3))

          if (present(ene_bin)) then
            call get_phasespace_charge( phasespace%ps_type, charge, species, l, lp, ene_bin )
          else
            call get_phasespace_charge( phasespace%ps_type, charge, species, l, lp )
          endif

          ! deposit on slab
          call deposit_3D( buffer, params%n, &
            xp1, xp2, xp3, charge, np, &
            params%min, params%max, &
            (/ 1, 1, lchunk /), (/ params%n(1), params%n(2), uchunk /) )

          l = l + p_par_buf_size
        enddo
      endif
    enddo
    species => spec_arr(1)%s

    ! accumulate data on node 0
    call reduce_array_size( no_co, buffer, bsize, p_sum )

    ! write chunk
    if ( root( no_co ) ) then
      ! normalize data
      call normalize_3D( buffer, (/ 1, 1, lchunk /), (/ params%n(1), params%n(2), uchunk /), &
                          phasespace, params, &
                          species%dx, species%coordinates )

      ! chunk start coordinates are 0 indexed
      chunk % start( 1:3 ) = [0, 0, lchunk-1]
      chunk % count( 1:3 ) = [params%n(1), params%n(2),uchunk - lchunk + 1]
      chunk % stride( 1:3 ) = 1

      chunk % data = c_loc( buffer )

      call diagFile % write_cdset( dset, chunk )

    endif

    ! Process next chunk
    lchunk = lchunk + chunkSize
  enddo

end subroutine save_phasespace_3D_tiles
!-----------------------------------------------------------------------------------------

end module m_species_phasespace_tiles
