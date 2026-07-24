# 1 "emf/solver/os-emf-solver-ck.f03"
# 1 "<built-in>" 1
# 1 "<built-in>" 3
# 467 "<built-in>" 3
# 1 "<command line>" 1
# 1 "<built-in>" 2
# 1 "emf/solver/os-emf-solver-ck.f03" 2
!-----------------------------------------------------------------------------------------
! Cole-Karkkainen Solver
! - M. Karkkainen, et.al., "Low-Dispersion Wake Field Calculation Tools", Proceedings of
! ICAP 2006, Chamonix, France, vol. 2, pp. 35-40, Jan. 2007.
! - J. B. Cole, "A high-accuracy Yee algorithm based on nonstandard finite differences:
! new developments and verifications", IEEE Trans. Antennas Prop., vol. 50, no. 9,
! Sept 2002, pp. 1185–1191
!-----------------------------------------------------------------------------------------

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
# 11 "emf/solver/os-emf-solver-ck.f03" 2
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
# 12 "emf/solver/os-emf-solver-ck.f03" 2

module m_emf_solver_ck

use m_parameters
use m_emf_define
use m_vdf_define
use m_input_file

private

type, extends( t_emf_solver ), public :: t_emf_solver_ck

    procedure(dedt_ck), pointer :: dedt => null()
    procedure(dbdt_ck), pointer :: dbdt => null()

    ! Solver parameters
    real(p_double) :: k1, k2, k3

contains

    procedure :: init => init_ck
    procedure :: advance => advance_ck
    procedure :: min_gc => min_gc_ck
    procedure :: name => name_ck
    procedure :: read_input => read_input_ck
    procedure :: test_stability => test_stability_ck


end type t_emf_solver_ck

interface
    subroutine dedt_ck( this, e, b, jay, dt )
        import t_emf_solver_ck, t_vdf, p_double
        class( t_emf_solver_ck ), intent(in) :: this
        type( t_vdf ), intent(inout) :: e
        type( t_vdf ), intent(in) :: b, jay
        real(p_double), intent(in) :: dt
    end subroutine
end interface

interface
    subroutine dbdt_ck( this, b, e, dt )
        import t_emf_solver_ck, t_vdf, p_double
        class( t_emf_solver_ck ), intent(in) :: this
        type( t_vdf ), intent(inout) :: b
        type( t_vdf ), intent(in) :: e
        real(p_double), intent(in) :: dt
    end subroutine
end interface

contains

function name_ck( this )

    implicit none

    class( t_emf_solver_ck ), intent(in) :: this
    character( len = p_maxlen_emf_solver_name ) :: name_ck

    name_ck = "Cole-Karkkainen"

end function name_ck


subroutine dbdt_2d_ck( this, b, e, dt )

    implicit none

    class( t_emf_solver_ck ), intent(in) :: this
    type( t_vdf ), intent( inout ) :: b
    type( t_vdf ), intent(in) :: e
    real(p_double), intent(in) :: dt

    real(p_k_fld) :: dtdx1, dtdx2
    real(p_k_fld) :: k1, k2
    integer :: i1, i2

    dtdx1 = real(dt/b%dx_(1), p_k_fld) !*c*c !c=1
    dtdx2 = real(dt/b%dx_(2), p_k_fld) !*c*c !c=1

    k1 = real( this % k1, p_k_fld )
    k2 = real( this % k2, p_k_fld )

    !$omp parallel do private(i1)
    do i2 = -1, b%nx_(2)+1
        do i1 = -1, b%nx_(1)+1

            !B1
            b%f2( 1, i1, i2 ) = b%f2( 1, i1, i2 ) &
                            - dtdx2 * ( k1 * ( e%f2( 3, i1 , i2+1 ) - e%f2( 3, i1 , i2 )) &
                                    + k2 * ( e%f2( 3, i1+1, i2+1 ) - e%f2( 3, i1+1, i2 ) &
                                            + e%f2( 3, i1-1, i2+1 ) - e%f2( 3, i1-1, i2 )))

            !B2
            b%f2( 2, i1, i2 ) = b%f2( 2, i1, i2 ) &
                            + dtdx1 * ( k1 * ( e%f2( 3, i1+1, i2 ) - e%f2( 3, i1, i2 )) &
                                    + k2 * ( e%f2( 3, i1+1, i2+1 ) - e%f2( 3, i1, i2+1 ) &
                                            + e%f2( 3, i1+1, i2-1 ) - e%f2( 3, i1, i2-1 )))

            !B3
            b%f2( 3, i1, i2 ) = b%f2( 3, i1, i2 ) &
                                - dtdx1 * ( k1 * ( e%f2( 2, i1+1, i2 ) - e%f2( 2, i1, i2 )) &
                                        + k2 * ( e%f2( 2, i1+1, i2+1 ) - e%f2( 2, i1, i2+1 ) &
                                                + e%f2( 2, i1+1, i2-1 ) - e%f2( 2, i1, i2-1 ))) &
                                + dtdx2 * ( k1 * ( e%f2( 1, i1 , i2+1 ) - e%f2( 1, i1 , i2 )) &
                                        + k2 * ( e%f2( 1, i1+1, i2+1 ) - e%f2( 1, i1+1, i2 ) &
                                                + e%f2( 1, i1-1, i2+1 ) - e%f2( 1, i1-1, i2 )))

        enddo
    enddo
    !$omp end parallel do

end subroutine dbdt_2d_ck


subroutine dbdt_3d_ck( this, b, e, dt )

    implicit none

    class( t_emf_solver_ck ), intent(in) :: this
    type( t_vdf ), intent( inout ) :: b
    type( t_vdf ), intent(in) :: e
    real(p_double), intent(in) :: dt

    ! local variables
    real(p_k_fld) :: dtdx1, dtdx2, dtdx3
    integer :: i1, i2, i3

    real(p_k_fld) :: k1, k2, k3

    dtdx1 = real( dt/b%dx_(1), p_k_fld )
    dtdx2 = real( dt/b%dx_(2), p_k_fld )
    dtdx3 = real( dt/b%dx_(3), p_k_fld )

    k1 = real( this % k1, p_k_fld )
    k2 = real( this % k2, p_k_fld )
    k3 = real( this % k3, p_k_fld )

    !$omp parallel do private(i2,i1)
    do i3 = -1, b%nx_(3)+1
        do i2 = -1, b%nx_(2)+1
            do i1 = -1, b%nx_(1)+1


            !B1
            b%f3( 1, i1, i2, i3 ) = b%f3( 1, i1, i2, i3 ) &
                    - dtdx2 * ( k1 * ( e%f3( 3, i1 , i2+1, i3 ) - e%f3( 3, i1 , i2 , i3 )) &
                            + k2 * ( e%f3( 3, i1+1, i2+1, i3 ) - e%f3( 3, i1+1, i2 , i3 ) &
                                    + e%f3( 3, i1-1, i2+1, i3 ) - e%f3( 3, i1-1, i2 , i3 ) &
                                    + e%f3( 3, i1 , i2+1, i3+1 ) - e%f3( 3, i1 , i2 , i3+1 ) &
                                    + e%f3( 3, i1 , i2+1, i3-1 ) - e%f3( 3, i1 , i2 , i3-1 )) &
                            + k3 * ( e%f3( 3, i1+1, i2+1, i3+1 ) - e%f3( 3, i1+1, i2 , i3+1 ) &
                                    + e%f3( 3, i1-1, i2+1, i3-1 ) - e%f3( 3, i1-1, i2 , i3-1 ) &
                                    + e%f3( 3, i1+1, i2+1, i3-1 ) - e%f3( 3, i1+1, i2 , i3-1 ) &
                                    + e%f3( 3, i1-1, i2+1, i3+1 ) - e%f3( 3, i1-1, i2 , i3+1 ))) &
                    + dtdx3 * ( k1 * ( e%f3( 2, i1 , i2 , i3+1 ) - e%f3( 2, i1 , i2 , i3 )) &
                            + k2 * ( e%f3( 2, i1 , i2+1, i3+1 ) - e%f3( 2, i1 , i2+1, i3 ) &
                                    + e%f3( 2, i1 , i2-1, i3+1 ) - e%f3( 2, i1 , i2-1, i3 ) &
                                    + e%f3( 2, i1+1, i2 , i3+1 ) - e%f3( 2, i1+1, i2 , i3 ) &
                                    + e%f3( 2, i1-1, i2 , i3+1 ) - e%f3( 2, i1-1, i2 , i3 )) &
                            + k3 * ( e%f3( 2, i1+1, i2+1, i3+1 ) - e%f3( 2, i1+1, i2+1, i3 ) &
                                    + e%f3( 2, i1-1, i2-1, i3+1 ) - e%f3( 2, i1-1, i2-1, i3 ) &
                                    + e%f3( 2, i1+1, i2-1, i3+1 ) - e%f3( 2, i1+1, i2-1, i3 ) &
                                    + e%f3( 2, i1-1, i2+1, i3+1 ) - e%f3( 2, i1-1, i2+1, i3 )))

            !B2
            b%f3( 2, i1, i2, i3 ) = b%f3( 2, i1, i2, i3 ) &
                    + dtdx1 * ( k1 * ( e%f3( 3, i1+1, i2 , i3 ) - e%f3( 3, i1 , i2 , i3 )) &
                            + k2 * ( e%f3( 3, i1+1, i2+1, i3 ) - e%f3( 3, i1 , i2+1, i3 ) &
                                    + e%f3( 3, i1+1, i2-1, i3 ) - e%f3( 3, i1 , i2-1, i3 ) &
                                    + e%f3( 3, i1+1, i2 , i3+1 ) - e%f3( 3, i1 , i2 , i3+1 ) &
                                    + e%f3( 3, i1+1, i2 , i3-1 ) - e%f3( 3, i1 , i2 , i3-1 )) &
                            + k3 * ( e%f3( 3, i1+1, i2+1, i3+1 ) - e%f3( 3, i1 , i2+1, i3+1 ) &
                                    + e%f3( 3, i1+1, i2-1, i3-1 ) - e%f3( 3, i1 , i2-1, i3-1 ) &
                                    + e%f3( 3, i1+1, i2+1, i3-1 ) - e%f3( 3, i1 , i2+1, i3-1 ) &
                                    + e%f3( 3, i1+1, i2-1, i3+1 ) - e%f3( 3, i1 , i2-1, i3+1 ))) &
                    - dtdx3 * ( k1 * ( e%f3( 1, i1 , i2 , i3+1 ) - e%f3( 1, i1 , i2 , i3 )) &
                            + k2 * ( e%f3( 1, i1+1, i2 , i3+1 ) - e%f3( 1, i1+1, i2 , i3 ) &
                                    + e%f3( 1, i1-1, i2 , i3+1 ) - e%f3( 1, i1-1, i2 , i3 ) &
                                    + e%f3( 1, i1 , i2+1, i3+1 ) - e%f3( 1, i1 , i2+1, i3 ) &
                                    + e%f3( 1, i1 , i2-1, i3+1 ) - e%f3( 1, i1 , i2-1, i3 )) &
                            + k3 * ( e%f3( 1, i1+1, i2+1, i3+1 ) - e%f3( 1, i1+1, i2+1, i3 ) &
                                    + e%f3( 1, i1-1, i2-1, i3+1 ) - e%f3( 1, i1-1, i2-1, i3 ) &
                                    + e%f3( 1, i1+1, i2-1, i3+1 ) - e%f3( 1, i1+1, i2-1, i3 ) &
                                    + e%f3( 1, i1-1, i2+1, i3+1 ) - e%f3( 1, i1-1, i2+1, i3 )))

            !B3
            b%f3( 3, i1, i2, i3 ) = b%f3( 3, i1, i2, i3 ) &
                    - dtdx1 * ( k1 * ( e%f3( 2, i1+1, i2 , i3 ) - e%f3( 2, i1 , i2 , i3 )) &
                            + k2 * ( e%f3( 2, i1+1, i2+1, i3 ) - e%f3( 2, i1 , i2+1, i3 ) &
                                    + e%f3( 2, i1+1, i2-1, i3 ) - e%f3( 2, i1 , i2-1, i3 ) &
                                    + e%f3( 2, i1+1, i2 , i3+1 ) - e%f3( 2, i1 , i2 , i3+1 ) &
                                    + e%f3( 2, i1+1, i2 , i3-1 ) - e%f3( 2, i1 , i2 , i3-1 )) &
                            + k3 * ( e%f3( 2, i1+1, i2+1, i3+1 ) - e%f3( 2, i1 , i2+1, i3+1 ) &
                                    + e%f3( 2, i1+1, i2-1, i3-1 ) - e%f3( 2, i1 , i2-1, i3-1 ) &
                                    + e%f3( 2, i1+1, i2+1, i3-1 ) - e%f3( 2, i1 , i2+1, i3-1 ) &
                                    + e%f3( 2, i1+1, i2-1, i3+1 ) - e%f3( 2, i1 , i2-1, i3+1 ))) &
                    + dtdx2 * ( k1 * ( e%f3( 1, i1 , i2+1, i3 ) - e%f3( 1, i1 , i2 , i3 )) &
                            + k2 * ( e%f3( 1, i1+1, i2+1, i3 ) - e%f3( 1, i1+1, i2 , i3 ) &
                                    + e%f3( 1, i1-1, i2+1, i3 ) - e%f3( 1, i1-1, i2 , i3 ) &
                                    + e%f3( 1, i1 , i2+1, i3+1 ) - e%f3( 1, i1 , i2 , i3+1 ) &
                                    + e%f3( 1, i1 , i2+1, i3-1 ) - e%f3( 1, i1 , i2 , i3-1 )) &
                            + k3 * ( e%f3( 1, i1+1, i2+1, i3+1 ) - e%f3( 1, i1+1, i2 , i3+1 ) &
                                    + e%f3( 1, i1-1, i2+1, i3-1 ) - e%f3( 1, i1-1, i2 , i3-1 ) &
                                    + e%f3( 1, i1+1, i2+1, i3-1 ) - e%f3( 1, i1+1, i2 , i3-1 ) &
                                    + e%f3( 1, i1-1, i2+1, i3+1 ) - e%f3( 1, i1-1, i2 , i3+1 )))

            enddo
        enddo
    enddo
    !$omp end parallel do

end subroutine dbdt_3d_ck

subroutine dedt_2d_ck( this, e, b, jay, dt )

    implicit none

    class( t_emf_solver_ck ), intent(in) :: this
    type( t_vdf ), intent(inout) :: e
    type( t_vdf ), intent(in) :: b, jay
    real(p_double), intent(in) :: dt

    integer :: i1, i2
    real(p_k_fld) :: dtdx1, dtdx2, dtif

    dtdx1 = real( dt/e%dx_(1), p_k_fld )
    dtdx2 = real( dt/e%dx_(2), p_k_fld )
    dtif = real( dt, p_k_fld )

    !$omp parallel do private(i1)
    do i2 = 0, e%nx_(2)+1
        do i1 = 0, e%nx_(1)+1
            ! E1
            e%f2(1, i1, i2) = e%f2(1, i1, i2) &
                            - dtif * jay%f2(1, i1, i2) &
                            + dtdx2 * ( b%f2(3, i1, i2) - b%f2(3, i1, i2-1) )

            ! E2
            e%f2(2, i1, i2) = e%f2(2, i1, i2) &
                            - dtif * jay%f2(2, i1, i2) &
                            - dtdx1 * ( b%f2(3, i1, i2) - b%f2(3, i1-1, i2) )

            ! E3
            e%f2(3, i1, i2) = e%f2(3, i1, i2) &
                            - dtif * jay%f2(3, i1, i2) &
                            + dtdx1 * ( b%f2(2, i1, i2) - b%f2(2, i1-1, i2)) &
                            - dtdx2 * ( b%f2(1, i1, i2) - b%f2(1, i1, i2-1))
        enddo
    enddo
    !$omp end parallel do

end subroutine dedt_2d_ck

subroutine dedt_3d_ck( this, e, b, jay, dt )

    implicit none

    class( t_emf_solver_ck ), intent(in) :: this
    type( t_vdf ), intent(inout) :: e
    type( t_vdf ), intent(in) :: b, jay

    real(p_double), intent(in) :: dt

    real(p_k_fld) :: dtdx1, dtdx2, dtdx3, dtif
    integer :: i1, i2, i3

    dtdx1 = real( dt/e%dx_(1), p_k_fld )
    dtdx2 = real( dt/e%dx_(2), p_k_fld )
    dtdx3 = real( dt/e%dx_(3), p_k_fld )
    dtif = real( dt, p_k_fld )

    !$omp parallel do private(i2,i1)
    do i3 = 0, e%nx_(3)+1
        do i2 = 0, e%nx_(2)+1
            do i1 = 0, e%nx_(1)+1

            e%f3( 1, i1, i2, i3 ) = e%f3( 1, i1, i2, i3 ) &
                - dtif * jay%f3( 1, i1, i2, i3 ) &
                + dtdx2 * ( b%f3( 3, i1, i2, i3 ) - b%f3( 3, i1, i2-1, i3 ) ) &
                - dtdx3 * ( b%f3( 2, i1, i2, i3 ) - b%f3( 2, i1, i2, i3-1 ) )

            e%f3( 2, i1, i2, i3 ) = e%f3( 2, i1, i2, i3 ) &
                - dtif * jay%f3( 2, i1, i2, i3 ) &
                + dtdx3 * ( b%f3( 1, i1, i2, i3 ) - b%f3( 1, i1, i2, i3-1 ) ) &
                - dtdx1 * ( b%f3( 3, i1, i2, i3 ) - b%f3( 3, i1-1, i2, i3 ) )

            e%f3( 3, i1, i2, i3 ) = e%f3( 3, i1, i2, i3 ) &
                - dtif * jay%f3( 3, i1, i2, i3 ) &
                + dtdx1 * ( b%f3( 2, i1, i2, i3 ) - b%f3( 2, i1-1, i2, i3 ) ) &
                - dtdx2 * ( b%f3( 1, i1, i2, i3 ) - b%f3( 1, i1, i2-1, i3 ) )

            enddo
        enddo
    enddo
    !$omp end parallel do

end subroutine dedt_3d_ck

subroutine test_stability_ck( this, dt, dx )

    implicit none

    class( t_emf_solver_ck ), intent(in) :: this
    real(p_double), intent(in) :: dt
    real(p_double), dimension(:), intent(in) :: dx

    real(p_double) :: cour
    integer :: i

    cour = dx(1)

    do i = 2, p_x_dim
        if ( dx(1) /= dx(i) ) then
            if ( mpi_node() == 0 ) then
                print *, '(*error*) C-K EMF solver stability condition violated, aborting'
                print *, 'Please ensure that dx1 = dx2 = dx3.'
            endif
            call abort_program(-32)
        endif
    enddo


    if (dt > cour) then
        if ( mpi_node() == 0 ) then
           print *, '(*error*) C-K EMF solver stability condition violated, aborting'
           print *, 'dx   = ', dx(1:p_x_dim)
           print *, 'dt   = ', dt, ' max(dt) = ', cour
        endif
        call abort_program(p_err_invalid)
    endif

end subroutine test_stability_ck

subroutine min_gc_ck( this, min_gc )

    implicit none

    class( t_emf_solver_ck ), intent(in) :: this
    integer, dimension(:,:), intent(inout) :: min_gc

    min_gc(p_lower,:) = 3
    min_gc(p_upper,:) = 2

end subroutine


subroutine init_ck( this, gix_pos, coords )

    implicit none

    class( t_emf_solver_ck ), intent(inout) :: this
    integer, intent(in), dimension(:) :: gix_pos
    integer, intent(in) :: coords

    select case (p_x_dim)
    case(1)
        if ( mpi_node() == 0 ) then
            write(0,*) "(*error*) CK EMF solver not implemented in 1D"
            write(0,*) "(*error*) aborting..."
        endif
        call abort_program( p_err_notimplemented )
    case(2)
        if ( coords == p_cylindrical_b ) then
            if ( mpi_node() == 0 ) then
                write(0,*) "(*error*) CK EMF solver not implemented in 2D cylindrical geometry"
                write(0,*) "(*error*) aborting..."
            endif
            call abort_program( p_err_notimplemented )
        else
            this % dedt => dedt_2d_ck
            this % dbdt => dbdt_2d_ck
            this % k1 = 0.583333333_p_k_fld
            this % k2 = 0.20833333_p_k_fld
         endif
    case(3)
        this % dedt => dedt_3d_ck
        this % dbdt => dbdt_3d_ck
        this % k1 = 0.583333333_p_k_fld
        this % k2 = 0.083333333_p_k_fld
        this % k3 = 0.020833333_p_k_fld

    end select

end subroutine init_ck

subroutine read_input_ck( this, input_file, dx, dt )

    implicit none

    class( t_emf_solver_ck ), intent(inout) :: this
    class( t_input_file ), intent(inout) :: input_file
    real(p_double), dimension(:), intent(in) :: dx
    real(p_double), intent(in) :: dt

    real :: dummy_
    integer :: ierr

    namelist /nl_emf_solver/ dummy_

    dummy_ = huge(1.0)

    call get_namelist( input_file, "nl_emf_solver", ierr )

    if ( ierr == 0 ) then
        read (input_file%nml_text, nml = nl_emf_solver, iostat = ierr)
        if (ierr /= 0) then
            if ( mpi_node() == 0 ) then
                write(0,*) ""
                write(0,*) "   Error reading C-K EMF solver parameters"
                write(0,*) "   aborting..."
            endif
            stop
        endif
    endif

    if ( dummy_ /= huge(1.0) ) then
        if ( mpi_node() == 0 ) then
            write(0,*) ""
            write(0,*) "   Error reading C-K EMF solver parameters"
            write(0,*) "   This solver does not accept any input parameters."
            write(0,*) "   aborting..."
        endif
        stop
    endif

end subroutine read_input_ck

subroutine advance_ck( this, e, b, jay, dt, bnd_con )

    implicit none

    class( t_emf_solver_ck ), intent(inout) :: this
    class( t_vdf ), intent( inout ) :: e, b, jay
    class( t_emf_bound ), intent(inout) :: bnd_con
    real(p_double), intent(in) :: dt

    real(p_double) :: dt_b, dt_e

    dt_b = dt / 2.0_p_double
    dt_e = dt

    ! Advance B half time step
    call this % dbdt( b, e, dt_b )
    call bnd_con%update_boundary_b( e, b, step = 1 )

    ! Advance E one full time step
    call this % dedt( e, b, jay, dt_e )
    call bnd_con%update_boundary_e( e, b )

    ! Advance B another half time step
    call this % dbdt( b, e, dt_b )
    call bnd_con%update_boundary_b( e, b, step = 2 )

end subroutine advance_ck

end module m_emf_solver_ck
