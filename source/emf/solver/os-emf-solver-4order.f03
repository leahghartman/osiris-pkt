!-----------------------------------------------------------------------------------------
! 4th order accurate field solver
! A. Taflove & S.C.Hagness, "Computational Electrodynamics: The Finite-Difference
!   Time-Domain Method", Chapter 4, Section 4.9.2, pp. 143
!-----------------------------------------------------------------------------------------

#include "os-config.h"
#include "os-preprocess.fpp"

module m_emf_solver_4order

use m_parameters
use m_restart
use m_emf_define
use m_vdf_define
use m_input_file
private  

type, extends( t_emf_solver ), public :: t_emf_solver_4order

    procedure(dedt_4order), pointer :: dedt => null()
    procedure(dbdt_4order), pointer :: dbdt => null()

contains

    procedure :: init => init_4order
    procedure :: advance => advance_4order
    procedure :: min_gc => min_gc_4order
    procedure :: name => name_4order
    procedure :: read_input => read_input_4order
    procedure :: test_stability => test_stability_4order

end type t_emf_solver_4order

interface
    subroutine dedt_4order( this, e, b, jay, dt )
        import t_emf_solver_4order, t_vdf, p_double
        class( t_emf_solver_4order ), intent(in) :: this
        class( t_vdf ), intent(inout) :: e
        class( t_vdf ), intent(in) :: b, jay
        real(p_double), intent(in) :: dt
    end subroutine
end interface

interface
    subroutine dbdt_4order( this, b, e, dt )
        import t_emf_solver_4order, t_vdf, p_double
        class( t_emf_solver_4order ), intent(in) :: this
        class( t_vdf ), intent(inout) :: b
        class( t_vdf ), intent(in) :: e
        real(p_double), intent(in) :: dt
    end subroutine
end interface

contains

function name_4order( this )
  
    implicit none
  
    class( t_emf_solver_4order ), intent(in) :: this
    character( len = p_maxlen_emf_solver_name ) :: name_4order
  
    name_4order = "4th order accurate"
  
end function name_4order

subroutine dbdt_1d_4order( this, b, e, dt )

    implicit none

    class( t_emf_solver_4order ), intent(in) :: this
    class( t_vdf ), intent( inout ) :: b
    class( t_vdf ), intent(in) :: e
    real(p_double), intent(in) :: dt

    real(p_k_fld) :: dtdx1
    integer :: i1
    dtdx1 = real( dt/(24 * b%dx_(1)), p_k_fld)

    !$omp parallel do
    do i1 = -2, b%nx_(1)+3
    !B1
    ! b%f1( 1, i1 ) = b%f1( 1, i1 )

    !B2
    b%f1( 2, i1 ) = b%f1( 2, i1 ) &
                        + dtdx1 * ( 27 * ( e%f1( 3, i1+1 ) - e%f1( 3, i1   ) ) &
                                       - ( e%f1( 3, i1+2 ) - e%f1( 3, i1-1 ) ) )

    !B3
    b%f1( 3, i1 ) = b%f1( 3, i1 ) &
                        - dtdx1 * ( 27 * ( e%f1( 2, i1+1 ) - e%f1( 2, i1   ) ) &
                                       - ( e%f1( 2, i1+2 ) - e%f1( 2, i1-1 ) ) )
    enddo
    !$omp end parallel do

end subroutine dbdt_1d_4order
    
subroutine dbdt_2d_4order( this, b, e, dt )

    implicit none

    class( t_emf_solver_4order ), intent(in) :: this
    class( t_vdf ), intent( inout ) :: b
    class( t_vdf ), intent(in) :: e
    real(p_double), intent(in) :: dt

    real(p_k_fld) :: dtdx1, dtdx2
    integer :: i1, i2

    dtdx1 = real( dt/(24 * b%dx_(1)), p_k_fld)
    dtdx2 = real( dt/(24 * b%dx_(2)), p_k_fld)

    !$omp parallel do private(i1)
    do i2 = -2, b%nx_(2)+3
        do i1 = -2, b%nx_(1)+3

            !B1
            b%f2( 1, i1, i2 ) = b%f2( 1, i1, i2 ) &
                                - dtdx2 * ( 27 * ( e%f2( 3, i1, i2+1 ) - e%f2( 3, i1, i2 ) ) &
                                                - ( e%f2( 3, i1, i2+2 ) - e%f2( 3, i1, i2-1 ) ) )

            !B2
            b%f2( 2, i1, i2 ) = b%f2( 2, i1, i2 ) &
                                + dtdx1 * ( 27 * ( e%f2( 3, i1+1, i2 ) - e%f2( 3, i1  , i2 ) )  &
                                                - ( e%f2( 3, i1+2, i2 ) - e%f2( 3, i1-1, i2 ) ) )

            !B3
            b%f2( 3, i1, i2 ) = b%f2( 3, i1, i2 ) &
                                - dtdx1 * (  27 * ( e%f2( 2, i1+1, i2 ) - e%f2( 2, i1, i2 ) ) &
                                                - ( e%f2( 2, i1+2, i2 ) - e%f2( 2, i1-1, i2 ) ) ) &
                                + dtdx2 * (  27 * ( e%f2( 1, i1, i2+1 ) - e%f2( 1, i1, i2 ) ) &
                                                - ( e%f2( 1, i1, i2+2 ) - e%f2( 1, i1, i2-1 ) ) )

        enddo
    enddo
    !$omp end parallel do

end subroutine dbdt_2d_4order

subroutine dbdt_cyl_2d_4order( this, b, e, dt )

    implicit none

    class( t_emf_solver_4order ), intent(in) :: this
    class( t_vdf ), intent( inout ) :: b
    class( t_vdf ), intent( in )    :: e
    real(p_double),  intent(in)    :: dt

    real(p_k_fld) :: dtdz, dtdr
    real(p_k_fld) :: tmp, rm, rp, rmm, rpp
    integer :: i1, i2, i2_0, gshift_i2

    dtdz = real( dt/(24 * b%dx_(1)), p_k_fld )
    dtdr = real( dt/(24 * b%dx_(2)), p_k_fld )

    gshift_i2 = this%gir_pos - 2

    if ( gshift_i2 < 0 ) then
        ! This node contains the cylindrical axis so start the solver in cell 2
        i2_0 = 2
    else
        i2_0 = 0
    endif

    do i2 = i2_0, b%nx_(2)+3

        rpp  = i2 + gshift_i2 + 1.5     ! r position of the corner of cell i2 + 2
        rp   = i2 + gshift_i2 + 0.5     ! r position of the corner of cell i2 + 1
        rm   = i2 + gshift_i2 - 0.5     ! r position of the corner of cell i2
        rmm  = i2 + gshift_i2 - 1.5     ! r position of the corner of cell i2 - 1

        tmp  = dtdr/(i2 + gshift_i2)    ! (dt/dr) / position of the middle of cell i2

        do i1 = -2, b%nx_(1)+3
            !B1
            b%f2( 1, i1, i2 ) = b%f2(1, i1, i2) - tmp * &
                                    (     - rpp * e%f2(3, i1, i2+2 ) + &
                                        27 * rp * e%f2(3, i1, i2+1 ) - &
                                        27 * rm * e%f2(3, i1, i2   ) + &
                                            rmm * e%f2(3, i1, i2-1 ) )

            !B2
            b%f2( 2, i1, i2 ) = b%f2( 2, i1, i2 ) + dtdz * &
                                    (     - e%f2( 3, i1+2, i2 ) + &
                                        27 * e%f2( 3, i1+1, i2 ) - &
                                        27 * e%f2( 3, i1  , i2 ) + &
                                            e%f2( 3, i1-1, i2 ) )

            !B3
            b%f2( 3, i1, i2 ) = b%f2( 3, i1, i2 ) &
                                    - dtdz * (    - e%f2( 2, i1+2, i2 ) + &
                                                27 * e%f2( 2, i1+1, i2 ) - &
                                                27 * e%f2( 2, i1  , i2 ) + &
                                                    e%f2( 2, i1-1, i2 ) )  &
                                    + dtdr * (    - e%f2( 1, i1, i2+2 ) + &
                                                27 * e%f2( 1, i1, i2+1 ) - &
                                                27 * e%f2( 1, i1, i2   ) + &
                                                    e%f2( 1, i1, i2-1 ) )

        enddo
    enddo

    if ( gshift_i2 < 0 ) then

        ! guard cells (i2 = -2, -1, 0)
        do i2 = -2, 0
            do i1 = -2, b%nx_(1)+3
            b%f2( 1, i1, i2 ) =   b%f2( 1, i1, 2 - i2 )
            b%f2( 2, i1, i2 ) = - b%f2( 2, i1, 3 - i2 )
            b%f2( 3, i1, i2 ) = - b%f2( 3, i1, 2 - i2 )
            enddo
        enddo

        ! axial cell (i2 = 1)
        do i1 = -2, b%nx_(1)+3
            ! rmid = 0 in this cell so the rot B component is calculated differently
            ! this is not consistent with the 4th order derivative.
            b%f2( 1, i1, 1 ) =   b%f2(1, i1, 1) - 4 * dtdr * e%f2(3, i1, 2)

            ! B2 is not defined on axis so just reflect the corresponding value inside the box
            b%f2( 2, i1, 1 ) = - b%f2( 2, i1, 2 )

            ! B3 is zero on axis
            b%f2( 3, i1, 1 ) = 0.0_p_k_fld
        enddo

    endif

end subroutine dbdt_cyl_2d_4order
    
subroutine dbdt_3d_4order( this, b, e, dt )

    implicit none

    class( t_emf_solver_4order ), intent(in) :: this
    class( t_vdf ), intent( inout ) :: b
    class( t_vdf ), intent(in) :: e
    real(p_double), intent(in) :: dt

    real(p_k_fld) :: dtdx1, dtdx2, dtdx3
    integer :: i1, i2, i3

    dtdx1 = real( dt/(24 * b%dx_(1)), p_k_fld)
    dtdx2 = real( dt/(24 * b%dx_(2)), p_k_fld)
    dtdx3 = real( dt/(24 * b%dx_(3)), p_k_fld)

    !$omp parallel do private(i1,i2)
    do i3 = -2, b%nx_(3)+3
        do i2 = -2, b%nx_(2)+3
            do i1 = -2, b%nx_(1)+3

            !B1
            b%f3( 1, i1, i2, i3 ) = b%f3( 1, i1, i2, i3 ) &
                                - dtdx2 * ( 27 * ( e%f3( 3, i1, i2+1, i3 ) - e%f3( 3, i1, i2  , i3 ) ) &
                                            - ( e%f3( 3, i1, i2+2, i3 ) - e%f3( 3, i1, i2-1, i3 ) ) ) &
                                + dtdx3 * ( 27 * ( e%f3( 2, i1, i2, i3+1 ) - e%f3( 2, i1, i2  , i3 ) ) &
                                            - ( e%f3( 2, i1, i2, i3+2 ) - e%f3( 2, i1, i2, i3-1 ) ) )

            !B2
            b%f3( 2, i1, i2, i3 ) = b%f3( 2, i1, i2, i3 ) &
                                + dtdx1 * ( 27 * ( e%f3( 3, i1+1, i2, i3 ) - e%f3( 3, i1  , i2, i3 ) )  &
                                            - ( e%f3( 3, i1+2, i2, i3 ) - e%f3( 3, i1-1, i2, i3 ) ) ) &
                                - dtdx3 * ( 27 * ( e%f3( 1, i1, i2, i3+1 ) - e%f3( 1, i1, i2, i3   ) )  &
                                            - ( e%f3( 1, i1, i2, i3+2 ) - e%f3( 1, i1, i2, i3-1 ) ) )

            !B3
            b%f3( 3, i1, i2, i3 ) = b%f3( 3, i1, i2, i3 ) &
                                - dtdx1 * ( 27 * ( e%f3( 2, i1+1, i2, i3 ) - e%f3( 2, i1  , i2, i3 ) ) &
                                            - ( e%f3( 2, i1+2, i2, i3 ) - e%f3( 2, i1-1, i2, i3 ) ) ) &
                                + dtdx2 * ( 27 * ( e%f3( 1, i1, i2+1, i3 ) - e%f3( 1, i1, i2  , i3 ) ) &
                                            - ( e%f3( 1, i1, i2+2, i3 ) - e%f3( 1, i1, i2-1, i3 ) ) )

            enddo
        enddo
    enddo
    !$omp end parallel do

end subroutine dbdt_3d_4order
    
    
subroutine dedt_1d_4order( this, e, b, jay, dt )

    implicit none

    class( t_emf_solver_4order ), intent(in) :: this
    class( t_vdf ), intent(inout) :: e
    class( t_vdf ), intent(in) :: b, jay
    real(p_double), intent(in) :: dt

    real(p_k_fld) :: dtdx1, dtif
    integer :: i1

    dtdx1 = real( dt/(24 * e%dx_(1)), p_k_fld )
    dtif = real( dt, p_k_fld )

    ! version 1, advance 1 cell at a time

    !$omp parallel do
    do i1 = 0, e%nx_(1)+2
        ! E1
        e%f1(1, i1) = e%f1(1, i1) - dtif * jay%f1(1, i1)

        ! E2
        e%f1(2, i1) = e%f1(2, i1) &
                        - dtif * jay%f1(2, i1) &
                        - dtdx1 * ( 27 * ( b%f1(3, i1  ) - b%f1(3, i1-1) ) &
                                        - ( b%f1(3, i1+1) - b%f1(3, i1-2) ) )

        ! E3
        e%f1(3, i1) = e%f1(3, i1) &
                        - dtif * jay%f1(3, i1) &
                        + dtdx1 * ( 27 * ( b%f1(2, i1  ) - b%f1(2, i1-1) ) &
                                        - ( b%f1(2, i1+1) - b%f1(2, i1-2) ) )
    enddo
    !$omp end parallel do

end subroutine dedt_1d_4order
    
subroutine dedt_2d_4order( this, e, b, jay, dt )

    implicit none

    class( t_emf_solver_4order ), intent(in) :: this
    class( t_vdf ), intent(inout) :: e
    class( t_vdf ), intent(in) :: b, jay
    real(p_double), intent(in) :: dt

    integer :: i1, i2
    real(p_k_fld) :: dtdx1, dtdx2, dtif

    dtdx1 = real( dt/(24 * e%dx_(1)), p_k_fld )
    dtdx2 = real( dt/(24 * e%dx_(2)), p_k_fld )
    dtif = real( dt, p_k_fld )

    !$omp parallel do private(i1)
    do i2 = 0, e%nx_(2)+2
        do i1 = 0, e%nx_(1)+2
            ! E1
            e%f2(1, i1, i2) = e%f2(1, i1, i2) &
                            - dtif * jay%f2(1, i1, i2) &
                            + dtdx2 * ( 27 * ( b%f2(3, i1, i2  ) - b%f2(3, i1, i2-1) ) &
                                            - ( b%f2(3, i1, i2+1) - b%f2(3, i1, i2-2) ) )

            ! E2
            e%f2(2, i1, i2) = e%f2(2, i1, i2) &
                            - dtif * jay%f2(2, i1, i2) &
                            - dtdx1 * ( 27 * ( b%f2(3, i1  , i2) - b%f2(3, i1-1, i2) ) &
                                            - ( b%f2(3, i1+1, i2) - b%f2(3, i1-2, i2) ) )

            ! E3
            e%f2(3, i1, i2) = e%f2(3, i1, i2) &
                            - dtif * jay%f2(3, i1, i2) &
                            + dtdx1 * ( 27 * ( b%f2(2, i1  , i2) - b%f2(2, i1-1, i2) ) &
                                            - ( b%f2(2, i1+1, i2) - b%f2(2, i1-2, i2) ) ) &
                            - dtdx2 * ( 27 * ( b%f2(1, i1,   i2) - b%f2(1, i1, i2-1) ) &
                                            - ( b%f2(1, i1, i2+1) - b%f2(1, i1, i2-2) ) )
        enddo
    enddo
    !$omp end parallel do

end subroutine dedt_2d_4order
    
subroutine dedt_3d_4order( this, e, b, jay, dt )

    implicit none

    class( t_emf_solver_4order ), intent(in) :: this
    class( t_vdf ), intent(inout) :: e
    class( t_vdf ), intent(in) :: b, jay
    real(p_double), intent(in) :: dt

    integer :: i1, i2, i3
    real(p_k_fld) :: dtdx1, dtdx2, dtdx3, dtif

    dtdx1 = real( dt/(24 * e%dx_(1)), p_k_fld )
    dtdx2 = real( dt/(24 * e%dx_(2)), p_k_fld )
    dtdx3 = real( dt/(24 * e%dx_(3)), p_k_fld )
    dtif = real( dt, p_k_fld )

    !$omp parallel do private(i1,i2)
    do i3 = 0, e%nx_(3)+2
        do i2 = 0, e%nx_(2)+2
            do i1 = 0, e%nx_(1)+2
            ! E1
            e%f3(1, i1, i2, i3) = e%f3(1, i1, i2, i3) &
                                - dtif * jay%f3(1, i1, i2, i3) &
                                + dtdx2 * ( 27 * ( b%f3(3, i1, i2  , i3) - b%f3(3, i1, i2-1, i3) ) &
                                                - ( b%f3(3, i1, i2+1, i3) - b%f3(3, i1, i2-2, i3) ) ) &
                                - dtdx3 * ( 27 * ( b%f3(2, i1, i2, i3  ) - b%f3(2, i1, i2, i3-1) ) &
                                                - ( b%f3(2, i1, i2, i3+1) - b%f3(2, i1, i2, i3-2) ) )


            ! E2
            e%f3(2, i1, i2, i3) = e%f3(2, i1, i2, i3) &
                                - dtif * jay%f3(2, i1, i2, i3) &
                                - dtdx1 * ( 27 * ( b%f3(3, i1  , i2, i3) - b%f3(3, i1-1, i2, i3) ) &
                                                - ( b%f3(3, i1+1, i2, i3) - b%f3(3, i1-2, i2, i3) ) ) &
                                + dtdx3 * ( 27 * ( b%f3(1, i1, i2, i3  ) - b%f3(1, i1, i2, i3-1) ) &
                                                - ( b%f3(1, i1, i2, i3+1) - b%f3(1, i1, i2, i3-2) ) )


            ! E3
            e%f3(3, i1, i2, i3) = e%f3(3, i1, i2, i3) &
                                - dtif * jay%f3(3, i1, i2, i3) &
                                + dtdx1 * ( 27 * ( b%f3(2, i1  , i2, i3) - b%f3(2, i1-1, i2, i3) ) &
                                                - ( b%f3(2, i1+1, i2, i3) - b%f3(2, i1-2, i2, i3) ) ) &
                                - dtdx2 * ( 27 * ( b%f3(1, i1,   i2, i3) - b%f3(1, i1, i2-1, i3) ) &
                                                - ( b%f3(1, i1, i2+1, i3) - b%f3(1, i1, i2-2, i3) ) )
            enddo
        enddo
    enddo
    !$omp end parallel do

end subroutine dedt_3d_4order
    
subroutine dedt_cyl_2d_4order( this, e, b, jay, dt )

    implicit none

    class( t_emf_solver_4order ), intent(in) :: this
    class( t_vdf ), intent(inout) :: e
    class( t_vdf ), intent(in) :: b, jay
    real(p_double),  intent(in)    :: dt

    real(p_k_fld) :: dtdz, dtdr, dtif
    real(p_k_fld) :: tmp, rcp, rcpp, rcm, rcmm
    integer :: i1, i2, gshift_i2

    dtdz = real( dt/(24 * b%dx_(1)), p_k_fld )
    dtdr = real( dt/(24 * b%dx_(2)), p_k_fld )
    dtif = real( dt, p_k_fld )

    gshift_i2 = this % gir_pos - 2

    do i2 = 0, e%nx_(2)+2

        tmp  = dtdr / (( i2 + gshift_i2 ) - 0.5)
        rcpp = i2 + gshift_i2 + 1
        rcp  = i2 + gshift_i2
        rcm  = i2 + gshift_i2 - 1
        rcmm = i2 + gshift_i2 - 2

        do i1 = 0, e%nx_(1)+2
            ! E1
            e%f2(1, i1, i2) = e%f2(1, i1, i2) &
                            - dtif * jay%f2(1, i1, i2) &
                            + tmp * (    - rcpp * b%f2(3, i1, i2+1) + &
                                        27 * rcp * b%f2(3, i1, i2  ) - &
                                        27 * rcm * b%f2(3, i1, i2-1) + &
                                            rcmm * b%f2(3, i1, i2-2) )

            ! E2
            e%f2(2, i1, i2) = e%f2(2, i1, i2) &
                            - dtif * jay%f2(2, i1, i2) &
                            - dtdz * (     - b%f2(3, i1+1, i2) + &
                                        27 * b%f2(3, i1,   i2) - &
                                        27 * b%f2(3, i1-1, i2) + &
                                                b%f2(3, i1-2, i2) )

            ! E3
            e%f2(3, i1, i2) = e%f2(3, i1, i2) &
                            - dtif * jay%f2(3, i1, i2) &
                            + dtdz * (     - b%f2(2, i1+1, i2) + &
                                        27 * b%f2(2, i1  , i2) - &
                                        27 * b%f2(2, i1-1, i2) + &
                                                b%f2(2, i1-2, i2) ) &
                            - dtdr * (     - b%f2(1, i1, i2+1) + &
                                        27 * b%f2(1, i1, i2  ) - &
                                        27 * b%f2(1, i1, i2-1) + &
                                                b%f2(1, i1, i2-2) )
        enddo
    enddo

    ! Correct values on axial boundary. (there is a little redundancy between this and update_boundary
    ! it needs to be checked)
    if ( gshift_i2 < 0 ) then

        ! guard cell 1 (i2 = 0)
        do i1 = 0, e%nx_(1)+2
            e%f2( 1, i1, 0 ) =   e%f2( 1, i1, 1 )
            e%f2( 2, i1, 0 ) = - e%f2( 2, i1, 2 )
            e%f2( 3, i1, 0 ) = - e%f2( 3, i1, 1 )
        enddo

        ! axial cell (i2 = 1)
        do i1 = 0, e%nx_(1)+2
            e%f2(1, i1, 1) =   e%f2(1, i1, 2)
            e%f2(2, i1, 1) =   0
            e%f2(3, i1, 1) = - e%f2(3, i1, 2)
        enddo
    endif


end subroutine dedt_cyl_2d_4order

subroutine test_stability_4order( this, dt, dx )
    
    implicit none

    class( t_emf_solver_4order ), intent(in) :: this
    real(p_double), intent(in)  ::  dt
    real(p_double), dimension(:), intent(in)  ::  dx

    real(p_double) :: cour
    integer :: i

    cour = 0.0
    do i = 1, p_x_dim
        cour = cour + 1.0/(dx(i))**2
    enddo
    cour = sqrt(1.0/cour)

    cour = cour * 6.0_p_double / 7.0_p_double

    if (dt > cour) then
        if ( mpi_node() == 0 ) then
            print *, '(*error*) 4th order EMF solver stability condition violated, aborting'
            print *, 'dx   = ', dx(1:p_x_dim)
            print *, 'dt   = ', dt, ' max(dt) = ', cour
         endif
         call abort_program(p_err_invalid)
    endif

end subroutine test_stability_4order

subroutine read_input_4order( this, input_file, dx, dt )

    implicit none

    class( t_emf_solver_4order ), intent(inout) :: this
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
                write(0,*) "   Error reading 4th order EMF solver parameters"
                write(0,*) "   aborting..."
            endif
            stop
        endif
    endif

    if ( dummy_ /= huge(1.0) ) then
        if ( mpi_node() == 0 ) then
            write(0,*) ""
            write(0,*) "   Error reading 4th order EMF solver parameters"
            write(0,*) "   This solver does not accept any input parameters."
            write(0,*) "   aborting..."
        endif
        stop
    endif

end subroutine read_input_4order

subroutine min_gc_4order( this, min_gc )

    implicit none

    class( t_emf_solver_4order ), intent(in) :: this
    integer, dimension(:,:), intent(inout) :: min_gc

    min_gc(p_lower,:) = 4
    min_gc(p_upper,:) = 5

end subroutine

subroutine init_4order( this, gix_pos, coords )

    implicit none

    class( t_emf_solver_4order ), intent(inout) :: this
    integer, intent(in), dimension(:) :: gix_pos
    integer, intent(in) :: coords

    select case (p_x_dim)
    case(1)
        this % dedt => dedt_1d_4order
        this % dbdt => dbdt_1d_4order
    case(2)
        if ( coords == p_cylindrical_b ) then
            ! Cylindrical (r-z) coordinates solvers
            this % dedt => dedt_cyl_2d_4order
            this % dbdt => dbdt_cyl_2d_4order
            ! Radial position of lower cell
            this % gir_pos = gix_pos(2)
        else
            this % dedt => dedt_2d_4order
            this % dbdt => dbdt_2d_4order
        endif
    case(3)
        this % dedt => dedt_3d_4order
        this % dbdt => dbdt_3d_4order
    end select

end subroutine init_4order

subroutine advance_4order( this, e, b, jay, dt, bnd_con )
  
    implicit none
  
    class( t_emf_solver_4order ), intent(inout) :: this
    class( t_vdf ), intent( inout )  ::  e, b, jay
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
    
end subroutine advance_4order

end module m_emf_solver_4order