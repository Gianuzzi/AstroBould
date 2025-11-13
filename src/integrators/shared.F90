!> Module with shared variables and routines for integrators
module shared
    implicit none

    ! Constants
    real(kind=8), parameter :: ZERO = 0.0
    real(kind=8), parameter :: ONE = 1.d0
    real(kind=8), parameter :: TWO = 2.d0
    real(kind=8), parameter :: THREE = 3.d0
    real(kind=8), parameter :: FOUR = 4.d0
    real(kind=8), parameter :: FIVE = 5.d0
    real(kind=8), parameter :: C1_2 = 0.5d0
    real(kind=8), parameter :: C1_3 = ONE/3.d0
    real(kind=8), parameter :: C1_4 = 0.25d0
    real(kind=8), parameter :: C1_5 = 0.2d0
    real(kind=8), parameter :: C1_6 = ONE/6.d0
    real(kind=8), parameter :: C1_7 = ONE/7.d0
    real(kind=8), parameter :: C1_8 = 0.125d0
    real(kind=8), parameter :: C1_9 = ONE/9.d0
    real(kind=8), parameter :: C1_12 = ONE/12.d0
    real(kind=8), parameter :: C1_15 = ONE/15.d0
    real(kind=8), parameter :: C1_18 = ONE/18.d0
    real(kind=8), parameter :: C1_48 = ONE/48.d0
    real(kind=8), parameter :: C1_840 = ONE/840.d0
    real(kind=8), parameter :: C2_3 = TWO/3.d0
    real(kind=8), parameter :: C2_5 = TWO/5.d0
    real(kind=8), parameter :: C2_27 = TWO/27.d0
    real(kind=8), parameter :: C2_45 = TWO/45.d0
    real(kind=8), parameter :: C3_2 = 1.5d0
    real(kind=8), parameter :: C3_4 = 0.75d0
    real(kind=8), parameter :: C3_8 = 0.375d0
    real(kind=8), parameter :: C3_16 = THREE/16.d0
    real(kind=8), parameter :: C4_15 = FOUR/15.d0
    real(kind=8), parameter :: C5_6 = FIVE/6.d0
    real(kind=8), parameter :: C5_12 = FIVE/12.d0
    real(kind=8), parameter :: C8_9 = 8.d0/9.d0
    real(kind=8), parameter :: SQ3_6 = sqrt(3.d0)/6.d0

    real(kind=8), parameter :: SAFE_LOW = 1d-30

    ! BASIC CONFIGURATION INPUT
    integer(kind=4) :: NDIM = 1  ! Number of dimensions
    real(kind=8) :: E_TOL = 1.d-14 ! Error tolerance
    real(kind=8) :: BETA = 0.15d0  ! Learning rate
    real(kind=8) :: DT_MIN = 1.d-6  ! Minimum dt
    integer(kind=4) :: EXTRA = 0  ! Amount of extra variables (without der) that are not POS
    !! derived
    integer(kind=4) :: NDIM2
    integer(kind=4) :: EXTRA2

    ! Global
    real(kind=8), allocatable :: der(:) ! Used in callers

    abstract interface

        !!!! TEMPLATES

        ! Keep going template
        function function_check_keep_tem (y) result(keep_going)
            implicit none
            real(kind=8), dimension(:), intent(in) :: y
            logical :: keep_going
        end function function_check_keep_tem

        ! f__ (t, y__) = (f_i (t, y__), ..., f_n (t, y__)) = der__
        function dydt_tem (t, y) result (derivate)
            implicit none
            real(kind=8), intent(in)               :: t
            real(kind=8), dimension(:), intent(in) :: y
            real(kind=8), dimension(size (y))      :: derivate
        end function dydt_tem

    end interface
    
    contains

        pure function get_index(i) result(idx)
            implicit none
            integer(kind=4), intent(in) :: i
            integer(kind=4) :: idx
            idx = EXTRA2 + NDIM2 * (i - 1) + 1
        end function get_index

    
end module shared