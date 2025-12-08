!> Module with shared variables and routines for integrators
module shared
    use iso_fortran_env, only: real64, real32  ! or int32, etc.
    implicit none

    ! Default precision
#ifdef WP
    integer, parameter :: wp = WP
#else
    integer, parameter :: wp = real64
#endif

    ! Constants
    real(wp), parameter :: ZERO = 0.0e0_wp
    real(wp), parameter :: ONE = 1.e0_wp
    real(wp), parameter :: TWO = 2.e0_wp
    real(wp), parameter :: THREE = 3.e0_wp
    real(wp), parameter :: FOUR = 4.e0_wp
    real(wp), parameter :: FIVE = 5.e0_wp
    real(wp), parameter :: C1_2 = ONE/TWO
    real(wp), parameter :: C1_3 = ONE/3.e0_wp
    real(wp), parameter :: C1_4 = ONE/FOUR
    real(wp), parameter :: C1_5 = ONE/FIVE
    real(wp), parameter :: C1_6 = ONE/6.e0_wp
    real(wp), parameter :: C1_7 = ONE/7.e0_wp
    real(wp), parameter :: C1_8 = ONE/8.e0_wp
    real(wp), parameter :: C1_9 = ONE/9.e0_wp
    real(wp), parameter :: C1_12 = ONE/12.e0_wp
    real(wp), parameter :: C1_15 = ONE/15.e0_wp
    real(wp), parameter :: C1_18 = ONE/18.e0_wp
    real(wp), parameter :: C1_25 = ONE/25.e0_wp
    real(wp), parameter :: C1_48 = ONE/48.e0_wp
    real(wp), parameter :: C1_192 = ONE/192.e0_wp
    real(wp), parameter :: C1_840 = ONE/840.e0_wp
    real(wp), parameter :: C2_3 = TWO/3.e0_wp
    real(wp), parameter :: C2_5 = TWO/5.e0_wp
    real(wp), parameter :: C2_7 = TWO/7.e0_wp
    real(wp), parameter :: C2_9 = TWO/9.e0_wp
    real(wp), parameter :: C2_15 = TWO/15.e0_wp
    real(wp), parameter :: C2_27 = TWO/27.e0_wp
    real(wp), parameter :: C2_45 = TWO/45.e0_wp
    real(wp), parameter :: C3_2 = THREE/TWO
    real(wp), parameter :: C3_4 = THREE/FOUR
    real(wp), parameter :: C3_7 = THREE/7.e0_wp
    real(wp), parameter :: C3_8 = THREE/8.e0_wp
    real(wp), parameter :: C3_16 = THREE/16.e0_wp
    real(wp), parameter :: C4_3 = FOUR/THREE
    real(wp), parameter :: C4_5 = FOUR/FIVE
    real(wp), parameter :: C4_7 = FOUR/7.e0_wp
    real(wp), parameter :: C4_15 = FOUR/15.e0_wp
    real(wp), parameter :: C5_6 = FIVE/6.e0_wp
    real(wp), parameter :: C5_12 = FIVE/12.e0_wp
    real(wp), parameter :: C5_36 = FIVE/36.e0_wp
    real(wp), parameter :: C6_7 = 6.e0_wp/7.e0_wp
    real(wp), parameter :: C7_90 = 7.e0_wp/90.e0_wp
    real(wp), parameter :: C8_7 = 8.e0_wp/7.e0_wp
    real(wp), parameter :: C8_9 = 8.e0_wp/9.e0_wp
    real(wp), parameter :: C8_75 = 8.e0_wp/75.e0_wp
    real(wp), parameter :: C12_7 = 12.e0_wp/7.e0_wp
    real(wp), parameter :: C12_90 = 12.e0_wp/90.e0_wp
    real(wp), parameter :: C15_4 = 15.e0_wp/FOUR
    real(wp), parameter :: C16_45 = 16.e0_wp/45.e0_wp
    real(wp), parameter :: C16_90 = 16.e0_wp/90.e0_wp
    real(wp), parameter :: C32_90 = 32.e0_wp/90.e0_wp
    real(wp), parameter :: SQ3_6 = sqrt(3.e0_wp)/6.e0_wp
    real(wp), parameter :: SQ15_5 = sqrt(15.e0_wp)/5.e0_wp
    real(wp), parameter :: SQ15_15 = sqrt(15.e0_wp)/15.e0_wp
    real(wp), parameter :: SQ15_24 = sqrt(15.e0_wp)/24.e0_wp
    real(wp), parameter :: SQ15_30 = sqrt(15.e0_wp)/30.e0_wp

    real(wp), parameter :: SAFE_LOW = 1e-30_wp

    ! BASIC CONFIGURATION INPUT
    integer(kind=4) :: NDIM = 1  ! Number of dimensions
    real(wp) :: E_TOL = 1.e-14_wp ! Error tolerance
    real(wp) :: BETA = 0.15e0_wp  ! Learning rate
    real(wp) :: DT_MIN = 1.e-6_wp  ! Minimum dt
    integer(kind=4) :: EXTRA = 0  ! Amount of extra variables (without der) that are not POS
    !! derived
    integer(kind=4) :: NDIM2 = 2
    integer(kind=4) :: EXTRA2 = 0
    real(wp) :: DT_MIN_NOW = 1.e-6_wp  ! Minimum dt can be lower if dt is low

    ! Global
    real(wp), allocatable :: der(:) ! Used in callers

    abstract interface

    !!!! TEMPLATES

        ! Keep going template
        function function_check_keep_tem(y) result(keep_going)
            import :: wp
            implicit none
            real(wp), dimension(:), intent(in) :: y
            logical :: keep_going
        end function function_check_keep_tem

        ! f__ (t, y__) = (f_i (t, y__), ..., f_n (t, y__)) = der__
        function dydt_tem(t, y) result(derivate)
            import :: wp
            implicit none
            real(wp), intent(in)               :: t
            real(wp), dimension(:), intent(in) :: y
            real(wp), dimension(size(y))      :: derivate
        end function dydt_tem

    end interface

contains

    pure function get_index(i) result(idx)
        implicit none
        integer(kind=4), intent(in) :: i
        integer(kind=4) :: idx
        idx = EXTRA2 + NDIM2*(i - 1) + 1
    end function get_index

end module shared
