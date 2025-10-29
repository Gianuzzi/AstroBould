!> Module with different integrators.
module integrators
    implicit none
    !! Changes WILL be made, following:
    !!! - https://github.com/jacobwilliams/rklib.git
    !!! - https://github.com/jacobwilliams/Fortran-Astrodynamics-Toolkit.git

    ! Define dydt_i pointer
    type dydt_i
        procedure(dydt_i_tem), pointer, nopass :: f_i => null ()
    end type dydt_i

    ! For adaptive step and implicit (might be overwritten)
    integer(kind=4)         :: MAX_N_ITER = 100
    real(kind=8), parameter :: MAX_DT_FAC = 5.0d0, SAFE_LOW = 1.d-40
    real(kind=8)            :: E_TOL = 1.d-13

    ! Aux Constants
    real(kind=8), parameter :: ZERO   = 0.0d0
    real(kind=8), parameter :: ONE    = 1.0d0
    real(kind=8), parameter :: TWO    = 2.0d0
    real(kind=8), parameter :: THREE  = 3.0d0
    real(kind=8), parameter :: FOUR   = 4.0d0
    real(kind=8), parameter :: FIVE   = 5.0d0
    real(kind=8), parameter :: C1_2   = ONE/TWO, C3_2 = C1_2 * THREE
    real(kind=8), parameter :: C1_3   = ONE/THREE, C2_3 = C1_3 * TWO
    real(kind=8), parameter :: C1_4   = ONE/FOUR, C3_4 = C1_4 * THREE
    real(kind=8), parameter :: C1_5   = ONE/FIVE, C2_5 = C1_5 * TWO
    real(kind=8), parameter :: C1_6   = ONE/6.d0, C5_6 = C1_6 * FIVE
    real(kind=8), parameter :: C1_7   = ONE/7.d0
    real(kind=8), parameter :: C1_8   = ONE/8.d0, C3_8 = C1_8 * THREE
    real(kind=8), parameter :: C1_9   = ONE/9.d0, C2_9 = C1_9 * TWO, C4_9 = C1_9 * FOUR, C8_9 = C1_9 * 8.d0
    real(kind=8), parameter :: C1_12  = ONE/12.d0, C5_12 = C1_12 * FIVE
    real(kind=8), parameter :: C1_15  = ONE/15.d0, C4_15 = C1_15 * FOUR
    real(kind=8), parameter :: C1_18  = ONE/18.d0, C5_18 = C1_18 * FIVE
    real(kind=8), parameter :: C2_27  = TWO/27.d0
    real(kind=8), parameter :: C3_16  = THREE/16.d0
    real(kind=8), parameter :: C5_36  = FIVE/36.d0
    real(kind=8), parameter :: C2_45  = TWO/45.d0
    real(kind=8), parameter :: C1_48  = ONE/48.d0
    real(kind=8), parameter :: C1_840 = ONE/84.d1, C41_840 = C1_840 * 41.d0
    real(kind=8), parameter :: SQ3    = sqrt(THREE)
    real(kind=8), parameter :: SQ15   = sqrt(15.d0)
    real(kind=8), parameter :: SQ3_6  = SQ3 * C1_6, SQ3_2 = SQ3_6 * TWO
    real(kind=8), parameter :: SQ15_5 = SQ15 * C1_5, SQ15_10 = SQ15_5 * TWO, SQ15_24 = SQ15/24.d0
    real(kind=8), parameter :: SQ15_30 = SQ15/3.d1, SQ15_15 = SQ15_30 * TWO

    abstract interface
    
        function function_check_keep_tem (y) result(keep_going)
            implicit none
            real(kind=8), dimension(:), intent(in) :: y
            logical :: keep_going
        end function function_check_keep_tem

        !---------------------------------------------------------------------------------------------
        ! ND -> ND
        !---------------------------------------------------------------------------------------------

        ! f_i (t, y__) = der
        real(kind=8) function dydt_i_tem (t, y) result (der)
            implicit none
            real(kind=8), intent(in)               :: t
            real(kind=8), dimension(:), intent(in) :: y
        end function dydt_i_tem

        ! f__ (t, y__) = (f_i (t, y__), ..., f_n (t, y__)) = der__
        function dydt_tem (t, y) result (der)
            implicit none
            real(kind=8), intent(in)               :: t
            real(kind=8), dimension(:), intent(in) :: y
            real(kind=8), dimension(size (y))      :: der
            ! Here must be every f_i defined explicitly
        end function dydt_tem

        ! in (t, y__, dt, f__, ynew__) -> ynew__
        ! Remember that, in this case,
        !  f__ == (f_1, ..., f_N) must be
        !  pre-defined explicitly
        subroutine integ_tem (t, y, der, dt, dydt, ynew)
            import :: dydt_tem
            implicit none
            real(kind=8), intent(in)                       :: t, dt
            real(kind=8), dimension(:), intent(in)         :: y
            real(kind=8), dimension(size (y)), intent(in)  :: der
            procedure(dydt_tem)                            :: dydt
            real(kind=8), dimension(size (y)), intent(out) :: ynew
        end subroutine integ_tem
        
        ! in (t, y__, dt, f__,  osol, oerr, yaux__, ynew__) -> osol, oerr, yaux__, ynew__
        ! Remember that, in this case,
        !  f__ == (f_1, ..., f_N) must be
        !  pre-defined explicitly
        subroutine embedded_tem (t, y, der, dt, dydt, osol, oaux, yaux, ynew)
            import :: dydt_tem
            implicit none
            real(kind=8), intent(in)                       :: t, dt
            real(kind=8), dimension(:), intent(in)         :: y
            real(kind=8), dimension(size (y)), intent(in)  :: der
            procedure(dydt_tem)                            :: dydt
            integer(kind=4), intent(out)                   :: osol, oaux
            real(kind=8), dimension(size (y)), intent(out) :: ynew, yaux
        end subroutine embedded_tem

        subroutine leapfrog_tem (t, y, der, dt, dydt, ynew)
            import :: dydt_tem
            implicit none
            real(kind=8), intent(in)                       :: t, dt
            real(kind=8), dimension(:), intent(in)         :: y
            real(kind=8), dimension(size (y)), intent(in)  :: der
            procedure(dydt_tem)                            :: dydt
            real(kind=8), dimension(size (y)), intent(out) :: ynew
        end subroutine leapfrog_tem

        subroutine calc_rk (t, y, dt, dydt, kin, kout)
            import :: dydt_tem
            implicit none
            real(kind=8), intent(in)                                          :: t, dt
            procedure(dydt_tem)                                               :: dydt
            real(kind=8), dimension(:, :), intent(in)                         :: kin
            real(kind=8), dimension(size (kin, 1)), intent(in)                :: y
            real(kind=8), dimension(size (kin, 1), size (kin,2)), intent(out) :: kout
        end subroutine calc_rk

    end interface

    contains
    
        !---------------------------------------------------------------------------------------------
        ! ND -> ND
        !---------------------------------------------------------------------------------------------

        !------------------------------------------ SOLVERS ------------------------------------------

        !! Implicit methods solvers

        subroutine solve_1k_implicit (t, y, dt, dydt, kprev, cte, k)
            implicit none
            real(kind=8), intent(in)                       :: t, dt, cte
            real(kind=8), dimension(:), intent(in)         :: y
            real(kind=8), dimension(size (y)), intent(in)  :: kprev
            real(kind=8), dimension(size (y)), intent(out) :: k
            real(kind=8), dimension(size (y))              :: kaux
            procedure(dydt_tem)                            :: dydt
            integer(kind=4)                                :: i
            
            k = dydt (t, y + dt * kprev)
            do i = 1, MAX_N_ITER
                kaux = k
                k    = dydt (t, y + dt * (kprev + cte * k))
                if (maxval( abs ((kaux - k) / (kaux + SAFE_LOW))) .le. E_TOL) then                    
                    exit
                end if
            end do
        end subroutine solve_1k_implicit

        subroutine solve_rk_implicit (t, y, dt, dydt, solver, rk)
            implicit none
            real(kind=8), intent(in)                          :: t, dt
            procedure(dydt_tem)                               :: dydt
            real(kind=8), dimension(:, :), intent(inout)      :: rk
            real(kind=8), dimension(size (rk,1), size (rk,2)) :: rkold
            real(kind=8), dimension(size (rk,1)), intent(in)  :: y
            procedure(calc_rk)                                :: solver
            integer(kind=4)                                   :: i

            do i = 1, MAX_N_ITER
                rkold = rk
                call solver (t, y, dt, dydt, rkold, rk)
                if (maxval( abs ((rkold - rk) / (rkold + SAFE_LOW))) .le. E_TOL) then
                    exit
                end if
            end do
        end subroutine solve_rk_implicit

        !! Runge Kutta methods solvers

        subroutine get_rks (t, y, dt, dydt, m, rk) ! Unused in the current version
            implicit none
            real(kind=8), intent(in)                                         :: t, dt
            real(kind=8), dimension(:), intent(in)                           :: y
            real(kind=8), dimension(:,:), intent(in)                         :: m
            real(kind=8), dimension((size (m,1) - 1), size (y)), intent(out) :: rk ! In columns
            real(kind=8), dimension(size (y))                                :: rkaux, added
            procedure(dydt_tem)                                              :: dydt
            integer(kind=4)                                                  :: i, j
            
            do i = 1, size (m, 1) - 1 ! Rows
                rkaux = ZERO
                do j = 1, i ! Cols (<i bc its inf triang)
                    added = m(i+1,j) * rk(j,:)
                    rkaux = rkaux + added
                end do
                rk(i, :) = dydt (t + m(1,i) * dt, y + dt * rkaux)
            end do
        end subroutine get_rks
        
        subroutine solve_rk (t, y, dt, dydt, m, ynew) ! Unused in the current version
            implicit none
            real(kind=8), intent(in)                            :: t, dt
            real(kind=8), dimension(:), intent(in)              :: y
            real(kind=8), dimension(:,:), intent(in)            :: m
            real(kind=8), dimension((size (m,1) - 1), size (y)) :: rk
            procedure(dydt_tem)                                 :: dydt
            real(kind=8), dimension(size (y)), intent(out)      :: ynew
            integer(kind=4)                                     :: i 
            
            call get_rks (t, y, dt, dydt, m, rk)
            
            do i = 1, size (y)
                ynew(i) = y(i) + dt * dot_product ((/m(2:, size (m, 1))/), rk(:,i))
            end do
        end subroutine solve_rk

        subroutine solve_rk_embed (t, y, dt, dydt, m, maux, yaux, ynew) ! Unused in the current version
            implicit none
            real(kind=8), intent(in)                              :: t, dt
            real(kind=8), dimension(:), intent(in)                :: y
            real(kind=8), dimension(:,:), intent(in)              :: m
            real(kind=8), dimension((size (m,1) - 1), size (y))   :: rk
            real(kind=8), dimension((size (m,1) - 1)), intent(in) :: maux
            procedure(dydt_tem)                                   :: dydt
            real(kind=8), dimension(size (y)), intent(out)        :: yaux, ynew
            integer(kind=4)                                       :: i 
            
            call get_rks (t, y, dt, dydt, m, rk)
            
            do i = 1, size (y)
                yaux(i) = y(i) + dt * dot_product (maux, rk(:,i))
                ynew(i) = y(i) + dt * dot_product ((/m(2:, size (m, 1))/), rk(:,i))
            end do
        end subroutine solve_rk_embed

        !! Embedded methods solver

        recursive subroutine solve_embed (t, y, der, dt_adap, dydt, integ, e_tol, beta, dt_min, dt_used, ynew)
            implicit none
            real(kind=8), intent(in)                            :: t, e_tol, beta, dt_min
            real(kind=8), intent(inout)                         :: dt_adap, dt_used
            real(kind=8), dimension(:), intent(in)              :: y
            real(kind=8), dimension(size (y)), intent(in)       :: der
            integer(kind=4), save                               :: osol, oaux, iter = 0
            real(kind=8), dimension(size (y))                   :: yaux, yscal
            procedure(dydt_tem)                                 :: dydt
            procedure(embedded_tem)                             :: integ
            real(kind=8), dimension(size (y)), intent(out)      :: ynew            
            real(kind=8)                                        :: e_calc, ratio

            iter    = iter + 1
            dt_adap = max (dt_adap, dt_min)            
            
            call integ (t, y, der, dt_adap, dydt, osol, oaux, yaux, ynew)

            yscal = abs (y) + abs (dt_adap * der) + SAFE_LOW
            
            e_calc = max (maxval (abs ((ynew - yaux) / yscal )), SAFE_LOW)
            ratio  = e_tol / e_calc
            if (ratio > ONE) then
                dt_used = dt_adap
                dt_adap = dt_adap * min (beta * ratio**(ONE / osol), MAX_DT_FAC)
                iter = 0
            else
                if (abs(dt_adap - dt_min) .le. E_TOL) then !E_TOL?
                    dt_used = dt_min
                    iter = 0
                else
                    dt_adap = dt_adap * min (beta * ratio**(ONE / oaux), MAX_DT_FAC)
                    if ((dt_adap .ne. dt_adap) .or. (dt_adap < dt_min) .or. (iter .eq. MAX_N_ITER)) then
                        dt_adap = dt_min
                        dt_used = dt_min
                        call integ (t, y, der, dt_adap, dydt, osol, oaux, yaux, ynew)
                        iter = 0
                    else
                        call solve_embed (t, y, der, dt_adap, dydt, integ, e_tol, beta, dt_min, dt_used, ynew)
                    end if 
                end if
            end if
        end subroutine solve_embed

        !! Runge Kutta half_step solver

        recursive subroutine solve_rk_half_step (t, y, der, dt_adap, dydt, integ, ord, e_tol, beta, dt_min, dt_used, ynew)
            implicit none
            integer(kind=4), intent(in)                    :: ord
            real(kind=8), intent(in)                       :: t, e_tol, beta, dt_min
            procedure(integ_tem)                           :: integ
            real(kind=8), dimension(:), intent(in)         :: y
            real(kind=8), dimension(size (y)), intent(in)  :: der
            real(kind=8), dimension(size (y))              :: derhalf
            real(kind=8), dimension(size (y))              :: yhalf, yaux, yscal
            procedure(dydt_tem)                            :: dydt
            real(kind=8), dimension(size (y)), intent(out) :: ynew
            real(kind=8), intent(inout)                    :: dt_adap, dt_used
            real(kind=8)                                   :: e_calc, ratio, hdt_adap
            integer(kind=4)                                :: iter = 0
            real(kind=8), save                             :: twotordmi1, ordplus1, ordr

            ordr       = ord * ONE
            ordplus1   = ordr + ONE
            twotordmi1 = TWO**ord - ONE
            iter       = iter + 1
            dt_adap    = max (dt_adap, dt_min)
            hdt_adap   = C1_2 * dt_adap

            call integ (t, y, der,  dt_adap, dydt,  ynew)
            call integ (t, y, der, hdt_adap, dydt, yhalf)
            derhalf = dydt (t + hdt_adap, yhalf)
            call integ (t + hdt_adap, yhalf, derhalf, hdt_adap, dydt,  yaux)

            yscal = abs (y) + abs (dt_adap * der) + SAFE_LOW
            
            e_calc = max (maxval (abs ((ynew - yaux) / yscal) / twotordmi1), SAFE_LOW)
            ratio  = e_tol / e_calc
            if (ratio > ONE) then
                dt_used = dt_adap
                dt_adap = dt_adap * min (beta * ratio**(ONE / ordplus1), MAX_DT_FAC)
                iter    = 0
            else
                if (abs(dt_adap-dt_min) .le. E_TOL) then !E_TOL?
                    dt_used = dt_min
                    iter = 0
                else
                    dt_adap = dt_adap * min (beta * ratio**(ONE / ordr), MAX_DT_FAC)
                    if ((dt_adap .ne. dt_adap) .or. (dt_adap .le. dt_min) .or. (iter .eq. MAX_N_ITER)) then
                        dt_used = dt_min
                        dt_adap = dt_min
                        call integ (t, y, der, dt_adap, dydt, ynew)
                        iter = 0
                    else
                        call solve_rk_half_step (t, y, der, dt_adap, dydt, integ, ord, e_tol, beta, dt_min, dt_used, ynew)
                    end if
                end if
            end if
        end subroutine solve_rk_half_step

        !! Adaptive Leap Frog Solver

        recursive subroutine solve_leapfrog (t, y, der, dt_adap, dydt, leapfrog, e_tol, beta, dt_min, dt_used, ynew)
            implicit none
            real(kind=8), intent(in)                            :: t, e_tol, beta, dt_min
            real(kind=8), intent(inout)                         :: dt_adap, dt_used
            real(kind=8), dimension(:), intent(in)              :: y
            real(kind=8), dimension(size (y)), intent(in)       :: der
            integer(kind=4), save                               :: iter = 0
            real(kind=8), dimension(size (y))                   :: derhalf, yhalf, yaux, yscal
            procedure(dydt_tem)                                 :: dydt
            procedure(leapfrog_tem)                             :: leapfrog
            real(kind=8), dimension(size (y)), intent(out)      :: ynew            
            real(kind=8)                                        :: e_calc, ratio, dt_half

            iter    = iter + 1
            dt_adap = max (dt_adap, dt_min)
            dt_half = C1_2 * dt_adap
            
            call leapfrog (t, y, der,  dt_adap, dydt, ynew)
            call leapfrog (t, y, der, dt_half, dydt, yhalf)
            derhalf = dydt (t + dt_half, yhalf)
            call leapfrog (t + dt_half, yhalf, derhalf, dt_half, dydt, yaux)

            yscal = abs (y) + abs (dt_adap * der) + SAFE_LOW
            
            e_calc = max (maxval (abs ((ynew - yaux) / yscal)), SAFE_LOW)
            ratio  = e_tol / e_calc

            if (ratio > ONE) then
                dt_used = dt_adap
                dt_adap = dt_adap * min (beta * ratio**C1_3, MAX_DT_FAC) ! 1/3 porque es 1/(O(2) + 1)
                iter    = 0
            else
                if (abs (dt_adap - dt_min) .le. E_TOL) then !E_TOL?
                    dt_used = dt_min
                    iter = 0
                else
                    dt_adap = dt_adap * min (beta * ratio**C1_2, MAX_DT_FAC)
                    if ((dt_adap .ne. dt_adap) .or. (dt_adap .le. dt_min) .or. (iter .eq. MAX_N_ITER)) then
                        dt_used = dt_min
                        dt_adap = dt_min
                        call leapfrog (t, y, der, dt_adap, dydt, ynew)
                        iter = 0
                    else
                        call solve_leapfrog (t, y, der, dt_adap, dydt, leapfrog, e_tol, beta, dt_min, dt_used, ynew)
                    end if
                end if
            end if
        end subroutine solve_leapfrog


        !---------------------------------------- INTEGRATORS ----------------------------------------

        !! Runge Kutta methods (implicit and explicit)
        
        subroutine Euler1 (t, y, der, dt, dydt, ynew)
            implicit none
            real(kind=8), intent(in)                       :: t, dt
            real(kind=8), dimension(:), intent(in)         :: y
            real(kind=8), dimension(size (y)), intent(in)  :: der
            procedure(dydt_tem)                            :: dydt
            real(kind=8), dimension(size (y)), intent(out) :: ynew
            
            ! k1 = dydt (t, y)
            ynew = y + dt * der
        end subroutine Euler1
        
        subroutine Euler_back1 (t, y, der, dt, dydt, ynew) ! Implicit
            implicit none
            real(kind=8), intent(in)                       :: t, dt
            real(kind=8), dimension(:), intent(in)         :: y
            real(kind=8), dimension(size (y)), intent(in)  :: der ! Unused as it is implicit
            procedure(dydt_tem)                            :: dydt
            real(kind=8), dimension(size (y)), intent(out) :: ynew
            real(kind=8), dimension(size (y))              :: kaux, k1
            real(kind=8)                                   :: aux

            kaux = ZERO
            aux  = ONE
            call solve_1k_implicit (t + dt, y, dt, dydt, kaux, aux, k1)

            ynew = y + dt * k1
        end subroutine Euler_back1

        subroutine Euler_center2 (t, y, der, dt, dydt, ynew) ! Implicit
            implicit none
            real(kind=8), intent(in)                       :: t, dt
            real(kind=8), dimension(:), intent(in)         :: y
            real(kind=8), dimension(size (y)), intent(in)  :: der ! Unused as it is implicit
            procedure(dydt_tem)                            :: dydt
            real(kind=8), dimension(size (y)), intent(out) :: ynew
            real(kind=8), dimension(size (y))              :: kaux, k1
            real(kind=8)                                   :: aux

            kaux = ZERO
            aux  = C1_2
            call solve_1k_implicit (t + dt * C1_2, y, dt, dydt, kaux, aux, k1)

            ynew = y + dt * k1
        end subroutine Euler_center2
        
        subroutine Crank_Nicolson2 (t, y, der, dt, dydt, ynew) ! Implicit
            implicit none
            real(kind=8), intent(in)                       :: t, dt
            real(kind=8), dimension(:), intent(in)         :: y
            real(kind=8), dimension(size (y)), intent(in)  :: der
            procedure(dydt_tem)                            :: dydt
            real(kind=8), dimension(size (y)), intent(out) :: ynew
            real(kind=8), dimension(size (y))              :: kaux, k2
            real(kind=8)                                   :: aux

            ! k1   = dydt (t, y)
            kaux = der * C1_2
            aux  = C1_2
            call solve_1k_implicit (t + dt, y, dt, dydt, kaux, aux, k2)

            ynew = y + dt * (der + k2) * C1_2
        end subroutine Crank_Nicolson2

        subroutine Heun2 (t, y, der, dt, dydt, ynew)
            implicit none
            real(kind=8), intent(in)                       :: t, dt
            real(kind=8), dimension(:), intent(in)         :: y
            real(kind=8), dimension(size (y)), intent(in)  :: der
            procedure(dydt_tem)                            :: dydt
            real(kind=8), dimension(size (y)), intent(out) :: ynew              
            real(kind=8), dimension(size (y))              :: k2
            
            ! k1 = dydt (t,      y)
            k2 = dydt (t + dt, y + dt * der)
            
            ynew = y + dt * (der + k2) * C1_2
        end subroutine Heun2

        subroutine midpoint2 (t, y, der, dt, dydt, ynew)
            implicit none
            real(kind=8), intent(in)                       :: t, dt
            real(kind=8), dimension(:), intent(in)         :: y
            real(kind=8), dimension(size (y)), intent(in)  :: der
            procedure(dydt_tem)                            :: dydt
            real(kind=8), dimension(size (y)), intent(out) :: ynew
            real(kind=8), dimension(size (y))              :: k2
            
            ! k1 = dydt (t,            y)
            k2 = dydt (t + dt * C1_2, y + dt * der * C1_2)
            
            ynew = y + dt * k2
        end subroutine midpoint2

        subroutine strange2 (t, y, der, dt, dydt, ynew)
            implicit none
            real(kind=8), intent(in)                       :: t, dt
            real(kind=8), dimension(:), intent(in)         :: y
            real(kind=8), dimension(size (y)), intent(in)  :: der
            procedure(dydt_tem)                            :: dydt
            real(kind=8), dimension(size (y)), intent(out) :: ynew
            real(kind=8), dimension(size (y))              :: k2
            
            ! k1   = dydt (t,             y)
            k2   = dydt (t + dt * C3_4, y + dt * der * C3_4)
            
            ynew = y + dt * (der + k2 * TWO) * C1_3
        end subroutine strange2

        subroutine Ralston2 (t, y, der, dt, dydt, ynew)
            implicit none
            real(kind=8), intent(in)                       :: t, dt
            real(kind=8), dimension(:), intent(in)         :: y
            real(kind=8), dimension(size (y)), intent(in)  :: der
            procedure(dydt_tem)                            :: dydt
            real(kind=8), dimension(size (y)), intent(out) :: ynew
            real(kind=8), dimension(size (y))              :: k2
            
            ! k1 = dydt (t,             y)
            k2 = dydt (t + dt * C2_3, y + dt * der * C2_3)
            
            ynew = y + dt * (der + k2 * THREE) * C1_4
        end subroutine Ralston2

        subroutine Hammer_Hollingsworth2 (t, y, der, dt, dydt, ynew) ! Implicit
            implicit none
            real(kind=8), intent(in)                       :: t, dt
            real(kind=8), dimension(:), intent(in)         :: y
            real(kind=8), dimension(size (y)), intent(in)  :: der
            procedure(dydt_tem)                            :: dydt
            real(kind=8), dimension(size (y)), intent(out) :: ynew
            real(kind=8), dimension(size (y))              :: kaux, k2
            real(kind=8)                                   :: aux

            ! k1 = dydt (t, y)
            kaux = der * C1_3
            aux  = C1_3
            call solve_1k_implicit (t + dt * C2_3, y, dt, dydt, kaux, aux, k2)

            ynew = y + dt * (der + k2 * THREE) * C1_4
        end subroutine Hammer_Hollingsworth2
        
        subroutine Kraaijevanger_Spijker2 (t, y, der, dt, dydt, ynew) ! Implicit
            implicit none
            real(kind=8), intent(in)                       :: t, dt
            real(kind=8), dimension(:), intent(in)         :: y
            real(kind=8), dimension(size (y)), intent(in)  :: der ! Unused as it is implicit
            procedure(dydt_tem)                            :: dydt
            real(kind=8), dimension(size (y)), intent(out) :: ynew
            real(kind=8), dimension(size (y))              :: kaux, k1, k2
            real(kind=8)                                   :: aux

            kaux = ZERO
            aux  = C1_2
            call solve_1k_implicit (t + dt * C1_2, y, dt, dydt, kaux, aux, k1)
            kaux = - k1 * C1_2
            aux  = TWO
            call solve_1k_implicit (t + dt * C3_2, y, dt, dydt, kaux, aux, k2)

            ynew = y + dt * (- k1 + k2 * THREE) * C1_2
        end subroutine Kraaijevanger_Spijker2
        
        subroutine Qin_Zhang2 (t, y, der, dt, dydt, ynew) ! Implicit
            implicit none
            real(kind=8), intent(in)                       :: t, dt
            real(kind=8), dimension(:), intent(in)         :: y
            real(kind=8), dimension(size (y)), intent(in)  :: der ! Unused as it is implicit
            procedure(dydt_tem)                            :: dydt
            real(kind=8), dimension(size (y)), intent(out) :: ynew
            real(kind=8), dimension(size (y))              :: kaux, k1, k2
            real(kind=8), parameter                        :: aux = C1_4

            kaux = ZERO
            call solve_1k_implicit (t + dt * C1_4, y, dt, dydt, kaux, aux, k1)
            kaux = k1 * C1_2
            call solve_1k_implicit (t + dt * C3_4, y, dt, dydt, kaux, aux, k2)

            ynew = y + dt * (k1 + k2) * C1_2
        end subroutine Qin_Zhang2

        subroutine Runge_Kutta3 (t, y, der, dt, dydt, ynew)
            implicit none
            real(kind=8), intent(in)                       :: t, dt
            real(kind=8), dimension(:), intent(in)         :: y
            real(kind=8), dimension(size (y)), intent(in)  :: der
            procedure(dydt_tem)                            :: dydt
            real(kind=8), dimension(size (y)), intent(out) :: ynew
            real(kind=8), dimension(size (y))              :: k2, k3
            
            ! k1 = dydt (t,             y)
            k2 = dydt (t + dt * C1_2, y + dt * der * C1_2)
            k3 = dydt (t + dt,        y + dt * (- der + k2 * 2))
            
            ynew = y + dt * (der + k2 * FOUR + k3) * C1_6
        end subroutine Runge_Kutta3

        subroutine Heun3 (t, y, der, dt, dydt, ynew)
            implicit none
            real(kind=8), intent(in)                       :: t, dt
            real(kind=8), dimension(:), intent(in)         :: y
            real(kind=8), dimension(size (y)), intent(in)  :: der
            procedure(dydt_tem)                            :: dydt
            real(kind=8), dimension(size (y)), intent(out) :: ynew
            real(kind=8), dimension(size (y))              :: k2, k3
            
            ! k1 = dydt (t,             y)
            k2 = dydt (t + dt * C1_3, y + dt * der * C1_3)
            k3 = dydt (t + dt * C2_3, y + dt * k2 * C2_3)
            
            ynew = y + dt * (der + k3 * THREE) * C1_4
        end subroutine Heun3

        subroutine Ralston3 (t, y, der, dt, dydt, ynew)
            implicit none
            real(kind=8), intent(in)                       :: t, dt
            real(kind=8), dimension(:), intent(in)         :: y
            real(kind=8), dimension(size (y)), intent(in)  :: der
            procedure(dydt_tem)                            :: dydt
            real(kind=8), dimension(size (y)), intent(out) :: ynew
            real(kind=8), dimension(size (y))              :: k2, k3

            ! k1 = dydt (t,             y)
            k2 = dydt (t + dt * C1_2, y + dt * der * C1_2)
            k3 = dydt (t + dt * C3_4, y + dt * k2 * C3_4)

            ynew = y + dt * (der * TWO + k2 * THREE + k3 * FOUR) * C1_9
        end subroutine Ralston3

        subroutine SSPRK3 (t, y, der, dt, dydt, ynew)
            implicit none
            real(kind=8), intent(in)                       :: t, dt
            real(kind=8), dimension(:), intent(in)         :: y
            real(kind=8), dimension(size (y)), intent(in)  :: der
            procedure(dydt_tem)                            :: dydt
            real(kind=8), dimension(size (y)), intent(out) :: ynew
            real(kind=8), dimension(size (y))              :: k2, k3

            ! k1 = dydt (t,            y)
            k2 = dydt (t + dt,       y + dt * der)
            k3 = dydt (t + dt * C1_2, y + dt * (der + k2) * C1_4)

            ynew = y + dt * (der + k2 + k3 * FOUR) * C1_6
        end subroutine SSPRK3

        subroutine Crouzeix3 (t, y, der, dt, dydt, ynew) ! Implicit
            implicit none
            real(kind=8), intent(in)                       :: t, dt
            real(kind=8), dimension(:), intent(in)         :: y
            real(kind=8), dimension(size (y)), intent(in)  :: der ! Unused as it is implicit
            procedure(dydt_tem)                            :: dydt
            real(kind=8), dimension(size (y)), intent(out) :: ynew
            real(kind=8), dimension(size (y))              :: kaux, k1, k2
            real(kind=8)                                   :: aux

            kaux = ZERO
            aux  = C1_2 + SQ3_6
            call solve_1k_implicit (t + dt * aux, y, dt, dydt, kaux, aux, k1)
            kaux = - k1 * SQ3_6 * 2
            call solve_1k_implicit (t + dt * (C1_2 - SQ3_6), y, dt, dydt, kaux, aux, k2)

            ynew = y + dt * (k1 + k2) * C1_2
        end subroutine Crouzeix3

        subroutine Runge_Kutta_implicit3 (t, y, der, dt, dydt, ynew) ! Implicit
            implicit none
            real(kind=8), intent(in)                       :: t, dt
            real(kind=8), dimension(:), intent(in)         :: y
            real(kind=8), dimension(size (y)), intent(in)  :: der ! Unused as it is implicit
            procedure(dydt_tem)                            :: dydt
            real(kind=8), dimension(size (y)), intent(out) :: ynew
            real(kind=8), dimension(size (y))              :: kaux, k1, k2, k3, k4
            real(kind=8), parameter                        :: aux = C1_2

            kaux = ZERO
            call solve_1k_implicit (t + dt * C1_2, y, dt, dydt, kaux, aux, k1)
            kaux = k1 * C1_6
            call solve_1k_implicit (t + dt * C2_3, y, dt, dydt, kaux, aux, k2)
            kaux = (- k1 + k2) * C1_2
            call solve_1k_implicit (t + dt * C1_2, y, dt, dydt, kaux, aux, k3)
            kaux = ((k1 - k2) * 3 + k3) * C1_2
            call solve_1k_implicit (t + dt, y, dt, dydt, kaux, aux, k4)

            ynew = y + dt * (kaux + k4 * C1_2)
        end subroutine Runge_Kutta_implicit3

        subroutine Ralston4 (t, y, der, dt, dydt, ynew)
            implicit none
            real(kind=8), intent(in)                       :: t, dt
            real(kind=8), dimension(:), intent(in)         :: y
            real(kind=8), dimension(size (y)), intent(in)  :: der
            procedure(dydt_tem)                            :: dydt
            real(kind=8), dimension(size (y)), intent(out) :: ynew
            real(kind=8), dimension(size (y))              :: k2, k3, k4

            ! k1 = dydt (t,                     y)
            k2 = dydt (t + dt * C2_5,         y + dt * der * C2_5)
            k3 = dydt (t + dt * 0.45573725d0, y + dt * (der * 0.29697761d0 + k2 * 0.15875964d0))
            k4 = dydt (t + dt,                y + dt * (der * 0.21810040d0 - k2 * 3.05096516d0 + k3 * 3.83286476d0))

            ynew = y + dt * (der * 0.17476028d0 - k2 * 0.55148066d0 + k3 * 1.20553560d0 + k4 * 0.17118478d0)
        end subroutine Ralston4

        ! subroutine Lobatto4 (t, y, der, dt, dydt, ynew) ! Implicit
        !     implicit none
        !     real(kind=8), intent(in)                       :: t, dt
        !     real(kind=8), dimension(:), intent(in)         :: y
        !     real(kind=8), dimension(size (y)), intent(in)  :: der
        !     procedure(dydt_tem)                            :: dydt
        !     real(kind=8), dimension(size (y)), intent(out) :: ynew
        !     real(kind=8), dimension(size (y))              :: kaux, k2, k3
        !     real(kind=8), parameter                        :: aux = C1_4

        !     ! k1   = dydt (t, y)
        !     kaux = der * C1_4
        !     call solve_1k_implicit (t + dt * C1_2, y, dt, dydt, kaux, aux, k2)
        !     k3 = dydt (t + dt, y + dt * k2)

        !     ynew = y + dt * (der + k2 * FOUR + k3) * C1_6
        ! end subroutine Lobatto4

        subroutine Runge_Kutta4 (t, y, der, dt, dydt, ynew)
            implicit none
            real(kind=8), intent(in)                       :: t, dt
            real(kind=8), dimension(:), intent(in)         :: y
            real(kind=8), dimension(size (y)), intent(in)  :: der
            procedure(dydt_tem)                            :: dydt
            real(kind=8), dimension(size (y)), intent(out) :: ynew
            real(kind=8), dimension(size (y))              :: k2, k3, k4
            real(kind=8)                                   :: aut
              
            aut = dt * C1_2

            ! k1 = dydt (t,       y)
            k2 = dydt (t + aut, y + dt * der * C1_2)
            k3 = dydt (t + aut, y + dt * k2 * C1_2)
            k4 = dydt (t + dt,  y + dt * k3)
        
            ynew = y + dt * (der + (k2 + k3) * TWO + k4) * C1_6
        end subroutine Runge_Kutta4

        subroutine Gauss_Legendre4 (t, y, der, dt, dydt, ynew) ! Implicit
            implicit none
            real(kind=8), intent(in)                       :: t, dt
            real(kind=8), dimension(:), intent(in)         :: y
            real(kind=8), dimension(size (y)), intent(in)  :: der ! Unused as it is implicit
            procedure(dydt_tem)                            :: dydt
            real(kind=8), dimension(size (y)), intent(out) :: ynew
            real(kind=8), dimension(size (y), 2)           :: rk

            rk = ZERO
            call solve_rk_implicit (t, y, dt, dydt, FunK_GL4, rk)

            ynew = y + dt * (rk(:,1) + rk(:,2)) * C1_2

            contains
                subroutine FunK_GL4 (t, y, dt, dydt, kin, kout) !! Funk for Gauss_Legendre4
                    implicit none
                    real(kind=8), intent(in)                                          :: t, dt
                    procedure(dydt_tem)                                               :: dydt
                    real(kind=8), dimension(:, :), intent(in)                         :: kin
                    real(kind=8), dimension(size (kin, 1)), intent(in)                :: y
                    real(kind=8), dimension(size (kin, 1), size (kin,2)), intent(out) :: kout
        
                    kout(:,1) = dydt (t + dt * (C1_2 - SQ3_6), y + dt * ( &
                        & kin(:,1) * C1_4 + &
                        & kin(:,2) * (C1_4 - SQ3_6)))
                    kout(:,2) = dydt (t + dt * (C1_2 + SQ3_6), y + dt * ( &
                        & kin(:,1) * (C1_4 + SQ3_6) + &
                        & kin(:,2) * C1_4))
                end subroutine FunK_GL4
        end subroutine Gauss_Legendre4 

        subroutine Runge_Kutta_four_oct4 (t, y, der, dt, dydt, ynew)
            implicit none
            real(kind=8), intent(in)                       :: t, dt
            real(kind=8), dimension(:), intent(in)         :: y
            real(kind=8), dimension(size (y)), intent(in)  :: der
            procedure(dydt_tem)                            :: dydt
            real(kind=8), dimension(size (y)), intent(out) :: ynew
            real(kind=8), dimension(size (y))              :: k2, k3, k4, kaux
              
            ! k1   = dydt (t,             y)
            kaux = der * C1_3
            k2   = dydt (t + dt * C1_3, y + dt * kaux)
            k3   = dydt (t + dt * C2_3, y + dt * (- kaux + k2))
            k4   = dydt (t + dt,        y + dt * (der - k2 + k3))
        
            ynew = y + dt * (der + (k2 + k3) * THREE + k4) * C1_8  
        end subroutine Runge_Kutta_four_oct4

        subroutine Runge_Kutta5 (t, y, der, dt, dydt, ynew)
            implicit none
            real(kind=8), intent(in)                       :: t, dt
            real(kind=8), dimension(:), intent(in)         :: y
            real(kind=8), dimension(size (y)), intent(in)  :: der
            procedure(dydt_tem)                            :: dydt
            real(kind=8), dimension(size (y)), intent(out) :: ynew
            real(kind=8), dimension(size (y))              :: k2, k3, k4, k5, k6
            real(kind=8)                                   :: aut
              
            aut = dt * C1_4
              
            ! k1 = dydt (t,             y)
            k2 = dydt (t + aut,       y + dt * der * C1_4)
            k3 = dydt (t + aut,       y + dt * (der + k2) * C1_8)
            k4 = dydt (t + dt * C1_2, y + dt * (- k2 * C1_2 + k3))
            k5 = dydt (t + dt * C3_4, y + dt * (der + k4 * THREE) * C3_16)
            k6 = dydt (t + dt,        y + dt * (- der * THREE + k2 * TWO + (k3 - k4) * 12.d0 + k5 * 8.d0) * C1_7)
        
            ynew = y + dt * ((der + k6) * 7.d0 + (k3 + k5) * 32.d0 + k4 * 12.d0)/9.d1
        end subroutine Runge_Kutta5

        subroutine Gauss_Legendre6 (t, y, der, dt, dydt, ynew) ! Implicit
            implicit none
            real(kind=8), intent(in)                       :: t, dt
            real(kind=8), dimension(:), intent(in)         :: y
            real(kind=8), dimension(size (y)), intent(in)  :: der ! Unused as it is implicit
            procedure(dydt_tem)                            :: dydt
            real(kind=8), dimension(size (y)), intent(out) :: ynew
            real(kind=8), dimension(size (y), 3)           :: rk

            rk = ZERO
            call solve_rk_implicit (t, y, dt, dydt, FunK_GL6, rk)

            ynew = y + dt * ((rk(:,1) + rk(:,3)) * FIVE + rk(:,2) * 8.d0) * C1_18

            contains
                subroutine FunK_GL6 (t, y, dt, dydt, kin, kout)  !!! Funk for Gauss_Legendre6
                    implicit none
                    real(kind=8), intent(in)                                          :: t, dt
                    procedure(dydt_tem)                                               :: dydt
                    real(kind=8), dimension(:, :), intent(in)                         :: kin
                    real(kind=8), dimension(size (kin, 1)), intent(in)                :: y
                    real(kind=8), dimension(size (kin, 1), size (kin,2)), intent(out) :: kout
        
                    kout(:,1) = dydt (t + dt * (ONE - SQ15_5) * C1_2,  y + dt * (&
                            & kin(:,1) * C5_36 + &
                            & kin(:,2) * (C2_3 - SQ15_5) * C1_3 + &
                            & kin(:,3) * (C5_36 - SQ15_30)))
                    kout(:,2) = dydt (t + dt * C1_2, y + dt * (&
                            & kin(:,1) * (C5_36 + SQ15_24)+ &
                            & kin(:,2) * C2_9 + &
                            & kin(:,3) * (C5_36 - SQ15_24)))
                    kout(:,3) = dydt (t + dt * (ONE + SQ15_5) * C1_2, y + dt * (&
                            & kin(:,1) * (C5_36 + SQ15_30) + &
                            & kin(:,2) * (C2_9 + SQ15_15) + &
                            & kin(:,3) * C5_36))
                end subroutine FunK_GL6
        end subroutine Gauss_Legendre6 

        subroutine Runge_Kutta6 (t, y, der, dt, dydt, ynew)
            implicit none
            real(kind=8), intent(in)                       :: t, dt
            real(kind=8), dimension(:), intent(in)         :: y
            real(kind=8), dimension(size (y)), intent(in)  :: der
            procedure(dydt_tem)                            :: dydt
            real(kind=8), dimension(size (y)), intent(out) :: ynew
            real(kind=8), dimension(size (y))              :: k2, k3, k4, k5, k6, k7

            ! k1 = dydt(t, y)
            k2 = dydt(t + dt * C1_3, y + dt * C1_3 * der)
            k3 = dydt(t + dt * C1_3, y + dt * C1_3 * k2)
            k4 = dydt(t + dt * C1_2, y + dt * C1_2 * k3)
            k5 = dydt(t + dt * C2_3, y + dt * C2_3 * k4)
            k6 = dydt(t + dt, y + dt * k5)
            k7 = dydt(t + dt, y + dt * k6)

            ynew = y + dt * (C1_12 * der + C1_3 * k4 + C1_3 * k5 + C1_12 * k7)
        end subroutine Runge_Kutta6

        subroutine Abbas6 (t, y, der, dt, dydt, ynew)
            implicit none
            real(kind=8), intent(in)                       :: t, dt
            real(kind=8), dimension(:), intent(in)         :: y
            real(kind=8), dimension(size (y)), intent(in)  :: der
            procedure(dydt_tem)                            :: dydt
            real(kind=8), dimension(size (y)), intent(out) :: ynew
            real(kind=8), dimension(size (y))              :: k2, k3, k4, k5, k6, k7
              
            ! k1 = dydt (t,             y)
            k2 = dydt (t + dt * C1_3, y + dt * der * C1_3)
            k3 = dydt (t + dt * C2_3, y + dt * k2 * C2_3)
            k4 = dydt (t + dt * C1_3, y + dt * (der + k2 * FOUR - k3) * C1_12)
            k5 = dydt (t + dt * C5_6, y + dt * (der * 25.d0 - k2 * 1.1d2 + k3 * 35.d0 + k4 * 90.d0)/48.d0)
            k6 = dydt (t + dt * C1_6, y + dt * (der * 0.15d0 - k2 * 0.55d0 - k3 * C1_8 + k4 * C1_2 + k5 * 0.1d0))
            k7 = dydt (t + dt,        y + dt * (- der * 195.75d0 + k2 * 495.d0 + k3 * 53.75d0 - k4 * 590.d0 + k5 * 32.d0 + &
            & k6 * 4.d2)/195.d0)
        
            ynew = y + dt * ((der + k7) * 13.d0 + (k3 + k4) * 55.d0 + (k5 + k6) * 32.d0) * 5.d-3
        end subroutine Abbas6

        !!

        !! Embedded

        subroutine Heun_Euler2_1 (t, y, der, dt, dydt, osol, oaux, yaux, ynew)
            implicit none
            real(kind=8), intent(in)                       :: t, dt
            real(kind=8), dimension(:), intent(in)         :: y
            real(kind=8), dimension(size (y)), intent(in)  :: der
            procedure(dydt_tem)                            :: dydt
            integer(kind=4), intent(out)                   :: osol, oaux
            real(kind=8), dimension(size (y)), intent(out) :: ynew, yaux
            real(kind=8), dimension(size (y))              :: k2
            
            osol = 2
            oaux = 1
            
            ! k1 = dydt (t,      y)
            k2 = dydt (t + dt, y + dt * der)
            
            ynew = y + dt * (der + k2) * C1_2
            yaux = y + dt * der
        end subroutine Heun_Euler2_1

        subroutine Fehlberg1_2 (t, y, der, dt, dydt, osol, oaux, yaux, ynew)
            implicit none
            real(kind=8), intent(in)                       :: t, dt
            real(kind=8), dimension(:), intent(in)         :: y
            real(kind=8), dimension(size (y)), intent(in)  :: der
            procedure(dydt_tem)                            :: dydt
            integer(kind=4), intent(out)                   :: osol, oaux
            real(kind=8), dimension(size (y)), intent(out) :: ynew, yaux
            real(kind=8), dimension(size (y))              :: k2, k3
            
            osol = 1
            oaux = 2
            
            ! k1 = dydt (t,             y)
            k2 = dydt (t + dt * C1_2, y + dt * der * C1_2)

            ynew = y + dt * (der + k2 * 255.d0) * 0.00390625d0

            k3 = dydt (t + dt, ynew)
            
            yaux = y + dt * (der + k2 * 500.d0 + k3) * 0.001953125d0
        end subroutine Fehlberg1_2

        subroutine Bogacki_Shampine3_2 (t, y, der, dt, dydt, osol, oaux, yaux, ynew)
            implicit none
            real(kind=8), intent(in)                       :: t, dt
            real(kind=8), dimension(:), intent(in)         :: y
            real(kind=8), dimension(size (y)), intent(in)  :: der
            procedure(dydt_tem)                            :: dydt
            integer(kind=4), intent(out)                   :: osol, oaux
            real(kind=8), dimension(size (y)), intent(out) :: ynew, yaux
            real(kind=8), dimension(size (y))              :: k2, k3, k4
            
            osol = 3
            oaux = 2
            
            ! k1 = dydt (t,             y)
            k2 = dydt (t + dt * C1_2, y + dt * der * C1_2)
            k3 = dydt (t + dt * C3_4, y + dt * k2 * C3_4)

            ynew = y + dt * (der * TWO + k2 * 3 + k3 * FOUR) * C1_9

            k4 = dydt (t + dt, ynew)           
            
            yaux = y + dt * (der * 7.d0/24.d0 + k2 * C1_4 + k3 * C1_3 + k4 * C1_8)
        end subroutine Bogacki_Shampine3_2

        subroutine Zonneveld4_3 (t, y, der, dt, dydt, osol, oaux, yaux, ynew)
            implicit none
            real(kind=8), intent(in)                       :: t, dt
            real(kind=8), dimension(:), intent(in)         :: y
            real(kind=8), dimension(size (y)), intent(in)  :: der
            procedure(dydt_tem)                            :: dydt
            integer(kind=4), intent(out)                   :: osol, oaux
            real(kind=8), dimension(size (y)), intent(out) :: ynew, yaux
            real(kind=8), dimension(size (y))              :: k2, k3, k4, k5, kaux
            real(kind=8)                                   :: aut
              
            aut = dt * C1_2
            
            osol = 4
            oaux = 3
            
            ! k1 = dydt (t,             y)
            k2 = dydt (t + aut,       y + dt * der * C1_2)
            k3 = dydt (t + aut,       y + dt * k2 * C1_2)
            k4 = dydt (t + dt,        y + dt * k3)
            k5 = dydt (t + dt * C3_4, y + dt * (der * FIVE + k2 * 7.d0 + k3 * 13.d0 - k4) * 0.03125d0)
            
            kaux = (k2 + k3)

            ynew = y + dt * (der + kaux * TWO + k4) * C1_6
            yaux = y + dt * (- der * THREE + kaux * 14.d0 + k4 * 13.d0 - k5 * 32.d0) * C1_6
        end subroutine Zonneveld4_3

        subroutine Merson4_3 (t, y, der, dt, dydt, osol, oaux, yaux, ynew)
            implicit none
            real(kind=8), intent(in)                       :: t, dt
            real(kind=8), dimension(:), intent(in)         :: y
            real(kind=8), dimension(size (y)), intent(in)  :: der
            procedure(dydt_tem)                            :: dydt
            integer(kind=4), intent(out)                   :: osol, oaux
            real(kind=8), dimension(size (y)), intent(out) :: ynew, yaux
            real(kind=8), dimension(size (y))              :: k2, k3, k4, k5
            
            osol = 4
            oaux = 3

            ! k1 = dydt (t,             y)
            k2 = dydt (t + dt * C1_3, y + dt * der * C1_3)
            k3 = dydt (t + dt * C1_3, y + dt * (der + k2) * C1_6)
            k4 = dydt (t + dt * C1_2, y + dt * (der + k3 * THREE) * C1_8)
            k5 = dydt (t + dt,        y + dt * (der - k3 * THREE + k4 * FOUR) * C1_2)
            
            ynew = y + dt * (der + k4 * FOUR + k5) * C1_6
            yaux = y + dt * (der + k3 * THREE + k4 * FOUR + k5 * TWO) * 0.1d0
        end subroutine Merson4_3

        subroutine Fehlberg4_5 (t, y, der, dt, dydt, osol, oaux, yaux, ynew)
            implicit none
            real(kind=8), intent(in)                       :: t, dt
            real(kind=8), dimension(:), intent(in)         :: y
            real(kind=8), dimension(size (y)), intent(in)  :: der
            procedure(dydt_tem)                            :: dydt
            integer(kind=4), intent(out)                   :: osol, oaux
            real(kind=8), dimension(size (y)), intent(out) :: ynew, yaux
            real(kind=8), dimension(size (y))              :: k2, k3, k4, k5, k6
            
            osol = 4
            oaux = 5
            
            ! k1 = dydt (t,                    y)
            k2 = dydt (t + dt * C1_4,        y + dt * der * C1_4)
            k3 = dydt (t + dt * C3_8,        y + dt * (der * THREE + k2 * 9.d0) * 0.03125d0)
            k4 = dydt (t + dt * 12.d0/13.d0, y + dt * (der * 1932.d0 - k2 * 7200.d0 + k3 * 7296.d0)/2197.d0)
            k5 = dydt (t + dt,               y + dt * ((der * 8341.d0 + k3 * 29440.d0 - k4 * 845.d0)/4104.d0 - k2 * 8.d0))
            k6 = dydt (t + dt * C1_2,        y + dt * ((- der * 1216.d0 + k4 * 1859.d0)/4104.d0 + k2 * TWO - &
                                               & k3 * 3544.d0/2565.d0 - k5 * 0.275d0))
            
            ynew = y + dt * ((der * 475.d0 + k4 * 2197.d0)/4104.d0 + k3 * 1408.d0/2565.d0 - k5 * C1_5)
            yaux = y + dt * ((der * 6688.d0 + k4 * 28561.d0 + k6 * 2052.d0)/56430.d0 + k3 * 6656.d0/12825.d0 - k5 * 0.18d0)
        end subroutine Fehlberg4_5

        subroutine Cash_Karp5_4 (t, y, der, dt, dydt, osol, oaux, yaux, ynew)
            implicit none
            real(kind=8), intent(in)                       :: t, dt
            real(kind=8), dimension(:), intent(in)         :: y
            real(kind=8), dimension(size (y)), intent(in)  :: der
            procedure(dydt_tem)                            :: dydt
            integer(kind=4), intent(out)                   :: osol, oaux
            real(kind=8), dimension(size (y)), intent(out) :: ynew, yaux
            real(kind=8), dimension(size (y))              :: k2, k3, k4, k5, k6
            
            osol = 5
            oaux = 4
            
            ! k1 = dydt (t,                y)
            k2 = dydt (t + dt * C1_5,    y + dt * der * C1_5)
            k3 = dydt (t + dt * 0.3d0,   y + dt * (der * THREE + k2 * 9.d0) * 0.025d0)
            k4 = dydt (t + dt * 0.6d0,   y + dt * (der * THREE - k2 * 9.d0 + k3 * 12.d0) * 0.1d0)
            k5 = dydt (t + dt,           y + dt * (- der * 11.d0 + k2 * 135.d0 - k3 * 140.d0 + k4 * 70.d0)/54.d0)
            k6 = dydt (t + dt * 0.875d0, y + dt * (der * 3262.d0 + k2 * 37800.d0 + k3 * 4600.d0 + k4 * 44275.d0 + &
                                           & k5 * 6831.d0)/110592.d0)
            
            ynew = y + dt * (der * 37.d0/378.d0 + k3 * 250.d0/621.d0 + k4 * 125.d0/594.d0 + k6 * 512.d0/1771.d0)
            yaux = y + dt * ((der * 5650.d0 + k4 * 13525.d0)/55296.d0 + k3 * 18575.d0/48384.d0 + k5 * 277.d0/14336.d0 + k6 * C1_4)
        end subroutine Cash_Karp5_4

        subroutine Dormand_Prince5_4 (t, y, der, dt, dydt, osol, oaux, yaux, ynew)
            implicit none
            real(kind=8), intent(in)                       :: t, dt
            real(kind=8), dimension(:), intent(in)         :: y
            real(kind=8), dimension(size (y)), intent(in)  :: der
            procedure(dydt_tem)                            :: dydt
            integer(kind=4), intent(out)                   :: osol, oaux
            real(kind=8), dimension(size (y)), intent(out) :: ynew, yaux
            real(kind=8), dimension(size (y))              :: k2, k3, k4, k5, k6, k7

            osol = 5
            oaux = 4

            ! k1 = dydt (t,              y)
            k2 = dydt (t + dt * C1_5,  y + dt * der * C1_5)
            k3 = dydt (t + dt * 0.3d0, y + dt * (der * 0.075d0 + k2 * 0.225d0))
            k4 = dydt (t + dt * 0.8d0, y + dt * (der * 44.d0 - k2 * 168.d0 + k3 * 160.d0)/45.d0)
            k5 = dydt (t + dt * C8_9,  y + dt * (der * 19372.d0 - k2 * 76080.d0 + k3 * 64448.d0 - k4 * 1908.d0)/6561.d0)
            k6 = dydt (t + dt,         y + dt * ((der * 9017.d0 - k2 * 34080.d0 + k4 * 882.d0)/3168.d0 + &
                                         & k3 * 46732.d0/5247.d0 - k5 * 5103.d0/18656.d0))
            
            ynew = y + dt * ((der * 35.d0 + k4 * 250.d0)/384.d0 + k3 * 500.d0/1113.d0 - k5 * 2187.d0/6784.d0 + k6 * 11.d0/84.d0)
            
            k7 = dydt (t + dt, ynew)

            yaux = y + dt * ((der * 5179.d0 + k4 * 35370.d0)/57600.d0 + k3 * 7571.d0/16695.d0 - k5 * 92097.d0/339200.d0 + &
                     & k6 * 187.d0/2100.d0 + k7 * 0.025d0)
        end subroutine Dormand_Prince5_4

        subroutine Verner6_5 (t, y, der, dt, dydt, osol, oaux, yaux, ynew)
            implicit none
            real(kind=8), intent(in)                       :: t, dt
            real(kind=8), dimension(:), intent(in)         :: y
            real(kind=8), dimension(size (y)), intent(in)  :: der
            procedure(dydt_tem)                            :: dydt
            integer(kind=4), intent(out)                   :: osol, oaux
            real(kind=8), dimension(size (y)), intent(out) :: ynew, yaux
            real(kind=8), dimension(size (y))              :: k2, k3, k4, k5, k6, k7, k8

            osol = 6
            oaux = 5
            
            ! k1 = dydt (t,              y)
            k2 = dydt (t + dt * C1_6,  y + dt * der * C1_6)
            k3 = dydt (t + dt * C4_15, y + dt * (der * FOUR + k2 * 16.d0)/75.d0)
            k4 = dydt (t + dt * C2_3,  y + dt * (der * FIVE - k2 * 16.d0 + k3 * 15.d0) * C1_6)
            k5 = dydt (t + dt * C5_6,  y + dt * ((- der * 165.d0 - k3 * 425.d0)/64.d0 + (k2 * 880.d0 + k4 * 85.d0)/96.d0))
            k6 = dydt (t + dt,         y + dt * ((der * 612.d0 + k5 * 88.d0)/255.d0 - k2 * 8.d0 + (k3 * 4015.d0 - &
                                         & k4 * 187.d0)/612.d0))
            k7 = dydt (t + dt * C1_15, y + dt * ((- der * 8263.d0 + k2 * 24800.d0)/15000.d0 - k3 * 643.d0/680.d0 - k4 * 0.324d0 + &
                                         & k5 * 2484.d0/10625.d0))
            k8 = dydt (t + dt,         y + dt * (der * 3501.d0/1720.d0 + (297275.d0 * k3 - 367200 * k2)/52632.d0 - &
                                         & k4 * 319.d0/2322.d0 + k5 * 24068.d0/84065.d0 + k7 * 3850.d0/26703.d0))
            
            ynew = y + dt * (der * 0.075d0 + k3 * 875.d0/2244.d0 + (k4 * 3703.d0 + k7 * 125.d0)/11592.d0 + k5 * 264.d0/1955.d0 + &
                     & k8 * 43.d0/616.d0)
            yaux = y + dt * ((der * 13.d0 + k4 * 50.d0)/160.d0 + (k3 * 2375.d0 + k6 * 408.d0)/5984.d0 + k5 * 12.d0/85.d0)
        end subroutine Verner6_5

        subroutine Fehlberg7_8 (t, y, der, dt, dydt, osol, oaux, yaux, ynew)
            implicit none
            real(kind=8), intent(in)                       :: t, dt
            real(kind=8), dimension(:), intent(in)         :: y
            real(kind=8), dimension(size (y)), intent(in)  :: der
            procedure(dydt_tem)                            :: dydt
            integer(kind=4), intent(out)                   :: osol, oaux
            real(kind=8), dimension(size (y)), intent(out) :: ynew, yaux
            real(kind=8), dimension(size (y))              :: k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12, k13, kaux
            
            osol = 7
            oaux = 8
            
            ! k1   = dydt (t,              y)
            k2   = dydt (t + dt * C2_27, y + dt * der * C2_27)
            k3   = dydt (t + dt * C1_9,  y + dt * (der + k2 * THREE)/36.d0)
            k4   = dydt (t + dt * C1_6,  y + dt * (der + k3 * THREE)/24.d0)
            k5   = dydt (t + dt * C5_12, y + dt * (der * C5_12 + (k4 - k3) * 1.5625d0))
            k6   = dydt (t + dt * C1_2,  y + dt * (der + k4 * FIVE + k5 * FOUR) * 0.02d0)
            k7   = dydt (t + dt * C5_6,  y + dt * (- der * 25.d0 + k4 * 125.d0 - k5 * 260.d0 + k6 * 250.d0)/108.d0)
            k8   = dydt (t + dt * C1_6,  y + dt * (der * 93.d0 + k5 * 244.d0 - k6 * 200.d0 + k7 * 13.d0)/900.d0)
            k9   = dydt (t + dt * C2_3,  y + dt * (der * TWO + (- k4 * 795.d0 + k5 * 1408.d0 - k6 * 1070.d0 +&
                                           & k7 * 67.d0)/9.d1 + k8 * THREE))
            k10  = dydt (t + dt * C1_3,  y + dt * ((- der * 91.d0 + k4 * 23.d0  + k6 * 622.d0)/108.d0 - k5 * 976.d0/135.d0 + &
                                           & (- k7 * 19.d0 + k8 * 170.d0 - k9 * 5.d0)/60.d0))
            kaux = - k4 * 8525 + k5 * 17984
            k11  = dydt (t + dt,         y + dt * (der * 2383.d0 + kaux - k6 * 15050.d0 + k7 * 2133.d0 + &
                                           & k8 * 2250.d0 + k9 * 1125.d0 + k10 * 1800.d0)/4100.d0)
            k12  = dydt (t,              y + dt * (((der - k7)/205.d0 + ((- k6 + k10) * TWO + (- k8 + k9))/41.d0) * THREE))
            k13  = dydt (t + dt,         y + dt * ((- der * 1777.d0 + kaux - k6 * 14450.d0 + k7 * 2193.d0 + &
                                           & k8 * 2550.d0 + k9 * 825.d0 + k10 * 1900.d0)/4100.d0 + k12))

            kaux = k6 * 272.d0 + (k7 + k8) * 216.d0 + (k9 + k10) * 27.d0

            ynew = y + dt * ((der + k11) * 41.d0 + kaux) * C1_840
            yaux = y + dt * (kaux + (k12 + k13) * 41.d0) * C1_840
        end subroutine Fehlberg7_8
        
        subroutine Dormand_Prince8_7 (t, y, der, dt, dydt, osol, oaux, yaux, ynew)
            implicit none
            real(kind=8), intent(in)                       :: t, dt
            real(kind=8), dimension(:), intent(in)         :: y
            real(kind=8), dimension(size (y)), intent(in)  :: der
            procedure(dydt_tem)                            :: dydt
            integer(kind=4), intent(out)                   :: osol, oaux
            real(kind=8), dimension(size (y)), intent(out) :: ynew, yaux
            real(kind=8), dimension(size (y))              :: k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12, k13

            osol = 8
            oaux = 7

            ! k1  = dydt (t, y)
            k2  = dydt (t + dt * C1_18, y + dt * (C1_18 * der))
            k3  = dydt (t + dt * C1_12, y + dt * (C1_48 * der + 0.0625d0 * k2))
            k4  = dydt (t + dt * C1_8, y + dt * (0.03125d0 * der + 0.09375d0 * k3))
            k5  = dydt (t + dt * 0.3125d0, y + dt * (0.3125d0 * der - 1.171875d0 * k3 +  1.171875d0 * k4))
            k6  = dydt (t + dt * C3_8, y + dt * (0.0375d0 * der + 0.1875d0 * k4 + 0.15d0 * k5))
            k7  = dydt (t + dt * 0.1475d0, y + dt * (29443841.d0/614563906.d0 * der + &
            & 77736538.d0/692538347.d0 * k4 - 28693883.d0/1125000000.d0 * k5 + 23124283.d0/1800000000.d0 * k6))
            k8  = dydt (t + dt * 0.465d0, y + dt * (16016141.d0/946692911.d0 * der + &
            & 61564180.d0/158732637.d0 * k4 + 22789713.d0/633445777.d0 * k5 + 545815736.d0/2771057229.d0 * k6 - &
            &   180193667.d0/1043307555.d0 * k7))
            k9  = dydt (t + dt * 5490023248.d0/9719169821.d0, y + dt * (39632708.d0/573591083.d0 * der - &
            & 433636366.d0/683701615.d0 * k4 - 421739975.d0/2616292301.d0 * k5 + 100302831.d0/723423059.d0 * k6 + &
            & 790204164.d0/839813087.d0 * k7 + 800635310.d0/3783071287.d0 * k8))
            k10 = dydt (t + dt * 0.65d0, y + dt * (246121993.d0/1340847787.d0 * der - &
            & 37695042795.d0/15268766246.d0 * k4 - 309121744.d0/1061227803.d0 * k5 - 12992083.d0/490766935.d0 * k6 + &
            &  6005943493.d0/2108947869.d0 * k7 + 393006217.d0/1396673457.d0 * k8 + 123872331.d0/1001029789.d0 * k9))
            k11 = dydt (t + dt * 1201146811.d0/1299019798.d0, y + dt * (-1028468189.d0/846180014.d0 * der + &
            &    8478235783.d0/508512852.d0 * k4 + 1311729495.d0/1432422823.d0 * k5 - 10304129995.d0/1701304382.d0 * k6 - &
            & 48777925059.d0/3047939560.d0 * k7 + 15336726248.d0/1032824649.d0 * k8 - 45442868181.d0/3398467696.d0 * k9 + &
            &  3065993473.d0/597172653.d0 * k10))
            k12 = dydt (t + dt, y + dt * (185892177.d0/718116043.d0 * der - &
            & 3185094517.d0/667107341.d0 * k4 - 477755414.d0/1098053517.d0 * k5 - 703635378.d0/230739211.d0 * k6 + &
            & 5731566787.d0/1027545527.d0 * k7 + 5232866602.d0/850066563.d0 * k8 - 4093664535.d0/808688257.d0 * k9 + &
            & 3962137247.d0/1805957418.d0 * k10 + 65686358.d0/487910083.d0 * k11))
            k13 = dydt (t + dt, y + dt * (403863854.d0/491063109.d0 * der - &
            & 5068492393.d0/434740067.d0 * k4 - 411421997.d0/543043805.d0 * k5 + 652783627.d0/914296604.d0 * k6 + &
            & 11173962825.d0/925320556.d0 * k7 - 13158990841.d0/6184727034.d0 * k8 +  3936647629.d0/1978049680.d0 * k9 - &
            & 160528059.d0/685178525.d0 * k10 + 248638103.d0/1413531060.d0 * k11))

            ynew  = y + dt * (13451932.d0/455176623.d0 * der - 808719846.d0/976000145.d0 * k6 + &
            & 1757004468.d0/5645159321.d0 * k7 + 656045339.d0/265891186.d0 * k8 - 3867574721.d0/1518517206.d0 * k9 + &
            & 465885868.d0/322736535.d0 * k10 + 53011238.d0/667516719.d0 * k11 + C2_45 * k12)
            yaux  = y + dt * (14005451.d0/335480064.d0 * der - 59238493.d0/1068277825.d0 * k6 + &
            & 181606767.d0/758867731.d0 * k7 + 561292985.d0/797845732.d0 * k8 - 1041891430.d0/1371343529.d0 * k9 + &
            & 760417239.d0/1151165299.d0 * k10 +  118820643.d0/751138087.d0 * k11 - 528747749.d0/2220607170.d0 * k12 + C1_4 * k13)
        end subroutine Dormand_Prince8_7

        !!
        
        !! Bulirsch Stoer main integrator
    
        subroutine Bulirsch_Stoer (t, y, dt_adap, dydt, e_tol, dt_used, ynew)
            implicit none
            real(kind=8), intent(in)                       :: t, e_tol
            real(kind=8), intent(inout)                    :: dt_adap, dt_used
            real(kind=8), dimension(:), intent(in)         :: y
            procedure(dydt_tem)                            :: dydt
            real(kind=8), dimension(size (y)), intent(out) :: ynew
            real(kind=8)                                   :: time, dtry
            integer(kind=4), save                          :: sizey
            
            sizey = size (y)
            time  = t
            dtry  = dt_adap
            call bstep (ynew, dydt, sizey, time, dtry, e_tol, dt_used, dt_adap)
        end subroutine Bulirsch_Stoer

        !!!! Auxiliar subroutines for Bulirsch_Stoer

        subroutine bstep (y, dydt, sizey, x, htry, eps, hdid, hnext)
            implicit none
            integer(kind=4), intent(in) :: sizey
            real(kind=8), dimension(sizey), intent(inout) :: y
            real(kind=8), intent(inout) :: x
            real(kind=8), intent(in)    :: htry, eps
            real(kind=8), intent(out)   :: hdid, hnext
            procedure(dydt_tem)         :: dydt
            integer(kind=4), parameter :: imax = 9, kmaxx = imax - 1
            real(kind=8), parameter    :: safe1 = .25d0, safe2 = .7d0
            real(kind=8), parameter    :: redmax = 1.d-5, redmin = .7d0
            real(kind=8), parameter    :: tini = 1.d-30, scalmx = .1d0
            integer(kind=4), parameter, dimension(imax) :: nseq = (/2, 4, 6, 8, 10, 12, 14, 16, 18/)
            integer(kind=4), save                       :: kmax, kopt
            real(kind=8), dimension(kmaxx, kmaxx), save :: alf
            real(kind=8), dimension(kmaxx)              :: err
            real(kind=8), dimension(imax), save         :: arr ! a in NR F90
            real(kind=8), save                          :: xnew, epsold = -1.d0
            real(kind=8)                   :: eps1, errmax, fact, h, red, scala, wrkmin, xest
            real(kind=8), dimension(sizey) :: yerr, ysav, yseq, der, yscal
            logical       :: reduct
            logical, save :: first = .True.

            real(kind=8), dimension(kmaxx * 2)        :: xpz
            real(kind=8), dimension(sizey, kmaxx * 2) :: qcolpz
            real(kind=8)    :: work
            integer(kind=4) :: k, iq, i ,km, kk

            der = dydt (x, y)
            yscal = abs (y) + abs (htry * der) + SAFE_LOW
            
            if (abs(eps - epsold) > SAFE_LOW) then !E_TOL? ! A new tolerance, so reinitialize.
                hnext = -1.0d29 ! “Impossible” values.
                xnew  = -1.0d29
                eps1  = safe1 * eps
                arr(1) = nseq(1) + 1
                do k = 1, kmaxx
                    arr(k + 1) = arr(k) + nseq(k + 1)
                end do
                ! Compute α(k, q):
                do iq = 2, kmaxx
                    do k = 1, iq - 1
                        alf(k, iq) = eps1**((arr(k + 1) - arr(iq + 1)) / ((arr(iq + 1) - arr(1) + ONE) * (TWO * k + 1)))
                    end do
                end do
                epsold = eps
                do kopt = 2, kmaxx - 1 ! Determine optimal row number for convergence.
                    if (arr(kopt + 1) > arr(kopt) * alf(kopt - 1, kopt)) exit
                end do
                kmax = kopt
            end if
            h = htry
            ysav = y
            if ((abs(h - hnext) > SAFE_LOW) .or. (abs(x - xnew) > SAFE_LOW)) then !E_TOL? ! A new stepsize or a new integration: Reestablish the order window.
                first = .True.
                kopt = kmax
            end if
            reduct = .False.
            main_loop: do
                do k = 1, kmax ! Evaluate the sequence of modiﬁed midpoint integrations.
                    xnew = x + h
                    if (abs(xnew - x) < SAFE_LOW) then !E_TOL?
                        print*, "Step size underflow in bstep at ", x
                        exit main_loop ! Luckily, hard_exit will handle it
                        ! stop 2
                        ! return
                    end if
                    call mmid (ysav, der, sizey, x, h, nseq(k), yseq, dydt)
                    yscal = abs (y) + abs (h * der) + SAFE_LOW
                    xest = (h / nseq(k))**2 ! Squared, since error series is even.
                    call pzextr (k, xest, yseq, y, yerr, sizey, qcolpz, xpz) ! Perform extrapolation.
                    if (k .ne. 1) then ! Compute normalized error estimate eps(k).
                        errmax = tini
                        do i = 1, sizey
                            errmax = max (errmax, abs(yerr(i)/yscal(i))) ! FALTA SCALADO? !! Ya no?
                        end do
                        errmax = errmax / eps
                        km = k - 1
                        err(km) = (errmax / safe1)**(ONE / (TWO * km + 1))
                    end if
                    if (k .ne. 1 .and. (k .ge. kopt - 1 .or. first)) then ! In order window.
                        if (errmax < ONE) exit main_loop ! Converged.
                        if (k .eq. kmax .or. k .eq. kopt + 1) then ! Check for possible stepsize reduction.
                            red = safe2 / err(km)
                            exit
                        else if (k .eq. kopt) then
                            if (alf(kopt - 1, kopt) < err(km)) then
                                red = ONE / err(km)
                                exit
                            end if
                        else if (kopt .eq. kmax) then
                            if (alf(km, kmax - 1) < err(km)) then
                                red = alf(km, kmax - 1) * safe2 / err(km)
                                exit
                            end if
                        else if (alf(km, kopt) < err(km)) then
                            red = alf(km, kopt - 1) / err(km)
                            exit
                        end if
                    end if
                end do
                red = max(min(red, redmin), redmax)
                h = h * red
                reduct = .True.
            end do main_loop
            x = xnew
            hdid = h
            first = .False.
            wrkmin = 1.d35
            do kk = 1, km
                fact = max (err(kk), scalmx)
                work = fact * arr(kk + 1)
                if (work < wrkmin) then
                    scala = fact
                    wrkmin = work
                    kopt = kk + 1
                end if
            end do
            hnext = h / scala
            if (kopt .ge. k .and. kopt .ne. kmax .and. .not. reduct) then! Check for possible order increase, but not if stepsize was just reduced.
                fact = max (scala / alf(kopt - 1, kopt), scalmx)
                if (arr(kopt + 1) * fact .le. wrkmin) then
                    hnext = h / fact
                    kopt = kopt + 1
                end if
            end if
        end subroutine bstep

        subroutine mmid (y, dydx, sizey, xs, htot, nstep, yout, dydt)
            integer(kind=4), intent(in) :: sizey, nstep
            procedure(dydt_tem) :: dydt
            real(kind=8), intent(in) :: xs, htot
            real(kind=8), dimension(sizey) :: ym, yn
            real(kind=8), dimension(sizey), intent(in) :: y, dydx
            real(kind=8), dimension(sizey), intent(out) :: yout
            real(kind=8) :: x, h, h2, swap
            integer(kind=4) :: i, n 

            h  = htot / (nstep * ONE) ! Stepsize this trip.
            ym = y
            yn = y + h * dydx ! First step.
            x  = xs + h
            yout = dydt(x, yn) ! Will use yout for temporary storage of derivatives.
            h2 = TWO * h
            do n = 2, nstep ! General step.
                do i = 1, sizey
                    swap = ym(i) + h2 * yout(i)
                    ym(i) = yn(i)
                    yn(i) = swap
                end do
                x = x + h
                yout = dydt(x, yn)
            end do
            yout = C1_2 * (ym + yn + h * yout) ! Last step.
        end subroutine mmid

        subroutine pzextr (iest, xest, yest, yz, dy, sizey, qcol, x)
            integer(kind=4), intent(in) :: iest, sizey
            real(kind=8), intent(in)    :: xest
            real(kind=8), dimension(sizey)              :: d
            real(kind=8), dimension(sizey), intent(in)  :: yest
            real(kind=8), dimension(sizey), intent(out) :: dy, yz
            integer(kind=4), parameter                  :: IEST_MAX=16
            real(kind=8), dimension(IEST_MAX), intent(inout)        :: x
            real(kind=8), dimension(sizey, IEST_MAX), intent(inout) :: qcol
            integer(kind=4) :: j, k1
            real(kind=8)    :: delta, f1, f2, q

            x(iest) = xest  ! Save current independent variable.
            dy = yest
            yz = yest
            if (iest .eq. 1) then  ! Store ﬁrst estimate in ﬁrst column.
                qcol(:, 1) = yest(:)
            else
                d = yest
                do k1 = 1, iest - 1
                    delta = ONE / (x(iest - k1) - xest)
                    f1 = xest * delta
                    f2 = x(iest - k1) * delta
                    do j = 1, sizey ! Propagate tableau 1 diagonal more.
                        q = qcol(j, k1)
                        qcol(j, k1) = dy(j)
                        delta = d(j) - q
                        dy(j) = f1 * delta
                        d(j)  = f2 * delta
                        yz(j) = yz(j) + dy(j)
                    end do  
                end do
                qcol(:, iest) = dy(:)
            end if
        end subroutine pzextr

        !!!! End auxiliar subroutines for Bulirsch_Stoer

        !!
        
        !! Bulirsch Stoer 2 integrator (only useful for [x, y, vx, vy,...] systems)

        subroutine Bulirsch_Stoer2 (t, y, dt_adap, dydt, e_tol, dt_used, ynew)
            implicit none
            real(kind=8), intent(in) :: t, e_tol
            real(kind=8), intent(inout) :: dt_adap, dt_used
            real(kind=8), dimension(:), intent(in) :: y
            procedure(dydt_tem) :: dydt
            real(kind=8), dimension(size (y)), intent(out) :: ynew
            real(kind=8) :: time, dtry
            integer(kind=4), save :: sizex
            
            sizex = int(size(y) / 4, 4)
            time  = t
            dtry  = dt_adap
            call bstep2 (ynew, dydt, sizex, t, dtry, e_tol, dt_used, dt_adap)
        end subroutine Bulirsch_Stoer2

        !!!! Auxiliar subroutine for Bulirsch_Stoer2

        subroutine bstep2 (y, dydt, sizex, t, htry, eps, hdid, hnext)
            implicit none
            procedure(dydt_tem) :: dydt
            integer(kind=4), intent(in) :: sizex
            real(kind=8), intent(in) :: htry, t, eps
            real(kind=8), intent(inout) :: y(:)
            real(kind=8), intent(out) :: hdid, hnext
            real(kind=8), parameter :: SHRINK = 0.55d0
            real(kind=8), parameter :: GROW   = 1.3d0
            real(kind=8) :: der(sizex*4), dt, tt
            integer(kind=4) :: i, j, j1, k, ns, ndo
            real(kind=8) :: tmp0, tmp1, tmp2, errmax, h, hx2, h2(8)
            real(kind=8), dimension(2,sizex) :: x, v, x0, v0, a0, xend, vend, a
            real(kind=8) :: xscal(sizex), vscal(sizex)
            real(kind=8) :: dx(2,sizex,8), dv(2,sizex,8)

            !------------------------------------------------------------------------------
            ! Calculate accelerations at the start of the step
            der = dydt(t, y)

            do i = 1, sizex
                x0(1,i) = y(4*i - 3)
                x0(2,i) = y(4*i - 2)
                v0(1,i) = y(4*i - 1)
                v0(2,i) = y(4*i)
                a0(1,i) = der(4*i - 1)
                a0(2,i) = der(4*i)
                ! Arrays used to scale the relative error (1 / absolute distance or velocity)
                xscal(i) = ONE / sqrt(x0(1,i)*x0(1,i) + x0(2,i)*x0(2,i))
                vscal(i) = ONE / sqrt(v0(1,i)*v0(1,i) + v0(2,i)*v0(2,i))
            end do
            
            dt = htry
            ! Try to do the step using current stepsize
            do ndo = 0, MAX_N_ITER
                tt = t
                ! For each value of NS, do a modified-midpoint integration with 2N substeps
                do ns = 1, 8
                    h = dt / (TWO * dble(ns))
                    h2(ns) = C1_4 / (ns * ns)
                    hx2 = h * TWO
                    x(:,1:sizex) = x0(:,1:sizex) + h * v0(:,1:sizex)
                    v(:,1:sizex) = v0(:,1:sizex) + h * a0(:,1:sizex)
                    tt = tt + h

                    a = get_a(tt, dydt, x, v)
                    xend(:,1:sizex) = x0(:,1:sizex) + hx2 * v(:,1:sizex)
                    vend(:,1:sizex) = v0(:,1:sizex) + hx2 * a(:,1:sizex)
                    tt = tt + hx2
                    
                    do j = 2, ns
                        a  = get_a(tt, dydt, xend, vend)
                        x(:,1:sizex) = x(:,1:sizex) + hx2 * vend(:,1:sizex)
                        v(:,1:sizex) = v(:,1:sizex) + hx2 * a(:,1:sizex)
                        tt = tt + hx2
                        
                        a  = get_a(tt, dydt, x, v)
                        xend(:,1:sizex) = xend(:,1:sizex) + hx2 * v(:,1:sizex)
                        vend(:,1:sizex) = vend(:,1:sizex) + hx2 * a(:,1:sizex)
                        tt = tt + hx2
                    end do
                    
                    a = get_a(tt, dydt, xend, vend)
                    dx(:,1:sizex,ns) = C1_2 * (xend(:,1:sizex) + x(:,1:sizex) + h * vend(:,1:sizex))
                    dv(:,1:sizex,ns) = C1_2 * (vend(:,1:sizex) + v(:,1:sizex) + h * a(:,1:sizex))
            
                    ! Update the DX and DV arrays used for polynomial extrapolation
                    do j = ns - 1, 1, -1
                        j1 = j + 1
                        tmp0 = ONE / (h2(j) - h2(ns))
                        tmp1 = tmp0 * h2(j1)
                        tmp2 = tmp0 * h2(ns)
                        dx(:,1:sizex,j) = tmp1 * dx(:,1:sizex,j1) - tmp2 * dx(:,1:sizex,j)
                        dv(:,1:sizex,j) = tmp1 * dv(:,1:sizex,j1) - tmp2 * dv(:,1:sizex,j)
                    end do
            
                    ! After several integrations, test the relative error on extrapolated values
                    if (ns > 3) then
                        ! Maximum relative position and velocity errors (last D terms added)
                        errmax = ZERO
                        do k = 1, sizex
                            tmp1 = sqrt(dx(1,k,1)*dx(1,k,1) + dx(2,k,1)*dx(2,k,1))
                            tmp2 = sqrt(dv(1,k,1)*dv(1,k,1) + dv(2,k,1)*dv(2,k,1))
                            errmax = max(errmax, tmp1*xscal(k), tmp2*vscal(k))
                        end do
                        ! If error is smaller than TOL, update position and velocity arrays, and exit
                        if (errmax <= eps) then ! * eps?
                            x0 = ZERO
                            v0 = ZERO
                            do j = 1, ns
                                x0(:,1:sizex) = x0(:,1:sizex) + dx(:,1:sizex,j)
                                v0(:,1:sizex) = v0(:,1:sizex) + dv(:,1:sizex,j)
                            end do
                            ! Recommend a stepsize for the next call to this subroutine
                            if (ns >= 8) hnext = dt * SHRINK
                            if (ns < 7)  hnext = dt * GROW
                            ! Update
                            hdid = dt
                            do i = 1, sizex
                                y(4*i - 3) = x0(1,i)
                                y(4*i - 2) = x0(2,i)
                                y(4*i - 1) = v0(1,i)
                                y(4*i)     = v0(2,i)
                            end do
                            !end
                            return
                        end if
                    end if
                end do
                ! If errors were too large, redo the step with half the previous step size.
                dt = dt * C1_2
            end do

            ! Check max iter
            if (ndo == MAX_N_ITER) then
                print*,  "Warning, eps not reached in bstep2"
                hdid  = dt * TWO
                hnext = hdid
            end if

            ! Update
            hdid = dt
            do i = 1, sizex
                y(4*i - 3) = x0(1,i)
                y(4*i - 2) = x0(2,i)
                y(4*i - 1) = v0(1,i)
                y(4*i)     = v0(2,i)
            end do
            !end

            contains
                function get_a (t, dydt, x, v) result(a)
                    implicit none
                    real(kind=8), intent(in) :: t
                    procedure(dydt_tem) :: dydt
                    real(kind=8), dimension(2,sizex), intent(in) :: x, v
                    real(kind=8), dimension(2,sizex) :: a
                    real(kind=8), dimension(sizex*4) :: y, der
                    integer(kind=4) :: i
                    do i = 1, sizex
                        y(4*i - 3) = x(1,i)
                        y(4*i - 2) = x(2,i)
                        y(4*i - 1) = v(1,i)
                        y(4*i)     = v(2,i)
                    end do
                    der = dydt(t, y)
                    do i = 1, sizex
                        a(1,i) = der(4*i - 1)
                        a(2,i) = der(4*i)
                    end do
                end function get_a
        end subroutine bstep2

        !!

        !! Leap Frog Integrators (only useful for [x, y, vx, vy,...] systems)

        !!! Leap Frog (KDK)

        subroutine leapfrog_KDK (t, y, der, dt, dydt, ynew)
            implicit none
            real(kind=8), intent(in) :: t, dt
            real(kind=8), dimension(:), intent(in) :: y
            real(kind=8), dimension(size (y)), intent(in) :: der
            procedure(dydt_tem) :: dydt
            real(kind=8), dimension(size (y)), intent(out) :: ynew
            real(kind=8), dimension(2,size (y)) :: x0, v0, x1, v1, v05, a0
            integer(kind=4) :: i, sizey

            sizey = int(size(y) / 4, 4)

            !------------------------------------------------------------------------------
            ! Get x0, v0, and a0
            do i = 1, sizey
                x0(1,i) = y(4*i - 3)
                x0(2,i) = y(4*i - 2)
                v0(1,i) = y(4*i - 1)
                v0(2,i) = y(4*i)
                a0(1,i) = der(4*i - 1)
                a0(2,i) = der(4*i)
            end do
            ! Calculate v05
            v05 = v0 + a0 * dt * C1_2
            ! Calculate x1 using v05
            x1  = x0 + v05 * dt
            ! Calculate v1 using accelerations at x1 and v05
            v1  = v05 + get_a(t + dt, dydt, x1, v05, sizey) * dt * C1_2

            ! Update ynew
            do i = 1, sizey
                ynew(4*i - 3) = x1(1,i)
                ynew(4*i - 2) = x1(2,i)
                ynew(4*i - 1) = v1(1,i)
                ynew(4*i)     = v1(2,i)
            end do

            contains
                function get_a (t, dydt, x, v, sizey) result(a)
                    implicit none
                    real(kind=8), intent(in)                     :: t
                    procedure(dydt_tem)                          :: dydt
                    integer(kind=4), intent(in)                  :: sizey
                    real(kind=8), dimension(2,sizey), intent(in) :: x, v
                    real(kind=8), dimension(2,sizey)             :: a
                    real(kind=8), dimension(sizey*4) :: y, der
                    integer(kind=4)                  :: i
                    do i = 1, sizey
                        y(4*i - 3) = x(1,i)
                        y(4*i - 2) = x(2,i)
                        y(4*i - 1) = v(1,i)
                        y(4*i)     = v(2,i)
                    end do
                    der = dydt(t, y)
                    do i = 1, sizey
                        a(1,i) = der(4*i - 1)
                        a(2,i) = der(4*i)
                    end do
                end function get_a
        end subroutine leapfrog_KDK

        !!! Leap Frog (DKD)

        subroutine leapfrog_DKD (t, y, der, dt, dydt, ynew)
            implicit none
            real(kind=8), intent(in) :: t, dt
            real(kind=8), dimension(:), intent(in) :: y
            real(kind=8), dimension(size (y)), intent(in) :: der
            procedure(dydt_tem) :: dydt
            real(kind=8), dimension(size (y)), intent(out) :: ynew
            real(kind=8), dimension(2,size (y)) :: x0, v0, x1, v1, x05
            integer(kind=4) :: i, sizey

            sizey = int(size(y) / 4, 4)

            !------------------------------------------------------------------------------
            ! Get x0, v0
            do i = 1, sizey
                x0(1,i) = y(4*i - 3)
                x0(2,i) = y(4*i - 2)
                v0(1,i) = y(4*i - 1)
                v0(2,i) = y(4*i)
            end do

            ! Calculate x05 at the start of the step
            x05  = x0 + v0 * dt * C1_2
            ! Calculate v1 using accelerations at x05 and v0
            v1  = v0 + get_a(t + dt, dydt, x05, v0, sizey) * dt
            ! Calculate x1 using x05 and v1
            x1  = x05 + v1 * dt * C1_2

            ! Update ynew
            do i = 1, sizey
                ynew(4*i - 3) = x1(1,i)
                ynew(4*i - 2) = x1(2,i)
                ynew(4*i - 1) = v1(1,i)
                ynew(4*i)     = v1(2,i)
            end do

            contains
                function get_a (t, dydt, x, v, sizey) result(a)
                    implicit none
                    real(kind=8), intent(in) :: t
                    procedure(dydt_tem) :: dydt
                    integer(kind=4), intent(in) :: sizey
                    real(kind=8), dimension(2,sizey), intent(in) :: x, v
                    real(kind=8), dimension(2,sizey) :: a
                    real(kind=8), dimension(sizey*4) :: y, der
                    integer(kind=4) :: i
                    do i = 1, sizey
                        y(4*i - 3) = x(1,i)
                        y(4*i - 2) = x(2,i)
                        y(4*i - 1) = v(1,i)
                        y(4*i)     = v(2,i)
                    end do
                    der = dydt(t, y)
                    do i = 1, sizey
                        a(1,i) = der(4*i - 1)
                        a(2,i) = der(4*i)
                    end do
                end function get_a
        end subroutine leapfrog_DKD

        !---------------------------------------------------------------------------------------------
        ! CALLERS:
        !---------------------------------------------------------------------------------------------

        !---------------------------------
        !   Call an RK (not embedded) integrator
        !---------------------------------

        subroutine integ_caller (t, y, dt_min, dydt, integ, dt, ynew, check_fun)
            implicit none
            real(kind=8), intent(in) :: t, dt, dt_min
            real(kind=8), dimension(:), intent(in) :: y
            procedure(function_check_keep_tem), optional :: check_fun
            procedure(dydt_tem) :: dydt
            procedure(integ_tem) :: integ
            real(kind=8), dimension(size (y)), intent(out) :: ynew
            real(kind=8), dimension(size (y)) :: yaux, der
            real(kind=8) :: time, t_end, dt_used
            logical :: keep = .True.

            ynew    = y
            time    = t
            t_end   = time + dt
            dt_used = min (dt_min, dt)
            do while (time < t_end)
                if (present(check_fun)) then ! If Check Continue function present
                    keep = check_fun(ynew)
                    if (.not. keep) return ! Exit subroutine
                end if
                yaux = ynew
                der  = dydt (time, yaux)
                call integ (time, yaux, der, dt_used, dydt, ynew)
                time = time + dt_used
            end do
        end subroutine integ_caller

        !---------------------------------
        !   Call an embedded integrator
        !---------------------------------
        
        subroutine embedded_caller (t, y, dt_adap, dydt, integ, e_tol, beta, dt_min, dt, ynew, check_fun)
            implicit none
            real(kind=8), intent(in) :: t, e_tol, beta, dt_min, dt
            real(kind=8), dimension(:), intent(in) :: y
            procedure(function_check_keep_tem), optional :: check_fun
            real(kind=8), intent(inout) :: dt_adap
            procedure(dydt_tem) :: dydt
            procedure(embedded_tem) :: integ
            real(kind=8), dimension(size (y)), intent(out) :: ynew
            real(kind=8), dimension(size (y)) :: yaux, der
            real(kind=8) :: time, t_end, dtmin, dtused
            logical :: keep = .True.

            ynew  = y
            time  = t
            t_end = time + dt
            dtmin = min (dt_min, dt)
            do while (time < t_end)
                if (present(check_fun)) then ! If Check Continue function present
                    keep = check_fun(ynew)
                    if (.not. keep) then ! If Hard Exit is True
                        dt_adap = time - t ! Replace dt_adap with actual dt used
                        return ! Exit subroutine
                    end if
                end if
                yaux  = ynew
                dt_adap = min (dt_adap, t_end - time)
                der = dydt (time, yaux)
                call solve_embed (time, yaux, der, dt_adap, dydt, integ, e_tol, beta, dtmin, dtused, ynew)
                time = time + dtused
            end do
        end subroutine embedded_caller

        !---------------------------------
        !   Call an RK integrator with adaptive timestep
        !---------------------------------

        subroutine solve_rk_half_step_caller (t, y, dt_adap, dydt, integ, p, e_tol, beta, dt_min, dt, ynew, check_fun)
            implicit none
            real(kind=8), intent(in) :: t, e_tol, beta, dt, dt_min
            real(kind=8), dimension(:), intent(in) :: y
            procedure(function_check_keep_tem), optional :: check_fun
            real(kind=8), intent(inout) :: dt_adap
            procedure(dydt_tem) :: dydt
            procedure(integ_tem) :: integ
            integer(kind=4), intent(in) :: p
            real(kind=8), dimension(size (y)), intent(out) :: ynew
            real(kind=8), dimension(size (y)) :: yaux, der
            real(kind=8) :: time, t_end, dtmin, dtused
            logical :: keep = .True.

            ynew  = y
            time  = t
            t_end = time + dt
            dtmin = dt_min
            do while (time < t_end)
                if (present(check_fun)) then ! If Check Continue function present
                    keep = check_fun(ynew)
                    if (.not. keep) then ! If Hard Exit is True
                        dt_adap = time - t ! Replace dt_adap with actual dt used
                        return ! Exit subroutine
                    end if
                end if
                yaux = ynew
                dt_adap = min (dt_adap, t_end - time)
                der = dydt (time, yaux)
                call solve_rk_half_step (time, yaux, der, dt_adap, dydt, integ, p, e_tol, beta, dtmin, dtused, ynew)
                time = time + dtused
            end do
        end subroutine solve_rk_half_step_caller

        !---------------------------------
        !   Call Bulirsch_Stoer
        !---------------------------------
        
        subroutine BStoer_caller (t, y, dt_adap, dydt, e_tol, dt, ynew, check_fun) !, dt_min
            implicit none
            real(kind=8), intent(in) :: t, e_tol, dt!, dt_min
            real(kind=8), dimension(:), intent(in) :: y
            procedure(function_check_keep_tem), optional :: check_fun
            real(kind=8), intent(inout) :: dt_adap
            procedure(dydt_tem) :: dydt
            real(kind=8), dimension(size (y)), intent(out) :: ynew
            real(kind=8), dimension(size (y)) :: yaux
            real(kind=8) :: time, t_end, dtused! , dtmin
            logical :: keep = .True.

            ynew  = y
            time  = t
            t_end = time + dt
            do while (time < t_end)
                if (present(check_fun)) then ! If Check Continue function present
                    keep = check_fun(ynew)
                    if (.not. keep) then ! If Hard Exit is True
                        dt_adap = time - t ! Replace dt_adap with actual dt used
                        return ! Exit subroutine
                    end if
                end if
                yaux = ynew
                dt_adap = min(dt_adap, t_end - time)
                call Bulirsch_Stoer (time, yaux, dt_adap, dydt, e_tol, dtused, ynew)
                time = time + dtused    
            end do
        end subroutine BStoer_caller

        !---------------------------------
        !   Call Bulirsch_Stoer2 ! Only useful for [x, y, vx, vy,...] systems
        !---------------------------------

        subroutine BStoer_caller2 (t, y, dt_adap, dydt, e_tol, dt, ynew, check_fun) !,dt_min
            implicit none
            real(kind=8), intent(in) :: t, e_tol, dt!, dt_min
            real(kind=8), dimension(:), intent(in) :: y
            procedure(function_check_keep_tem), optional :: check_fun
            real(kind=8), intent(inout) :: dt_adap
            procedure(dydt_tem) :: dydt
            real(kind=8), dimension(size (y)), intent(out) :: ynew
            real(kind=8), dimension(size (y)) :: yaux
            real(kind=8) :: time, t_end, dtused!, dtmin
            logical :: keep = .True.

            ynew  = y
            time  = t
            t_end = time + dt
            ! dtmin = min (dt_min, dt)
            do while (time < t_end)
                if (present(check_fun)) then ! If Check Continue function present
                    keep = check_fun(ynew)
                    if (.not. keep) then ! If Hard Exit is True
                        dt_adap = time - t ! Replace dt_adap with actual dt used
                        return ! Exit subroutine
                    end if
                end if
                yaux  = ynew
                ! dt_adap = min(max(dt_adap, dt_min), t_end - time)
                dt_adap = min(dt_adap, t_end - time)
                call Bulirsch_Stoer2 (time, yaux, dt_adap, dydt, e_tol, dtused, ynew)
                time = time + dtused
            end do
        end subroutine BStoer_caller2

        !------------------------------------------------
        !  Call leapfrog integrator (adaptive timestep) ! Only useful for [x, y, vx, vy,...] systems
        !------------------------------------------------

        subroutine leapfrog_caller (t, y, dt_adap, dydt, leapfrog, e_tol, beta, dt_min, dt, ynew, check_fun)
            implicit none
            real(kind=8), intent(in) :: t, e_tol, beta, dt, dt_min
            real(kind=8), dimension(:), intent(in) :: y
            procedure(function_check_keep_tem), optional :: check_fun
            real(kind=8), intent(inout) :: dt_adap
            procedure(dydt_tem) :: dydt
            procedure(leapfrog_tem) :: leapfrog
            real(kind=8), dimension(size (y)), intent(out) :: ynew
            real(kind=8), dimension(size (y)) :: yaux, der
            real(kind=8) :: time, t_end, dtmin, dtused
            logical :: keep = .True.

            ynew  = y
            time  = t
            t_end = time + dt
            dtmin = dt_min
            do while (time < t_end)
                if (present(check_fun)) then ! If Check Continue function present
                    keep = check_fun(ynew)
                    if (.not. keep) then ! If Hard Exit is True
                        dt_adap = time - t ! Replace dt_adap with actual dt used
                        return ! Exit subroutine
                    end if
                end if
                yaux = ynew
                dt_adap = min (dt_adap, t_end - time)
                der = dydt (time, yaux)
                call solve_leapfrog (time, yaux, der, dt_adap, dydt, leapfrog, e_tol, beta, dtmin, dtused, ynew)
                time = time + dtused
            end do
        end subroutine leapfrog_caller

        !---------------------------------------------------------------------------------------------
end module integrators
