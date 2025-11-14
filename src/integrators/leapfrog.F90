!> Module with LeapFrog integrators (only useful for [X, V] systems)
module leapfrog
    use shared
    implicit none
    private
    public :: init_leapfrog, free_leapfrog, leapfrog_caller, leapfrog_fixed_caller

    ! Workspace arrays
    real(kind=8), allocatable :: y05(:), der05(:)  ! integrator
    real(kind=8), allocatable :: yscal(:)  ! solver
    real(kind=8), allocatable :: yaux(:)  ! solver
    real(kind=8), allocatable :: ycaller(:)  ! caller

    ! Constants
    integer(kind=4), parameter :: MAX_N_ITER = 350
    real(kind=8), parameter :: MAX_DT_FACTOR = 5.0d0


    ! Pointer to leapfrog used
    procedure (leapfrog_tem), pointer :: leapfrog_ptr => null ()

    abstract interface

        subroutine leapfrog_tem (sizey, y, dydt, t, dt, deri, ynew)
            import :: dydt_tem
            implicit none
            integer(kind=4), intent(in) :: sizey
            real(kind=8), dimension(sizey), intent(in) :: y
            procedure(dydt_tem) :: dydt
            real(kind=8), intent(in) :: t, dt
            real(kind=8), dimension(sizey), intent(in) :: deri
            real(kind=8), dimension(sizey), intent(out) :: ynew
        end subroutine leapfrog_tem
    
    end interface
    
    contains
    
        !!!! HANDLER

        subroutine init_leapfrog(sizey, which)
            implicit none
            integer(kind=4), intent(in) :: sizey
            integer(kind=4), intent(in) :: which

            if (MOD(sizey, 2) > 0) then
                print*, "ERROR: Can not use LeapFrog with uneven array."
                stop 1
            end if
            
            ! integrator
            allocate(y05(sizey), der05(sizey))

            ! solver
            allocate(yscal(sizey))
            allocate(yaux(sizey))

            ! caller
            allocate(ycaller(sizey))

            if (which > 0) then
                leapfrog_ptr => leapfrog_KDK
            else
                leapfrog_ptr => leapfrog_DKD
            end if
        end subroutine init_leapfrog

        subroutine free_leapfrog()
            implicit none
            if (allocated(y05)) deallocate(y05)
            if (allocated(der05)) deallocate(der05)
            if (allocated(yscal)) deallocate(yscal)
            if (allocated(yaux)) deallocate(yaux)
            if (allocated(ycaller)) deallocate(ycaller)
            nullify(leapfrog_ptr)
        end subroutine free_leapfrog

        
        !!!! Main subroutines for LeapFrog        

        subroutine leapfrog_KDK (sizey, y, dydt, t, dt, deri, ynew)
            implicit none
            integer(kind=4), intent(in) :: sizey
            real(kind=8), dimension(sizey), intent(in) :: y
            procedure(dydt_tem) :: dydt
            real(kind=8), intent(in) :: t, dt
            real(kind=8), dimension(sizey), intent(in) :: deri
            real(kind=8), dimension(sizey), intent(out) :: ynew

            integer(kind=4) :: i
            
            ! Initial set
            do i = 1, sizey, 2
                ! Calculate v05
                ynew(i+1) = y(i+1) + deri(i+1) * dt * C1_2
                ! Calculate x1 (y05) using v05
                ynew(i) = y(i) + ynew(i+1) * dt
            end do
            
            ! Calculate a_aux using accelerations at x1 and v05
            der05(:sizey) = dydt(t + dt, ynew)

            ! Update v1 using a_aux
            do i = 2, sizey, 2
                ynew(i) = ynew(i) + der05(i) * dt * C1_2
            end do
        end subroutine leapfrog_KDK

        !!! Leap Frog (DKD)

        subroutine leapfrog_DKD (sizey, y, dydt, t, dt, deri, ynew)
            implicit none
            integer(kind=4), intent(in) :: sizey
            real(kind=8), dimension(sizey), intent(in) :: y
            procedure(dydt_tem) :: dydt
            real(kind=8), intent(in) :: t, dt
            real(kind=8), dimension(sizey), intent(in) :: deri
            real(kind=8), dimension(sizey), intent(out) :: ynew

            integer(kind=4) :: i
            
            ! Initial set
            ynew = y
            
            ! Calculate x05 at the start of the step
            do i = 1, sizey, 2
                ynew(i) = ynew(i) + ynew(i+1) * dt * C1_2
            end do
            
            ! Calculate a_aux using accelerations at x05 and v0
            der05(:sizey) = dydt(t + dt, ynew)

            ! Update v1 and x1 (with x05 and v1)
            do i = 1, sizey, 2
                ynew(i+1) = ynew(i+1) + der05(i+1) * dt
                ynew(i) = ynew(i) + ynew(i+1) * dt * C1_2
            end do
        end subroutine leapfrog_DKD
        

        !------------------------------------------------
        !  Solver LeapFrog (adaptive timestep)
        !------------------------------------------------

        recursive subroutine solve_leapfrog (sizey, y, dydt, t, dt_adap, dt_used, deri, leapfrog, ynew)
            implicit none
            integer(kind=4), intent(in) :: sizey
            real(kind=8), dimension(sizey), intent(in) :: y
            procedure(dydt_tem) :: dydt
            real(kind=8), intent(in) :: t
            real(kind=8), intent(inout) :: dt_adap
            real(kind=8), intent(out) :: dt_used
            real(kind=8), dimension(sizey), intent(in) :: deri
            procedure (leapfrog_tem), pointer :: leapfrog
            real(kind=8), dimension(sizey), intent(out) :: ynew
            
            integer(kind=4), save :: iter = 0
            real(kind=8) :: e_calc, ratio, dt_half

            iter = iter + 1
            dt_adap = max (dt_adap, DT_MIN_NOW)
            dt_half = C1_2 * dt_adap

            ! yscal
            yscal(:sizey) = abs (y) + abs (dt_adap * deri) + SAFE_LOW

            ! y(t, dt) -> ynew
            call leapfrog (sizey, y, dydt, t, dt_adap, deri, ynew)
            
            ! y(t, dt/2) -> y05
            call leapfrog (sizey, y, dydt, t, dt_half, deri, y05(:sizey))

            ! d[y05(t + dt/2, dt/2)] / dt -> der05
            der05(:sizey) = dydt (t + dt_half, y05(:sizey))

            ! y05(t + dt/2, dt/2) -> yaux
            call leapfrog (sizey, y05(:sizey), dydt, t + dt_half, dt_half, der05(:sizey), yaux(:sizey))

            ! Error
            e_calc = max (maxval (abs ((ynew - yaux(:sizey)) / yscal(:sizey))), SAFE_LOW)
            ratio = E_TOL / e_calc

            if (ratio > ONE) then
                dt_used = dt_adap
                dt_adap = dt_adap * min (BETA * ratio**C1_3, MAX_DT_FACTOR) ! 1/3 porque es 1/(O(2) + 1)
                iter = 0

            else
                if (abs (dt_adap - DT_MIN_NOW) .le. E_TOL) then !E_TOL?
                    dt_used = DT_MIN_NOW
                    iter = 0

                else
                    dt_adap = dt_adap * min (BETA * ratio**C1_2, MAX_DT_FACTOR)

                    if ((dt_adap /= dt_adap) .or. (dt_adap .le. DT_MIN_NOW) .or. (iter == MAX_N_ITER)) then
                        dt_used = DT_MIN_NOW
                        dt_adap = DT_MIN_NOW

                        call leapfrog (sizey, y, dydt, t, dt_adap, deri, ynew)
                        iter = 0

                    else
                        call solve_leapfrog (sizey, y, dydt, t, dt_adap, dt_used, deri, leapfrog, ynew)
                                            
                    end if

                end if

            end if

        end subroutine solve_leapfrog


        ! ------------------------------------------------
        !  Call leapfrog integrator (adaptive timestep)
        ! ------------------------------------------------

        subroutine leapfrog_caller (t, y, dt_adap, dydt, dt, ynew, check_fun)
            implicit none
            real(kind=8), intent(in) :: t
            real(kind=8), dimension(:), intent(in) :: y
            real(kind=8), intent(inout) :: dt_adap  ! This try and also next try
            procedure(dydt_tem) :: dydt
            real(kind=8), intent(in) :: dt
            real(kind=8), dimension(size (y)), intent(out) :: ynew
            procedure(function_check_keep_tem), optional :: check_fun
            
            integer(kind=4) :: sizey
            real(kind=8) :: time, t_end, dt_used
            logical :: keep = .True.
            logical, save :: has_check = .False.

            sizey = size (y)
            has_check = present(check_fun)

            ynew = y
            time = t
            t_end = time + dt
            
            do while (time < t_end)

                if (has_check) then ! If Check Continue function present
                    keep = check_fun(ynew)
                    if (.not. keep) then ! If Hard Exit is True
                        dt_adap = time - t ! Replace dt_adap with actual dt used
                        return ! Exit subroutine
                    end if
                end if

                ycaller(:sizey) = ynew

                dt_adap = min (dt_adap, t_end - time)
                DT_MIN_NOW = min (DT_MIN, dt_adap)

                der(:sizey) = dydt (time, ycaller(:sizey))
                                    
                call solve_leapfrog (sizey, ycaller(:sizey), dydt, time, dt_adap, dt_used, &
                                    & der(:sizey), leapfrog_ptr, ynew)

                time = time + dt_used
            end do

        end subroutine leapfrog_caller


        ! ------------------------------------------------
        !  Call leapfrog integrator (fixed timestep)
        ! ------------------------------------------------

        subroutine leapfrog_fixed_caller (t, y, dt_adap, dydt, dt, ynew, check_fun)
            implicit none
            real(kind=8), intent(in) :: t
            real(kind=8), dimension(:), intent(in) :: y
            real(kind=8), intent(inout) :: dt_adap  ! This is each sub-step
            procedure(dydt_tem) :: dydt
            real(kind=8), intent(in) :: dt ! This is full step
            real(kind=8), dimension(size (y)), intent(out) :: ynew
            procedure(function_check_keep_tem), optional :: check_fun
            
            integer(kind=4) :: sizey
            real(kind=8) :: time, t_end
            logical :: keep = .True.
            logical, save :: has_check = .False.

            sizey = size (y)
            has_check = present(check_fun)

            ynew = y
            time = t
            t_end = time + dt
            
            do while (time < t_end)

                if (has_check) then ! If Check Continue function present
                    keep = check_fun(ynew)
                    if (.not. keep) then ! If Hard Exit is True
                        dt_adap = time - t ! Replace dt_adap with actual dt used
                        return ! Exit subroutine
                    end if
                end if

                ycaller(:sizey) = ynew

                dt_adap = min (dt_adap, t_end - time)
                DT_MIN_NOW = min (DT_MIN, dt_adap)

                der(:sizey) = dydt (time, ycaller(:sizey))

                call leapfrog_ptr (sizey, ycaller(:sizey), dydt, time, dt_adap, der(:sizey), ynew)

                time = time + dt_adap
            end do

        end subroutine leapfrog_fixed_caller


    
end module leapfrog