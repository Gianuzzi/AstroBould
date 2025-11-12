!> Module with LeapFrog integrators (only useful for [X, V] systems)
module leapfrog
    use shared
    implicit none
    private

    ! Workspace arrays
    real(kind=8), allocatable :: x0(:,:), v0(:,:), a0(:,:)
    real(kind=8), allocatable :: x1(:,:), v1(:,:), aux05(:,:)
    real(kind=8), allocatable :: a(:,:)
    real(kind=8), allocatable :: der05(:)
    real(kind=8), allocatable :: y05(:), yaux(:), yscal(:)

    ! Constants
    integer(kind=4), parameter :: MAX_N_ITER = 350
    real(kind=8), parameter :: MAX_DT_FACTOR = 5.0d0


    ! Pointer to leapfrog used
    procedure (leapfrog_tem), pointer :: leapfrog_ptr => null ()

    abstract interface

        subroutine leapfrog_tem (sizey, y, dydt, size_nd, t, dt, deri, ynew)
            import :: dydt_tem
            implicit none
            integer(kind=4), intent(in) :: sizey
            real(kind=8), dimension(sizey), intent(in) :: y
            procedure(dydt_tem) :: dydt
            integer(kind=4), intent(in) :: size_nd
            real(kind=8), intent(in) :: t, dt
            real(kind=8), dimension(sizey), intent(in) :: deri
            real(kind=8), dimension(sizey), intent(out) :: ynew
        end subroutine leapfrog_tem
    
    end interface
    
    contains
    
        !!!! HANDLER

        subroutine init_leapfrof(sizey, which)
            implicit none
            integer(kind=4), intent(in) :: sizey
            integer(kind=4), intent(in) :: which
            integer(kind=4) :: size_nd

            size_nd = int(sizey / ndim2, 4)

            allocate(x0(ndim,size_nd), v0(ndim,size_nd), a0(ndim,size_nd))
            allocate(x1(ndim,size_nd), v1(ndim,size_nd), aux05(ndim,size_nd))
            allocate(a(ndim,size_nd))
            allocate(der05(sizey))
            allocate(y05(sizey), yaux(sizey), yscal(sizey))

            if (which > 0) then
                leapfrog_ptr => leapfrog_KDK
            else
                leapfrog_ptr => leapfrog_DKD
            end if
        end subroutine init_leapfrof

        
        !!!! Main subroutines for LeapFrog        

        subroutine leapfrog_KDK (sizey, y, dydt, size_nd, t, dt, deri, ynew)
            implicit none
            integer(kind=4), intent(in) :: sizey
            real(kind=8), dimension(sizey), intent(in) :: y
            procedure(dydt_tem) :: dydt
            integer(kind=4), intent(in) :: size_nd
            real(kind=8), intent(in) :: t, dt
            real(kind=8), dimension(sizey), intent(in) :: deri
            real(kind=8), dimension(sizey), intent(out) :: ynew

            integer(kind=4) :: auxi, auxidx
            integer(kind=4) :: idim, i
            
            ! Initial set
            ynew = y

            
            ! Get x0, v0, and a0
            do i = 1, size_nd
                auxi = get_index(i)
                do idim = 1, ndim
                    auxidx = auxi + 2 * (idim - 1)
                    x0(idim,i) = y(auxidx)
                    v0(idim,i) = y(auxidx + 1)
                    a0(idim,i) = deri(auxidx + 1)
                end do
            end do

            ! Calculate v05
            aux05(:,:size_nd) = v0(:,:size_nd) + a0(:,:size_nd) * dt * C12

            ! Calculate x1 using v05
            x1(:,:size_nd) = x0(:,:size_nd) + aux05(:,:size_nd) * dt
            
            ! Calculate a and v1 using accelerations at x1 and v05
            call get_a (sizey, size_nd, t + dt, x1(:,:size_nd), aux05(:,:size_nd), dydt, a(:,:size_nd))
            v1(:,:size_nd) = aux05(:,:size_nd) + a(:,:size_nd) * dt * C12

            ! Update ynew
            do i = 1, size_nd
                auxi = get_index(i)
                do idim = 1, ndim
                    auxidx = auxi + 2 * (idim - 1)
                    ynew(auxidx) = x1(idim,i)
                    ynew(auxidx + 1) = v1(idim,i)
                end do
            end do

        end subroutine leapfrog_KDK

        !!! Leap Frog (DKD)

        subroutine leapfrog_DKD (sizey, y, dydt, size_nd, t, dt, deri, ynew)
            implicit none
            integer(kind=4), intent(in) :: sizey
            real(kind=8), dimension(sizey), intent(in) :: y
            procedure(dydt_tem) :: dydt
            integer(kind=4), intent(in) :: size_nd
            real(kind=8), intent(in) :: t, dt
            real(kind=8), dimension(sizey), intent(in) :: deri
            real(kind=8), dimension(sizey), intent(out) :: ynew

            integer(kind=4) :: auxi, auxidx
            integer(kind=4) :: idim, i
            
            ! Initial set
            ynew = y

            
            ! Get x0, v0
            do i = 1, size_nd
                auxi = get_index(i)
                do idim = 1, ndim
                    auxidx = auxi + 2 * (idim - 1)
                    x0(idim,i) = y(auxidx)
                    v0(idim,i) = y(auxidx + 1)
                    a0(idim,i) = deri(auxidx + 1)
                end do
            end do

            ! Calculate x05 at the start of the step
            aux05(:,:size_nd) = x0(:,:size_nd) + v0(:,:size_nd) * dt * C12
            
            ! Calculate a and v1 using accelerations at x05 and v0
            call get_a (sizey, size_nd, t + dt, aux05(:,:size_nd), v0(:,:size_nd), dydt, a(:,:size_nd))
            v1(:,:size_nd) = v0(:,:size_nd) + a(:,:size_nd) * dt

            ! Calculate x1 using x05 and v1
            x1(:,:size_nd) = aux05(:,:size_nd) + v1(:,:size_nd) * dt * C12

            ! Update ynew
            do i = 1, size_nd
                auxi = get_index(i)
                do idim = 1, ndim
                    auxidx = auxi + 2 * (idim - 1)
                    ynew(auxidx) = x1(idim,i)
                    ynew(auxidx + 1) = v1(idim,i)
                end do
            end do

        end subroutine leapfrog_DKD
        

        !------------------------------------------------
        !  Solver LeapFrog (adaptive timestep)
        !------------------------------------------------

        recursive subroutine solve_leapfrog (sizey, y, dydt, size_nd, t, dt_adap, dt_used, deri, leapfrog, ynew)
            implicit none
            integer(kind=4), intent(in)                         :: sizey
            real(kind=8), dimension(sizey), intent(in)          :: y
            procedure(dydt_tem)                                 :: dydt
            integer(kind=4), intent(in)                         :: size_nd
            real(kind=8), intent(in)                            :: t
            real(kind=8), intent(inout)                         :: dt_adap
            real(kind=8), intent(out)                           :: dt_used
            real(kind=8), dimension(sizey), intent(in)          :: deri
            procedure (leapfrog_tem), pointer                   :: leapfrog
            real(kind=8), dimension(sizey), intent(out)         :: ynew
            
            integer(kind=4), save                               :: iter = 0
            real(kind=8)                                        :: e_calc, ratio, dt_half

            iter = iter + 1
            dt_adap = max (dt_adap, DT_MIN)
            dt_half = C12 * dt_adap

            ! yscal
            yscal(:sizey) = abs (y) + abs (dt_adap * deri) + SAFE_LOW

            ! ynew
            call leapfrog (sizey, y, dydt, size_nd, t, dt_adap, deri, ynew)

            ! y05
            call leapfrog (sizey, y, dydt, size_nd, t, dt_half, deri, y05(:sizey))

            ! der05
            der05(:sizey) = dydt (t + dt_half, y05(:sizey))

            ! yaux
            call leapfrog (sizey, y, dydt, size_nd, t + dt_half, dt_half, der05(:sizey), yaux(:sizey))

            
            e_calc = max (maxval (abs ((ynew(:sizey) - yaux(:sizey)) / yscal(:sizey))), SAFE_LOW)
            ratio  = E_TOL / e_calc

            if (ratio > ONE) then
                dt_used = dt_adap
                dt_adap = dt_adap * min (BETA * ratio**C13, MAX_DT_FACTOR) ! 1/3 porque es 1/(O(2) + 1)
                iter = 0

            else

                if (abs (dt_adap - DT_MIN) .le. E_TOL) then !E_TOL?
                    dt_used = DT_MIN
                    iter = 0

                else
                    dt_adap = dt_adap * min (BETA * ratio**C12, MAX_DT_FACTOR)

                    if ((dt_adap /= dt_adap) .or. (dt_adap .le. DT_MIN) .or. (iter == MAX_N_ITER)) then
                        dt_used = DT_MIN
                        dt_adap = DT_MIN

                        ! ynew; with new dt_adap
                        call leapfrog (sizey, y, dydt, size_nd, t, dt_adap, deri, ynew)
                        iter = 0

                    else
                        call solve_leapfrog (sizey, y, dydt, size_nd, t, dt_adap, dt_used, deri, leapfrog, ynew)
                                            
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
            
            integer(kind=4) :: sizey, size_nd
            real(kind=8) :: time, t_end, dt_used
            logical :: keep = .True.
            logical, save :: has_check = .False.

            sizey = size (y)
            size_nd = int(sizey / ndim2, 4)
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

                yaux(:sizey) = ynew
                dt_adap = min (dt_adap, t_end - time)
                der(:sizey) = dydt (time, yaux(:sizey))
                                    
                call solve_leapfrog (sizey, yaux(:sizey), dydt, size_nd, time, dt_adap, dt_used, &
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
            
            integer(kind=4) :: sizey, size_nd
            real(kind=8) :: time, t_end
            logical :: keep = .True.
            logical, save :: has_check = .False.

            sizey = size (y)
            size_nd = int(sizey / ndim2, 4)
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

                yaux(:sizey) = ynew
                dt_adap = min (dt_adap, t_end - time)
                der(:sizey) = dydt (time, yaux(:sizey))
                call leapfrog_ptr (sizey, y, dydt, size_nd, time, dt_adap, der(:sizey), ynew)

                time = time + dt_adap
            end do

        end subroutine leapfrog_fixed_caller


    
end module leapfrog