!> Module with BStoer integrator
module bstoer
    use shared
    implicit none
    private
    public :: init_BS, free_BS, BStoer_caller

    ! Workspace arrays
    real(kind=8), allocatable :: ysav(:), yseq(:), yerr(:), yscal(:) ! bstep
    real(kind=8), allocatable :: qcolpz(:,:)  ! bstep
    real(kind=8), allocatable :: mmid_ym(:), mmid_yn(:)  ! mmid
    real(kind=8), allocatable :: pz_d(:)  ! pzextr
    
    contains
    
        !!!! HANDLERS

        subroutine init_BS(sizey)
            implicit none
            integer(kind=4), intent(in) :: sizey
            
            ! bstep
            allocate(ysav(sizey), yseq(sizey), yerr(sizey), yscal(sizey))
            allocate(qcolpz(sizey,16))

            ! mmid
            allocate(mmid_ym(sizey), mmid_yn(sizey))

            ! pzextr
            allocate(pz_d(sizey))
        end subroutine init_BS

        subroutine free_BS()
            implicit none
            
            ! bstep
            if (allocated(ysav)) deallocate(ysav)
            if (allocated(yseq)) deallocate(yseq)
            if (allocated(yerr)) deallocate(yerr)
            if (allocated(yscal)) deallocate(yscal)
            if (allocated(qcolpz)) deallocate(qcolpz)

            ! mmid
            if (allocated(mmid_ym)) deallocate(mmid_ym)
            if (allocated(mmid_yn)) deallocate(mmid_yn)

            ! pzextr
            if (allocated(pz_d)) deallocate(pz_d)
        end subroutine free_BS


        !!!! Auxiliar subroutines for Bulirsch_Stoer

        subroutine bstep (sizey, y, dydt, x, htry, hdid, hnext)
            implicit none
            integer(kind=4), intent(in) :: sizey
            real(kind=8), dimension(sizey), intent(inout) :: y
            procedure(dydt_tem) :: dydt
            real(kind=8), intent(in) :: x
            real(kind=8), intent(in) :: htry
            real(kind=8), intent(out) :: hdid, hnext
            
            real(kind=8), parameter :: safe1 = .25d0, safe2 = .7d0
            real(kind=8), parameter :: redmax = 1.d-5, redmin = .7d0
            real(kind=8), parameter :: tini = 1.d-30, scalmx = .1d0
            integer(kind=4), parameter, dimension(9) :: nseq = (/2, 4, 6, 8, 10, 12, 14, 16, 18/)
            integer(kind=4), save :: kmax, kopt
            real(kind=8), dimension(8, 8), save :: alf
            real(kind=8), dimension(8) :: err
            integer(kind=4), dimension(9), save :: arr ! a in NR F90
            real(kind=8), save :: xnew, epsold = -1.d0
            real(kind=8) :: eps1, errmax, fact, h, red, scala, wrkmin, xest
            logical :: reduct
            logical, save :: first = .True.
            real(kind=8), dimension(16) :: xpz
            real(kind=8) :: work
            integer(kind=4) :: k, iq, i ,km, kk

            der(:sizey) = dydt (x, y)
            yscal(:sizey) = abs (y) + abs (htry * der(:sizey)) + SAFE_LOW
            
            if (abs(E_TOL - epsold) > SAFE_LOW) then !E_TOL? ! A new tolerance, so reinitialize.
                hnext = -1.0d29 ! “Impossible” values.
                xnew  = -1.0d29
                eps1  = safe1 * E_TOL
                arr(1) = nseq(1) + 1
                do k = 1, 8
                    arr(k + 1) = arr(k) + nseq(k + 1)
                end do
                ! Compute α(k, q):
                do iq = 2, 8
                    do k = 1, iq - 1
                        alf(k, iq) = eps1**((arr(k + 1) - arr(iq + 1)) / ((arr(iq + 1) - arr(1) + ONE) * (TWO * k + 1)))
                    end do
                end do
                epsold = E_TOL
                do kopt = 2, 7 ! Determine optimal row number for convergence.
                    if (arr(kopt + 1) > arr(kopt) * alf(kopt - 1, kopt)) exit
                end do
                kmax = kopt
            end if
            h = htry
            ysav(:sizey) = y
            if ((abs(h - hnext) > SAFE_LOW) .or. (abs(x - xnew) > SAFE_LOW)) then !E_TOL?
                ! A new stepsize or a new integration: Reestablish the order window.
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
                    call mmid (sizey, ysav(:sizey), der(:sizey), x, h, nseq(k), yseq(:sizey), dydt)
                    yscal(:sizey) = abs (y) + abs (h * der(:sizey)) + SAFE_LOW
                    xest = (h / dble(nseq(k)))**2 ! Squared, since error series is even.
                    call pzextr (sizey, k, xest, yseq(:sizey), y, yerr(:sizey), qcolpz(1:sizey, :), xpz) ! Perform extrapolation.
                    if (k /= 1) then ! Compute normalized error estimate eps(k).
                        errmax = tini
                        do i = 1, sizey
                            errmax = max (errmax, abs(yerr(i)/yscal(i))) ! FALTA SCALADO? !! Ya no?
                        end do
                        errmax = errmax / E_TOL
                        km = k - 1
                        err(km) = (errmax / safe1)**(ONE / (TWO * km + 1))
                    end if
                    if (k /= 1 .and. (k .ge. kopt - 1 .or. first)) then ! In order window.
                        if (errmax < ONE) exit main_loop ! Converged.
                        if (k == kmax .or. k == kopt + 1) then ! Check for possible stepsize reduction.
                            red = safe2 / err(km)
                            exit
                        else if (k == kopt) then
                            if (alf(kopt - 1, kopt) < err(km)) then
                                red = ONE / err(km)
                                exit
                            end if
                        else if (kopt == kmax) then
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
            ! x = xnew  !! x now is intent(in)
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
            ! Check for possible order increase, but not if stepsize was just reduced.
            if (kopt .ge. k .and. kopt /= kmax .and. .not. reduct) then
                fact = max (scala / alf(kopt - 1, kopt), scalmx)
                if (arr(kopt + 1) * fact .le. wrkmin) then
                    hnext = h / fact
                    kopt = kopt + 1
                end if
            end if
        end subroutine bstep

        subroutine mmid (sizey, y, dydx, xs, htot, nstep, yout, dydt)
            integer(kind=4), intent(in) :: sizey
            real(kind=8), dimension(sizey), intent(in) :: y, dydx
            real(kind=8), intent(in) :: xs, htot
            integer(kind=4), intent(in) :: nstep
            real(kind=8), dimension(sizey), intent(out) :: yout            
            procedure(dydt_tem) :: dydt

            real(kind=8) :: x, h, h2, swap
            integer(kind=4) :: i, n 

            h = htot / (nstep * ONE) ! Stepsize this trip.
            mmid_ym(1:sizey) = y
            mmid_yn(1:sizey) = y + h * dydx ! First step.
            x = xs + h
            yout = dydt(x, mmid_yn) ! Will use yout for temporary storage of derivatives.
            h2 = TWO * h
            do n = 2, nstep ! General step.
                do i = 1, sizey
                    swap = mmid_ym(i) + h2 * yout(i)
                    mmid_ym(i) = mmid_yn(i)
                    mmid_yn(i) = swap
                end do
                x = x + h
                yout = dydt(x, mmid_yn(1:sizey))
            end do
            yout = C1_2 * (mmid_ym(1:sizey) + mmid_yn(1:sizey) + h * yout) ! Last step.
        end subroutine mmid

        subroutine pzextr (sizey, iest, xest, yest, yz, dy, qcol, x)
            implicit none
            integer(kind=4), intent(in) :: sizey
            integer(kind=4), intent(in) :: iest
            real(kind=8), intent(in) :: xest
            real(kind=8), dimension(sizey), intent(in) :: yest
            real(kind=8), dimension(sizey), intent(out) :: yz, dy
            real(kind=8), dimension(sizey,16), intent(inout) :: qcol
            real(kind=8), dimension(16), intent(inout) :: x
            
            integer(kind=4) :: j, k1
            real(kind=8)    :: delta, f1, f2, q

            x(iest) = xest  ! Save current independent variable.
            dy = yest
            yz = yest
            if (iest == 1) then  ! Store ﬁrst estimate in ﬁrst column.
                qcol(:, 1) = yest(:)
            else
                pz_d = yest
                do k1 = 1, iest - 1
                    delta = ONE / (x(iest - k1) - xest)
                    f1 = xest * delta
                    f2 = x(iest - k1) * delta
                    do j = 1, sizey ! Propagate tableau 1 diagonal more.
                        q = qcol(j, k1)
                        qcol(j, k1) = dy(j)
                        delta = pz_d(j) - q
                        dy(j) = f1 * delta
                        pz_d(j)  = f2 * delta
                        yz(j) = yz(j) + dy(j)
                    end do  
                end do
                qcol(:, iest) = dy(:)
            end if
        end subroutine pzextr


        !---------------------------------
        !   Call Bulirsch_Stoer
        !---------------------------------
        
        subroutine BStoer_caller (t, y, dt_adap, dydt, dt, ynew, check_fun)
            implicit none
            real(kind=8), intent(in) :: t
            real(kind=8), dimension(:), intent(in) :: y
            real(kind=8), intent(inout) :: dt_adap  ! This try and also next try
            procedure(dydt_tem) :: dydt
            real(kind=8), intent(in) :: dt
            real(kind=8), dimension(size (y)), intent(out) :: ynew
            procedure(function_check_keep_tem), optional :: check_fun
            
            integer(kind=4) :: sizey
            real(kind=8) :: time, t_end, dt_try, dt_used
            logical :: keep = .True.
            logical :: has_check = .False.

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

                dt_adap = min(dt_adap, t_end - time)                
                dt_try = dt_adap

                call bstep (sizey, ynew, dydt, time, dt_try, dt_used, dt_adap)

                time = time + dt_used

            end do

        end subroutine BStoer_caller

    
end module bstoer