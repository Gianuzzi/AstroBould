!> Module with BStoer2 integrator (only useful for [x0, v0, x1, v1, ...] systems)
module bstoer2
    use shared
    implicit none
    private

    ! Workspace arrays
    real(kind=8), allocatable :: y0(:), der0(:)
    real(kind=8), allocatable :: arr_scal(:) 
    real(kind=8), allocatable :: yaux(:), deraux(:)
    real(kind=8), allocatable :: yend(:)
    real(kind=8), allocatable :: dy(:,:)

    ! Constants
    integer(kind=4), parameter :: MAX_N_ITER = 350
    
    contains
    
        !!!! HANDLER

        subroutine init_BS2(sizey)
            implicit none
            integer(kind=4), intent(in) :: sizey
            integer(kind=4) :: n_x2
            if (MOD(sizey, 2) > 0) then
                print*, "ERROR: Can not use BStoer2 with uneven array."
                stop 1
            end if

            if (allocated(y0)) deallocate(y0)
            if (allocated(der0)) deallocate(der0)
            if (allocated(arr_scal)) deallocate(arr_scal)
            if (allocated(yaux)) deallocate(yaux)
            if (allocated(deraux)) deallocate(deraux)
            if (allocated(yend)) deallocate(yend)
            if (allocated(dy)) deallocate(dy)

            allocate(y0(sizey), der0(sizey))

            ! arr_scal layout: EXTRA2 entries + 2*n_x (X and V scalings per body)
            n_x2 = int((sizey - EXTRA2) / NDIM2, 4) * 2
            allocate(arr_scal(EXTRA2 + n_x2))  ! |x0|, |v0|, |x1|, |v1|,..., |xA|, |vA|, |xB|, |vB|

            allocate(yaux(sizey), deraux(sizey))
            allocate(yend(sizey))
            allocate(dy(sizey,8))

        end subroutine init_BS2


        !!!! Auxiliar subroutines for Bulirsch_Stoer 2        

        subroutine bstep (sizey, y, dydt, t, htry, hdid, hnext)
            implicit none
            integer(kind=4), intent(in) :: sizey
            real(kind=8), dimension(sizey), intent(inout) :: y
            procedure(dydt_tem) :: dydt
            real(kind=8), intent(in) :: t, htry
            real(kind=8), intent(out) :: hdid, hnext

            real(kind=8), parameter :: SHRINK = 0.55d0
            real(kind=8), parameter :: GROW   = 1.3d0

            integer(kind=4) :: n_x
            integer(kind=4) :: auxi
            integer(kind=4) :: idim2
            integer(kind=4) :: i, j, j1, k, ns, ndo
            real(kind=8) :: tmp0, tmp1, tmp2
            real(kind=8) :: errmax
            real(kind=8) :: h, hx2, h2(8)
            real(kind=8) :: dt, tt

            ! NX2
            n_x = int((sizey - EXTRA2) / NDIM2, 4)

            !------------------------------------------------------------------------------
            ! Calculate accelerations at the start of the step
            der(:sizey) = dydt(t, y)

            y0(:sizey) = y
            der0(:sizey) = der

            ! Calculate the scale array
            do i = 1, EXTRA2
                arr_scal(i) = ONE / (abs(y0(i)) + SAFE_LOW)
            end do
            do i = 1, n_x
                auxi = get_index(i)
                tmp1 = ZERO
                tmp2 = ZERO
                do idim2 = 0, NDIM2 - 2, 2
                    tmp1 = tmp1 + y0(auxi + idim2)**2  ! X
                    tmp2 = tmp2 + y0(auxi + idim2 + 1)**2  ! V
                end do
                k = i + EXTRA2
                arr_scal(k) = ONE / sqrt(tmp1 + SAFE_LOW)  ! X
                arr_scal(k+1) = ONE / sqrt(tmp2 + SAFE_LOW)  ! V
            end do
            
            dt = htry
            ! Try to do the step using current stepsize
            do ndo = 0, MAX_N_ITER
                tt = t

                ! For each value of NS, do a modified-midpoint integration with 2N substeps
                do ns = 1, 8
                    h = dt * C12 / dble(ns)
                    h2(ns) = C14 / (ns * dble(ns))
                    hx2 = h * TWO
                    do i = 1, sizey, 2
                        yaux(i) = y0(i) + h * y0(i+1)
                        yaux(i+1) = y0(i+1) + h * der0(i+1)
                    end do
                    tt = tt + h

                    deraux(:sizey) = dydt(tt, yaux(:sizey))
                    do i = 1, sizey, 2
                        yend(i) = y0(i) + hx2 * yaux(i+1)
                        yend(i+1) = y0(i+1) + hx2 * deraux(i+1)
                    end do
                    tt = tt + hx2
                    
                    do j = 2, ns
                        deraux(:sizey) = dydt(tt, yend(:sizey))
                        do i = 1, sizey, 2
                            yaux(i) = yaux(i) + hx2 * yend(i+1)
                            yaux(i+1) = yaux(i+1) + hx2 * deraux(i+1)
                        end do
                        tt = tt + hx2
                        
                        deraux(:sizey) = dydt(tt, yaux(:sizey))
                        do i = 1, sizey, 2
                            yend(i) = yend(i) + hx2 * yaux(i+1)
                            yend(i+1) = yend(i+1) + hx2 * deraux(i+1)
                        end do
                        tt = tt + hx2
                    end do
                    
                    deraux(:sizey) = dydt(tt, yend(:sizey))
                    do i = 1, sizey, 2
                        dy(i,ns) = C12 * (yend(i) + yaux(i) + h * yend(i+1))
                        dy(i+1,ns) = C12 * (yend(i+1) + yaux(i+1) + h * deraux(i+1))
                    end do
            
                    ! Update the DX and DV arrays used for polynomial extrapolation
                    do j = ns - 1, 1, -1
                        j1 = j + 1
                        tmp0 = ONE / (h2(j) - h2(ns))
                        tmp1 = tmp0 * h2(j1)
                        tmp2 = tmp0 * h2(ns)
                        dy(:sizey,j) = tmp1 * dy(:sizey,j1) - tmp2 * dy(:sizey,j)
                    end do
            
                    ! After several integrations, test the relative error on extrapolated values
                    if (ns > 3) then

                        ! Maximum relative position and velocity errors (last D terms added)
                        errmax = ZERO

                        do i = 1, EXTRA2
                            errmax = MAX(errmax, abs(dy(i,1)) * arr_scal(i))
                        end do
                        do i = 1, n_x
                            auxi = get_index(i)
                            tmp1 = ZERO
                            tmp2 = ZERO
                            do idim2 = 0, NDIM2 - 2, 2
                                tmp1 = tmp1 + dy(auxi + idim2, 1)**2
                                tmp2 = tmp2 + dy(auxi + idim2 + 1, 1)**2
                            end do
                            k = i + EXTRA2
                            errmax = MAX(errmax, sqrt(tmp1) * arr_scal(k), sqrt(tmp2) * arr_scal(k+1))
                        end do

                        ! If error is smaller than TOL, update position and velocity arrays, and exit
                        if (errmax <= E_TOL) then ! * E_TOL?
                            y0(:sizey) = ZERO

                            do j = 1, ns
                                y0(:sizey) = y0(:sizey) + dy(:sizey,j)
                            end do

                            ! Recommend a stepsize for the next call to this subroutine
                            if (ns >= 8) hnext = dt * SHRINK
                            if (ns < 7)  hnext = dt * GROW

                            ! Update
                            hdid = dt
                            y = y0(:sizey)
                            
                            return

                        end if

                    end if

                end do
                ! If errors were too large, redo the step with half the previous step size.
                dt = dt * C12

            end do

            ! Check max iter
            if (ndo == MAX_N_ITER) then
                print*, "Warning, E_TOL not reached in bstep2"
                hdid = dt
                hnext = hdid
            end if

            ! Update
            hdid = dt
            y = y0(:sizey)

        end subroutine bstep



        !---------------------------------
        !   Call Bulirsch_Stoer 2
        !---------------------------------
        
        subroutine BStoer2_caller (t, y, dt_adap, dydt, dt, ynew, check_fun)
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

                dt_adap = min(dt_adap, t_end - time)                
                dt_try = dt_adap

                call bstep (sizey, ynew, dydt, time, dt_try, dt_used, dt_adap)

                time = time + dt_used

            end do

        end subroutine BStoer2_caller

    
end module bstoer2