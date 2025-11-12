!> Module with BStoer2 integrator (only useful for [X, V] systems)
module bstoer2
    use shared
    implicit none
    private

    ! Workspace arrays
    real(kind=8), allocatable :: x(:,:), v(:,:), a(:,:)
    real(kind=8), allocatable :: x0(:,:), v0(:,:), a0(:,:)
    real(kind=8), allocatable :: xend(:,:), vend(:,:)
    real(kind=8), allocatable :: xscal(:), vscal(:)
    real(kind=8), allocatable :: dx(:,:,:), dv(:,:,:)

    ! Constants
    integer(kind=4), parameter :: MAX_N_ITER = 350
    
    contains
    
        !!!! HANDLER

        subroutine init_BS2(sizey)
            implicit none
            integer(kind=4), intent(in) :: sizey
            integer(kind=4) :: size_nd

            size_nd = int(sizey / ndim2, 4)

            allocate(x(ndim,size_nd), v(ndim,size_nd), a(ndim,size_nd))
            allocate(x0(ndim,size_nd), v0(ndim,size_nd), a0(ndim,size_nd))
            allocate(xend(ndim,size_nd), vend(ndim,size_nd))
            allocate(xscal(size_nd), vscal(size_nd))
            allocate(dx(ndim,size_nd,8), dv(ndim,size_nd,8))

        end subroutine init_BS2


        !!!! Auxiliar subroutines for Bulirsch_Stoer 2        

        subroutine bstep (sizey, y, dydt, size_nd, t, htry, hdid, hnext)
            implicit none
            integer(kind=4), intent(in) :: sizey
            real(kind=8), dimension(sizey), intent(inout) :: y
            procedure(dydt_tem) :: dydt
            integer(kind=4), intent(in) :: size_nd
            real(kind=8), intent(in) :: t, htry
            real(kind=8), intent(out) :: hdid, hnext

            real(kind=8), parameter :: SHRINK = 0.55d0
            real(kind=8), parameter :: GROW   = 1.3d0

            integer(kind=4) :: auxi, auxidx
            integer(kind=4) :: idim, i, j, j1, k, ns, ndo
            real(kind=8) :: tmp0, tmp1, tmp2
            real(kind=8) :: errmax
            real(kind=8) :: h, hx2, h2(8)
            real(kind=8) :: dt, tt

            !------------------------------------------------------------------------------
            ! Calculate accelerations at the start of the step
            der(:sizey) = dydt(t, y)

            do i = 1, size_nd
                auxi = get_index(i)
                do idim = 1, ndim
                    auxidx = auxi + 2 * (idim - 1)
                    x0(idim,i) = y(auxidx)
                    v0(idim,i) = y(auxidx + 1)
                    a0(idim,i) = der(auxidx + 1)
                end do
                ! Arrays used to scale the relative error (1 / absolute distance or velocity)
                xscal(i) = ONE / sqrt(sum(x0(:,i)**2))
                vscal(i) = ONE / sqrt(sum(v0(:,i)**2))
            end do
            
            dt = htry
            ! Try to do the step using current stepsize
            do ndo = 0, MAX_N_ITER
                tt = t

                ! For each value of NS, do a modified-midpoint integration with 2N substeps
                do ns = 1, 8
                    h = dt / (TWO * dble(ns))
                    h2(ns) = C14 / (ns * dble(ns))
                    hx2 = h * TWO
                    x(:,:size_nd) = x0(:,:size_nd) + h * v0(:,:size_nd)
                    v(:,:size_nd) = v0(:,:size_nd) + h * a0(:,:size_nd)
                    tt = tt + h

                    call get_a(sizey, size_nd, tt, x(:,:size_nd), v(:,:size_nd), dydt, a(:,:size_nd))
                    xend(:,:size_nd) = x0(:,:size_nd) + hx2 * v(:,:size_nd)
                    vend(:,:size_nd) = v0(:,:size_nd) + hx2 * a(:,:size_nd)
                    tt = tt + hx2
                    
                    do j = 2, ns
                        call get_a(sizey, size_nd, tt, xend(:,:size_nd), vend(:,:size_nd), dydt, a(:,:size_nd))
                        x(:,:size_nd) = x(:,:size_nd) + hx2 * vend(:,:size_nd)
                        v(:,:size_nd) = v(:,:size_nd) + hx2 * a(:,:size_nd)
                        tt = tt + hx2
                        
                        call get_a(sizey, size_nd, tt, x(:,:size_nd), v(:,:size_nd), dydt, a(:,:size_nd))
                        xend(:,:size_nd) = xend(:,:size_nd) + hx2 * v(:,:size_nd)
                        vend(:,:size_nd) = vend(:,:size_nd) + hx2 * a(:,:size_nd)
                        tt = tt + hx2
                    end do
                    
                    call get_a(sizey, size_nd, tt, xend(:,:size_nd), vend(:,:size_nd), dydt, a(:,:size_nd))
                    dx(:,:size_nd,ns) = C12 * (xend(:,:size_nd) + x(:,:size_nd) + h * vend(:,:size_nd))
                    dv(:,:size_nd,ns) = C12 * (vend(:,:size_nd) + v(:,:size_nd) + h * a(:,:size_nd))
            
                    ! Update the DX and DV arrays used for polynomial extrapolation
                    do j = ns - 1, 1, -1
                        j1 = j + 1
                        tmp0 = ONE / (h2(j) - h2(ns))
                        tmp1 = tmp0 * h2(j1)
                        tmp2 = tmp0 * h2(ns)
                        dx(:,:size_nd,j) = tmp1 * dx(:,:size_nd,j1) - tmp2 * dx(:,:size_nd,j)
                        dv(:,:size_nd,j) = tmp1 * dv(:,:size_nd,j1) - tmp2 * dv(:,:size_nd,j)
                    end do
            
                    ! After several integrations, test the relative error on extrapolated values
                    if (ns > 3) then

                        ! Maximum relative position and velocity errors (last D terms added)
                        errmax = ZERO
                        do k = 1, size_nd
                            tmp1 = sqrt(sum(dx(:,k,1)**2))
                            tmp2 = sqrt(sum(dv(:,k,1)**2))
                            errmax = max(errmax, tmp1*xscal(k), tmp2*vscal(k))
                        end do

                        ! If error is smaller than TOL, update position and velocity arrays, and exit
                        if (errmax <= E_TOL) then ! * E_TOL?
                            x0(:,:size_nd) = ZERO
                            v0(:,:size_nd) = ZERO

                            do j = 1, ns
                                x0(:,:size_nd) = x0(:,:size_nd) + dx(:,:size_nd,j)
                                v0(:,:size_nd) = v0(:,:size_nd) + dv(:,:size_nd,j)
                            end do

                            ! Recommend a stepsize for the next call to this subroutine
                            if (ns >= 8) hnext = dt * SHRINK
                            if (ns < 7)  hnext = dt * GROW

                            ! Update
                            hdid = dt
                            do i = 1, size_nd
                                auxi = get_index(i)
                                do idim = 1, ndim
                                    auxidx = auxi + 2 * (idim - 1)
                                    y(auxidx) = x0(idim,i)
                                    y(auxidx + 1) = v0(idim,i)
                                end do
                            end do
                            
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
                hdid = dt * TWO
                hnext = hdid
            end if

            ! Update
            hdid = dt
            do i = 1, size_nd
                auxi = get_index(i)
                do idim = 1, ndim
                    auxidx = auxi + 2 * (idim - 1)
                    y(auxidx) = x0(idim,i)
                    y(auxidx + 1) = v0(idim,i)
                end do
            end do

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
            
            integer(kind=4) :: sizey, size_nd
            real(kind=8) :: time, t_end, dt_try, dt_used
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

                dt_adap = min(dt_adap, t_end - time)                
                dt_try = dt_adap

                call bstep (sizey, ynew, dydt, size_nd, time, dt_try, dt_used, dt_adap)

                time = time + dt_used

            end do

        end subroutine BStoer2_caller

    
end module bstoer2