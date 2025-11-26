!> Module with BStoer2 integrator (standard only useful for [x0, v0, x1, v1, ...] systems)
module bstoer2
    use shared
    implicit none
    private
    public :: init_BS2, free_BS2, BStoer2_caller

    ! Workspace arrays
    real(wp), allocatable :: y0(:), der0(:)
    real(wp), allocatable :: arr_scal(:) 
    real(wp), allocatable :: yaux(:), deraux(:)
    real(wp), allocatable :: yend(:)
    real(wp), allocatable :: dy(:,:)

    ! Constants
    integer(kind=4), parameter :: MAX_N_ITER = 350


    ! Pointer to bstep used
    procedure (bstep_tem), pointer :: bstep_ptr => null ()

    abstract interface

        subroutine bstep_tem (sizey, y, dydt, t, htry, hdid, hnext)
            import :: wp
            import :: dydt_tem
            implicit none
            integer(kind=4), intent(in) :: sizey
            real(wp), dimension(sizey), intent(inout) :: y
            procedure(dydt_tem) :: dydt
            real(wp), intent(in) :: t, htry
            real(wp), intent(out) :: hdid, hnext
        end subroutine bstep_tem
    
    end interface
    
    contains
    
        !!!! HANDLER

        subroutine init_BS2(sizey, standard)
            implicit none
            integer(kind=4), intent(in) :: sizey
            logical, intent(in), optional :: standard
            logical :: is_std  = .False.
            integer(kind=4) :: n_x2

            if (present(standard)) is_std = standard

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
            if (n_x2 .le. 0) then
                print*, "ERROR: Bad n_x2:", n_x2, "Check NEXTRA and NDIM."
                stop 1
            end if
            allocate(arr_scal(EXTRA2 + n_x2))  ! |x0|, |v0|, |x1|, |v1|,..., |xA|, |vA|, |xB|, |vB|

            allocate(yaux(sizey), deraux(sizey))
            allocate(yend(sizey))
            allocate(dy(sizey,8))

            ! Set pointer
            if (is_std) then
                bstep_ptr => bstep_std
            else
                bstep_ptr => bstep
            end if

        end subroutine init_BS2

        subroutine free_BS2()
            implicit none
            if (allocated(y0)) deallocate(y0)
            if (allocated(der0)) deallocate(der0)
            if (allocated(arr_scal)) deallocate(arr_scal)
            if (allocated(yaux)) deallocate(yaux)
            if (allocated(yend)) deallocate(yend)
            if (allocated(dy)) deallocate(dy)
        end subroutine free_BS2


        !!!! Auxiliar subroutines for Bulirsch_Stoer 2
        
        subroutine bstep(sizey, y, dydt, t, htry, hdid, hnext)
            implicit none
            integer(kind=4), intent(in) :: sizey
            real(wp), dimension(sizey), intent(inout) :: y
            procedure(dydt_tem) :: dydt
            real(wp), intent(in) :: t, htry
            real(wp), intent(out) :: hdid, hnext

            real(wp), parameter :: SHRINK = 0.55e0_wp
            real(wp), parameter :: GROW   = 1.3e0_wp

            integer(kind=4) :: n_x
            integer(kind=4) :: auxi
            integer(kind=4) :: i, j, j1, k, ns, ndo, dim_idx
            real(wp) :: tmp0, tmp1, tmp2
            real(wp) :: errmax
            real(wp) :: h, hx2, h2(8)
            real(wp) :: dt, tt

            ! NX2 - number of NDIM objects
            n_x = int((sizey - EXTRA2) / NDIM2, 4)

            !------------------------------------------------------------------------------
            ! Calculate accelerations at the start of the step
            der(:sizey) = dydt(t, y)

            y0(:sizey) = y
            der0(:sizey) = der

            ! Calculate the scale array
            ! 1D variables
            do i = 1, EXTRA2
                arr_scal(i) = ONE / (abs(y0(i)) + SAFE_LOW)
            end do
            
            ! NDIM variables
            do i = 1, n_x
                auxi = EXTRA2 + (i-1)*NDIM2 + 1  ! Starting index for this NDIM object
                
                ! Calculate norms for positions and velocities
                tmp1 = ZERO
                tmp2 = ZERO
                do dim_idx = 0, NDIM-1
                    tmp1 = tmp1 + y0(auxi + dim_idx)**2  ! Positions
                    tmp2 = tmp2 + y0(auxi + NDIM + dim_idx)**2  ! Velocities
                end do
                
                ! Store scales for this NDIM object
                k = EXTRA2 + (i-1)*2 + 1  ! Index in arr_scal for this object
                arr_scal(k) = ONE / (sqrt(tmp1) + SAFE_LOW)  ! Position scale
                arr_scal(k+1) = ONE / (sqrt(tmp2) + SAFE_LOW)  ! Velocity scale
            end do
            
            dt = htry
            ! Try to do the step using current stepsize
            do ndo = 0, MAX_N_ITER
                tt = t

                ! For each value of NS, do a modified-midpoint integration with 2N substeps
                do ns = 1, 8
                    h = dt * C1_2 / real(ns, kind=wp)
                    h2(ns) = C1_4 / (ns * real(ns, kind=wp))
                    hx2 = h * TWO
                    
                    ! Initialize yaux for modified midpoint method
                    ! 1D variables
                    do i = 1, EXTRA2, 2
                        yaux(i) = y0(i) + h * y0(i+1)        ! Position update
                        yaux(i+1) = y0(i+1) + h * der0(i+1)  ! Velocity update
                    end do
                    
                    ! NDIM variables
                    do i = 1, n_x
                        auxi = EXTRA2 + (i-1)*NDIM2 + 1
                        do dim_idx = 0, NDIM-1
                            ! Position update
                            yaux(auxi + dim_idx) = y0(auxi + dim_idx) + h * y0(auxi + NDIM + dim_idx)
                            ! Velocity update  
                            yaux(auxi + NDIM + dim_idx) = y0(auxi + NDIM + dim_idx) + h * der0(auxi + NDIM + dim_idx)
                        end do
                    end do
                    
                    tt = tt + h
                    deraux(:sizey) = dydt(tt, yaux(:sizey))
                    
                    ! First step of modified midpoint
                    ! 1D variables
                    do i = 1, EXTRA2, 2
                        yend(i) = y0(i) + hx2 * yaux(i+1)        ! Position
                        yend(i+1) = y0(i+1) + hx2 * deraux(i+1)  ! Velocity
                    end do
                    
                    ! NDIM variables
                    do i = 1, n_x
                        auxi = EXTRA2 + (i-1)*NDIM2 + 1
                        do dim_idx = 0, NDIM-1
                            ! Position update
                            yend(auxi + dim_idx) = y0(auxi + dim_idx) + hx2 * yaux(auxi + NDIM + dim_idx)
                            ! Velocity update
                            yend(auxi + NDIM + dim_idx) = y0(auxi + NDIM + dim_idx) + hx2 * deraux(auxi + NDIM + dim_idx)
                        end do
                    end do
                    
                    tt = tt + hx2
                    
                    ! Subsequent steps of modified midpoint
                    do j = 2, ns
                        deraux(:sizey) = dydt(tt, yend(:sizey))
                        
                        ! 1D variables
                        do i = 1, EXTRA2, 2
                            yaux(i) = yaux(i) + hx2 * yend(i+1)        ! Position
                            yaux(i+1) = yaux(i+1) + hx2 * deraux(i+1)  ! Velocity
                        end do
                        
                        ! NDIM variables
                        do i = 1, n_x
                            auxi = EXTRA2 + (i-1)*NDIM2 + 1
                            do dim_idx = 0, NDIM-1
                                ! Position update
                                yaux(auxi + dim_idx) = yaux(auxi + dim_idx) + hx2 * yend(auxi + NDIM + dim_idx)
                                ! Velocity update
                                yaux(auxi + NDIM + dim_idx) = yaux(auxi + NDIM + dim_idx) + hx2 * deraux(auxi + NDIM + dim_idx)
                            end do
                        end do
                        
                        tt = tt + hx2
                        
                        deraux(:sizey) = dydt(tt, yaux(:sizey))
                        
                        ! 1D variables
                        do i = 1, EXTRA2, 2
                            yend(i) = yend(i) + hx2 * yaux(i+1)        ! Position
                            yend(i+1) = yend(i+1) + hx2 * deraux(i+1)  ! Velocity
                        end do
                        
                        ! NDIM variables
                        do i = 1, n_x
                            auxi = EXTRA2 + (i-1)*NDIM2 + 1
                            do dim_idx = 0, NDIM-1
                                ! Position update
                                yend(auxi + dim_idx) = yend(auxi + dim_idx) + hx2 * yaux(auxi + NDIM + dim_idx)
                                ! Velocity update
                                yend(auxi + NDIM + dim_idx) = yend(auxi + NDIM + dim_idx) + hx2 * deraux(auxi + NDIM + dim_idx)
                            end do
                        end do
                        
                        tt = tt + hx2
                    end do
                    
                    deraux(:sizey) = dydt(tt, yend(:sizey))
                    
                    ! Final combination for Richardson extrapolation
                    ! 1D variables
                    do i = 1, EXTRA2, 2
                        dy(i,ns) = C1_2 * (yend(i) + yaux(i) + h * yend(i+1))
                        dy(i+1,ns) = C1_2 * (yend(i+1) + yaux(i+1) + h * deraux(i+1))
                    end do
                    
                    ! NDIM variables
                    do i = 1, n_x
                        auxi = EXTRA2 + (i-1)*NDIM2 + 1
                        do dim_idx = 0, NDIM-1
                            ! Position
                            dy(auxi + dim_idx, ns) = C1_2 * (yend(auxi + dim_idx) + yaux(auxi + dim_idx) + &
                                                    h * yend(auxi + NDIM + dim_idx))
                            ! Velocity
                            dy(auxi + NDIM + dim_idx, ns) = C1_2 * (yend(auxi + NDIM + dim_idx) + &
                                                        yaux(auxi + NDIM + dim_idx) + h * deraux(auxi + NDIM + dim_idx))
                        end do
                    end do
            
                    ! Update the DX and DV arrays used for polynomial extrapolation
                    do j = ns - 1, 1, -1
                        j1 = j + 1
                        tmp0 = ONE / (h2(j) - h2(ns))
                        tmp1 = tmp0 * h2(j1)
                        tmp2 = tmp0 * h2(ns)
                        dy(1:sizey,j) = tmp1 * dy(1:sizey,j1) - tmp2 * dy(1:sizey,j)
                    end do
            
                    ! After several integrations, test the relative error on extrapolated values
                    if (ns > 3) then
                        ! Maximum relative position and velocity errors (last D terms added)
                        errmax = ZERO

                        ! 1D variables
                        do i = 1, EXTRA2
                            errmax = MAX(errmax, abs(dy(i,1)) * arr_scal(i))
                        end do
                        
                        ! NDIM variables
                        do i = 1, n_x
                            auxi = EXTRA2 + (i-1)*NDIM2 + 1
                            k = EXTRA2 + (i-1)*2 + 1  ! Index in arr_scal for this object
                            
                            ! Calculate norms for position and velocity errors
                            tmp1 = ZERO
                            tmp2 = ZERO
                            do dim_idx = 0, NDIM-1
                                tmp1 = tmp1 + dy(auxi + dim_idx, 1)**2
                                tmp2 = tmp2 + dy(auxi + NDIM + dim_idx, 1)**2
                            end do
                            
                            errmax = MAX(errmax, sqrt(tmp1) * arr_scal(k), sqrt(tmp2) * arr_scal(k+1))
                        end do

                        ! If error is smaller than TOL, update position and velocity arrays, and exit
                        if (errmax <= E_TOL) then
                            y0(:sizey) = ZERO

                            do j = 1, ns
                                y0(:sizey) = y0(:sizey) + dy(1:sizey,j)
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
                dt = dt * C1_2
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



        !!!! STANDARD subroutines for Bulirsch_Stoer 2. Apply when arr = /x0, v0, x1, v1,.../

        subroutine bstep_std (sizey, y, dydt, t, htry, hdid, hnext)
            implicit none
            integer(kind=4), intent(in) :: sizey
            real(wp), dimension(sizey), intent(inout) :: y
            procedure(dydt_tem) :: dydt
            real(wp), intent(in) :: t, htry
            real(wp), intent(out) :: hdid, hnext

            real(wp), parameter :: SHRINK = 0.55e0_wp
            real(wp), parameter :: GROW   = 1.3e0_wp

            integer(kind=4) :: n_x
            integer(kind=4) :: auxi
            integer(kind=4) :: idim2
            integer(kind=4) :: i, j, j1, k, ns, ndo
            real(wp) :: tmp0, tmp1, tmp2
            real(wp) :: errmax
            real(wp) :: h, hx2, h2(8)
            real(wp) :: dt, tt

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
                    h = dt * C1_2 / real(ns, kind=wp)
                    h2(ns) = C1_4 / (ns * real(ns, kind=wp))
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
                        dy(i,ns) = C1_2 * (yend(i) + yaux(i) + h * yend(i+1))
                        dy(i+1,ns) = C1_2 * (yend(i+1) + yaux(i+1) + h * deraux(i+1))
                    end do
            
                    ! Update the DX and DV arrays used for polynomial extrapolation
                    do j = ns - 1, 1, -1
                        j1 = j + 1
                        tmp0 = ONE / (h2(j) - h2(ns))
                        tmp1 = tmp0 * h2(j1)
                        tmp2 = tmp0 * h2(ns)
                        dy(1:sizey,j) = tmp1 * dy(1:sizey,j1) - tmp2 * dy(1:sizey,j)
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
                                y0(:sizey) = y0(:sizey) + dy(1:sizey,j)
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
                dt = dt * C1_2

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

        end subroutine bstep_std


        !---------------------------------
        !   Call Bulirsch_Stoer 2
        !---------------------------------
        
        subroutine BStoer2_caller (t, y, dt_adap, dydt, dt, ynew, check_fun)
            implicit none
            real(wp), intent(in) :: t
            real(wp), dimension(:), intent(in) :: y
            real(wp), intent(inout) :: dt_adap  ! This try and also next try
            procedure(dydt_tem) :: dydt
            real(wp), intent(in) :: dt
            real(wp), dimension(size (y)), intent(out) :: ynew
            procedure(function_check_keep_tem), optional :: check_fun
            
            integer(kind=4) :: sizey
            real(wp) :: time, t_end, dt_try, dt_used
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

                call bstep_ptr (sizey, ynew, dydt, time, dt_try, dt_used, dt_adap)

                time = time + dt_used

            end do

        end subroutine BStoer2_caller

    
end module bstoer2