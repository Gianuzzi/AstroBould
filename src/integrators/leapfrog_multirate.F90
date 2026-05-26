!> Multi-rate KDK Leapfrog integrator with collision sub-stepping.
module leapfrog_multirate
    use shared

    implicit none
    private
    public :: init_leapfrog_multirate, free_leapfrog_multirate, leapfrog_multirate_caller

    ! Workspace arrays
    real(wp), allocatable :: ycaller(:) ! caller
    real(wp), allocatable :: yaux(:) ! solver

contains

    !!!! HANDLER

    subroutine init_leapfrog_multirate(sizey, dummy)
        implicit none
        integer(kind=4), intent(in) :: sizey
        integer(kind=4), intent(in) :: dummy
        allocate(ycaller(sizey))
        allocate(yaux(sizey))
    end subroutine init_leapfrog_multirate

    subroutine free_leapfrog_multirate()
        implicit none
        if (allocated(ycaller)) deallocate(ycaller)
        if (allocated(yaux)) deallocate(yaux)
    end subroutine free_leapfrog_multirate

    ! =========================================================================
    !  Drift: pos += dt * vel, in-place on y
    ! =========================================================================
    subroutine drift_step(y, dt)
        implicit none
        real(wp), dimension(:), intent(inout) :: y
        real(wp), intent(in) :: dt

        integer(kind=4) :: i, j

        ! ── 1D extra variable pairs (pos at odd, vel at even index) ──
        do i = 1, EXTRA2, 2
            y(i) = modulo(y(i) + dt * y(i + 1), DOSPI)
        end do

        ! ── NDIM body variables ───
        do i = EXTRA2 + 1, size(y), NDIM2
            do j = 0, NDIM - 1
                y(i + j) = y(i + j) + dt * y(i + NDIM + j)
            end do
        end do

    end subroutine drift_step


    ! =========================================================================
    !  Core: one large KDK step
    ! =========================================================================
    subroutine leapfrog_KDK_multirate(sizey, y, dydt_slow, dydt_fast, t, dt_large, dt_small, deri, ynew)
        implicit none
        integer(kind=4), intent(in) :: sizey
        real(wp), dimension(sizey), intent(in) :: y
        procedure(dydt_tem) :: dydt_slow, dydt_fast
        real(wp), intent(in) :: t, dt_large, dt_small
        real(wp), dimension(sizey), intent(in) :: deri
        real(wp), dimension(sizey), intent(out) :: ynew

        integer(kind=4) :: N_sub, k
        real(wp) :: dt_small_eff, dt_half_large, t_sub

        ! ── N_sub: nearest integer, minimum 1 ────
        N_sub = max(1, nint(dt_large / dt_small))
        dt_small_eff = dt_large / real(N_sub, wp)
        dt_half_large = C1_2 * dt_large

        ynew = y

        ! ── 1. Half slow kick ─────
        ynew = ynew + dt_half_large * deri
        ! Only velocity slots of der (accelerations) are non-zero here;
        ! position slots of ynew are unchanged.

        ! ── 2. Collision sub-steps ──────
        t_sub = t

        do k = 1, N_sub

            ! 2a. Half fast kick — velocity update only
            yaux(:sizey) = dydt_fast(t_sub, ynew)
            ynew = ynew + C1_2 * dt_small_eff * yaux(:sizey)

            ! 2b. Drift — use velocity slots of ynew directly
            call drift_step(ynew, dt_small_eff)

            ! 2c. Half fast kick — velocity update only
            yaux(:sizey) = dydt_fast(t_sub + dt_small_eff * C1_2, ynew)  ! * C1_2 is to have Verlet see pos at ~ t + dt_small_eff/2
            ynew = ynew + C1_2 * dt_small_eff * yaux(:sizey)

            t_sub = t_sub + dt_small_eff

        end do

        ! ── 3. Half slow kick ────
        ! Positions are now at t + dt_large.
        yaux(:sizey) = dydt_slow(t + dt_large, ynew)
        ynew = ynew + dt_half_large * yaux(:sizey)

    end subroutine leapfrog_KDK_multirate


    ! =========================================================================
    !  Caller: advance y from t to t+dt
    ! =========================================================================
    subroutine leapfrog_multirate_caller(t, y, dt_large, dydt_slow, dt_small, dydt_fast, dt, ynew, check_fun)
        implicit none
        real(wp), intent(in) :: t
        real(wp), dimension(:), intent(in) :: y
        real(wp), intent(inout) :: dt_large
        procedure(dydt_tem) :: dydt_slow
        real(wp), intent(in) :: dt_small
        procedure(dydt_tem) :: dydt_fast
        real(wp), intent(in) :: dt
        real(wp), dimension(size(y)), intent(out) :: ynew
        procedure(function_check_keep_tem), optional :: check_fun

        integer(kind=4) :: sizey
        real(wp) :: time, t_end
        logical :: keep = .True.
        logical :: has_check = .False.

        sizey = size(y)
        has_check = present(check_fun)

        ynew = y
        time = t
        t_end = time + dt

        do while (time < t_end)

            if (has_check) then ! If Check Continue function present
                keep = check_fun(ynew)
                if (.not. keep) then ! If Hard Exit is True
                    dt_large = time - t ! Replace dt_large with actual dt used
                    return ! Exit subroutine
                end if
            end if

            ycaller(:sizey) = ynew
            dt_large = min(dt_large, t_end - time)
            DT_MIN_NOW = min(DT_MIN, dt_large)

            der(:sizey) = dydt_slow(time, ycaller(:sizey))

            call leapfrog_KDK_multirate(sizey, ycaller(:sizey), dydt_slow, dydt_fast, time, dt_large, dt_small, der(:sizey), ynew)

            time = time + dt_large
        end do

    end subroutine leapfrog_multirate_caller

end module leapfrog_multirate
