!> Multi-rate KDK Leapfrog integrator: gravity kicks with dt_large; the drift (with the contact
!> forces resolved by sub-steps for the bodies that may touch) is delegated to drift_fast.
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
    !  Core: one large KDK step
    ! =========================================================================
    subroutine leapfrog_KDK_multirate(sizey, y, dydt_slow, drift_fast, t, dt_large, dt_small, deri, ynew)
        implicit none
        integer(kind=4), intent(in) :: sizey
        real(wp), dimension(sizey), intent(in) :: y
        procedure(dydt_tem) :: dydt_slow
        procedure(drift_fast_tem) :: drift_fast
        real(wp), intent(in) :: t, dt_large, dt_small
        real(wp), dimension(sizey), intent(in) :: deri
        real(wp), dimension(sizey), intent(out) :: ynew

        integer(kind=4) :: i
        real(wp) :: dt_half_large

        dt_half_large = C1_2 * dt_large

        ynew = y

        ! ── 1. Half slow kick ─────
        ynew = ynew + dt_half_large * deri
        ! Only velocity slots of der (accelerations) are non-zero here;
        ! position slots of ynew are unchanged.

        ! ── 2. Drift, with the contact forces resolved by sub-steps ──────
        ! 1D extra variable pairs (pos at odd, vel at even index)
        do i = 1, EXTRA2, 2
            ynew(i) = modulo(ynew(i) + dt_large * ynew(i + 1), DOSPI)
        end do
        ! Bodies: only those that may touch take the sub-steps (see collisions_drift)
        call drift_fast(t, ynew, dt_large, dt_small)

        ! ── 3. Half slow kick ────
        ! Positions are now at t + dt_large.
        yaux(:sizey) = dydt_slow(t + dt_large, ynew)
        ynew = ynew + dt_half_large * yaux(:sizey)

    end subroutine leapfrog_KDK_multirate


    ! =========================================================================
    !  Caller: advance y from t to t+dt
    ! =========================================================================
    subroutine leapfrog_multirate_caller(t, y, dt_large, dydt_slow, dt_small, drift_fast, dt, ynew, check_fun)
        implicit none
        real(wp), intent(in) :: t
        real(wp), dimension(:), intent(in) :: y
        real(wp), intent(inout) :: dt_large
        procedure(dydt_tem) :: dydt_slow
        real(wp), intent(in) :: dt_small
        procedure(drift_fast_tem) :: drift_fast
        real(wp), intent(in) :: dt
        real(wp), dimension(size(y)), intent(out) :: ynew
        procedure(function_check_keep_tem), optional :: check_fun

        integer(kind=4) :: sizey
        real(wp) :: time, t_end, dt_step
        logical :: keep
        logical :: has_check

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

            ! This is a fixed-step scheme: dt_large is the nominal step and is NOT shrunk here,
            ! or the short step that closes an interval would be inherited by every later call.
            ! The last step absorbs the round-off remainder instead of taking a ~0 extra step.
            dt_step = min(dt_large, t_end - time)
            if ((t_end - time) - dt_step < 1.e-10_wp * dt_large) dt_step = t_end - time
            DT_MIN_NOW = min(DT_MIN, dt_step)

            der(:sizey) = dydt_slow(time, ycaller(:sizey))

            call leapfrog_KDK_multirate(sizey, ycaller(:sizey), dydt_slow, drift_fast, time, dt_step, dt_small, der(:sizey), ynew)

            if (dt_step >= t_end - time) then
                time = t_end
            else
                time = time + dt_step
            end if
        end do

    end subroutine leapfrog_multirate_caller

end module leapfrog_multirate
