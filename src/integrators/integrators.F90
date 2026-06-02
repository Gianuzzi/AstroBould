!> Module handler for integrators
module integrators
    use shared
    use bstoer
    use bstoer2
    use leapfrog
    use embedded
    use runge_kutta
    use leapfrog_multirate
    
    implicit none
    private
    public :: init_integrator, free_integrator, integrate, integrate_substeps, check_integrator, fallback_integrator

    ! Pointer to leapfrog used
    procedure(integrator_caller), pointer :: integrate => null()
    procedure(integrator_substep_caller), pointer :: integrate_substeps => null()


    abstract interface

        subroutine integrator_caller(t, y, dt_adap, dydt, dt, ynew, check_fun)
            import :: wp
            import :: dydt_tem
            import :: function_check_keep_tem
            implicit none
            real(wp), intent(in) :: t
            real(wp), dimension(:), intent(in) :: y
            real(wp), intent(inout) :: dt_adap  ! This is each sub-step
            procedure(dydt_tem) :: dydt
            real(wp), intent(in) :: dt ! This is full step
            real(wp), dimension(size(y)), intent(out) :: ynew
            procedure(function_check_keep_tem), optional :: check_fun
        end subroutine integrator_caller

        subroutine integrator_substep_caller(t, y, dt_large, dydt_slow, dt_small, dydt_fast, dt, ynew, check_fun)
            import :: wp
            import :: dydt_tem
            import :: function_check_keep_tem
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
        end subroutine integrator_substep_caller

    end interface

contains

    subroutine init_integrator(integrator, sizey, n_dimensions, n_extra, min_dt, err_tol, learning_beta, fix_dt, standard)
        implicit none
        integer(kind=4), intent(in) :: integrator, sizey, n_dimensions, n_extra
        real(wp), intent(in) :: min_dt, err_tol, learning_beta
        logical, intent(in) :: fix_dt
        logical, intent(in), optional :: standard
        logical :: is_std = .False.

        ! Size of y
        SIZEY0 = sizey

        ! Fixed dt
        FIXED_DT = fix_dt

        ! Allocate der, used in callers
        allocate (der(SIZEY0))

        ! NDimensions
        NDIM = n_dimensions
        NDIM2 = NDIM*2

        ! Extra parameters
        EXTRA = n_extra
        EXTRA2 = EXTRA*2

        ! Minimum dt
        DT_MIN = min_dt
        DT_MIN_NOW = DT_MIN

        ! Total error
        E_TOL = err_tol

        ! Learning rate
        BETA = learning_beta

        ! Check if std
        if (present(standard)) is_std = standard

        ! Set the integrator
        if (integrator == -30) then
            call init_leapfrog_multirate(SIZEY0, 0)
            integrate_substeps => leapfrog_multirate_caller
            print*, "ALL SET"
        
        else if (integrator > -30 .and. integrator < -3) then
            call init_runge_kutta(SIZEY0, abs(integrator + 3))
            if (FIXED_DT) then
                integrate => runge_kutta_fixed_caller
            else
                integrate => runge_kutta_caller
            end if

        else if (integrator == -3) then
            call init_BS2(SIZEY0)
            integrate => BStoer2_caller

        else if (integrator == -2) then
            call init_leapfrog(SIZEY0, 0, is_std)
            if (FIXED_DT) then
                integrate => leapfrog_fixed_caller
            else
                integrate => leapfrog_caller
            end if

        else if (integrator == -1) then
            call init_leapfrog(SIZEY0, 1, is_std)
            if (FIXED_DT) then
                integrate => leapfrog_fixed_caller
            else
                integrate => leapfrog_caller
            end if

        else if (integrator == 0) then
            call init_BS(SIZEY0)
            integrate => BStoer_caller

        else if (integrator >= 1) then
            call init_embedded(SIZEY0, integrator)
            if (FIXED_DT) then
                integrate => embedded_fixed_caller
            else
                integrate => embedded_caller
            end if

        end if
    end subroutine init_integrator

    subroutine free_integrator()
        implicit none
        if (allocated(der)) deallocate (der)
        call free_BS2()
        call free_leapfrog()
        call free_BS()
        call free_embedded()
        call free_runge_kutta()
    end subroutine free_integrator

    subroutine fallback_integrator()
        implicit none
        integer(kind=4) :: which

        call free_integrator()

        if (FIXED_DT) then
            which = 12 ! Fallback to Dopri8
        else
            which = 6 ! Fallback to RK6
        end if

        call init_integrator(which, SIZEY0, NDIM, EXTRA, DT_MIN, E_TOL, BETA, FIXED_DT)

    end subroutine fallback_integrator

    subroutine check_integrator(aux_logical)
        implicit none
        logical, intent(out) :: aux_logical
        ! we only have to check BStoer
        if (reached_underflow) then
            ! print*, "Underflow detected in BStoer. Switching to fallback integrator."
            call fallback_integrator()
            reached_underflow = .False. ! Reset the flag for future checks
            aux_logical = .True.
        else
            aux_logical = .False.
        end if
    end subroutine check_integrator

end module integrators
