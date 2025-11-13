!> Module handler for integrators
module integrators
    use shared
    use bstoer
    use bstoer2
    use leapfrog
    use embedded
    implicit none
    private
    public :: init_integrator, free_integrator, integrate

    ! Pointer to leapfrog used
    procedure (integrator_caller), pointer :: integrate => null ()

    abstract interface

        subroutine integrator_caller (t, y, dt_adap, dydt, dt, ynew, check_fun)
            import :: dydt_tem
            import :: function_check_keep_tem
            implicit none
            real(kind=8), intent(in) :: t
            real(kind=8), dimension(:), intent(in) :: y
            real(kind=8), intent(inout) :: dt_adap  ! This is each sub-step
            procedure(dydt_tem) :: dydt
            real(kind=8), intent(in) :: dt ! This is full step
            real(kind=8), dimension(size (y)), intent(out) :: ynew
            procedure(function_check_keep_tem), optional :: check_fun
        end subroutine integrator_caller

    end interface
    
    
    contains

            subroutine init_integrator (integrator, sizey, ndimensions, n_extra, min_dt, err_tol, learning_beta, fix_dt)
                implicit none
                integer(kind=4), intent(in) :: integrator, sizey, ndimensions, n_extra
                real(kind=8), intent(in) :: min_dt, err_tol, learning_beta
                logical, intent(in), optional :: fix_dt
                logical :: is_fix_dt = .False.

                if (present(fix_dt)) is_fix_dt = fix_dt
                
                ! Set the integrator
                if (integrator == -3) then
                    call init_BS2(sizey)
                    integrate => BStoer2_caller
                else if (integrator == -2) then
                    call init_leapfrog(sizey,0)
                    if (is_fix_dt) then
                        integrate => leapfrog_fixed_caller
                    else
                        integrate => leapfrog_caller
                    end if
                else if (integrator == -1) then
                    call init_leapfrog(sizey,1)
                    if (is_fix_dt) then
                        integrate => leapfrog_fixed_caller
                    else
                        integrate => leapfrog_caller
                    end if
                else if (integrator == 0) then
                    call init_BS(sizey)
                    integrate => BStoer_caller
                else if (integrator >= 1) then
                    call init_embedded(sizey,integrator)
                    if (is_fix_dt) then
                        integrate => embedded_fixed_caller
                    else
                        integrate => embedded_caller
                    end if
                end if

                ! Allocate der, used in callers
                allocate(der(sizey))

                ! NDimensions
                NDIM = ndimensions
                NDIM2 = NDIM * 2

                ! Extra parameters
                EXTRA = n_extra
                EXTRA2 = EXTRA * 2

                ! Minimum dt
                DT_MIN = min_dt

                ! Total error
                E_TOL = err_tol

                ! Learning rate
                BETA = learning_beta
            end subroutine init_integrator

            subroutine free_integrator()
                if (allocated(der)) deallocate(der)
                call free_BS2()
                call free_leapfrog()
                call free_BS()
                call free_embedded()
            end subroutine free_integrator
    
end module integrators