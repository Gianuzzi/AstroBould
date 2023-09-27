module run
    
    implicit none
    integer(kind=4)                         :: n_iter                         ! N_iterations
    integer(kind=4)                         :: n_points                       ! N_outputs
    real(kind=8)                            :: t, t0, tf                      ! Times
    real(kind=8)                            :: dt, dt_adap, dt_min, dt_out    ! dTimes
    real(kind=8)                            :: t_add, logt, fixt              ! Times [output]
    real(kind=8), dimension(:), allocatable :: t_out                          ! Times [output]
    logical, parameter                      :: TRUE = .true., FALSE = .false. ! True | False
    logical                                 :: logsp                          ! Logarithmic spacing
    real(kind=8)                            :: start_time, final_time         ! Excecution start and elapsed time
    real(kind=8), parameter                 :: MIN_VAL = 1.d-15                ! Avoid extremely low values
    ! real(kind=8), parameter                 :: eps   = 1.0d-11                ! Precisión

    contains

        subroutine get_t_outs(t0, tf, n_points, dt_out, logsp, t_out)
            implicit none
            real(kind=8), intent(in)                 :: t0, tf
            real(kind=8), intent(inout)              :: dt_out
            real(kind=8), dimension(:), allocatable  :: t_out
            integer(kind=4), intent(inout)           :: n_points
            logical                                  :: logsp
            real(kind=8)                             :: t_aux, t_add
            integer(kind=4)                          :: i
            real(kind=8)                             :: npointsr

            if (dt_out > 0.d0) then
                if (dt_out > tf - t0) then
                    write(*,*) "dt_out should be less than tf - t0."
                    stop 1
                end if
                n_points = int((tf - t0) / dt_out) + 1
                npointsr = n_points * 1.d0
            else
                npointsr = n_points * 1.d0
                dt_out = (tf - t0) / (npointsr - 1.d0)
            end if
            if (n_points < 2) then
                write(*,*) "n_points should be greater than 2."
                stop 1
            end if
            allocate (t_out(n_points))
            t_out(1) = t0
            if (logsp .eqv. .TRUE.) then
                t_aux = exp (log (tf - t0 + 1.) / npointsr)
                t_add = t_aux
                do i = 2, n_points - 1
                    t_add    = t_add * t_aux
                    t_out(i) = t0 + t_add - 1.
                end do
            else
                t_aux = (tf - t0) / (npointsr - 1.d0)
                do i = 1, n_points - 1
                    t_out(i + 1) = t0 + t_aux * i
                end do
            end if
            t_out(n_points) = tf
        end subroutine get_t_outs

end module run