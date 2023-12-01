module run
    use const, only: cero, inf
    
    implicit none
    integer(kind=4) :: n_iter = 0                                              ! N_iterations
    integer(kind=4) :: n_points = 0                                            ! N_outputs
    real(kind=8) :: t = cero, t0 = inf, tf = cero                              ! Times
    real(kind=8) :: dt = cero, dt_adap = cero, dt_min = cero, dt_out = cero    ! dTimes
    real(kind=8) :: t_add = cero, logt = cero, fixt = cero                     ! Times [output]
    real(kind=8), dimension(:), allocatable :: t_out, omega_out                ! Times [output]
    logical :: logsp = .False.                                                 ! Logarithmic spacing

    contains

        subroutine get_t_outs(t0, tf, n_points, dt_out, logsp, t_out, omega_out, file_tout)
            implicit none
            real(kind=8), intent(in) :: t0, tf
            integer(kind=4), intent(inout) :: n_points
            real(kind=8), intent(inout) :: dt_out
            logical, intent(in) :: logsp
            real(kind=8), dimension(:), allocatable, intent(out) :: t_out, omega_out
            character(LEN=*), optional :: file_tout
            real(kind=8) :: t_aux, t_add
            integer(kind=4) :: i, io
            real(kind=8) :: npointsr

            if (t0 >= tf) then
                write(*,*) "ERROR: t0 >= tf"
                stop 1
            end if
            if (present(file_tout)) then
                n_points = 2
                open(unit=10, file=file_tout, status="old", action="read")
                do
                    read(10, *, iostat=io) t_aux
                    if ((io /= 0) .or. (t_aux > tf)) exit
                    n_points = n_points + 1
                end do
                close(10)
            else 
                if (dt_out > 0.d0) then
                    if (dt_out > tf - t0) then
                        write(*,*) "ERROR: dt_out >= (tf - t0)"
                        stop 1
                    end if
                    n_points = int((tf - t0) / dt_out) + 1
                    npointsr = n_points * 1.d0
                else
                    npointsr = n_points * 1.d0
                    dt_out = (tf - t0) / (npointsr - 1.d0)
                end if
            end if
            if (n_points < 2) then
                write(*,*) "ERROR: n_points < 2"
                stop 1
            end if
            allocate (t_out(n_points))
            t_out(1) = t0
            if (.not. present(file_tout)) then
                if (logsp .eqv. .True.) then
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
            else
                allocate (omega_out(n_points))
                open(unit=10, file=file_tout, status="old", action="read")
                do i = 2, n_points - 1 
                    read(10, *) t_out(i), omega_out(i)
                end do
                close(10)
                omega_out(n_points) = omega_out(n_points - 1)
            end if
            t_out(n_points) = tf
        end subroutine get_t_outs

end module run