!> Module with time and checkpoint calculation routines
module times
    use constants, only: uno, uno2, tini
    use auxiliar, only: quicksort
    implicit none
    real(kind=8), dimension(:), allocatable :: output_times  ! Vector con tiempos de salida
    real(kind=8), dimension(:), allocatable :: checkpoint_times  ! Vector con checkpoints
    
    contains

        ! Set the output times
        subroutine set_output_times(t0, tf, n_out, dt_out, case_dist, t_out)
            implicit none
            real(kind=8), intent(in) :: t0, tf
            integer(kind=4), intent(inout) :: n_out
            real(kind=8), intent(inout) :: dt_out
            integer(kind=4), intent(in) :: case_dist
            real(kind=8), dimension(:), allocatable, intent(out) :: t_out
            real(kind=8) :: t_aux, t_add
            integer(kind=4) :: i
            real(kind=8) :: npointsr

            if (dt_out > 0.d0) then
                if (dt_out > tf - t0) then
                    write (*,*) "ERROR: dt_out >= (tf - t0)"
                    stop 1
                end if
                n_out = int((tf - t0) / dt_out) + 1
                npointsr = n_out * 1.d0
            else
                npointsr = n_out * 1.d0
                dt_out = (tf - t0) / (npointsr - 1.d0)
            end if
            if (n_out < 2) then
                write (*,*) "ERROR: n_out < 2"
                stop 1
            end if
            allocate (t_out(n_out))  ! Here is allocated output_times
            select case (case_dist)
            case(0) ! Lineal
                t_aux = (tf - t0) / (npointsr - uno)
                do i = 2, n_out - 1
                    t_out(i) = t0 + t_aux * (i-1)
                end do
            case(1) ! Logarítmico
                t_aux = exp(log(tf - t0 + 1.) / npointsr)
                t_add = t_aux
                do i = 2, n_out - 1
                    t_add    = t_add * t_aux
                    t_out(i) = t0 + t_add - 1.
                end do
            case(2) ! Combinado
                !! Teniendo en cuenta que 0 es t0 y n_out es tf, si se pide
                !!  una cantidad impar n_out de puntos, se van a hacer 
                !!  lin = (n_out-1)/2 puntos lineales y log = lin-1 puntos logarítmicos
                !!  Si se pide una cantidad par, lin = log = n_out/2
                ! Primero lineal
                npointsr = uno * (ceiling(npointsr * uno2 - tini) + 1) ! - tini para evitar errores de redondeo
                t_aux = (tf - t0) / (npointsr - uno)
                do i = 2, ceiling(npointsr - tini) - 1 ! - tini para evitar errores de redondeo
                    t_out(i) = t0 + t_aux * (i-1)
                end do
                ! Luego logarítmico
                if (mod(n_out, 2) .ne. 0) npointsr = npointsr - uno
                t_aux = exp(log(tf - t0 + 1.) / npointsr)
                t_add = t_aux
                do i = n_out - floor(npointsr + tini) + 2, n_out - 1 ! + tini para evitar errores de redondeo
                    t_add    = t_add * t_aux
                    t_out(i) = t0 + t_add - 1.
                end do
                ! Y al final ordenamos
                call quicksort(t_out, 2, n_out-1)
            end select
            t_out(1) = t0
            t_out(n_out) = tf
        end subroutine set_output_times

    
        ! Free the time arrays
        subroutine free_times_arrays()
            implicit none

            if (allocated(output_times)) deallocate(output_times)
            if (allocated(checkpoint_times)) deallocate(checkpoint_times)
        end subroutine free_times_arrays

end module times