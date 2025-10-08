!> Module with time and checkpoint calculation routines
module times
    use constants, only: cero, uno, uno2, tini
    use auxiliary, only: quicksort, merge_sort_and_unique
    implicit none
    
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

            if (dt_out > cero) then
                if (dt_out > tf - t0) then
                    write(*,*) "ERROR: dt_out >= (tf - t0)"
                    stop 1
                end if
                n_out = int((tf - t0) / dt_out) + 1
                npointsr = dble(n_out)
            else
                npointsr = dble(n_out)
                dt_out = (tf - t0) / (npointsr - uno)
            end if
            if (n_out < 2) then
                write(*,*) "ERROR: n_out < 2"
                stop 1
            end if
            if (allocated(t_out)) deallocate(t_out)  ! Deallocate if needed
            allocate (t_out(n_out))  ! Here is allocated output_times
            select case (case_dist)
            case(0) ! Lineal
                t_aux = (tf - t0) / (npointsr - uno)
                do i = 2, n_out - 1
                    t_out(i) = t0 + t_aux * (i - 1)
                end do
            case(1) ! Logarítmico
                t_aux = exp(log(tf - t0 + uno) / npointsr)
                t_add = t_aux
                do i = 2, n_out - 1
                    t_add = t_add * t_aux
                    t_out(i) = t0 + t_add - uno
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
                    t_out(i) = t0 + t_aux * (i - 1)
                end do
                ! Luego logarítmico
                if (mod(n_out, 2) .ne. 0) npointsr = npointsr - uno
                t_aux = exp(log(tf - t0 + uno) / npointsr)
                t_add = t_aux
                do i = n_out - floor(npointsr + tini) + 2, n_out - 1 ! + tini para evitar errores de redondeo
                    t_add = t_add * t_aux
                    t_out(i) = t0 + t_add - uno
                end do
                ! Y al final ordenamos
                call quicksort(t_out, 2, n_out-1)
            end select
            t_out(1) = t0
            t_out(n_out) = tf
        end subroutine set_output_times

        ! Expand the checkpoints
        subroutine expand_checkpoints(extra_checkp_times, checkp_times, checkp_is_a, checkp_is_b, new_checkp_num)
            implicit none
            ! --- Arguments ---
            real(kind=8), intent(in)  :: extra_checkp_times(:)
            real(kind=8), allocatable, intent(inout) :: checkp_times(:)
            logical, allocatable, intent(inout) :: checkp_is_a(:)
            logical, allocatable, intent(inout) :: checkp_is_b(:)
            integer(kind=4), intent(inout) :: new_checkp_num

            ! --- Locals ---
            real(kind=8), allocatable :: merged_times(:)
            logical, allocatable :: ina(:), inb(:)
            logical, allocatable :: merged_a(:), merged_b(:)
            integer(kind=4) :: kfin, i, i_a
            integer(kind=4) :: n_old, n_new

            n_new = size(extra_checkp_times)
            if (allocated(checkp_times)) then
                n_old = size(checkp_times)
            else
                write(*,*) "ERROR: Old checkpoint times not initialized yet."
                stop 1
            end if

            ! --- Case 1: No old checkpoints yet ---
            if (n_old == 0) then
                write(*,*) "ERROR: Old checkpoint times has 0 values."
                stop 1
            end if

            ! --- Case 2: Merge with existing ---
            call merge_sort_and_unique(checkp_times, extra_checkp_times, ina, inb, merged_times, kfin)

            allocate(merged_a(kfin), merged_b(kfin))
            merged_a = .False.
            merged_b = .False.

            i_a = 0

            do i = 1, kfin
                if (ina(i)) then
                    i_a = i_a + 1
                    merged_a(i) = merged_a(i) .or. checkp_is_a(i_a)
                    merged_b(i) = merged_b(i) .or. checkp_is_b(i_a)
                end if
            end do

            ! --- Replace old arrays safely ---
            call move_alloc(merged_times, checkp_times)
            call move_alloc(merged_a, checkp_is_a)
            call move_alloc(merged_b, checkp_is_b)

            new_checkp_num = kfin
        end subroutine expand_checkpoints

end module times