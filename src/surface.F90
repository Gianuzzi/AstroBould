!> Module with surface section calculation routines.
module surface
    use constants, only: wp, cero, dos, myepsilon, unit_dist, unit_vel, G

    implicit none
    private
    public :: section_st, init_section, crossed_section, get_jacobi_constant

    type section_st
        logical :: active = .False.            ! Flag is the surface is active
        integer(kind=4) :: idx_surface = -1    ! index that defines the section (e.g. yP = 8)
        real(wp) :: valor = cero               ! value of the section (usually 0.0d0)
        integer(kind=4) :: direction = 0       ! +1: upward crossing, -1: downward crossing, 0: any direction
        integer(kind=4) :: idx_condition = -1  ! optional extra condition index (e.g. vyP = 10)
        character(2) :: name = "NO"
        real(wp) :: condition_min = cero       ! extra condition: y(idx_condition) > condition_min
        logical :: use_condition = .False.     ! Flag if uses extra condition
        logical :: surface_is_r = .False.      ! Flag is the surface is r
        logical :: condition_is_r = .False.    ! Flag is the extra condition is r
        logical :: surface_is_vr = .False.     ! Flag is the surface is vr
        logical :: condition_is_vr = .False.   ! Flag is the extra condition is vr
        character(2) :: condition_name = "NO"
    end type section_st

    integer(kind=4), parameter :: offset = 6
    integer(kind=4), parameter :: idx_x = offset + 1
    integer(kind=4), parameter :: idx_y = offset + 2
    integer(kind=4), parameter :: idx_vx = offset + 3
    integer(kind=4), parameter :: idx_vy = offset + 4
    character(len=2), dimension(6), parameter :: names = (/"x ", "y ", "vx", "vy", "r ", "vr"/) 

contains

    pure subroutine init_section(sec, idx_surface, valor, direction, idx_condition, condition_min)
        implicit none
        type(section_st), intent(inout) :: sec
        integer(kind=4), intent(in) :: idx_surface
        real(wp), intent(in) :: valor
        integer(kind=4), intent(in) :: direction
        integer(kind=4), intent(in), optional :: idx_condition
        real(wp), intent(in), optional :: condition_min

        if ((idx_surface < 1) .or. (idx_surface > 6)) then
            sec%active = .False.
            return  ! No surface
        end if
        
        ! x, y, vx, vy, r, vr
        ! 1, 2,  3,  4, 5,  6

        ! Check if r or vr
        sec%surface_is_r = idx_surface == 5
        sec%surface_is_vr = idx_surface == 6

        ! Add offset
        sec%idx_surface = idx_surface + offset

        ! Set with units
        if ((idx_surface == 1) .or. (idx_surface == 2) .or. (idx_surface == 5)) then
            sec%valor = valor * unit_dist
        else
            sec%valor = valor * unit_vel
        end if
        sec%direction = direction

        ! Set name
        sec%name = names(idx_surface)

        ! Extra condition
        sec%use_condition = .False.
        if (present(idx_condition)) then
            
            if ((idx_condition >= 1) .and. (idx_condition <= 6)) then
                
                ! Check if r or vr
                sec%condition_is_r = idx_condition == 5
                sec%condition_is_vr = idx_condition == 6

                ! Add offset
                sec%idx_condition = idx_condition + offset

                ! Set with units
                if ((idx_condition == 1) .or. (idx_condition == 2) .or. (idx_condition == 5)) then
                    sec%condition_min = condition_min * unit_dist
                else
                    sec%condition_min = condition_min * unit_vel
                end if

                ! Set name
                sec%condition_name = names(idx_condition)

                sec%use_condition = .True.
            
            end if

        end if

        sec%active = .True.

    end subroutine

    pure function compute_r(y) result(r)
        implicit none
        real(wp), intent(in) :: y(:)
        real(wp) :: r

        r = sqrt(y(idx_x)**2 + y(idx_y)**2)

    end function compute_r

    pure function compute_vr(y) result(vr)
        implicit none
        real(wp), intent(in) :: y(:)
        real(wp) :: r
        real(wp) :: vr

        r = compute_r(y)

        if (r < myepsilon) then
            vr = cero
        else
            vr = (y(idx_x)*y(idx_vx) + y(idx_y)*y(idx_vy)) / r
        end if

    end function compute_vr

    pure subroutine crossed_section(sec, y_old, y_new, f_old, f_new, has_crossed)
        implicit none
        type(section_st), intent(in) :: sec
        real(wp), intent(in) :: y_old(:), y_new(:)
        real(wp), intent(inout) :: f_old, f_new
        logical, intent(inout) :: has_crossed
        real(wp) :: val_old, val_new
        logical :: direction_ok, condition_ok

        ! Default
        val_old = cero
        val_new = cero
        has_crossed = .False.

        ! No surface
        if (.not. sec%active) return        

         ! ---- Surface evaluation ----
        if (sec%surface_is_r) then
            val_old = compute_r(y_old)
            val_new = compute_r(y_new)
        else if (sec%surface_is_vr) then
            val_old = compute_vr(y_old)
            val_new = compute_vr(y_new)
        else
            val_old = y_old(sec%idx_surface)
            val_new = y_new(sec%idx_surface)
        end if

        ! Get values
        f_old = val_old - sec%valor
        f_new = val_new - sec%valor

        ! Must change sign
        if (f_old * f_new > cero) return

        ! ---- Direction control ----
        select case (sec%direction)
            case (1)
                direction_ok = (f_old <= cero) .and. (f_new > cero)
            case (-1)
                direction_ok = (f_old >= cero) .and. (f_new < cero)
            case default
                direction_ok = .True.
        end select

        if (.not. direction_ok) return

        ! Optional extra condition (evaluated at new step)
        if (sec%use_condition) then
            if (sec%condition_is_r) then
                condition_ok = compute_r(y_new) > sec%condition_min
            else if (sec%condition_is_vr) then
                condition_ok = compute_vr(y_new) > sec%condition_min
            else
                condition_ok = y_new(sec%idx_condition) > sec%condition_min
            end if

            if (.not. condition_ok) return
        end if

        has_crossed = .True.

    end subroutine crossed_section


    !----------------------------------------------------------
    ! Gravitational potential
    !----------------------------------------------------------
    pure function get_gravitational_potential(x, y, masses, pos) result(U)
        implicit none
        real(wp), intent(in) :: x, y
        real(wp), intent(in) :: masses(:)   ! masses(N)
        real(wp), intent(in) :: pos(:, :)   ! pos(N, 2)
        real(wp) :: U, r
        integer(kind=4) :: i, N

        U = cero
        N = size(masses)

        do i = 1, N
            r = sqrt( (x - pos(i,1))**2 &
                    + (y - pos(i,2))**2 )

            U = U + G*masses(i) / r
        end do

    end function get_gravitational_potential


    !----------------------------------------------------------
    ! Jacobi constant
    !----------------------------------------------------------
    pure function get_jacobi_constant(x, y, vx, vy, masses, pos, lam) result(J)
        implicit none

        real(wp), intent(in) :: x, y
        real(wp), intent(in) :: vx, vy
        real(wp), intent(in) :: lam
        real(wp), intent(in) :: masses(:)
        real(wp), intent(in) :: pos(:, :)
        real(wp) :: J
        real(wp) :: U

        U = get_gravitational_potential(x, y, masses, pos)

        J = lam**2 * (x**2 + y**2) &
            + dos * U &
            - (vx**2 + vy**2)

    end function get_jacobi_constant

end module surface