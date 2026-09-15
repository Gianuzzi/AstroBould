!> Module with surface section calculation routines.
module surface
    use constants, only: wp, cero, dos, uno2, myepsilon, infinito, unit_dist, unit_vel, G

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

        ! Set index
        sec%idx_surface = idx_surface

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

                ! Set condition index
                sec%idx_condition = idx_condition

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

        r = sqrt(y(1)**2 + y(2)**2)

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
            vr = (y(1)*y(3) + y(2)*y(4)) / r
        end if

    end function compute_vr

    pure subroutine crossed_section(sec, y_old, y_new, alpha, error, has_crossed)
        implicit none
        type(section_st), intent(in) :: sec
        real(wp), intent(in) :: y_old(:), y_new(:)
        real(wp), intent(inout) :: alpha, error
        logical, intent(inout) :: has_crossed(:)
        real(wp) :: val_old, val_new, f_old, f_new
        logical :: direction_ok, condition_ok
        real(wp) :: y_old_this(4), y_new_this(4)
        integer(kind=4) :: ntotal, idx
        
        ntotal = size(has_crossed)

        ! Default
        alpha = infinito
        error = infinito
        do idx = 1, ntotal
            has_crossed(idx) = .False.
        end do

        ! No surface
        if (.not. sec%active) return

        ! Loop over all particles
        do idx = 1, ntotal

            y_old_this = y_old(4*idx + 3 : 4*idx + 6) ! [TODO : This is a bit ugly, but works]
            y_new_this = y_new(4*idx + 3 : 4*idx + 6) ! [TODO : This is a bit ugly, but works]

            ! ---- Surface evaluation ----
            if (sec%surface_is_r) then
                val_old = compute_r(y_old_this)
                val_new = compute_r(y_new_this)
            else if (sec%surface_is_vr) then
                val_old = compute_vr(y_old_this)
                val_new = compute_vr(y_new_this)
            else
                val_old = y_old_this(sec%idx_surface)
                val_new = y_new_this(sec%idx_surface)
            end if

            ! Get values
            f_old = val_old - sec%valor
            f_new = val_new - sec%valor

            ! Must change sign
            if (f_old * f_new > cero) cycle  ! This not crossing

            ! ---- Direction control ----
            select case (sec%direction)
                case (1)
                    direction_ok = (f_old <= cero) .and. (f_new > cero)
                case (-1)
                    direction_ok = (f_old >= cero) .and. (f_new < cero)
                case default
                    direction_ok = .True.
            end select

            if (.not. direction_ok) cycle  ! This not crossing in the right direction

            ! Optional extra condition (evaluated at new step)
            if (sec%use_condition) then
                if (sec%condition_is_r) then
                    condition_ok = compute_r(y_new_this) > sec%condition_min
                else if (sec%condition_is_vr) then
                    condition_ok = compute_vr(y_new_this) > sec%condition_min
                else
                    condition_ok = y_new_this(sec%idx_condition) > sec%condition_min
                end if

                if (.not. condition_ok) cycle  ! This not crossing because extra condition not satisfied
            end if


            ! Crossed!
            has_crossed(idx) = .True.
            ! abs just in case, but should be positive since they have different signs
            alpha = min(alpha, abs(-f_old / (f_new - f_old)))
            ! Estimate of the error in the crossing point (the smaller of the two values) 
            error = min(error, min(abs(f_old), abs(f_new)))  
        
        end do

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