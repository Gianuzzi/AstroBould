!> Module with main collision routines.
module collisions
    use iso_fortran_env, only: int64
    use omp_lib
    use constants, only: wp, G, cero, uno, uno2, dos, twopi, tini
    use auxiliary, only: rotate2D
    use parameters, only: sim, m_arr, R_arr, get_index

    private
    public :: init_collisions, collisions_brute, collisions_grid, collisions_verlet, verlet_rebuilds

    abstract interface
        subroutine get_xy_rotated_tem(xy, r, dt, mu, omega)
            import :: wp
            implicit none
            real(wp), dimension(2), intent(inout) :: xy
            real(wp), intent(in) :: r, dt, mu, omega
        end subroutine get_xy_rotated_tem
    end interface

    procedure(get_xy_rotated_tem), pointer :: get_xy_rotated => null()

    integer(kind=4), save :: vlist_n
    real(wp), save        :: vrcut
    integer(kind=4), save, allocatable :: vlist(:, :)
    real(wp), save, allocatable        :: vlist_pos(:, :)
    real(wp), save, allocatable        :: vlist_r(:)
    real(wp), save                     :: vlist_time
    logical,  save                     :: vlist_built = .False.
    logical,  save, allocatable        :: list_collided(:)

    integer(int64), save :: verlet_rebuilds = 0_int64

contains

    subroutine init_collisions(N_total, sinodic)
        implicit none
        integer(kind=4), intent(in) :: N_total
        logical,         intent(in) :: sinodic

        if (allocated(list_collided)) deallocate(list_collided)
        allocate(list_collided(N_total))
        list_collided = .False.

        if (sinodic) then
            get_xy_rotated => get_xy_rotated_sinodic
        else
            get_xy_rotated => get_xy_rotated_inertial
        end if

    end subroutine init_collisions

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!! SOFT-SPHERE COLLISION ROUTINES !!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !> Soft-sphere force kernel with ATOMIC writes on all der and list_collided.
    !> Safe to call from any parallel region.
    subroutine soft_sphere_force(y, der, i, j, gamma_n, gamma_t, are_moons)
        implicit none
        real(wp), dimension(:), intent(in)    :: y
        real(wp), dimension(:), intent(inout) :: der
        integer(kind=4),        intent(in)    :: i, j
        real(wp),               intent(in)    :: gamma_n, gamma_t
        logical,                intent(in)    :: are_moons

        integer(kind=4) :: idx, jdx
        real(wp) :: dr_vec(2), dr2, dr, rcoll, overlap
        real(wp) :: dr_ver(2), dt_ver(2)
        real(wp) :: dv_vec(2), dvr, dvt
        real(wp) :: F_n, F_t, F_vec(2)
        real(wp) :: gamma_n_pair, gamma_t_pair
        real(wp) :: mi, mj
        real(wp) :: rand1, rn
        real(wp) :: aux_real

        idx = get_index(i)
        jdx = get_index(j)

        dr_vec = y(jdx:jdx+1) - y(idx:idx+1)
        dr2    = dr_vec(1)*dr_vec(1) + dr_vec(2)*dr_vec(2)
        rcoll  = R_arr(i) + R_arr(j)

        if (dr2 >= rcoll*rcoll) return

        if (dr2 < tini) then
            dv_vec = y(jdx+2:jdx+3) - y(idx+2:idx+3)
            rn = dv_vec(1)*dv_vec(1) + dv_vec(2)*dv_vec(2)
            if (rn > tini) then
                rn     = sqrt(rn)
                dr_ver = dv_vec / rn
            else
                call random_number(rand1)
                rand1  = twopi * rand1
                dr_ver = [cos(rand1), sin(rand1)]
            end if
            dr      = tini
            overlap = rcoll
        else
            dr      = sqrt(dr2)
            overlap = rcoll - dr
            dr_ver  = dr_vec / dr
        end if

        dt_ver = [-dr_ver(2), dr_ver(1)]

        if (are_moons) then
            mi       = m_arr(i)
            mj       = m_arr(j)
            aux_real = dos * sqrt(sim%kappa_col_moon * (mi * mj / (mi + mj)))
            gamma_n_pair = min(sim%gamma_col_moon_n, uno) * aux_real
            gamma_t_pair = min(sim%gamma_col_moon_t, uno) * aux_real
        else
            mi           = uno
            mj           = uno
            gamma_n_pair = gamma_n
            gamma_t_pair = gamma_t
        end if

        dv_vec = y(jdx+2:jdx+3) - y(idx+2:idx+3)
        dvr    = dv_vec(1)*dr_ver(1) + dv_vec(2)*dr_ver(2)
        dvt    = dv_vec(1)*dt_ver(1) + dv_vec(2)*dt_ver(2)

        F_n = sim%kappa_col_part * overlap
        if (gamma_n_pair > cero) F_n = F_n - gamma_n_pair * min(dvr, cero)
        F_n = max(F_n, cero)

        F_t = cero
        if ((gamma_t_pair > cero) .and. (sim%coulomb_mu_col > cero) .and. (F_n > cero)) then
            F_t = -gamma_t_pair * dvt
            F_t = sign(min(abs(F_t), sim%coulomb_mu_col * F_n), F_t)
        end if

        F_vec = F_n * dr_ver + F_t * dt_ver

        ! ── Newton's 3rd law — ATOMIC on all four velocity slots ────────────
        !$OMP ATOMIC UPDATE
        der(jdx+2) = der(jdx+2) + F_vec(1) / mj
        !$OMP ATOMIC UPDATE
        der(jdx+3) = der(jdx+3) + F_vec(2) / mj
        !$OMP ATOMIC UPDATE
        der(idx+2) = der(idx+2) - F_vec(1) / mi
        !$OMP ATOMIC UPDATE
        der(idx+3) = der(idx+3) - F_vec(2) / mi

        !$OMP ATOMIC WRITE
        list_collided(i) = .True.
        !$OMP ATOMIC WRITE
        list_collided(j) = .True.

    end subroutine soft_sphere_force


    !> Brute-force O(N²).
    !> Outer i loop is PARALLEL — atomics inside soft_sphere_force handle conflicts.
    !> DYNAMIC schedule because the triangular inner loop makes load uneven.
    subroutine collisions_brute(y, der, first_index, last_index, gamma_n, gamma_t, are_moons)
        implicit none
        real(wp), dimension(:), intent(in)    :: y
        real(wp), dimension(:), intent(inout) :: der
        integer(kind=4),        intent(in)    :: first_index, last_index
        real(wp),               intent(in)    :: gamma_n, gamma_t
        logical,                intent(in)    :: are_moons

        integer(kind=4) :: i, j

        !$OMP PARALLEL DO                  &
        !$OMP   DEFAULT(SHARED)            &
        !$OMP   PRIVATE(i, j)              &
        !$OMP   SCHEDULE(DYNAMIC, 400)
        do i = first_index, last_index - 1
            do j = i + 1, last_index
                call soft_sphere_force(y, der, i, j, gamma_n, gamma_t, are_moons)
            end do
        end do
        !$OMP END PARALLEL DO

    end subroutine collisions_brute


    !> Cell-list O(N).
    !> The linked-list traversal (head/next) is not thread-safe to build in
    !> parallel, so this routine stays serial. The force kernel itself is called
    !> serially here — if this loop ever becomes the bottleneck, use colouring.
    subroutine collisions_grid(y, der, first_index, last_index, gamma_n, gamma_t, are_moons)
        implicit none
        real(wp), dimension(:), intent(in)    :: y
        real(wp), dimension(:), intent(inout) :: der
        integer(kind=4),        intent(in)    :: first_index, last_index
        real(wp),               intent(in)    :: gamma_n, gamma_t
        logical,                intent(in)    :: are_moons

        real(wp) :: L_cell, x_min, x_max, y_min, y_max
        integer(int64) :: Ncx, Ncy, Nc
        integer(int64) :: cx, cy, ci, cx2, cy2, ci2
        integer(kind=4) :: dcx, dcy
        integer(kind=4) :: i, j, idx
        integer(kind=4), allocatable :: head(:), next(:)

        if (are_moons) then
            L_cell = dos * maxval(R_arr(first_index:last_index))
        else
            L_cell = dos * R_arr(first_index)
        end if
        L_cell = max(L_cell, sim%grid_col_min_cell_size)

        x_min = y(get_index(first_index))
        x_max = x_min
        y_min = y(get_index(first_index)+1)
        y_max = y_min

        do i = first_index, last_index
            idx   = get_index(i)
            x_min = min(x_min, y(idx))
            x_max = max(x_max, y(idx))
            y_min = min(y_min, y(idx+1))
            y_max = max(y_max, y(idx+1))
        end do

        x_min = x_min - L_cell
        x_max = x_max + L_cell
        y_min = y_min - L_cell
        y_max = y_max + L_cell

        Ncx = max(1, ceiling((x_max - x_min) / L_cell))
        Ncy = max(1, ceiling((y_max - y_min) / L_cell))

        if (Ncx > sim%grid_col_max_cells / Ncy) then
            call collisions_brute(y, der, first_index, last_index, gamma_n, gamma_t, are_moons)
            return
        end if

        Nc = Ncx * Ncy
        allocate(head(0:Nc-1), next(first_index:last_index))
        head = -1

        do i = first_index, last_index
            idx = get_index(i)
            cx  = min(int((y(idx)   - x_min) / L_cell, int64), Ncx - 1_int64)
            cy  = min(int((y(idx+1) - y_min) / L_cell, int64), Ncy - 1_int64)
            cx  = max(cx, 0_int64)
            cy  = max(cy, 0_int64)
            ci  = cx + Ncx * cy
            next(i)  = head(ci)
            head(ci) = i
        end do

        do i = first_index, last_index
            idx = get_index(i)
            cx  = min(max(int((y(idx)   - x_min) / L_cell, int64), 0_int64), Ncx - 1_int64)
            cy  = min(max(int((y(idx+1) - y_min) / L_cell, int64), 0_int64), Ncy - 1_int64)

            do dcy = -1, 1
                do dcx = -1, 1
                    cx2 = cx + dcx
                    cy2 = cy + dcy
                    if (cx2 < 0_int64 .or. cx2 >= Ncx) cycle
                    if (cy2 < 0_int64 .or. cy2 >= Ncy) cycle
                    ci2 = cx2 + Ncx * cy2
                    j   = head(ci2)
                    do while (j /= -1)
                        if (j > i) call soft_sphere_force(y, der, i, j, gamma_n, gamma_t, are_moons)
                        j = next(j)
                    end do
                end do
            end do
        end do

        deallocate(head, next)

    end subroutine collisions_grid


    !> Build Verlet neighbour list (serial — modifies global saved state).
    subroutine build_verlet_list(time, y, first_index, last_index, are_moons)
        implicit none
        real(wp),               intent(in) :: time
        real(wp), dimension(:), intent(in) :: y
        integer(kind=4),        intent(in) :: first_index, last_index
        logical,                intent(in) :: are_moons

        integer(kind=4) :: i, j, idx, jdx, npairs, N_particles
        real(wp) :: r_cut, r_cut2, dr_vec(2), dr2
        integer(kind=4), allocatable :: tmp_list(:, :)

        real(wp) :: x_min, x_max, y_min, y_max
        integer(int64) :: Ncx, Ncy, Nc, cx, cy, ci, cx2, cy2, ci2
        integer(kind=4) :: dcx, dcy
        integer(kind=4), allocatable :: head(:), next(:)

        N_particles = last_index - first_index + 1

        if (are_moons) then
            vrcut = dos * maxval(R_arr(first_index:last_index))
        else
            vrcut = dos * R_arr(first_index)
        end if
        vrcut  = max(vrcut, sim%grid_col_min_cell_size)
        r_cut  = vrcut * (uno + sim%verlet_skin_factor)
        r_cut2 = r_cut * r_cut

        allocate(tmp_list(2, N_particles*(N_particles-1)/2))
        npairs = 0

        x_min = y(get_index(first_index))
        x_max = x_min
        y_min = y(get_index(first_index)+1)
        y_max = y_min

        do i = first_index, last_index
            idx   = get_index(i)
            x_min = min(x_min, y(idx))
            x_max = max(x_max, y(idx))
            y_min = min(y_min, y(idx+1))
            y_max = max(y_max, y(idx+1))
        end do

        x_min = x_min - vrcut
        x_max = x_max + vrcut
        y_min = y_min - vrcut
        y_max = y_max + vrcut

        Ncx = max(1, ceiling((x_max - x_min) / vrcut))
        Ncy = max(1, ceiling((y_max - y_min) / vrcut))

        if (Ncx > sim%grid_col_max_cells / Ncy) then
            do i = first_index, last_index - 1
                idx = get_index(i)
                do j = i + 1, last_index
                    jdx    = get_index(j)
                    dr_vec = y(jdx:jdx+1) - y(idx:idx+1)
                    dr2    = dr_vec(1)*dr_vec(1) + dr_vec(2)*dr_vec(2)
                    if (dr2 < r_cut2) then
                        npairs = npairs + 1
                        tmp_list(1, npairs) = i
                        tmp_list(2, npairs) = j
                    end if
                end do
            end do
        else
            Nc = Ncx * Ncy
            allocate(head(0:Nc-1), next(first_index:last_index))
            head = -1

            do i = first_index, last_index
                idx = get_index(i)
                cx  = min(int((y(idx)   - x_min) / vrcut, int64), Ncx - 1_int64)
                cy  = min(int((y(idx+1) - y_min) / vrcut, int64), Ncy - 1_int64)
                cx  = max(cx, 0_int64)
                cy  = max(cy, 0_int64)
                ci  = cx + Ncx * cy
                next(i)  = head(ci)
                head(ci) = i
            end do

            do i = first_index, last_index
                idx = get_index(i)
                cx  = min(max(int((y(idx)   - x_min) / vrcut, int64), 0_int64), Ncx - 1_int64)
                cy  = min(max(int((y(idx+1) - y_min) / vrcut, int64), 0_int64), Ncy - 1_int64)

                do dcy = -1, 1
                    do dcx = -1, 1
                        cx2 = cx + dcx
                        cy2 = cy + dcy
                        if (cx2 < 0_int64 .or. cx2 >= Ncx) cycle
                        if (cy2 < 0_int64 .or. cy2 >= Ncy) cycle
                        ci2 = cx2 + Ncx * cy2
                        j   = head(ci2)
                        do while (j /= -1)
                            if (j > i) then
                                jdx    = get_index(j)
                                dr_vec = y(jdx:jdx+1) - y(idx:idx+1)
                                dr2    = dr_vec(1)*dr_vec(1) + dr_vec(2)*dr_vec(2)
                                if (dr2 < r_cut2) then
                                    npairs = npairs + 1
                                    tmp_list(1, npairs) = i
                                    tmp_list(2, npairs) = j
                                end if
                            end if
                            j = next(j)
                        end do
                    end do
                end do
            end do

            deallocate(head, next)
        end if

        if (allocated(vlist))     deallocate(vlist)
        if (allocated(vlist_pos)) deallocate(vlist_pos)
        if (allocated(vlist_r))   deallocate(vlist_r)

        allocate(vlist(2, npairs))
        vlist(:, 1:npairs) = tmp_list(:, 1:npairs)
        vlist_n = npairs
        deallocate(tmp_list)

        allocate(vlist_pos(2, first_index:last_index))
        allocate(vlist_r(first_index:last_index))
        do i = first_index, last_index
            idx           = get_index(i)
            vlist_pos(1, i) = y(idx)
            vlist_pos(2, i) = y(idx+1)
            vlist_r(i)      = sqrt(y(idx)**2 + y(idx+1)**2)
        end do

        vlist_built = .True.
        vlist_time  = time

        !$OMP ATOMIC UPDATE
        verlet_rebuilds = verlet_rebuilds + 1_int64

    end subroutine build_verlet_list


    !> Verlet pair loop — PARALLEL over k (pairs).
    !> Each pair calls soft_sphere_force which uses ATOMIC writes internally.
    !> Rebuild check is serial (modifies global saved state).
    subroutine collisions_verlet(time, y, der, first_index, last_index, gamma_n, gamma_t, are_moons)
        implicit none
        real(wp),               intent(in)    :: time
        real(wp), dimension(:), intent(in)    :: y
        real(wp), dimension(:), intent(inout) :: der
        integer(kind=4),        intent(in)    :: first_index, last_index
        real(wp),               intent(in)    :: gamma_n, gamma_t
        logical,                intent(in)    :: are_moons

        integer(kind=4) :: i, j, k, idx
        real(wp)        :: r_skin2, dx, dy, drift2, dt
        real(wp), dimension(2) :: xyr

        ! ── Rebuild check (serial — touches global saved state) ─────────────
        if (.not. vlist_built) then
            call build_verlet_list(time, y, first_index, last_index, are_moons)
        else if (size(vlist_pos, 2) /= last_index - first_index + 1) then
            call build_verlet_list(time, y, first_index, last_index, are_moons)
        else
            r_skin2 = (vrcut * sim%verlet_skin_factor * uno2)**2
            dt      = time - vlist_time

            do i = first_index, last_index
                idx    = get_index(i)
                xyr    = vlist_pos(:, i)
                call get_xy_rotated(xyr, vlist_r(i), dt, G*m_arr(1), y(2))
                dx     = y(idx)   - xyr(1)
                dy     = y(idx+1) - xyr(2)
                drift2 = dx*dx + dy*dy
                if (drift2 > r_skin2) then
                    call build_verlet_list(time, y, first_index, last_index, are_moons)
                    exit
                end if
            end do
        end if

        ! ── Pair force loop (PARALLEL over k) ──────
        !$OMP PARALLEL DO              &
        !$OMP   DEFAULT(SHARED)        &
        !$OMP   PRIVATE(k, i, j)       &
        !$OMP   SCHEDULE(DYNAMIC, 640)
        do k = 1, vlist_n
            i = vlist(1, k)
            j = vlist(2, k)
            call soft_sphere_force(y, der, i, j, gamma_n, gamma_t, are_moons)
        end do
        !$OMP END PARALLEL DO

    end subroutine collisions_verlet


    pure subroutine get_xy_rotated_inertial(xy, r, dt, GM, dummy)
        implicit none
        real(wp), dimension(2), intent(inout) :: xy
        real(wp),               intent(in)    :: r, dt, GM, dummy
        real(wp) :: n, dphi
        n    = sqrt(GM / r**3)
        dphi = n * dt
        xy   = rotate2D(xy, dphi)
    end subroutine get_xy_rotated_inertial

    pure subroutine get_xy_rotated_sinodic(xy, r, dt, GM, omega)
        implicit none
        real(wp), dimension(2), intent(inout) :: xy
        real(wp),               intent(in)    :: r, dt, GM, omega
        real(wp) :: n, dphi
        n    = sqrt(GM / r**3)
        dphi = (n - omega) * dt
        xy   = rotate2D(xy, dphi)
    end subroutine get_xy_rotated_sinodic

end module collisions
