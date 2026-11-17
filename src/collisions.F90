!> Module with main collision routines.
module collisions
    use iso_fortran_env, only: int64
    use omp_lib
    use constants, only: wp, G, cero, uno, uno2, dos, pi, twopi, tini
    use auxiliary, only: rotate2D
    use parameters, only: sim, m_arr, R_arr

    private
    public :: init_collisions, set_coll_parameters, &
            & collisions_brute, collisions_grid, collisions_verlet, &
            & verlet_rebuilds, verlet_caches, list_collided, &
            & collisions_drift, mr_steps, mr_active_sum, mr_pairs_sum, mr_missed, mr_rebuilds


    integer(kind=4), save :: vlist_n
    real(wp), save :: vskin  ! Verlet skin: pairs are listed up to (contact distance + vskin)
    logical, save :: is_sinodic = .False.  ! Coordinates given in the frame rotating with the asteroid
    integer(kind=4), save, allocatable :: vlist(:, :)
    real(wp), save, allocatable :: vlist_pos(:, :)
    real(wp), save, allocatable :: vlist_r(:)
    real(wp), save, allocatable :: vlist_n_kep(:)
    real(wp), save :: vlist_mean_n_kep
    real(wp), save :: vlist_time
    logical, save :: vlist_built = .False.
    logical, save, allocatable :: list_collided(:)

    integer(int64), save :: verlet_rebuilds = 0_int64
    integer(int64), save :: verlet_caches = 0_int64

    ! ── Multi-rate collision stepping (collisions_drift) ──────────────────
    logical, save, allocatable :: mr_is_active(:)        ! Body takes collision sub-steps in this step
    integer(kind=4), save, allocatable :: mr_act(:)      ! Active bodies
    integer(kind=4), save :: mr_n_act = 0
    integer(kind=4), save, allocatable :: mr_pairs(:, :) ! Pairs that may touch during the step (moon pairs first)
    integer(kind=4), save :: mr_np = 0, mr_np_moon = 0
    real(wp), save, allocatable :: mr_acc(:, :)          ! Contact accelerations, per body
    integer(int64), save, allocatable :: mr_prev_keys(:) ! Sorted pair keys of the previous step
    integer(kind=4), save :: mr_n_prev = 0, mr_prev_Ntotal = -1
    ! Sorted cell search (mr_pairs_grid)
    integer(int64), save, allocatable :: gkey(:), gktmp(:), gcx(:), gcy(:), sk(:)
    integer(kind=4), save, allocatable :: gord(:), gtmp(:), sid(:)
    real(wp), save, allocatable :: gu(:), sx(:, :), srs(:)
    integer(kind=4), save, allocatable :: mr_cand(:, :)  ! Candidate pairs for the look-ahead horizon (moon pairs first)
    integer(kind=4), save :: mr_nc = 0, mr_nc_moon = 0
    real(wp), save :: mr_horizon_left = cero             ! Time left before the candidates must be rebuilt
    real(wp), save :: mr_u_typ = cero                    ! Typical random speed at the last rebuild
    integer(kind=4), save :: mr_cand_Ntotal = -1
    integer(int64), save :: mr_rebuilds = 0_int64
    real(wp), parameter :: MR_MAX_NT = 0.07_wp           ! Look-ahead horizon: n_max * T <= this
    integer(kind=4), parameter :: MR_MAX_K = 8           ! ... and at most this many steps
    real(wp), parameter :: MR_FIXED_MARGIN = 0.05_wp     ! Fixed part of the margin, in contact distances
    ! Diagnostics
    integer(int64), save :: mr_steps = 0_int64       ! Large steps taken
    integer(int64), save :: mr_active_sum = 0_int64  ! Sum over steps of active bodies
    integer(int64), save :: mr_pairs_sum = 0_int64   ! Sum over steps of candidate pairs
    integer(int64), save :: mr_missed = 0_int64      ! Contacts found already started without the pair being listed

contains

    !> Index of body i in the state array (local copy so that the compiler can inline it:
    !> procedures from another module are not inlined without link-time optimization)
    pure elemental function get_index(i) result(idx)
        implicit none
        integer(kind=4), intent(in) :: i
        integer(kind=4) :: idx
        idx = 4*i - 1
    end function get_index

    subroutine init_collisions(N_total, sinodic)
        implicit none
        integer(kind=4), intent(in) :: N_total
        logical,         intent(in) :: sinodic

        if (allocated(list_collided)) deallocate(list_collided)
        allocate(list_collided(N_total))
        list_collided = .False.

        ! In the sinodic frame the bodies rotate with (n - omega): used by the Verlet drift check
        is_sinodic = sinodic

        ! Any previous Verlet list is no longer valid
        vlist_built = .False.

    end subroutine init_collisions

    subroutine set_coll_parameters(kappa_n, gamma_n, epsilon, Tdur, gamma_t, beta, explicit)
        implicit none
        real(wp), intent(inout) :: kappa_n, gamma_n
        real(wp), intent(inout) :: epsilon, Tdur
        real(wp), intent(inout) :: gamma_t
        real(wp), intent(in) :: beta
        logical, intent(in) :: explicit

        real(wp) :: dzeta, eta

        if (explicit) then
            if (kappa_n <= cero) return  ! No spring: nothing to derive (caller disables the collisions)
            dzeta = min(gamma_n / (dos * sqrt(kappa_n)), uno - tini)
            eta = dzeta * pi / sqrt(uno - dzeta*dzeta)
            epsilon = exp(-eta)
            Tdur = sqrt((pi*pi + eta*eta) / kappa_n)
        else
            ! epsilon = 0 would need an infinite damping, and Tdur = 0 an infinite stiffness
            if (epsilon <= cero .or. epsilon > uno .or. Tdur <= cero) then
                kappa_n = cero
                gamma_n = cero
                gamma_t = cero
                return  ! Caller disables the collisions (kappa = 0)
            end if
            eta = -log(epsilon)
            kappa_n = (pi*pi + eta*eta) / Tdur**2
            gamma_n = dos * eta / Tdur
        end if

        gamma_t = beta * gamma_n

    end subroutine set_coll_parameters

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!! SOFT-SPHERE COLLISION ROUTINES !!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !> Soft-sphere force kernel with ATOMIC writes on all der and list_collided.
    !> Safe to call from any parallel region.
    !> Soft-sphere contact force between bodies i and j (no side effects).
    !> touching = .False. if they are not in contact (F_vec is then zero).
    !> The force on j is +F_vec and on i is -F_vec; mi and mj are the masses to divide by.
    pure subroutine pair_force(y, i, j, are_moons, F_vec, mi, mj, touching)
        implicit none
        real(wp), dimension(:), intent(in) :: y
        integer(kind=4),        intent(in) :: i, j
        logical,                intent(in) :: are_moons
        real(wp),               intent(out) :: F_vec(2), mi, mj
        logical,                intent(out) :: touching

        integer(kind=4) :: idx, jdx
        real(wp) :: dr_vec(2), dr2, dr, rcoll, overlap
        real(wp) :: dr_ver(2), dt_ver(2)
        real(wp) :: dv_vec(2), dvr, dvt
        real(wp) :: F_n, F_t
        real(wp) :: kappa_pair, gamma_n_pair, gamma_t_pair
        real(wp) :: meff
        real(wp) :: rand1, rn

        F_vec = cero
        mi = uno
        mj = uno
        touching = .False.

        idx = get_index(i)
        jdx = get_index(j)

        dr_vec(1) = y(jdx) - y(idx)
        dr_vec(2) = y(jdx+1) - y(idx+1)

        dr2 = dr_vec(1)*dr_vec(1) + dr_vec(2)*dr_vec(2)
        rcoll = R_arr(i) + R_arr(j)

        if (dr2 >= rcoll*rcoll) return
        touching = .True.

        dv_vec(1) = y(jdx+2) - y(idx+2)
        dv_vec(2) = y(jdx+3) - y(idx+3)

        if (dr2 < tini) then

            rn = dv_vec(1)*dv_vec(1) + dv_vec(2)*dv_vec(2)
            if (rn > tini) then
                rn = sqrt(rn)
                dr_ver = dv_vec / rn
            else
                ! Coincident and at relative rest: any direction works. Take one from the pair
                ! indices (golden-ratio sequence) so that runs stay deterministic, also with OpenMP
                rand1 = twopi * modulo(0.6180339887498949_wp * real(i + 7*j, wp), uno)
                dr_ver(1) = cos(rand1)
                dr_ver(2) = sin(rand1)
            end if

            dr = tini
            overlap = rcoll

        else

            dr = sqrt(dr2)
            dr_ver = dr_vec / dr
            overlap = rcoll - dr 

        end if

        ! Minimum effective separation: dr_eff = max(dr, dr_factor * rcoll), which caps the overlap
        ! (and so the force) for deeply interpenetrated bodies
        overlap = min(overlap, (uno - sim%dr_factor_col) * rcoll)
        dt_ver(1) = -dr_ver(2)
        dt_ver(2) = dr_ver(1)

        if (are_moons) then
            mi = m_arr(i)
            mj = m_arr(j)
            meff = mi * mj / (mi + mj)
            kappa_pair = sim%kappa_col_moon
            gamma_n_pair = sim%gamma_col_moon_n
            gamma_t_pair = sim%gamma_col_moon_t
        else
            mi = uno
            mj = uno
            meff = uno2
            kappa_pair = sim%kappa_col_part
            gamma_n_pair = sim%gamma_col_part_n
            gamma_t_pair = sim%gamma_col_part_t
        end if

        dvr = dv_vec(1)*dr_ver(1) + dv_vec(2)*dr_ver(2)
        dvt = dv_vec(1)*dt_ver(1) + dv_vec(2)*dt_ver(2)

        ! In the rotating frame v_inertial = v + omega z x r: the relative velocity gains omega*dr along
        ! the tangent (the normal part is unchanged). Friction must act on the inertial sliding velocity,
        ! or it would drag the pair towards co-rotation with the asteroid.
        if (is_sinodic) dvt = dvt + y(2) * dr

        ! Spring + dashpot. While separating the dashpot pulls the bodies together:
        ! if only repulsion is wanted, the normal force is clipped at zero
        F_n = kappa_pair * overlap - gamma_n_pair * dvr
        if (sim%use_only_repulsive_col) F_n = max(F_n, cero)

        F_t = cero
        if ((gamma_t_pair > cero) .and. (sim%coulomb_mu_col > cero) .and. (F_n > cero)) then
            F_t = -gamma_t_pair * dvt
            F_t = sign(min(abs(F_t), sim%coulomb_mu_col * F_n), F_t)
        end if

        F_vec = meff * (F_n * dr_ver + F_t * dt_ver)

    end subroutine pair_force


    !> Soft-sphere force kernel with ATOMIC writes on all der and list_collided.
    !> Safe to call from any parallel region.
    subroutine soft_sphere_force(y, der, i, j, are_moons)
        implicit none
        real(wp), dimension(:), intent(in) :: y
        real(wp), dimension(:), intent(inout) :: der
        integer(kind=4),        intent(in) :: i, j
        logical,                intent(in) :: are_moons

        integer(kind=4) :: idx, jdx
        real(wp) :: F_vec(2), mi, mj
        logical :: touching

        call pair_force(y, i, j, are_moons, F_vec, mi, mj, touching)
        if (.not. touching) return

        idx = get_index(i)
        jdx = get_index(j)

        ! ── Newton's 3rd law — ATOMIC on all four velocity slots ────────────
        !$OMP ATOMIC UPDATE
        der(jdx+2) = der(jdx+2) + F_vec(1) / mj
        !$OMP ATOMIC UPDATE
        der(jdx+3) = der(jdx+3) + F_vec(2) / mj
        !$OMP ATOMIC UPDATE
        der(idx+2) = der(idx+2) - F_vec(1) / mi
        !$OMP ATOMIC UPDATE
        der(idx+3) = der(idx+3) - F_vec(2) / mi

        ! Not need to  be ATOMIC bc only changes to TRUE
        list_collided(i) = .True.
        list_collided(j) = .True.

    end subroutine soft_sphere_force


    !> Brute-force O(N²).
    !> Outer i loop is PARALLEL — atomics inside soft_sphere_force handle conflicts.
    !> DYNAMIC schedule because the triangular inner loop makes load uneven.
    subroutine collisions_brute(y, der, first_index, last_index, are_moons)
        implicit none
        real(wp), dimension(:), intent(in) :: y
        real(wp), dimension(:), intent(inout) :: der
        integer(kind=4),        intent(in) :: first_index, last_index
        logical,                intent(in) :: are_moons

        integer(kind=4) :: i, j

        !$OMP PARALLEL DO                &
        !$OMP DEFAULT(SHARED)            &
        !$OMP PRIVATE(i, j)              &
        !$OMP SCHEDULE(GUIDED,64)
        do i = first_index, last_index - 1
            do j = i + 1, last_index
                call soft_sphere_force(y, der, i, j, are_moons)
            end do
        end do
        !$OMP END PARALLEL DO

    end subroutine collisions_brute


    !> Cell-list O(N).
    !> The linked-list traversal (head/next) is not thread-safe to build in
    !> parallel, so this routine stays serial. The force kernel itself is called
    !> serially here — if this loop ever becomes the bottleneck, use colouring.
    subroutine collisions_grid(y, der, first_index, last_index, are_moons)
        implicit none
        real(wp), dimension(:), intent(in) :: y
        real(wp), dimension(:), intent(inout) :: der
        integer(kind=4),        intent(in) :: first_index, last_index
        logical,                intent(in) :: are_moons

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

        !$OMP PARALLEL DO &
        !$OMP DEFAULT(SHARED) &
        !$OMP PRIVATE(i, idx) &
        !$OMP REDUCTION(min:x_min,y_min) &
        !$OMP REDUCTION(max:x_max,y_max)
        do i = first_index, last_index
            idx = get_index(i)

            x_min = min(x_min, y(idx))
            x_max = max(x_max, y(idx))
            y_min = min(y_min, y(idx+1))
            y_max = max(y_max, y(idx+1))
        end do
        !$OMP END PARALLEL DO

        x_min = x_min - L_cell
        x_max = x_max + L_cell
        y_min = y_min - L_cell
        y_max = y_max + L_cell

        Ncx = max(1_int64, ceiling((x_max - x_min) / L_cell, int64))
        Ncy = max(1_int64, ceiling((y_max - y_min) / L_cell, int64))

        ! Too many cells: enlarge them instead of falling back to an O(N^2) search.
        ! Cells wider than the cut-off still find every pair (only with more candidates per cell).
        do while (Ncx > sim%grid_col_max_cells / Ncy)
            L_cell = max(1.1_wp * L_cell, sqrt((x_max - x_min) * (y_max - y_min) / real(sim%grid_col_max_cells, wp)))
            Ncx = max(1_int64, ceiling((x_max - x_min) / L_cell, int64))
            Ncy = max(1_int64, ceiling((y_max - y_min) / L_cell, int64))
        end do


        Nc = Ncx * Ncy
        allocate(head(0:Nc-1), next(first_index:last_index))
        head = -1

        ! Create linked lists serially — if this becomes a bottleneck, use colouring to parallelise.
        do i = first_index, last_index
            idx = get_index(i)
            cx = min(int((y(idx)   - x_min) / L_cell, int64), Ncx - 1_int64)
            cy = min(int((y(idx+1) - y_min) / L_cell, int64), Ncy - 1_int64)
            cx = max(cx, 0_int64)
            cy = max(cy, 0_int64)
            ci = cx + Ncx * cy
            next(i) = head(ci)
            head(ci) = i
        end do

        !$OMP PARALLEL DO                                         &
        !$OMP DEFAULT(SHARED)                                     &
        !$OMP PRIVATE(i, idx, cx, cy, cx2, cy2, ci2, j, dcx, dcy) &
        !$OMP SCHEDULE(DYNAMIC, 32)
        do i = first_index, last_index
            idx = get_index(i)
            cx = min(max(int((y(idx)   - x_min) / L_cell, int64), 0_int64), Ncx - 1_int64)
            cy = min(max(int((y(idx+1) - y_min) / L_cell, int64), 0_int64), Ncy - 1_int64)

            do dcy = -1, 1
                do dcx = -1, 1
                    cx2 = cx + dcx
                    cy2 = cy + dcy
                    if (cx2 < 0_int64 .or. cx2 >= Ncx) cycle
                    if (cy2 < 0_int64 .or. cy2 >= Ncy) cycle
                    ci2 = cx2 + Ncx * cy2
                    j = head(ci2)
                    do while (j /= -1)
                        if (j > i) then
                            call soft_sphere_force(y, der, i, j, are_moons)
                        end if
                        j = next(j)
                    end do
                end do
            end do
        end do
        !$OMP END PARALLEL DO

        deallocate(head, next)

    end subroutine collisions_grid


    !> Append a pair to a growable (2, :) list.
    subroutine push_pair(list, npairs, i, j)
        implicit none
        integer(kind=4), allocatable, intent(inout) :: list(:, :)
        integer(kind=4), intent(inout) :: npairs
        integer(kind=4), intent(in) :: i, j
        integer(kind=4), allocatable :: tmp(:, :)

        if (npairs == size(list, 2)) then
            allocate(tmp(2, 2*size(list, 2)))
            tmp(:, 1:npairs) = list(:, 1:npairs)
            call move_alloc(tmp, list)
        end if
        npairs = npairs + 1
        list(1, npairs) = i
        list(2, npairs) = j
    end subroutine push_pair


    !> Build Verlet neighbour list (serial — modifies global saved state).
    !> Pairs closer than r_cut = contact + skin are stored, where contact is the largest
    !> possible contact distance (2 Rmax) and skin = verlet_skin_factor * contact.
    !> The list stays valid while no body moves more than skin/2 (see collisions_verlet).
    subroutine build_verlet_list(time, y, first_index, last_index, are_moons)
        implicit none
        real(wp),               intent(in) :: time
        real(wp), dimension(:), intent(in) :: y
        integer(kind=4),        intent(in) :: first_index, last_index
        logical,                intent(in) :: are_moons

        integer(kind=4) :: i, j, idx, jdx, npairs, N_particles
        real(wp) :: contact, r_cut, r_cut2, dr_vec(2), dr2, L_cell
        integer(kind=4), allocatable :: tmp_list(:, :)

        real(wp) :: x_min, x_max, y_min, y_max
        integer(int64) :: Ncx, Ncy, Nc, cx, cy, ci, cx2, cy2, ci2
        integer(kind=4) :: dcx, dcy
        integer(kind=4), allocatable :: head(:), next(:)
        real(wp) :: sq_GM

        sq_GM = sqrt(G * m_arr(1))

        N_particles = last_index - first_index + 1

        if (are_moons) then
            contact = dos * maxval(R_arr(first_index:last_index))
        else
            contact = dos * R_arr(first_index)
        end if
        vskin = contact * sim%verlet_skin_factor
        r_cut = contact + vskin
        r_cut2 = r_cut * r_cut

        ! The cells must be at least r_cut wide: only the 8 surrounding cells are searched
        L_cell = max(r_cut, sim%grid_col_min_cell_size)

        ! Growable pair list (a full N(N-1)/2 allocation is neither needed nor affordable)
        allocate(tmp_list(2, max(16, 8*N_particles)))
        npairs = 0

        x_min = y(get_index(first_index))
        x_max = x_min
        y_min = y(get_index(first_index)+1)
        y_max = y_min

        !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,idx) &
        !$OMP REDUCTION(min:x_min,y_min) &
        !$OMP REDUCTION(max:x_max,y_max)
        do i = first_index, last_index
            idx = get_index(i)

            x_min = min(x_min, y(idx))
            x_max = max(x_max, y(idx))
            y_min = min(y_min, y(idx+1))
            y_max = max(y_max, y(idx+1))
        end do
        !$OMP END PARALLEL DO

        x_min = x_min - L_cell
        x_max = x_max + L_cell
        y_min = y_min - L_cell
        y_max = y_max + L_cell

        Ncx = max(1_int64, ceiling((x_max - x_min) / L_cell, int64))
        Ncy = max(1_int64, ceiling((y_max - y_min) / L_cell, int64))

        ! Too many cells: enlarge them instead of falling back to an O(N^2) search.
        ! Cells wider than the cut-off still find every pair (only with more candidates per cell).
        do while (Ncx > sim%grid_col_max_cells / Ncy)
            L_cell = max(1.1_wp * L_cell, sqrt((x_max - x_min) * (y_max - y_min) / real(sim%grid_col_max_cells, wp)))
            Ncx = max(1_int64, ceiling((x_max - x_min) / L_cell, int64))
            Ncy = max(1_int64, ceiling((y_max - y_min) / L_cell, int64))
        end do

        if (Ncx > sim%grid_col_max_cells / Ncy) then
            ! Not reached (cells are enlarged above); kept as a safe fallback
            do i = first_index, last_index - 1
                idx = get_index(i)
                do j = i + 1, last_index
                    jdx = get_index(j)
                    dr_vec = y(jdx:jdx+1) - y(idx:idx+1)
                    dr2 = dr_vec(1)*dr_vec(1) + dr_vec(2)*dr_vec(2)
                    if (dr2 < r_cut2) call push_pair(tmp_list, npairs, i, j)
                end do
            end do
        else
            Nc = Ncx * Ncy
            allocate(head(0:Nc-1), next(first_index:last_index))
            head = -1

            ! Create linked lists serially — if this becomes a bottleneck, use colouring to parallelise.
            do i = first_index, last_index
                idx = get_index(i)
                cx = min(int((y(idx)   - x_min) / L_cell, int64), Ncx - 1_int64)
                cy = min(int((y(idx+1) - y_min) / L_cell, int64), Ncy - 1_int64)
                cx = max(cx, 0_int64)
                cy = max(cy, 0_int64)
                ci = cx + Ncx * cy
                next(i) = head(ci)
                head(ci) = i
            end do

            do i = first_index, last_index
                idx = get_index(i)
                cx = min(max(int((y(idx)   - x_min) / L_cell, int64), 0_int64), Ncx - 1_int64)
                cy = min(max(int((y(idx+1) - y_min) / L_cell, int64), 0_int64), Ncy - 1_int64)

                do dcy = -1, 1
                    do dcx = -1, 1
                        cx2 = cx + dcx
                        cy2 = cy + dcy
                        if (cx2 < 0_int64 .or. cx2 >= Ncx) cycle
                        if (cy2 < 0_int64 .or. cy2 >= Ncy) cycle
                        ci2 = cx2 + Ncx * cy2
                        j = head(ci2)
                        do while (j /= -1)
                            if (j > i) then
                                jdx = get_index(j)
                                dr_vec = y(jdx:jdx+1) - y(idx:idx+1)
                                dr2 = dr_vec(1)*dr_vec(1) + dr_vec(2)*dr_vec(2)
                                if (dr2 < r_cut2) call push_pair(tmp_list, npairs, i, j)
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
        if (allocated(vlist_n_kep)) deallocate(vlist_n_kep)

        allocate(vlist(2, npairs))
        vlist(:, 1:npairs) = tmp_list(:, 1:npairs)
        vlist_n = npairs
        deallocate(tmp_list)

        allocate(vlist_pos(2, first_index:last_index))
        allocate(vlist_r(first_index:last_index))
        allocate(vlist_n_kep(first_index:last_index))

        !$OMP PARALLEL DO           &
        !$OMP DEFAULT(SHARED)       &
        !$OMP PRIVATE(i, idx)       &
        !$OMP SCHEDULE(STATIC)
        do i = first_index, last_index
            idx = get_index(i)
            vlist_pos(1, i) = y(idx)
            vlist_pos(2, i) = y(idx+1)
            vlist_r(i) = max(sqrt(y(idx)**2 + y(idx+1)**2), tini)
            vlist_n_kep(i) = sq_GM / vlist_r(i)**(1.5_wp)
        end do
        !$OMP END PARALLEL DO

        vlist_mean_n_kep = sum(vlist_n_kep(first_index:last_index)) / real(N_particles, wp)

        vlist_built = .True.
        vlist_time = time

        verlet_rebuilds = verlet_rebuilds + 1_int64

    end subroutine build_verlet_list


    !> Verlet pair loop — PARALLEL over k (pairs).
    !> Each pair calls soft_sphere_force which uses ATOMIC writes internally.
    !> Rebuild check is serial (modifies global saved state).
    !>
    !> The list is valid while no body has moved more than skin/2 from where it was when
    !> the list was built. Since a rigid rotation of the whole reference set keeps every
    !> pair distance, the reference positions are first rotated with the mean orbital
    !> motion, so that only the motion relative to it counts as drift.
    subroutine collisions_verlet(time, y, der, first_index, last_index, are_moons)
        implicit none
        real(wp),               intent(in) :: time
        real(wp), dimension(:), intent(in) :: y
        real(wp), dimension(:), intent(inout) :: der
        integer(kind=4),        intent(in) :: first_index, last_index
        logical,                intent(in) :: are_moons

        integer(kind=4) :: i, j, k, idx
        real(wp) :: r_skin2, dx, dy, dt, drift2_max
        real(wp) :: xr, yr
        real(wp) :: dphi, cos_dphi, sin_dphi
        logical :: rebuild

        ! ── Rebuild check (serial — touches global saved state) ─────────────
        rebuild = .not. vlist_built
        if (.not. rebuild) rebuild = (lbound(vlist_pos, 2) /= first_index) .or. (ubound(vlist_pos, 2) /= last_index)

        if (.not. rebuild) then
            r_skin2 = (vskin * uno2)**2
            dt = time - vlist_time
            if (is_sinodic) then
                dphi = (vlist_mean_n_kep - y(2)) * dt  ! Rotating frame: (n - omega)
            else
                dphi = vlist_mean_n_kep * dt
            end if
            cos_dphi = cos(dphi)
            sin_dphi = sin(dphi)
            drift2_max = cero

            !$OMP PARALLEL DO                  &
            !$OMP DEFAULT(SHARED)              &
            !$OMP PRIVATE(i, idx, xr, yr, dx, dy) &
            !$OMP REDUCTION(max:drift2_max)    &
            !$OMP SCHEDULE(STATIC)
            do i = first_index, last_index
                idx = get_index(i)
                xr = vlist_pos(1, i) * cos_dphi - vlist_pos(2, i) * sin_dphi
                yr = vlist_pos(1, i) * sin_dphi + vlist_pos(2, i) * cos_dphi
                dx = y(idx)   - xr
                dy = y(idx+1) - yr
                drift2_max = max(drift2_max, dx*dx + dy*dy)
            end do
            !$OMP END PARALLEL DO

            rebuild = drift2_max > r_skin2
        end if

        if (rebuild) call build_verlet_list(time, y, first_index, last_index, are_moons)

        ! ── Pair force loop (PARALLEL over k) ──────
        !$OMP PARALLEL DO            &
        !$OMP DEFAULT(SHARED)        &
        !$OMP PRIVATE(k, i, j)       &
        !$OMP SCHEDULE(STATIC)
        do k = 1, vlist_n
            i = vlist(1, k)
            j = vlist(2, k)
            call soft_sphere_force(y, der, i, j, are_moons)
        end do
        !$OMP END PARALLEL DO

        verlet_caches = verlet_caches + 1_int64

    end subroutine collisions_verlet


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!! MULTI-RATE COLLISION STEPPING !!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !> Can bodies i and j come into contact during a step dt?
    !> Free (straight-line) motion is assumed; the margin covers the skin and the velocity
    !> changes that collisions with third bodies can produce during the step.
    pure logical function may_touch(y, i, j, dt, u_typ)
        implicit none
        real(wp), dimension(:), intent(in) :: y
        integer(kind=4),        intent(in) :: i, j
        real(wp),               intent(in) :: dt, u_typ
        integer(kind=4) :: idx, jdx
        real(wp) :: dx(2), dv(2), dv2, ts, dmin(2), rcoll, margin, vrel

        idx = get_index(i)
        jdx = get_index(j)
        dx(1) = y(jdx) - y(idx)
        dx(2) = y(jdx+1) - y(idx+1)
        dv(1) = y(jdx+2) - y(idx+2)
        dv(2) = y(jdx+3) - y(idx+3)
        dv2 = dv(1)*dv(1) + dv(2)*dv(2)
        vrel = sqrt(dv2)
        rcoll = R_arr(i) + R_arr(j)
        ! Velocity changes during the step (a collision with a third body changes the velocity by
        ! at most about its relative speed, ~ u_typ), plus a small fixed fraction for round-off
        margin = MR_FIXED_MARGIN*rcoll + dt*(vrel + dos*u_typ)
        ts = cero
        if (dv2 > tini) ts = min(max(-(dx(1)*dv(1) + dx(2)*dv(2))/dv2, cero), dt)
        dmin = dx + dv*ts
        may_touch = (dmin(1)*dmin(1) + dmin(2)*dmin(2)) < (rcoll + margin)**2
    end function may_touch


    !> Append a pair to the look-ahead candidates
    subroutine mr_cand_add(i, j)
        implicit none
        integer(kind=4), intent(in) :: i, j
        integer(kind=4), allocatable :: tmp(:, :)
        if (mr_nc == size(mr_cand, 2)) then
            allocate(tmp(2, 2*size(mr_cand, 2)))
            tmp(:, 1:mr_nc) = mr_cand(:, 1:mr_nc)
            call move_alloc(tmp, mr_cand)
        end if
        mr_nc = mr_nc + 1
        mr_cand(1, mr_nc) = min(i, j)
        mr_cand(2, mr_nc) = max(i, j)
    end subroutine mr_cand_add


    !> Register a candidate pair; count it as a missed contact if the bodies are already
    !> overlapping although the pair was not listed in the previous step.
    subroutine mr_add_pair(y, i, j, N_total)
        implicit none
        real(wp), dimension(:), intent(in) :: y
        integer(kind=4), intent(in) :: i, j, N_total
        integer(kind=4), allocatable :: tmp(:, :)
        integer(kind=4) :: idx, jdx, a, b
        real(wp) :: dx, dy, rcoll
        integer(int64) :: key

        if (mr_np == size(mr_pairs, 2)) then
            allocate(tmp(2, 2*size(mr_pairs, 2)))
            tmp(:, 1:mr_np) = mr_pairs(:, 1:mr_np)
            call move_alloc(tmp, mr_pairs)
        end if
        mr_np = mr_np + 1
        mr_pairs(1, mr_np) = min(i, j)
        mr_pairs(2, mr_np) = max(i, j)

        if (mr_n_prev > 0) then
            idx = get_index(i)
            jdx = get_index(j)
            dx = y(jdx) - y(idx)
            dy = y(jdx+1) - y(idx+1)
            rcoll = R_arr(i) + R_arr(j)
            if (dx*dx + dy*dy < rcoll*rcoll) then
                key = int(min(i, j), int64)*int(N_total + 1, int64) + int(max(i, j), int64)
                ! Binary search in the previous (sorted) keys
                a = 1
                b = mr_n_prev
                do while (a < b)
                    if (mr_prev_keys((a + b)/2) < key) then
                        a = (a + b)/2 + 1
                    else
                        b = (a + b)/2
                    end if
                end do
                if (mr_prev_keys(a) /= key) mr_missed = mr_missed + 1_int64
            end if
        end if
    end subroutine mr_add_pair


    !> Candidate pairs among bodies first..last, by brute force (few bodies: moons)
    subroutine mr_pairs_brute(y, dt, first, last)
        implicit none
        real(wp), dimension(:), intent(in) :: y
        real(wp), intent(in) :: dt
        integer(kind=4), intent(in) :: first, last
        integer(kind=4) :: i, j
        do i = first, last - 1
            do j = i + 1, last
                if (may_touch(y, i, j, dt, cero)) call mr_cand_add(i, j)
            end do
        end do
    end subroutine mr_pairs_brute


    !> Candidate pairs among bodies first..last, with cells sorted in memory.
    !> Each body searches a radius that bounds how close any partner can come during dt
    !> (random speed relative to the circular flow, plus the Keplerian shear).
    !> The bodies are sorted by cell (row-major) with a radix sort and their data copied in that
    !> order, so each neighbouring row of cells is a contiguous range found by running pointers:
    !> the search is linear in the number of bodies and reads memory sequentially.
    subroutine mr_pairs_grid(y, dt, first, last)
        implicit none
        real(wp), dimension(:), intent(in) :: y
        real(wp), intent(in) :: dt
        integer(kind=4), intent(in) :: first, last
        integer(kind=4) :: i, j, idx, k, nb, p, q, dcy, row, lo, hi, mid
        integer(kind=4) :: ptr(-1:1)
        integer(int64) :: cxmin, cymin, cxmax, cymax, Wx, kp, kt, klo, khi
        real(wp) :: GMa, xa(4), rx, ry, r, vc, ux, uy, u2sum, u_typ, nmax, base, shear, L_cell
        real(wp) :: dx(2), dv(2), dv2, ts, dmin(2), rcoll, margin, rsi

        nb = last - first + 1
        if (nb < 2) return
        call mr_grid_alloc(nb)

        ! Random speed of each body relative to the circular flow around the asteroid
        GMa = G*m_arr(1)
        xa = y(3:6)
        u2sum = cero
        nmax = cero
        do i = first, last
            idx = get_index(i)
            rx = y(idx) - xa(1)
            ry = y(idx+1) - xa(2)
            r = max(sqrt(rx*rx + ry*ry), tini)
            vc = sqrt(GMa/r)
            ux = y(idx+2) - xa(3) + vc*ry/r
            uy = y(idx+3) - xa(4) - vc*rx/r
            gu(i - first + 1) = sqrt(ux*ux + uy*uy)
            u2sum = u2sum + gu(i - first + 1)**2
            nmax = max(nmax, vc/r)
        end do
        u_typ = sqrt(u2sum/real(nb, wp))
        mr_u_typ = u_typ

        ! |dv_ij| <= u_i + u_j + (velocity difference of the circular flow) <= u_i + u_j + 2.5 n |dx|
        base = dos*maxval(R_arr(first:last))*(uno + MR_FIXED_MARGIN)
        shear = max(uno - 5.0_wp*nmax*dt, uno2)
        ! Cells sized for a body with twice the typical random speed
        L_cell = max((base + dt*(10.0_wp*u_typ))/shear, sim%grid_col_min_cell_size)

        ! Cell keys (row-major)
        cxmin = huge(1_int64); cymin = huge(1_int64); cxmax = -huge(1_int64); cymax = -huge(1_int64)
        do i = first, last
            idx = get_index(i)
            p = i - first + 1
            gcx(p) = floor(y(idx)/L_cell, int64)
            gcy(p) = floor(y(idx+1)/L_cell, int64)
            cxmin = min(cxmin, gcx(p)); cxmax = max(cxmax, gcx(p))
            cymin = min(cymin, gcy(p)); cymax = max(cymax, gcy(p))
        end do
        Wx = cxmax - cxmin + 3_int64   ! Spare column on each side: neighbour keys never wrap rows
        do p = 1, nb
            gkey(p) = (gcy(p) - cymin + 1_int64)*Wx + (gcx(p) - cxmin + 1_int64)
            gord(p) = p
        end do
        call radix_sort_keys(nb)

        ! Data in sorted order
        do q = 1, nb
            p = gord(q)
            i = p + first - 1
            idx = get_index(i)
            sid(q) = i
            sk(q) = gkey(q)   ! gkey is already sorted (gord holds the permutation)
            sx(1, q) = y(idx);   sx(2, q) = y(idx+1)
            sx(3, q) = y(idx+2); sx(4, q) = y(idx+3)
            srs(q) = (base + dt*(4.0_wp*gu(p) + dos*u_typ))/shear
        end do

        ptr = 1
        do q = 1, nb
            i = sid(q)
            rsi = srs(q)
            k = max(1, ceiling(rsi/L_cell))
            kp = sk(q)
            do dcy = -k, k
                klo = kp + int(dcy, int64)*Wx - int(k, int64)
                khi = kp + int(dcy, int64)*Wx + int(k, int64)
                if (k == 1) then
                    ! Running pointer: klo grows with q, so the range start only moves forward
                    row = dcy
                    do while (ptr(row) <= nb)
                        if (sk(ptr(row)) >= klo) exit
                        ptr(row) = ptr(row) + 1
                    end do
                    lo = ptr(row)
                else
                    ! Fast body: binary search
                    lo = 1
                    hi = nb + 1
                    do while (lo < hi)
                        mid = (lo + hi)/2
                        if (sk(mid) < klo) then
                            lo = mid + 1
                        else
                            hi = mid
                        end if
                    end do
                end if
                do p = lo, nb
                    kt = sk(p)
                    if (kt > khi) exit
                    if (p == q) cycle
                    ! Each pair once, from the body with the larger search radius
                    j = sid(p)
                    if (.not. ((srs(p) < rsi) .or. ((srs(p) == rsi) .and. (j > i)))) cycle
                    ! Cheap rejection: no partner can come closer than its search radius allows
                    dx = sx(1:2, p) - sx(1:2, q)
                    if (dx(1)*dx(1) + dx(2)*dx(2) >= rsi*rsi) cycle
                    ! May they touch during dt? (same test as may_touch, on the sorted copies)
                    dv = sx(3:4, p) - sx(3:4, q)
                    dv2 = dv(1)*dv(1) + dv(2)*dv(2)
                    rcoll = R_arr(i) + R_arr(j)
                    margin = MR_FIXED_MARGIN*rcoll + dt*(sqrt(dv2) + dos*u_typ)
                    ts = cero
                    if (dv2 > tini) ts = min(max(-(dx(1)*dv(1) + dx(2)*dv(2))/dv2, cero), dt)
                    dmin = dx + dv*ts
                    if ((dmin(1)*dmin(1) + dmin(2)*dmin(2)) < (rcoll + margin)**2) call mr_cand_add(i, j)
                end do
            end do
        end do

    end subroutine mr_pairs_grid


    !> Work arrays of the sorted cell search
    subroutine mr_grid_alloc(nb)
        implicit none
        integer(kind=4), intent(in) :: nb
        if (allocated(gkey)) then
            if (size(gkey) >= nb) return
            deallocate(gkey, gord, gtmp, gktmp, gcx, gcy, gu, sid, sk, sx, srs)
        end if
        allocate(gkey(nb), gord(nb), gtmp(nb), gktmp(nb), gcx(nb), gcy(nb), gu(nb), &
               & sid(nb), sk(nb), sx(4, nb), srs(nb))
    end subroutine mr_grid_alloc


    !> LSD radix sort of gkey(1:n) (non-negative), carrying the permutation gord
    subroutine radix_sort_keys(n)
        implicit none
        integer(kind=4), intent(in) :: n
        integer(kind=4), parameter :: BITS = 11, NBUCK = 2**BITS
        integer(kind=4) :: cnt(0:NBUCK-1), i, b, shift
        integer(int64) :: kmax

        kmax = maxval(gkey(1:n))
        shift = 0
        do while (shiftr(kmax, shift) > 0_int64)
            cnt = 0
            do i = 1, n
                b = int(iand(shiftr(gkey(i), shift), int(NBUCK - 1, int64)))
                cnt(b) = cnt(b) + 1
            end do
            do b = 1, NBUCK - 1
                cnt(b) = cnt(b) + cnt(b - 1)
            end do
            do i = n, 1, -1
                b = int(iand(shiftr(gkey(i), shift), int(NBUCK - 1, int64)))
                gktmp(cnt(b)) = gkey(i)
                gtmp(cnt(b)) = gord(i)
                cnt(b) = cnt(b) - 1
            end do
            gkey(1:n) = gktmp(1:n)
            gord(1:n) = gtmp(1:n)
            shift = shift + BITS
        end do
    end subroutine radix_sort_keys


    !> Contact accelerations on the active bodies (mr_acc), from the candidate pairs
    subroutine mr_forces(y)
        implicit none
        real(wp), dimension(:), intent(in) :: y
        integer(kind=4) :: a, p, i, j
        real(wp) :: F_vec(2), mi, mj
        logical :: touching

        do a = 1, mr_n_act
            i = mr_act(a)
            mr_acc(1, i) = cero
            mr_acc(2, i) = cero
        end do
        do p = 1, mr_np
            i = mr_pairs(1, p)
            j = mr_pairs(2, p)
            call pair_force(y, i, j, p <= mr_np_moon, F_vec, mi, mj, touching)
            if (.not. touching) cycle
            mr_acc(1, j) = mr_acc(1, j) + F_vec(1)/mj
            mr_acc(2, j) = mr_acc(2, j) + F_vec(2)/mj
            mr_acc(1, i) = mr_acc(1, i) - F_vec(1)/mi
            mr_acc(2, i) = mr_acc(2, i) - F_vec(2)/mi
            list_collided(i) = .True.
            list_collided(j) = .True.
        end do
    end subroutine mr_forces


    !> In-place heap sort of int64 keys
    subroutine sort_keys(a, n)
        implicit none
        integer(int64), intent(inout) :: a(:)
        integer(kind=4), intent(in) :: n
        integer(kind=4) :: i, last
        integer(int64) :: tmp
        do i = n/2, 1, -1
            call sift(i, n)
        end do
        do last = n, 2, -1
            tmp = a(1); a(1) = a(last); a(last) = tmp
            call sift(1, last - 1)
        end do
    contains
        subroutine sift(start, nn)
            integer(kind=4), intent(in) :: start, nn
            integer(kind=4) :: root, child
            integer(int64) :: t
            root = start
            do while (2*root <= nn)
                child = 2*root
                if (child < nn) then
                    if (a(child) < a(child + 1)) child = child + 1
                end if
                if (a(root) < a(child)) then
                    t = a(root); a(root) = a(child); a(child) = t
                    root = child
                else
                    return
                end if
            end do
        end subroutine sift
    end subroutine sort_keys


    !> Steps covered by one candidate search.
    !> The straight-line relative motion assumed by the search must stay accurate (n_max T small).
    !> Within that, the search radius grows as base + T v (base ~ contact distance, v ~ random
    !> speeds), so the search cost per step goes as (base + T v)^2 / T: smallest at T v ~ base.
    integer(kind=4) function mr_horizon_steps(y, dt, last_moon, N_total)
        implicit none
        real(wp), dimension(:), intent(in) :: y
        real(wp), intent(in) :: dt
        integer(kind=4), intent(in) :: last_moon, N_total
        integer(kind=4) :: i, idx, k_shear, k_cost
        real(wp) :: r2min, rx, ry, nmax, base

        r2min = huge(cero)
        do i = 2, N_total
            idx = get_index(i)
            rx = y(idx) - y(3)
            ry = y(idx+1) - y(4)
            r2min = min(r2min, rx*rx + ry*ry)
        end do
        nmax = sqrt(G*m_arr(1)/max(r2min, tini)**1.5_wp)
        k_shear = max(1, min(MR_MAX_K, int(MR_MAX_NT/(nmax*dt))))

        k_cost = MR_MAX_K
        if ((N_total > last_moon) .and. (mr_u_typ > cero)) then
            base = dos*maxval(R_arr(last_moon+1:N_total))*(uno + MR_FIXED_MARGIN)
            k_cost = max(1, nint(base/(dt*10.0_wp*mr_u_typ)))  ! Cells are sized with 10 u_typ T
        end if
        mr_horizon_steps = min(k_shear, k_cost)
    end function mr_horizon_steps


    !> Advance the positions of all bodies by dt (the drift of a KDK step), resolving soft-sphere
    !> contacts with sub-steps of ~dt_small. Only the bodies that may touch during dt take the
    !> sub-steps (kick-drift-kick with the contact forces); all the others drift once.
    !> This is the same scheme as sub-stepping every body (their contact kicks would be zero),
    !> at a cost proportional to the bodies in or near contact instead of to all bodies.
    subroutine collisions_drift(y, dt, dt_small, last_moon, N_total)
        implicit none
        real(wp), dimension(:), intent(inout) :: y
        real(wp), intent(in) :: dt, dt_small
        integer(kind=4), intent(in) :: last_moon, N_total
        integer(kind=4) :: a, p, i, idx, N_sub, k
        real(wp) :: h, horizon

        ! Work arrays
        if (allocated(mr_is_active)) then
            if (size(mr_is_active) < N_total) deallocate(mr_is_active, mr_act, mr_acc)
        end if
        if (.not. allocated(mr_is_active)) then
            allocate(mr_is_active(N_total), mr_act(N_total), mr_acc(2, N_total))
            mr_is_active = .False.
            mr_acc = cero
            mr_n_act = 0
        end if
        if (.not. allocated(mr_pairs)) allocate(mr_pairs(2, 1024))
        if (.not. allocated(mr_cand)) allocate(mr_cand(2, 1024))
        if (.not. allocated(mr_prev_keys)) allocate(mr_prev_keys(1024))

        ! Body indices change when bodies are removed: the previous pair list is then meaningless
        if (N_total /= mr_prev_Ntotal) mr_n_prev = 0

        do a = 1, mr_n_act
            mr_is_active(mr_act(a)) = .False.
        end do
        mr_n_act = 0
        mr_np = 0
        list_collided = .False.

        ! ── Candidates: pairs that may touch within a look-ahead horizon of several steps ──
        ! (rebuilt when the horizon is used up, or when bodies were removed)
        if ((mr_horizon_left < dt*(uno - 1.e-9_wp)) .or. (N_total /= mr_cand_Ntotal)) then
            horizon = dt*real(mr_horizon_steps(y, dt, last_moon, N_total), wp)
            mr_nc = 0
            if (sim%use_moon_soft_sphere_col .and. (last_moon > 2)) call mr_pairs_brute(y, horizon, 2, last_moon)
            mr_nc_moon = mr_nc
            if (sim%use_part_soft_sphere_col) call mr_pairs_grid(y, horizon, last_moon + 1, N_total)
            mr_horizon_left = horizon
            mr_cand_Ntotal = N_total
            mr_rebuilds = mr_rebuilds + 1_int64
        end if
        mr_horizon_left = mr_horizon_left - dt

        ! ── Pairs that may touch during this step (from the candidates) ──
        do p = 1, mr_nc
            if (p == mr_nc_moon + 1) mr_np_moon = mr_np
            if (may_touch(y, mr_cand(1, p), mr_cand(2, p), dt, merge(cero, mr_u_typ, p <= mr_nc_moon))) &
                & call mr_add_pair(y, mr_cand(1, p), mr_cand(2, p), N_total)
        end do
        if (mr_nc_moon == mr_nc) mr_np_moon = mr_np

        do p = 1, mr_np
            do k = 1, 2
                i = mr_pairs(k, p)
                if (.not. mr_is_active(i)) then
                    mr_is_active(i) = .True.
                    mr_n_act = mr_n_act + 1
                    mr_act(mr_n_act) = i
                end if
            end do
        end do

        ! Keys of this step's pairs, for the missed-contact check of the next step
        if (size(mr_prev_keys) < mr_np) then
            deallocate(mr_prev_keys)
            allocate(mr_prev_keys(2*mr_np))
        end if
        do p = 1, mr_np
            mr_prev_keys(p) = int(mr_pairs(1, p), int64)*int(N_total + 1, int64) + int(mr_pairs(2, p), int64)
        end do
        mr_n_prev = mr_np
        if (mr_np > 1) call sort_keys(mr_prev_keys, mr_np)
        mr_prev_Ntotal = N_total
        if (mr_np == 0) mr_n_prev = 0

        ! ── Bodies that cannot touch: one drift ──
        do i = 1, N_total
            if (mr_is_active(i)) cycle
            idx = get_index(i)
            y(idx)   = y(idx)   + dt*y(idx+2)
            y(idx+1) = y(idx+1) + dt*y(idx+3)
        end do

        ! ── Bodies that may touch: KDK sub-steps with the contact forces ──
        if (mr_n_act > 0) then
            N_sub = max(1, nint(dt/dt_small))
            h = dt/real(N_sub, wp)
            do k = 1, N_sub
                call mr_forces(y)
                do a = 1, mr_n_act
                    i = mr_act(a)
                    idx = get_index(i)
                    y(idx+2:idx+3) = y(idx+2:idx+3) + uno2*h*mr_acc(:, i)
                    y(idx:idx+1) = y(idx:idx+1) + h*y(idx+2:idx+3)
                end do
                call mr_forces(y)
                do a = 1, mr_n_act
                    i = mr_act(a)
                    idx = get_index(i)
                    y(idx+2:idx+3) = y(idx+2:idx+3) + uno2*h*mr_acc(:, i)
                end do
            end do
        end if

        mr_steps = mr_steps + 1_int64
        mr_active_sum = mr_active_sum + int(mr_n_act, int64)
        mr_pairs_sum = mr_pairs_sum + int(mr_np, int64)

    end subroutine collisions_drift

end module collisions
