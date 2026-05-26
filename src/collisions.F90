!> Module with main collision routines
module collisions
    use iso_fortran_env, only: int64
    use constants, only: wp, G, cero, uno, uno2, dos, twopi, tini
    use auxiliary, only: rotate2D
    use parameters, only: sim, m_arr, R_arr, get_index

    private
    public :: init_collisions, collisions_brute, collisions_grid, collisions_verlet, verlet_rebuilds

    abstract interface
        ! Here must be every f_i defined explicitly
        subroutine get_xy_rotated_tem(xy, r, dt, mu, omega)
            import :: wp
            implicit none
            real(wp), dimension(2), intent(inout) :: xy
            real(wp), intent(in) :: r, dt, mu, omega
        end subroutine get_xy_rotated_tem

    end interface

    procedure(get_xy_rotated_tem), pointer :: get_xy_rotated => null()

    ! ── Verlet neighbour list (persistent across sub-steps) ── 
    ! r_skin: extra shell beyond 2R — pairs within (2R + r_skin) are listed.
    ! Rebuild triggered when any particle moves more than r_skin/2 since last build.
    integer(kind=4), save :: vlist_n                    ! number of pairs
    real(wp), save :: vrcut                             ! current Verlet cutoff (2R + skin)
    integer(kind=4), save, allocatable :: vlist(:, :)   ! (2, vlist_n) pair indicesx0
    real(wp), save, allocatable :: vlist_pos(:, :)      ! (2, N) positions at last build
    real(wp), save, allocatable :: vlist_r(:)           ! (vlist_n) pairwise distances at last build
    real(wp), save :: vlist_time                        ! current Verlet time
    logical, save :: vlist_built = .False.
    logical, save, allocatable :: list_collided(:)      ! Track collisions for diagnostics (size N, indexed by particle ID)
    
    integer(int64), save :: verlet_rebuilds = 0_int64

    contains

        subroutine init_collisions(N_total, sinodic)
            implicit none
            integer(kind=4), intent(in) :: N_total
            logical, intent(in) :: sinodic

            if (allocated(list_collided)) deallocate(list_collided)
            allocate(list_collided(N_total))
            list_collided = .False.

            if (sinodic) then
                get_xy_rotated => get_xy_rotated_sinodic
            else
                get_xy_rotated => get_xy_rotated_inertial
            end if

        end subroutine init_collisions

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!! SOFT-SPHERE COLLISION ROUTINES !!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
        ! ── Shared force kernel ───
        ! Computes the soft-sphere force between two overlapping particles and
        ! accumulates it into der (Newton's 3rd law).
        !
        ! Normal force:
        !   F_n = kappa * delta - gamma_n * dvr    (clamped to >= 0)
        !
        ! Tangential force (Coulomb-limited friction):
        !   F_t = -gamma_t * dv_tan               (viscous sliding damping)
        !   |F_t| <= coulomb_mu * |F_n|           (Coulomb cap)
        !   kappa_t * delta_t not included: static tangential spring needs history
        !
        ! Called by both brute-force and grid routines.
        subroutine soft_sphere_force(y, der, i, j, gamma_n, gamma_t, are_moons)
            implicit none
            real(wp), dimension(:), intent(in) :: y
            real(wp), dimension(:), intent(inout) :: der
            integer(kind=4), intent(in) :: i, j
            real(wp), intent(in) :: gamma_n, gamma_t
            logical, intent(in) :: are_moons

            integer(kind=4) :: idx, jdx
            real(wp) :: dr_vec(2), dr2, dr, rcoll, overlap
            real(wp) :: dr_ver(2), dt_ver(2)          ! normal and tangential unit vectors
            real(wp) :: dv_vec(2), dvr, dvt           ! relative vel components
            real(wp) :: F_n, F_t                      ! force magnitudes
            real(wp) :: F_vec(2)
            real(wp) :: gamma_n_pair, gamma_t_pair
            real(wp) :: mi, mj
            real(wp) :: rand1, rn                     ! for random fallback direction
            real(wp) :: aux_real

            idx = get_index(i)
            jdx = get_index(j)

            ! ── Relative position ───
            dr_vec = y(jdx:jdx+1) - y(idx:idx+1)
            dr2 = dr_vec(1)*dr_vec(1) + dr_vec(2)*dr_vec(2)
            rcoll = R_arr(i) + R_arr(j)

            if (dr2 >= rcoll*rcoll) return  ! no overlap

            ! ── Zero-distance guard: random direction from previous relative vel ──
            if (dr2 < tini) then
                ! Try to use relative velocity direction as contact normal
                dv_vec = y(jdx+2:jdx+3) - y(idx+2:idx+3)
                rn = dv_vec(1)*dv_vec(1) + dv_vec(2)*dv_vec(2)
                if (rn > tini) then
                    rn = sqrt(rn)
                    dr_ver = dv_vec / rn          ! direction of approach
                else
                    ! Truly degenerate: random unit vector
                    call random_number(rand1)
                    rand1 = twopi * rand1  ! uniform angle in [0, 2pi)
                    dr_ver = [cos(rand1), sin(rand1)]
                end if
                dr = tini
                overlap = rcoll
            else
                dr = sqrt(dr2)
                overlap = rcoll - dr
                dr_ver = dr_vec / dr
            end if
            

            ! ── Tangential unit vector (90° rotation of normal) ───
            ! dt_ver points in the direction of tangential sliding.
            dt_ver = [-dr_ver(2), dr_ver(1)]

            ! ── Masses ────
            if (are_moons) then
                mi = m_arr(i)
                mj = m_arr(j)
                aux_real = dos * sqrt(sim%kappa_col_moon * (mi * mj / (mi + mj)))
                gamma_n_pair = min(sim%gamma_col_moon_n, uno) * aux_real
                gamma_t_pair = min(sim%gamma_col_moon_t, uno) * aux_real
            else
                mi = uno
                mj = uno
                gamma_n_pair = gamma_n
                gamma_t_pair = gamma_t
            end if

            ! ── Relative velocity components ───
            dv_vec = y(jdx+2:jdx+3) - y(idx+2:idx+3)
            dvr = dv_vec(1)*dr_ver(1) + dv_vec(2)*dr_ver(2)  ! normal component
            dvt = dv_vec(1)*dt_ver(1) + dv_vec(2)*dt_ver(2)  ! tangential component

            ! ── Normal force (spring + damping, clamped >= 0) ───
            ! Clamp ensures force is always repulsive — never attractive.
            F_n = sim%kappa_col_part * overlap
            if (gamma_n_pair > cero) then
                F_n = F_n - gamma_n_pair * min(dvr, cero)
            end if
            F_n = max(F_n, cero) ! ← clamp: F_n cannot go negative

            ! ── Tangential force (viscous, Coulomb-limited) ───
            ! Viscous sliding: F_t = -gamma_t * dv_tan
            ! Coulomb cap: |F_t| <= mu * F_n  (no friction beyond this)
            F_t = cero
            if ((gamma_t_pair > cero) .and. (sim%coulomb_mu_col > cero) .and. (F_n > cero)) then
                F_t = -gamma_t_pair * dvt
                ! Coulomb cap
                F_t = sign(min(abs(F_t), sim%coulomb_mu_col * F_n), F_t)
            end if

            ! ── Assemble total force on j ───
            F_vec = F_n * dr_ver + F_t * dt_ver

            ! ── Newton's 3rd law ────
            der(jdx+2:jdx+3) = der(jdx+2:jdx+3) + F_vec / mj
            der(idx+2:idx+3) = der(idx+2:idx+3) - F_vec / mi

            list_collided(i) = .True.
            list_collided(j) = .True.

        end subroutine soft_sphere_force
    

        ! ── Brute-force O(N²) — used when N_particles <= sim%grid_col_min_bodies ───
        subroutine collisions_brute(y, der, first_index, last_index, gamma_n, gamma_t, are_moons)
            implicit none
            real(wp), dimension(:), intent(in) :: y
            real(wp), dimension(:), intent(inout) :: der
            integer(kind=4), intent(in) :: first_index, last_index
            real(wp), intent(in) :: gamma_n, gamma_t
            logical, intent(in) :: are_moons
    
            integer(kind=4) :: i, j
    
            do i = first_index, last_index - 1
                do j = i + 1, last_index
                    call soft_sphere_force(y, der, i, j, gamma_n, gamma_t, are_moons)
                end do
            end do
    
        end subroutine collisions_brute
    

        ! ── Cell-list O(N) — used when N_particles > sim%grid_col_min_bodies ──
        ! Cell size ~ 2R (collision diameter): only the 9 surrounding cells
        ! need to be checked for each particle.
        subroutine collisions_grid(y, der, first_index, last_index, gamma_n, gamma_t, are_moons)
            implicit none
            real(wp), dimension(:), intent(in) :: y
            real(wp), dimension(:), intent(inout) :: der
            integer(kind=4), intent(in) :: first_index, last_index
            real(wp), intent(in) :: gamma_n, gamma_t
            logical, intent(in) :: are_moons
    
            real(wp) :: L_cell, x_min, x_max, y_min, y_max
            integer(int64) :: Ncx, Ncy, Nc
            integer(int64) :: cx, cy, ci
            integer(int64) :: cx2, cy2, ci2
            integer(kind=4) :: dcx, dcy
            integer(kind=4) :: i, j, idx
            integer(kind=4), allocatable :: head(:), next(:)
    
            ! ── Cell size ~ 2R (collision diameter): only the 9 surrounding cells
            if (are_moons) then
                L_cell = dos * maxval(R_arr(first_index:last_index))
            else
                L_cell = dos * R_arr(first_index)
            end if

            L_cell = max(L_cell, sim%grid_col_min_cell_size)  ! Avoid too small cells (too much overhead)

            ! ── Domain bounds from particle positions ───
            x_min = y(get_index(first_index))
            x_max = x_min
            y_min = y(get_index(first_index)+1)
            y_max = y_min

            do i = first_index, last_index
                idx = get_index(i)
                x_min = min(x_min, y(idx))
                x_max = max(x_max, y(idx))
                y_min = min(y_min, y(idx+1))
                y_max = max(y_max, y(idx+1))
            end do

            ! Add one cell of padding so boundary particles aren't clipped
            x_min = x_min - L_cell
            x_max = x_max + L_cell
            y_min = y_min - L_cell
            y_max = y_max + L_cell
    
            Ncx = max(1, ceiling((x_max - x_min) / L_cell))
            Ncy = max(1, ceiling((y_max - y_min) / L_cell))

            ! ── Safety: fall back to brute-force if grid is too large ──
            if (Ncx > sim%grid_col_max_cells / Ncy) then
                call collisions_brute(y, der, first_index, last_index, gamma_n, gamma_t, are_moons)
                return
            end if

            Nc = Ncx * Ncy
    
            ! ── Build linked-list cell structure ───
            allocate(head(0:Nc-1), next(first_index:last_index))
            head = -1  ! -1 = empty cell
    
            do i = first_index, last_index
                idx = get_index(i)
                cx = min(int((y(idx)   - x_min) / L_cell, int64), Ncx - 1_int64)
                cy = min(int((y(idx+1) - y_min) / L_cell, int64), Ncy - 1_int64)
                ! Clamp particles outside domain to boundary cells
                cx = max(cx, 0_int64)
                cy = max(cy, 0_int64)
                ci = cx + Ncx * cy
                next(i)  = head(ci)
                head(ci) = i
            end do
    
            ! ── Check only 9-cell neighbourhood ───
            do i = first_index, last_index
                idx = get_index(i)
                cx = min(max(int((y(idx) - x_min) / L_cell, int64), 0_int64), Ncx - 1_int64)
                cy = min(max(int((y(idx+1) - y_min) / L_cell, int64), 0_int64), Ncy - 1_int64)
    
                do dcy = -1, 1
                    do dcx = -1, 1
                        cx2 = cx + dcx
                        cy2 = cy + dcy
                        if (cx2 < 0_int64 .or. cx2 >= Ncx) cycle
                        if (cy2 < 0_int64 .or. cy2 >= Ncy) cycle
        
                        ci2 = cx2 + Ncx * cy2
                        j   = head(ci2)
                        do while (j /= -1)
                            if (j > i) then  ! avoid double-counting
                                call soft_sphere_force(y, der, i, j, gamma_n, gamma_t, are_moons)
                            end if
                            j = next(j)
                        end do
                    end do
                end do
            end do
    
            deallocate(head, next)
    
        end subroutine collisions_grid


        ! ── Verlet neighbour list ────
        ! Builds a list of all pairs within (2R + skin). The list is reused across
        ! sub-steps and only rebuilt when any particle has moved > skin/2 since
        ! the last build. This amortizes the O(N²) build cost over many sub-steps.
        subroutine build_verlet_list(time, y, first_index, last_index, are_moons)
            implicit none
            real(wp), intent(in) :: time
            real(wp), dimension(:), intent(in) :: y
            integer(kind=4), intent(in) :: first_index, last_index
            logical, intent(in) :: are_moons

            integer(kind=4) :: i, j, idx, jdx, npairs, N_particles
            real(wp) :: r_cut, r_cut2, dr_vec(2), dr2
            integer(kind=4), allocatable :: tmp_list(:, :)

            real(wp) :: x_min, x_max, y_min, y_max
            integer(int64) :: Ncx, Ncy, Nc
            integer(int64) :: cx, cy, ci
            integer(int64) :: cx2, cy2, ci2
            integer(kind=4) :: dcx, dcy
            integer(kind=4), allocatable :: head(:), next(:)

            N_particles = last_index - first_index + 1

            ! ── Cell size ~ 2R (collision diameter): only the 9 surrounding cells
            if (are_moons) then
                vrcut = dos * maxval(R_arr(first_index:last_index))
            else
                vrcut = dos * R_arr(first_index)
            end if

            vrcut = max(vrcut, sim%grid_col_min_cell_size)  ! Avoid too small cells (too much overhead)

            r_cut = vrcut * (uno + sim%verlet_skin_factor)  ! = 2R + skin
            r_cut2 = r_cut * r_cut

            ! Temporary storage — worst case N*(N-1)/2 pairs
            allocate(tmp_list(2, N_particles*(N_particles-1)/2))
            npairs = 0            

            ! ── Domain bounds from particle positions ───
            x_min = y(get_index(first_index))
            x_max = x_min
            y_min = y(get_index(first_index)+1)
            y_max = y_min

            do i = first_index, last_index
                idx = get_index(i)
                x_min = min(x_min, y(idx))
                x_max = max(x_max, y(idx))
                y_min = min(y_min, y(idx+1))
                y_max = max(y_max, y(idx+1))
            end do

            ! Add one cell of padding so boundary particles aren't clipped
            x_min = x_min - vrcut
            x_max = x_max + vrcut
            y_min = y_min - vrcut
            y_max = y_max + vrcut
    
            Ncx = max(1, ceiling((x_max - x_min) / vrcut))
            Ncy = max(1, ceiling((y_max - y_min) / vrcut))

            ! ── Safety: fall back to brute-force if grid is too large ──
            if (Ncx > sim%grid_col_max_cells / Ncy) then
                do i = first_index, last_index - 1
                    idx = get_index(i)
                    do j = i + 1, last_index
                        jdx = get_index(j)
                        dr_vec = y(jdx:jdx+1) - y(idx:idx+1)
                        dr2 = dr_vec(1)*dr_vec(1) + dr_vec(2)*dr_vec(2)
                        if (dr2 < r_cut2) then
                            npairs = npairs + 1
                            tmp_list(1, npairs) = i
                            tmp_list(2, npairs) = j
                        end if
                    end do
                end do
            else

                Nc = Ncx * Ncy
        
                ! ── Build linked-list cell structure ───
                allocate(head(0:Nc-1), next(first_index:last_index))
                head = -1  ! -1 = empty cell
        
                do i = first_index, last_index
                    idx = get_index(i)
                    cx = min(int((y(idx)   - x_min) / vrcut, int64), Ncx - 1_int64)
                    cy = min(int((y(idx+1) - y_min) / vrcut, int64), Ncy - 1_int64)
                    ! Clamp particles outside domain to boundary cells
                    cx = max(cx, 0_int64)
                    cy = max(cy, 0_int64)
                    ci = cx + Ncx * cy
                    next(i)  = head(ci)
                    head(ci) = i
                end do
        
                ! ── Check only 9-cell neighbourhood ───
                do i = first_index, last_index
                    idx = get_index(i)
                    cx = min(max(int((y(idx) - x_min) / vrcut, int64), 0_int64), Ncx - 1_int64)
                    cy = min(max(int((y(idx+1) - y_min) / vrcut, int64), 0_int64), Ncy - 1_int64)
        
                    do dcy = -1, 1
                        do dcx = -1, 1
                            cx2 = cx + dcx
                            cy2 = cy + dcy
                            if (cx2 < 0_int64 .or. cx2 >= Ncx) cycle
                            if (cy2 < 0_int64 .or. cy2 >= Ncy) cycle
            
                            ci2 = cx2 + Ncx * cy2
                            j = head(ci2)
                            do while (j /= -1)
                                if (j > i) then  ! avoid double-counting
                                    jdx = get_index(j)
                                    dr_vec = y(jdx:jdx+1) - y(idx:idx+1)
                                    dr2 = dr_vec(1)*dr_vec(1) + dr_vec(2)*dr_vec(2)
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

            ! Store compacted list
            if (allocated(vlist)) deallocate(vlist)
            allocate(vlist(2, npairs))
            vlist(:, 1:npairs) = tmp_list(:, 1:npairs)
            vlist_n = npairs
            deallocate(tmp_list)

            ! Store positions at build time (for drift check)
            if (allocated(vlist_pos)) deallocate(vlist_pos)
            if (allocated(vlist_r)) deallocate(vlist_r)
            allocate(vlist_pos(2, first_index:last_index))
            allocate(vlist_r(first_index:last_index))
            do i = first_index, last_index
                idx = get_index(i)
                vlist_pos(1, i) = y(idx)
                vlist_pos(2, i) = y(idx+1)
                vlist_r(i) = sqrt(y(idx)**2 + y(idx+1)**2)
            end do

            vlist_built = .True.
            vlist_time = time

            verlet_rebuilds = verlet_rebuilds + 1

        end subroutine build_verlet_list

        ! ── Verlet collision loop ───
        ! Checks if list needs rebuilding (any particle drifted > skin/2),
        ! then loops only over listed pairs.
        subroutine collisions_verlet(time, y, der, first_index, last_index, gamma_n, gamma_t, are_moons)
            implicit none
            real(wp), intent(in) :: time
            real(wp), dimension(:), intent(in)  :: y
            real(wp), dimension(:), intent(inout) :: der
            integer(kind=4), intent(in) :: first_index, last_index
            real(wp), intent(in) :: gamma_n, gamma_t
            logical, intent(in) :: are_moons

            integer(kind=4) :: i, j, k, idx 
            real(wp) :: r_skin2, dx, dy, drift2
            real(wp) :: dt
            real(wp), dimension(2) :: xyr

            ! ── Rebuild check ───
            ! Rebuild if: never built, size changed, or any particle drifted > skin/2
            if (.not. vlist_built) then
                call build_verlet_list(time, y, first_index, last_index, are_moons)
            else if (size(vlist_pos, 2) /= last_index - first_index + 1) then
                call build_verlet_list(time, y, first_index, last_index, are_moons)
            else
                r_skin2 = (vrcut * sim%verlet_skin_factor * uno2)**2  ! (skin/2)²

                dt = time - vlist_time
                
                do i = first_index, last_index
                    idx = get_index(i)
                    xyr = vlist_pos(:, i)

                    call get_xy_rotated(xyr, vlist_r(i), dt, G*m_arr(1), y(2))
                    dx = y(idx) - xyr(1)
                    dy = y(idx+1) - xyr(2)

                    drift2 = dx*dx + dy*dy
                    if (drift2 > r_skin2) then
                        call build_verlet_list(time, y, first_index, last_index, are_moons)
                        exit  ! rebuilt — no need to check further
                    end if
                end do
            end if

            ! ── Apply forces for all listed pairs ──
            do k = 1, vlist_n
                i = vlist(1, k)
                j = vlist(2, k)
                call soft_sphere_force(y, der, i, j, gamma_n, gamma_t, are_moons)
            end do

        end subroutine collisions_verlet

        pure subroutine get_xy_rotated_inertial(xy, r, dt, GM, dummy)
            implicit none
            real(wp), dimension(2), intent(inout) :: xy
            real(wp), intent(in) :: r, dt, GM, dummy

            real(wp) :: n
            real(wp) :: dphi
            
            n = sqrt(GM / r**3)
            dphi = n * dt

            xy = rotate2D(xy, dphi)

        end subroutine get_xy_rotated_inertial

        pure subroutine get_xy_rotated_sinodic(xy, r, dt, GM, omega)
            implicit none
            real(wp), dimension(2), intent(inout) :: xy
            real(wp), intent(in) :: r, dt, GM, omega

            real(wp) :: n
            real(wp) :: dphi
            
            n = sqrt(GM / r**3)
            dphi = (n - omega) * dt   ! subtract synodic rotation to get inertial frame angle change

            xy = rotate2D(xy, dphi)
            
        end subroutine get_xy_rotated_sinodic

end module collisions
