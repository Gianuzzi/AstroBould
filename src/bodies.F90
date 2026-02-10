!> Module with System, Asteroid, Moons, and Particles structures and routines
module bodies
    use constants, only: wp, cero, uno, uno2, uno3, dos, G, pi, twopi, myepsilon, sqepsilon, tini, infinito, &
                         & unit_mass, unit_time, unit_dist, unit_vel, unit_ener, unit_angm, radian
    use celestial, only: get_a_corot, get_acc_and_pot_single, get_acc_and_pot_single_with_J2, elem, coord, coord2geom
    use auxiliary, only: quickargsort, quickargsort_int, rotate2D

    implicit none

    type :: sphere_st
        integer(kind=4) :: id = -1  ! Identifier
        real(wp) :: mu_to_primary = cero ! Mass ratio to primary body [config]
        real(wp) :: mu_to_asteroid = cero ! Mass ratio to asteroid
        real(wp) :: mass = cero  ! Mass
        real(wp) :: radius = cero  ! Radius  [config]
        real(wp) :: initial_theta ! Initial angle from X (asteroid CM)
        real(wp) :: theta_from_primary ! Initial angle from primary(X)  [config]
        real(wp) :: dist_to_asteroid = cero  ! Distance to asteroid CM
        real(wp), dimension(4) :: coordinates_CM = cero  ! Coordinates from asteroid CM (non rotated)
    end type sphere_st

    type :: primary_st
        logical :: is_sphere = .True. ! Whether the primary is a sphere, or an ellipsoid
        integer(kind=4) :: id = -1  ! Identifier ! Will be 0
        real(wp) :: mu_to_asteroid = cero ! Mass ratio to asteroid
        real(wp) :: mass = cero  ! Mass
        real(wp) :: radius = cero  ! Radius  [config]
        real(wp) :: initial_theta ! Initial angle from X (asteroid CM)
        real(wp) :: dist_to_asteroid = cero  ! Distance to asteroid CM
        real(wp), dimension(4) :: coordinates_CM = cero  ! Coordinates from asteroid CM (non rotated)
        real(wp), dimension(3) :: semi_axis = cero  ! a, b, c
        real(wp) :: C20 = cero
        real(wp) :: C22 = cero
    end type primary_st

    type :: boulder_z_st
        logical :: active  ! Whether this is active or not
        real(wp) :: mu_to_primary = cero  ! Mass ratio to primary
        real(wp) :: mu_to_asteroid = cero ! Mass ratio to asteroid
        real(wp) :: mass = cero  ! Mass
        real(wp) :: dist_to_asteroid = cero  ! Distance in Z
    end type boulder_z_st

    type :: asteroid_st
        integer(kind=4) :: Nbodies = 0  ! Amount of boulders + primary
        integer(kind=4) :: Nboulders = 0  ! Amount of boulders
        type(primary_st) :: primary  ! Primary object
        type(sphere_st), allocatable :: boulders(:)  ! Without primary (from 1)
        type(boulder_z_st) :: boulder_z
        real(wp) :: mass = cero  ! Mass  [config -]
        real(wp) :: radius = cero  ! Radius
        real(wp) :: theta ! Rotational angle from X [dynamic]
        real(wp) :: omega = cero  ! Spin  [config] [dynamic]
        real(wp) :: a_corotation = cero  ! Corotation a, for massless
        real(wp) :: rotational_period = cero  ! Period of rotation  [config]
        real(wp) :: omega_kep = cero  ! Keplerian omega boulders would have
        real(wp) :: lambda_kep = cero  ! Ratio from omega to omega_kep  [config] [dynamic]
        real(wp) :: inertia = cero  ! Inertia moment
        real(wp) :: ang_mom_rot = cero  ! Rotational angular momentum [dynamic]
        real(wp) :: ang_mom_orb = cero  ! Orbital angular momentum [dynamic]
        real(wp) :: dist_to_cm = cero  ! Distance to origin [dynamic]
        real(wp), dimension(4) :: elements = cero  ! a, e, M, w
        real(wp), dimension(4) :: coordinates = cero  ! x, y, vx, vy
        real(wp) :: C20 = cero  ! C20 of the asteroid as a whole
        real(wp) :: e_rot = cero ! Rotational energy [dynamic]
        real(wp) :: e_kin = cero ! Kinetic energy [dynamic]
        real(wp) :: chaos_a(2) = (/infinito, cero/) ! (a_min, a_max)
        real(wp) :: chaos_e(2) = (/infinito, cero/) ! (e_min, e_max)
    end type asteroid_st

    type :: particle_st
        logical :: active  ! Whether this particle is active (orbiting) or not
        integer(kind=4) :: id = -1  ! Identifier
        real(wp) :: dist_to_cm = cero  ! Distance to origin
        real(wp), dimension(4) :: coordinates = cero  ! x, y vx, vy
        real(wp), dimension(4) :: elements = cero  ! a, e, M, w  [config]
        real(wp), dimension(4) :: geometric = cero  ! a_geom, e_geom, M_geom, w_geom [config]
        real(wp) :: mmr = cero  ! mean motion ratio to asteroid  [config]
        real(wp) :: mmr_geom = cero  ! geometric mean motion ratio to asteroid
        real(wp) :: chaos_a(2) = (/infinito, cero/) ! (a_min, a_max)
        real(wp) :: chaos_e(2) = (/infinito, cero/) ! (e_min, e_max)
        real(wp) :: chaos_a_geom(2) = (/infinito, cero/) ! (a_geom_min, a_geom_max)
        real(wp) :: chaos_e_geom(2) = (/infinito, cero/) ! (e_geom_min, e_geom_max)
        real(wp) :: variational(4) = cero ! dx, dy, dvx, dvy  ! For Megno
        real(wp) :: var_dist = cero ! variational distance
        real(wp) :: lyapunov = uno ! lyapunov value
        real(wp) :: megno = cero ! megno value
        real(wp) :: tmax = cero ! max time integrated
        integer(kind=4) :: merged_to = -1  ! IDX of body it was merged to
    end type particle_st

    type, extends(particle_st) :: moon_st
        real(wp) :: mu_to_asteroid = cero ! Mass ratio to asteroid  [config]
        real(wp) :: mass = cero  ! Masa
        real(wp) :: radius = cero  ! Radio
        real(wp) :: inertia = cero  ! Inertia moment
        real(wp) :: ang_mom_orb = cero  ! Orbital angular momentum [dynamic]
        real(wp) :: ang_mom_rot = cero  ! Rotational angular momentum [dynamic] (Just for conservation)
        real(wp) :: e_kin = cero ! Kinetic energy [dynamic]
        real(wp) :: e_rot = cero ! Rotational energy [dynamic] (Just for conservation)
    end type moon_st


    ! This struct will store the CM (system) properties and objects
    type :: system_st
        real(wp) :: time = cero
        real(wp) :: mass = cero
        real(wp) :: energy = cero
        real(wp) :: ang_mom = cero
        real(wp) :: inertia = cero
        real(wp) :: eta_col = uno
        real(wp) :: f_col = uno
        type(asteroid_st) :: asteroid
        integer(kind=4) :: Nmoons = 0
        integer(kind=4) :: Nmoons_active = 0
        type(moon_st), allocatable :: moons(:)
        integer(kind=4) :: Nparticles = 0
        integer(kind=4) :: Nparticles_active = 0
        type(particle_st), allocatable :: particles(:)
        real(wp) :: megno_epsilon = cero
        integer(kind=4) :: reference_frame = 0  ! Just to have an idea. 0:B, 1:J, 2:A
    end type system_st

    !!!!!!!!!!!   Used I/O formats   !!!!!!!!!
    character(30), parameter :: i3r23 = "(I7, I7, I7, 23(1X, 1PE22.15))"  ! 3 int y 23 real
    character(30), parameter :: i3r24 = "(I7, I7, I7, 24(1X, 1PE22.15))"  ! 3 int y 24 real
    character(26), parameter :: i2r15 = "(I7, I7, 15(1X, 1PE22.15))"  ! 2 int y 15 real
    character(26), parameter :: i2r16 = "(I7, I7, 16(1X, 1PE22.15))"  ! 2 int y 15 real
    character(25), parameter :: i2r9 = "(I7, I7, 9(1X, 1PE22.15))"  ! 3 int y 9 real
    character(18), parameter :: s1i5x5 = "(5(A, 1X, I5, 1X))"
    character(18), parameter :: r13 = "(13(1X, 1PE22.15))"

contains
    !  ------------------------   MEMORY    -------------------------------

    subroutine allocate_asteroid(self, Nboulders)
        implicit none
        type(asteroid_st), intent(inout) :: self
        integer(kind=4), intent(in) :: Nboulders

        if (allocated(self%boulders)) then
            write (*, *) "ERROR: Number of boulder already defined."
            STOP 1
        end if

        allocate (self%boulders(Nboulders))  ! FROM 1 to Nboulders
        self%Nbodies = Nboulders + 1  ! +1 to include primary
        self%Nboulders = Nboulders
    end subroutine allocate_asteroid

    subroutine allocate_moons(self, Nmoons)
        implicit none
        type(moon_st), allocatable, intent(inout) :: self(:)
        integer(kind=4), intent(in) :: Nmoons

        if (Nmoons > 0) then
            if (allocated(self)) then
                write (*, *) "ERROR: Number of moons already defined."
                write (*, *) "        Use moons file or config file; not both."
                STOP 1
            end if

            allocate (self(Nmoons))
        end if
    end subroutine allocate_moons

    subroutine allocate_particles(self, Nparticles)
        implicit none
        type(particle_st), allocatable, intent(inout) :: self(:)
        integer(kind=4), intent(in) :: Nparticles

        if (Nparticles > 0) then
            if (allocated(self)) then
                write (*, *) "ERROR: Number of particles already defined."
                write (*, *) "        Use particles file or config file; not both."
                STOP 1
            end if

            allocate (self(Nparticles))
        end if
    end subroutine allocate_particles

    pure subroutine free_asteroid(self)
        implicit none
        type(asteroid_st), intent(inout) :: self

        if (allocated(self%boulders)) deallocate (self%boulders)
    end subroutine free_asteroid

    pure subroutine free_system(self)
        implicit none
        type(system_st), intent(inout) :: self

        call free_asteroid(self%asteroid)
        if (allocated(self%moons)) deallocate (self%moons)
        if (allocated(self%particles)) deallocate (self%particles)
    end subroutine free_system

    !  -----------------------   ADD OBJECTS    ---------------------------

    subroutine add_primary(self, mass_primary_or_ast, semi_a_primary, semi_b_primary, semi_c_primary)
        implicit none
        type(asteroid_st), intent(inout) :: self
        real(wp), intent(in) :: mass_primary_or_ast, semi_a_primary, semi_b_primary, semi_c_primary
        real(wp) :: aux_radius = cero, radius_primary = cero

        if (self%primary%id /= -1) then
            write (*, *) "ERROR: Primary already set."
            stop 1
        end if

        ! Ensure not anti-radius
        aux_radius = semi_a_primary*semi_b_primary*semi_c_primary
        if (aux_radius < myepsilon) then
            write (*, *) "ERROR: Primary can not have zero or negative radius."
            stop 1
        end if
        radius_primary = (aux_radius)**uno3

        self%primary%is_sphere = abs(radius_primary - semi_c_primary) < myepsilon  ! c == radius => Sphere

        self%primary%semi_axis = (/semi_a_primary, semi_b_primary, semi_c_primary/)
        self%primary%mass = mass_primary_or_ast
        self%primary%radius = radius_primary
        self%primary%id = 0

        self%primary%C20 = (dos*semi_c_primary**2 - semi_a_primary**2 - semi_b_primary**2) &
                            & /(10.e0_wp*radius_primary**2)
        self%primary%C22 = (semi_a_primary**2 - semi_b_primary**2)/(20.e0_wp*radius_primary**2)

    end subroutine add_primary

    subroutine add_boulder_z(self, mu_boulder_z_to_primary)
        implicit none
        type(asteroid_st), intent(inout) :: self
        real(wp), intent(in) :: mu_boulder_z_to_primary
        
        self%boulder_z%mu_to_primary = mu_boulder_z_to_primary
        self%boulder_z%active = abs(mu_boulder_z_to_primary) > myepsilon
        
    end subroutine add_boulder_z

    subroutine add_boulder(self, mu_to_primary, radius, theta_from_primary)
        implicit none
        type(asteroid_st), intent(inout) :: self
        real(wp), intent(in) :: mu_to_primary, radius, theta_from_primary
        integer(kind=4) :: i
        logical :: slot_found

        slot_found = .False. ! Default

        ! Ensure not anti-radius
        if (radius < cero) then
            write (*, *) "ERROR: Boulders can not have negative radius. Use negative mass if needed."
            stop 1
        end if

        ! Look for first empty boulder slot (id = -1)
        do i = 1, self%Nboulders
            if (self%boulders(i)%id == -1) then
                self%boulders(i)%id = i
                self%boulders(i)%mu_to_primary = mu_to_primary
                self%boulders(i)%radius = radius
                self%boulders(i)%theta_from_primary = theta_from_primary
                slot_found = .True.
                exit
            end if
        end do

        if (.not. slot_found) then
            write (*, *) "ERROR: No available slot to add a new boulder."
            stop 1
        end if
    end subroutine add_boulder

    subroutine add_moon(self, mu_to_asteroid, ele_a, ele_e, ele_M, ele_w, mmr, radius, id_start)
        implicit none
        type(moon_st), dimension(:), intent(inout) :: self
        real(wp), intent(in) :: mu_to_asteroid, ele_a, ele_e, ele_M, ele_w, mmr, radius
        integer(kind=4), intent(in), optional :: id_start
        integer(kind=4) :: i
        integer(kind=4) :: id0 = 0  ! Where to start using IDs from (first is 1)
        logical :: slot_found

        slot_found = .False. ! Default

        ! Ensure not massless
        if (mu_to_asteroid < myepsilon) then
            write (*, *) "ERROR: Moons can not have zero or negative mass (I think)..."
            stop 1
        end if

        ! Ensure not anti-radius
        if (radius < cero) then
            write (*, *) "ERROR: Moons can not have negative radius."
            stop 1
        end if

        ! Set starting ID
        if (present(id_start)) id0 = id_start

        ! Look for first empty moon slot (id = -1)
        do i = 1, size(self)
            if (self(i)%id == -1) then
                self(i)%id = i + id0
                self(i)%active = .True.
                self(i)%mu_to_asteroid = mu_to_asteroid
                self(i)%elements(1) = ele_a
                self(i)%elements(2) = ele_e
                self(i)%elements(3) = ele_M
                self(i)%elements(4) = ele_w
                self(i)%mmr = mmr
                self(i)%radius = radius
                slot_found = .True.
                exit
            end if
        end do

        if (.not. slot_found) then
            write (*, *) "ERROR: No available slot to add a new moon."
            stop 1
        end if
    end subroutine add_moon

    subroutine add_particle(self, ele_a, ele_e, ele_M, ele_w, mmr, id_start)
        implicit none
        type(particle_st), dimension(:), intent(inout) :: self
        real(wp), intent(in) :: ele_a, ele_e, ele_M, ele_w, mmr
        integer(kind=4), intent(in), optional :: id_start
        integer(kind=4) :: i
        integer(kind=4) :: id0 = 0  ! Where to start using IDs from (first is 1)
        logical :: slot_found

        slot_found = .False. ! Default

        ! Set starting ID
        if (present(id_start)) id0 = id_start

        ! Look for first empty particle slot (id = -1)
        do i = 1, size(self)
            if (self(i)%id == -1) then
                self(i)%active = .True.
                self(i)%id = i + id0
                self(i)%elements(1) = ele_a
                self(i)%elements(2) = ele_e
                self(i)%elements(3) = ele_M
                self(i)%elements(4) = ele_w
                self(i)%mmr = mmr
                slot_found = .True.
                exit
            end if
        end do

        if (.not. slot_found) then
            write (*, *) "ERROR: No available slot to add a new particle."
            stop 1
        end if
    end subroutine add_particle

    !  -----------------------   INIT OBJECTS   ---------------------------

    ! Init asteroid parameters from config
    subroutine init_asteroid_params(self, lambda_kep, rotational_period)
        implicit none
        type(asteroid_st), intent(inout) :: self
        real(wp), intent(in) :: lambda_kep, rotational_period
        real(wp), allocatable :: coords_from_primary(:, :)
        real(wp) :: coords_cm_from_primary(4)
        real(wp) :: aux_real, aux_real4(4)
        integer(kind=4) :: i

        ! Ensure primary
        if (self%primary%id == -1) then
            write (*, *) "ERROR: Primary not loaded."
            stop 1
        end if

        ! Ensure all boulders
        do i = 1, self%Nboulders
            if (self%boulders(i)%id == -1) then
                write (*, *) "ERROR: Not all boulders loaded. Missing:", i
                stop 1
            end if
        end do

        ! Correct primary Mass
        if (self%primary%mass < cero) then
            aux_real = self%boulder_z%mu_to_primary + sum(self%boulders%mu_to_primary)
            self%primary%mass = abs(self%primary%mass)/(uno + aux_real)
        end if

        ! Set boulders Masses
        self%boulder_z%mass = self%boulder_z%mu_to_primary*self%primary%mass
        do i = 1, self%Nboulders
            self%boulders(i)%mass = self%boulders(i)%mu_to_primary*self%primary%mass
        end do

        ! Set asteroid Mass and boulders Ratios to primary
        self%mass = self%primary%mass + self%boulder_z%mass + sum(self%boulders%mass)  ! Mass
        if (self%mass <= cero) then
            write (*, *) "ERROR: Asteroid mass is zero or negative"
            stop 1
        end if

        self%primary%mu_to_asteroid = self%primary%mass/self%mass  ! Mass ratio of primary to asteroid
        self%boulder_z%mu_to_asteroid = self%boulder_z%mass/self%mass  ! Mass ratio of boulder z to asteroid
        self%boulders%mu_to_asteroid = self%boulders%mass/self%mass  ! Mass ratio of boulders to asteroid

        ! Set Rotations
        self%theta = cero  ! Initial angle with X axis
        self%omega_kep = sqrt(G*self%mass/self%primary%radius**3) ! Keplerian mean motion boulders would have
        if (abs(lambda_kep) > myepsilon) then
            self%lambda_kep = lambda_kep
            self%omega = self%omega_kep*self%lambda_kep  ! Angular velocity of primary
            self%rotational_period = twopi/abs(self%omega) ! Rotational period of primary
        else if (abs(rotational_period) > myepsilon) then
            self%rotational_period = abs(rotational_period)
            self%omega = (twopi/rotational_period) ! Angular velocity of primary
            self%lambda_kep = self%omega/self%omega_kep ! Ratio of omegas
        else  ! No rotation at all
            self%rotational_period = cero
            self%omega = cero ! Angular velocity of primary
            self%lambda_kep = cero ! Ratio of omegas
        end if
        self%a_corotation = get_a_corot(self%mass, self%omega)

        ! Calculate Coordinates and angle to CM from boulders
        coords_cm_from_primary = cero
        if (self%Nboulders > 0) then
            allocate (coords_from_primary(self%Nboulders, 4))
            coords_from_primary = cero
            do i = 1, self%Nboulders
                coords_from_primary(i, 1) = self%primary%radius*cos(self%boulders(i)%theta_from_primary) ! X from primary
                coords_from_primary(i, 2) = self%primary%radius*sin(self%boulders(i)%theta_from_primary) ! Y from primary
                coords_from_primary(i, 3) = -self%omega*coords_from_primary(i, 2)  ! vX
                coords_from_primary(i, 4) = self%omega*coords_from_primary(i, 1)  ! vY
            end do
            !! Get coords of CM from primary
            coords_cm_from_primary(1) = dot_product(coords_from_primary(:, 1), self%boulders(1:)%mu_to_asteroid)
            coords_cm_from_primary(2) = dot_product(coords_from_primary(:, 2), self%boulders(1:)%mu_to_asteroid)
            coords_cm_from_primary(3) = dot_product(coords_from_primary(:, 3), self%boulders(1:)%mu_to_asteroid)
            coords_cm_from_primary(4) = dot_product(coords_from_primary(:, 4), self%boulders(1:)%mu_to_asteroid)
            if (abs(coords_cm_from_primary(1)) < myepsilon) coords_cm_from_primary(1) = cero  ! Little correction
            if (abs(coords_cm_from_primary(2)) < myepsilon) coords_cm_from_primary(2) = cero  ! Little correction
            if (abs(coords_cm_from_primary(3)) < myepsilon) coords_cm_from_primary(3) = cero  ! Little correction
            if (abs(coords_cm_from_primary(4)) < myepsilon) coords_cm_from_primary(4) = cero  ! Little correction
        end if

        !! Set Distances and Angles
        !!! Primary
        self%primary%coordinates_CM = -coords_cm_from_primary
        self%primary%dist_to_asteroid = sqrt(coords_cm_from_primary(1)*coords_cm_from_primary(1) + &
                                               & coords_cm_from_primary(2)*coords_cm_from_primary(2))
        if ((self%Nboulders > 0) .and. (self%primary%dist_to_asteroid > myepsilon)) then
            self%primary%initial_theta = modulo(atan2(coords_cm_from_primary(2), coords_cm_from_primary(1)) + pi, twopi)
        else
            self%primary%initial_theta = cero
        end if
        !!! Boulder z
        self%boulder_z%dist_to_asteroid = sqrt(self%primary%radius**2 - self%primary%dist_to_asteroid**2)
        !!! Boulders
        do i = 1, self%Nboulders
            aux_real4 = coords_from_primary(i, :) - coords_cm_from_primary
            self%boulders(i)%coordinates_CM = aux_real4
            self%boulders(i)%initial_theta = atan2(aux_real4(2), aux_real4(1))
            self%boulders(i)%dist_to_asteroid = sqrt(aux_real4(1)*aux_real4(1) + aux_real4(2)*aux_real4(2))
        end do

        ! Set Coordinates
        self%coordinates = cero  ! We begin at the origin

        ! Set Radius
        if (self%primary%is_sphere) then
            self%radius = self%primary%radius
            do i = 1, self%Nboulders
                self%radius = max(self%radius, self%boulders(i)%radius + self%boulders(i)%dist_to_asteroid) ! Radius
            end do
        else  ! If ellipsoid, set maximum semi-axis
            self%radius = self%primary%semi_axis(1)  ! a
        end if

        ! Set Inertia
        !! Iz of Ellipsoid, or sphere if equal
        self%inertia = 0.2e0_wp*self%primary%mass*(self%primary%semi_axis(1)**2 + self%primary%semi_axis(2)**2)
        do i = 1, self%Nboulders
            !! Inertia Sphere boulder
            aux_real = 0.4e0_wp*self%boulders(i)%mass*self%boulders(i)%radius**2
            !! Sphere + Steiner
            self%inertia = self%inertia + aux_real + self%boulders(i)%mass*self%boulders(i)%dist_to_asteroid**2
        end do

        ! Set Angular momentum
        self%ang_mom_rot = self%inertia*self%omega ! Rotational
        !! Can not set ang_mom_orb here.

        ! Set Energy
        self%e_rot = uno2*self%inertia*self%omega**2  ! Rotational
        !! Can not set e_kin here. We are assuming v=0

        !! Free memory
        if (allocated(coords_from_primary)) deallocate (coords_from_primary)
    end subroutine init_asteroid_params

    ! Init moon derived parameters
    subroutine init_moons_params(self, asteroid)
        implicit none
        type(moon_st), dimension(:), intent(inout) :: self
        type(asteroid_st), intent(in) :: asteroid
        real(wp) :: combined_mass, aux_real
        real(wp) :: aux_real6(6)
        integer(kind=4) :: i

        do i = 1, size(self)
            if (self(i)%id == -1) then  ! Ensure all moons
                write (*, *) "ERROR: Not all moons loaded. Missing:", i
                stop 1
            end if
            self(i)%mass = self(i)%mu_to_asteroid*asteroid%mass

            ! Define coordinates
            !! Elements are asteroid-centric relative
            combined_mass = asteroid%mass + self(i)%mass
            aux_real = get_a_corot(combined_mass, asteroid%omega)  ! Ask. Including mass?
            if (self(i)%mmr > myepsilon) then
                if (aux_real < myepsilon) then  ! Ensure rotating
                    write (*, *) "ERROR: Can not set moon MMR with non rotating asteroid."
                    stop 1
                end if
                self(i)%elements(1) = self(i)%mmr**(2.e0_wp/3.e0_wp)*aux_real  ! a
            else
                if (aux_real < myepsilon) then  ! Ensure rotating
                    self(i)%mmr = cero
                else
                    self(i)%mmr = (self(i)%elements(1)/aux_real)**(1.5e0_wp) ! MMR
                end if
            end if
            call coord(combined_mass, &
                        & self(i)%elements(1), &
                        & self(i)%elements(2), &
                        & cero, & ! inc
                        & self(i)%elements(3), &
                        & self(i)%elements(4), &
                        & cero, & ! Omega
                        & aux_real6)  ! This are asteroid-centric
            self(i)%coordinates = (/aux_real6(1:2), aux_real6(4:5)/)

            ! Inertia
            self%inertia = 0.4e0_wp*self(i)%mass*self(i)%radius*self(i)%radius

            ! Angular momentum
            self(i)%ang_mom_rot = cero
            self(i)%ang_mom_orb = self(i)%mass*(self(i)%coordinates(1)*self(i)%coordinates(4) - &
                                                & self(i)%coordinates(2)*self(i)%coordinates(3))  ! Z component

            ! Energy
            self(i)%e_rot = cero
            self(i)%e_kin = uno2*self(i)%mass*(self(i)%coordinates(3)*self(i)%coordinates(3) + &
                                                 & self(i)%coordinates(4)*self(i)%coordinates(4))  ! Kinetic
        end do
    end subroutine init_moons_params

    ! Init particles derived parameters
    subroutine init_particles_params(self, asteroid)
        implicit none
        type(particle_st), dimension(:), intent(inout) :: self
        type(asteroid_st), intent(in) :: asteroid
        integer(kind=4) :: i
        real(wp) :: aux_real6(6)

        do i = 1, size(self)
            if (self(i)%id == -1) then  ! Ensure all moons
                write (*, *) "ERROR: Not all particles loaded. Missing:", i
                stop 1
            end if

            ! Define coordinates
            !! Elements are asteroid-centric
            if (self(i)%mmr > myepsilon) then
                if (asteroid%a_corotation < myepsilon) then  ! Ensure rotating
                    write (*, *) "ERROR: Can not set particle MMR with non rotating asteroid."
                    stop 1
                end if
                self(i)%elements(1) = self(i)%mmr**(2.e0_wp/3.e0_wp)*asteroid%a_corotation  ! a
            else
                if (asteroid%a_corotation < myepsilon) then  ! Ensure rotating
                    self(i)%mmr = cero
                else
                    self(i)%mmr = (self(i)%elements(1)/asteroid%a_corotation)**(1.5e0_wp) ! MMR
                end if
            end if
            call coord(asteroid%mass, &
                        & self(i)%elements(1), &
                        & self(i)%elements(2), &
                        & cero, &  ! inc
                        & self(i)%elements(3), &
                        & self(i)%elements(4), &
                        & cero, &  ! Omega
                        & aux_real6)  ! This are asteroid-centric
            self(i)%coordinates = (/aux_real6(1:2), aux_real6(4:5)/)
        end do
    end subroutine init_particles_params

    !  --------------------   CHANGE PARAMETERS    ------------------------

    ! Update asteroid parameters according to new rotation
    pure subroutine spin_asteroid(self, theta, omega)
        implicit none
        type(asteroid_st), intent(inout) :: self
        real(wp), intent(in) :: theta
        real(wp), intent(in) :: omega
        integer(kind=4) :: i

        ! Update Rotation
        self%theta = modulo(theta, twopi)  ! Update THETA
        if (abs(omega) > myepsilon) then
            self%a_corotation = get_a_corot(self%mass, omega)  ! This is a bit odd now...
            self%rotational_period = twopi/omega ! Rotational period of primary
            self%lambda_kep = omega/self%omega_kep ! Ratio of omegas
            self%omega = omega  ! Update OMEGA
        else
            self%a_corotation = cero
            self%rotational_period = cero ! Rotational period of primary
            self%lambda_kep = cero  ! Ratio of omegas
            self%omega = cero ! Update OMEGA
        end if

        ! Update primary coordinates
        self%primary%coordinates_CM(1) = self%primary%dist_to_asteroid* &
                                         & cos(self%theta + self%primary%initial_theta)  ! x
        self%primary%coordinates_CM(2) = self%primary%dist_to_asteroid* &
                                         & sin(self%theta + self%primary%initial_theta)  ! y
        self%primary%coordinates_CM(3) = -self%omega*self%primary%coordinates_CM(2)  ! vx
        self%primary%coordinates_CM(4) = self%omega*self%primary%coordinates_CM(1)  ! vy

        ! Update boulders coordinates
        do i = 1, self%Nboulders
            self%boulders(i)%coordinates_CM(1) = self%boulders(i)%dist_to_asteroid* &
                                                 & cos(self%theta + self%boulders(i)%initial_theta)  ! x
            self%boulders(i)%coordinates_CM(2) = self%boulders(i)%dist_to_asteroid* &
                                                 & sin(self%theta + self%boulders(i)%initial_theta)  ! y
            self%boulders(i)%coordinates_CM(3) = -self%omega*self%boulders(i)%coordinates_CM(2)  ! vx
            self%boulders(i)%coordinates_CM(4) = self%omega*self%boulders(i)%coordinates_CM(1)  ! vy
        end do

        ! Update Angular momentum
        self%ang_mom_rot = self%inertia*omega ! Rotational

        ! Update Energy
        self%e_rot = uno2*self%inertia*omega**2  ! Rotational
    end subroutine spin_asteroid

    ! Update asteroid parameters according to new coordinates
    pure subroutine shift_asteroid(self, coordinates)
        implicit none
        type(asteroid_st), intent(inout) :: self
        real(wp), intent(in) :: coordinates(4)  ! Barycentric

        ! Update Coords and distance
        self%coordinates = coordinates ! Update COORDINATES
        self%dist_to_cm = sqrt(coordinates(1)*coordinates(1) + &
                             & coordinates(2)*coordinates(2))

        ! Update Angular momentum
        self%ang_mom_orb = self%mass*(coordinates(1)*coordinates(4) - &
                                      & coordinates(2)*coordinates(3))  ! Z component

        ! Update Energy
        self%e_kin = uno2*self%mass*(coordinates(3)*coordinates(3) + &
                                       & coordinates(4)*coordinates(4))  ! Kinetic
    end subroutine shift_asteroid

    ! Update asteroid parameters according to new coordinates
    pure subroutine grow_asteroid(self, mass_to_add)
        implicit none
        type(asteroid_st), intent(inout) :: self
        real(wp), intent(in) :: mass_to_add
        real(wp) :: growth
        integer(kind=4) :: i

        if (mass_to_add < myepsilon) return  ! No mass to add

        ! Get mass growth ratio
        growth = uno + (mass_to_add/self%mass)

        !! Update primary mass [mass ratios are kept the same]
        self%primary%mass = self%primary%mass*growth
        !! Update boulder z masses [mass ratios are kept the same]
        self%boulder_z%mass = self%boulder_z%mass*growth
        !! Update boulders masses [mass ratios are kept the same]
        do i = 1, self%Nboulders
            self%boulders(i)%mass = self%boulders(i)%mass*growth
        end do
        self%mass = self%mass*growth  ! Update MASS

        ! Update derived rotation parameters
        self%omega_kep = sqrt(G*self%mass/self%primary%radius**3)  ! Check
        self%a_corotation = get_a_corot(self%mass, self%omega)  ! Corotation a
        if (abs(self%omega) > cero) self%lambda_kep = self%omega/self%omega_kep ! Ratio of omegas

        ! Update Inertia
        self%inertia = self%inertia*growth

        ! Update angular momentum
        self%ang_mom_rot = self%ang_mom_rot*growth ! Rotational
        self%ang_mom_orb = self%ang_mom_orb*growth ! Orbital

        ! Update Energy
        self%e_rot = self%e_rot*growth  ! Rotational
        self%e_kin = self%e_kin*growth  ! Kinetic
    end subroutine grow_asteroid

    ! Update a single moon parameters according to new coordinates
    pure subroutine shift_single_moon(self, coordinates)
        implicit none
        type(moon_st), intent(inout) :: self
        real(wp), intent(in) :: coordinates(4)  ! These are barycentric now

        ! Update coordinates
        self%coordinates = coordinates  ! Update COORDINATES

        !! Distance to CM
        self%dist_to_cm = sqrt(coordinates(1)*coordinates(1) + &
                             & coordinates(2)*coordinates(2))

        ! Angular momentum
        self%ang_mom_orb = self%mass*(self%coordinates(1)*self%coordinates(4) - &
                                      & self%coordinates(2)*self%coordinates(3))  ! Z component

        ! Energy
        self%e_kin = uno2*self%mass*(coordinates(3)*coordinates(3) + &
                                       & coordinates(4)*coordinates(4))  ! Kinetic
    end subroutine shift_single_moon

    ! Update a single particle parameters according to new coordinates
    pure subroutine shift_single_particle(self, coordinates)
        implicit none
        type(particle_st), intent(inout) :: self
        real(wp), intent(in) :: coordinates(4)  ! These are barycentric now

        ! Update coordinates
        self%coordinates = coordinates  ! Update COORDINATES

        !! Distance to CM
        self%dist_to_cm = sqrt(coordinates(1)*coordinates(1) + &
                             & coordinates(2)*coordinates(2))
    end subroutine shift_single_particle

    ! Center system
    pure subroutine center_sytem(self)
        implicit none
        type(system_st), intent(inout) :: self
        real(wp) :: mcm, rv_cm(4)
        integer(kind=4) :: i

        ! Get cm
        call get_cm(self, mcm, rv_cm)

        ! Shif asteroid
        call shift_asteroid(self%asteroid, self%asteroid%coordinates - rv_cm)

        ! Shift moons
        do i = 1, self%Nmoons_active
            call shift_single_moon(self%moons(i), self%moons(i)%coordinates - rv_cm)
        end do

        ! Shift particles
        do i = 1, self%Nparticles_active
            call shift_single_particle(self%particles(i), self%particles(i)%coordinates - rv_cm)
        end do

    end subroutine center_sytem

    !  ---------------------   GET PARAMETERS    --------------------------

    ! (Re)Calculate system mass and CM
    pure subroutine get_cm(self, mass, coordinates)
        implicit none
        type(system_st), intent(in) :: self
        real(wp), dimension(4), intent(out) :: coordinates
        real(wp), intent(out) :: mass
        real(wp) :: aux_4(4)
        integer(kind=4) :: i

        ! Calculate cm with all contributing parts
        mass = self%asteroid%primary%mass
        coordinates = (self%asteroid%primary%coordinates_CM + self%asteroid%coordinates)*self%asteroid%primary%mass

        mass = mass + self%asteroid%boulder_z%mass
        coordinates = coordinates + self%asteroid%coordinates*self%asteroid%boulder_z%mass  ! Assume they are on the asteroid CM

        do i = 1, self%asteroid%Nboulders
            aux_4 = self%asteroid%boulders(i)%coordinates_CM + self%asteroid%coordinates
            coordinates = coordinates + aux_4*self%asteroid%boulders(i)%mass
            mass = mass + self%asteroid%boulders(i)%mass
        end do

        do i = 1, self%Nmoons_active
            coordinates = coordinates + self%moons(i)%coordinates*self%moons(i)%mass
            mass = mass + self%moons(i)%mass
        end do

        coordinates = coordinates/mass  ! rcm = sum(ri * mi) / mcm

    end subroutine get_cm

    ! Calculate system energy and angular momentum
    pure subroutine calculate_energy_ang_mom_inertia(self, energy, ang_mom, inertia)
        implicit none
        type(system_st), intent(in) :: self
        real(wp), intent(out) :: energy, ang_mom, inertia
        real(wp) :: e_pot
        integer(kind=4) :: i, j
        real(wp) :: dr(2), dist

        ! Asteroid energy and ang_mom are obtaines considering a rigid body

        ! Calculate energy and ang_mom
        ang_mom = self%asteroid%ang_mom_rot + self%asteroid%ang_mom_orb
        energy = self%asteroid%e_kin + self%asteroid%e_rot ! e_pot below
        inertia = self%asteroid%inertia

        !! e_pot
        e_pot = cero
        !! Asteroid (with moons)
        do i = 1, self%Nmoons_active
            dr = self%moons(i)%coordinates(1:2) - self%asteroid%coordinates(1:2)
            dist = sqrt(dr(1)*dr(1) + dr(2)*dr(2))
            inertia = inertia + self%moons(i)%mass*dist*dist
            if (dist < myepsilon) cycle
            e_pot = e_pot - (self%asteroid%mass*self%moons(i)%mass)/dist
        end do
        !! Moons (here, only moons with moons are added)
        do i = 1, self%Nmoons_active - 1
            do j = 2, self%Nmoons_active
                dr = self%moons(j)%coordinates(1:2) - self%moons(i)%coordinates(1:2)
                dist = sqrt(dr(1)*dr(1) + dr(2)*dr(2))
                if (dist < myepsilon) cycle
                e_pot = e_pot - (self%moons(i)%mass*self%moons(j)%mass)/dist
            end do
            !! Energy and Ang Mom (TBD: Last moon)
            energy = energy + self%moons(i)%e_kin + self%moons(i)%e_rot  ! Energy
            ang_mom = ang_mom + self%moons(i)%ang_mom_orb + self%moons(i)%ang_mom_rot  ! Angular Momentum
        end do

        energy = energy + e_pot*G  ! G

        !! Now the last Moon
        if (self%Nmoons_active > 0) then
            ang_mom = ang_mom + self%moons(self%Nmoons_active)%ang_mom_orb + self%moons(self%Nmoons_active)%ang_mom_rot
            energy = energy + self%moons(self%Nmoons_active)%e_kin + self%moons(self%Nmoons_active)%e_rot
        end if
    end subroutine calculate_energy_ang_mom_inertia

    ! Get amount of active bodies, including asteroid
    pure subroutine get_Nactive(self, Nactive)
        implicit none
        type(system_st), intent(in) :: self
        integer(kind=4), intent(inout) :: Nactive
        Nactive = 1 + self%Nmoons_active + self%Nparticles_active  ! Includes asteroid
    end subroutine get_Nactive

    ! Get amount values in array to generate (includes possible megno)
    pure function get_y_nvalues(self, new_Nactive) result(y_nvalues)
        implicit none
        type(system_st), intent(in) :: self
        integer(kind=4), intent(in), optional :: new_Nactive
        integer(kind=4) :: y_nvalues
        integer(kind=4) :: Nactive

        if (present(new_Nactive)) then
            Nactive = new_Nactive
        else
            call get_Nactive(self, Nactive)
        end if

        y_nvalues = 2 + 4*Nactive
        
        if (self%megno_epsilon > tini) y_nvalues = y_nvalues + 4*self%Nparticles_active
    end function get_y_nvalues

    ! Get moon by ID
    pure subroutine get_moon_by_ID(self, id, moon)
        implicit none
        type(system_st), intent(in) :: self
        integer(kind=4), intent(in) :: id
        type(moon_st), intent(out) :: moon
        integer(kind=4) :: i

        do i = 1, self%Nmoons
            if (self%moons(i)%id == id) then
                moon = self%moons(i)
                exit
            end if
        end do
    end subroutine get_moon_by_ID

    ! Get particle by ID
    pure subroutine get_particle_by_ID(self, id, particle)
        implicit none
        type(system_st), intent(in) :: self
        integer(kind=4), intent(in) :: id
        type(particle_st), intent(out) :: particle
        integer(kind=4) :: i

        do i = 1, self%Nparticles
            if (self%particles(i)%id == id) then
                particle = self%particles(i)
                exit
            end if
        end do
    end subroutine get_particle_by_ID

    !  -----------------   UPDATE INTERNAL PARAMETERS    ------------------

    ! Recalculate elements
    pure subroutine update_elements(self, reference_frame)
        implicit none
        type(system_st), intent(inout) :: self
        integer(kind=4), intent(in), optional :: reference_frame
        integer(kind=4) :: j, ref_used
        real(wp) :: a, e, i, M, w, O, a_corot
        real(wp) :: coords_shifted(4), coords_cm_in(4), mass_cm_in
        logical :: mmr_from_n

        ! Set reference frame
        if (present(reference_frame)) then 
            ref_used = reference_frame
        else
            ref_used = self%reference_frame
        end if

        ! This may not have too much of a meaning to be fair, depending on the reference frame used
        a_corot = get_a_corot(self%asteroid%mass, self%asteroid%omega)  ! Astrocentric here
        ! a_corot = get_a_corot(self%mass, self%ang_mom/self%inertia)  ! Something strange
        self%asteroid%a_corotation = a_corot
        mmr_from_n = (a_corot < myepsilon) .or. (self%Nmoons_active > 0)

        if ((ref_used .eq. 0) .or. (ref_used .eq. 1)) then  ! Baryc or Jacobi

            coords_shifted = cero
            coords_cm_in = self%asteroid%coordinates
            mass_cm_in = self%asteroid%mass
            
            ! ! Deprecated. Not used in reality
            ! if ((self%asteroid%dist_to_cm > myepsilon) .and. (self%Nmoons_active > 0)) then
            !     call elem(mass_cm_in, (/self%asteroid%coordinates(1:2), cero, self%asteroid%coordinates(3:4), cero/), &
            !             & a, e, i, M, w, O)
            !     self%asteroid%elements = (/a, e, M, w/)
            ! end if

            ! Assume mass ordered
            !!call sort_moons(self%moons)  ! If not ordered
            do j = 1, self%Nmoons_active
                coords_shifted = self%moons(j)%coordinates - coords_cm_in
                call elem(mass_cm_in + self%moons(j)%mass, &
                    & (/coords_shifted(1:2), cero, coords_shifted(3:4), cero/), &
                    a, e, i, M, w, O)
                self%moons(j)%elements = (/a, e, M, w/)
                self%moons(j)%mmr = self%asteroid%omega / (sqrt(G * (mass_cm_in + self%moons(j)%mass) / a**3)) ! MMR

                ! Update cm in
                coords_cm_in = mass_cm_in*coords_cm_in + self%moons(j)%mass*self%moons(j)%coordinates
                if (ref_used .eq. 0) then  ! If baryc
                    self%moons(j)%elements(1) = a * mass_cm_in/(mass_cm_in + self%moons(j)%mass)
                end if
                mass_cm_in = mass_cm_in + self%moons(j)%mass
                coords_cm_in = coords_cm_in/mass_cm_in
            end do
            
            ! Here, the cm is already defined
            do j = 1, self%Nparticles_active
                call elem(mass_cm_in, &
                    & (/self%particles(j)%coordinates(1:2), cero, &
                    & self%particles(j)%coordinates(3:4), cero/), &
                    a, e, i, M, w, O)
                self%particles(j)%elements = (/a, e, M, w/)
                if (mmr_from_n) then
                    self%particles(j)%mmr = self%asteroid%omega / (sqrt(G * mass_cm_in / a**3)) ! MMR
                else
                    self%particles(j)%mmr = (a/a_corot)**(1.5e0_wp) ! MMR
                end if
            end do

        else  ! Astrocentric

            coords_shifted = cero
            coords_cm_in = self%asteroid%coordinates
            mass_cm_in = self%asteroid%mass

            ! Moons
            do j = 1, self%Nmoons_active
                coords_shifted = self%moons(j)%coordinates - coords_cm_in
                call elem(mass_cm_in, &
                        & (/coords_shifted(1:2), cero, coords_shifted(3:4), cero/), &
                        a, e, i, M, w, O)
                self%moons(j)%elements = (/a, e, M, w/)
                self%moons(j)%mmr = self%asteroid%omega / (sqrt(G * (mass_cm_in + self%moons(j)%mass) / a**3)) ! MMR
            end do

            ! Particles
            do j = 1, self%Nparticles_active
                coords_shifted = self%particles(j)%coordinates - coords_cm_in
                call elem(mass_cm_in, &
                        & (/coords_shifted(1:2), cero, coords_shifted(3:4), cero/), &
                        a, e, i, M, w, O)
                self%particles(j)%elements = (/a, e, M, w/)
                if (mmr_from_n) then
                    self%particles(j)%mmr = self%asteroid%omega / (sqrt(G * mass_cm_in / a**3)) ! MMR
                else
                    self%particles(j)%mmr = (a/a_corot)**(1.5e0_wp) ! MMR
                end if
            end do
        end if
    end subroutine update_elements

    ! Recalculate geometric elements
    pure subroutine update_geometric(self)
        implicit none
        type(system_st), intent(inout) :: self
        integer(kind=4) :: j
        real(wp) :: a, e, i, M, w, O
        real(wp) :: n
        real(wp) :: coords_ast(4), ast_cm(4)
        real(wp) :: mass_ast, radius_ast, omega_ast
        real(wp) :: J2_ast

        coords_ast = cero
        mass_ast = self%asteroid%mass
        radius_ast = self%asteroid%radius
        omega_ast = self%asteroid%omega
        ast_cm = self%asteroid%coordinates
        J2_ast = -self%asteroid%C20

        ! Moons
        do j = 1, self%Nmoons_active
            coords_ast = self%moons(j)%coordinates - ast_cm
            call coord2geom(mass_ast + self%moons(j)%mass, radius_ast, J2_ast, &
                    & (/coords_ast(1:2), cero, coords_ast(3:4), cero/), &
                    a, e, i, M, w, O, n)
            M = modulo(M, twopi)
            w = modulo(w, twopi)
            self%moons(j)%geometric = (/a, e, M, w/)
            self%moons(j)%mmr_geom = omega_ast/n
        end do

        ! Particles
        do j = 1, self%Nparticles_active
            coords_ast = self%particles(j)%coordinates - ast_cm
            call coord2geom(mass_ast, radius_ast, J2_ast, &
                    & (/coords_ast(1:2), cero, coords_ast(3:4), cero/), &
                    a, e, i, M, w, O, n)
            M = modulo(M, twopi)
            w = modulo(w, twopi)
            self%particles(j)%geometric = (/a, e, M, w/)
            self%particles(j)%mmr_geom = omega_ast/n
        end do
    end subroutine update_geometric

    ! (Re)calculate all system main parameters
    pure subroutine recalculate_all(self)
        implicit none
        type(system_st), intent(inout) :: self
        real(wp) :: total_mass, energy, ang_mom, inertia
        real(wp) :: rvcm(4)
        integer(kind=4) :: i

        ! Get mass and then CM too
        call get_cm(self, total_mass, rvcm)

        self%mass = total_mass  ! Update MASS

        ! Set system at CM, shifting all objects.
        call shift_asteroid(self%asteroid, self%asteroid%coordinates - rvcm)
        do i = 1, self%Nmoons_active
            call shift_single_moon(self%moons(i), self%moons(i)%coordinates - rvcm)
        end do
        do i = 1, self%Nparticles_active
            call shift_single_particle(self%particles(i), self%particles(i)%coordinates - rvcm)
        end do

        ! Set Energy and Angular Momentum
        call calculate_energy_ang_mom_inertia(self, energy, ang_mom, inertia)
        self%ang_mom = ang_mom  !Set ANG_MOM
        self%energy = energy  !Set ENERGY
        self%inertia = inertia !Set INERTIA
    end subroutine recalculate_all

    !  -------------------   INIT MAIN SYSTEM   ---------------------------

    ! Init whole system
    subroutine init_system(self, asteroid, moons, particles, lambda_kep, rotational_period, reference_frame)
        implicit none
        type(system_st), intent(inout) :: self
        type(asteroid_st), intent(inout) :: asteroid
        type(moon_st), allocatable, intent(inout) :: moons(:)
        type(particle_st), allocatable, intent(inout) :: particles(:)
        real(wp), intent(in) :: lambda_kep, rotational_period
        integer(kind=4), intent(in) :: reference_frame
        integer(kind=4) :: i, j
        integer(kind=4), allocatable :: id_list(:)

        ! Set main reference frame
        self%reference_frame = reference_frame

        ! Initialize the objects and Create the system
        !! Asteroid
        call init_asteroid_params(asteroid, lambda_kep, rotational_period)
        self%asteroid = asteroid  ! CREATE ASTEROID
        !! Moons
        if (allocated(moons)) then
            call init_moons_params(moons, asteroid)
            self%Nmoons = size(moons)
            self%Nmoons_active = size(moons)
            self%moons = moons  ! CREATE MOONS
        else
            self%Nmoons = 0
            self%Nmoons_active = 0
        end if
        !! Particles
        if (allocated(particles)) then
            call init_particles_params(particles, asteroid)
            self%Nparticles = size(particles)
            self%Nparticles_active = size(particles)
            self%particles = particles  ! CREATE PARTICLES
        else
            self%Nparticles = 0
            self%Nparticles_active = 0
        end if

        ! Check unique ids
        allocate (id_list(self%Nmoons + self%Nparticles))
        do i = 1, self%Nmoons
            id_list(i) = self%moons(i)%id
            do j = i - 1, 1, -1  ! Back loop
                if (id_list(i) == id_list(j)) then
                    write (*, *) "WARNING: Moons IDs are not unique. Repeated:", id_list(i)
                    return
                end if
            end do
        end do
        do i = self%Nmoons + 1, self%Nmoons + self%Nparticles
            id_list(i) = self%particles(i - self%Nmoons)%id
            do j = i - 1, 1, -1  ! Back loop
                if (id_list(i) == id_list(j)) then
                    write (*, *) "WARNING: Particles + Moons IDs are not unique. Repeated:", id_list(i)
                    return
                end if
            end do
        end do

        ! Right now, all coordinates are AsteroidCentric
        call recalculate_all(self)

        ! Now, they are barycentric
        !! Sanity check for low values in case of no moons
        if ((self%Nmoons == 0) .and. self%asteroid%dist_to_cm > cero) then
            ! Manual shift
            self%asteroid%coordinates = cero
            self%asteroid%dist_to_cm = cero
            self%asteroid%ang_mom_orb = cero
            self%asteroid%e_kin = cero
        end if

    end subroutine init_system

    ! Set extra system parameters
    pure subroutine set_system_extra(self, time, eta_collision, f_collision, manual_J2, J2_from_primary)
        implicit none
        type(system_st), intent(inout) :: self
        real(wp), intent(in) :: time
        real(wp), intent(in) :: eta_collision, f_collision, manual_J2
        logical, intent(in) :: J2_from_primary

        ! Initial time. Set TIME
        self%time = time

        ! Initial collisional eta and f. Set eta_col and f_col
        self%eta_col = eta_collision
        self%f_col = f_collision
        if ((abs(manual_J2) > cero) .and. self%asteroid%primary%is_sphere) then
            if ((.not. J2_from_primary) .or. (self%asteroid%Nboulders .eq. 0)) then
                self%asteroid%C20 = -manual_J2  ! C20 = -J2
            else
                self%asteroid%primary%C20 = -manual_J2  ! C20 = -J2
            end if
        end if

    end subroutine set_system_extra

    ! ---------------------------   ADD MEGNO    ---------------------------

    subroutine init_megno(self, eps)
        implicit none
        type(system_st), intent(inout) :: self
        real(wp), intent(in) :: eps
        real(wp), dimension(4) :: rndm_arr
        real(wp) :: norma
        integer(kind=4) :: i

        if (self%Nparticles_active > 0) then

            do i = 1, self%Nparticles
                rndm_arr = uno
                call random_number(rndm_arr)
                rndm_arr = rndm_arr - 0.5
                norma = sqrt(sum(rndm_arr*rndm_arr))
                rndm_arr = rndm_arr / norma
                self%particles(i)%variational = rndm_arr * eps
                self%particles(i)%var_dist = eps
            end do

            self%megno_epsilon = eps

        else

            self%megno_epsilon = cero

        end if
        
    end subroutine init_megno

    !  ------------------------   OBJECT SWAP  -----------------------------

    ! Swap between two moons, using the positional index
    subroutine swap_moons(self, i, j)
        implicit none
        type(moon_st), intent(inout) :: self(:)
        integer(kind=4), intent(in) :: i, j
        type(moon_st) :: tmp_moon

        ! Safety check
        if (i < 1 .or. j < 1 .or. i > size(self) .or. j > size(self)) then
            write (*, *) "ERROR: indices out of range in swap_moons"
            stop 1
        end if

        ! Perform swap
        tmp_moon = self(i)
        self(i) = self(j)
        self(j) = tmp_moon
    end subroutine swap_moons

    ! Swap between two particles, using the positional index
    subroutine swap_particles(self, i, j)
        implicit none
        type(particle_st), intent(inout) :: self(:)
        integer(kind=4), intent(in) :: i, j
        type(particle_st) :: tmp_particles

        ! Safety check
        if (i < 1 .or. j < 1 .or. i > size(self) .or. j > size(self)) then
            write (*, *) "ERROR: indices out of range in swap_particles"
            stop 1
        end if

        ! Perform swap
        tmp_particles = self(i)
        self(i) = self(j)
        self(j) = tmp_particles
    end subroutine swap_particles

    ! Sort moons according to their mass
    pure subroutine sort_moons(self)
        implicit none    
        type(moon_st), intent(inout) :: self(:)
        integer(kind=4), dimension(size(self)) :: order
        type(moon_st), dimension(size(self)) :: tmp_moon
        integer(kind=4) :: i, Nmoons

        Nmoons = size(self)
        
        ! Get order
        call quickargsort(self%mass, order, 1, Nmoons)

        ! Get a copy
        tmp_moon = self

        ! Reorder
        do i = 1, Nmoons
            self(i) = tmp_moon(order(i))
        end do
    end subroutine sort_moons

    !  -------------------------   CHECK  ---------------------------------

    ! Check if syst is in CM = 0
    subroutine check_coordinates(self, error)
        implicit none
        type(system_st), intent(in) :: self
        integer(kind=4), intent(inout), optional :: error
        real(wp), dimension(4) :: coordinates
        real(wp) :: mass
        real(wp) :: aux_real

        ! Calculate cm
        call get_cm(self, mass, coordinates)

        ! Check MASS
        if (abs(mass - self%mass)/mass > myepsilon) then
            write (*, *) "WARNING: Total mass differs from CM mass. Error:", abs(mass - self%mass)/mass
            if (present(error)) error = error + 1
        end if

        ! Check POS
        aux_real = sqrt(coordinates(1)*coordinates(1) + &
                      & coordinates(2)*coordinates(2))
        if (aux_real > myepsilon) then
            !write (*,*) "WARNING: CM is not centered at 0.", aux_real
            if (present(error)) error = error + 1
        end if

        ! Check VEL
        aux_real = sqrt(coordinates(3)*coordinates(3) + &
                      & coordinates(4)*coordinates(4))
        if (aux_real > myepsilon) then
            !write (*,*) "WARNING: CM velocity is not 0.", aux_real
            if (present(error)) error = error + 1
        end if
    end subroutine check_coordinates

    !  -------------------   OBJECT DEACTIVATION  -------------------------

    ! Remove a moon from a system by deactivating and moving it to bottom
    subroutine deactivate_moon_i(self, i, keep_order, error)
        implicit none
        type(system_st), intent(inout) :: self
        integer(kind=4), intent(in) :: i
        logical, intent(in) :: keep_order
        type(moon_st), allocatable :: tmp_moons(:)  ! In case of keep order
        integer(kind=4), intent(inout), optional :: error
        integer(kind=4) :: j

        ! Check
        if (.not. self%moons(i)%active) then
            write (*, *) "WARNING: Can not deactivate moon already deactivated:", i
            if (present(error)) error = error + 1
            return
        end if

        ! Deactivate moon
        self%moons(i)%active = .False.

        ! Set tmax
        self%moons(i)%tmax = self%time

        ! If it is the last active, do not swap
        if (i < self%Nmoons_active) then
            if (keep_order) then
                allocate(tmp_moons(self%Nmoons_active))
                tmp_moons = self%moons(:self%Nmoons_active)
                do j = 1, i-1
                    self%moons(j) = tmp_moons(j)
                end do
                do j = i+1, self%Nmoons_active
                    self%moons(j-1) = tmp_moons(j)
                end do
                deallocate(tmp_moons)
            else
                ! Switch deactivated moon with last active
                do j = self%Nmoons_active, 1, -1
                    if (i == j) cycle
                    if (self%moons(j)%active) then
                        call swap_moons(self%moons, i, j)
                        exit
                    end if
                end do
            end if
        end if
        self%Nmoons_active = self%Nmoons_active - 1  ! Update NMOONS_ACTIVE
    end subroutine deactivate_moon_i

    ! Remove a particle from a system by deactivating and moving it to bottom
    subroutine deactivate_particle_i(self, i, error)
        implicit none
        type(system_st), intent(inout) :: self
        integer(kind=4), intent(in) :: i
        integer(kind=4), intent(inout), optional :: error
        integer(kind=4) :: j

        ! Check
        if (.not. self%particles(i)%active) then
            write (*, *) "WARNING: Can not deactivate particle already deactivated:", i
            if (present(error)) error = error + 1
        end if

        ! Deactivate particle
        self%particles(i)%active = .False.

        ! Set tmax
        self%particles(i)%tmax = self%time

        ! Switch deactivated particle with last active
        if (i < self%Nparticles_active) then
            do j = self%Nparticles_active, 1, -1
                if (i == j) cycle
                if (self%particles(j)%active) then
                    call swap_particles(self%particles, i, j)
                    exit
                end if
            end do
        end if
        self%Nparticles_active = self%Nparticles_active - 1  ! Update NPARTICLES_ACTIVE
    end subroutine deactivate_particle_i

    !  -------------------   UPDATE / GENERATION  -------------------------

    pure subroutine update_megno(self, der, rescale)
        implicit none
        type(system_st), intent(inout) :: self
        real(wp), dimension(:), intent(in) :: der
        logical, intent(in) :: rescale
        real(wp) :: dist2, prod2
        integer(kind=4) :: i, vdx

        ! Check if something to do
        if (self%megno_epsilon < tini) return
        
        vdx = 2 + (self%Nmoons_active + self%Nparticles_active) * 4 + 1
        do i = 1, self%Nparticles_active
            
            dist2 = max(sum(self%particles(i)%variational*self%particles(i)%variational), tini)

            prod2 = sum(self%particles(i)%variational*der(vdx:vdx+3))

            self%particles(i)%lyapunov = prod2 / dist2
            self%particles(i)%megno = self%particles(i)%lyapunov * self%time

            if (rescale) self%particles(i)%variational = self%particles(i)%variational * self%megno_epsilon / sqrt(dist2)

            vdx = vdx + 4

        end do    
        
    end subroutine update_megno

    ! Update bodies chaos and MENGO at system
    pure subroutine update_chaos(self, reference_frame)
        implicit none
        type(system_st), intent(inout) :: self
        integer(kind=4), intent(in) :: reference_frame
        real(wp) :: aux_real2(2)
        integer(kind=4) :: i

        call update_elements(self, reference_frame)

        ! Asteroid
        !! a
        aux_real2(1) = min(self%asteroid%chaos_a(1), self%asteroid%elements(1))
        aux_real2(2) = max(self%asteroid%chaos_a(2), self%asteroid%elements(1))
        self%asteroid%chaos_a = aux_real2
        !! e
        aux_real2(1) = min(self%asteroid%chaos_e(1), self%asteroid%elements(2))
        aux_real2(2) = max(self%asteroid%chaos_e(2), self%asteroid%elements(2))
        self%asteroid%chaos_e = aux_real2

        ! Moons
        do i = 1, self%Nmoons_active
            !! a
            aux_real2(1) = min(self%moons(i)%chaos_a(1), self%moons(i)%elements(1))
            aux_real2(2) = max(self%moons(i)%chaos_a(2), self%moons(i)%elements(1))
            self%moons(i)%chaos_a = aux_real2
            !! e
            aux_real2(1) = min(self%moons(i)%chaos_e(1), self%moons(i)%elements(2))
            aux_real2(2) = max(self%moons(i)%chaos_e(2), self%moons(i)%elements(2))
            self%moons(i)%chaos_e = aux_real2
            !! time
            self%moons(i)%tmax = self%time
        end do

        ! Particles
        do i = 1, self%Nparticles_active
            !! a
            aux_real2(1) = min(self%particles(i)%chaos_a(1), self%particles(i)%elements(1))
            aux_real2(2) = max(self%particles(i)%chaos_a(2), self%particles(i)%elements(1))
            self%particles(i)%chaos_a = aux_real2
            !! e
            aux_real2(1) = min(self%particles(i)%chaos_e(1), self%particles(i)%elements(2))
            aux_real2(2) = max(self%particles(i)%chaos_e(2), self%particles(i)%elements(2))
            self%particles(i)%chaos_e = aux_real2
            !! time
            self%particles(i)%tmax = self%time
        end do

    end subroutine update_chaos

    ! Update bodies chaos at system
    pure subroutine update_chaos_geometric(self)
        implicit none
        type(system_st), intent(inout) :: self
        real(wp) :: aux_real2(2)
        integer(kind=4) :: i

        call update_geometric(self)

        ! Moons
        do i = 1, self%Nmoons_active
            !! a
            aux_real2(1) = min(self%moons(i)%chaos_a_geom(1), self%moons(i)%geometric(1))
            aux_real2(2) = max(self%moons(i)%chaos_a_geom(2), self%moons(i)%geometric(1))
            self%moons(i)%chaos_a_geom = aux_real2
            !! e
            aux_real2(1) = min(self%moons(i)%chaos_e_geom(1), self%moons(i)%geometric(2))
            aux_real2(2) = max(self%moons(i)%chaos_e_geom(2), self%moons(i)%geometric(2))
            self%moons(i)%chaos_e_geom = aux_real2
            !! time
            self%moons(i)%tmax = self%time
        end do

        ! Particles
        do i = 1, self%Nparticles_active
            !! a
            aux_real2(1) = min(self%particles(i)%chaos_a_geom(1), self%particles(i)%geometric(1))
            aux_real2(2) = max(self%particles(i)%chaos_a_geom(2), self%particles(i)%geometric(1))
            self%particles(i)%chaos_a_geom = aux_real2
            !! e
            aux_real2(1) = min(self%particles(i)%chaos_e_geom(1), self%particles(i)%geometric(2))
            aux_real2(2) = max(self%particles(i)%chaos_e_geom(2), self%particles(i)%geometric(2))
            self%particles(i)%chaos_e_geom = aux_real2
            !! time
            self%particles(i)%tmax = self%time
        end do

    end subroutine update_chaos_geometric

    ! Update system state. ! Assumes theta, omega, x, y, vx, vy,... ! BARYCENTRIC
    subroutine update_system_from_array(self, time, array, error)
        implicit none
        type(system_st), intent(inout) :: self
        real(wp), intent(in) :: time
        real(wp), dimension(:), intent(in) :: array ! Includes theta and omega
        integer(kind=4), intent(inout), optional :: error
        real(wp), dimension(4) :: aux_coordinates
        real(wp) :: energy, ang_mom, inertia
        integer(kind=4) :: i, idx

        ! Update TIME
        self%time = time

        ! Update Asteroid rotation
        call spin_asteroid(self%asteroid, array(1), array(2))

        ! Update Asteroid Position and derived parameters
        aux_coordinates = array(3:6) !
        call shift_asteroid(self%asteroid, aux_coordinates)

        ! Update Moons Position and derived parameters
        idx = 7  ! First 2 + 4 were the asteroid
        do i = 1, self%Nmoons_active
            aux_coordinates = array(idx:idx + 3)
            call shift_single_moon(self%moons(i), aux_coordinates)
            idx = idx + 4
        end do

        ! Update Particles Position and derived parameters
        do i = 1, self%Nparticles_active
            aux_coordinates = array(idx:idx + 3)
            call shift_single_particle(self%particles(i), aux_coordinates)
            idx = idx + 4
        end do

        ! Check CM
        call check_coordinates(self, error)

        ! Update energy and ang_mom
        call calculate_energy_ang_mom_inertia(self, energy, ang_mom, inertia)
        self%ang_mom = ang_mom  !Set ANG_MOM
        self%energy = energy  !Set ENERGY
        self%inertia = inertia !Set INERTIA

        ! If MEGNO present, we have to add the extra bodies
        if (self%megno_epsilon > tini) then
            do i = 1, self%Nparticles_active
                self%particles(i)%variational = array(idx:idx + 3)
                idx = idx + 4
            end do
        end if

    end subroutine update_system_from_array

    ! Create mass and coordinates. STARTS FROM 1
    pure subroutine generate_arrays(self, mass_arr, radius_arr, array)
        implicit none
        type(system_st), intent(in) :: self
        real(wp), dimension(:), intent(inout) :: mass_arr
        real(wp), dimension(:), intent(inout) :: radius_arr
        real(wp), dimension(:), intent(inout) :: array  ! Includes Theta and Omega
        integer(kind=4) :: i, idx

        ! Set asteroid values
        mass_arr(1) = self%asteroid%mass
        radius_arr(1) = self%asteroid%radius
        array(1) = self%asteroid%theta
        array(2) = self%asteroid%omega
        array(3:6) = self%asteroid%coordinates

        ! Set Moons values
        idx = 7  ! First 4 were the asteroid
        do i = 1, self%Nmoons_active
            mass_arr(i + 1) = self%moons(i)%mass
            radius_arr(i + 1) = self%moons(i)%radius
            array(idx:idx + 3) = self%moons(i)%coordinates
            idx = idx + 4
        end do

        ! Set particles values
        do i = 1, self%Nparticles_active
            mass_arr(self%Nmoons_active + i + 1) = cero  ! Needed ?
            radius_arr(self%Nmoons_active + i + 1) = cero  ! Needed ?
            array(idx:idx + 3) = self%particles(i)%coordinates
            idx = idx + 4
        end do

        ! If MEGNO present, we have to add the extra equations
        if (self%megno_epsilon > tini) then
            do i = 1, self%Nparticles_active
                array(idx:idx + 3) = self%particles(i)%variational
                idx = idx + 4
            end do
        end if

    end subroutine generate_arrays

    ! Update system state. ! Assumes theta, omega, a, e, M, w,...
    subroutine update_system_from_elements(self, time, array, reference_frame, error)
        implicit none
        type(system_st), intent(inout) :: self
        real(wp), intent(in) :: time
        real(wp), dimension(:), intent(in) :: array ! Includes theta and omega
        integer(kind=4), intent(in) :: reference_frame
        integer(kind=4), intent(inout), optional :: error
        real(wp), dimension(size(array)) :: arr_coor
        real(wp), dimension(6) :: aux_coordinates
        real(wp), dimension(4) :: coords_in
        real(wp) :: mass_in, mass_add
        real(wp) :: a_effect
        real(wp) :: mass_cm
        real(wp), dimension(4) :: coords_cm
        integer(kind=4) :: nactive, i, idx

        arr_coor(1) = modulo(array(1), twopi)
        arr_coor(2) = array(2)

        ! Check if a asteroid too small
        if ((array(3) < myepsilon) .or. (reference_frame > 0)) then ! All 0
            arr_coor(3:6) = cero
        else
            call coord(self%mass, &
                    & array(3), &
                    & array(4), &
                    & cero, &
                    & array(5), &
                    & array(6), &
                    & cero, &
                    & aux_coordinates)
            arr_coor(3:6) = aux_coordinates((/1,2,4,5/))
        end if

        ! Get Nactive
        call get_Nactive(self, nactive)

        ! Get inner
        mass_in = self%asteroid%mass
        coords_in = cero

        if ((reference_frame .eq. 0) .or. (reference_frame .eq. 1)) then  ! Barycentric or Jacobi         

            ! Moons (massive)
            do i = 2, self%Nmoons_active + 1
                idx = i*4 - 1
                mass_add = self%moons(i-1)%mass
                if (reference_frame .eq. 0) then ! Barycenter
                    a_effect = array(idx) * (mass_in + mass_add) / mass_in
                else
                    a_effect = array(idx)
                end if
                call coord(mass_in + mass_add, &
                        & a_effect, &
                        & array(idx + 1), &
                        & cero, &
                        & array(idx + 2), &
                        & array(idx + 3), &
                        & cero, &
                        & aux_coordinates)
                ! Jacobi -> barycentric
                arr_coor(idx:idx+3) = aux_coordinates((/1,2,4,5/)) + coords_in
                ! Update barycenter
                coords_in = coords_in*mass_in + arr_coor(idx:idx+3)*mass_add
                mass_in = mass_in + mass_add
                coords_in = coords_in/mass_in
            end do
            
            ! Shift to CM
            do i = 1, self%Nmoons_active + 1
                idx = i*4 - 1
                arr_coor(idx:idx+3) = arr_coor(idx:idx+3) - coords_in
            end do

            ! Particles
            do i = self%Nmoons_active + 2, nactive
                idx = i*4 - 1
                call coord(mass_in, &
                        & array(idx), &
                        & array(idx + 1), &
                        & cero, &
                        & array(idx + 2), &
                        & array(idx + 3), &
                        & cero, &
                        & aux_coordinates)
                arr_coor(idx:idx+3) = aux_coordinates((/1,2,4,5/))
            end do
        
        else ! Astrocentric

            mass_cm = mass_in
            coords_cm = cero

            ! Moons (massive)
            do i = 2, self%Nmoons_active + 1
                idx = i*4 - 1
                call coord(mass_in, &
                        & array(idx), &
                        & array(idx + 1), &
                        & cero, &
                        & array(idx + 2), &
                        & array(idx + 3), &
                        & cero, &
                        & aux_coordinates)
                arr_coor(idx:idx+3) = aux_coordinates((/1,2,4,5/))
                coords_cm = coords_cm + arr_coor(idx:idx+3)*self%moons(i-1)%mass
                mass_cm = mass_cm + self%moons(i-1)%mass
            end do

            ! Particles
            do i = self%Nmoons_active + 2, nactive
                idx = i*4 - 1
                call coord(mass_in, &
                        & array(idx), &
                        & array(idx + 1), &
                        & cero, &
                        & array(idx + 2), &
                        & array(idx + 3), &
                        & cero, &
                        & aux_coordinates)
                arr_coor(idx:idx+3) = aux_coordinates((/1,2,4,5/))
            end do

            ! Now we find the CM and move everything
            coords_cm = coords_cm/mass_cm
            do i = 1, nactive
                idx = i*4 - 1
                arr_coor(idx:idx+3) = arr_coor(idx:idx+3) - coords_cm
            end do

        end if

        call update_system_from_array(self, time, arr_coor, error) ! THey must be barycentric here

    end subroutine update_system_from_elements

    !  -----------------   GRAVITY CALCULATION  ---------------------------

    ! Get acceleration and potential energy from primary
    subroutine get_acc_and_pot_asteroid(self, xy_target, acc, pot, inside)
        implicit none
        type(asteroid_st) :: self
        real(wp), intent(in) :: xy_target(2)
        real(wp), intent(inout) :: acc(2), pot
        logical, intent(inout), optional :: inside
        logical :: has_inside = .False.
        integer(kind=4) :: i
        real(wp) :: dx, dy, dr, dr2, dr4, dz2
        real(wp) :: xy_centered(2)
        real(wp) :: xy_rotated(2)
        real(wp) :: mu, R2, cos2th, sin2th  ! ellipsoid
        real(wp) :: J2K_coef

        ! Init
        has_inside = present(inside)

        ! Primary
        xy_centered = self%primary%coordinates_CM(1:2) + self%coordinates(1:2)

        !! Sphere
        if (self%primary%is_sphere) then
            if (abs(self%primary%C20) > 0) then
                J2K_coef = 1.5e0_wp*self%primary%radius**2*self%primary%C20 
                call get_acc_and_pot_single_with_J2 ( &
                    & self%primary%mass, &
                    & xy_centered, &
                    & xy_target, &
                    & self%primary%radius, &
                    & J2K_coef, &
                    & acc, pot, inside)
            else
                call get_acc_and_pot_single( &
                    & self%primary%mass, &
                    & xy_centered, &
                    & xy_target, &
                    & self%primary%radius, &
                    & acc, pot, inside)
            end if

        !! Ellipsoid
        else
            ! Anti-rotate target to check if inside
            xy_rotated = rotate2D(xy_target, -self%theta)
            dx = xy_rotated(1) - xy_centered(1)
            dy = xy_rotated(2) - xy_centered(2)

            if (has_inside) then
                inside = (dx/self%primary%semi_axis(1))**2 &
                     & + (dy/self%primary%semi_axis(2))**2 < uno
            end if

            dr2 = dx*dx + dy*dy
            dr = sqrt(dr2)

            dr4 = dr2*dr2
            mu = G*self%primary%mass
            R2 = self%primary%radius**2
            cos2th = cos(dos*self%theta)
            sin2th = sin(dos*self%theta)

            if (dr > myepsilon) then
                !! Now, calculate from non rotated
                dx = xy_target(1) - xy_centered(1)
                dy = xy_target(2) - xy_centered(2)

                pot = pot + mu*( &
                        & -dos*dr4 &
                        & + R2*( &
                        & dr2*self%primary%C20 &
                        & - 6*( &
                        & (dx - dy)*(dx + dy)*cos2th &
                        & + dos*dx*dy*sin2th) &
                        & *self%primary%C22) &
                        & )/(dos*dr**5)

                acc(1) = acc(1) + mu*( &
                            & -dos*dr4*dx &
                            & + 3*R2*( &
                            & dr2*dx*self%primary%C20 &
                            & + dos*( &
                            & dx*(dos*dr2 - 5*dx**2 + 5*dy**2)*cos2th &
                            & + dos*(dr2 - 5*dx**2)*dy*sin2th &
                            & )*self%primary%C22) &
                        & )/(dos*dr**7)

                acc(2) = acc(2) + mu*( &
                            & -dos*dr4*dy &
                            & + 3*R2*( &
                            & dr2*dy*self%primary%C20 &
                            & + dos*( &
                            & dy*(-dos*dr2 - 5*dx**2 + 5*dy**2)*cos2th &
                            & + dos*(dr2 - 5*dy**2)*dx*sin2th &
                            & )*self%primary%C22) &
                            & )/(dos*dr**7)
            end if

        end if

        ! Boulders
        do i = 1, self%Nboulders
            !! Skip if already inside
            if (has_inside) then
                if (inside) return
            end if

            ! Center and calculate
            xy_centered = self%boulders(i)%coordinates_CM(1:2) + self%coordinates(1:2)
            call get_acc_and_pot_single( &
                & self%boulders(i)%mass, &
                & xy_centered, &
                & xy_target, &
                & self%boulders(i)%radius, &
                & acc, pot, inside)
        end do

        ! Asteroid as a whole, if C20 from it
        if (abs(self%C20) > 0) then
            ! Center
            dx = xy_target(1) - self%coordinates(1)
            dy = xy_target(2) - self%coordinates(2)
            dr2 = dx*dx + dy*dy
            dr = sqrt(dr2)
            ! Add J2 contribution from CM
            J2K_coef = 1.5e0_wp*self%radius**2*self%C20 
            acc = acc + G*self%mass/(dr2*dr)*(/dx, dy/)*J2K_coef/dr2
            pot = pot + G*self%mass/dr*J2K_coef*uno3/dr2
        
        else if (self%boulder_z%active) then
            ! Center
            dx = xy_target(1) - self%coordinates(1)
            dy = xy_target(2) - self%coordinates(2)
            dz2 = self%boulder_z%dist_to_asteroid**2
            dr2 = dx*dx + dy*dy + dz2
            dr = sqrt(dr2)
            acc = acc - G*self%boulder_z%mass*(/dx, dy/)/(dr2*dr)
            pot = pot - G*self%boulder_z%mass/dr

        end if

    end subroutine get_acc_and_pot_asteroid

    ! Get acceleration and potential energy from single moon
    subroutine get_acc_and_pot_moon(self, xy_target, acc, pot, inside)
        implicit none
        type(moon_st) :: self
        real(wp), intent(in) :: xy_target(2)
        real(wp), intent(inout) :: acc(2), pot
        logical, intent(inout), optional :: inside

        call get_acc_and_pot_single(self%mass, self%coordinates(1:2), xy_target, self%radius, acc, pot, inside)
    end subroutine get_acc_and_pot_moon

    ! Calculate acceleration and potential at a given coordinate
    subroutine get_acc_and_pot_xy(self, xy_target, acc, pot, inside)
        implicit none
        type(system_st), intent(in) :: self
        real(wp), intent(in) :: xy_target(2)
        real(wp), intent(out) :: acc(2), pot
        logical, intent(inout), optional :: inside
        integer(kind=4) :: i
        logical :: has_inside, is_inside

        ! Init
        has_inside = present(inside)
        is_inside = .False.
        acc = cero
        pot = cero

        ! Asteroid
        call get_acc_and_pot_asteroid(self%asteroid, xy_target, acc, pot, is_inside)

        !! Check if inside
        if (is_inside) then
            pot = -G*self%asteroid%mass/self%asteroid%radius  ! Convenience
            acc = cero
            if (has_inside) inside = .True.
            return
        end if

        ! Moons
        do i = 1, self%Nmoons_active
            call get_acc_and_pot_moon(self%moons(i), &
                                      & xy_target, &
                                      & acc, pot, &
                                      is_inside)
            ! Check if inside boulder
            if (is_inside) then
                pot = -G*self%moons(i)%mass/self%moons(i)%radius
                acc = cero
                if (has_inside) inside = .True.
                return
            end if
        end do

    end subroutine get_acc_and_pot_xy

    !  ------------------   COLLISIONS/ESCAPES   --------------------------

    ! Merge a moon into asteroid. NO DEACTIVATION HERE
    subroutine merge_moon_into_ast(asteroid, moon, error)
        implicit none
        type(asteroid_st), intent(inout) :: asteroid
        type(moon_st), intent(inout) :: moon
        integer(kind=4), intent(inout), optional :: error
        real(wp) :: m_cm, rv_cm(4)
        real(wp) :: new_ang_mom_rot
        real(wp) :: aux_real24(2, 4)

        ! Check if moon is active
        if (.not. moon%active) then
            write (*, *) "WARNING: Merging an inactive moon into the asteroid:", moon%id
            if (present(error)) error = error + 1
        end if

        ! Get these 2 bodies CM properties
        m_cm = asteroid%mass + moon%mass
        rv_cm = (asteroid%mass*asteroid%coordinates + moon%mass*moon%coordinates)/m_cm

        ! Calculate L_rot from this CM
        aux_real24(1, :) = asteroid%coordinates - rv_cm
        aux_real24(2, :) = moon%coordinates - rv_cm
        new_ang_mom_rot = (aux_real24(1, 1)*aux_real24(1, 4) - aux_real24(1, 2)*aux_real24(1, 3))*asteroid%mass + &
                        & (aux_real24(2, 1)*aux_real24(2, 4) - aux_real24(2, 2)*aux_real24(2, 3))*moon%mass

        ! Add rot ang mom that the asteroid has before merging (and the moon too!)
        new_ang_mom_rot = new_ang_mom_rot + moon%ang_mom_rot + asteroid%ang_mom_rot

        ! Update MASS and derivates
        call grow_asteroid(asteroid, moon%mass)

        ! Update ROTATION
        call spin_asteroid(asteroid, asteroid%theta, new_ang_mom_rot/asteroid%inertia)

        ! Update COORDINATES and derivates
        !! Get total mass if given
        call shift_asteroid(asteroid, rv_cm)
    end subroutine merge_moon_into_ast

    ! Merge a moon(i) into asteroid, in a system. DEACTIVATION HERE.
    subroutine merge_moon_i_into_ast(self, i, error)
        implicit none
        type(system_st), intent(inout) :: self
        integer(kind=4), intent(in) :: i
        integer(kind=4), intent(inout), optional :: error
        real(wp) :: aux_real, ang_mom

        ! Perform the merge
        call merge_moon_into_ast(self%asteroid, self%moons(i), error)

        ! Add the flag
        self%moons(i)%merged_to = 0 ! Merged to asteroid

        ! Deactivate the moon
        call deactivate_moon_i(self, i, .False., error)  ! Possible reordering done outside

        ! Check CM
        call check_coordinates(self, error)

        ! Check ANGULAR MOMENTUM
        call calculate_energy_ang_mom_inertia(self, aux_real, ang_mom, aux_real)
        aux_real = max(sqepsilon, abs(ang_mom))
        if (abs(ang_mom - self%ang_mom)/aux_real > sqepsilon) then
            write (*, *) "WARNING: Total angular momentum relative error:", abs(ang_mom - self%ang_mom)/aux_real
            if (present(error)) error = error + 1
            return
        end if

        ! ENERGY loss ocurrs, and the amount is given by:
        !! dE = G m1 m2 / r - (m1 m2 / (m1 + m2)) (v2 - v1)**2 / 2
    end subroutine merge_moon_i_into_ast

    ! Eject a moon(i), in a system. DEACTIVATION HERE
    subroutine eject_moon_i(self, i, error)
        implicit none
        type(system_st), intent(inout) :: self
        integer(kind=4), intent(in) :: i
        integer(kind=4), intent(inout), optional :: error
        logical :: keep_order

        ! Add the flag
        self%moons(i)%merged_to = -2 ! Ejected

        ! Deactivate the moon
        keep_order = self%reference_frame .eq. 1 ! Only if Jacobi
        call deactivate_moon_i(self, i, keep_order, error)

        ! Update almost all
        call recalculate_all(self)
    end subroutine eject_moon_i

    ! Resolve escapes
    subroutine resolve_escapes(self, r_max, unit_file)
        implicit none
        type(system_st), intent(inout) :: self
        real(wp), intent(in) :: r_max
        integer(kind=4), intent(in), optional :: unit_file
        integer(kind=4) :: i
        integer(kind=4) :: m_act0
        logical :: again, do_write
        real(wp) :: r_max2, dist2_to_ast, ast_coord(2), aux_coord(2)

        if (r_max .le. cero) return  ! Nothing to do

        ! Init
        do_write = present(unit_file) .and. unit_file > 0   !!! -std=f08
        ast_coord = self%asteroid%coordinates(1:2)
        r_max2 = r_max*r_max

        again = .True.
        do while (again)
            m_act0 = self%Nmoons_active

            ! Moons
            do i = m_act0, 1, -1  ! Backwards loop
                aux_coord = self%moons(i)%coordinates(1:2) - ast_coord
                dist2_to_ast = aux_coord(1)*aux_coord(1) + aux_coord(2)*aux_coord(2)
                ! if (self%moons(i)%dist_to_cm > r_max) then
                if (dist2_to_ast > r_max2) then
                    if (do_write) write (unit_file, s1i5x5) "Moon", i, "(", self%moons(i)%id, ") escaped."
                    call eject_moon_i(self, i)
                end if
            end do

            ! Particles
            do i = self%Nparticles_active, 1, -1  ! Backwards loop
                aux_coord = self%particles(i)%coordinates(1:2) - ast_coord
                dist2_to_ast = aux_coord(1)*aux_coord(1) + aux_coord(2)*aux_coord(2)
                ! if (self%particles(i)%dist_to_cm > r_max) then
                if (dist2_to_ast > r_max2) then
                    self%particles(i)%merged_to = -2  ! Escape
                    if (do_write) write (unit_file, s1i5x5) "Particle", i, "(", self%particles(i)%id, ") escaped."
                    call deactivate_particle_i(self, i)
                end if
            end do
            again = (self%Nmoons_active .ge. 1) .and. (m_act0 > self%Nmoons_active)

        end do
    end subroutine resolve_escapes

    ! Collide a moon(j) with a moon(i), in a system
    subroutine collide_2_moons(self, i, j, outcome, error)
        implicit none
        type(system_st), intent(inout) :: self
        integer(kind=4), intent(in) :: i, j
        integer(kind=4), intent(out) :: outcome
        integer(kind=4), intent(inout), optional :: error
        real(wp) :: m_cm, rv_cm(4)
        real(wp) :: L_rot, L_orb
        real(wp) :: aux_real, ang_mom, aux_real24(2, 4)
        real(wp) :: mi, mj
        real(wp) :: ri(2), rj(2)
        real(wp) :: vi(2), vj(2)
        real(wp) :: dr_vec(2), dr, dr_ver(2)
        real(wp) :: dv_vec(2), dv2
        real(wp) :: dv_rad
        real(wp) :: Li_orb, Lj_orb
        real(wp) :: mu, Ekin, Epot, vesc2
        real(wp) :: vin, vit, vjn, vjt
        real(wp) :: vin_new, vjn_new
        real(wp) :: Li_tmp, Lj_tmp
        real(wp) :: L_tmp
        real(wp) :: denom
        real(wp) :: delta_vi, delta_vj
        real(wp) :: vi_col(2), vj_col(2)
        real(wp) :: overlap, corr1, corr2
        real(wp) :: new_rvi(4), new_rvj(4)

        ! Result
        outcome = 0

        ! Individual attrs
        ri = self%moons(i)%coordinates(1:2)
        rj = self%moons(j)%coordinates(1:2)

        ! Relative attrs
        dr_vec = ri - rj  ! Relative pos
        dr = sqrt(dr_vec(1)*dr_vec(1) + dr_vec(2)*dr_vec(2))  ! Distance

        ! Skip if not touching
        if (dr .gt. (self%moons(j)%radius + self%moons(i)%radius)) return ! Nothing to do

        ! Relative position versor
        dr_ver = dr_vec/dr

        ! Masses
        mi = self%moons(i)%mass
        mj = self%moons(j)%mass

        ! Relative attrs
        m_cm = mi + mj  ! Combined mass

        ! Velocities
        vi = self%moons(i)%coordinates(3:4)
        vj = self%moons(j)%coordinates(3:4)

        ! Relative attrs
        dv_vec = vi - vj  ! Relative vel
        dv_rad = dv_vec(1)*dr_ver(1) + dv_vec(2)*dr_ver(2)

        ! Are they already moving appart ?
        if (dv_rad > cero) then
            ! Push apart proportionally to masses
            overlap = self%moons(i)%radius + self%moons(j)%radius - dr
            if (overlap > cero) then
                corr1 = overlap*(mj/m_cm)
                corr2 = overlap*(mi/m_cm)

                new_rvi(1:2) = ri + corr1*dr_ver
                new_rvj(1:2) = rj - corr2*dr_ver

                new_rvi(3:4) = vi
                new_rvj(3:4) = vj

                ! Keep velocities unchanged
                call shift_single_moon(self%moons(i), new_rvi)
                call shift_single_moon(self%moons(j), new_rvj)

            end if

            return  ! Nothing else to do

        end if

        ! Energy
        dv2 = dv_vec(1)*dv_vec(1) + dv_vec(2)*dv_vec(2)  ! Squared relative velocity
        mu = mi*mj/m_cm  ! Reduced mass
        Ekin = uno2*mu*dv2
        Epot = -G*mi*mj/dr

        ! Escape velocity
        vesc2 = dos*G*m_cm/dr

        ! Check if moons are active
        if ((.not. self%moons(i)%active) .or. ((.not. self%moons(j)%active))) then
            if (.not. self%moons(j)%active) write (*, *) "WARNING: Colliding an inactive moon:", self%moons(j)%active
            if (.not. self%moons(i)%active) write (*, *) "WARNING: Colliding to an inactive moon:", self%moons(i)%active
            if (present(error)) error = error + 1
        end if

        ! Check if merge
        !! Merge if bounded or eta = 1
        if (((Ekin + Epot < -self%f_col*abs(Epot)) .and. (dv2 < vesc2*(uno - self%f_col))) .or. &
            & (self%eta_col .ge. uno)) then

            ! Result (merge)
            outcome = 1

            ! Get these 2 moons CM properties
            rv_cm = (mi*self%moons(i)%coordinates + mj*self%moons(j)%coordinates)/m_cm

            ! Calculate L_rot from this CM
            aux_real24(1, :) = self%moons(i)%coordinates - rv_cm
            aux_real24(2, :) = self%moons(j)%coordinates - rv_cm
            L_rot = (aux_real24(1, 1)*aux_real24(1, 4) - aux_real24(1, 2)*aux_real24(1, 3))*mi + &
                & (aux_real24(2, 1)*aux_real24(2, 4) - aux_real24(2, 2)*aux_real24(2, 3))*mj

            ! Update MASS and derivates
            self%moons(i)%mass = m_cm
            self%moons(i)%mu_to_asteroid = m_cm/self%asteroid%mass

            ! Update RADIUS
            self%moons(i)%radius = max(self%moons(i)%radius, self%moons(j)%radius)   ! HEURISTIC Maybe mean density ???

            ! Update ENERGY and ANGULAR MOMENTUM
            self%moons(i)%inertia = 0.4e0_wp*m_cm*self%moons(i)%radius*self%moons(i)%radius
            self%moons(i)%ang_mom_rot = self%moons(i)%ang_mom_rot + self%moons(j)%ang_mom_rot + L_rot
            self%moons(i)%e_rot = self%moons(i)%e_rot + self%moons(j)%e_rot + (uno2*L_rot*L_rot/self%moons(i)%inertia)

            ! Add the flag
            self%moons(j)%merged_to = self%moons(i)%id  ! Merged j into i

            ! Deactivate the moon
            call deactivate_moon_i(self, j, .False., error)  ! Possible reordering done below

            ! Update COORDINATES and derivates
            call shift_single_moon(self%moons(i), rv_cm)

        else  ! Not plastic

            ! Result (not plastic)
            outcome = 2

            ! Angular momentum BEFORE collision (orbital part)
            Li_orb = mi*(ri(1)*vi(2) - ri(2)*vi(1))
            Lj_orb = mj*(rj(1)*vj(2) - rj(2)*vj(1))
            L_orb = Li_orb + Lj_orb

            ! Project velocities
            vin = vi(1)*dr_ver(1) + vi(2)*dr_ver(2)
            vit = -vi(1)*dr_ver(2) + vi(2)*dr_ver(1)
            vjn = vj(1)*dr_ver(1) + vj(2)*dr_ver(2)
            vjt = -vj(1)*dr_ver(2) + vj(2)*dr_ver(1)

            ! Normal after collision (restitution handled by self%eta_col)
            vin_new = ((mi*vin + mj*vjn) - (uno - self%eta_col)*mj*(vin - vjn))/m_cm
            vjn_new = ((mi*vin + mj*vjn) + (uno - self%eta_col)*mi*(vin - vjn))/m_cm

            ! --- Position correction (push apart to avoid overlap) ---
            overlap = self%moons(i)%radius + self%moons(j)%radius - dr
            if (overlap > cero) then
                corr1 = overlap*(mj/m_cm)
                corr2 = overlap*(mi/m_cm)
                new_rvi(1) = ri(1) + corr1*dr_ver(1)
                new_rvi(2) = ri(2) + corr1*dr_ver(2)
                new_rvj(1) = rj(1) - corr2*dr_ver(1)
                new_rvj(2) = rj(2) - corr2*dr_ver(2)
            else
                new_rvi(1:2) = ri
                new_rvj(1:2) = rj
            end if

            ! Build the post-impulse cartesian velocities (before tangential tweak)
            ! vi_col and vj_col are the velocities after impulse using vin_new/vit etc.
            vi_col(1) = vin_new*dr_ver(1) - vit*dr_ver(2)
            vi_col(2) = vin_new*dr_ver(2) + vit*dr_ver(1)
            vj_col(1) = vjn_new*dr_ver(1) - vjt*dr_ver(2)
            vj_col(2) = vjn_new*dr_ver(2) + vjt*dr_ver(1)

            ! Angular momentum AFTER position correction but BEFORE tangential tweak
            Li_tmp = mi*(new_rvi(1)*vi_col(2) - new_rvi(2)*vi_col(1))
            Lj_tmp = mj*(new_rvj(1)*vj_col(2) - new_rvj(2)*vj_col(1))
            L_tmp = Li_tmp + Lj_tmp

            ! Denominator: lever arms for tangential correction (explicit symmetric form)
            denom = mi*((new_rvi(1) - new_rvj(1))*dr_ver(1) + (new_rvi(2) - new_rvj(2))*dr_ver(2))

            if (abs(denom) > myepsilon) then
                ! Solve for tangential velocity correction (delta on vi_col along t)
                delta_vi = (L_orb - L_tmp)/denom
                delta_vj = -(mi/mj)*delta_vi

                ! Apply corrections along tangent direction, preserving linear momentum
                vi_col(1) = vi_col(1) - delta_vi*dr_ver(2)
                vi_col(2) = vi_col(2) + delta_vi*dr_ver(1)
                vj_col(1) = vj_col(1) - delta_vj*dr_ver(2)
                vj_col(2) = vj_col(2) + delta_vj*dr_ver(1)
            end if

            ! Pack the new states to pass to shift_single_moon
            new_rvi(3) = vi_col(1)
            new_rvi(4) = vi_col(2)
            new_rvj(3) = vj_col(1)
            new_rvj(4) = vj_col(2)

            ! Update COORDINATES and VELOCITIES
            call shift_single_moon(self%moons(i), new_rvi)
            call shift_single_moon(self%moons(j), new_rvj)

        end if

        ! Check CM
        call check_coordinates(self, error)

        ! Check ANGULAR MOMENTUM
        call calculate_energy_ang_mom_inertia(self, aux_real, ang_mom, aux_real)
        aux_real = max(sqepsilon, abs(ang_mom))
        if (abs(ang_mom - self%ang_mom)/aux_real > sqepsilon) then
            write (*, *) "WARNING: Total angular momentum differs from previous. Err_rel:", abs(ang_mom - self%ang_mom)/aux_real
            if (present(error)) error = error + 1
        end if

        ! ENERGY loss may ocurr, and the amount is given by something like:
        !! dE = G m1 m2 / r - (m1 m2 / (m1 + m2)) (v2 - v1)**2 / 2
    end subroutine collide_2_moons

    ! Resolve collisions
    subroutine resolve_collisions(self, r_min, unit_file)
        implicit none
        type(system_st), intent(inout) :: self
        real(wp), intent(in) :: r_min
        integer(kind=4), intent(in), optional :: unit_file
        integer(kind=4) :: i, j
        integer(kind=4) :: m_act0
        integer(kind=4) :: m_act, p_act
        integer(kind=4) :: moon_id
        integer(kind=4) :: outcome
        real(wp) :: boul_pos(2), moon_pos(2)
        real(wp) :: boul_rad, moon_rad
        real(wp) :: dr_vec(2), dr
        real(wp) :: xy_rotated(2), dx, dy  ! For ellipsoid
        logical :: again, do_write, was_any_moon_something

        ! Init
        do_write = present(unit_file) .and. unit_file > 0   !!! -std=f08
        
        was_any_moon_something = .False.
        again = .True.
        do while (again)

            again = .False.
            m_act0 = self%Nmoons_active

            ! Moons -> Moons
            do j = m_act0, 2, -1  ! Backwards loop
                inner_loop: do i = j - 1, 1, -1  ! Backwards loop
                    call collide_2_moons(self, i, j, outcome)
                    if (outcome > 0) then
                        if (do_write) then
                            if (outcome == 1) then
                                write (unit_file, s1i5x5) "Merged moon ", j, &
                                                & "(", self%moons(j)%id, ") into moon ", &
                                                & i, "(", self%moons(i)%id, ")."
                            else
                                write (unit_file, s1i5x5) "Collision between moon ", j, &
                                                & "(", self%moons(j)%id, ") and moon ", &
                                                & i, "(", self%moons(i)%id, ")."
                            end if
                        end if
                        again = .True.
                        was_any_moon_something = .True.
                        exit inner_loop
                    end if
                end do inner_loop
            end do

            ! Moons -> Primary + Boulders
            m_act = self%Nmoons_active
            if (r_min .le. cero) then  ! Only if not r_min

                ! Primary
                boul_pos = self%asteroid%primary%coordinates_CM(1:2) + self%asteroid%coordinates(1:2)
                if (self%asteroid%primary%is_sphere) then  ! Sphere
                    boul_rad = self%asteroid%primary%radius
                    do i = m_act, 1, -1  ! Backwards loop
                        dr_vec = self%moons(i)%coordinates(1:2) - boul_pos  ! Relative pos
                        dr = sqrt(dr_vec(1)*dr_vec(1) + dr_vec(2)*dr_vec(2))  ! Distance
                        if (dr .le. boul_rad + self%moons(i)%radius) then
                            if (do_write) write (unit_file, s1i5x5) "Merged moon ", i, "(", self%moons(i)%id, ") into asteroid."
                            call merge_moon_i_into_ast(self, i)
                            was_any_moon_something = .True.
                        end if
                    end do
                else  ! Ellipsoid
                    do i = m_act, 1, -1  ! Backwards loop
                        xy_rotated = rotate2D(self%moons(i)%coordinates(1:2), -self%asteroid%theta)
                        dx = xy_rotated(1) - boul_pos(1) + self%moons(i)%radius
                        dy = xy_rotated(2) - boul_pos(2) + self%moons(i)%radius
                        if (((dx/self%asteroid%primary%semi_axis(1))**2 &
                         & + (dy/self%asteroid%primary%semi_axis(2))**2) < uno) then
                            if (do_write) write (unit_file, s1i5x5) "Merged moon ", i, "(", self%moons(i)%id, ") into asteroid."
                            call merge_moon_i_into_ast(self, i)
                            was_any_moon_something = .True.
                        end if
                    end do
                end if

                m_act = self%Nmoons_active
                ! Boulders
                do j = 1, self%asteroid%Nboulders
                    boul_pos = self%asteroid%boulders(j)%coordinates_CM(1:2) + self%asteroid%coordinates(1:2)
                    boul_rad = self%asteroid%boulders(j)%radius
                    do i = m_act, 1, -1  ! Backwards loop
                        dr_vec = self%moons(i)%coordinates(1:2) - boul_pos  ! Relative pos
                        dr = sqrt(dr_vec(1)*dr_vec(1) + dr_vec(2)*dr_vec(2))  ! Distance
                        if (dr .le. boul_rad + self%moons(i)%radius) then
                            if (do_write) write (unit_file, s1i5x5) "Merged moon ", i, "(", self%moons(i)%id, ") into asteroid."
                            call merge_moon_i_into_ast(self, i)
                        end if
                    end do
                end do

            else ! Check only if distance to asteroid is lower than r_min

                do i = m_act, 1, -1  ! Backwards loop
                    dr_vec = self%moons(i)%coordinates(1:2) - self%asteroid%coordinates(1:2)  ! Relative pos
                    dr = sqrt(dr_vec(1)*dr_vec(1) + dr_vec(2)*dr_vec(2))  ! Distance
                    if (dr .le. r_min + self%moons(i)%radius) then
                        if (do_write) write (unit_file, s1i5x5) "Merged moon ", i, "(", self%moons(i)%id, ") into asteroid."
                        call merge_moon_i_into_ast(self, i)
                        was_any_moon_something = .True.
                    end if
                end do

            end if

            again = again .or. ((self%Nmoons_active .ge. 1) .and. (m_act0 > self%Nmoons_active))

        end do

        ! Particles -> Moons
        m_act = self%Nmoons_active
        p_act = self%Nparticles_active
        do j = 1, m_act
            moon_pos = self%moons(j)%coordinates(1:2)
            moon_rad = self%moons(j)%radius
            moon_id = self%moons(j)%id
            do i = p_act, 1, -1  ! Backwards loop
                dr_vec = self%particles(i)%coordinates(1:2) - moon_pos  ! Relative pos
                dr = sqrt(dr_vec(1)*dr_vec(1) + dr_vec(2)*dr_vec(2))  ! Distance
                if (dr .le. moon_rad) then
                    if (do_write) write (unit_file, s1i5x5) "Merged particle ", i, "(", self%particles(i)%id, ") into moon ", &
                                    & j, "(", self%moons(j)%id, ")."
                    self%particles(i)%merged_to = moon_id  ! Set where merged to
                    call deactivate_particle_i(self, i)  ! Deactivate
                end if
            end do
        end do

        ! Particles -> Primary + Boulders (Only if not r_min)
        p_act = self%Nparticles_active
        if (r_min .le. cero) then

            ! Primary
            boul_pos = self%asteroid%primary%coordinates_CM(1:2) + self%asteroid%coordinates(1:2)
            if (self%asteroid%primary%is_sphere) then  ! Sphere
                boul_rad = self%asteroid%primary%radius
                do i = p_act, 1, -1  ! Backwards loop
                    dr_vec = self%particles(i)%coordinates(1:2) - boul_pos  ! Relative pos
                    dr = sqrt(dr_vec(1)*dr_vec(1) + dr_vec(2)*dr_vec(2))  ! Distance
                    if (dr .le. boul_rad) then
                        if (do_write) write (unit_file, s1i5x5) "Merged particle ", i, &
                                                            & "(", self%particles(i)%id, ") into asteroid."
                        self%particles(i)%merged_to = 0  ! Set where merged to (Asteroid)
                        call deactivate_particle_i(self, i)  ! Deactivate
                    end if
                end do
            else  ! Ellipsoid
                do i = p_act, 1, -1  ! Backwards loop
                    xy_rotated = rotate2D(self%particles(i)%coordinates(1:2), -self%asteroid%theta)
                    dx = xy_rotated(1) - boul_pos(1)
                    dy = xy_rotated(2) - boul_pos(2)
                    if (((dx/self%asteroid%primary%semi_axis(1))**2 &
                     & + (dy/self%asteroid%primary%semi_axis(2))**2) < uno) then
                        if (do_write) write (unit_file, s1i5x5) "Merged particle ", i, &
                                                            & "(", self%particles(i)%id, ") into asteroid."
                        self%particles(i)%merged_to = 0  ! Set where merged to (Asteroid)
                        call deactivate_particle_i(self, i)  ! Deactivate
                    end if
                end do
            end if

            p_act = self%Nparticles_active
            ! Boulders
            do j = 1, self%asteroid%Nboulders
                boul_pos = self%asteroid%boulders(j)%coordinates_CM(1:2) + self%asteroid%coordinates(1:2)
                boul_rad = self%asteroid%boulders(j)%radius
                do i = p_act, 1, -1  ! Backwards loop
                    dr_vec = self%particles(i)%coordinates(1:2) - boul_pos  ! Relative pos
                    dr = sqrt(dr_vec(1)*dr_vec(1) + dr_vec(2)*dr_vec(2))  ! Distance
                    if (dr .le. boul_rad) then
                        if (do_write) write (unit_file, s1i5x5) "Merged particle ", i, &
                                                            & "(", self%particles(i)%id, ") into asteroid."
                        self%particles(i)%merged_to = 0  ! Set where merged to (Asteroid)
                        call deactivate_particle_i(self, i)  ! Deactivate
                    end if
                end do
            end do

        else ! Check only if distance to asteroid is lower than r_min

            do i = p_act, 1, -1  ! Backwards loop
                dr_vec = self%particles(i)%coordinates(1:2) - self%asteroid%coordinates(1:2)  ! Relative pos
                dr = sqrt(dr_vec(1)*dr_vec(1) + dr_vec(2)*dr_vec(2))  ! Distance
                if (dr .le. r_min) then
                    if (do_write) write (unit_file, s1i5x5) "Merged particle ", i, "(", self%particles(i)%id, ") into asteroid."
                    self%particles(i)%merged_to = 0  ! Set where merged to (Asteroid)
                    call deactivate_particle_i(self, i)  ! Deactivate
                end if
            end do

        end if

        ! Reorder if needed
        if ((self%reference_frame .eq. 1 ) .and. was_any_moon_something) call sort_moons(self%moons)
    end subroutine resolve_collisions

    !  -------------------------   IO    ----------------------------------

    ! Write elements asteroid
    subroutine write_ast_elem(self, unit_file)
        implicit none
        type(system_st), intent(in) :: self
        integer(kind=4), intent(in) :: unit_file

        write (unit_file, i2r16) &
            & 0, &  ! ID
            & -1, &  ! type
            & self%time/unit_time, &  ! time
            & self%asteroid%theta/radian, &  ! theta
            & self%asteroid%omega*unit_time, &  ! omega
            & self%asteroid%elements(1)/unit_dist, &  ! a
            & self%asteroid%elements(2), &  ! e
            & self%asteroid%elements(3)/radian, &  ! M
            & self%asteroid%elements(4)/radian, &  ! w
            & cero, &   ! MMR
            & self%asteroid%mass/unit_mass, &  ! mass
            & self%asteroid%radius/unit_dist, &  ! radius
            & self%asteroid%dist_to_cm/unit_dist, &  ! distance
            & self%asteroid%chaos_a/unit_dist, &  ! da
            & self%asteroid%chaos_e, & ! de
            & cero  ! MEGNO
    end subroutine write_ast_elem

    ! Write elements moon i
    subroutine write_moon_i_elem(self, i, unit_file)
        implicit none
        type(system_st), intent(in) :: self
        integer(kind=4), intent(in) :: i, unit_file

        write (unit_file, i2r16) &
            & self%moons(i)%id, &  ! ID
            & 1, &  ! type
            & self%time/unit_time, &  ! time
            & self%asteroid%theta/radian, &  ! theta
            & self%asteroid%omega*unit_time, &  ! omega
            & self%moons(i)%elements(1)/unit_dist, &  ! a
            & self%moons(i)%elements(2), &  ! e
            & self%moons(i)%elements(3)/radian, &  ! M
            & self%moons(i)%elements(4)/radian, &  ! w
            & self%moons(i)%mmr, &  ! MMR
            & self%moons(i)%mass/unit_mass, &  ! mass
            & self%moons(i)%radius/unit_dist, &  ! radius
            & self%moons(i)%dist_to_cm/unit_dist, &  ! distance
            & self%moons(i)%chaos_a/unit_dist, &  ! da
            & self%moons(i)%chaos_e, & ! de
            & cero  ! MEGNO
    end subroutine write_moon_i_elem

    ! Write elements particle i
    subroutine write_particle_i_elem(self, i, unit_file)
        implicit none
        type(system_st), intent(in) :: self
        integer(kind=4), intent(in) :: i, unit_file
        real(kind=wp) :: megno
        if (self%time > tini) then
            megno = dos * self%particles(i)%megno / self%time
        else
            megno = cero
        end if

        write (unit_file, i2r16) &
            & self%particles(i)%id, &   ! ID
            & 2, &  ! type
            & self%time/unit_time, &  ! time
            & self%asteroid%theta/radian, &  ! theta
            & self%asteroid%omega*unit_time, &  ! omega
            & self%particles(i)%elements(1)/unit_dist, &  ! a
            & self%particles(i)%elements(2), &  ! e
            & self%particles(i)%elements(3)/radian, &  ! M
            & self%particles(i)%elements(4)/radian, &  ! w
            & self%particles(i)%mmr, &  ! MMR
            & cero, &  ! mass
            & cero, &  ! radius
            & self%particles(i)%dist_to_cm/unit_dist, &  ! distance
            & self%particles(i)%chaos_a/unit_dist, &  ! da
            & self%particles(i)%chaos_e, & ! de
            & megno  ! MEGNO
    end subroutine write_particle_i_elem

    ! Write elements ALL
    subroutine write_elem(self, unit_file)
        implicit none
        type(system_st), intent(in) :: self
        integer(kind=4), intent(in) :: unit_file
        integer(kind=4) :: i
        integer(kind=4), dimension(:), allocatable :: ids

        allocate (ids(max(self%Nmoons_active, self%Nparticles_active)))
        if (self%asteroid%chaos_a(2) > myepsilon) call write_ast_elem(self, unit_file)

        if (self%Nmoons_active > 1) then
            do i = 1, self%Nmoons_active
                ids(i) = i
            end do
            call quickargsort_int(self%moons%id, ids, 1, self%Nmoons_active)
            do i = 1, self%Nmoons_active
                call write_moon_i_elem(self, ids(i), unit_file)
            end do
        else if (self%Nmoons_active == 1) then
            call write_moon_i_elem(self, 1, unit_file)
        end if

        if (self%Nparticles_active > 1) then
            do i = 1, self%Nparticles_active
                ids(i) = i
            end do
            call quickargsort_int(self%particles%id, ids, 1, self%Nparticles_active)
            do i = 1, self%Nparticles_active
                call write_particle_i_elem(self, ids(i), unit_file)
            end do
        else if (self%Nparticles_active == 1) then
            call write_particle_i_elem(self, 1, unit_file)
        end if

        deallocate (ids)

    end subroutine write_elem

    ! Write coordinates asteroid
    subroutine write_ast_coor(self, unit_file)
        implicit none
        type(system_st), intent(in) :: self
        integer(kind=4), intent(in) :: unit_file
        integer(kind=4) :: i
        real(wp) :: coords(4)

        ! Asteroid
        write (unit_file, i2r9) &
                & - 1, &  ! ID
                & -1, &  ! type
                & self%time/unit_time, &  !time
                & self%asteroid%theta/radian, &  ! theta
                & self%asteroid%omega*unit_time, &  ! omega
                & self%asteroid%coordinates(1:2)/unit_dist, &  ! x y
                & self%asteroid%coordinates(3:4)/unit_vel, &  ! vx vy
                & self%asteroid%mass/unit_mass, &  ! mass
                & self%asteroid%radius/unit_dist  ! radius

        ! Else, only if at least 1 boulder
        if (self%asteroid%Nboulders == 0) return

        ! Primary
        coords = self%asteroid%primary%coordinates_CM + self%asteroid%coordinates
        write (unit_file, i2r9) &
            & 0, &  ! ID
            & 0, &  ! type
            & self%time/unit_time, &  !time
            & (self%asteroid%theta + self%asteroid%primary%initial_theta)/radian, &  ! theta
            & self%asteroid%omega*unit_time, &  ! omega
            & coords(1:2)/unit_dist, &  ! x y
            & coords(3:4)/unit_vel, &  ! vx vy
            & self%asteroid%primary%mass/unit_mass, &  ! mass
            & self%asteroid%primary%radius/unit_dist  ! radius

        ! Boulders
        do i = 1, self%asteroid%Nboulders
            coords = self%asteroid%boulders(i)%coordinates_CM + self%asteroid%coordinates
            write (unit_file, i2r9) &
                & i, &  ! ID
                & 0, &  ! type
                & self%time/unit_time, &  !time
                & (self%asteroid%theta + self%asteroid%boulders(i)%initial_theta)/radian, &  ! theta
                & self%asteroid%omega*unit_time, &  ! omega
                & coords(1:2)/unit_dist, &  ! x y
                & coords(3:4)/unit_vel, &  ! vx vy
                & self%asteroid%boulders(i)%mass/unit_mass, &  ! mass
                & self%asteroid%boulders(i)%radius/unit_dist  ! radius
        end do

    end subroutine write_ast_coor

    ! Write coordinates moon i
    subroutine write_moon_i_coor(self, i, unit_file)
        implicit none
        type(system_st), intent(in) :: self
        integer(kind=4), intent(in) :: i, unit_file

        write (unit_file, i2r9) &
            & self%moons(i)%id + self%asteroid%Nboulders, &  ! ID + Nbould to avoid duplicates
            & 1, &  ! type
            & self%time/unit_time, &  !time
            & self%asteroid%theta/radian, &  ! theta
            & self%asteroid%omega*unit_time, &  ! omega
            & self%moons(i)%coordinates(1:2)/unit_dist, &  ! x y
            & self%moons(i)%coordinates(3:4)/unit_vel, &  ! vx vy
            & self%moons(i)%mass/unit_mass, &  ! mass
            & self%moons(i)%radius/unit_dist  ! radius
    end subroutine write_moon_i_coor

    ! Write coordinates particle i
    subroutine write_particle_i_coor(self, i, unit_file)
        implicit none
        type(system_st), intent(in) :: self
        integer(kind=4), intent(in) :: i, unit_file

        write (unit_file, i2r9) &
            & self%particles(i)%id + self%asteroid%Nboulders, &  ! ID + Nbould to avoid duplicates
            & 2, &  ! type
            & self%time/unit_time, &  !time
            & self%asteroid%theta/radian, &  ! theta
            & self%asteroid%omega*unit_time, &  ! omega
            & self%particles(i)%coordinates(1:2)/unit_dist, &  ! x y
            & self%particles(i)%coordinates(3:4)/unit_vel, &  ! vx vy
            & cero, &  ! mass
            & cero  ! radius
    end subroutine write_particle_i_coor

    ! Write coordinates ALL
    subroutine write_coor(self, unit_file)
        implicit none
        type(system_st), intent(in) :: self
        integer(kind=4), intent(in) :: unit_file
        integer(kind=4) :: i
        integer(kind=4), dimension(:), allocatable :: ids

        allocate (ids(max(self%Nmoons_active, self%Nparticles_active)))
        call write_ast_coor(self, unit_file)

        if (self%Nmoons_active > 1) then
            do i = 1, self%Nmoons_active
                ids(i) = i
            end do
            call quickargsort_int(self%moons%id, ids, 1, self%Nmoons_active)
            do i = 1, self%Nmoons_active
                call write_moon_i_coor(self, ids(i), unit_file)
            end do
        else if (self%Nmoons_active == 1) then
            call write_moon_i_coor(self, 1, unit_file)
        end if

        if (self%Nparticles_active > 1) then
            do i = 1, self%Nparticles_active
                ids(i) = i
            end do
            call quickargsort_int(self%particles%id, ids, 1, self%Nparticles_active)
            do i = 1, self%Nparticles_active
                call write_particle_i_coor(self, ids(i), unit_file)
            end do
        else if (self%Nparticles_active == 1) then
            call write_particle_i_coor(self, 1, unit_file)
        end if

        deallocate (ids)

    end subroutine write_coor

    ! Write both: elements and coordinates. Unit_file and Unit_file + 1
    subroutine write_both(self, unit_file)
        implicit none
        type(system_st), intent(in) :: self
        integer(kind=4), intent(in) :: unit_file
        call write_elem(self, unit_file)
        call write_coor(self, unit_file + 1)
    end subroutine write_both

    ! Write chaos asteroid
    subroutine write_ast_chaos(initial, actual, unit_file)
        implicit none
        type(system_st), intent(in) :: initial
        type(system_st), intent(in) :: actual
        integer(kind=4), intent(in) :: unit_file

        write (unit_file, i3r24) &
            & 0, &  ! ID
            & -1, &  ! type
            & -1, &  ! merged_to
            & actual%time/unit_time, &  ! time
            & initial%asteroid%theta/radian, &  ! theta
            & initial%asteroid%omega*unit_time, &  ! omega
            & initial%asteroid%elements(1)/unit_dist, &  ! a
            & initial%asteroid%elements(2), &  ! e
            & initial%asteroid%elements(3)/radian, &  ! M
            & initial%asteroid%elements(4)/radian, &  ! w
            & cero, &   ! MMR
            & initial%asteroid%mass/unit_mass, &  ! mass
            & initial%asteroid%radius/unit_dist, &  ! radius
            & actual%asteroid%theta/radian, &  ! theta
            & actual%asteroid%omega*unit_time, &  ! omega
            & actual%asteroid%elements(1)/unit_dist, &  ! a
            & actual%asteroid%elements(2), &  ! e
            & actual%asteroid%elements(3)/radian, &  ! M
            & actual%asteroid%elements(4)/radian, &  ! w
            & cero, &   ! MMR
            & actual%asteroid%mass/unit_mass, &  ! mass
            & actual%asteroid%radius/unit_dist, &  ! radius
            & actual%asteroid%chaos_a/unit_dist, &  ! da
            & actual%asteroid%chaos_e, & ! de
            & cero  ! MEGNO
    end subroutine write_ast_chaos

    ! Write chaos moon i
    subroutine write_moon_i_chaos(initial, actual, i, unit_file)
        implicit none
        type(system_st), intent(in) :: initial
        type(system_st), intent(in) :: actual
        integer(kind=4), intent(in) :: i, unit_file
        integer(kind=4):: i_initial

        do i_initial = 1, initial%Nmoons
            if (actual%moons(i)%id == initial%moons(i_initial)%id) exit
        end do

        if (i_initial == -1) then
            write (*, *) "ERROR: Moon not found in initial system:", actual%moons(i)%id
            stop 2
        end if

        write (unit_file, i3r24) &
            & actual%moons(i)%id, &  ! ID
            & 1, &  ! type
            & actual%moons(i)%merged_to, &  ! merged_to
            & actual%moons(i)%tmax/unit_time, &  ! time
            & initial%asteroid%theta/radian, &  ! theta
            & initial%asteroid%omega*unit_time, &  ! omega
            & initial%moons(i_initial)%elements(1)/unit_dist, &  ! a
            & initial%moons(i_initial)%elements(2), &  ! e
            & initial%moons(i_initial)%elements(3)/radian, &  ! M
            & initial%moons(i_initial)%elements(4)/radian, &  ! w
            & initial%moons(i_initial)%mmr, &   ! MMR
            & initial%moons(i_initial)%mass/unit_mass, &  ! mass
            & initial%moons(i_initial)%radius/unit_dist, &  ! radius
            & actual%asteroid%theta/radian, &  ! theta
            & actual%asteroid%omega*unit_time, &  ! omega
            & actual%moons(i)%elements(1)/unit_dist, &  ! a
            & actual%moons(i)%elements(2), &  ! e
            & actual%moons(i)%elements(3)/radian, &  ! M
            & actual%moons(i)%elements(4)/radian, &  ! w
            & actual%moons(i)%mmr, &   ! MMR
            & actual%moons(i)%mass/unit_mass, &  ! mass
            & actual%moons(i)%radius/unit_dist, &  ! radius
            & actual%moons(i)%chaos_a/unit_dist, &  ! da
            & actual%moons(i)%chaos_e, & ! de
            & cero  ! MEGNO
    end subroutine write_moon_i_chaos

    ! Write chaos particle i
    subroutine write_particle_i_chaos(initial, actual, i, unit_file)
        implicit none
        type(system_st), intent(in) :: initial
        type(system_st), intent(in) :: actual
        integer(kind=4), intent(in) :: i, unit_file
        integer(kind=4) :: i_initial
        real(kind=wp) :: megno
        if (actual%time > tini) then
            megno = dos * actual%particles(i)%megno / actual%time
        else
            megno = cero
        end if

        do i_initial = 1, initial%Nparticles
            if (actual%particles(i)%id == initial%particles(i_initial)%id) exit
        end do

        if (i_initial == -1) then
            write (*, *) "ERROR: Particle not found in initial system:", actual%particles(i)%id
            stop 2
        end if

        write (unit_file, i3r24) &
            & actual%particles(i)%id, &  ! ID
            & 2, &  ! type
            & actual%particles(i)%merged_to, &  ! merged_to
            & actual%particles(i)%tmax/unit_time, &  ! time
            & initial%asteroid%theta/radian, &  ! theta
            & initial%asteroid%omega*unit_time, &  ! omega
            & initial%particles(i_initial)%elements(1)/unit_dist, &  ! a
            & initial%particles(i_initial)%elements(2), &  ! e
            & initial%particles(i_initial)%elements(3)/radian, &  ! M
            & initial%particles(i_initial)%elements(4)/radian, &  ! w
            & initial%particles(i_initial)%mmr, &   ! MMR
            & cero, &  ! mass
            & cero, &  ! radius
            & actual%asteroid%theta/radian, &  ! theta
            & actual%asteroid%omega*unit_time, &  ! omega
            & actual%particles(i)%elements(1)/unit_dist, &  ! a
            & actual%particles(i)%elements(2), &  ! e
            & actual%particles(i)%elements(3)/radian, &  ! M
            & actual%particles(i)%elements(4)/radian, &  ! w
            & actual%particles(i)%mmr, &   ! MMR
            & cero, &  ! mass
            & cero, &  ! radius
            & actual%particles(i)%chaos_a/unit_dist, &  ! da
            & actual%particles(i)%chaos_e, & ! de
            & megno  ! MEGNO
    end subroutine write_particle_i_chaos

    ! Write chaos ALL
    subroutine write_chaos(initial, actual, unit_file)
        implicit none
        type(system_st), intent(in) :: initial
        type(system_st), intent(in) :: actual
        integer(kind=4), intent(in) :: unit_file
        integer(kind=4) :: i
        integer(kind=4), dimension(:), allocatable :: ids

        allocate (ids(max(actual%Nmoons, actual%Nparticles)))

        if (actual%asteroid%chaos_a(2) > myepsilon) call write_ast_chaos(initial, actual, unit_file)

        if (actual%Nmoons > 1) then
            do i = 1, actual%Nmoons
                ids(i) = i
            end do
            call quickargsort_int(actual%moons%id, ids, 1, actual%Nmoons)
            do i = 1, actual%Nmoons
                call write_moon_i_chaos(initial, actual, ids(i), unit_file)
            end do
        else if (actual%Nmoons == 1) then
            call write_moon_i_chaos(initial, actual, 1, unit_file)
        end if

        if (actual%Nparticles > 1) then
            do i = 1, actual%Nparticles
                ids(i) = i
            end do
            call quickargsort_int(actual%particles%id, ids, 1, actual%Nparticles)
            do i = 1, actual%Nparticles
                call write_particle_i_chaos(initial, actual, ids(i), unit_file)
            end do
        else if (actual%Nparticles == 1) then
            call write_particle_i_chaos(initial, actual, 1, unit_file)
        end if

        deallocate (ids)

    end subroutine write_chaos

    ! Write geometric elements moon i
    subroutine write_moon_i_geom(self, i, unit_file)
        implicit none
        type(system_st), intent(in) :: self
        integer(kind=4), intent(in) :: i, unit_file

        write (unit_file, i2r16) &
            & self%moons(i)%id, &  ! ID
            & 1, &  ! type
            & self%time/unit_time, &  ! time
            & self%asteroid%theta/radian, &  ! theta
            & self%asteroid%omega*unit_time, &  ! omega
            & self%moons(i)%geometric(1)/unit_dist, &  ! a
            & self%moons(i)%geometric(2), &  ! e
            & self%moons(i)%geometric(3)/radian, &  ! M
            & self%moons(i)%geometric(4)/radian, &  ! w
            & self%moons(i)%mmr_geom, &  ! MMR
            & self%moons(i)%mass/unit_mass, &  ! mass
            & self%moons(i)%radius/unit_dist, &  ! radius
            & self%moons(i)%dist_to_cm/unit_dist, &  ! distance
            & self%moons(i)%chaos_a_geom/unit_dist, &  ! da
            & self%moons(i)%chaos_e_geom, & ! de
            & cero  ! MEGNO
    end subroutine write_moon_i_geom

    ! Write geometric elements particle i
    subroutine write_particle_i_geom(self, i, unit_file)
        implicit none
        type(system_st), intent(in) :: self
        integer(kind=4), intent(in) :: i, unit_file
        real(kind=wp) :: megno
        if (self%time > tini) then
            megno = dos * self%particles(i)%megno / self%time
        else
            megno = cero
        end if

        write (unit_file, i2r16) &
            & self%particles(i)%id, &   ! ID
            & 2, &  ! type
            & self%time/unit_time, &  ! time
            & self%asteroid%theta/radian, &  ! theta
            & self%asteroid%omega*unit_time, &  ! omega
            & self%particles(i)%geometric(1)/unit_dist, &  ! a
            & self%particles(i)%geometric(2), &  ! e
            & self%particles(i)%geometric(3)/radian, &  ! M
            & self%particles(i)%geometric(4)/radian, &  ! w
            & self%particles(i)%mmr_geom, &  ! MMR
            & cero, &  ! mass
            & cero, &  ! radius
            & self%particles(i)%dist_to_cm/unit_dist, &  ! distance
            & self%particles(i)%chaos_a_geom/unit_dist, &  ! da
            & self%particles(i)%chaos_e_geom, & ! de
            & megno  ! MEGNO
    end subroutine write_particle_i_geom

    ! Write geometric elements ALL
    subroutine write_geom(self, unit_file)
        implicit none
        type(system_st), intent(in) :: self
        integer(kind=4), intent(in) :: unit_file
        integer(kind=4) :: i
        integer(kind=4), dimension(:), allocatable :: ids

        allocate (ids(max(self%Nmoons_active, self%Nparticles_active)))

        if (self%Nmoons_active > 1) then
            do i = 1, self%Nmoons_active
                ids(i) = i
            end do
            call quickargsort_int(self%moons%id, ids, 1, self%Nmoons_active)
            do i = 1, self%Nmoons_active
                call write_moon_i_geom(self, ids(i), unit_file)
            end do
        else if (self%Nmoons_active == 1) then
            call write_moon_i_geom(self, 1, unit_file)
        end if

        if (self%Nparticles_active > 1) then
            do i = 1, self%Nparticles_active
                ids(i) = i
            end do
            call quickargsort_int(self%particles%id, ids, 1, self%Nparticles_active)
            do i = 1, self%Nparticles_active
                call write_particle_i_geom(self, ids(i), unit_file)
            end do
        else if (self%Nparticles_active == 1) then
            call write_particle_i_geom(self, 1, unit_file)
        end if

        deallocate (ids)

    end subroutine write_geom

    ! Write geometric chaos moon i
    subroutine write_moon_i_geomchaos(initial, actual, i, unit_file)
        implicit none
        type(system_st), intent(in) :: initial
        type(system_st), intent(in) :: actual
        integer(kind=4), intent(in) :: i, unit_file
        integer(kind=4):: i_initial

        do i_initial = 1, initial%Nmoons
            if (actual%moons(i)%id == initial%moons(i_initial)%id) exit
        end do

        if (i_initial == -1) then
            write (*, *) "ERROR: Moon not found in initial system:", actual%moons(i)%id
            stop 2
        end if

        write (unit_file, i3r24) &
            & actual%moons(i)%id, &  ! ID
            & 1, &  ! type
            & actual%moons(i)%merged_to, &  ! merged_to
            & actual%moons(i)%tmax/unit_time, &  ! time
            & initial%asteroid%theta/radian, &  ! theta
            & initial%asteroid%omega*unit_time, &  ! omega
            & initial%moons(i_initial)%geometric(1)/unit_dist, &  ! a
            & initial%moons(i_initial)%geometric(2), &  ! e
            & initial%moons(i_initial)%geometric(3)/radian, &  ! M
            & initial%moons(i_initial)%geometric(4)/radian, &  ! w
            & initial%moons(i_initial)%mmr_geom, &   ! MMR
            & initial%moons(i_initial)%mass/unit_mass, &  ! mass
            & initial%moons(i_initial)%radius/unit_dist, &  ! radius
            & actual%asteroid%theta/radian, &  ! theta
            & actual%asteroid%omega*unit_time, &  ! omega
            & actual%moons(i)%geometric(1)/unit_dist, &  ! a
            & actual%moons(i)%geometric(2), &  ! e
            & actual%moons(i)%geometric(3)/radian, &  ! M
            & actual%moons(i)%geometric(4)/radian, &  ! w
            & actual%moons(i)%mmr_geom, &   ! MMR
            & actual%moons(i)%mass/unit_mass, &  ! mass
            & actual%moons(i)%radius/unit_dist, &  ! radius
            & actual%moons(i)%chaos_a_geom/unit_dist, &  ! da
            & actual%moons(i)%chaos_e_geom, & ! de
            & cero  ! MEGNO
    end subroutine write_moon_i_geomchaos

    ! Write geometric chaos particle i
    subroutine write_particle_i_geomchaos(initial, actual, i, unit_file)
        implicit none
        type(system_st), intent(in) :: initial
        type(system_st), intent(in) :: actual
        integer(kind=4), intent(in) :: i, unit_file
        integer(kind=4) :: i_initial
        real(kind=wp) :: megno
        if (actual%time > tini) then
            megno = dos * actual%particles(i)%megno / actual%time
        else
            megno = cero
        end if

        do i_initial = 1, initial%Nparticles
            if (actual%particles(i)%id == initial%particles(i_initial)%id) exit
        end do

        if (i_initial == -1) then
            write (*, *) "ERROR: Particle not found in initial system:", actual%particles(i)%id
            stop 2
        end if

        write (unit_file, i3r24) &
            & actual%particles(i)%id, &  ! ID
            & 2, &  ! type
            & actual%particles(i)%merged_to, &  ! merged_to
            & actual%particles(i)%tmax/unit_time, &  ! time
            & initial%asteroid%theta/radian, &  ! theta
            & initial%asteroid%omega*unit_time, &  ! omega
            & initial%particles(i_initial)%geometric(1)/unit_dist, &  ! a
            & initial%particles(i_initial)%geometric(2), &  ! e
            & initial%particles(i_initial)%geometric(3)/radian, &  ! M
            & initial%particles(i_initial)%geometric(4)/radian, &  ! w
            & initial%particles(i_initial)%mmr_geom, &   ! MMR
            & cero, &  ! mass
            & cero, &  ! radius
            & actual%asteroid%theta/radian, &  ! theta
            & actual%asteroid%omega*unit_time, &  ! omega
            & actual%particles(i)%geometric(1)/unit_dist, &  ! a
            & actual%particles(i)%geometric(2), &  ! e
            & actual%particles(i)%geometric(3)/radian, &  ! M
            & actual%particles(i)%geometric(4)/radian, &  ! w
            & actual%particles(i)%mmr_geom, &   ! MMR
            & cero, &  ! mass
            & cero, &  ! radius
            & actual%particles(i)%chaos_a_geom/unit_dist, &  ! da
            & actual%particles(i)%chaos_e_geom, & ! de
            & megno ! MEGNO
    end subroutine write_particle_i_geomchaos

    ! Write geometric chaos ALL
    subroutine write_geomchaos(initial, actual, unit_file)
        implicit none
        type(system_st), intent(in) :: initial
        type(system_st), intent(in) :: actual
        integer(kind=4), intent(in) :: unit_file
        integer(kind=4) :: i
        integer(kind=4), dimension(:), allocatable :: ids

        allocate (ids(max(actual%Nmoons, actual%Nparticles)))

        if (actual%Nmoons > 1) then
            do i = 1, actual%Nmoons
                ids(i) = i
            end do
            call quickargsort_int(actual%moons%id, ids, 1, actual%Nmoons)
            do i = 1, actual%Nmoons
                call write_moon_i_geomchaos(initial, actual, ids(i), unit_file)
            end do
        else if (actual%Nmoons == 1) then
            call write_moon_i_geomchaos(initial, actual, 1, unit_file)
        end if

        if (actual%Nparticles > 1) then
            do i = 1, actual%Nparticles
                ids(i) = i
            end do
            call quickargsort_int(actual%particles%id, ids, 1, actual%Nparticles)
            do i = 1, actual%Nparticles
                call write_particle_i_geomchaos(initial, actual, ids(i), unit_file)
            end do
        else if (actual%Nparticles == 1) then
            call write_particle_i_geomchaos(initial, actual, 1, unit_file)
        end if

        deallocate (ids)

    end subroutine write_geomchaos

    ! Write diagnostic information ALL
    subroutine write_diagnostics(initial, actual, unit_file)
        type(system_st), intent(in) :: initial
        type(system_st), intent(in) :: actual
        integer(kind=4), intent(in), optional :: unit_file
        integer(kind=4) :: my_file = 6  ! StdOut
        real(wp) :: mass, COM(4), energy, ang_mom, inertia

        call get_cm(actual, mass, COM)
        call calculate_energy_ang_mom_inertia(actual, energy, ang_mom, inertia)
        if (present(unit_file)) my_file = unit_file

        write (my_file, r13) actual%time/unit_time, &  ! time
                         & mass/unit_mass, &  ! mass
                         & mass/initial%mass, &  ! mass / mass0
                         & COM(1:2)/unit_dist, &  ! x y
                         & COM(3:4)/unit_vel, &  ! vx vy
                         & actual%asteroid%omega*unit_time, &  ! omega
                         & actual%asteroid%omega/initial%asteroid%omega, &  ! omega / omega0
                         & energy/unit_ener, &  ! energy
                         & energy/initial%energy, &  ! energy / energy0
                         & ang_mom/unit_angm, &  ! ang_mom
                         & ang_mom/initial%ang_mom  ! ang_mom / ang_mom0
    end subroutine write_diagnostics

    !  --------------------------   COPY    -------------------------------

    ! Copy the objects from another system
    subroutine copy_objects(other, self)
        implicit none
        type(system_st), intent(in) :: other
        type(system_st), intent(inout) :: self

        self%asteroid = other%asteroid
        if (self%Nmoons > 0) then
            self%moons = other%moons
            self%Nmoons_active = other%Nmoons_active
        end if
        if (self%Nparticles > 0) then
            self%particles = other%particles
            self%Nparticles_active = other%Nparticles_active
        end if

        call recalculate_all(self)

    end subroutine copy_objects

end module bodies
