!> Module with Asteroid, Moons, and Particles structures and routines
module bodies
    use constants
    use celestial
    implicit none

    type :: sphere_st
        integer(kind=8) :: id  ! Identifier
        real(kind=8) :: mu_to_primary ! Mass ratio to primary body [config]
        real(kind=8) :: mu_to_asteroid ! Mass ratio to asteroid 
        real(kind=8) :: mass  ! Mass
        real(kind=8) :: radius  ! Radius  [config]
        real(kind=8) :: initial_theta ! Initial angle from X (asteroid CM)
        real(kind=8) :: theta_from_primary ! Initial angle from primary(X)  [config]
        real(kind=8) :: dist_to_asteroid  ! Distance to asteroid CM
    end type sphere_st

    type :: asteroid_st
        integer(kind=4) :: Nbodies  ! Amount of boulders + primary
        integer(kind=4) :: Nboulders  ! Amount of boulders
        type(sphere_st), allocatable :: boulders(:)  ! Includes primary  (from 0)
        real(kind=8) :: mass  ! Mass  [config -]
        real(kind=8) :: radius  ! Radius
        real(kind=8) :: theta ! Rotational angle from X [dynamic]
        real(kind=8) :: omega  ! Spin  [config] [dynamic]
        real(kind=8) :: a_corotation  ! Corotation a, for massless
        real(kind=8) :: rotational_period  ! Period of rotation  [config]
        real(kind=8) :: omega_kep  ! Keplerian omega boulders would have
        real(kind=8) :: lambda_kep  ! Ratio from omega to omega_kep  [config] [dynamic]
        real(kind=8) :: inertia  ! Inertia moment
        real(kind=8) :: ang_mom_rot  ! Rotational angular momentum [dynamic]
        real(kind=8) :: ang_mom_orb  ! Orbital angular momentum [dynamic]
        real(kind=8) :: dist_to_cm  ! Distance to origin [dynamic]
        real(kind=8), dimension(4) :: elements  ! a, e, M, w
        real(kind=8), dimension(4) :: coordinates  ! x, y, vx, vy
        real(kind=8) :: e_rot ! Rotational energy [dynamic]
        real(kind=8) :: e_kin ! Kinetic energy [dynamic]
        real(kind=8) :: chaos_a(2) ! (a_max, a_min)
        real(kind=8) :: chaos_e(2) ! (e_max, e_min)
    end type asteroid_st


    type :: particle_st
        logical :: active  ! Whether this particle is active (orbiting) or not
        integer(kind=8) :: id  ! Identifier
        real(kind=8) :: dist_to_cm  ! Distance to origin
        real(kind=8), dimension(4) :: elements  ! a, e, M, w  [config]
        real(kind=8), dimension(4) :: coordinates  ! x, y vx, vy
        real(kind=8) :: mmr  ! initial mean motion ratio to asteroid  [config] 
        real(kind=8) :: chaos_a(2) ! (a_max, a_min)
        real(kind=8) :: chaos_e(2) ! (e_max, e_min)
        real(kind=8) :: tmax ! max time integrated
    end type particle_st

    type, extends(particle_st) :: moon_st
        real(kind=8) :: mu_to_asteroid ! Mass ratio to asteroid  [config] 
        real(kind=8) :: mass  ! Masa
        real(kind=8) :: radius  ! Radio
        real(kind=8) :: ang_mom_orb  ! Orbital angular momentum [dynamic]
        real(kind=8) :: ang_mom_rot  ! Rotational angular momentum [dynamic] (Just for conservation)
        real(kind=8) :: e_kin ! Kinetic energy [dynamic]
        real(kind=8) :: e_rot ! Rotational energy [dynamic] (Just for conservation)
    end type moon_st

    ! This struct will store the CM (system) properties and objects
    type :: system_st
        real(kind=8) :: time
        real(kind=8) :: mass
        real(kind=8) :: energy
        real(kind=8) :: ang_mom
        type(asteroid_st) :: asteroid
        integer(kind=4) :: Nmoons
        integer(kind=4) :: Nmoons_active
        type(moon_st), allocatable :: moons(:)
        integer(kind=4) :: Nparticles
        integer(kind=4) :: Nparticles_active
        type(particle_st), allocatable :: particles(:)
    end type system_st

    
    contains

        subroutine allocate_asteroid(self, Nboulders)
            implicit none
            type(asteroid_st), intent(inout) :: self
            integer(kind=4), intent(in) :: Nboulders
            integer(kind=4) :: i

            if (allocated(self%boulders)) then
                write (*,*) "ERROR: Number of boulder already defined."
                STOP 1
            end if

            allocate(self%boulders(0:Nboulders))  ! FROM 0 to Nboulders
            self%Nbodies = Nboulders + 1
            self%Nboulders = Nboulders

            do i = 0, self%Nboulders
                self%boulders(i)%id = -1
            end do
        end subroutine allocate_asteroid

        subroutine allocate_moons(self, Nmoons)
            implicit none
            type(moon_st), allocatable, intent(inout) :: self(:)
            integer(kind=4), intent(in) :: Nmoons
            integer(kind=4) :: i
            
            if (Nmoons > 0) then
                if (allocated(self)) then
                    write (*,*) "ERROR: Number of moons already defined."
                    write (*,*) "        Use moons file or config file; not both."
                    STOP 1
                end if

                allocate(self(Nmoons))
                do i = 1, Nmoons
                    self(i)%id = -1
                end do
            end if
        end subroutine allocate_moons

        subroutine allocate_particles(self, Nparticles)
            implicit none
            type(particle_st), allocatable, intent(inout) :: self(:)
            integer(kind=4), intent(in) :: Nparticles
            integer(kind=4) :: i
            
            if (Nparticles > 0) then
                if (allocated(self)) then
                    write (*,*) "ERROR: Number of particles already defined."
                    write (*,*) "        Use particles file or config file; not both."
                    STOP 1
                end if

                allocate(self(Nparticles))
                do i = 1, Nparticles
                    self(i)%id = -1
                end do
            end if
        end subroutine allocate_particles

        subroutine add_primary(self, mass_primary_or_ast, radius_primary)
            implicit none
            type(asteroid_st), intent(inout) :: self
            real(kind=8), intent(in) :: mass_primary_or_ast, radius_primary
            
            if (self%boulders(0)%id .ne. -1) then
                write(*,*) "ERROR: Primary already set."
                stop 1
            end if

            ! Ensure not anti-radius
            if (radius_primary < epsilon) then
                write(*,*) "ERROR: Primary can not have zero or negative radius."
                stop 1
            end if

            self%boulders(0)%mu_to_primary = uno
            self%boulders(0)%mass = mass_primary_or_ast
            self%boulders(0)%radius = radius_primary
            self%boulders(0)%theta_from_primary = cero
            self%boulders(0)%id = 0
            
        end subroutine add_primary

        subroutine add_boulder(self, mu_to_primary, radius, theta_from_primary)
            implicit none
            type(asteroid_st), intent(inout) :: self
            real(kind=8), intent(in) :: mu_to_primary, radius, theta_from_primary
            integer(kind=4) :: i
            logical :: slot_found

            slot_found = .false.

            ! Ensure not anti-radius
            if (radius < cero) then
                write(*,*) "ERROR: Boulders can not have negative radius. Use negative mass if needed."
                stop 1
            end if

            ! Look for first empty boulder slot (id = -1)
            do i = 1, self%Nboulders
                if (self%boulders(i)%id == -1) then
                    self%boulders(i)%id = i
                    self%boulders(i)%mu_to_primary = mu_to_primary
                    self%boulders(i)%radius = radius
                    self%boulders(i)%theta_from_primary = theta_from_primary
                    slot_found = .true.
                    exit
                end if
            end do

            if (.not. slot_found) then
                write(*,*) "ERROR: No available slot to add a new boulder."
                stop 1
            end if
        end subroutine add_boulder

        subroutine add_moon(self, mu_to_asteroid, ele_a, ele_e, ele_M, ele_w, mmr, radius, id_start)
            implicit none
            type(moon_st), dimension(:), intent(inout) :: self
            real(kind=8), intent(in) :: mu_to_asteroid, ele_a, ele_e, ele_M, ele_w, mmr, radius
            integer(kind=4), intent(in), optional :: id_start
            integer(kind=4) :: i
            integer(kind=4) :: id0 = 0  ! Where to start using IDs from
            logical :: slot_found

            slot_found = .False.

            ! Ensure not massless
            if (mu_to_asteroid < epsilon) then
                write(*,*) "ERROR: Moons can not have zero or negative mass (I think)..."
                stop 1
            end if

            ! Ensure not anti-radius
            if (radius < cero) then
                write(*,*) "ERROR: Moons can not have negative radius."
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
                write(*,*) "ERROR: No available slot to add a new moon."
                stop 1
            end if
        end subroutine add_moon

        subroutine add_particle(self, ele_a, ele_e, ele_M, ele_w, mmr, id_start)
            implicit none
            type(particle_st), dimension(:), intent(inout) :: self
            real(kind=8), intent(in) :: ele_a, ele_e, ele_M, ele_w, mmr
            integer(kind=4), intent(in), optional :: id_start
            integer(kind=4) :: i
            integer(kind=4) :: id0 = 0  ! Where to start using IDs from
            logical :: slot_found

            slot_found = .False.
            
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
                write(*,*) "ERROR: No available slot to add a new particle."
                stop 1
            end if
        end subroutine add_particle

        ! Init asteroid parameters from config
        subroutine init_asteroid_params(self, lambda_kep, rotational_period)
            implicit none
            type(asteroid_st), intent(inout) :: self
            real(kind=8), intent(in) :: lambda_kep, rotational_period
            real(kind=8), allocatable :: pos_from_primary(:,:)
            real(kind=8) :: pos_cm_from_primary(2)
            real(kind=8) :: aux_real, aux_real2(2)
            integer(kind=4) :: i

            ! Ensure all boulders
            do i = 0, self%Nboulders
                if (self%boulders(i)%id == -1) then
                    write(*,*) "ERROR: Not all boulders loaded. Missing:", i
                    stop 1
                end if
            end do
            
            ! Correct Body 0 mass
            aux_real = self%boulders(0)%mass
            if (aux_real < cero) self%boulders(0)%mass = abs(aux_real) / (uno + sum(self%boulders(1:)%mu_to_primary))
            
            ! Boulders mass
            do i = 1, self%Nboulders
                self%boulders(i)%mass = self%boulders(i)%mu_to_primary * self%boulders(0)%mass
            end do

            ! Mass and ratios
            self%mass = sum(self%boulders(:)%mass)  ! Mass
            if (self%mass <= cero) then
                write(*,*) "ERROR: Asteroid mass is zero or negative"
                stop 1
            end if
            self%boulders(:)%mu_to_asteroid = self%boulders(:)%mass / self%mass  ! Mass ratio to asteroid

            ! Calculate distance and angle to CM from boulders
            allocate(pos_from_primary(self%Nboulders,2))
            pos_from_primary = cero
            do i = 1, self%Nboulders
                pos_from_primary(i,1) = self%boulders(0)%radius * cos(self%boulders(i)%theta_from_primary) ! X boulder from primary
                pos_from_primary(i,2) = self%boulders(0)%radius * sin(self%boulders(i)%theta_from_primary) ! Y boulder from primary
            end do
            pos_cm_from_primary(1) = dot_product(pos_from_primary(:,1), self%boulders(1:)%mu_to_asteroid)
            pos_cm_from_primary(2) = dot_product(pos_from_primary(:,2), self%boulders(1:)%mu_to_asteroid)

            !! Set angles and distances
            !!! Body 0
            if (self%Nboulders > 0) then
                self%boulders(0)%initial_theta = mod(atan2(pos_cm_from_primary(2), pos_cm_from_primary(1)) + pi, twopi)
            else 
                self%boulders(0)%initial_theta = cero
            end if
            self%boulders(0)%dist_to_asteroid = sqrt(pos_cm_from_primary(1) * pos_cm_from_primary(1) + &
                                                   & pos_cm_from_primary(2) * pos_cm_from_primary(2))
            !!! Boulders
            do i = 1, self%Nboulders
                aux_real2 = pos_from_primary(i,:) - pos_cm_from_primary
                self%boulders(i)%initial_theta = atan2(aux_real2(2), aux_real2(1))
                self%boulders(i)%dist_to_asteroid = sqrt(aux_real2(1) * aux_real2(1) + aux_real2(2) * aux_real2(2))
            end do

            ! Set Radius
            self%radius = cero
            do i = 0, self%Nboulders
                self%radius = max(self%radius, self%boulders(i)%radius + self%boulders(i)%dist_to_asteroid) ! Radius
            end do

            ! Set Rotations
            self%theta = cero  ! Initial angle with X axis
            self%omega_kep = sqrt(G * self%mass / self%boulders(0)%radius**3) ! Keplerian mean motion boulders would have
            if (abs(lambda_kep) > epsilon) then    
                self%lambda_kep = lambda_kep
                self%omega = self%omega_kep * self%lambda_kep  ! Angular velocity of body 0
                self%rotational_period = twopi / abs(self%omega) ! Rotational period of body 0
            else if (abs(rotational_period) > epsilon) then
                self%rotational_period = abs(rotational_period)
                self%omega = (twopi / rotational_period) ! Angular velocity of body 0
                self%lambda_kep = self%omega / self%omega_kep ! Ratio of omegas
            else  ! No rotation at all
                self%rotational_period = cero
                self%omega = cero ! Angular velocity of body 0
                self%lambda_kep = cero ! Ratio of omegas
            end if
            self%a_corotation = get_a_corot(self%mass, self%omega)

            ! Set Inertia
            self%inertia = cero
            do i = 0, self%Nboulders
                aux_real = 0.4d0 * self%boulders(i)%mass * self%boulders(i)%radius**2 ! Inertia Sphere
                self%inertia = self%inertia + aux_real + self%boulders(i)%mass * self%boulders(i)%dist_to_asteroid**2 ! Sphere + Steiner
            end do

            ! Set Angular momentum
            self%ang_mom_rot = self%inertia * self%omega ! Rotational
            !! Can not set ang_mom_orb here.

            ! Set Energy
            self%e_rot = uno2 * self%inertia * self%omega**2  ! Rotational
            !! Can not set e_kin here. We are assuming v=0


            !! Free memory
            deallocate(pos_from_primary)
        end subroutine init_asteroid_params

        ! Init moon derived parameters
        subroutine init_moons_params(self, asteroid)
            implicit none
            type(moon_st), dimension(:), intent(inout) :: self
            type(asteroid_st), intent(in) :: asteroid
            real(kind=8) :: combined_mass, aux_real 
            real(kind=8) :: aux_real6(6)
            integer(kind=4) :: i

            do i = 1, size(self)
                if (self(i)%id == -1) then  ! Ensure all moons
                    write(*,*) "ERROR: Not all moons loaded. Missing:", i
                    stop 1
                end if
                self(i)%mass = self(i)%mu_to_asteroid * asteroid%mass

                ! Define coordinates
                !! Elements are ateroid-centric
                combined_mass = asteroid%mass + self(i)%mass
                aux_real = get_a_corot(combined_mass, asteroid%omega)  ! Ask. Including mass?
                if (self(i)%mmr > epsilon) then
                    if (aux_real < epsilon) then  ! Ensure rotating
                        write(*,*) "ERROR: Can not set moon MMR with non rotating asteroid."
                        stop 1
                    end if
                    self(i)%elements(1) = self(i)%mmr**(2.d0/3.d0) * aux_real  ! a
                else
                    if (aux_real < epsilon) then  ! Ensure rotating
                        self(i)%mmr = cero
                    else
                        self(i)%mmr = (self(i)%elements(1) / aux_real)**(1.5d0) ! MMR
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
                
                ! Angular momentum
                self(i)%ang_mom_rot = cero
                self(i)%ang_mom_orb = self(i)%mass * (self(i)%coordinates(1) * self(i)%coordinates(4) - &
                                                    & self(i)%coordinates(2) * self(i)%coordinates(3))  ! Z component

                ! Energy
                self(i)%e_rot = cero
                self(i)%e_kin = uno2 * self(i)%mass * (self(i)%coordinates(3) * self(i)%coordinates(3) + &
                                                     & self(i)%coordinates(4) * self(i)%coordinates(4))  ! Kinetic
            end do
        end subroutine init_moons_params

        ! Init particles derived parameters
        subroutine init_particles_params(self, asteroid)
            implicit none
            type(particle_st), dimension(:), intent(inout) :: self
            type(asteroid_st), intent(in) :: asteroid
            integer(kind=4) :: i
            real(kind=8) :: aux_real6(6)
            do i = 1, size(self)
                if (self(i)%id == -1) then  ! Ensure all moons
                    write(*,*) "ERROR: Not all particles loaded. Missing:", i
                    stop 1
                end if

                ! Define coordinates
                !! Elements are ateroid-centric
                if (self(i)%mmr > epsilon) then
                    if (asteroid%a_corotation < epsilon) then  ! Ensure rotating
                        write(*,*) "ERROR: Can not set particle MMR with non rotating asteroid."
                        stop 1
                    end if
                    self(i)%elements(1) = self(i)%mmr**(2.d0/3.d0) * asteroid%a_corotation  ! a
                else
                    if (asteroid%a_corotation < epsilon) then  ! Ensure rotating
                        self(i)%mmr = cero
                    else
                        self(i)%mmr = (self(i)%elements(1) / asteroid%a_corotation)**(1.5d0) ! MMR 
                    end if
                end if
                call coord(asteroid%mass, &
                            & self(i)%elements(1), &
                            & self(i)%elements(2), &
                            & cero, &  ! inc
                            & self(i)%elements(3), &
                            & self(i)%elements(4), &
                            & cero, &  ! Omega
                            & aux_real6 )  ! This are asteroid-centric
                self(i)%coordinates = (/aux_real6(1:2), aux_real6(4:5)/)
            end do
        end subroutine init_particles_params

        ! Update asteroid parameters according to new rotation
        subroutine spin_asteroid(self, theta, omega)
            implicit none
            type(asteroid_st), intent(inout) :: self
            real(kind=8), intent(in) :: theta
            real(kind=8), intent(in) :: omega

            ! Update Rotation
            self%theta = theta  ! Update THETA
            if (abs(omega) > epsilon) then
                self%a_corotation = get_a_corot(self%mass, omega)
                self%rotational_period = twopi / omega ! Rotational period of body 0
                self%lambda_kep = omega / self%omega_kep ! Ratio of omegas
                self%omega = omega  ! Update OMEGA
            else 
                self%a_corotation = cero
                self%rotational_period = cero ! Rotational period of body 0
                self%lambda_kep = cero  ! Ratio of omegas
                self%omega = cero ! Update OMEGA
            end if

            ! Update Angular momentum
            self%ang_mom_rot = self%inertia * omega ! Rotational

            ! Update Energy
            self%e_rot = uno2 * self%inertia * omega**2  ! Rotational
        end subroutine spin_asteroid

        ! Update asteroid parameters according to new coordinates
        subroutine shift_asteroid(self, coordinates, total_mass)
            implicit none
            type(asteroid_st), intent(inout) :: self
            real(kind=8), intent(in) :: coordinates(4)  ! Barycentric
            real(kind=8), intent(in), optional :: total_mass
            real(kind=8) :: total_mass_
            real(kind=8) :: a, e, i, M, w, O

            ! Update Coords and distance
            self%coordinates = coordinates ! Update COORDINATES
            self%dist_to_cm = sqrt(coordinates(1) * coordinates(1) + &
                                 & coordinates(2) * coordinates(2))

            ! Update Elements
            !! Get total system mass
            if(present(total_mass)) then 
                if (total_mass < self%mass) then
                    write(*,*) "ERROR: Total mass can not be lower than asteroid mass."
                    stop 1
                end if
                total_mass_ = total_mass
            else 
                total_mass_ = self%mass  ! Assume only the asteroid
            end if
            if (self%dist_to_cm > epsilon) then
                call elem(total_mass_, (/coordinates(1:2), cero, coordinates(3:4), cero/), a, e, i, M, w, O)
                self%elements(1) = a
                self%elements(2) = e
                self%elements(3) = M
                self%elements(4) = w
            end if

            ! Update Angular momentum
            self%ang_mom_orb = self%mass * (coordinates(1) * coordinates(4) - &
                                          & coordinates(2) * coordinates(3))  ! Z component

            ! Update Energy
            self%e_kin = uno2 * self%mass * (coordinates(3) * coordinates(3) + &
                                           & coordinates(4) * coordinates(4))  ! Kinetic
        end subroutine shift_asteroid

        ! Update asteroid parameters according to new coordinates
        subroutine grow_asteroid(self, mass_to_add)
            implicit none
            type(asteroid_st), intent(inout) :: self
            real(kind=8), intent(in) :: mass_to_add
            real(kind=8) :: growth
            integer(kind=4) :: i

            if (mass_to_add < epsilon) then
                write(*,*) " No mass to add."
                return  ! No mass to add
            end if

            ! Get mass growth ratio
            growth = uno + (mass_to_add / self%mass)

            !! Update masses [mass ratios are kept the same]
            do i = 0, self%Nboulders
                self%boulders(i)%mass = self%boulders(i)%mass * growth
            end do
            self%mass = self%mass * growth  ! Update MASS
            
            ! Update derived rotation parameters
            self%omega_kep = sqrt(G * self%mass / self%boulders(0)%radius**3)  ! Check
            self%a_corotation = get_a_corot(self%mass, self%omega)  ! Corotation a
            if (abs(self%omega) > cero) self%lambda_kep = self%omega / self%omega_kep ! Ratio of omegas

            ! Update Inertia
            self%inertia = self%inertia * growth

            ! Update angular momentum
            self%ang_mom_rot = self%ang_mom_rot * growth ! Rotational
            self%ang_mom_orb = self%ang_mom_orb * growth ! Orbital

            ! Update Energy
            self%e_rot = self%e_rot * growth  ! Rotational
            self%e_kin = self%e_kin * growth  ! Kinetic
        end subroutine grow_asteroid

        ! Update a single moon parameters according to new coordinates
        subroutine shift_single_moon(self, coordinates, total_mass)
            implicit none
            type(moon_st), intent(inout) :: self
            real(kind=8), intent(in) :: coordinates(4)  ! These are barycentric now
            real(kind=8), intent(in) :: total_mass
            real(kind=8) :: a, e, i, M, w, O

            ! Update coordinates
            self%coordinates = coordinates  ! Update COORDINATES

            !! Distance to CM
            self%dist_to_cm = sqrt(coordinates(1) * coordinates(1) + &
                                 & coordinates(2) * coordinates(2))
            
            ! Get elements (barycentric)
            call elem(total_mass, (/coordinates(1:2), cero, coordinates(3:4), cero/), a, e, i, M, w, O)
            self%elements(1) = a
            self%elements(2) = e
            self%elements(3) = M
            self%elements(4) = w

            ! Angular momentum
            self%ang_mom_orb = self%mass * (self%coordinates(1) * self%coordinates(4) - &
                                          & self%coordinates(2) * self%coordinates(3))  ! Z component

            ! Energy
            self%e_kin = uno2 * self%mass * (coordinates(3) * coordinates(3) + &
                                           & coordinates(4) * coordinates(4))  ! Kinetic
        end subroutine shift_single_moon

        ! Update a single particle parameters according to new coordinates
        subroutine shift_single_particle(self, coordinates, total_mass)
            implicit none
            type(particle_st), intent(inout) :: self
            real(kind=8), intent(in) :: coordinates(4)  ! These are barycentric now
            real(kind=8), intent(in) :: total_mass
            real(kind=8) :: a, e, i, M, w, O

            ! Update coordinates
            self%coordinates = coordinates  ! Update COORDINATES

            !! Distance to CM
            self%dist_to_cm = sqrt(coordinates(1) * coordinates(1) + &
                                 & coordinates(2) * coordinates(2))
            
            ! Get elements (barycentric)
            call elem(total_mass, (/coordinates(1:2), cero, coordinates(3:4), cero/), a, e, i, M, w, O)
            self%elements(1) = a
            self%elements(2) = e
            self%elements(3) = M
            self%elements(4) = w
        end subroutine shift_single_particle

        ! Get barycentric coordinates (position and velocities) of boulders
        subroutine get_boulder_i_coord(self, i, coordinates)
            implicit none
            type(asteroid_st), intent(in) :: self
            integer(kind=4), intent(in) :: i
            real(kind=8), dimension(4), intent(out) :: coordinates
            
            coordinates(1) = self%boulders(i)%dist_to_asteroid * cos(self%theta + self%boulders(i)%initial_theta)  ! x
            coordinates(2) = self%boulders(i)%dist_to_asteroid * sin(self%theta + self%boulders(i)%initial_theta)  ! y
            coordinates(3) = -self%omega * coordinates(2)  ! vx
            coordinates(4) = self%omega * coordinates(1)  ! vy
            
            coordinates = coordinates + self%coordinates  ! Shift to asteroid
        end subroutine get_boulder_i_coord

        ! Recalculate system mass and CM
        subroutine get_cm(self, mass, coordinates)
            implicit none
            type(system_st), intent(in) :: self
            real(kind=8), dimension(4), intent(out) :: coordinates
            real(kind=8), intent(out) :: mass
            real(kind=8) :: aux_4(4) = cero
            integer(kind=4) :: i

            ! Calculate cm
            mass = cero
            coordinates = cero
            do i = 0, self%asteroid%Nboulders
                call get_boulder_i_coord(self%asteroid, i, aux_4)
                coordinates = coordinates + aux_4 * self%asteroid%boulders(i)%mass
                mass = mass + self%asteroid%boulders(i)%mass
            end do
            do i = 1, self%Nmoons_active
                coordinates = coordinates + self%moons(i)%coordinates * self%moons(i)%mass
                mass = mass + self%moons(i)%mass
            end do
            coordinates = coordinates / mass
        end subroutine get_cm

        ! Calculate system energy and angular momentum
        subroutine calculate_energy_and_ang_mom(self, energy, ang_mom)
            implicit none
            type(system_st), intent(in) :: self
            real(kind=8), intent(out) :: energy, ang_mom
            real(kind=8) :: e_pot
            integer(kind=4) :: i, j
            real(kind=8) :: dr(2), dist

            ! Asteroid energy and ang_mom is used as a rigid body

            ! Calculate energy and ang_mom
            ang_mom = self%asteroid%ang_mom_rot + self%asteroid%ang_mom_orb
            energy = self%asteroid%e_kin + self%asteroid%e_rot ! e_pot below
            print*, "Ast K | R:", self%asteroid%e_kin, self%asteroid%e_rot

            !! e_pot
            e_pot = cero
            !! Asteroid (with moons)
            do i = 1, self%Nmoons_active
                dr = self%moons(i)%coordinates(1:2) - self%asteroid%coordinates(1:2)
                dist = sqrt(dr(1) * dr(1) + dr(2) * dr(2))
                if (dist < epsilon) cycle
                e_pot = e_pot - (self%asteroid%mass * self%moons(i)%mass) / dist
            end do
            !! Moons (here, only moons with moons are added)
            do i = 1, self%Nmoons_active - 1
                do j = 2, self%Nmoons_active
                    dr = self%moons(j)%coordinates(1:2) - self%moons(i)%coordinates(1:2)
                    dist = sqrt(dr(1) * dr(1) + dr(2) * dr(2))
                    if (dist < epsilon) cycle
                    e_pot = e_pot - (self%moons(i)%mass * self%moons(j)%mass) / dist
                end do
                energy = energy + self%moons(i)%e_kin + self%moons(i)%e_rot  ! Energy (TBD: Last moon)
                ang_mom = ang_mom + self%moons(i)%ang_mom_orb + self%moons(i)%ang_mom_rot  ! Angular Momentum (TBD: Last moon)
                print*, "Moon K | R:", self%moons(i)%e_kin,  self%moons(i)%e_rot
            end do

            energy = energy + e_pot * G  ! G

            if (self%Nmoons_active > 0) then
                ang_mom = ang_mom + self%moons(self%Nmoons_active)%ang_mom_orb + self%moons(self%Nmoons_active)%ang_mom_rot  ! Last Moon
                energy = energy + self%moons(self%Nmoons_active)%e_kin + self%moons(self%Nmoons_active)%e_rot  ! Last Moon
                print*, "Moon K | R:", self%moons(self%Nmoons_active)%e_kin,  self%moons(self%Nmoons_active)%e_rot
            end if
            print*, "Ep:",  e_pot * G
            print*, "Tot:",  energy
        end subroutine calculate_energy_and_ang_mom

        ! (Re)calculate all system main parameters
        subroutine recalculate_all(self)
            implicit none
            type(system_st), intent(inout) :: self
            real(kind=8) :: total_mass, energy, ang_mom
            real(kind=8) :: rvcm(4), coords_shift(4)
            integer(kind=4) :: i

            ! Get mass and then CM too
            call get_cm(self, total_mass, rvcm)

            self%mass = total_mass  ! Update MASS

            ! Set system at CM, shifting all objects.
            coords_shift = self%asteroid%coordinates - rvcm
            call shift_asteroid(self%asteroid, coords_shift, total_mass)
            do i = 1, self%Nmoons_active
                coords_shift = self%moons(i)%coordinates - rvcm
                call shift_single_moon(self%moons(i), coords_shift, total_mass)
            end do
            do i = 1, self%Nparticles_active
                coords_shift = self%particles(i)%coordinates - rvcm
                call shift_single_particle(self%particles(i), coords_shift, total_mass)
            end do

            ! Set Energy and Angular Momentum
            call calculate_energy_and_ang_mom(self, energy, ang_mom)
            self%ang_mom = ang_mom  !Set ANG_MOM
            self%energy = energy  !Set ENERGY
        end subroutine recalculate_all

        ! Init whole system
        subroutine init_system(self, asteroid, moons, particles, lambda_kep, rotational_period, time)
            implicit none
            type(system_st), intent(inout) :: self
            type(asteroid_st), intent(inout) :: asteroid
            type(moon_st), allocatable, intent(inout) :: moons(:)
            type(particle_st), allocatable, intent(inout) :: particles(:)
            real(kind=8), intent(in), optional :: time
            real(kind=8), intent(in) :: lambda_kep, rotational_period

            ! Initial time. Set TIME
            self%time = cero
            if (present(time)) self%time = time

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

            call recalculate_all(self)
        end subroutine init_system

        ! Swap between two moons, using the positional index
        subroutine swap_moons(self, i, j)
            implicit none
            type(moon_st), intent(inout) :: self(:)
            integer, intent(in) :: i, j
            type(moon_st) :: tmp_moon

            ! Safety check
            if (i < 1 .or. j < 1 .or. i > size(self) .or. j > size(self)) then
                write(*,*) "Error: indices out of range in swap_moons"
                return
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
            integer, intent(in) :: i, j
            type(particle_st) :: tmp_particles

            ! Safety check
            if (i < 1 .or. j < 1 .or. i > size(self) .or. j > size(self)) then
                write(*,*) "Error: indices out of range in swap_particles"
                return
            end if

            ! Perform swap
            tmp_particles = self(i)
            self(i) = self(j)
            self(j) = tmp_particles
        end subroutine swap_particles

        ! Check if syst is in CM = 0
        subroutine check_coordinates(self)
            implicit none
            type(system_st), intent(in) :: self
            real(kind=8), dimension(4) :: coordinates
            real(kind=8) :: mass
            real(kind=8) :: aux_real

            ! Calculate cm
            call get_cm(self, mass, coordinates)

            ! Check MASS
            if (abs(mass - self%mass)/mass > epsilon) then
                write(*,*) "ERROR: Total mass differs from CM mass.", mass, self%mass !abs(mass - self%mass)/mass
                stop 1
            end if

            ! Check POS
            aux_real = sqrt(coordinates(1) * coordinates(1) + &
                          & coordinates(2) * coordinates(2))
            if (aux_real > epsilon) then
                write(*,*) "ERROR: CM is not centered at 0.", aux_real
                stop 1
            end if

            ! Check VEL
            aux_real = sqrt(coordinates(3) * coordinates(3) + &
                          & coordinates(4) * coordinates(4))
            if (aux_real > epsilon) then
                write(*,*) "ERROR: CM velocity is not 0.", aux_real
                stop 1
            end if
        end subroutine check_coordinates

        ! Remove a moon from a system by deactivating and moving it to bottom
        subroutine deactivate_moon_i(self, i)
            implicit none
            type(system_st), intent(inout) :: self
            integer(kind=4), intent(in) :: i
            integer(kind=4) :: j

            ! Check
            if (.not. self%moons(i)%active)  then
                write(*,*) "ERROR: Can not deactivate moon already deactivated:", i
                stop 1
            end if

            ! Deactivate moon
            self%moons(i)%active = .False.

            ! If it is the last active, do not swap
            if (.not. (i .eq. self%Nmoons_active)) then
                ! Switch deactivated moon with last active
                do j = self%Nmoons_active, 1, -1
                    if (i .eq. j) cycle
                    if (self%moons(j)%active) then
                        call swap_moons(self%moons, i, j)
                        exit
                    end if
                end do
            end if
            self%Nmoons_active = self%Nmoons_active - 1  ! Update NMOONS_ACTIVE
        end subroutine deactivate_moon_i

        ! Remove a particle from a system by deactivating and moving it to bottom
        subroutine deactivate_particle_i(self, i)
            implicit none
            type(system_st), intent(inout) :: self
            integer(kind=4), intent(in) :: i
            integer(kind=4) :: j

            ! Check
            if (.not. self%particles(i)%active)  then
                write(*,*) "ERROR: Can not deactivate particle already deactivated:", i
                stop 1
            end if

            ! Deactivate moon
            self%particles(i)%active = .False.

            ! Switch deactivated moon with last active
            do j = self%Nparticles_active, 1, -1
                if (i .eq. j) cycle
                if (self%particles(j)%active) then
                    call swap_particles(self%particles, i, j)
                    exit
                end if
            end do
            self%Nparticles_active = self%Nparticles_active - 1  ! Update NPARTICLES_ACTIVE
        end subroutine deactivate_particle_i

        ! Merge a moon into asteroid. NO DEACTIVATION HERE
        subroutine merge_moon_into_ast(asteroid, moon, total_mass)
            implicit none
            type(asteroid_st), intent(inout) :: asteroid
            type(moon_st), intent(inout) :: moon
            real(kind=8), intent(in), optional :: total_mass
            real(kind=8) :: total_mass_
            real(kind=8) :: m_cm, rv_cm(4)
            real(kind=8) :: l_rot
            real(kind=8) :: aux_real24(2,4)

            ! Check if moon is active
            if (.not. moon%active) then
                write(*,*) "ERROR: Can not merge an inactive moon into the asteroid:", moon%id
                stop 1
            end if

            ! Get these 2 bodies CM properties
            m_cm = asteroid%mass + moon%mass
            rv_cm = (asteroid%mass * asteroid%coordinates + moon%mass * moon%coordinates) / m_cm

            ! Calculate L_rot from this CM
            aux_real24(1,:) = asteroid%coordinates - rv_cm
            aux_real24(2,:) = moon%coordinates - rv_cm
            l_rot = (aux_real24(1,1) * aux_real24(1,4) - aux_real24(1,2) * aux_real24(1,3)) * asteroid%mass + &
                  & (aux_real24(2,1) * aux_real24(2,4) - aux_real24(2,2) * aux_real24(2,3)) * moon%mass

            ! Update MASS and derivates
            call grow_asteroid(asteroid, moon%mass)

            ! Update ROTATION  (Remember the moon can have ang_mom)
            l_rot = l_rot + moon%ang_mom_rot
            call spin_asteroid(asteroid, asteroid%theta, l_rot / asteroid%inertia)

            ! Update COORDINATES and derivates
            !! Get total mass if given
            total_mass_ = m_cm
            if (present(total_mass)) total_mass_ = total_mass
            call shift_asteroid(asteroid, rv_cm, total_mass_)
        end subroutine merge_moon_into_ast

        ! Merge a moon(i) into asteroid, in a system. DEACTIVATION HERE
        subroutine merge_moon_i_into_ast(self, i)
            implicit none
            type(system_st), intent(inout) :: self
            integer(kind=4), intent(in) :: i
            real(kind=8) :: aux_real, ang_mom
            
            ! Perform the merge
            call merge_moon_into_ast(self%asteroid, self%moons(i), self%mass)

            ! Deactivate the moon
            call deactivate_moon_i(self, i)

            ! Check CM
            call check_coordinates(self)

            ! Check ANGULAR MOMENTUM
            call calculate_energy_and_ang_mom(self, aux_real, ang_mom)
            aux_real = max(epsilon, abs(ang_mom))
            if (abs(ang_mom - self%ang_mom)/aux_real > epsilon) then
                write(*,*) "ERROR: Total angular momentum differs from previously calculated."
                stop 1
            end if

            ! ENERGY loss ocurrs, and the amount is given by:
            !! dE = G m1 m2 / r - (m1 m2 / (m1 + m2)) (v2 - v1)**2 / 2
        end subroutine merge_moon_i_into_ast

        ! Eject a moon(i), in a system. DEACTIVATION HERE
        subroutine eject_moon_i(self, i)
            implicit none
            type(system_st), intent(inout) :: self
            integer(kind=4), intent(in) :: i
            
            ! Deactivate the moon
            call deactivate_moon_i(self, i)

            ! Update almost all
            call recalculate_all(self)
        end subroutine eject_moon_i

        ! Merge a moon(i) with a moon(j), in a system. DEACTIVATION HERE
        subroutine merge_2_moons(self, i, j)
            implicit none
            type(system_st), intent(inout) :: self
            integer(kind=4), intent(in) :: i, j
            real(kind=8) :: m_cm, rv_cm(4)
            real(kind=8) :: l_rot, new_inert
            real(kind=8) :: aux_real, ang_mom, aux_real24(2,4)
            
            ! Check if moons ar active
            if ((.not. self%moons(i)%active) .or. ((.not. self%moons(j)%active))) then
                write(*,*) "ERROR: Can not merge moons if one is inactive. i, j:", self%moons(i)%active, self%moons(j)%active
                stop 1
            end if

            ! Get these 2 moons CM properties
            m_cm = self%moons(i)%mass + self%moons(j)%mass
            rv_cm = (self%moons(i)%mass * self%moons(i)%coordinates + &
                   & self%moons(j)%mass * self%moons(j)%coordinates) / m_cm

            ! Calculate L_rot from this CM
            aux_real24(1,:) = self%moons(i)%coordinates - rv_cm
            aux_real24(2,:) = self%moons(j)%coordinates - rv_cm
            l_rot = (aux_real24(1,1) * aux_real24(1,4) - aux_real24(1,2) * aux_real24(1,3)) * self%moons(i)%mass + &
                  & (aux_real24(2,1) * aux_real24(2,4) - aux_real24(2,2) * aux_real24(2,3)) * self%moons(j)%mass

            ! Update MASS and derivates
            self%moons(i)%mass = m_cm
            self%moons(i)%mu_to_asteroid = m_cm / self%asteroid%mass

            ! Update RADIUS
            self%moons(i)%radius = max(self%moons(i)%radius, self%moons(j)%radius)

            ! Update ENERGY and ANGULAR MOMENTUM
            new_inert = 0.4d0 * m_cm * self%moons(i)%radius * self%moons(i)%radius
            self%moons(i)%ang_mom_rot = self%moons(i)%ang_mom_rot + l_rot
            self%moons(i)%e_rot = self%moons(i)%e_rot + (uno2 * l_rot * l_rot / new_inert)
            
            ! Deactivate the moon
            call deactivate_moon_i(self, j)
            
            ! Update COORDINATES and derivates
            call shift_single_moon(self%moons(i), rv_cm, self%mass)

            ! Check CM
            call check_coordinates(self)

            ! Check ANGULAR MOMENTUM
            call calculate_energy_and_ang_mom(self, aux_real, ang_mom)
            aux_real = max(epsilon, abs(ang_mom))
            if (abs(ang_mom - self%ang_mom)/aux_real > epsilon) then
                write(*,*) "ERROR: Total angular momentum differs from previously calculated."
                stop 1
            end if

            ! ENERGY loss ocurrs, and the amount is given by:
            !! dE = G m1 m2 / r - (m1 m2 / (m1 + m2)) (v2 - v1)**2 / 2
        end subroutine merge_2_moons

        ! Update bodies chaos at system
        subroutine update_chaos(self)
            implicit none
            type(system_st), intent(inout) :: self
            real(kind=8) :: aux_real2(2)
            integer(kind=4) :: i

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

        ! Update system state. ! Assumes x,y,vx,vy 
        subroutine update_system_from_pos(self, time, theta, omega, cordinates_arr)
            implicit none
            type(system_st), intent(inout) :: self
            real(kind=8), intent(in) :: time
            real(kind=8), intent(in) :: theta, omega
            real(kind=8), dimension(:), intent(in) :: cordinates_arr
            real(kind=8), dimension(4) :: aux_coordinates
            real(kind=8) :: energy, ang_mom
            integer(kind=4) :: i, idx

            ! Update TIME
            self%time = time

            ! Update Asteroid rotation
            call spin_asteroid(self%asteroid, theta, omega)

            ! Update Asteroid Position and derived parameters
            aux_coordinates = cordinates_arr(1:4) !
            call shift_asteroid(self%asteroid, aux_coordinates, self%mass)

            ! Update Moons Position and derived parameters
            idx = 4  ! First 4 were the asteroid
            do i = 1, self%Nmoons_active
                aux_coordinates = cordinates_arr(idx:idx+3)
                call shift_single_moon(self%moons(i), aux_coordinates, self%mass)
                idx = idx + 4
            end do

            ! Update Particles Position and derived parameters
            do i = 1, self%Nparticles_active
                aux_coordinates = cordinates_arr(idx:idx+3)
                call shift_single_particle(self%particles(i), aux_coordinates, self%mass)
                idx = idx + 4
            end do

            ! Update CHAOS
            call update_chaos(self)
            
            ! Check CM
            call check_coordinates(self)

            ! Update energy and ang_mom
            call calculate_energy_and_ang_mom(self, energy, ang_mom)
            self%energy = energy
            self%ang_mom = ang_mom
        end subroutine update_system_from_pos

        ! Create mass and coordinates
        subroutine generate_mass_and_coordinates_arrays(self, mass_arr, coordinates_arr)
            implicit none
            type(system_st), intent(in) :: self
            real(kind=8), dimension(:), intent(inout) :: mass_arr
            real(kind=8), dimension(:), intent(inout) :: coordinates_arr
            integer(kind=4) :: i, idx

            ! Set asteroid values
            mass_arr(1) = self%asteroid%mass
            coordinates_arr(1:4) = self%asteroid%coordinates

            ! Set Moons values
            idx = 5  ! First 4 were the asteroid
            do i = 1, self%Nmoons_active
                mass_arr(i+1) = self%moons(i)%mass
                coordinates_arr(idx:idx+3) = self%moons(i)%coordinates
                idx = idx + 4
            end do

            ! Set particles values
            do i = 1, self%Nparticles
                coordinates_arr(idx:idx+3) = self%particles(i)%coordinates
                idx = idx + 4
            end do
        end subroutine generate_mass_and_coordinates_arrays

        ! Calculate acceleration and potential at a given coordinate  (Missing G)
        subroutine get_acc_and_pot_xy(self, xy_target, acc, pot, surface)
            implicit none
            type(system_st), intent(in) :: self
            real(kind=8), intent(in) :: xy_target(2)
            real(kind=8), intent(out) :: acc(2), pot
            real(kind=8), intent(in), optional :: surface
            real(kind=8) :: aux_real4(4) = cero
            integer(kind=4) :: i
            logical :: has_surface = .False., inside = .False.
            acc = cero
            pot = cero

            has_surface = present(surface)
            ! Asteroid
            do i = 0, self%asteroid%Nboulders
                inside = .False.
                call get_boulder_i_coord(self%asteroid, i, aux_real4)
                call get_acc_and_pot_single(self%asteroid%boulders(i)%mass, &
                                          & (/aux_real4(1), aux_real4(2)/), &
                                          & xy_target, self%asteroid%boulders(i)%radius, &
                                          & acc, pot, &
                                          inside)
                if (inside .and. has_surface) then
                    acc = (/cero, cero/)
                    pot = surface
                    return
                end if
            end do

            ! Moons
            do i = 1, self%Nmoons_active
                call get_acc_and_pot_single(self%moons(i)%mass, &
                                          & (/self%moons(i)%coordinates(1), self%moons(i)%coordinates(2)/), &
                                          & xy_target, self%moons(i)%radius, &
                                          & acc, pot, &
                                          inside)
                if (inside .and. has_surface) then
                    pot = cero
                    call get_acc_and_pot_single(self%moons(i)%mass, &
                                          & (/self%moons(i)%coordinates(1), self%moons(i)%coordinates(2)/), &
                                          & (/self%moons(i)%coordinates(1) + self%moons(i)%radius, &
                                          &   self%moons(i)%coordinates(2)/), &
                                          & tini, &
                                          & acc, pot)
                    acc = cero
                    return
                end if
            end do
        end subroutine get_acc_and_pot_xy

        subroutine free_asteroid(self)
            implicit none
            type(asteroid_st), intent(inout) :: self
            if (allocated(self%boulders)) deallocate(self%boulders)
        end subroutine free_asteroid

        subroutine free_system(self)
            implicit none
            type(system_st), intent(inout) :: self
            call free_asteroid(self%asteroid)
            if (allocated(self%moons)) deallocate(self%moons)
            if (allocated(self%particles)) deallocate(self%particles)
        end subroutine free_system
end module bodies