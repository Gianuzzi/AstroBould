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
        real(kind=8) :: theta ! Angle from X (asteroid CM)
        real(kind=8) :: theta_from_primary ! Initial angle from primary(X)  [config]
        real(kind=8) :: dist_to_asteroid  ! Distance to asteroid CM
    end type sphere_st

    type :: asteroid_st
        integer(kind=4) :: Nbodies  ! Amount of boulders + primary
        integer(kind=4) :: Nboulders  ! Amount of boulders
        type(sphere_st), allocatable :: boulders(:)  ! Includes primary  (from 0)
        real(kind=8) :: mass  ! Mass  [config -]
        real(kind=8) :: radius  ! Radius
        real(kind=8) :: theta ! Angle from X [dynamic]
        real(kind=8) :: omega  ! Spin  [config] [dynamic]
        real(kind=8) :: inertia  ! Inertia moment
        real(kind=8) :: ang_mom  ! Angular momentum [dynamic]
        real(kind=8) :: a_corotation  ! Corotation a, for massless
        real(kind=8) :: rotational_period  ! Period of rotation  [config]
        real(kind=8) :: omega_kep  ! Keplerian omega boulders would have
        real(kind=8) :: lambda_kep  ! Ratio from omega to omega_kep  [config] [dynamic]
        real(kind=8) :: dist_to_cm  ! Distance to origin [dynamic]
        real(kind=8), dimension(6) :: elements  ! a, e, i, M, w, O
        real(kind=8), dimension(6) :: coordinates  ! x, y, z, vx, vy, vz
    end type asteroid_st

    type(asteroid_st) :: asteroid

    type :: moon_st
        integer(kind=8) :: id  ! Identifier
        real(kind=8) :: mu_to_asteroid ! Mass ratio to asteroid  [config] 
        real(kind=8) :: mass  ! Masa
        ! real(kind=8) :: radius  ! Radio
        real(kind=8) :: dist_to_cm  ! Distance to origin [dynamic]
        real(kind=8) :: mmr  ! initial mean motion ratio to asteroid  [config] 
        real(kind=8), dimension(6) :: elements  ! a, e, i, M, w, O [config]
        real(kind=8), dimension(6) :: coordinates  ! x, y, z, vx, vy, vz
    end type moon_st
    
    integer(kind=4) :: Nmoons
    type(moon_st), allocatable :: moons(:)
    logical :: use_moons

    type :: particle_st
        integer(kind=8) :: id  ! Identifier
        real(kind=8) :: dist_to_cm  ! Distance to origin
        real(kind=8) :: mmr  ! initial mean motion ratio to asteroid  [config] 
        real(kind=8), dimension(6) :: elements  ! a, e, i, M, w, O  [config]
        real(kind=8), dimension(6) :: coordinates  ! x, y, z, vx, vy, vz
    end type particle_st

    integer(kind=4) :: Nparticles
    type(particle_st), allocatable :: particles(:)
    logical :: use_particles
    
    contains

        subroutine allocate_asteriod(self, Nboulders)
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
        end subroutine allocate_asteriod

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
                use_moons = .True.
            else 
                use_moons = .False.
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
                use_particles = .True.
            else 
                use_particles = .False.
            end if
        end subroutine allocate_particles

        subroutine add_primary(self, mass_primary_or_ast, radius_primary)
            implicit none
            type(asteroid_st), intent(inout) :: self
            real(kind=8), intent(in) :: mass_primary_or_ast, radius_primary
            
            if (self%boulders(0)%id == -1) then
                write(*,*) "ERROR: Primary already set."
                stop 1
            end if

            self%boulders(0)%mass = mass_primary_or_ast
            self%boulders(0)%radius = radius_primary
            self%boulders(0)%mu_to_primary = uno
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

        subroutine add_moon(self, mu_to_asteroid, ele_a, ele_e, ele_M, ele_w, mmr)
            implicit none
            type(moon_st), dimension(:), intent(inout) :: self
            real(kind=8), intent(in) :: mu_to_asteroid, ele_a, ele_e, ele_M, ele_w, mmr
            integer(kind=4) :: i
            logical :: slot_found

            slot_found = .false.

            ! Ensure not massless
            if (mu_to_asteroid < tini) then
                write(*,*) "ERROR: Moons can not have zero or negative mass."
                stop 1
            end if

            ! Look for first empty moon slot (id = -1)
            do i = 1, size(self)
                if (self(i)%id == -1) then
                    self(i)%id = i
                    self(i)%mu_to_asteroid = mu_to_asteroid
                    self(i)%elements(1) = ele_a
                    self(i)%elements(2) = ele_e
                    self(i)%elements(3) = cero  ! i
                    self(i)%elements(4) = ele_M
                    self(i)%elements(5) = ele_w
                    self(i)%elements(6) = cero  ! O
                    self(i)%mmr = mmr
                    slot_found = .true.
                    exit
                end if
            end do

            if (.not. slot_found) then
                write(*,*) "ERROR: No available slot to add a new moon."
                stop 1
            end if
        end subroutine add_moon

        subroutine add_particle(self, ele_a, ele_e, ele_M, ele_w, mmr)
            implicit none
            type(particle_st), dimension(:), intent(inout) :: self
            real(kind=8), intent(in) :: ele_a, ele_e, ele_M, ele_w, mmr
            integer(kind=4) :: i
            logical :: slot_found

            slot_found = .false.

            ! Look for first empty particle slot (id = -1)
            do i = 1, size(self)
                if (self(i)%id == -1) then
                    self(i)%id = i
                    self(i)%elements(1) = ele_a
                    self(i)%elements(2) = ele_e
                    self(i)%elements(3) = cero  ! i
                    self(i)%elements(4) = ele_M
                    self(i)%elements(5) = ele_w
                    self(i)%elements(6) = cero  ! O
                    self(i)%mmr = mmr
                    slot_found = .true.
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
            allocate(pos_from_primary(self%Nbodies,2))
            pos_from_primary = cero
            do i = 1, self%Nboulders
                pos_from_primary(i,1) = self%boulders(0)%radius * cos(self%boulders(i)%theta_from_primary) ! X boulder from primary
                pos_from_primary(i,2) = self%boulders(0)%radius * sin(self%boulders(i)%theta_from_primary) ! Y boulder from primary
            end do
            pos_cm_from_primary(1) = dot_product(pos_from_primary(:,1), self%boulders(:)%mu_to_asteroid)
            pos_cm_from_primary(2) = dot_product(pos_from_primary(:,2), self%boulders(:)%mu_to_asteroid)

            !! Set angles and distances
            !!! Body 0
            if (self%Nboulders > 0) then
                self%boulders(0)%theta = mod(atan2(pos_cm_from_primary(2), pos_cm_from_primary(1)) + pi, twopi)
            else 
                self%boulders(0)%theta = cero
            end if
            self%boulders(0)%dist_to_asteroid = sqrt(pos_cm_from_primary(1) * pos_cm_from_primary(1) + pos_cm_from_primary(2) * pos_cm_from_primary(2))
            !!! Boulders
            do i = 1, self%Nboulders
                aux_real2 = pos_from_primary(i,:) - pos_cm_from_primary
                self%boulders(i)%theta = atan2(aux_real2(2), aux_real2(1))
                self%boulders(i)%dist_to_asteroid = sqrt(aux_real2(1) * aux_real2(1) + aux_real2(2) * aux_real2(2))
            end do

            ! Radius
            self%radius = cero
            do i = 0, self%Nboulders
                self%radius = max(self%radius, self%boulders(i)%radius + self%boulders(i)%dist_to_asteroid) ! Radius
            end do

            ! Rotations
            self%omega_kep = sqrt(G * self%mass / self%boulders(0)%radius**3) ! Keplerian mean motion boulders would have
            if (lambda_kep > tini) then    
                self%lambda_kep = lambda_kep
                self%omega = self%omega_kep * self%lambda_kep  ! Angular velocity of body 0 [rad/day]
                self%rotational_period = twopi / self%omega ! Periodo de rotación del cuerpo 0
            else
                self%rotational_period = rotational_period
                self%omega = (twopi / rotational_period) ! Angular velocity of body 0
                self%lambda_kep = self%omega / self%omega_kep ! Ratio of omegas
            end if
            self%a_corotation = get_a_corot(self%mass, self%omega)

            ! Inertia and angular momentum
            self%inertia = cero
            do i = 0, self%Nboulders
                aux_real = 0.4d0 * self%boulders(i)%mass * self%boulders(i)%radius**2 ! Inertia Sphere
                self%inertia = self%inertia + aux_real + self%boulders(i)%mass * self%boulders(i)%dist_to_asteroid**2 ! Sphere + Steiner
            end do
            self%ang_mom = self%inertia * self%omega ! Rotational


            !! Free memory
            deallocate(pos_from_primary)
        end subroutine init_asteroid_params

        ! Init moon derived params
        subroutine init_moons_params(self, asteroid)
            implicit none
            type(moon_st), dimension(:), intent(inout) :: self
            type(asteroid_st), intent(in) :: asteroid
            real(kind=8) :: combined_mass, aux_real 
            integer(kind=4) :: i

            do i = 1, size(self)
                if (self(i)%id == -1) then  ! Ensure all moons
                    write(*,*) "ERROR: Not all boulders loaded. Missing:", i
                    stop 1
                end if
                self(i)%mass = self(i)%mu_to_asteroid * asteroid%mass

                ! Define coordinates
                !! Elements are ateroid-centric
                combined_mass = asteroid%mass + self(i)%mass
                aux_real = get_a_corot(combined_mass, asteroid%omega)  ! Ask. Including mass?
                if (self(i)%mmr > tini) then
                    self(i)%elements(1) = self(i)%mmr**(2.d0/3.d0) * aux_real  ! a
                else
                    self(i)%mmr = (self(i)%elements(1) / aux_real)**(1.5d0) ! MMR 
                end if
                call coord(combined_mass, &
                            & self(i)%elements(1), &
                            & self(i)%elements(2), &
                            & self(i)%elements(3), &
                            & self(i)%elements(4), &
                            & self(i)%elements(5), &
                            & self(i)%elements(6), &
                            self(i)%coordinates )  ! This are asteroid-centric
            end do
        end subroutine init_moons_params

        ! Init particles derived params
        subroutine init_particles_params(self, asteroid)
            implicit none
            type(particle_st), dimension(:), intent(inout) :: self
            type(asteroid_st), intent(in) :: asteroid
            integer(kind=4) :: i
            do i = 1, size(self)
                if (self(i)%id == -1) then  ! Ensure all moons
                    write(*,*) "ERROR: Not all particles loaded. Missing:", i
                    stop 1
                end if

                ! Define coordinates
                !! Elements are ateroid-centric
                if (self(i)%mmr > tini) then
                    self(i)%elements(1) = self(i)%mmr**(2.d0/3.d0) * asteroid%a_corotation  ! a
                else
                    self(i)%mmr = (self(i)%elements(1) / asteroid%a_corotation)**(1.5d0) ! MMR 
                end if
                call coord(asteroid%mass, &
                            & self(i)%elements(1), &
                            & self(i)%elements(2), &
                            & self(i)%elements(3), &
                            & self(i)%elements(4), &
                            & self(i)%elements(5), &
                            & self(i)%elements(6), &
                            self(i)%coordinates )  ! This are asteroid-centric
            end do
        end subroutine init_particles_params

        subroutine free_asteroid()
            implicit none
            if (allocated(asteroid%boulders)) deallocate(asteroid%boulders)
        end subroutine free_asteroid

        subroutine free_moons()
            implicit none
            if (allocated(moons)) deallocate(moons)
        end subroutine free_moons

        subroutine free_particles()
            implicit none
            if (allocated(particles)) deallocate(particles)
        end subroutine free_particles

    
end module bodies