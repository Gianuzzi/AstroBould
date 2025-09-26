program main
    use parameters
    use integrators
    use celestial, only: get_Period
    use bodies
    use times
    use derivates
    implicit none
    external :: create_map   ! Declare the external function for map
    logical :: only_potential_map = .False.
    integer(kind=4) :: aux_int
    character(1) :: aux_character1
    character(20) :: aux_character20
    character(30) :: aux_character30
    logical :: aux_logical
    real(kind=8) :: aux_real
    real(kind=8), dimension(:,:), allocatable :: aux_2D  ! To read files
    integer(kind=4) :: i, j  ! Loops
    integer(kind=4) :: new_Nactive = 0,  y_nvalues = 0
    integer(kind=4) :: unit_file = -1  ! Where to write escapes/collisions
    
    use_configfile = .True. ! Usar archivo de configuración
    arguments_number = command_argument_count()
    do i = 1, arguments_number
        call get_command_argument(i, aux_character30)
        if (trim(aux_character30) .eq. "--noconfig") then
            use_configfile = .False.
            exit
        end if
    end do

    if (use_configfile) call read_config_file(input, "config.ini", configfile_exists) ! Leemos archivo de configuración

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!! Editar aquí abajo de ser necesario !!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (.not. use_configfile) then ! Usaremos parámetros por defecto
        
        ! Asteroide central
        !! Primary mass
        input%mass_primary = 6.3d18 ! Masa del cuerpo 0 [kg] ! -x =>  mAst = x

        !! Primary shape

        !!! Tri_Axial
        input%use_triaxial = .False.  ! Logical para determinar si se usa triax, o Boulders
        input%triax_a_primary = cero  ! Semieje a
        input%triax_b_primary = cero  ! Semieje b
        input%triax_c_primary = cero  ! Semieje c

        !!! Explicit primary radius
        input%radius_primary = 129.d0 ! Radio del cuerpo 0 [km]

        !! ROTACIÓN
        input%asteroid_rotational_period = 7.004d0/24.d0  ! Periodo de rotación [day]
        !lambda_kep = 0.471d0      ! Cociente omega/wk

        !! Boulders        
        input%Nboulders = 1 ! Número de boulders

        if (input%Nboulders > 0) then
            call allocate_params_asteroid(input%Nboulders)!! Alocatamos (No tocar)

            boulders_in(1,1) = 1.d-1   ! Cociente de masas entre boulder 1 y primary
            boulders_in(1,2) = cero    ! Ángulo de fase del boulder 1 [deg]
            boulders_in(1,3) = 2.5d0   ! Radio del boulder 1 [km]       

        end if

        !! Moons
        input%Nmoons = 0 ! Número de boulders

        if (input%Nmoons > 0) then
            call allocate_params_moons(input%Nmoons)!! Alocatamos (No tocar)

            moons_in(1,1) = 1.d-6   ! Cociente de masas entre luna 1 y asteroide
            moons_in(1,2) = cero    ! Semieje [km]
            moons_in(1,3) = cero    ! Eccentricidad
            moons_in(1,4) = cero    ! M [deg]
            moons_in(1,5) = cero    ! w [deg]
            moons_in(1,6) = 11.1d0  ! MMR
            moons_in(1,7) = cero    ! radius [km]

        end if

        !! Particles
        input%Nparticles = 0 ! Número de boulders

        if (input%Nparticles > 0) then
            call allocate_params_particles(input%Nparticles)!! Alocatamos (No tocar)

            particles_in(1,1) = cero   ! Semieje [km]
            particles_in(1,2) = cero   ! Eccentricidad
            particles_in(1,3) = cero   ! M [deg]
            particles_in(1,4) = cero   ! w [deg]
            particles_in(1,5) = 8.1d0  ! MMR
            
        end if

        !!!! Merges
        input%use_merge_part_mass = .True. ! Merge particles into asteroid
        input%use_merge_massive = .True. ! Merge massive bodies

        !!!! Stops
        input%use_stop_no_moon_left = .True. ! Stop if no more moons left
        input%use_stop_no_part_left = .True. ! Stop if no more particles left
        
        !!!! Stokes
        input%use_stokes = .False.
        input%stokes_a_damping_time = infinity                           ! [day]
        input%stokes_e_damping_time = input%stokes_a_damping_time / 1.d2 ! [day]
        input%stokes_active_time = cero                                  ! [day] Tiempo que actúa stokes

        !!!! Naive-Stokes (drag)
        input%use_drag = .False.
        input%drag_coefficient = cero  ! Eta
        input%drag_active_time = cero  ! [day] Tiempo que actúa drag


        !!! Parámetros corrida

        !!!! Tiempos
        input%initial_time = cero      ! Initial time [day]
        input%final_time = 2.d3        ! Final time [day]
        input%case_output_type = 0     ! 0: Linear ; 1: Logarithmic ; 2: Combination
        input%output_timestep = cero   ! Output timestep [day] (used if case_output_type != 1)
        input%output_number = 10000    ! Number of outputs (if output_timestep = 0)

        !!!! Error
        input%learning_rate = 0.85d0   ! [For adaptive step integrators] Learning rate
        input%error_digits = 12        ! [For adaptive step integrators] Digits for relative error

        !!!! Colision y escape
        input%min_distance = -1.5d0    ! Min distance before impact [km] ! 0 => R0 + max(Rboul)
        input%max_distance = -1.d2     ! Max distance before escape [km] ! -x => R0 * x

        !!! Map
        input%use_potential_map = .False.
        input%mapfile = ""
        input%map_grid_size_x = 500
        input%map_grid_size_y = 500
        input%map_min_x = -500.d0
        input%map_max_x = 500.d0
        input%map_min_y = -500.d0
        input%map_max_y = 500.d0

        !!! Output: "" or "no", if not used
        input%datafile = ""
        input%chaosfile = ""
        input%multfile = ""

        !!!!! Screeen
        input%use_screen = .True. ! Print info in screen
        input%use_datascreen = .True. ! Print data in screen

        !!! Input: "" or "no", if not used
        input%tomfile = ""
        input%moonsfile = ""
        input%particlesfile = ""

        !!! Parallel
        input%use_parallel = .False.
        input%requested_threads = 1 ! Number of threads to use !! -1 => all available

    end if


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!! No tocar de aquí a abajo !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    ! Load command line arguments
    call load_command_line_arguments(input, use_configfile)
    
    ! Message
    if ((arguments_number .le. 1) .and. (use_configfile .and. (.not. configfile_exists))) then  ! Global variable
        print*, "WARNING: No command line arguments."
        print*, "¿Do you want to run with in-code set parameters? [y/N]"
        read (*,*) aux_character1
        if (index("YySs", aux_character1) == 0) then
            print*, "Exiting."
            stop
        end if
    end if
    
    ! Init derived parameters
    sim%input_params_st = input ! Create sim with input parameters
    call set_derived_parameters(sim) ! Inicializamos parámetros derivados
    if (sim%use_screen) unit_file = 6  ! Unit_file to std out ( for collisions and escapes )

    
    ! <> Check Output

    ! Mensaje
    if (sim%use_screen .and. .not. configfile_exists) then
        write (*,*) "WARNING: Configurations file could not be read: config.ini"
        write (*,*) "         Using in-code parameters."
    end if

    if (.not. any((/sim%use_datascreen, sim%use_datafile, sim%use_chaosfile, &
                  & sim%use_potential_map, sim%use_multiple_outputs/))) then ! No tiene sentido hacer nada
        if (.not. sim%use_screen) then
            write (*,*) ACHAR(10)
            write (*,*)  "EXITING: No output to be generated."
            stop 1
        else 
            write (*,*)  "WARNING: No integration to be done. Just printing information on screen."
        end if
    end if

    
    ! <> Check Bodies

    !! Read Moons or Particles if needed

    !!! Moons
    if (sim%use_command_body .and. sim%cl_body_is_moon) then
        if (allocated(moons_in)) then  ! Deallocate if necessary
            if (sim%use_screen) then
                write (*,*) ACHAR(10)
                write (*,*) "WARNING: Command line moon used. Ignoring moons in config file."
            end if
            deallocate(moons_in)
        else if (sim%use_moonsfile .and. sim%use_screen) then
            write (*,*) ACHAR(10)
            write (*,*) "WARNING: Command line moon used. Ignoring moons file: ", trim(sim%moonsfile)
        end if
        call allocate_params_moons(1)
        ! Fill the array
        moons_in(1,:) = cl_body_in
    else if (sim%use_moonsfile) then
        if (sim%Nmoons > 0) then
            write (*,*) ACHAR(10)
            write (*,*) "ERROR: Can not read moons from moons file and config file."
            stop 1
        end if
        if (sim%use_screen) write (*,*) "Reading moons from file: ", trim(sim%moonsfile)
        !!!! Formato: [mu, a, e, M, w, MMR, radius]
        call read_columns_file(sim%moonsfile, aux_2D)
        sim%Nmoons = size(aux_2D, 1)
        if (sim%Nmoons .eq. 0) then
            write (*,*) ACHAR(10)
            write (*,*) "ERROR: Could not read moons from moons file."
            stop 1
        end if
        sim%use_moons = .True.
        aux_int = size(aux_2D, 2)
        if (sim%use_screen) then
            write (*,s1i1) "  Read ", sim%Nmoons, " rows."
            write (*,s1i1) "  Read ", aux_int, " columns."
        end if
        call allocate_params_moons(sim%Nmoons)
        ! Fill the arrays
        moons_in = cero
        do i = 1, min(aux_int,7)
            moons_in(:sim%Nmoons,i) = aux_2D(:,i)
        end do
        deallocate(aux_2D)
    end if

    !!! Particles
    if (sim%use_command_body .and. (.not. sim%cl_body_is_moon)) then
        if (allocated(particles_in)) then  ! Deallocate if necessary
            if (sim%use_screen) then
                write (*,*) ACHAR(10)
                write (*,*) "WARNING: Command line particle used. Ignoring particles in config file."
            end if
            deallocate(particles_in)
        else if (sim%use_particlesfile .and. sim%use_screen) then
            write (*,*) ACHAR(10)
            write (*,*) "WARNING: Command line particle used. Ignoring particles file: ", trim(sim%particlesfile)
        end if
        call allocate_params_particles(1)
        ! Fill the array
        particles_in(1,:) = cl_body_in(2:6)
    else if (sim%use_particlesfile) then
        if (sim%Nparticles > 0) then
            write (*,*) ACHAR(10)
            write (*,*) "ERROR: Can not read particles from particles file and config file."
            stop 1
        end if
        if (sim%use_screen) write (*,*) "Reading particles from file: ", trim(sim%particlesfile)
        !!!! Formato: [mu, a, e, M, w, MMR]
        call read_columns_file(sim%particlesfile, aux_2D)
        sim%Nparticles = size(aux_2D, 1)
        if (sim%Nparticles .eq. 0) then
            write (*,*) ACHAR(10)
            write (*,*) "ERROR: Could not read particles from particles file."
            stop 1
        end if
        sim%use_particles = .True.
        aux_int = size(aux_2D, 2)
        if (sim%use_screen) then
            write (*,s1i1) "  Read ", sim%Nparticles, " rows."
            write (*,s1i1) "  Read ", aux_int, " columns."
        end if
        call allocate_params_particles(sim%Nparticles)
        ! Fill the arrays
        particles_in = cero
        do i = 1, min(aux_int,5)
            particles_in(:sim%Nparticles,i) = aux_2D(:,i)
        end do
        deallocate(aux_2D)
    end if


    ! <> Set and check numbers
    !! Set
    sim%Ntotal = 1 + sim%Nmoons + sim%Nparticles
    call update_sim_Nactive(sim, sim%Nparticles, sim%Nmoons)

    !! Check
    if (sim%Ntotal - 1 == 0) then ! No Hay partículas para integrar?
        if (.not. sim%use_potential_map) then
            write (*,*) ACHAR(10)
            write (*,*)  "EXITING: No moon/particle to integrate."
            stop 1
        end if
        only_potential_map = .True.  ! Global
        if (sim%use_screen) then
            write (*,*)  "WARNING: Performing only potential map calculation."
        end if
    end if

    if ((.not. sim%use_boulders) .and. (.not. sim%use_triaxial)) then
        write (*,*) "WARNING: No boulders to integrate. Rotation effects dissabled."
    end if
    

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  Parallel  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (sim%use_screen) then
        write (*,*) ACHAR(5)
        if (sim%use_parallel) then
            write (*,*) "---------- Parallel integration----------"
            write (*,s1i1) "  Available threads:", available_threads
            write (*,s1i1) "  Requested threads:", sim%requested_threads
            write (*,s1i1) "  Used threads     :", my_threads
        else
            write (*,*) " Serial integration (No parallel)"
        end if
        write (*,*) ACHAR(5)
    end if
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!! COMIENZO DE CÁLCULOS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (sim%use_screen) then
        write (*,s1i1) "Starting simulation number: ", simulation_number
        write (*,*) ACHAR(5)
        write (*,*) "---------- Initial parameters ----------"
        write (*,*) ACHAR(5)
    end if


    
    !!!!!!!!!!!!!!!!!!!!!!!!!! Bodies !!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !  ----->   REMEMBER TO SET UNITS BELOW   <---------

    !!!!!!!!!!!!!! Asteroid !!!!!!!!!!!!!!
    
    call allocate_asteroid(asteroid, sim%Nboulders)  ! Alocate
    call add_primary(asteroid, &   ! Create primary 
                    & sim%mass_primary * unit_mass, &  ! mass
                    & sim%triax_a_primary * unit_dist, &  ! Semi-axis a
                    & sim%triax_b_primary * unit_dist, &  ! Semi-axis b
                    & sim%triax_c_primary * unit_dist)    ! Semi-axis c
    do i = 1, sim%Nboulders ! Add boulders
        call add_boulder(asteroid, &
                        & boulders_in(i,1), &              ! mu
                        & boulders_in(i,2) * unit_dist, &  ! radius
                        & boulders_in(i,3) * radian)       ! theta
    end do

    !!!!!!!!!!!!!! Moons !!!!!!!!!!!!!!

    call allocate_moons(moons_arr, sim%Nmoons)
    do i = 1, sim%Nmoons  ! Add boulders
        call add_moon(moons_arr, &
                     & moons_in(i,1), &              ! mu
                     & moons_in(i,2) * unit_dist, &  ! a
                     & moons_in(i,3), &              ! e
                     & moons_in(i,4) * radian, &    ! M
                     & moons_in(i,5) * radian, &    ! w
                     & moons_in(i,6), &              ! MMR
                     & moons_in(i,7) * unit_dist)    ! radius
    end do

    !!!!!!!!!!!!!! Particles !!!!!!!!!!!!!

    ! Init particles
    call allocate_particles(particles_arr, sim%Nparticles)
    do i = 1, sim%Nparticles  ! Add particles
        call add_particle(particles_arr, &
                         & particles_in(i,1) * unit_dist, &  ! a
                         & particles_in(i,2), &              ! e
                         & particles_in(i,3) * radian, &     ! M
                         & particles_in(i,4) * radian, &     ! w
                         & particles_in(i,5), &              ! MMR
                         sim%Nmoons)                         ! ID to star from
    end do
        
    !!!!!!!!!!!!!! SYSTEM !!!!!!!!!!!!!!

    call init_system(system, asteroid, moons_arr, particles_arr, &
                    & sim%lambda_kep, &                             ! keplerian omega
                    & sim%asteroid_rotational_period * unit_time, & ! asteroid period
                    & sim%initial_time * unit_time)

    !! <> Messages
    
    ! Asteroid
    if (sim%use_screen) then
        write (*,*) "Asteroid Rotation:"
        write (*,s1r1) "  a_corot               :", system%asteroid%a_corotation / unit_dist, "[km]"
        write (*,s1r1) "  omega                 :", system%asteroid%omega * unit_time, "[rad day⁻¹]"
        write (*,s1r1) "  omega_kep             :", system%asteroid%omega_kep * unit_time, "[rad day⁻¹]"
        write (*,s1r1) "  lambda_kep            :", system%asteroid%lambda_kep
        write (*,s1r1) "  Period                :", system%asteroid%rotational_period * 24 * unit_time, "[hs]"
        write (*,s1r1) "  Inertial momentum     :", system%asteroid%inertia / (unit_mass * unit_dist**2), "[kg km²]"
        write (*,s1r1) "  Angular momentum (rot):", system%asteroid%ang_mom_rot / unit_mass / &
                                    & (unit_dist**2) * unit_time, "[kg km² day⁻¹]"
        write (*,*) ACHAR(5)
        if ((.not. sim%use_boulders) .and. (.not. sim%use_triaxial)) then
            write(*,*) "Rotation is considered just to define initial parameters."
            write (*,*) ACHAR(5)
        end if
    end if

    ! Individual boulder positions
    allocate(boulders_coords(0:sim%Nboulders, 4))
    if (sim%use_screen .and. sim%use_boulders) then
        write (*,*) "   Asteroid-centric boulders: [x, y, vx, vy, mass, distance, radius]"
        do i = 0, sim%Nboulders
            call get_boulder_i_coord(system%asteroid, i, boulders_coords(i,:))
            write (*,i1r7) i, &
                     & boulders_coords(i,1:2) / unit_dist, &
                     & boulders_coords(i,3:4) / unit_vel, &
                     & system%asteroid%boulders(i)%mass / unit_mass, &
                     & system%asteroid%boulders(i)%dist_to_asteroid / unit_dist, &
                     & system%asteroid%boulders(i)%radius / unit_dist
        end do
        write (*,*) ACHAR(5)
    end if

    ! Moons positions
    if (sim%use_screen .and. sim%use_moons) then
        write (*,*) "   Barycentric moons: [x, y, vx, vy, mass, distance]"
        do i = 1, sim%Nmoons
            write (*,i1r7) i, &
                     & system%moons(i)%coordinates(1:2) / unit_dist, &
                     & system%moons(i)%coordinates(3:4) / unit_vel, &
                     & system%moons(i)%mass / unit_mass, &
                     & system%moons(i)%dist_to_cm / unit_dist
        end do
        write (*,*) ACHAR(5)
    end if

    ! Particles positions
    if (sim%use_screen .and. sim%use_particles) then
        write (*,*) "   Barycentric particles: [x, y, vx, vy, distance]"
        do i = 1, sim%Nparticles
            write (*,i1r7) i, &
                     & system%particles(i)%coordinates(1:2) / unit_dist, &
                     & system%particles(i)%coordinates(3:4) / unit_vel, &
                     & system%particles(i)%dist_to_cm / unit_dist
        end do
        write (*,*) ACHAR(5)
    end if

    ! Energy and ang_mom
    if (sim%use_screen) then
        write (*,s1r1) "   Mass            : ", system%mass / unit_mass, "[kg]"
        write (*,s1r1) "   Energy          : ", system%energy / unit_ener * (metro**2 / segundo), "[J]"
        write (*,s1r1) "   Angular Momentum: ", system%ang_mom / unit_angm * (metro**2 / segundo**2), "[J s⁻¹]"
        write (*,*) ACHAR(5)
    end if

    ! <<<< Save initial data >>>>
    initial_system = system
    allocate(boulders_data(0:sim%Nboulders, 4))  ! To use in dydt
    do i = 0, sim%Nboulders
        call get_boulder_i_coord(system%asteroid, i, boulders_coords(i,:))
        boulders_data(i,1) = system%asteroid%boulders(i)%mass
        boulders_data(i,2) = system%asteroid%boulders(i)%radius
        boulders_data(i,3) = system%asteroid%boulders(i)%initial_theta
        boulders_data(i,4) = system%asteroid%boulders(i)%dist_to_asteroid
    end do


    ! <> Parallel revisted
    ! Redefine OMP threads if necessary. Lower equal to particle number
    if ((sim%use_parallel) .and. ((sim%Ntotal - 1) < my_threads)) then
        my_threads = sim%Ntotal - 1
        if (sim%use_screen) then
            write (*,s1i1) "WARNING: Threads reduced to ", my_threads
            write (*,*) ACHAR(5)
        end if
        !$ call omp_set_num_threads(my_threads)
    end if


    !!!!!!!!!!!!!!!!!!!!!!!!!!!! EXTERNA EFFECTS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (sim%use_screen) then
        write (*,*) ACHAR(5)
        write (*,*) ("---------- External/Internal forces ----------")
        write (*,*) ACHAR(5)
    end if

    !! Tri-axial gravity
    if (sim%use_triaxial) then
        call init_ellipsoid(sim%triax_a_primary * unit_dist, sim%triax_b_primary * unit_dist, sim%triax_c_primary * unit_dist)
        if (sim%use_screen) then
            write (*,*) "Tri-axial potential"
            write (*,s1r1) "  semi-axis-a :", sim%triax_a_primary, "[km]"
            write (*,s1r1) "  semi-axis-b :", sim%triax_b_primary, "[km]"
            write (*,s1r1) "  semi-axis-c :", sim%triax_c_primary, "[km]"
            write (*,s1r1) "   Effective Radius :", sim%Reffective, "[km]"
            write (*,s1r1) "   C_20:", sim%C20
            write (*,s1r1) "   C_22:", sim%C22
            write (*,*) ACHAR(5)
        end if
    end if


    !! Moons gravity
    if (sim%use_moon_gravity .and. sim%use_screen) then
        write (*,*) "Gravity between moons DEACTIVATED."
        write (*,*) ACHAR(5)
    end if

    !! Omega Damping
    if (sim%use_omega_damping .and. ((.not. sim%use_boulders) .and. (.not. sim%use_triaxial))) then
        if (sim%use_screen) write (*,*) "WARNING: Skipping omega damping without boulders."
        sim%use_omega_damping = .False.
    else if (sim%use_omega_damping) then
        if (sim%use_lin_omega_damp) then
            aux_real = - system%asteroid%omega / (sim%omega_lin_damping_time * unit_time - sim%initial_time)
            call init_damping(aux_real, &
                            & cero, &
                            & sim%omega_damp_active_time * unit_time, &
                            & 1)
            if (sim%use_screen) then
                write (*,*) "Omega Damping: linear"
                write (*,s1r1) " tau_o :", sim%omega_lin_damping_time / (system%asteroid%rotational_period / unit_time), "[Prot]"
                write (*,*) ACHAR(5)
            end if
        else if (sim%use_exp_omega_damp) then
            call init_damping(sim%omega_exp_damping_time * unit_time, &
                            & cero, &
                            & sim%omega_damp_active_time * unit_time, &
                            & 2)
            if (sim%use_screen) then
                write (*,*) "Omega Damping: exponential"
                write (*,s1r1) " tau_o :", sim%omega_exp_damping_time / (system%asteroid%rotational_period / unit_time), "[Prot]"
                write (*,*) ACHAR(5)
            end if
        else if (sim%use_poly_omega_damp) then
            call init_damping(sim%omega_exp_damp_poly_A, &
                            & sim%omega_exp_damp_poly_A, &
                            & sim%omega_damp_active_time * unit_time, &
                            & 3)
            
            if (sim%use_screen) then
                write (*,*) "Omega Damping: poly-exponential"
                write (*,s1r1) "  A :", sim%omega_exp_damp_poly_A
                write (*,s1r1) "  B :", sim%omega_exp_damp_poly_B
                write (*,*) ACHAR(5)
            end if
        end if
    end if

    !! Stokes
    if (sim%use_stokes) then
        call init_stokes(sim%stokes_a_damping_time * unit_time, &
                       & sim%stokes_e_damping_time * unit_time, &
                       & sim%stokes_active_time * unit_time)
        if (sim%use_screen) then
            write (*,*) "Stokes"
            write (*,s1r1) " tau_a   : ", sim%stokes_a_damping_time, "[days]"
            write (*,s1r1) " tau_e   : ", sim%stokes_e_damping_time, "[days]"
            write (*,s1r1) " t_stokes: ", sim%stokes_active_time, "[days]"
            write (*,*) ACHAR(5)
        end if
    end if

    !! Naive-Stokes (Drag)
    if (sim%use_drag) then
        call init_drag(sim%drag_coefficient, sim%drag_active_time * unit_time)
        if (sim%use_screen) then
            write (*,*) "Drag"
            write (*,s1r1) " eta   : ", sim%drag_coefficient
            write (*,s1r1) " t_drag: ", sim%drag_active_time, "[days]"
            write (*,*) ACHAR(5)
        end if
    end if

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Escape/Colisión
    if (sim%min_distance < cero) then
        sim%min_distance = system%asteroid%radius * abs(sim%min_distance)
    else if (sim%min_distance < tini) then
        sim%min_distance = system%asteroid%radius
    else
        sim%min_distance = sim%min_distance * unit_dist
    end if
    if (sim%max_distance < cero) then
        sim%max_distance = system%asteroid%radius * abs(sim%max_distance)
    else 
        sim%max_distance = sim%max_distance * unit_dist
    end if
    if ((sim%max_distance > cero) .and. (sim%max_distance <= sim%min_distance)) then
        write (*,*) ACHAR(10)
        write (*,*) "ERROR: rmax <= rmin"
        stop 1
    end if
    if (sim%use_screen) then
        write (*,*) "Conditions for escape/collision"
        write (*,s1r1) "  rmin : ", sim%min_distance / unit_dist, "[km] =", sim%min_distance / system%asteroid%radius, "[Rast]"
        if (sim%max_distance > cero) then
            write (*,s1r1) "  rmax : ", sim%max_distance / unit_dist, "[km] =", sim%max_distance / system%asteroid%radius, "[Rast]"
        else
            write (*,*) "  rmax : Infinity"
        end if
        if (sim%use_any_merge) then
            if (sim%use_merge_part_mass) write (*,*) " Colliding particles into massive bodies will be removed."
            if (sim%use_merge_massive) write (*,*) " Colliding massive bodies will be merged."
        else
            write (*,*) " Collisions will not be solved."
        end if
        if (sim%use_any_stop) then
            if (sim%use_stop_no_part_left) write (*,*) " Simulation will stop if no more particles are left."
            if (sim%use_stop_no_moon_left) write (*,*) " Simulation will stop if no more moons are left."
        else
            write (*,*) " Simulation will stop if no more particles and moons are left."
        end if
        write (*,*) ACHAR(5)
    end if


    !!!!!!!!!!!!!!!!!!!!! MAPA Potencial y Aceleraciones !!!!!!!!!!!!!!!!!!!!!!
    if (sim%use_potential_map) then
        if (sim%use_screen) then
            write (*,*) ACHAR(5)
            write (*,*) "---------- MAP (Gravitational Energy and Acceleration) ----------"
            write (*,*) ACHAR(5)
            write (*,*) "Creating maps..."
        end if
        call create_map(sim, system)
        if (sim%use_screen) then
            write (*,*) "Maps saved to file: ", trim(sim%mapfile)
            write (*,*) ACHAR(5)
        end if
    end if
    if (only_potential_map) then  ! Global
        if (sim%use_screen) then
            write (*,*) "END: Only the Map was requested."
            write (*,*) ACHAR(5)
        end if
        stop 0
    end if

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!! FIN DE CÁLCULOS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!! INTEGRACIÓN !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Mensaje
    if (sim%use_screen) then
        write (*,*) ACHAR(5)
        write (*,*) "---------- Preparing integration ----------"
        write (*,*) ACHAR(5)
    end if


    !!!!!!!! TIEMPOS !!!!!!!

    ! Tiempos de integración
    sim%initial_time = sim%initial_time * unit_time
    if (sim%final_time < cero) then
        if (abs(system%asteroid%rotational_period) < tini) then
            write (*,*) "ERROR: Can not calculate total integration time with non-rotating asteroid."
        end if
        sim%final_time = abs(sim%final_time) * system%asteroid%rotational_period
    else
        sim%final_time = sim%final_time * unit_time
    end if
    if (sim%initial_time >= sim%final_time) then
        write (*,*) ACHAR(10)
        write (*,*) "ERROR: t0 >= tf"
        stop 1
    end if

    !! Timesteps
    sim%output_timestep = sim%output_timestep * unit_time
    min_timestep = min_timestep * unit_time ! Global

    !! Output times
    call set_output_times(sim%initial_time, sim%final_time, sim%output_number, &
                        & sim%output_timestep, sim%case_output_type, output_times)

    !! TOMFILE
    if (sim%use_tomfile) then
        !! En este caso, leeremos los tiempos desde un archivo
        if (sim%use_screen) write (*,*) "Reading times from TOM file: ", trim(sim%tomfile)
        call read_tomfile(sim%initial_time / unit_time, sim%final_time / unit_time, & 
                        & tom_times, tom_deltaomega, tom_deltamass, sim%tomfile) ! Read LOOP checkpoints
        tom_total_number = size(tom_times, 1)
        !!! Unidades
        tom_times = tom_times * unit_time
        if (allocated(tom_deltaomega)) tom_deltaomega = tom_deltaomega / unit_time
        if (allocated(tom_deltamass)) tom_deltamass = tom_deltamass * unit_mass
        !!! Condicion inicial (y final)
        tom_times(1) = sim%initial_time
        tom_times(tom_total_number) = sim%final_time
        if (allocated(tom_deltaomega)) then
            tom_deltaomega(1) = cero
            tom_deltaomega(tom_total_number) = cero
        end if
        if (allocated(tom_deltamass)) then
            tom_deltamass(1) = cero
            tom_deltamass(tom_total_number) = cero
        end if
        !!! Mensaje !
        if (sim%use_screen) then
            if (allocated(tom_deltamass)) then
                write (*,*) "  - 3 columns read: t, Delta_omega, Delta_m"
            else if (allocated(tom_deltaomega)) then
                write (*,*) "  - 2 columns read: t, Delta_omega"
            else
                write (*,*) "  - 1 column read: t"
            end if
        end if
        !! Ahora debemos combinar los tiempos de TOM con los tiempos de Output, y crear un nuevo vector de tiempos Checkpoints
        call merge_sort_and_unique(tom_times, output_times, &
                                   & checkpoint_is_tom, checkpoint_is_output, &
                                   & checkpoint_times, checkpoint_number)
        if (allocated(tom_deltaomega) .and. ((.not. sim%use_boulders) .and. (.not. sim%use_triaxial))) then
            if (sim%use_screen) write (*,*) "WARNING: Delta_omega has no sense without boulders. It will be ignored."
            deallocate(tom_deltaomega)            
        end if
    else
        !! En este caso, los tiempos de check son los mismos que los de LOOP
        checkpoint_number = sim%output_number
        allocate(checkpoint_times(checkpoint_number))
        allocate(checkpoint_is_output(checkpoint_number))
        allocate(checkpoint_is_tom(checkpoint_number))
        checkpoint_times = output_times
        checkpoint_is_output = .True.
        checkpoint_is_tom = .False.
        tom_total_number = 0
    end if

    ! Variable temporal (para integración)
    time = sim%initial_time ! Tiempo actual
    timestep = output_times(1) - sim%initial_time ! Paso de tiempo inicial
    min_timestep = max(min(min_timestep, sim%output_timestep), tini) ! Paso de tiempo mínimo
    if (system%asteroid%rotational_period > tini) then
        adaptive_timestep = system%asteroid%rotational_period * 0.01d0 ! Paso de tiempo adaptativo inicial: 1% del periodo de rotación
    else
        adaptive_timestep = infinity
        do i = 1, system%Nmoons_active
            adaptive_timestep = min(adaptive_timestep, &
                                  & get_Period(system%asteroid%mass + system%moons(i)%mass, system%moons(i)%elements(1)))
        end do
        do i = 1, system%Nparticles_active
            adaptive_timestep = min(adaptive_timestep, &
                                  & get_Period(system%asteroid%mass, system%particles(i)%elements(1)))
        end do
        adaptive_timestep = adaptive_timestep * unit_time
    end if

    !! Mensaje
    if (sim%use_screen) then
        write (*,*) "Times:"
        if (system%asteroid%rotational_period > tini) then
            write (*,s1r1) "    t0    : ", sim%initial_time / system%asteroid%rotational_period, "[Prot] = ", &
            & sim%initial_time / unit_time, "[day]"
            write (*,s1r1) "    tf    : ", sim%final_time / system%asteroid%rotational_period, "[Prot] = ", &
            & sim%final_time / unit_time, "[day]"
            write (*,s1r1) "    dt_out: ", sim%output_timestep / system%asteroid%rotational_period, "[Prot] = ", &
            & sim%output_timestep / unit_time, "[day]"
            write (*,s1r1) "    dt_min: ", min_timestep / system%asteroid%rotational_period, "[Prot] = ", &
            & min_timestep / unit_time, "[day]"
        else 
            write (*,s1r1) "    t0    : ", sim%initial_time / unit_time, "[day]"
            write (*,s1r1) "    tf    : ", sim%final_time / unit_time, "[day]"
            write (*,s1r1) "    dt_out: ", sim%output_timestep / unit_time, "[day]"
            write (*,s1r1) "    dt_min: ", min_timestep / unit_time, "[day]"
        end if
        write (*,s1i1) "    n_out : ", sim%output_number
        write (*,*) ACHAR(5)
    end if



    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!! Integration Arrays !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Allocate
    allocate(m_arr(sim%Ntotal))  ! Masses
    allocate(R_arr(sim%Ntotal))  ! Radii
    allocate(y_arr(2 + sim%Ntotal * 4))     ! Theta, Omega, Positions
    allocate(y_arr_new(2 + sim%Ntotal * 4)) ! Theta, Omega, Positions
    allocate(y_der(2 + sim%Ntotal * 4))     ! derivate(Theta, Omega, Positions)

    ! Init to 0
    m_arr = cero
    R_arr = cero
    y_arr = cero
    y_arr_new = cero
    y_der = cero

    ! get amount of values of Y to use
    y_nvalues = get_index(sim%Nactive) + 3  ! Update nvalues to use in y

    ! Get Arrays to integrate
    call generate_arrays(system, m_arr, R_arr, y_arr)

    ! CHECKS
    if (sim%use_screen) then
        write (*,*) ACHAR(5)
        write (*,*) "---------- FINAL CHECKS ----------"
        write (*,*) ACHAR(5) 
    end if

    !!! Check workers
    if (sim%use_parallel) then
        aux_int = MIN(MAX(sim%Ntotal, 3), my_threads)
        if (aux_int .ne. my_threads) then
            my_threads = aux_int
            !$ call OMP_SET_NUM_THREADS(my_threads)
            if (sim%use_screen) write (*,*) "WARNING: Se usarán ", my_threads, " hilos para la integración."
        end if
    end if


    !!!! Mensaje
    if (sim%use_screen .and. (sim%error_tolerance <= 1.d-16)) write (*,*) " WARNING: e_tol might be too low (<= 10⁻¹⁶)"


    !!!! Ahora si, si no hay que hacer más nada entonces terminamos.
    if (.not. any((/sim%use_datascreen, sim%use_datafile, sim%use_chaosfile, & 
                  & sim%use_potential_map, sim%use_multiple_outputs/))) then ! No tiene sentido hacer nada
        write (*,*) ACHAR(10)
        write (*,*) "As nothing will be stored, there is no integration."
        write (*,*) "Exiting."
        stop 1
    end if

    !!! Check if not too much multiple files
    if (sim%use_multiple_outputs) then
        if (sim%Ntotal > 1000) then
            write (*,*) ACHAR(5)
            write (*,*) "WARNING: More than 1000 output files will be created."
            write (*,*) "         This could generate an issue."
            write (*,*) "Do you wish to continue? (y/[N])"
            read (*,*) aux_character30
            aux_character1 = trim(adjustl(aux_character30))
            if ((aux_character1 == "y") .or. (aux_character1 == "s")) then
                if (sim%use_screen) then
                    write (*,*) "Resuming..."
                    write (*,*) ACHAR(5)
                end if
            else
                write (*,*) "Exiting."
                stop 1
            end if
        end if
    end if

    !!! Mensaje
    if (sim%use_screen) then
        if (sim%use_datafile) write (*,*) "General output file: ", trim(sim%datafile)
        if (sim%use_multiple_outputs) write (*,*) "Individuals output files: ", trim(sim%multfile) // "_*"
        if (sim%use_datascreen) then
            write (*,*) "Data will be printed on screen."
        else 
            write (*,*) "Data will not be printed on screen."
        end if
        write (*,*) ACHAR(5)
        write (*,*) "FINAL CKECHS OK"
        write (*,*) ACHAR(5)
    end if

    ! Free initial arrays
    call free_initial_arays()

    !!!! Mensaje de inicio
    if (sim%use_screen) then
        write (*,*) ACHAR(5)
        write (*,*) "---------- INTEGRATING ----------"
        write (*,*) ACHAR(5) 
    end if

    ! Inicializamos Punteros
    call define_pointers(sim)

    
    ! ABRIMOS ARCHIVOS
    !! Archivo de salida general
    if (sim%use_datafile) open (unit=20, file=trim(sim%datafile), status='replace', action='write', position="append")
    !! Archivos individuales
    if (sim%use_multiple_outputs) then
        do i = 0, sim%Ntotal  ! 0 is the asteroid
            write (aux_character20, *) i
            open (unit=200+i, file=trim(sim%multfile) // "_" // trim(adjustl(aux_character20)), &
                & status='replace', action='write', position="append")
        end do
    end if
    !! Chaos File
    if (sim%use_chaosfile) open (unit=40, file=trim(sim%chaosfile), status='replace', action='readwrite', position="append")

    !!!!!! MAIN LOOP INTEGRATION !!!!!!!
    keep_integrating = .True.  ! Init flag

    ! CHECK INITIAL CONDITIONS
    ! Apply colissions/escapes and checks
    call resolve_collisions(system, sim%min_distance, unit_file)
    call check_after_col(system, sim, keep_integrating)
    call resolve_escapes(system, sim%max_distance, unit_file)
    call check_after_esc(system, sim, keep_integrating)

    ! Update Nactive and y_arr if necessary
    call get_Nactive(system, new_Nactive)
    if (new_Nactive < sim%Nactive) then
        call update_sim_Nactive(sim, system%Nparticles_active, system%Nmoons_active)  ! Update sim Nactive
        y_nvalues = get_index(new_Nactive) + 3  ! Update nvalues to use in y
        call generate_arrays(system, m_arr, R_arr, y_arr)  ! Regenerate arrays
    end if
    
    ! Update Chaos
    call update_chaos(system, sim%use_baryc_output)

    !! Reduce threads if necessary...
    if (sim%use_parallel) then
        aux_int = MIN(MAX(sim%Nactive, 3), my_threads)
        if (aux_int .ne. my_threads) then
            my_threads = aux_int
            !$ call OMP_SET_NUM_THREADS(my_threads)
            if (sim%use_screen) then
                write (*,s1i1) " WARNING: Paralelism reduced to ", my_threads, " threads for the integration."
                write (*,*) ACHAR(5)
            end if
        end if
    end if

    ! Output initial conditions
    call write_a_to_individual(system, 200)
    !$OMP PARALLEL DEFAULT(SHARED) &
    !$OMP PRIVATE(i)
    !$OMP DO
    do i = 1, system%Nmoons_active
        call write_m_to_individual(system, i, 200+system%moons(i)%id)
    end do 
    !$OMP END DO NOWAIT
    !$OMP DO
    do i = 1, system%Nparticles_active
        call write_p_to_individual(system, i, 200+system%particles(i)%id)
    end do 
    !$OMP END DO NOWAIT
    !$OMP SECTIONS
    !$OMP SECTION
    call write_to_screen(system, 6)
    !$OMP SECTION
    call write_to_general(system, 20)
    call flush_output(20)
    !$OMP SECTION
    call flush_chaos(40) ! Update chaos
    !$OMP END SECTIONS
    !$OMP END PARALLEL

    ! >>>>>>>>>>>>>>>>>>----------------- MAIN LOOP  ---------------<<<<<<<<<<<<<<<<<<<<<<<<
    tom_index_number = 2 !!!! Inicializamos en 2 porque el primer checkpoint es el IC (t0)
    j = 2! From 2 because 1 is the IC (t0) !! +1 por si hay HardExit en el último
    main_loop: do while (.True.)
    
        hard_exit = .False. ! Resetear HardExit 
        is_premature_exit = .False. ! Resetear premature
        
        ! Check if all done
        !! Time end
        if (j == checkpoint_number + 1) keep_integrating = .False.

        !! Check if Particles/Moons left
        if (sim%Nactive == 1 .and. keep_integrating) then
            if (sim%use_screen) then
                write (*,*) ACHAR(5) 
                write (*,*) " No more active particles/moons left."
            end if
            keep_integrating = .False.
        end if

        ! Keep going?
        if (.not. keep_integrating) then
            if (sim%use_screen) then
                write (*,*) ACHAR(5) 
                write (*,s1r1) "Integration finished at time = ", time / unit_time, "[days]"
            end if
            exit main_loop
        end if

        
        ! Update dt
        timestep = checkpoint_times(j) - time


        !! Check TOM
        if (checkpoint_is_tom(j) .and. (tom_index_number <= tom_total_number)) then
            if (allocated(tom_deltaomega)) call spin_asteroid(system%asteroid, system%asteroid%theta, &
                                                            & system%asteroid%omega + tom_deltaomega(tom_index_number))
            if (allocated(tom_deltamass)) call grow_asteroid(system%asteroid, tom_deltamass(tom_index_number))
            tom_index_number = tom_index_number + 1
            if (sim%use_screen .and. (allocated(tom_deltaomega) .or. allocated(tom_deltamass))) then
                write (*,*) ACHAR(5)
                write (*,s1r1) "Updated asteroid Omega and mass following TOM data, at time = ", time / unit_time, "[days]"
                write (*,*) ACHAR(5)
            end if

            ! Apply colissions/escapes and checks
            call resolve_collisions(system, sim%min_distance, unit_file)
            call check_after_col(system, sim, keep_integrating)
            call resolve_escapes(system, sim%max_distance, unit_file)
            call check_after_esc(system, sim, keep_integrating)

            ! Update Nactive and y_arr if necessary
            call get_Nactive(system, new_Nactive)
            if (new_Nactive < sim%Nactive) then
                call update_sim_Nactive(sim, system%Nparticles_active, system%Nmoons_active)  ! Update sim Nactive
                y_nvalues = get_index(new_Nactive) + 3  ! Update nvalues to use in y
            end if

            call generate_arrays(system, m_arr, R_arr, y_arr)  ! Mandatory bc of new asteroid
        end if

        !!! Execute an integration method (uncomment/edit one of these)
        ! call integ_caller (time, y_arr(:y_nvalues), adaptive_timestep, dydt, &
        !     & Ralston4, timestep, y_arr_new(:y_nvalues), check_func)
        ! call rk_half_step_caller (time, y_arr(:y_nvalues), adaptive_timestep, dydt, &
        !     & Runge_Kutta5, 5, sim%error_tolerance, learning_rate, min_timestep, timestep, y_arr_new(:y_nvalues), check_func)
        ! call embedded_caller (time, y_arr(:y_nvalues), adaptive_timestep, dydt, Dormand_Prince8_7, &
        !    & sim%error_tolerance, learning_rate, min_timestep, timestep, y_arr_new(:y_nvalues), check_func)
        call BStoer_caller (time, y_arr(:y_nvalues), adaptive_timestep, dydt, &
            & sim%error_tolerance, timestep, y_arr_new(:y_nvalues), check_func)
        ! call BStoer_caller2 (time, y_arr(:y_nvalues), adaptive_timestep, dydt, &
        !     & sim%error_tolerance, timestep, y_arr_new(:y_nvalues), check_func)
        ! call leapfrog_caller (time, y_arr(:y_nvalues), adaptive_timestep, dydt, &
        !     & leapfrof_KDK, sim%error_tolerance, y_nvalues, min_timestep, timestep, y_arr_new(:y_nvalues), check_func)
    
        ! Check if it might be hard_exit
        if (hard_exit) then
            !! If so, the dt used is in dt_adap
            timestep = adaptive_timestep

            !! Check if premature_exit: If exit at less than 1% of finishing timestep
            if (timestep > cero) then
                is_premature_exit = (checkpoint_times(j) - (time + timestep)) / timestep > 0.01d0
            else
                is_premature_exit = checkpoint_times(j) - time < tini
            end if
        end if

        ! Update parameters
        time = time + timestep
        y_arr = y_arr_new

        y_arr(1) = modulo(y_arr(1), twopi)  ! Modulate theta

        ! Update from y_new
        call update_system_from_array(system, time, y_arr)

        ! Apply colissions and check
        call resolve_collisions(system, sim%min_distance, unit_file)
        call check_after_col(system, sim, keep_integrating)

        ! Apply escapes and check
        call resolve_escapes(system, sim%max_distance, unit_file)
        call check_after_esc(system, sim, keep_integrating)

        ! Update Nactive and y_arr if necessary
        call get_Nactive(system, new_Nactive)
        if (new_Nactive < sim%Nactive) then
            call update_sim_Nactive(sim, system%Nparticles_active, system%Nmoons_active)  ! Update sim Nactive
            y_nvalues = get_index(new_Nactive) + 3  ! Update nvalues to use in y
            call generate_arrays(system, m_arr, R_arr, y_arr)  ! Regenerate arrays
        end if
        
        ! Update Chaos
        call update_chaos(system, sim%use_baryc_output)
        
        ! Output
        if ((checkpoint_is_output(j)) .and. (.not. is_premature_exit)) then
            call write_a_to_individual(system, 200)
            !$OMP PARALLEL DEFAULT(SHARED) &
            !$OMP PRIVATE(i)
            !$OMP DO
            do i = 1, system%Nmoons_active
                call write_m_to_individual(system, i, 200+system%moons(i)%id)
            end do 
            !$OMP END DO NOWAIT
            !$OMP DO
            do i = 1, system%Nparticles_active
                call write_p_to_individual(system, i, 200+system%particles(i)%id)
            end do 
            !$OMP END DO NOWAIT
            !$OMP SECTIONS
            !$OMP SECTION
            call write_to_screen(system, 6)
            if (sim%use_diagnostics) call write_diagnotics(initial_system, system, 6)
            !$OMP SECTION
            call write_to_general(system, 20)
            call flush_output(20)
            !$OMP SECTION
            call write_to_chaos(system, initial_system, 40)
            call flush_chaos(40) ! Update chaos
            !$OMP END SECTIONS
            !$OMP END PARALLEL
        end if

        if (sim%use_percentage) call percentage(time, sim%final_time)

        ! Update j; only if not HardExit
        if (.not. is_premature_exit) j = j + 1
        
    end do main_loop

    !! Porcentaje final
    if (sim%use_percentage .and. time .ne. sim%final_time) call percentage(sim%final_time + uno, sim%final_time)

    !! Cerrar archivo de salida
    if (sim%use_datafile) then
        close (2)
        if (sim%use_screen) then
            write (*,*) ACHAR(10)
            write (*,*) "Output data saved to file: ", trim(sim%datafile)
        end if
    end if

    !! Cerrar archivos individuales
    if (sim%use_multiple_outputs) then
        do i = 0, sim%Ntotal
            close (200+i)
        end do
        if (sim%use_screen) then
            write (*,*) ACHAR(10)
            write (*,*) "Individual output data saved to files: ", trim(sim%multfile) // "_*"
        end if
    end if

    !! Final chaos
    if (sim%use_chaos) then
        if (sim%use_chaosfile) then
            !! Chaosfile
            inquire (unit=40, opened=aux_logical)
            if (.not. aux_logical) then
                open (unit=40, file=trim(sim%chaosfile), status='unknown', action='write')
            else
                rewind(40)
            end if
            call write_chaos(initial_system, system, 40)
            close (40)
            !! Mensaje
            if (sim%use_screen) then 
                write (*,*) ACHAR(10)
                write (*,*) "Chaos data saved into: ", trim(sim%chaosfile)
            end if
        end if
        if (sim%use_datascreen) then
            write (*,*) ACHAR(10)        
            write (*,*) "Chaos:"
            call write_chaos(initial_system, system, 6)
        end if
    end if

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!! LIBERACION MEMORIA!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (sim%use_screen) then
        write (*,*) ACHAR(10)
        write (*,*) "Freeing memory..."
        write (*,*) ACHAR(5)
    end if

    !!!!!!!!! TIMES !!!!!!!!
    call free_times_arrays()

    !!!!!!!!! PARAMS !!!!!!!!
    call free_parameters_arays()
    
    !!!!!!!!! POINTERS !!!!!!!!
    call nullify_pointers()    

    !!!!!!!!! BODIES !!!!!!!!
    call free_asteroid(asteroid)
    call free_system(system)
    call free_system(initial_system)

    if (sim%use_screen) then 
        write (*,*) ACHAR(5)
        write (*,*) "End of program"
        write (*,*) ACHAR(5)
    end if
    
end program main

subroutine create_map(sim, system)
    use constants, only: cero
    use bodies, only: system_st, get_acc_and_pot_xy
    use parameters, only: sim_params_st
    implicit none
    type(sim_params_st), intent(in) :: sim
    type(system_st), intent(in) :: system
    real(kind=8), dimension(:,:), allocatable :: pot
    real(kind=8), dimension(:,:,:), allocatable :: acc
    real(kind=8) :: dx_nx, dy_ny, rb(2)
    integer(kind=4) :: i, j

    !!! Check
    if ((sim%map_grid_size_x < 2) .or. (sim%map_grid_size_y < 2)) then
        write (*,*) "ERROR: Check that (ngx > 2) and (ngy > 2)"
        stop 1
    end if
    if ((sim%map_min_x >= sim%map_max_x) .or. (sim%map_min_y >= sim%map_max_y)) then
        write (*,*) "ERROR: Check that (xmin < xmax) and (ymin < ymax)"
        stop 1
    end if

    ! Allocate
    allocate(pot(sim%map_grid_size_x, sim%map_grid_size_y))
    allocate(acc(sim%map_grid_size_x, sim%map_grid_size_y, 2))
    
    pot = cero
    acc = cero
    dx_nx = (sim%map_max_x - sim%map_min_x) / sim%map_grid_size_x
    dy_ny = (sim%map_max_y - sim%map_min_y) / sim%map_grid_size_y

    !$OMP PARALLEL DEFAULT(SHARED) &
    !$OMP PRIVATE(i,j,rb)
    !$OMP DO SCHEDULE (STATIC)
    do i = 1, sim%map_grid_size_x
        do j = 1, sim%map_grid_size_y
            rb(1) = sim%map_min_x + i * dx_nx
            rb(2) = sim%map_min_y + j * dy_ny
            call get_acc_and_pot_xy(system, rb, acc(i,j,:), pot(i,j))
        end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL
    
    if (sim%use_screen) write (*,*) "Potential calculated. Writing file..."
    open (unit=50, file=trim(sim%mapfile), status='replace', action='write')
    do i = 1, sim%map_grid_size_x
        do j = 1, sim%map_grid_size_y
            write (50,*) sim%map_min_x + i * dx_nx, sim%map_min_y + j * dy_ny, pot(i,j), acc(i,j,:)
        end do
    end do
    close (50)

    ! Deallocate
    deallocate(pot)
    deallocate(acc)
end subroutine create_map

! ! module iso_fortran_env

! !     ! Nonintrinsic version for Lahey/Fujitsu Fortran for Linux. 
! !     ! See Subclause 13.8.2 of the Fortran 2003 standard. 
  
! !     implicit none 
! !     public 
  
! !     integer, parameter :: Character_Storage_Size = 8 
! !     integer, parameter :: Error_Unit = 0 
! !     integer, parameter :: File_Storage_Size = 8 
! !     integer, parameter :: Input_Unit = 5 
! !     integer, parameter :: IOSTAT_END = -1 
! !     integer, parameter :: IOSTAT_EOR = -2 
! !     integer, parameter :: Numeric_Storage_Size = 32 
! !     integer, parameter :: Output_Unit = 6 
  
! !   end module iso_fortran_env

