program main
    use parameters
    use integrators
    use celestial, only: get_Period, coord2geom
    use bodies
    use times
    use derivates
    use filtering

    implicit none

    external :: create_map   ! Declare the external function for map
    logical :: only_potential_map = .False.  ! Change this if only map desired

    ! Auxiliar variables
    integer(kind=4) :: aux_int
    character(1) :: aux_character1
    character(20) :: aux_character20
    logical :: aux_logical
    real(kind=8) :: aux_real
    real(kind=8), dimension(:), allocatable :: aux_1D
    real(kind=8), dimension(:,:), allocatable :: aux_2D  ! To read files

    integer(kind=4) :: i, j  ! Loops
    integer(kind=4) :: new_Nactive = 0  ! When escapes/collisions occur
    integer(kind=4) :: unit_file = -1  ! Where to write escapes/collisions

    use_configfile = .True. ! Manual override if needed
    if (use_configfile) call read_config_file(input, "config.ini", configfile_exists) ! Leemos archivo de configuración

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!! Edit from here to change configuration without configfile !!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (.not. use_configfile) then ! Usaremos parámetros por defecto

        ! Just print? (True) Or integrate? (False)
        input%only_print = .False.
        
        ! Bodies

        !! Central Asteroid
        !!! Primary mass
        input%mass_primary = 6.3d18 ! Masa del cuerpo 0 [kg] ! -x =>  mAst = x

        !!! Primary shape

        !!!! Tri_Axial
        input%use_triaxial = .False.  ! Logical para determinar si se usa triax, o Boulders
        input%triax_a_primary = cero  ! Semieje a
        input%triax_b_primary = cero  ! Semieje b
        input%triax_c_primary = cero  ! Semieje c

        !!!! Explicit primary radius
        input%radius_primary = 129.d0 ! Radio del cuerpo 0 [km]

        !!! Rotation
        input%asteroid_rotational_period = 7.004d0/24.d0  ! Periodo de rotación [day]
        !lambda_kep = 0.471d0      ! Cociente omega/wk

        !! Boulders        
        input%Nboulders = 1 ! Número de boulders

        if (input%Nboulders > 0) then  ! Specify boulders properties
            call allocate_params_asteroid(input%Nboulders)!! Alocatamos (No tocar)

            boulders_in(1,1) = 1.d-1   ! Cociente de masas entre boulder 1 y primary
            boulders_in(1,2) = cero    ! Ángulo de fase del boulder 1 [deg]
            boulders_in(1,3) = 2.5d0   ! Radio del boulder 1 [km]       

        end if

        !! Moons
        input%Nmoons = 0 ! Número de boulders

        if (input%Nmoons > 0) then  ! Specify moons properties
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

        if (input%Nparticles > 0) then ! Specify particles properties
            call allocate_params_particles(input%Nparticles)!! Alocatamos (No tocar)

            particles_in(1,1) = cero   ! Semieje [km]
            particles_in(1,2) = cero   ! Eccentricidad
            particles_in(1,3) = cero   ! M [deg]
            particles_in(1,4) = cero   ! w [deg]
            particles_in(1,5) = 8.1d0  ! MMR
            
        end if

        ! Interactions

        !! Collision y escapes
        input%min_distance = -1.5d0    ! Min distance before impact [km] ! 0 => R0 + max(Rboul)
        input%max_distance = -1.d2     ! Max distance before escape [km] ! -x => R0 * x

        !! Merges
        input%use_merge_part_mass = .True. ! Merge particles into asteroid
        input%use_merge_massive = .True. ! Merge massive bodies

        !!! Collisions factors
        input%eta_col = uno  ! 0: Elastic, 1: Plastic
        input%f_col = uno  ! Bounded: Etot < -f |Epot|

        !! Stops
        input%use_stop_no_moon_left = .True. ! Stop if no more moons left
        input%use_stop_no_part_left = .True. ! Stop if no more particles left

        ! Additional forces

        !! Manual J2 (for primary)
        input%use_manual_J2 = .False.
        input%manual_J2 = cero
        
        !! Stokes
        input%use_stokes = .False.
        input%stokes_a_damping_time = infinity                           ! [day]
        input%stokes_e_damping_time = input%stokes_a_damping_time / 1.d2 ! [day]
        input%stokes_active_time = cero                                  ! [day] Tiempo que actúa stokes

        !! Naive-Stokes (drag)
        input%use_drag = .False.
        input%drag_coefficient = cero  ! Eta
        input%drag_active_time = cero  ! [day] Tiempo que actúa drag


        ! Filter
        input%use_filter = .False.
        input%filter_dt = cero     ! dt to cutoff
        input%filter_nsamples = 0  ! Amount of integration steps of dt ~ df_fitered, per window
        input%filter_nwindows = 0  ! Amount of windows to compute the filtering
        input%filter_prefix = ""   ! empty => filt


        ! Simulation parameters

        !! Integrator
        input%integrator_ID = 0        ! Integrator to use. (0 = BS; see src/integrators/README.md)

        !! Times
        input%final_time = 2.d3        ! Final time [day]
        input%case_output_type = 0     ! 0: Linear ; 1: Logarithmic ; 2: Combination
        input%output_timestep = cero   ! Output timestep [day] (used if case_output_type != 1)
        input%output_number = 100      ! Number of outputs (if output_timestep = 0)
        input%extra_checkpoints = 2000 ! Number of extra checkpoints for chaos calculations

        !! Error and Tiemstepping
        input%use_adaptive = .True.    ! Whether to use an adaptive timestepping method
        input%dt_min = -1.d-2          ! Minimum dt to use. (if negative, |dt| * min_period)
        input%learning_rate = 0.85d0   ! [For adaptive step integrators] Learning rate
        input%error_digits = 12        ! [For adaptive step integrators] Digits for relative error

        ! Map
        input%use_potential_map = .False.
        input%mapfile = ""
        input%map_grid_size_x = 500
        input%map_grid_size_y = 500
        input%map_min_x = -500.d0
        input%map_max_x = 500.d0
        input%map_min_y = -500.d0
        input%map_max_y = 500.d0

        ! Output: "" or "no", if not used
        input%datafile = ""
        input%chaosfile = ""
        input%multfile = ""
        input%geometricfile = ""

        !! Screen
        input%use_screen = .True. ! Print info in screen
        input%use_datascreen = .True. ! Print data in screen

        ! Input: "" or "no", if not used
        input%tomfile = ""
        input%moonsfile = ""
        input%particlesfile = ""

        ! Parallel
        input%use_parallel = .False.
        input%requested_threads = 1 ! Number of threads to use !! -1 => all available

    end if


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!! No tocar de aquí a abajo !!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!! EXTRA ARGUMENTS AND DERIVED !!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    

    ! Load command line arguments
    call load_command_line_arguments(input, use_configfile)
    
    ! Message
    if ((arguments_number .le. 1) .and. (use_configfile .and. (.not. configfile_exists))) then  ! Global variable
        print*, "WARNING: No command line arguments."
        print*, "¿Do you want to run with in-code set parameters? [y/N]"
        read(*,*) aux_character20
        aux_character20 = adjustl(aux_character20)
        aux_character1 = aux_character20(1:1)   
        if (index("YySs", aux_character1) == 0) then
            write(*,*) "Exiting."
            stop 1
        else if (sim%use_screen) then
            write(*,*) "Resuming..."
            write(*,*) ACHAR(5)
        end if
    end if
    
    ! Init derived parameters
    sim%input_params_st = input ! Create sim with input parameters
    call set_derived_parameters(sim) ! Inicializamos parámetros derivados
    if (sim%use_screen) unit_file = 6  ! Unit_file to std out (for collisions and escapes)
    
    ! <> Check Output

    ! Mensaje
    if (sim%use_screen .and. .not. configfile_exists) then
        write(*,*) "WARNING: Configurations file could not be read: config.ini"
        write(*,*) "         Using in-code parameters."
    end if

    if (.not. any((/sim%use_datascreen, sim%use_datafile, sim%use_chaosfile, sim%use_geometricfile, &
                  & sim%use_potential_map, sim%use_multiple_outputs/))) then ! No tiene sentido hacer nada
        if (.not. sim%use_screen) then
            write(*,*) ACHAR(10)
            write(*,*)  "EXITING: No output to be generated."
            stop 1
        else 
            write(*,*)  "WARNING: No integration to be done. Just printing information on screen."
            ! Re-set only_print if not set
            sim%only_print = .True.
        end if
    end if

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!! READ BODIES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    
    ! <> Check Bodies

    !! Read Moons or Particles if needed

    !!! Moons
    if (sim%use_command_body .and. sim%cl_body_is_moon) then
        if (allocated(moons_in)) then  ! Deallocate if necessary
            if (sim%use_screen) then
                write(*,*) ACHAR(10)
                write(*,*) "WARNING: Command line moon used. Ignoring moons in config file."
            end if
            deallocate(moons_in)
        else if (sim%use_moonsfile .and. sim%use_screen) then
            write(*,*) ACHAR(10)
            write(*,*) "WARNING: Command line moon used. Ignoring moons file: ", trim(sim%moonsfile)
        end if
        call allocate_params_moons(1)
        ! Fill the array
        moons_in(1,:) = cl_body_in
    else if (sim%use_moonsfile) then
        if (sim%Nmoons > 0) then
            write(*,*) ACHAR(10)
            write(*,*) "ERROR: Can not read moons from moons file and config file."
            stop 1
        end if
        if (sim%use_screen) write(*,*) "Reading moons from file: ", trim(sim%moonsfile)
        !!!! Formato: [mu, a, e, M, w, MMR, radius]
        call read_columns_file(sim%moonsfile, aux_2D)
        sim%Nmoons = size(aux_2D, 1)
        if (sim%Nmoons == 0) then
            write(*,*) ACHAR(10)
            write(*,*) "ERROR: Could not read moons from moons file."
            stop 1
        end if
        sim%use_moons = .True.
        aux_int = size(aux_2D, 2)
        if (sim%use_screen) then
            write(*,s1i1) "  Read ", sim%Nmoons, " rows."
            write(*,s1i1) "  Read ", aux_int, " columns."
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
                write(*,*) ACHAR(10)
                write(*,*) "WARNING: Command line particle used. Ignoring particles in config file."
            end if
            deallocate(particles_in)
        else if (sim%use_particlesfile .and. sim%use_screen) then
            write(*,*) ACHAR(10)
            write(*,*) "WARNING: Command line particle used. Ignoring particles file: ", trim(sim%particlesfile)
        end if
        call allocate_params_particles(1)
        ! Fill the array
        particles_in(1,:) = cl_body_in(2:6)
    else if (sim%use_particlesfile) then
        if (sim%Nparticles > 0) then
            write(*,*) ACHAR(10)
            write(*,*) "ERROR: Can not read particles from particles file and config file."
            stop 1
        end if
        if (sim%use_screen) write(*,*) "Reading particles from file: ", trim(sim%particlesfile)
        !!!! Formato: [mu, a, e, M, w, MMR]
        call read_columns_file(sim%particlesfile, aux_2D)
        sim%Nparticles = size(aux_2D, 1)
        if (sim%Nparticles == 0) then
            write(*,*) ACHAR(10)
            write(*,*) "ERROR: Could not read particles from particles file."
            stop 1
        end if
        sim%use_particles = .True.
        aux_int = size(aux_2D, 2)
        if (sim%use_screen) then
            write(*,s1i1) "  Read ", sim%Nparticles, " rows."
            write(*,s1i1) "  Read ", aux_int, " columns."
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
        if (.not. (sim%use_potential_map .or. sim%only_print)) then
            write(*,*) ACHAR(10)
            write(*,*)  "EXITING: No moons/particles to integrate."
            stop 1
        else if (sim%use_potential_map .and. .not. sim%only_print) then
            only_potential_map = .True.  ! Global
            if (sim%use_screen) then
                write(*,*)  "WARNING: Performing only potential map calculation."
            write(*,*) ACHAR(5)
            end if
        end if
    end if

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!  Parallel  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (sim%use_screen) then
        write(*,*) ACHAR(5)
        if (sim%use_parallel) then
            write(*,*) "---------- Parallel integration----------"
            write(*,s1i1) "  Available threads:", available_threads
            write(*,s1i1) "  Requested threads:", sim%requested_threads
            write(*,s1i1) "  Used threads     :", my_threads
        else
            write(*,*) "Serial integration (No parallel)"
        end if
        write(*,*) ACHAR(5)
    end if
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!! CONFIGURATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (sim%use_screen) then
        write(*,s1i1) "Configurating simulation number: ", simulation_number
        write(*,*) ACHAR(5)
        write(*,*) "---------- Initial parameters ----------"
        write(*,*) ACHAR(5)
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
                    & sim%lambda_kep, &                           ! keplerian omega
                    & sim%asteroid_rotational_period * unit_time) ! asteroid period

    call set_system_extra(system, cero, sim%eta_col, sim%f_col, sim%manual_J2)  ! Extra parameters


    ! <<<< Save initial data >>>>
    initial_system = system

    ! <<<< Fill auxiliar arrays >>>>
    allocate(boulders_data(0:sim%Nboulders, 4))   ! To use in dydt
    allocate(boulders_coords(0:sim%Nboulders, 4)) ! To use in dydt

    !! Primary
    boulders_coords(0,:) = system%asteroid%primary%coordinates_CM + system%asteroid%coordinates
    boulders_data(0,1) = system%asteroid%primary%mass
    boulders_data(0,2) = system%asteroid%primary%radius
    boulders_data(0,3) = system%asteroid%primary%initial_theta
    boulders_data(0,4) = system%asteroid%primary%dist_to_asteroid

    !! Boulders
    do i = 1, sim%Nboulders
        boulders_coords(i,:) = system%asteroid%boulders(i)%coordinates_CM + system%asteroid%coordinates
        boulders_data(i,1) = system%asteroid%boulders(i)%mass
        boulders_data(i,2) = system%asteroid%boulders(i)%radius
        boulders_data(i,3) = system%asteroid%boulders(i)%initial_theta
        boulders_data(i,4) = system%asteroid%boulders(i)%dist_to_asteroid
    end do


    !! <<<< Messages >>>>
    
    ! Asteroid
    if (sim%use_screen) then
        if ((.not. sim%use_boulders) .and. (.not. sim%use_triaxial)) then
            write(*,*) "WARNING: No boulders or triaxial shape to integrate. "
            write(*,*) "         Rotation is considered just to define initial parameters."
            write(*,*) ACHAR(5)
        end if
        write(*,*) "Asteroid Rotation:"
        write(*,s1r1) "  a_corot               :", system%asteroid%a_corotation / unit_dist, "[km]"
        write(*,s1r1) "  omega                 :", system%asteroid%omega * unit_time, "[rad day⁻¹]"
        write(*,s1r1) "  omega_kep             :", system%asteroid%omega_kep * unit_time, "[rad day⁻¹]"
        write(*,s1r1) "  lambda_kep            :", system%asteroid%lambda_kep
        write(*,s1r1) "  Period                :", system%asteroid%rotational_period * 24 * unit_time, "[hs]"
        write(*,s1r1) "  Inertial momentum     :", system%asteroid%inertia / (unit_mass * unit_dist**2), "[kg km²]"
        write(*,s1r1) "  Angular momentum (rot):", system%asteroid%ang_mom_rot / unit_angm, "[kg km² day⁻¹]"
        write(*,*) ACHAR(5)
        write(*,s1r1) " Asteroid mass: ", asteroid%mass / unit_mass, "[kg]"
        write(*,*) ACHAR(5)
        
        !! Tri-axial
        if (sim%use_triaxial) then
            write(*,*) "Asteroid ellipsoid:"
            write(*,s1r1) "  semi-axis-a :", system%asteroid%primary%semi_axis(1), "[km]"
            write(*,s1r1) "  semi-axis-b :", system%asteroid%primary%semi_axis(2), "[km]"
            write(*,s1r1) "  semi-axis-c :", system%asteroid%primary%semi_axis(3), "[km]"
            write(*,s1r1) "   Effective Radius :", system%asteroid%primary%radius, "[km]"
            write(*,s1r1) "   C_20:", system%asteroid%primary%C20
            write(*,s1r1) "   C_22:", system%asteroid%primary%C22
            write(*,*) ACHAR(5)
        else
            write(*,s1r1) " Asteroid radius :", system%asteroid%radius, "[km]"
            write(*,*) ACHAR(5)
        end if

    end if


    ! < Messages >
    if (sim%use_screen) then
        write(*,*) "System configuration:"
        write(*,*) ACHAR(5)
    end if


    ! Individual boulder positions
    if (sim%use_screen .and. sim%use_boulders) then
        write(*,*) "   Asteroid-centric boulders: [x, y, vx, vy, mass, distance, radius]"
        do i = 0, sim%Nboulders  ! aux_arrays include primary
            write(*,i1r7) i, &
                & boulders_coords(i,1:2) / unit_dist, &
                & boulders_coords(i,3:4) / unit_vel, &
                & boulders_data(i,1) / unit_mass, &
                & boulders_data(i,4) / unit_dist, &
                & boulders_data(i,2) / unit_dist
        end do
        write(*,*) ACHAR(5)
    end if


    ! Moons positions
    if (sim%use_screen .and. sim%use_moons) then
        write(*,*) "   Barycentric moons: [x, y, vx, vy, mass, distance]"
        do i = 1, sim%Nmoons
            write(*,i1r7) i, &
                & system%moons(i)%coordinates(1:2) / unit_dist, &
                & system%moons(i)%coordinates(3:4) / unit_vel, &
                & system%moons(i)%mass / unit_mass, &
                & system%moons(i)%dist_to_cm / unit_dist
        end do
        write(*,*) ACHAR(5)
    end if


    ! Particles positions
    if (sim%use_screen .and. sim%use_particles) then
        write(*,*) "   Barycentric particles: [x, y, vx, vy, distance]"
        do i = 1, sim%Nparticles
            write(*,i1r7) i, &
                     & system%particles(i)%coordinates(1:2) / unit_dist, &
                     & system%particles(i)%coordinates(3:4) / unit_vel, &
                     & system%particles(i)%dist_to_cm / unit_dist
        end do
        write(*,*) ACHAR(5)
    end if


    ! Energy and ang_mom
    if (sim%use_screen) then
        write(*,s1r1) "  Total mass            : ", system%mass / unit_mass, "[kg]"
        write(*,s1r1) "  Total energy          : ", system%energy / unit_ener * (metro**2 / segundo), "[J]"
        write(*,s1r1) "  Total angular momentum: ", system%ang_mom / unit_angm * (metro**2 / segundo**2), "[J s⁻¹]"
        write(*,*) ACHAR(5)
    end if



    ! <> Parallel revisted
    ! Redefine OMP threads if necessary. Lower equal to particle number
    if ((sim%use_parallel) .and. ((sim%Ntotal - 1) < my_threads)) then
        my_threads = sim%Ntotal - 1
        if (sim%use_screen) then
            write(*,s1i1) "WARNING: Threads reduced to ", my_threads
            write(*,*) ACHAR(5)
        end if
        !$ call omp_set_num_threads(my_threads)
    end if


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!! EXTRA EFFECTS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !! Initial message
    if (sim%use_screen) then
        write(*,*) ACHAR(5)
        write(*,*) ("---------- External / Internal forces ----------")
        write(*,*) ACHAR(5)
    end if


    !! <<<< Tri-axial gravity >>>>
    if (sim%use_triaxial) then
        call init_ellipsoid(sim%triax_a_primary * unit_dist, sim%triax_b_primary * unit_dist, sim%triax_c_primary * unit_dist)
    end if

    !! <<<< Manual J2 >>>>
    if (sim%use_manual_J2) then
        call init_manual_J2(sim%manual_J2, system%asteroid%primary%radius)
        if (sim%use_screen) then
            write(*,*) "Manual J2"
            write(*,s1r1) "  J2 : ", sim%manual_J2
            write(*,*) ACHAR(5)
        end if
        any_extra_effect = .True. 
    end if


    !! <<<< Moons gravity >>>>
    if ((.not. sim%use_moon_gravity) .and. sim%Nmoons > 0) then
        if (sim%use_screen) then
            write(*,*) "Gravity between moons DEACTIVATED."
            write(*,*) ACHAR(5)
        end if
        any_extra_effect = .True.
    end if


    !! <<<< Omega Damping >>>>
    if (sim%use_omega_damping .and. ((.not. sim%use_boulders) .and. (.not. sim%use_triaxial))) then
        if (sim%use_screen) write(*,*) "WARNING: Skipping omega damping without boulders."
        sim%use_omega_damping = .False.
    else if (sim%use_omega_damping) then
        if (sim%use_lin_omega_damp) then
            aux_real = - system%asteroid%omega / (sim%omega_lin_damping_time * unit_time)
            call init_damping(aux_real, &
                            & cero, &
                            & sim%omega_damp_active_time * unit_time, &
                            & 1)
            if (sim%use_screen) then
                write(*,*) "Omega Damping: linear"
                write(*,s1r1) "  tau_o :", sim%omega_lin_damping_time / (system%asteroid%rotational_period / unit_time), "[Prot]"
                write(*,*) ACHAR(5)
            end if
        else if (sim%use_exp_omega_damp) then
            call init_damping(sim%omega_exp_damping_time * unit_time, &
                            & cero, &
                            & sim%omega_damp_active_time * unit_time, &
                            & 2)
            if (sim%use_screen) then
                write(*,*) "Omega Damping: exponential"
                write(*,s1r1) "  tau_o :", sim%omega_exp_damping_time / (system%asteroid%rotational_period / unit_time), "[Prot]"
                write(*,*) ACHAR(5)
            end if
        else if (sim%use_poly_omega_damp) then
            call init_damping(sim%omega_exp_damp_poly_A, &
                            & sim%omega_exp_damp_poly_A, &
                            & sim%omega_damp_active_time * unit_time, &
                            & 3)
            
            if (sim%use_screen) then
                write(*,*) "Omega Damping: poly-exponential"
                write(*,s1r1) "  A :", sim%omega_exp_damp_poly_A
                write(*,s1r1) "  B :", sim%omega_exp_damp_poly_B
                write(*,*) ACHAR(5)
            end if
        else 
            write(*,*) ACHAR(10)
            write(*,*) "ERROR: Can not identify omega damping model."
            stop 1
        end if
        any_extra_effect = .True.        
    end if


    !! <<<< Stokes >>>>
    if (sim%use_stokes) then
        call init_stokes(sim%stokes_a_damping_time * unit_time, &
                       & sim%stokes_e_damping_time * unit_time, &
                       & sim%stokes_active_time * unit_time)
        if (sim%use_screen) then
            write(*,*) "Stokes"
            write(*,s1r1) "  tau_a   : ", sim%stokes_a_damping_time, "[days]"
            write(*,s1r1) "  tau_e   : ", sim%stokes_e_damping_time, "[days]"
            write(*,s1r1) "  t_stokes: ", sim%stokes_active_time, "[days]"
            write(*,*) ACHAR(5)
        end if
        any_extra_effect = .True. 
    end if


    !! <<<< Naive-Stokes (Drag) >>>>
    if (sim%use_drag) then
        call init_drag(sim%drag_coefficient, sim%drag_active_time * unit_time)
        if (sim%use_screen) then
            write(*,*) "Drag"
            write(*,s1r1) "  eta   : ", sim%drag_coefficient
            write(*,s1r1) "  t_drag: ", sim%drag_active_time, "[days]"
            write(*,*) ACHAR(5)
        end if
        any_extra_effect = .True. 
    end if


    
    ! Final message
    if (sim%use_screen .and. .not. any_extra_effect) then
        write(*,*) "No extra internal / external effects activated."
        write(*,*) ACHAR(5)
    end if


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!! ESCAPES /COLLISIONS !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !! Initial message
    if (sim%use_screen) then
        write(*,*) ACHAR(5)
        write(*,*) ("---------- Escapes / Collisions parameters ----------")
        write(*,*) ACHAR(5)
    end if


    ! <<<< Escape/Colisión >>>>
    if (sim%min_distance < cero) then
        !! Check
        if (abs(sim%min_distance) < uno) then
            write(*,*) ACHAR(10)
            write(*,s1r1) "ERROR: rmin can not be lower than asteroid radius. Ratio:", abs(sim%min_distance)
            stop 1
        end if
        sim%min_distance = system%asteroid%radius * abs(sim%min_distance)
    else if (sim%min_distance > tini) then
        sim%min_distance = max(sim%min_distance * unit_dist, system%asteroid%radius)
    else  ! Assume its cero
        sim%min_distance = cero
    end if
    
    if (sim%max_distance < cero) then
        sim%max_distance = system%asteroid%radius * abs(sim%max_distance)
    else 
        sim%max_distance = sim%max_distance * unit_dist
    end if
    if ((sim%max_distance > cero) .and. (sim%max_distance <= sim%min_distance)) then
        write(*,*) ACHAR(10)
        write(*,*) "ERROR: rmax <= rmin"
        stop 1
    end if
    if (sim%use_screen) then
        write(*,*) "Conditions for escape / collision"
        if (sim%min_distance > tini) then
            write(*,s1r1) "   rmin : ", sim%min_distance / unit_dist, "[km] =", sim%min_distance / system%asteroid%radius, "[Rast]"
        else
            write(*,s1r1) "   rmin : Asteroid surface ~", system%asteroid%radius / unit_dist, "[km]"
        end if
        if (sim%max_distance > cero) then
            write(*,s1r1) "   rmax : ", sim%max_distance / unit_dist, "[km] =", sim%max_distance / system%asteroid%radius, "[Rast]"
        else
            write(*,*) "   rmax : Infinity"
        end if
        write(*,*) ACHAR(5)
        if (sim%use_particles) then
            if (sim%use_merge_part_mass) then
                write(*,*) "Colliding particles into massive bodies will be removed."
            else
                write(*,*) "Colliding particles into massive bodies will stop the integration."
            end if
            write(*,*) ACHAR(5)
        end if
        if (sim%use_moons) then
            if (sim%eta_col == uno) then
                if (sim%use_merge_massive) then
                    write(*,*) "Colliding massive bodies will be merged."
                else
                    write(*,*) "Colliding massive bodies will stop the integration."
                end if
            else
                write(*,s1r1) " Massive bodies collisional factor:", sim%eta_col
            end if
            write(*,*) ACHAR(5)
        end if
        if (sim%use_any_stop) then
            if (sim%use_stop_no_part_left .and. sim%use_particles) then
                write(*,*) " Simulation will stop if no more particles are left."
            end if
            if (sim%use_stop_no_moon_left .and. sim%use_moons) then
                write(*,*) " Simulation will stop if no more moons are left."
            end if
        else
            write(*,*) "Simulation will stop if no more particles and moons are left."
        end if
        write(*,*) ACHAR(5)
    end if


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! MAPS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !! Initial message and configuration
    if (sim%use_potential_map .and. .not. sim%only_print) then
        if (sim%use_screen) then
            write(*,*) ACHAR(5)
            write(*,*) "---------- MAP (Gravitational Energy and Acceleration) ----------"
            write(*,*) ACHAR(5)
            write(*,*) "Creating maps..."
        end if
        call create_map(sim, system)
        if (sim%use_screen) then
            write(*,*) "Maps saved to file: ", trim(sim%mapfile)
            write(*,*) ACHAR(5)
        end if
    end if

    !! Global configuration
    if (only_potential_map .and. .not. sim%only_print) then  ! Global
        if (sim%use_screen) then
            write(*,*) "END: Only the Map was requested."
            write(*,*) ACHAR(5)
        end if
        stop 0
    end if


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!! INTEGRACIÓN !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Initial message
    if (sim%use_screen) then
        write(*,*) ACHAR(5)
        write(*,*) "---------- Times ----------"
        write(*,*) ACHAR(5)
    end if


    !!!!!!!! TIEMPOS !!!!!!!

    ! <<<< Integration times >>>>
    if (sim%final_time < cero) then
        if (abs(system%asteroid%rotational_period) < tini) then
            write(*,*) "WARNING: Setting total integration time with non-rotating asteroid."
            write(*,*) ACHAR(5)
        end if
        sim%final_time = abs(sim%final_time) * system%asteroid%rotational_period
    else
        sim%final_time = sim%final_time * unit_time
    end if

    ! Check final time
    if (sim%final_time .le. cero) then
        write(*,*) ACHAR(10)
        write(*,*) "ERROR: tf < 0"
        stop 1
    end if


    !! <<<< Timesteps >>>>
    sim%output_timestep = sim%output_timestep * unit_time


    !! <<<< Output times >>>>
    call set_output_times(cero, sim%final_time, sim%output_number, &
                        & sim%output_timestep, sim%case_output_type, output_times)


    !! <<<< TOMFILE >>>>
    if (sim%use_tomfile) then

        call setup_TOM(tom, sim%tomfile, sim%final_time / unit_time, sim%use_screen)
        ! Warning if delta_omega with no rotation
        if (tom%use_domega .and. ((.not. sim%use_boulders) .and. (.not. sim%use_triaxial))) then
            if (sim%use_screen) then
                write(*,*) "WARNING: Delta_omega has no sense without boulders or triaxial elipsoid. It will be ignored."
            end if
            deallocate(tom%deltaomega)  ! Deallocate
            tom%use_domega = .False.
        end if
        ! SET UNITS
        if (tom%total_number > 0) tom%times = tom%times * unit_time
        if (tom%use_domega) tom%deltaomega = tom%deltaomega / unit_time
        if (tom%use_dmass) tom%deltamass = tom%deltamass * unit_mass

        !! Ahora debemos combinar los tiempos de TOM con los tiempos de Output, y crear un nuevo vector de tiempos Checkpoints
        call merge_sort_and_unique(tom%times, output_times, &
                                   & checkpoint_is_tom, checkpoint_is_output, &
                                   & checkpoint_times, sim%checkpoint_number)

    else

        !! En este caso, los tiempos de check son los mismos que los de LOOP
        sim%checkpoint_number = sim%output_number

        allocate(checkpoint_times(sim%checkpoint_number))
        allocate(checkpoint_is_output(sim%checkpoint_number))
        allocate(checkpoint_is_tom(sim%checkpoint_number))

        checkpoint_times = output_times

        checkpoint_is_output = .True.
        checkpoint_is_tom = .False.

    end if


    !! <<<< Extra checkpoint times >>>>
    if (sim%extra_checkpoints > 0) then

        !! Ensure extra_checkpoints is not the same as output_number
        aux_int = sim%extra_checkpoints + 2   ! (+2 to account fo initial and final times, and +1 for ???)
        if (sim%output_number - aux_int == 1) aux_int = aux_int + 2  ! Add 2 just to make a different array
        if (MOD(aux_int, sim%output_number) == 0) aux_int = aux_int + 1  ! Add 1 just to make a different array
        if (aux_int - sim%output_number == 1) aux_int = aux_int + 1  ! Add 1 just to make a different array

        !! Create more checkpoint times than only the output+tom ones
        call set_output_times(cero, sim%final_time, aux_int, aux_real, sim%case_output_type, aux_1D)

        !! Ahora debemos combinar los checkpoints previos con los nuevos, y crear un nuevo vector de tiempos Checkpoints
        call expand_checkpoints(aux_1D, checkpoint_times, checkpoint_is_output, checkpoint_is_tom, sim%checkpoint_number)

        ! Message
        if (sim%use_screen) write(*,s1i1) "Checkpoint times amount expanded. Total:", sim%checkpoint_number

    end if


    ! <<<<  Variables temporales (para integración) >>>>
    time = cero ! Tiempo actual
    timestep = checkpoint_times(1) ! Paso de tiempo inicial

    ! Get minimum period (already has units)
    sim%min_period = infinity
    if (system%asteroid%rotational_period > tini) sim%min_period = system%asteroid%rotational_period
    do i = 1, system%Nmoons_active
        sim%min_period = min(sim%min_period, get_Period(system%asteroid%mass + system%moons(i)%mass, system%moons(i)%elements(1)))
    end do
    do i = 1, system%Nparticles_active
        sim%min_period = min(sim%min_period, get_Period(system%asteroid%mass, system%particles(i)%elements(1)))
    end do

    ! Set minumum timestep  (already has units)
    if (sim%dt_min < cero) then
        sim%dt_min = sim%min_period * abs(sim%dt_min)
    else if (sim%dt_min > cero) then
        sim%dt_min = sim%dt_min * unit_time
    end if
    if (sim%dt_min < tini) sim%dt_min = sim%min_period * 1d-5  ! Heuristic

    ! Set initial adaptive and fixed timestep  (already has units)
    adaptive_timestep = sim%dt_min
    fixed_timestep = sim%dt_min

    !! Mensaje
    if (sim%use_screen) then
        write(*,*) "Times:"
        if (system%asteroid%rotational_period > tini) then
            write(*,s1r1) "  tf    : ", sim%final_time / system%asteroid%rotational_period, "[Prot] = ", &
            & sim%final_time / unit_time, "[day]"
            write(*,s1r1) "  dt_out: ", sim%output_timestep / system%asteroid%rotational_period, "[Prot] = ", &
            & sim%output_timestep / unit_time, "[day]"
            write(*,s1r1) "  dt_min: ", sim%dt_min / system%asteroid%rotational_period, "[Prot] = ", &
            & sim%dt_min / unit_time, "[day]"
        else 
            write(*,s1r1) "  tf    : ", sim%final_time / unit_time, "[day]"
            write(*,s1r1) "  dt_out: ", sim%output_timestep / unit_time, "[day]"
            write(*,s1r1) "  dt_min: ", sim%dt_min / unit_time, "[day]"
        end if
        write(*,*) ACHAR(5)
        write(*,s1i1) "  n_out         : ", sim%output_number
        write(*,s1i1) "  n_checkpoints : ", sim%checkpoint_number
        write(*,*) ACHAR(5)
    end if

    !! Possible Warning
    if (sim%dt_min > sim%output_timestep) then
        write(*,*) ACHAR(5)
        write(*,*) " WARNING: Min dt is greater than output dt."
        write(*,*) ACHAR(5)
    end if


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!! Integration Arrays !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    ! <<<< Allocate >>>>
    allocate(m_arr(sim%Ntotal))  ! Masses
    allocate(R_arr(sim%Ntotal))  ! Radii
    allocate(y_arr(2 + sim%Ntotal * 4))     ! Theta, Omega, Positions
    allocate(y_arr_new(2 + sim%Ntotal * 4)) ! Theta, Omega, Positions
    allocate(y_der(2 + sim%Ntotal * 4))     ! derivate(Theta, Omega, Positions)

    ! <<<< Init to 0 >>>>
    m_arr = cero
    R_arr = cero
    y_arr = cero
    y_arr_new = cero
    y_der = cero

    ! <<<< Amount of values of Y to use >>>>
    y_nvalues = get_index(sim%Nactive) + 3  ! Update nvalues to use in y

    ! <<<< Arrays to integrate >>>>
    call center_sytem(system)
    call generate_arrays(system, m_arr, R_arr, y_arr)


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! FILTER !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Initial message
    if (sim%use_screen) then
        write(*,*) ACHAR(5)
        write(*,*) "---------- Filter settings ----------"
        write(*,*) ACHAR(5)
    end if

    ! <<<< SET UP FILTER >>>>
    if (sim%use_filter) then
        ! Create aux system
        system_filtered = system

        ! Allocate arrays
        allocate(y_pre_filter(size(y_arr)))
        allocate(elem_filtered(size(y_arr)))
        allocate(tmp_y_arr(size(y_arr)))
        
        ! Set paramters
        if (sim%filter_dt < cero) sim%filter_dt = abs(sim%filter_dt) * twopi / system%asteroid%omega / unit_time
        call setup_filter(filter, &
                & sim%filter_dt * unit_time, &
                & sim%filter_nsamples, &
                & sim%filter_nwindows, &
                & .True., &  ! LOW PASS
                & sim%filter_model, &
                & y_nvalues)
        !! Update sim parameters
        sim%filter_dt = filter%dt
        sim%filter_nsamples = filter%n_samples
        sim%filter_nwindows = filter%n_windows
        
        ! Create filter file
        open (unit=30, file=trim(trim(sim%filter_prefix)), status='replace', action='write', position="append")
        write(30,*) "dt ", "n_samples ", "n_windows ", "size ", "total_dt"
        write(30,*) filter%dt, filter%n_samples, filter%n_windows, filter%size, filter%dt * filter%size
        write(30,*) filter%kernel
        close(30)

        ! CHECK
        if (filter%dt * filter%size > sim%final_time) then
            write(*,*) ACHAR(10)
            write(*,*) "ERROR: Filter has more window timespan than simulation final time."
            stop 1
        end if

        if (sim%use_screen) then
            write(*,*) "Filter ON"
            write(*,s1r1) "  dt    :", filter%dt / system%asteroid%omega, "[Prot] = ", &
                          & filter%dt / unit_time, "[day]"
            write(*,s1i1) "  size  :", filter%size
            write(*,s1r1) "  width :", filter%dt * filter%size / system%asteroid%omega, "[Prot] = ", &
                          & filter%dt * filter%size / unit_time, "[day]"
            write(*,s1r1) "  cut_off freq:", filter%omega_pass * radian / unit_time, "[day⁻¹] =>", &
                          & uno / (filter%omega_pass * radian / unit_time), "[day]"
            write(*,*) ACHAR(5)
        end if
    
    else if (sim%use_screen) then  ! NO FILTER
        write(*,*) "Filter OFF"
        write(*,*) ACHAR(5)
    end if


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! OUTPUT  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Initial message
    if (sim%use_screen) then
        write(*,*) ACHAR(5)
        write(*,*) "---------- OUTPUT ----------"
        write(*,*) ACHAR(5) 
    end if

    !!!! Ahora si, si no hay que hacer más nada entonces terminamos.
    if (.not. any((/sim%use_datascreen, sim%use_datafile, sim%use_chaosfile, sim%use_geometricfile, &
                  & sim%use_potential_map, sim%use_multiple_outputs/))) then ! No tiene sentido hacer nada
        write(*,*) ACHAR(10)
        write(*,*) "As nothing will be stored, there is no integration."
        write(*,*) "Exiting."
        stop 1
    end if


    !!! Check if not too much multiple files
    if (sim%use_multiple_outputs) then
        if (sim%Ntotal > 1000) then
            write(*,*) ACHAR(5)
            write(*,*) "WARNING: More than 1000 output files will be created."
            write(*,*) "         This could generate an issue."
            write(*,*) "Do you wish to continue? (y/[N])"
            read(*,*) aux_character20
            aux_character20 = adjustl(aux_character20)
            aux_character1 = aux_character20(1:1)   
            if (index("YySs", aux_character1) == 0) then
                write(*,*) "Exiting."
                stop 1
            else if (sim%use_screen) then
                write(*,*) "Resuming..."
                write(*,*) ACHAR(5)
            end if
        end if
    end if


    !!! Final message
    if (sim%use_screen) then
        if (sim%use_datafile) write(*,*) "General output file: ", trim(sim%datafile)
        if (sim%use_multiple_outputs) write(*,*) "Individuals output files: ", trim(sim%multfile) // "_*"
        if (sim%use_geometricfile) write(*,*) "Geometric output file: ", trim(sim%geometricfile)
        write(*,*) ACHAR(5)
        if (sim%use_datascreen) then
            write(*,*) "Data will be printed on screen."
        else 
            write(*,*) "Data will not be printed on screen."
        end if
        write(*,*) ACHAR(5)
    end if


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!! FINAL CHECKS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Initial message
    if (sim%use_screen) then
        write(*,*) ACHAR(5)
        write(*,*) "---------- FINAL CHECKS ----------"
        write(*,*) ACHAR(5) 
    end if

    !!! <<<< Check workers >>>>
    if (sim%use_parallel) then
        aux_int = MIN(MAX(sim%Ntotal, 3), my_threads)
        if (aux_int /= my_threads) then
            my_threads = aux_int
            !$ call OMP_SET_NUM_THREADS(my_threads)
            if (sim%use_screen) write(*,*) "WARNING: ", my_threads, " threads will be used for integration."
        end if
    end if


    !!!! Final message
    if (sim%use_screen) then
        if (sim%error_tolerance <= 1.d-16) write(*,*) " WARNING: e_tol might be too low (<= 10⁻¹⁶)"
        write(*,*) ACHAR(5)
        write(*,*) "FINAL CKECHS OK"
        write(*,*) ACHAR(5)
    end if



    ! <<<< Free initial arrays >>>>
    call free_initial_arays()



    !! <<<< CHECK IF GO ON >>>>
    if (sim%only_print) then
        if (sim%use_screen) then
            write(*,*) ACHAR(5)
            write(*,*) "Only printing was requested. No integration to perform."
            write(*,*) "Exiting."
            write(*,*) ACHAR(5)
        end if
        stop 0
    end if


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! SET INTEGRATOR
    call init_integrator(sim%integrator_ID, size(y_arr), 2, 1, sim%dt_min, sim%error_tolerance, sim%learning_rate, &
                        & .not. sim%use_adaptive)


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Integration  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!! Initial message
    if (sim%use_screen) then
        write(*,*) ACHAR(5)
        write(*,*) "---------- INTEGRATING ----------"
        write(*,*) ACHAR(5) 
    end if

    ! <<<< Inicializamos Punteros >>>>
    call define_writing_pointers(sim)

    
    ! <<<< ABRIMOS ARCHIVOS >>>>
    !     10:  config file
    !     11:  initial moons/particles (columns) file 
    !     12:  TOM file
    !     21:  data
    !     22:  chaos
    !     23:  geometric
    !     24:  chaos geometric
    !     30:  filter summary
    !     31:  filtered data
    !     32:  filtered chaos
    !     50:  potential map
    !  100+i:  individual bodies data 
    ! 3100+i:  filtered individual bodies data 


    !! Archivo de salida general
    if (sim%use_datafile) open (unit=21, file=trim(sim%datafile), status='replace', action='write', position="append")

    !! Chaos File
    if (sim%use_chaosfile) open (unit=22, file=trim(sim%chaosfile), status='replace', action='readwrite', position="append")

    !! Geometric file
    if (sim%use_geometricfile) open (unit=23, file=trim(sim%geometricfile),status='replace', action='write', position="append")

    !! Geometric chaos file
    if (sim%use_geomchaosfile) open (unit=24, file=trim(sim%geomchaosfile), status='replace', action='readwrite', position="append")

    !! Filter Files
    if (sim%use_filter) then

        !! Archivo de salida general
        if (sim%use_datafile) then
            open (unit=31, &
                & file=trim(sim%filter_prefix) // trim(sim%datafile), &
                & status='replace', action='write', position="append")
        end if

        !! Chaos File
        if (sim%use_chaosfile) then
            open (unit=32, &
                & file=trim(sim%filter_prefix) // trim(sim%chaosfile), &
                & status='replace', action='readwrite', position="append")
        end if

        ! !! Geometric File
        ! if (sim%use_geometricfile) then
        !     open (unit=33, &
        !         & file=trim(sim%filter_prefix) // trim(sim%geometricfile), &
        !         & status='replace', action='readwrite', position="append")
        ! end if

        !! Archivos individuales
        if (sim%use_multiple_outputs) then
            do i = 0, sim%Ntotal  ! 0 is the asteroid
                write (aux_character20, *) i
                open (unit=3100+i, &
                    & file=trim(sim%filter_prefix) // trim(sim%multfile) // "_" // trim(adjustl(aux_character20)), &
                    & status='replace', action='write', position="append")
            end do
        end if

    end if

    !! Archivos individuales
    if (sim%use_multiple_outputs) then
        do i = 0, sim%Ntotal  ! 0 is the asteroid
            write (aux_character20, *) i
            open (unit=100+i, file=trim(sim%multfile) // "_" // trim(adjustl(aux_character20)), &
                & status='replace', action='write', position="append")
        end do
    end if


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    !!!!!! MAIN LOOP INTEGRATION !!!!!!!
    keep_integrating = .True.  ! Init flag

    ! CHECK INITIAL CONDITIONS
    ! Apply colissions/escapes and checks
    call check_esc_and_col(system, unit_file)

    ! Update Nactive and y_arr if necessary
    call get_Nactive(system, new_Nactive)
    if (new_Nactive < sim%Nactive) then
        call update_sim_Nactive(sim, system%Nparticles_active, system%Nmoons_active)  ! Update sim Nactive
        y_nvalues = get_index(new_Nactive) + 3  ! Update nvalues to use in y
        call generate_arrays(system, m_arr, R_arr, y_arr)  ! Regenerate arrays
    end if
    
    ! Update Chaos (triggers update elements)
    call update_chaos(system, sim%use_baryc_output)

    ! Update geometric chaos (and elements) if needed
    if (sim%use_geometricfile) call update_chaos_geometric(system)

    !! Reduce threads if necessary...
    if (sim%use_parallel) then
        aux_int = MIN(MAX(sim%Nactive, 3), my_threads)
        if (aux_int /= my_threads) then
            my_threads = aux_int
            !$ call OMP_SET_NUM_THREADS(my_threads)
            if (sim%use_screen) then
                write(*,s1i1) " WARNING: Paralelism reduced to ", my_threads, " threads for the integration."
                write(*,*) ACHAR(5)
            end if
        end if
    end if

    ! Output initial conditions
    call generate_output(system)

    if (sim%use_filter) then
        call copy_objects(system, system_filtered)
        call generate_output(system_filtered, .True.)
    end if

    ! >>>>>>>>>>>>>>>>>>----------------- MAIN LOOP  ---------------<<<<<<<<<<<<<<<<<<<<<<<<
    tom%index_number = 2  !!!! Inicializamos en 2 porque el primer checkpoint es el IC (t0)
    j = 2  ! From 2 because 1 is the IC (t0) !! +1 por si hay HardExit en el último

    ! Check if filter used
    if (sim%use_filter) then

        ! Identify which will have filter, and which don't
        do first_idx_yes_filter = 2, sim%checkpoint_number
            ! If too short, cycle to next
            if (checkpoint_is_output(first_idx_yes_filter) .and. &
             & (checkpoint_times(first_idx_yes_filter) .ge. filter%half_width)) exit
        end do

        ! This is first initial time related to filtering
        next_t_filt = checkpoint_times(first_idx_yes_filter) - filter%half_width

        ! Find the last checkpoint before first filtering process
        do last_idx_no_filter = first_idx_yes_filter - 1, 1, -1  ! Backwards look
            if (checkpoint_is_output(first_idx_yes_filter) .and. &
             & (checkpoint_times(last_idx_no_filter) .le. next_t_filt)) exit
        end do

        ! Define last output done
        last_output = 1  ! The initial condition

        ! -----------------------------------------------------------------------------
        ! --------> LOOP WITH NO FILTER, up to last checkpoint without filter <--------
        loop_pre_filter: do while (j .le. last_idx_no_filter)
        
            hard_exit = .False. ! Resetear HardExit 
            is_premature_exit = .False. ! Resetear premature

            ! Check if Particles/Moons left
            if (sim%Nactive == 1) then

                ! Set flag
                keep_integrating = .False.

                if (sim%use_screen) then
                    write(*,*) ACHAR(5) 
                    write(*,*) " No more active particles/moons left."
                    write(*,*) ACHAR(5) 
                    write(*,s1r1) "Integration finished at time = ", time / unit_time, "[days]"
                end if

                exit loop_pre_filter

            end if
            
            ! Update dt
            next_time = checkpoint_times(j)  ! Get next time
            timestep = next_time - time

            ! INTEGRATE
            call integrate (time, y_arr(:y_nvalues), adaptive_timestep, dydt, timestep, y_arr_new(:y_nvalues), check_func)
        
            ! Check if it might be hard_exit
            if (hard_exit) then

                !! If so, the dt used is in dt_adap
                timestep = adaptive_timestep

                !! Check if premature_exit: If exit at less than 1% of finishing timestep
                if (timestep > cero) then
                    is_premature_exit = (next_time - (time + timestep)) / timestep > 0.01d0
                else  ! Weird case of timestep 0
                    is_premature_exit = next_time - time < tini
                end if

            end if

            ! Update parameters
            time = time + timestep
            y_arr(:y_nvalues) = y_arr_new(:y_nvalues)
            y_arr(1) = modulo(y_arr(1), twopi)  ! Modulate theta

            ! Update from y_new
            call update_system_from_array(system, time, y_arr)

            ! Apply escapes and colissions and check
            call check_esc_and_col(system, unit_file)

            ! Update Nactive and y_arr if necessary
            call get_Nactive(system, new_Nactive)
            if (new_Nactive < sim%Nactive) then
                call update_sim_Nactive(sim, system%Nparticles_active, system%Nmoons_active)  ! Update sim Nactive
                y_nvalues = get_index(new_Nactive) + 3  ! Update nvalues to use in y
                call generate_arrays(system, m_arr, R_arr, y_arr)  ! Regenerate arrays
            end if
        
            ! Update Chaos (triggers update elements)
            call update_chaos(system, sim%use_baryc_output)

            ! Update geometric Chaos (triggers update geometric elements)
            if (sim%use_geometricfile) call update_chaos_geometric(system)
            
            ! Output and Update j; only if not premature
            if (.not. is_premature_exit) then
                if (checkpoint_is_output(j) .and. (j > last_output)) then
                    call generate_output(system)
                    last_output = j
                end if
                j = j + 1
            end if

            ! Reset adaptive if needed
            if (.not. sim%use_adaptive) adaptive_timestep = fixed_timestep

            ! Percentage output
            if (sim%use_percentage) call percentage(time, sim%final_time)

        end do loop_pre_filter
        ! -----------------------------------------------------------------------------

        ! -----------------------------------------------------------------------------
        ! --------> LOOP WITH NO FILTER, up to first filter time <--------
        next_time = checkpoint_times(first_idx_yes_filter) - filter%half_width
        loop_until_first_filter: do while (keep_integrating)
        
            hard_exit = .False. ! Resetear HardExit 
            is_premature_exit = .False. ! Resetear premature

            ! Check if Particles/Moons left
            if (sim%Nactive == 1) then

                ! Set flag
                keep_integrating = .False.

                if (sim%use_screen) then
                    write(*,*) ACHAR(5) 
                    write(*,*) " No more active particles/moons left."
                    write(*,*) ACHAR(5) 
                    write(*,s1r1) "Integration finished at time = ", time / unit_time, "[days]"
                end if

                exit loop_until_first_filter

            end if
            
            ! Update dt
            timestep = next_time - time
            if (abs(timestep) < tini) exit loop_until_first_filter !! Manual CHECK

            ! INTEGRATE
            call integrate (time, y_arr(:y_nvalues), adaptive_timestep, dydt, timestep, y_arr_new(:y_nvalues), check_func)
        
            ! Check if it might be hard_exit
            if (hard_exit) then
                !! If so, the dt used is in dt_adap
                timestep = adaptive_timestep

                !! Check if premature_exit: If exit at less than 1% of finishing timestep
                if (timestep > cero) then
                    is_premature_exit = (next_time - (time + timestep)) / timestep > 0.01d0
                else  ! Weird case of timestep 0
                    is_premature_exit = next_time - time < tini
                end if

            end if

            ! Update parameters
            time = time + timestep
            y_arr(:y_nvalues) = y_arr_new(:y_nvalues)

            y_arr(1) = modulo(y_arr(1), twopi)  ! Modulate theta

            ! Update from y_new
            call update_system_from_array(system, time, y_arr)

            ! Apply escapes and colissions and check
            call check_esc_and_col(system, unit_file)

            ! Update Nactive and y_arr if necessary
            call get_Nactive(system, new_Nactive)
            if (new_Nactive < sim%Nactive) then

                call update_sim_Nactive(sim, system%Nparticles_active, system%Nmoons_active)  ! Update sim Nactive
                y_nvalues = get_index(new_Nactive) + 3  ! Update nvalues to use in y
                call generate_arrays(system, m_arr, R_arr, y_arr)  ! Regenerate arrays

            end if
            
            ! Update Chaos (triggers update elements)
            call update_chaos(system, sim%use_baryc_output)

            ! Update geometric Chaos (triggers update geometric elements)
            if (sim%use_geometricfile) call update_chaos_geometric(system)

            ! Reset adaptive if needed
            if (.not. sim%use_adaptive) adaptive_timestep = fixed_timestep

            ! Percentage output
            if (sim%use_percentage) call percentage(time, sim%final_time)            

        end do loop_until_first_filter
        ! -----------------------------------------------------------------------------

        ! SET TMP
        tmp_j = j
        tmp_sim = sim
        tmp_time = time
        tmp_system = system        
        tmp_y_nvalues = y_nvalues
        tmp_adaptive_timestep = adaptive_timestep
        tmp_y_arr(:y_nvalues) = y_arr(:y_nvalues)

        ! -----------------------------------------------------------------------------
        ! --------> WORK ON FIRST FILTER <--------
        ! Now we are at the first filtering data time needed
        ! This is exceptional bc may include checkpoints that may not have filters
        ! The proceeding is split in 3 steps:
        !  First: Integrate the filter using a parallel 'y_pre_filter' array
        !  Second: Integrate the checkpoints up to the principal (including). j is updated here
        !  Third: Set the array and time to next filter first time.

        ! =======> FISRT STEP <==========
        time_filt = time
        adaptive_timestep_filt = adaptive_timestep

        ! Get system and y_arr filtering
        call copy_objects(system, system_filtered)  ! From system to filtered; but no chaos (only objects)
        y_pre_filter(:y_nvalues) = y_arr(:y_nvalues)

        !! Store it
        call store_to_filter(filter, time_filt, y_pre_filter, y_nvalues, 1)

        ! Integrate and store filtered values
        do i = 2, filter%size

            !! Reset adaptive timestep if needed
            if (.not. sim%use_adaptive) adaptive_timestep_filt = fixed_timestep

            ! INTEGRATE (without check)
            call integrate (time_filt, y_pre_filter(:y_nvalues), adaptive_timestep_filt, dydt, filter%dt, y_arr_new(:y_nvalues))

            !! Update time
            time_filt = time_filt + filter%dt

            y_arr_new(1) = modulo(y_arr_new(1), twopi)  ! Modulate theta

            !! Store
            call store_to_filter(filter, time_filt, y_arr_new, y_nvalues, i)

            !! Update y
            y_pre_filter(:y_nvalues) = y_arr_new(:y_nvalues)
        
        end do

        ! Create filtered and do post-processing
        call apply_filter(y_nvalues, system, system_filtered, elem_filtered(:y_nvalues))
        call update_system_from_elements(system_filtered, &
                                        & checkpoint_times(first_idx_yes_filter), &
                                        & elem_filtered, &
                                        & sim%use_baryc_output)
        ! Update Chaos (triggers update elements)
        call update_chaos(system_filtered, sim%use_baryc_output)
        ! Filtered output
        call generate_output(system_filtered, .True.)

        ! =======> SECOND STEP <==========

        ! j is where I WILL be
        loop_missing_chkp: do while (j .le. first_idx_yes_filter)
        
            hard_exit = .False. ! Resetear HardExit 
            is_premature_exit = .False. ! Resetear premature

            ! Check if Particles/Moons left
            if (sim%Nactive == 1) then

                ! Set flag
                keep_integrating = .False.

                if (sim%use_screen) then
                    write(*,*) ACHAR(5) 
                    write(*,*) " No more active particles/moons left."
                    write(*,*) ACHAR(5) 
                    write(*,s1r1) "Integration finished at time = ", time / unit_time, "[days]"
                end if

                exit loop_missing_chkp

            end if
            
            ! Update dt
            next_time = checkpoint_times(j)  ! Get next time
            timestep = next_time - time

            ! INTEGRATE
            call integrate (time, y_arr(:y_nvalues), adaptive_timestep, dydt, timestep, y_arr_new(:y_nvalues), check_func)
        
            ! Check if it might be hard_exit
            if (hard_exit) then

                !! If so, the dt used is in dt_adap
                timestep = adaptive_timestep

                !! Check if premature_exit: If exit at less than 1% of finishing timestep
                if (timestep > cero) then
                    is_premature_exit = (next_time - (time + timestep)) / timestep > 0.01d0
                else  ! Weird case of timestep 0
                    is_premature_exit = next_time - time < tini
                end if

            end if

            ! Update parameters
            time = time + timestep
            y_arr(:y_nvalues) = y_arr_new(:y_nvalues)
            y_arr(1) = modulo(y_arr(1), twopi)  ! Modulate theta

            ! Update from y_new
            call update_system_from_array(system, time, y_arr)

            ! Apply escapes and colissions and check
            call check_esc_and_col(system, unit_file)

            ! Update Nactive and y_arr if necessary
            call get_Nactive(system, new_Nactive)
            if (new_Nactive < sim%Nactive) then
                call update_sim_Nactive(sim, system%Nparticles_active, system%Nmoons_active)  ! Update sim Nactive
                y_nvalues = get_index(new_Nactive) + 3  ! Update nvalues to use in y
                call generate_arrays(system, m_arr, R_arr, y_arr)  ! Regenerate arrays
            end if
        
            ! Update Chaos (triggers update elements)
            call update_chaos(system, sim%use_baryc_output)

            ! Update geometric Chaos (triggers update geometric elements)
            if (sim%use_geometricfile) call update_chaos_geometric(system)
            
            ! Output and Update j; only if not premature
            if (.not. is_premature_exit) then
                if (checkpoint_is_output(j) .and. (j > last_output)) then
                    call generate_output(system)
                    last_output = j
                end if
                j = j + 1
            end if

            ! Reset adaptive if needed
            if (.not. sim%use_adaptive) adaptive_timestep = fixed_timestep

            ! Percentage output
            if (sim%use_percentage) call percentage(time, sim%final_time)

        end do loop_missing_chkp

        ! =======> THIRD STEP <==========
        
        keep_integrating = keep_integrating .and. (j .le. sim%checkpoint_number)

        ! If not needed, skip
        if (keep_integrating) then

            ! -- Find the next output checkpoint index --
            next_output = j
            do i = j, sim%checkpoint_number
                if (checkpoint_is_output(i)) then
                    next_output = i
                    exit
                end if
            end do

            ! -- Compute the time at which filtering should start --
            next_t_filt = checkpoint_times(next_output) - filter%half_width

            ! -- Find the next checkpoint after that filtering time --
            aux_int = next_checkpoint
            do i = aux_int, next_output
                if (checkpoint_times(i) > next_t_filt) then
                    next_checkpoint = i
                    exit
                end if
            end do

            !! In this case, then next filter time has already passed
            if (time > next_t_filt) then

                ! We reset the system up to the tmp values
                j = tmp_j
                sim = tmp_sim
                time = tmp_time
                system = tmp_system
                y_nvalues = tmp_y_nvalues
                adaptive_timestep = tmp_adaptive_timestep
                y_arr(:y_nvalues)  = tmp_y_arr(:y_nvalues)
                aux_int = unit_file  ! For coll/escp writing
            
            else
                aux_int = -1  ! For coll/escp writing

            end if
            
            !! Now we just need to integrate the normal system up to next filter time

            !!! First, we go until last checkpoint before the filter time
            !!! No output here. At most, only extra checkpoints
            loop_until_next_checkpoint_minus_1: do while ((j < next_checkpoint) .and. keep_integrating)
            
                hard_exit = .False. ! Resetear HardExit 
                is_premature_exit = .False. ! Resetear premature

                ! Check if Particles/Moons left
                if (sim%Nactive == 1) then

                    ! Set flag
                    keep_integrating = .False.

                    if (sim%use_screen) then
                        write(*,*) ACHAR(5) 
                        write(*,*) " No more active particles/moons left."
                        write(*,*) ACHAR(5) 
                        write(*,s1r1) "Integration finished at time = ", time / unit_time, "[days]"
                    end if

                    exit loop_until_next_checkpoint_minus_1

                end if
                
                ! Update dt
                timestep = checkpoint_times(j) - time

                ! INTEGRATE
                call integrate (time, y_arr(:y_nvalues), adaptive_timestep, dydt, timestep, y_arr_new(:y_nvalues), check_func)
            
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
                y_arr(:y_nvalues) = y_arr_new(:y_nvalues)
                y_arr(1) = modulo(y_arr(1), twopi)  ! Modulate theta

                ! Update from y_new
                call update_system_from_array(system, time, y_arr)

                ! Apply escapes and colissions and check
                call check_esc_and_col(system, aux_int)

                ! Update Nactive and y_arr if necessary
                call get_Nactive(system, new_Nactive)
                if (new_Nactive < sim%Nactive) then
                    call update_sim_Nactive(sim, system%Nparticles_active, system%Nmoons_active)  ! Update sim Nactive
                    y_nvalues = get_index(new_Nactive) + 3  ! Update nvalues to use in y
                    call generate_arrays(system, m_arr, R_arr, y_arr)  ! Regenerate arrays
                end if
            
                ! Update Chaos (triggers update elements)
                call update_chaos(system, sim%use_baryc_output)

                ! Update geometric Chaos (triggers update geometric elements)
                if (sim%use_geometricfile) call update_chaos_geometric(system)

                ! Update j; only if not premature
                if (.not. is_premature_exit) j = j + 1

                ! Reset adaptive if needed
                if (.not. sim%use_adaptive) adaptive_timestep = fixed_timestep

            end do loop_until_next_checkpoint_minus_1

            !!! Now, we integrate up to next_filter
            next_time = next_t_filt
            loop_until_second_filter: do while (keep_integrating)
            
                hard_exit = .False. ! Resetear HardExit 
                is_premature_exit = .False. ! Resetear premature

                ! Check if Particles/Moons left
                if (sim%Nactive == 1) then

                    ! Set flag
                    keep_integrating = .False.

                    if (sim%use_screen) then
                        write(*,*) ACHAR(5) 
                        write(*,*) " No more active particles/moons left."
                        write(*,*) ACHAR(5) 
                        write(*,s1r1) "Integration finished at time = ", time / unit_time, "[days]"
                    end if

                    exit loop_until_second_filter

                end if
                
                ! Update dt
                timestep = next_time - time
                if (abs(timestep) < tini) exit loop_until_second_filter !! Manual CHECK

                ! INTEGRATE
                call integrate (time, y_arr(:y_nvalues), adaptive_timestep, dydt, timestep, y_arr_new(:y_nvalues), check_func)
            
                ! Check if it might be hard_exit
                if (hard_exit) then
                    !! If so, the dt used is in dt_adap
                    timestep = adaptive_timestep

                    !! Check if premature_exit: If exit at less than 1% of finishing timestep
                    if (timestep > cero) then
                        is_premature_exit = (next_time - (time + timestep)) / timestep > 0.01d0
                    else  ! Weird case of timestep 0
                        is_premature_exit = next_time - time < tini
                    end if

                end if

                ! Update parameters
                time = time + timestep
                y_arr(:y_nvalues) = y_arr_new(:y_nvalues)

                y_arr(1) = modulo(y_arr(1), twopi)  ! Modulate theta

                ! Update from y_new
                call update_system_from_array(system, time, y_arr)

                ! Apply escapes and colissions and check
                call check_esc_and_col(system, unit_file)

                ! Update Nactive and y_arr if necessary
                call get_Nactive(system, new_Nactive)
                if (new_Nactive < sim%Nactive) then

                    call update_sim_Nactive(sim, system%Nparticles_active, system%Nmoons_active)  ! Update sim Nactive
                    y_nvalues = get_index(new_Nactive) + 3  ! Update nvalues to use in y
                    call generate_arrays(system, m_arr, R_arr, y_arr)  ! Regenerate arrays

                end if
                
                ! Update Chaos (triggers update elements)
                call update_chaos(system, sim%use_baryc_output)

                ! Update geometric Chaos (triggers update geometric elements)
                if (sim%use_geometricfile) call update_chaos_geometric(system)

                ! Reset adaptive if needed
                if (.not. sim%use_adaptive) adaptive_timestep = fixed_timestep

                ! Percentage output
                if (sim%use_percentage) call percentage(time, sim%final_time)            

            end do loop_until_second_filter

        end if

        !!!! Now we are in the second filter (j = next_checkpoint)

        ! ===============================
        ! -----------------------------------------------------------------------------

        ! SET TMP
        tmp_j = j
        tmp_sim = sim
        tmp_time = time
        tmp_system = system
        tmp_y_nvalues = y_nvalues
        tmp_adaptive_timestep = adaptive_timestep
        tmp_y_arr(:y_nvalues) = y_arr(:y_nvalues)

        ! From now on, all following filtering methods are the same
        ! There are 2 arrays (systems): y_pre_filter and y_arr. They are asyncronous

        ! -----------------------------------------------------------------------------
        ! --------> LOOP WITH FILTER <--------
        ! j is where I WILL be
        loop_filter: do while (.True.)

            ! Check if all done
            !! Time end
            if (j > sim%checkpoint_number) keep_integrating = .False.

            ! Check if Particles/Moons left
            if (sim%Nactive == 1) then

                if (sim%use_screen) then
                    write(*,*) ACHAR(5) 
                    write(*,*) " No more active particles/moons left."
                end if

                keep_integrating = .False.

            end if

            ! Keep going?
            if (.not. keep_integrating) then

                if (sim%use_screen) then
                    write(*,*) ACHAR(5) 
                    write(*,s1r1) "Integration finished at time = ", time / unit_time, "[days]"
                end if

                exit loop_filter

            end if

            ! =============> Integrate and store filtered values <=============
            time_filt = time
            adaptive_timestep_filt = adaptive_timestep

            ! Get system and y_arr filtering
            call copy_objects(system, system_filtered)  ! From system to filtered; but no chaos (only objects)
            y_pre_filter(:y_nvalues) = y_arr(:y_nvalues)

            !! Store it
            call store_to_filter(filter, time_filt, y_pre_filter, y_nvalues, 1)

            ! Integrate and store filtered values
            do i = 2, filter%size  ! i = 2 bc i = 1 is where y_pre_filter already is

                !! Reset adaptive timestep if needed
                if (.not. sim%use_adaptive) adaptive_timestep_filt = fixed_timestep

                !! INTEGRATE
                call integrate (time_filt, y_pre_filter(:y_nvalues), adaptive_timestep_filt, dydt, filter%dt,&
                                &  y_arr_new(:y_nvalues))  ! No checks here

                !! Update time
                time_filt = time_filt + filter%dt

                !! Store
                call store_to_filter(filter, time_filt, y_arr_new, y_nvalues, i)

                !! Update y
                y_pre_filter(:y_nvalues) = y_arr_new(:y_nvalues)

            end do

            ! Create filtered
            call apply_filter(y_nvalues, system, system_filtered, elem_filtered(:y_nvalues))
            call update_system_from_elements(system_filtered, &
                                        & checkpoint_times(next_output), &
                                        & elem_filtered, &
                                        & sim%use_baryc_output)
            ! Update Chaos (triggers update elements)
            call update_chaos(system_filtered, sim%use_baryc_output)
            ! Filtered output
            call generate_output(system_filtered, .True.)

            ! =============> Get to Checkpoint j (with output) and process <=============
            
            ! LOOP WITH NO FILTER
            loop_until_checkp: do while ((j .le. next_output) .and. keep_integrating)
            
                hard_exit = .False. ! Resetear HardExit 
                is_premature_exit = .False. ! Resetear premature
                
                ! Check if all done
                !! Time end
                if (j == sim%checkpoint_number + 1) keep_integrating = .False.

                ! Check if Particles/Moons left
                if (sim%Nactive == 1) then

                    ! Set flag
                    keep_integrating = .False.

                    exit loop_until_checkp
                end if

                ! Set timestep
                timestep = checkpoint_times(j) - time
                
                ! INTEGRATE
                call integrate (time, y_arr(:y_nvalues), adaptive_timestep, dydt, timestep, y_arr_new(:y_nvalues), check_func)
            
            
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
                y_arr(:y_nvalues) = y_arr_new(:y_nvalues)

                y_arr(1) = modulo(y_arr(1), twopi)  ! Modulate theta

                ! Update from y_new
                call update_system_from_array(system, time, y_arr)

                ! Apply escapes and colissions and check
                call check_esc_and_col(system, unit_file)

                ! Update Nactive and y_arr if necessary
                call get_Nactive(system, new_Nactive)
                if (new_Nactive < sim%Nactive) then
                    call update_sim_Nactive(sim, system%Nparticles_active, system%Nmoons_active)  ! Update sim Nactive
                    y_nvalues = get_index(new_Nactive) + 3  ! Update nvalues to use in y
                    call generate_arrays(system, m_arr, R_arr, y_arr)  ! Regenerate arrays
                end if

                ! Update Chaos (triggers update elements)
                call update_chaos(system, sim%use_baryc_output)

                ! Update geometric Chaos (triggers update geometric elements)
                if (sim%use_geometricfile) call update_chaos_geometric(system)

                ! Update j and Output
                if (.not. is_premature_exit) then
                    if (checkpoint_is_output(j) .and. (j > last_output)) then
                        call generate_output(system)
                        last_output = j
                    end if
                    j = j + 1
                end if

                ! Reset adaptive if needed
                if (.not. sim%use_adaptive) adaptive_timestep = fixed_timestep

                ! Percentage output
                if (sim%use_percentage) call percentage(time, sim%final_time)

            end do loop_until_checkp

            ! =============> Get to next filter time <=============
            ! Now, j is already updated

            ! Check if next checkpoint exits
            if (j > sim%checkpoint_number) cycle  ! Will exit above

            ! -- Find the next output checkpoint index --
            next_output = j
            do i = j, sim%checkpoint_number
                if (checkpoint_is_output(i)) then
                    next_output = i
                    exit
                end if
            end do

            ! -- Compute the time at which filtering should start --
            next_t_filt = checkpoint_times(next_output) - filter%half_width

            ! -- Find the next checkpoint after that filtering time --
            aux_int = next_checkpoint
            do i = aux_int, next_output
                if (checkpoint_times(i) > next_t_filt) then
                    next_checkpoint = i
                    exit
                end if
            end do

            !! In this case, then next filter time has already passed
            if (time > next_t_filt) then
                ! We reset the system up to the tmp values
                j = tmp_j
                sim = tmp_sim
                time = tmp_time
                system = tmp_system
                y_nvalues = tmp_y_nvalues
                adaptive_timestep = tmp_adaptive_timestep
                y_arr(:y_nvalues)  = tmp_y_arr(:y_nvalues)
                aux_int = unit_file  ! For coll/escp writing
            
            else
                aux_int = -1  ! For coll/escp writing

            end if
                        
            !! Now we just need to integrate the normal system up to next filter time

            !!! First, we go until next checkpoint - 1
            !!! No output here. At most, only extra checkpoints
            loop2_until_next_checkpoint_minus_1: do while ((j < next_checkpoint) .and. keep_integrating)
            
                hard_exit = .False. ! Resetear HardExit 
                is_premature_exit = .False. ! Resetear premature

                ! Check if Particles/Moons left
                if (sim%Nactive == 1) then

                    ! Set flag
                    keep_integrating = .False.

                    exit loop2_until_next_checkpoint_minus_1  ! main loop

                end if
                
                ! Update dt
                timestep = checkpoint_times(j) - time

                ! INTEGRATE
                call integrate (time, y_arr(:y_nvalues), adaptive_timestep, dydt, timestep, y_arr_new(:y_nvalues), check_func)
            
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
                y_arr(:y_nvalues) = y_arr_new(:y_nvalues)
                y_arr(1) = modulo(y_arr(1), twopi)  ! Modulate theta

                ! Update from y_new
                call update_system_from_array(system, time, y_arr)

                ! Apply escapes and colissions and check
                call check_esc_and_col(system, aux_int)

                ! Update Nactive and y_arr if necessary
                call get_Nactive(system, new_Nactive)
                if (new_Nactive < sim%Nactive) then
                    call update_sim_Nactive(sim, system%Nparticles_active, system%Nmoons_active)  ! Update sim Nactive
                    y_nvalues = get_index(new_Nactive) + 3  ! Update nvalues to use in y
                    call generate_arrays(system, m_arr, R_arr, y_arr)  ! Regenerate arrays
                end if
            
                ! Update Chaos (triggers update elements)
                call update_chaos(system, sim%use_baryc_output)

                ! Update geometric Chaos (triggers update geometric elements)
                if (sim%use_geometricfile) call update_chaos_geometric(system)

                ! Update j; only if not premature
                if (.not. is_premature_exit) j = j + 1

                ! Reset adaptive if needed
                if (.not. sim%use_adaptive) adaptive_timestep = fixed_timestep

            end do loop2_until_next_checkpoint_minus_1

            !!! Now, we integrate up to next_filter
            next_time = next_t_filt
            loop2_until_second_filter: do while (keep_integrating)
            
                hard_exit = .False. ! Resetear HardExit 
                is_premature_exit = .False. ! Resetear premature

                ! Check if Particles/Moons left
                if (sim%Nactive == 1) then

                    ! Set flag
                    keep_integrating = .False.
                    
                    exit loop2_until_second_filter
                end if
                
                ! Update dt
                timestep = next_time - time
                if (abs(timestep) < tini) exit loop2_until_second_filter !! Manual EXIT

                ! INTEGRATE
                call integrate (time, y_arr(:y_nvalues), adaptive_timestep, dydt, timestep, y_arr_new(:y_nvalues), check_func)
            
                ! Check if it might be hard_exit
                if (hard_exit) then
                    !! If so, the dt used is in dt_adap
                    timestep = adaptive_timestep

                    !! Check if premature_exit: If exit at less than 1% of finishing timestep
                    if (timestep > cero) then
                        is_premature_exit = (next_time - (time + timestep)) / timestep > 0.01d0
                    else  ! Weird case of timestep 0
                        is_premature_exit = next_time - time < tini
                    end if

                end if

                ! Update parameters
                time = time + timestep
                y_arr(:y_nvalues) = y_arr_new(:y_nvalues)

                y_arr(1) = modulo(y_arr(1), twopi)  ! Modulate theta

                ! Update from y_new
                call update_system_from_array(system, time, y_arr)

                ! Apply escapes and colissions and check
                call check_esc_and_col(system, unit_file)

                ! Update Nactive and y_arr if necessary
                call get_Nactive(system, new_Nactive)
                if (new_Nactive < sim%Nactive) then

                    call update_sim_Nactive(sim, system%Nparticles_active, system%Nmoons_active)  ! Update sim Nactive
                    y_nvalues = get_index(new_Nactive) + 3  ! Update nvalues to use in y
                    call generate_arrays(system, m_arr, R_arr, y_arr)  ! Regenerate arrays

                end if
                
                ! Update Chaos (triggers update elements)
                call update_chaos(system, sim%use_baryc_output)

                ! Update geometric Chaos (triggers update geometric elements)
                if (sim%use_geometricfile) call update_chaos_geometric(system)

                ! Reset adaptive if needed
                if (.not. sim%use_adaptive) adaptive_timestep = fixed_timestep

                ! Percentage output
                if (sim%use_percentage) call percentage(time, sim%final_time)            

            end do loop2_until_second_filter

            ! Percentage output
            if (keep_integrating .and. sim%use_percentage) call percentage(time, sim%final_time)

            ! Update TMP
            tmp_j = j
            tmp_sim = sim
            tmp_time = time
            tmp_system = system
            tmp_y_nvalues = y_nvalues
            tmp_adaptive_timestep = adaptive_timestep
            tmp_y_arr(:y_nvalues) = y_arr(:y_nvalues)
        
        end do loop_filter
        ! -----------------------------------------------------------------------------

        ! Final time restoration
        if (j .eq. sim%checkpoint_number + 1) time = sim%final_time

    else

        ! LOOP WITH NO FILTER
        main_loop_no_filter: do while (.True.)
        
            hard_exit = .False. ! Resetear HardExit 
            is_premature_exit = .False. ! Resetear premature
            
            ! Check if all done
            !! Time end
            if (j == sim%checkpoint_number + 1) keep_integrating = .False.


            ! Check if Particles/Moons left
            if (sim%Nactive == 1) then
                if (sim%use_screen) then
                    write(*,*) ACHAR(5) 
                    write(*,*) " No more active particles/moons left."
                end if
                keep_integrating = .False.
            end if


            ! Keep going?
            if (.not. keep_integrating) then
                if (sim%use_screen) then
                    write(*,*) ACHAR(5) 
                    write(*,s1r1) "Integration finished at time = ", time / unit_time, "[days]"
                end if
                exit main_loop_no_filter
            end if

            
            ! Update dt
            timestep = checkpoint_times(j) - time
            
            ! INTEGRATE
            call integrate (time, y_arr(:y_nvalues), adaptive_timestep, dydt, timestep, y_arr_new(:y_nvalues), check_func)
        
        
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
            y_arr(:y_nvalues) = y_arr_new(:y_nvalues)

            y_arr(1) = modulo(y_arr(1), twopi)  ! Modulate theta

            ! Update from y_new
            call update_system_from_array(system, time, y_arr)

            ! Apply escapes and colissions and check
            call check_esc_and_col(system, unit_file)

            ! Update Nactive and y_arr if necessary
            call get_Nactive(system, new_Nactive)
            if (new_Nactive < sim%Nactive) then
                call update_sim_Nactive(sim, system%Nparticles_active, system%Nmoons_active)  ! Update sim Nactive
                y_nvalues = get_index(new_Nactive) + 3  ! Update nvalues to use in y
                call generate_arrays(system, m_arr, R_arr, y_arr)  ! Regenerate arrays
            end if

            ! Update Chaos (triggers update elements)
            call update_chaos(system, sim%use_baryc_output)

            ! Update geometric Chaos (triggers update geometric elements)
            if (sim%use_geometricfile) call update_chaos_geometric(system)

            ! Output
            if ((checkpoint_is_output(j)) .and. (.not. is_premature_exit)) call generate_output(system)

            ! Percentage output
            if (sim%use_percentage) call percentage(time, sim%final_time)

            ! Check update from TOM
            if (checkpoint_is_tom(j) .and. (tom%index_number <= tom%total_number)) then
                if (tom%use_domega) call spin_asteroid(system%asteroid, &
                                                     & system%asteroid%theta, &
                                                     & system%asteroid%omega + tom%deltaomega(tom%index_number))
                if (tom%use_dmass) call grow_asteroid(system%asteroid, tom%deltamass(tom%index_number))
                tom%index_number = tom%index_number + 1
                if (sim%use_screen .and. (tom%use_dmass .or. tom%use_domega)) then
                    write(*,*) ACHAR(5)
                    write(*,s1r1) "Updated asteroid Omega and mass following TOM data, at time = ", time / unit_time, "[days]"
                    write(*,*) ACHAR(5)
                end if

                ! Apply colissions/escapes and checks
                call check_esc_and_col(system, unit_file)

                ! Update Nactive and y_arr if necessary
                call get_Nactive(system, new_Nactive)
                if (new_Nactive < sim%Nactive) then
                    call update_sim_Nactive(sim, system%Nparticles_active, system%Nmoons_active)  ! Update sim Nactive
                    y_nvalues = get_index(new_Nactive) + 3  ! Update nvalues to use in y
                end if

                call generate_arrays(system, m_arr, R_arr, y_arr)  ! Mandatory bc of new asteroid
            end if

            ! Reset adaptive if needed
            if (.not. sim%use_adaptive) adaptive_timestep = fixed_timestep

            ! Update j; only if not premature
            if (.not. is_premature_exit) j = j + 1
            

        end do main_loop_no_filter
    
    end if

    !! Porcentaje final
    if (sim%use_percentage .and. time /= sim%final_time) call percentage(sim%final_time + uno, sim%final_time)

    !! Cerrar archivo de salida
    if (sim%use_datafile) then
        close (21)
        if (sim%use_screen) then
            write(*,*) ACHAR(10)
            write(*,*) "Output data saved to file: ", trim(sim%datafile)
        end if
    end if

    !! Final chaos y cierre
    if (sim%use_chaos) then
        if (sim%use_chaosfile) then
            !! Chaosfile
            inquire (unit=22, opened=aux_logical)
            if (.not. aux_logical) then
                open (unit=22, file=trim(sim%chaosfile), status='unknown', action='write')
            else
                rewind(22)
            end if
            call write_chaos(initial_system, system, 22)
            close (22)
            !! Mensaje
            if (sim%use_screen) then 
                write(*,*) ACHAR(10)
                write(*,*) "Chaos data saved into: ", trim(sim%chaosfile)
            end if
        end if
        if (sim%use_datascreen) then
            write(*,*) ACHAR(10)        
            write(*,*) "Chaos:"
            call write_chaos(initial_system, system, 6)
        end if
    end if

    !! Cerrar archivo de geométricos
    if (sim%use_geometricfile) then
        close (23)
        if (sim%use_screen) then
            write(*,*) ACHAR(10)
            write(*,*) "Geometric output data saved to file: ", trim(sim%geometricfile)
        end if
    end if

    !! Final geometric chaos y cierre
    if (sim%use_chaos) then
        if (sim%use_geomchaosfile) then
            !! Chaosfile
            inquire (unit=24, opened=aux_logical)
            if (.not. aux_logical) then
                open (unit=24, file=trim(sim%geomchaosfile), status='unknown', action='write')
            else
                rewind(24)
            end if
            call write_to_geomchaos(initial_system, system, 24)
            close (24)
            !! Mensaje
            if (sim%use_screen) then 
                write(*,*) ACHAR(10)
                write(*,*) "Geometric chaos data saved into: ", trim(sim%geomchaosfile)
            end if
            if (sim%use_datascreen) then
                write(*,*) ACHAR(10)        
                write(*,*) "Geometric chaos:"
                call write_geomchaos(initial_system, system, 6)
            end if
        end if
    end if

    !! Cerrar archivos de filtrado
    if (sim%use_filter) then

        !! Cerrar archivo de salida filtrado
        if (sim%use_datafile) then
            close (31)
            if (sim%use_screen) then
                write(*,*) ACHAR(10)
                write(*,*) "Filterd output data saved to file: ", trim(sim%filter_prefix) // trim(sim%datafile)
            end if
        end if

        !! Final chaos filtrados
        if (sim%use_chaos) then
            if (sim%use_chaosfile) then
                !! Chaosfile
                inquire (unit=32, opened=aux_logical)
                if (.not. aux_logical) then
                    open (unit=32, file=trim(sim%chaosfile), status='unknown', action='write')
                else
                    rewind(32)
                end if
                call write_chaos(initial_system, system_filtered, 32)
                close (32)
                !! Mensaje
                if (sim%use_screen) then 
                    write(*,*) ACHAR(10)
                    write(*,*) "Filterd chaos data saved into: ", "fil" // trim(sim%chaosfile)
                end if
            end if
            if (sim%use_datascreen) then
                write(*,*) ACHAR(10)        
                write(*,*) "Chaos filtered:"
                call write_chaos(initial_system, system_filtered, 6)
            end if
        end if

        ! !! Cerrar archivo de geométricos filtrado
        ! if (sim%use_geometricfile) then
        !     close (33)
        !     if (sim%use_screen) then
        !         write(*,*) ACHAR(10)
        !         write(*,*) "Filterd geometric data saved to file: ", trim(sim%filter_prefix) // trim(sim%geometricfile)
        !     end if
        ! end if

        !! Cerrar archivos individuales filtrados
        if (sim%use_multiple_outputs) then
            do i = 0, sim%Ntotal  ! 0 is the asteroid
                close (3100+i)
            end do
            if (sim%use_screen) then
                write(*,*) ACHAR(10)
                write(*,*) "Filterd individual output data saved to files: ", trim(sim%filter_prefix) // trim(sim%multfile) // "_*"
            end if
        end if

    end if

    !! Cerrar archivos individuales
    if (sim%use_multiple_outputs) then
        do i = 0, sim%Ntotal  ! 0 is the asteroid
            close (100+i)
        end do
        if (sim%use_screen) then
            write(*,*) ACHAR(10)
            write(*,*) "Individual output data saved to files: ", trim(sim%multfile) // "_*"
        end if
    end if


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!! LIBERACION MEMORIA!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    if (sim%use_screen) then
        write(*,*) ACHAR(10)
        write(*,*) "Freeing memory..."
        write(*,*) ACHAR(5)
    end if

    !!!!!!!!! PARAMS !!!!!!!!
    call free_parameters_arays()

    !!!!!!!!! TOM !!!!!!!!
    call free_tom(tom)
    
    !!!!!!!!! POINTERS !!!!!!!!
    call nullify_pointers()    

    !!!!!!!!! BODIES !!!!!!!!
    call free_asteroid(asteroid)
    call free_system(system)
    call free_system(initial_system)
    call free_system(tmp_system)

    !!!!!!!!! FILTER !!!!!!!!
    call free_filter(filter)
    call free_system(system_filtered)

    !!!!!!!!! INTEGRATORS !!!!!!!!
    call free_integrator()

    if (sim%use_screen) then 
        write(*,*) ACHAR(5)
        write(*,*) "End of program"
        write(*,*) ACHAR(5)
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
        write(*,*) ACHAR(10)
        write(*,*) "ERROR: Check that (ngx > 2) and (ngy > 2)"
        stop 1
    end if
    if ((sim%map_min_x >= sim%map_max_x) .or. (sim%map_min_y >= sim%map_max_y)) then
        write(*,*) ACHAR(10)
        write(*,*) "ERROR: Check that (xmin < xmax) and (ymin < ymax)"
        stop 1
    end if

    ! Allocate
    allocate(pot(sim%map_grid_size_x, sim%map_grid_size_y))
    allocate(acc(sim%map_grid_size_x, sim%map_grid_size_y, 2))
    
    pot = cero
    acc = cero
    dx_nx = (sim%map_max_x - sim%map_min_x) / (sim%map_grid_size_x - 1)
    dy_ny = (sim%map_max_y - sim%map_min_y) / (sim%map_grid_size_y - 1)

    !$OMP PARALLEL DEFAULT(SHARED) &
    !$OMP PRIVATE(i,j,rb)
    !$OMP DO SCHEDULE (STATIC)
    do i = 1, sim%map_grid_size_x
        do j = 1, sim%map_grid_size_y
            rb(1) = sim%map_min_x + (i - 1) * dx_nx
            rb(2) = sim%map_min_y + (j - 1) * dy_ny
            call get_acc_and_pot_xy(system, rb, acc(i,j,:), pot(i,j))
        end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL
    
    if (sim%use_screen) write(*,*) "Potential calculated. Writing file..."
    open (unit=50, file=trim(sim%mapfile), status='replace', action='write')
    do i = 1, sim%map_grid_size_x
        do j = 1, sim%map_grid_size_y
            write (50,'(5(1PE22.15,1X))') sim%map_min_x + (i - 1) * dx_nx, sim%map_min_y + (j - 1) * dy_ny, pot(i,j), acc(i,j,:)
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

