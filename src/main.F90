program main
    use forces
    use integrators
    
    implicit none

    call init_default_parameters()  ! Inicializamos parámetros

    existe_configfile = .True. ! Existe archivo de configuración
    arguments_number = command_argument_count()
    do i = 1, arguments_number
        call get_command_argument(i, aux_character30)
        if (trim(aux_character30) .eq. "--noconfig") then
            existe_configfile = .False.
            exit
        end if
    end do

    call read_config_file("config.ini", existe_configfile) ! Leemos archivo de configuración

    if (.not. existe_configfile) then ! Usaremos parámetros por defecto
        
        ! Asteroide central
        !! Primary
        mass_primary = 6.3d18 ! Masa del cuerpo 0 [kg]
        radius_primary = 129.d0 ! Radio del cuerpo 0 [km]

        !! ROTACIÓN
        asteroid_rotational_period = 7.004d0/24.d0  ! Periodo de rotación [day]
        !lambda_kep = 0.471d0      ! Cociente omega/wk

        !! Boulders        
        Nboulders = 1 ! Número de boulders

        if (Nboulders > 0) then
        
            !! Alocamos (No tocar)
            call allocate_asteriod_arrays(Nboulders)

            !!!!! COCIENTE DE MASAS
            mu_from_primary(1) = 1.d-1        ! Cociente de masas entre boulder 1 y primary
            ! mu_from_primary(2) = 1.d-6        ! Cociente de masas entre boulder 2 y primary
            ! mu_from_primary(3) = 1.d-4        ! Cociente de masas entre boulder 3 y primary
            ! mu_from_primary(4) = 1.d-7        ! Cociente de masas entre boulder 4 y primary

            !!! ÁNGULOS DE FASE RESPECTO AL ORIGEN [deg]
            theta_from_primary(1) = cero   ! Ángulo de fase del boulder 1[deg]
            ! theta_from_primary(2) = 90     ! Ángulo de fase del boulder 2[deg]
            ! theta_from_primary(3) = 180    ! Ángulo de fase del boulder 3[deg]
            ! theta_from_primary(4) = 270    ! Ángulo de fase del boulder 4[deg]

            !!! RADIOS
            radius_ast_arr(1) = 2.5d0   ! Radio del boulder 1 [km]
            ! radius_ast_arr(2) = 1.5d0   ! Radio del boulder 2 [km]
            ! radius_ast_arr(3) = 3.5d0   ! Radio del boulder 3 [km]
            ! radius_ast_arr(4) = 2.5d0   ! Radio del boulder 3 [km]
        
        end if

        !!!! Merges (colliding particles into asteroid)
        use_merge = .True. ! Merge particles into asteroid

        !!!! Torque
        use_torque = .False. ! Cant be used with explicit method, or exponential omega damping
        
        !!!! Stokes
        use_stokes = .False.
        stokes_a_damping_time = infinity     ! [day]
        stokes_e_damping_time = stokes_a_damping_time / 1.d2 ! [day]
        stokes_charac_time = cero                   ! [day] Tiempo que actua stokes

        !!!! Naive-Stokes (drag)
        use_naive_stokes = .False.
        drag_coefficient = 0.0d0

        !!!! Geo-Potential (J2)
        use_J2 = .False.
        J2_coefficient = cero

        !!! Parámetros corrida
        !!!! Tiempos
        initial_time = cero      ! Initial time [day]
        final_time = 2.d3      ! Final time [day]
        min_timestep = cero      ! Min timestep [day] ! Almost unused
        use_logspaced_output = .False.   ! LogSpaced outputs
        output_number = 10000      ! Number of outputs (if logsp=.True. or dt_out=0)
        output_timestep = cero      ! Output timestep [day] (if logsp = .False.)
        !!!! Error
        learning_rate = 0.85d0    ! [For adaptive step integrators] Learning rate
        error_digits = 12        ! [For adaptive step integrators] Digits for relative error
        !!!! Colision y escape
        min_distance = -1.5d0    ! Min distance before impact [km] ! 0 => R0 + max(Rboul)
        max_distance = -1.d2     ! Max distance before escape [km] ! -x => R0 * x
        
        !!!! Particle elements
        use_single_particle = .False.          ! Single particle? (THIS ONE BELOW)
        single_part_mass = 6.3e14              ! Mass of the particle [kg]
        single_part_elem_a = cero              ! Element a of the particle [km]
        single_part_elem_e = cero !0.1d0       ! ecc
        single_part_elem_M = cero              ! Mean anomaly (deg)
        single_part_elem_w = cero              ! Pericenter argument (deg)
        single_part_MMR = 5.2d0 !1.001d0 !3.3d0 ! resonancia nominal correspondiente

        !!! Output: "" or "no", if not used
        datafile = ""
        chaosfile = ""
        mapfile = ""
        multfile = ""
        !!!!! Screeen
        use_datascreen = .True. ! Print data in screen
        use_screen = .True. ! Print info in screen
        !!! Input: "" or "no", if not used
        tomfile = ""
        particlesfile = ""

        !!! Parallel
        use_parallel = .False.
        requested_threads = 1 ! Number of threads to use !! -1 => all available
    end if

    call load_command_line_arguments()
    
    call set_derived_parameters() ! Inicializamos parámetros derivados


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!! No tocar de aquí a abajo !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    !!!!!!!!    OUTPUT     !!!!!!!!!!

    if ((trim(datafile) /= "") .and. (trim(datafile) /= "no")) then
        use_datafile = .True.
    else 
        use_datafile = .False.
    end if
    if ((trim(chaosfile) /= "") .and. (trim(chaosfile) /= "no")) then
        use_chaosfile = .True.
        use_chaos = .True.
    else 
        use_chaosfile = .False.
    end if
    if ((trim(multfile) /= "") .and. (trim(multfile) /= "no")) then
        use_multiple_outputs = .True.
    else 
        use_multiple_outputs = .False.
    end if
    if ((trim(mapfile) /= "") .and. (trim(mapfile) /= "no")) then
        use_potential_map = .True.
    else 
        use_potential_map = .False.
    end if
    if ((trim(tomfile) /= "") .and. (trim(tomfile) /= "no")) then
        use_tomfile = .True.
    else 
        use_tomfile = .False.
    end if
    if ((trim(particlesfile) /= "") .and. (trim(particlesfile) /= "no")) then
        use_particlesfile = .True.
    else 
        use_particlesfile = .False.
    end if

    if (.not. any((/use_datascreen, use_datafile, use_chaosfile, use_potential_map, use_multiple_outputs/))) then ! No tiene sentido hacer nada
        if (.not. use_screen) then
            write (*,*) ACHAR(10)
            write (*,*)  "EXITING: No se guardará ni imprimirá ninguna salida."
            stop 1
        else 
            write (*,*)  "WARNING: No se integrará. Solo se imprimirá en pantalla información."
        end if
    end if

    if (.not. any((/use_single_particle, use_particlesfile/))) then ! No tiene sentido hacer nada
        write (*,*) ACHAR(10)
        write (*,*)  "EXITING: No se integrará ninguna partícula."
        stop 1
    end if
    
    ! Mensaje
    if (use_screen .and. .not. existe_configfile) then
        write (*,*) "WARNING: No se pudo leer el archivo de configuración: config.ini"
        write (*,*) "         Se utilizan los parámetros dados en el código."
    end if


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  Parallel  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (use_screen) then
        write (*,*) ACHAR(5)
        if (use_parallel) then
            write (*,*) "---------- Parallel integration----------"
            write (*,*) "  Available threads:", available_threads
            write (*,*) "  Requested threads:", requested_threads
            write (*,*) "  Used threads     :", my_threads
        else
            write (*,*) "---------- Serial integration ----------"
        end if
        write (*,*) ACHAR(5)
    end if
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!! COMIENZO DE CÁLCULOS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (use_screen) then
        write (*,*) ACHAR(5)
        write (*,*) "Comenzando simulación: ", simulation_number
        write (*,*) ACHAR(5)
        write (*,*) "---------- Parámetros iniciales ----------"
        write (*,*) ACHAR(5)
    end if
    
    !!!!!!!!!!!!!!!!!!!!!!!!!! Asteroide y Boulders !!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Radios (unidades)
    radius_primary = radius_primary * unit_dist ! Unidades según G
    radius_ast_arr = radius_ast_arr * unit_dist ! Unidades según G  
    radius_ast_arr(0) = radius_primary ! Radio del cuerpo 0  



    ! Masas
    mass_primary = mass_primary * unit_mass ! Unidades según G
    mass_ast_arr(0) = mass_primary ! Masa del cuerpo 0
    do i = 1, Nboulders
        mass_ast_arr(i) = mu_from_primary(i) * mass_primary ! Masa del boulder i
    end do
    asteroid_mass = sum(mass_ast_arr)  ! Masa del sistema
    mu_ast_arr = mass_ast_arr / asteroid_mass ! Mues respecto a mcm
    Gmass_ast_arr = G * mass_ast_arr
    Gasteroid_mass = G * asteroid_mass
    


    ! Rotaciones
    theta_from_primary = theta_from_primary * radian ! Ángulos de fase de los boulders [rad]
    omega_kep = sqrt(Gasteroid_mass / radius_primary**3) / unit_time ! Movimiento medio (kepleriano) que tendrían los boulder
    if (lambda_kep > tini) then    
        asteroid_omega = omega_kep * lambda_kep  ! Velocidad angular del cuerpo 0 [rad/day]
        asteroid_rotational_period = twopi / asteroid_omega ! Periodo de rotación del cuerpo 0 [day]
    else
        asteroid_omega = (twopi / asteroid_rotational_period) / unit_time ! Velocidad angular del cuerpo 0 [rad/day]
        lambda_kep = asteroid_omega / omega_kep ! Cociente de velocidades
    end if
    asteroid_omega2 = asteroid_omega * asteroid_omega
    asteroid_a_corot = (Gasteroid_mass / asteroid_omega2)**(1/3.)

    !! Mensaje
    if (use_screen) then
        write (*,*) "Rotación:"
        write (*,*) "  a_corot   :", asteroid_a_corot / unit_dist, "[km]"
        write (*,*) "  omega     :", asteroid_omega * unit_time, "[rad/day]"
        write (*,*) "  omega_kep :", omega_kep * unit_time, "[rad/day]"
        write (*,*) "  lambda_kep:", lambda_kep
        write (*,*) "  Period    :", asteroid_rotational_period * 24 * unit_time, "[hs]"
        write (*,*) ACHAR(5)
        write (*,*) "Masa cuerpo central    :", mass_primary / unit_mass, "[kg]"
        write (*,*) "Radio de cuerpo central:", radius_primary / unit_dist, "[km]"
        write (*,*) ACHAR(5)
    end if



    ! Coordenadas

    !! Astrocentricas
    do i = 1, Nboulders
        pos_from_primary(i,1) = radius_primary * cos(theta_from_primary(i)) ! Posición x del boulder i
        pos_from_primary(i,2) = radius_primary * sin(theta_from_primary(i)) ! Posición x del boulder i
        vel_from_primary(i,1) = - asteroid_omega * pos_from_primary(i,2) ! Velocidad x del boulder i
        vel_from_primary(i,2) =   asteroid_omega * pos_from_primary(i,1) ! Velocidad x del boulder i
        acc_from_primary(i,1) = - asteroid_omega2 * pos_from_primary(i,1) ! Aceleración x del boulder i
        acc_from_primary(i,2) = - asteroid_omega2 * pos_from_primary(i,2) ! Aceleración x del boulder i
    end do
    !!! Mensaje
    if (use_screen) then
        write (*,*) "Coordenadas:"
        write (*,*) "   m0-centricas: [x, y, vx, vy, ax, ay, mu, theta]"
        do i = 1, Nboulders
            write (*,*) i, &
            & pos_from_primary(i,:) / unit_dist, &
            & pos_from_primary(i,:) / unit_vel, &
            & acc_from_primary(i,:) / unit_acc, &
            & mu_from_primary(i), &
            & theta_from_primary(i) / radian
        end do
        write (*,*) ACHAR(5)
    end if


    !! Centro de masas, desde primary
    pos_ast_from_primary = cero
    vel_ast_from_primary = cero
    acc_ast_from_primary = cero
    theta_ast_from_primary = cero
    do i = 1, Nboulders
        pos_ast_from_primary = pos_ast_from_primary + mu_ast_arr(i) * pos_from_primary(i,:) ! Posición del asteroide
        vel_ast_from_primary = vel_ast_from_primary + mu_ast_arr(i) * vel_from_primary(i,:) ! Velocidad del asteroide
        acc_ast_from_primary = acc_ast_from_primary + mu_ast_arr(i) * acc_from_primary(i,:) ! Aceleración del asteroide
    end do
    if (Nboulders > 0) theta_ast_from_primary = atan2(pos_ast_from_primary(2), pos_ast_from_primary(1)) ! Ángulo del asteroide
    !!! Mensaje
    if (use_screen) then
        write (*,*) "   Centro de masas desde m0: [x, y, vx, vy, ax, ay, mass]"
        write (*,*) "         CM", pos_ast_from_primary / unit_dist, &
                                & vel_ast_from_primary / unit_vel, &
                                & acc_ast_from_primary / unit_acc, &
                                & asteroid_mass / unit_mass
        write (*,*) ACHAR(5)
    end if


    !! Baricentricas
    !!! Inicializamos asteroide
    pos_ast_arr(0,:) = - pos_ast_from_primary
    vel_ast_arr(0,:) = - vel_ast_from_primary
    acc_ast_arr(0,:) = - acc_ast_from_primary
    theta_ast_arr(0) = mod(theta_ast_from_primary + pi, twopi)
    dist_ast_arr(0) = sqrt(pos_ast_arr(0,1)*pos_ast_arr(0,1) + pos_ast_arr(0,2)*pos_ast_arr(0,2))
    do i = 1, Nboulders
        pos_ast_arr(i,:) = pos_from_primary(i,:) - pos_ast_from_primary
        vel_ast_arr(i,:) = vel_from_primary(i,:) - vel_ast_from_primary
        acc_ast_arr(i,:) = acc_from_primary(i,:) - acc_ast_from_primary
        theta_ast_arr(i) = atan2(pos_ast_arr(i,2), pos_ast_arr(i,1))
        dist_ast_arr(i) = sqrt(pos_ast_arr(i,1)*pos_ast_arr(i,1) + pos_ast_arr(i,2)*pos_ast_arr(i,2))
    end do
    !!! Mensaje
    if (use_screen) then
        write (*,*) "   Baricentricas: [x, y, vx, vy, ax, ay, mass, distance, theta]"
        do i = 0, Nboulders
            write (*,*) i, &
                     & pos_ast_arr(i,:) / unit_dist, &
                     & vel_ast_arr(i,:) / unit_vel, &
                     & acc_ast_arr(i,:) / unit_acc, &
                     & mass_ast_arr(i) / unit_mass, &
                     & dist_ast_arr(i) / unit_dist, &
                     & theta_ast_arr(i) / radian
        end do
        write (*,*) ACHAR(5)
    end if

    

    ! Angular momentum and Inertia
    asteroid_angmom = cero
    asteroid_inertia = cero
    do i = 0, Nboulders
        asteroid_angmom = asteroid_angmom + mass_ast_arr(i) * cross2D(pos_ast_arr(i,:), vel_ast_arr(i,:)) ! Traslacional
        ! asteroid_angmom = asteroid_angmom + 0.4d0 * mass_ast_arr(i) * radius_ast_arr(i)**2 * asteroid_omega ! Rotacional
        asteroid_inertia = asteroid_inertia + mass_ast_arr(i) * dist_ast_arr(i)**2 ! Inercia
    end do
    !! Mensaje
    if (use_screen) then
        write (*,*) "Momento angular total:", asteroid_angmom / unit_mass / (unit_dist**2) * unit_time, "[kg km^2 / day]"
        write (*,*) "Momento inercia asteroide", asteroid_inertia / (unit_mass * unit_dist**2), "[kg km^2]"
        write (*,*) "Check:"
        write (*,*) "    L / (I * omega):", asteroid_angmom / asteroid_inertia / asteroid_omega
        write (*,*) "    Relative error :", abs(asteroid_omega - (asteroid_angmom / asteroid_inertia)) / asteroid_omega
        write (*,*) ACHAR(5)
    end if



    ! Asteroide, desde centro del sistema
    asteroid_pos = cero
    asteroid_vel = cero
    asteroid_acc = cero
    asteroid_theta = cero
    asteroid_theta_correction = cero ! Corrección angular
    do i = 0, Nboulders
        asteroid_pos = asteroid_pos + mass_ast_arr(i) * pos_ast_arr(i,:)
        asteroid_vel = asteroid_vel + mass_ast_arr(i) * vel_ast_arr(i,:)
    end do
    asteroid_pos = asteroid_pos / asteroid_mass ! rcm = sum_i m_i * r_i / M
    asteroid_vel = asteroid_vel / asteroid_mass ! vcm = sum_i m_i * v_i / M
    !!!! Mensaje
    if (use_screen) then
        write (*,*) "Centro de masas del asteroide:"
        write (*,*) "  [x, y]  : ", asteroid_pos / unit_dist, "[km]"
        write (*,*) "  [vx, vy]: ", asteroid_vel / unit_vel, "[km/day]"
        write (*,*) ACHAR(5)
    end if
    !! Check y Centrado
    if (hard_center) then
        if (any(abs(asteroid_pos) > cero)) then
            if (use_screen) write (*,*) "  - Se centrarán las posiciones para que el asteroide esté en el origen."
            do i = 0, Nboulders
                pos_ast_arr(i,:) = pos_ast_arr(i,:) - asteroid_pos
            end do        
            asteroid_pos = cero
        end if
        if (any(abs(asteroid_vel) > cero)) then
            if (use_screen) write (*,*) "  - Se centrarán las velocidades para que el asteroide esté en el origen."
            do i = 0, Nboulders
                vel_ast_arr(i,:) = vel_ast_arr(i,:) - asteroid_vel
            end do
            asteroid_vel = cero
        end if
    end if

    

    ! EFECTOS EXTERNOS

    if (use_screen) then
        write (*,*) ACHAR(5)
        write (*,*) ("---------- Efectos ex/internos ----------")
        write (*,*) ACHAR(5)
    end if

    !! Variación en masa o en omega (funciones del tiempo)
    mass_exp_damping_time = cero !!! Deprecado
    
    if ((abs(omega_linear_damping_time) > tini) .and. (omega_linear_damping_time < infinity)) then
        if (use_explicit_method) then
            write (*,*) ACHAR(10)
            write (*,*) "ERROR: No se puede usar el método explícito con tau_o finito."
            write (*,*) "tau_o [Prot]:", omega_linear_damping_time * unit_time / asteroid_rotational_period
            stop 1
        end if
        omega_linear_damping_time = omega_linear_damping_time * unit_time
        omega_linear_damping_slope = - asteroid_omega / (omega_linear_damping_time - initial_time)
        write (*,*) "Explicit omega linear damping"
        write (*,*) "    tau_o :", omega_linear_damping_time / asteroid_rotational_period, "[Prot]"
    else
        omega_linear_damping_time = infinity
        omega_linear_damping_slope = cero
    end if
    if ((abs(omega_exp_damping_time) > tini) .and. (omega_exp_damping_time < infinity)) then
        if (use_explicit_method) then
            write (*,*) ACHAR(10)
            write (*,*) "ERROR: No se puede usar el método explícito con tau_o finito."
            write (*,*) "tau_o [Prot]:", omega_exp_damping_time * unit_time / asteroid_rotational_period
            stop 1
        end if
        omega_exp_damping_time = omega_exp_damping_time * unit_time
        write (*,*) "Explicit omega exponential damping"
        write (*,*) "    tau_o :", omega_exp_damping_time / asteroid_rotational_period, "[Prot]"
    else
        omega_exp_damping_time = infinity
    end if

    if ((abs(mass_exp_damping_time) > tini) .and. (mass_exp_damping_time < infinity)) then
        if (use_explicit_method) then
            write (*,*) ACHAR(10)
            write (*,*) "ERROR: No se puede usar el método explícito con tau_m finito."
            write (*,*) "tau_m [Prot]:", mass_exp_damping_time * unit_time / asteroid_rotational_period
            stop 1
        end if
        mass_exp_damping_time = mass_exp_damping_time * unit_time
        write (*,*) "Explicit mass exponential damping"
        write (*,*) "    tau_m :", mass_exp_damping_time / asteroid_rotational_period, "[Prot]"
    else
        mass_exp_damping_time = infinity
    end if

    !!! Check not both omega dampings
    if (omega_exp_damping_time < infinity .and. omega_linear_damping_time < infinity) then
        write (*,*) ACHAR(10)
        write (*,*) "ERROR: No se puede usar ambos decaimientos de omega (linear and exp) finitos."
        stop 1
    end if

    !!! Check que no esté torque también
    if (use_torque .and. ((omega_exp_damping_time < infinity) .or. (omega_linear_damping_time < infinity))) then
        write (*,*) ACHAR(10)
        write (*,*) "ERROR: No se puede usar torque y tau_o finito."
        stop 1
    end if
    !!!

    !! Stokes
    if (use_stokes) then
        call set_stokes_C_and_alpha(stokes_a_damping_time, stokes_e_damping_time, stokes_C, stokes_alpha)
        stokes_a_damping_time = stokes_a_damping_time * unit_time
        stokes_e_damping_time = stokes_e_damping_time * unit_time
        stokes_charac_time = stokes_charac_time * unit_time
        !!! Mensaje
        if (use_screen) then
            write (*,*) "Stokes"
            write (*,*) "    t_stokes: ", stokes_charac_time / asteroid_rotational_period, "[Prot]"
            write (*,*) "    C       : ", stokes_C
            write (*,*) "    alpha   : ", stokes_alpha
        end if
    else ! Just to be sure
        stokes_a_damping_time = infinity
        stokes_e_damping_time = infinity
        stokes_charac_time = cero
    end if


    !! Mensaje !
    if (use_screen) then
        if (use_naive_stokes) then !!!! Naive-Stokes
            write (*,*) "Naive-Stokes (drag radial)"
            write (*,*) "    Eta :", drag_coefficient
        end if
        if (use_J2) then
            write (*,*) "Geo-Potential (J2)"
            write (*,*) "    J2 :", J2_coefficient !!!! Geo-Potential (J2)
        end if
        if (use_torque) write (*,*) "Torque from particles to asteroid ACTIVATED"
        
        if (.not. any((/use_stokes, use_naive_stokes, use_J2, use_torque, &
            & omega_linear_damping_time < infinity, &
            & omega_exp_damping_time < infinity, &
            & mass_exp_damping_time < infinity/))) then
            write (*,*) "No se aplicarán efectos ex/internos."
        end if
    end if


    !!!!!!!!!!!!!!!!!!!!! MAPA Potencial y Aceleraciones !!!!!!!!!!!!!!!!!!!!!!
    if (use_potential_map) then
        if (use_screen) then
            write (*,*) ACHAR(5)
            write (*,*) "---------- MAPA potencial (y aceleracion) ----------"
            write (*,*) ACHAR(5)
            write (*,*) "Creando mapa de potencial..."
        end if
        call create_map(Nboulders, &
                      & map_grid_size_x, map_grid_size_y, &
                      & map_min_x, map_max_x, &
                      & map_min_y, map_max_y, &
                      & mass_ast_arr,pos_ast_arr,vel_ast_arr, &
                      & mapfile)
        if (use_screen) then
            write (*,*) "Guardado en el archivo: ", trim(mapfile)
            write (*,*) ACHAR(5)
        end if
    end if
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    
    if (use_screen) then
        write (*,*) ACHAR(5)
        write (*,*) "---------- Preparando partículas ----------"
        write (*,*) ACHAR(5)
    end if



    ! Particle(s)
    Nparticles = 0
    !! Partiles?
    if (use_single_particle) Nparticles = Nparticles + 1
    if (use_particlesfile) then
        if (use_screen) write (*,*) "Leyendo partículas desde archivo: ", trim(particlesfile)
        !!!! Formato: [mass, a, e, M, w, (MMR)]
        call read_columns_file(particlesfile, aux_particles_arr)
        Narticles_in_partfile = size(aux_particles_arr, 1)
        if (use_screen) then
            write (*,*) "  Se han leído ", Narticles_in_partfile, " filas."
            write (*,*) "  Se han leído ", size(aux_particles_arr, 2), " columnas."
        end if
        ! [Remember to free memory later...]
    else 
        Narticles_in_partfile = 0
    end if

    !!! Sumamos
    Nparticles = Nparticles + Narticles_in_partfile

    !!!! Redefine Chaos (needed to allocate arrays)
    use_chaos = (use_chaosfile .or. (use_datascreen .and. (Nparticles .eq. 1))) .or. use_chaos

    !! Allocate
    if (Nparticles > 0) then
        call allocate_particles(Nparticles)
    else 
        write (*,*) ACHAR(10)
        write (*,*) "ERROR: No hay partículas para integrar."
        stop 1
    end if

    !! Redefine OMP threads if necessary. Lower equal to particle number
    if ((use_parallel) .and. (Nparticles < my_threads)) then
        my_threads = Nparticles
        if (use_screen) then
            write (*,*) "WARNING: Se reducen los threads a ", my_threads
            write (*,*) ACHAR(5)
        end if
        !$ call omp_set_num_threads(my_threads)
    end if

    
    !! Create Index
    do i = 1, Nparticles
        particles_index(i) = i
        sorted_particles_index(i) = i
    end do
    if (use_particlesfile) then
        particles_mass(:Narticles_in_partfile) = aux_particles_arr(:,1) * unit_mass   ! Masa
        particles_elem(:Narticles_in_partfile,1) = aux_particles_arr(:,2) * unit_dist ! a
        particles_elem(:Narticles_in_partfile,2) = aux_particles_arr(:,3)             ! e
        particles_elem(:Narticles_in_partfile,3) = aux_particles_arr(:,4) * radian    ! M
        particles_elem(:Narticles_in_partfile,4) = aux_particles_arr(:,5) * radian    ! w
        if (size(aux_particles_arr, 2) > 5) then
            !$OMP PARALLEL IF(my_threads > 1 .AND. Narticles_in_partfile > 10) DEFAULT(SHARED) &
            !$OMP PRIVATE(i)
            !$OMP DO SCHEDULE (STATIC)
            do i = 1, Narticles_in_partfile
                if (aux_particles_arr(i,6) > tini) then
                    particles_MMR(i) = aux_particles_arr(i,6) ! MMR
                    particles_elem(i,1) = particles_MMR(i)**(2/3.) * asteroid_a_corot ! a
                else
                    particles_MMR(i) = (particles_elem(i,1) / asteroid_a_corot)**(1.5) ! MMR
                end if
            end do
            !$OMP END DO
            !$OMP END PARALLEL
        else 
            particles_MMR(:Narticles_in_partfile) = (particles_elem(:Narticles_in_partfile,1) / asteroid_a_corot)**(1.5)
        end if
        !! Deallocatamos
        deallocate(aux_particles_arr)
    end if
    
    if (use_single_particle) then
        !! Masa
        particles_mass(Nparticles) = single_part_mass * unit_mass     ! Masa
        !! Elementos orbitales
        particles_elem(Nparticles,1) = single_part_elem_a * unit_dist ! a
        particles_elem(Nparticles,2) = single_part_elem_e             ! e
        particles_elem(Nparticles,3) = single_part_elem_M * radian    ! M
        particles_elem(Nparticles,4) = single_part_elem_w * radian    ! w
        if (single_part_MMR > tini) then
            particles_elem(Nparticles,1) = single_part_MMR**(2/3.) * asteroid_a_corot
            particles_MMR(Nparticles) = single_part_MMR
        else
            particles_MMR(Nparticles) = (particles_elem(Nparticles,1) / asteroid_a_corot)**(1.5)
        end if
    end if


    !! Crear vector de coordenadas
    aux_real_arr6 = cero
    !$OMP PARALLEL IF(my_threads > 1 .AND. Nparticles > 10) DEFAULT(NONE) &
    !$OMP PRIVATE(i,aux_real_arr6) &
    !$OMP SHARED(Nparticles,asteroid_mass,particles_mass,particles_elem,particles_coord)
    !$OMP DO SCHEDULE (STATIC)
    do i = 1, Nparticles
        !!! Coordenadas y vectores
        call coord(asteroid_mass + particles_mass(i), & ! Masa
                  & particles_elem(i,1), particles_elem(i,2), cero, & ! a, e, i
                  & particles_elem(i,3), particles_elem(i,4), cero, & ! M, w, Omega
                  & aux_real_arr6)
        particles_coord(i,:) = (/aux_real_arr6(1:2), aux_real_arr6(4:5)/) ! [x, y, vx, vy]
    end do
    !$OMP END DO
    !$OMP END PARALLEL
    
    !! Crear vector distancias
    particles_dist = sqrt(sum(particles_coord(:,1:2) * particles_coord(:,1:2), 2)) ! sqrt(x^2 + y^2)

    !! Chaos and Outcome
    particles_outcome = 0
    if (use_chaos) then
        particles_max_a = cero
        particles_min_a = infinity
        particles_max_e = cero
        particles_min_e = infinity
        particles_times = cero
    end if
    particles_hexit = 0


    !!! Redefine logicals
    if (all(particles_mass .eq. cero)) use_merge = .False. ! No merge if no mass
    use_single_particle = Nparticles .eq. 1
    if (use_single_particle) use_elements = .True.
    
    !! Mensaje
    if (use_screen) then
        if (use_single_particle) then
            write (*,*) "Partícula simple:"
            write (*,*) "    Masa:", particles_mass(Nparticles) / unit_mass, "[kg]"
            write (*,*) "    Elementos orbitales:"
            write (*,*) "        a   :", particles_elem(Nparticles,1) / unit_dist, "[km]"
            write (*,*) "        e   :", particles_elem(Nparticles,2)
            write (*,*) "        M   :", particles_elem(Nparticles,3) / radian, "[deg]"
            write (*,*) "        w   :", particles_elem(Nparticles,4) / radian, "[deg]"
            write (*,*) "        MMR :", particles_MMR(Nparticles)
            write (*,*) "    Coordenadas:"
            write (*,*) "        x   :", particles_coord(Nparticles,1) / unit_dist, "[km]"
            write (*,*) "        y   :", particles_coord(Nparticles,2) / unit_dist, "[km]"
            write (*,*) "        vx  :", particles_coord(Nparticles,3) / unit_vel, "[km/day]"
            write (*,*) "        vy  :", particles_coord(Nparticles,4) / unit_vel, "[km/day]"
            write (*,*) "        dist:", particles_dist(Nparticles) / unit_dist, "[km]"
        else
            write (*,*) "Cantidad total de partículas:", Nparticles
        end if
        write (*,*) ACHAR(5)
    end if
    

    !! Definimos Tipo 1 o Tipo 2, dependiendo si tienen masa
    Nparticles_type1 = count(particles_mass > cero)
    Nparticles_type2 = count(particles_mass <= cero)

    ! ! ! Reordenamos las particulas, dependiendo si tienen masa o no [NO NECESSARY]
    ! ! if ((.not. use_single_particle) .and. (Nparticles_type1 > 1))  then
    ! !     call argsort(particles_mass, sorted_particles_index)
    ! !     do i = 1, Nparticles
    ! !         particles_mass(i) = particles_mass(sorted_particles_index(i))
    ! !         particles_elem(i,:) = particles_elem(sorted_particles_index(i),:)
    ! !         particles_coord(i,:) = particles_coord(sorted_particles_index(i),:)
    ! !         particles_dist(i) = particles_dist(sorted_particles_index(i))
    ! !         particles_MMR(i) = particles_MMR(sorted_particles_index(i))
    ! !     end do
    ! ! else 
    ! !     sorted_particles_index = particles_index
    ! ! end if



    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!! FIN DE CÁLCULOS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    Ntotal = 1 + Nboulders + Nparticles ! Número total de cuerpos a integra
    Nactive = Nparticles ! Número de cuerpos activos (no colisionados o escapados)


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!! INTEGRACIÓN !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Mensaje
    if (use_screen) then
        write (*,*) ACHAR(5)
        write (*,*) "---------- Preparando integración ----------"
        write (*,*) ACHAR(5)
    end if

    
    ! Escape/Colisión
    if (min_distance < cero) then
        min_distance = radius_primary * abs(min_distance)
    else if (min_distance < tini) then
        min_distance = radius_primary + maxval(radius_ast_arr(1:))
    else
        min_distance = min_distance * unit_dist
    end if
    if (max_distance < cero) then
        max_distance = radius_primary * abs(max_distance)
    else 
        max_distance = max_distance * unit_dist
    end if
    if (max_distance <= min_distance) then
        write (*,*) ACHAR(10)
        write (*,*) "ERROR: rmax <= rmin"
        stop 1
    end if
    if (use_screen) then
        write (*,*) "Condición escape/colisión"
        write (*,*) "    rmin : ", min_distance / unit_dist, "[km] =", min_distance / radius_primary, "[R0]"
        write (*,*) "    rmax : ", max_distance / unit_dist, "[km] =", max_distance / radius_primary, "[R0]"
        if (use_merge) then
            write (*,*) "  Las colisiones se resolverán como mergers al asteroide."
        else
            write (*,*) "  Las colisiones NO se resolverán."
        end if
        write (*,*) ACHAR(5)
    end if

    !!!!!!!! TIEMPOS !!!!!!!

    ! Tiempos de integración
    initial_time = initial_time * unit_time
    if (final_time < cero) then
        final_time = abs(final_time) * asteroid_rotational_period
    else
        final_time = final_time * unit_time
    end if
    if (initial_time >= final_time) then
        write (*,*) ACHAR(10)
        write (*,*) "ERROR: t0 >= tf"
        stop 1
    end if

    !! Timesteps
    output_timestep = output_timestep * unit_time
    min_timestep = min_timestep * unit_time

    !! Output times
    call set_output_times(initial_time, final_time, output_number, output_timestep, use_logspaced_output, output_times)

    !! TOMFILE
    if (use_tomfile) then
        !! En este caso, leeremos los tiempos desde un archivo
        if (use_screen) write (*,*) "Leyendo tiempos desde el archivo: ", trim(tomfile)
        call read_tomfile(initial_time / unit_time, final_time / unit_time, tom_times, tom_deltaomega, tom_deltamass, tomfile) ! Read LOOP checkpoints
        tom_total_number = size(tom_times, 1)
        !!! Unidades
        tom_times = tom_times * unit_time
        if (allocated(tom_deltaomega)) tom_deltaomega = tom_deltaomega / unit_time
        if (allocated(tom_deltamass)) tom_deltamass = tom_deltamass * unit_mass
        !!! Condicion inicial (y final)
        tom_times(1) = initial_time
        tom_times(tom_total_number) = final_time
        if (allocated(tom_deltaomega)) then
            tom_deltaomega(1) = cero
            tom_deltaomega(tom_total_number) = cero
        end if
        if (allocated(tom_deltamass)) then
            tom_deltamass(1) = cero
            tom_deltamass(tom_total_number) = cero
        end if
        !!! Mensaje !
        if (use_screen) then
            if (allocated(tom_deltamass)) then
                write (*,*) "  Se han leído 3 columnas: t, omega, m"
            else if (allocated(tom_deltaomega)) then
                write (*,*) "  Se han leído 2 columnas: t, omega"
            else
                write (*,*) "  Se ha leído 1 columna: t"
            end if
        end if
        !! Ahora debemos combinar los tiempos de TOM con los tiempos de Output, y crear un nuevo vector de tiempos Checkpoints
        call merge_sort_and_unique(tom_times, output_times, &
                                   & checkpoint_is_tom, checkpoint_is_output, &
                                   & checkpoint_times, checkpoint_number)
    else
        !! En este caso, los tiempos de check son los mismos que los de LOOP
        checkpoint_number = output_number
        allocate(checkpoint_times(checkpoint_number))
        allocate(checkpoint_is_output(checkpoint_number))
        allocate(checkpoint_is_tom(checkpoint_number))
        checkpoint_times = output_times
        checkpoint_is_output = .True.
        checkpoint_is_tom = .False.
        tom_total_number = 0
    end if

    ! Variable temporal (para integración)
    time = initial_time ! Tiempo actual
    timestep = output_times(1) - initial_time ! Paso de tiempo inicial
    min_timestep = max(min(min_timestep, output_timestep), tini) ! Paso de tiempo mínimo
    adaptive_timestep = min_timestep ! Paso de tiempo adaptativo inicial

    !! Mensaje
    if (use_screen) then
        write (*,*) "Tiempos:"
        write (*,*) "    t0    : ", initial_time / asteroid_rotational_period, "[Prot]", &
            & " = ", initial_time / unit_time, "[day]"
        write (*,*) "    tf    : ", final_time / asteroid_rotational_period, "[Prot]", &
            & " = ", final_time / unit_time, "[day]"
        write (*,*) "    dt_out: ", output_timestep / asteroid_rotational_period, "[Prot]", &
            & " = ", output_timestep / unit_time, "[day]"
        write (*,*) "    dt_min: ", min_timestep / asteroid_rotational_period, "[Prot]", &
            & " = ", min_timestep / unit_time, "[day]"
        write (*,*) "    n_out : ", output_number
        write (*,*) ACHAR(5)
    end if


    
    !!!!!!!! VECTOR A INTEGRAR !!!!!!!
    !!!! Nos pararemos en el centro de masas del sistema (asteroide).
    !! Calculate the center of mass of the asteroid
    !!!! Inicializamos vector
    if (use_version_1) then
        !!!!! Version 1: [x0, y0, vx0, vy0, x1, y1, vx1, vy1, ...]
        allocate(parameters_arr(4 * Ntotal))
        allocate(parameters_arr_new(4 * Ntotal))
        do i = 0, Nboulders
            parameters_arr(1 + 4 * i) = pos_ast_arr(i,1)
            parameters_arr(2 + 4 * i) = pos_ast_arr(i,2)
            parameters_arr(3 + 4 * i) = vel_ast_arr(i,1)
            parameters_arr(4 + 4 * i) = vel_ast_arr(i,2)
        end do
        !$OMP PARALLEL IF(my_threads > 1 .AND. Nparticles > 10) DEFAULT(SHARED) &
        !$OMP PRIVATE(i)
        !$OMP DO SCHEDULE (STATIC)
        do i = 1, Nparticles
            parameters_arr(1 + 4 * (i + Nboulders)) = particles_coord(i,1) ! xP
            parameters_arr(2 + 4 * (i + Nboulders)) = particles_coord(i,2) ! yP
            parameters_arr(3 + 4 * (i + Nboulders)) = particles_coord(i,3) ! vPx
            parameters_arr(4 + 4 * (i + Nboulders)) = particles_coord(i,4) ! vPy
        end do
        !$OMP END DO
        !$OMP END PARALLEL
        first_particle = 1 + 4 * (Nboulders + 1)
    else 
        !!!!! Version 2: [theta, omega, xA, yA, vxA, vyA, Part...]
        allocate(parameters_arr(6 + 4 * Nparticles))
        allocate(parameters_arr_new(6 + 4 * Nparticles))
        parameters_arr(1) = asteroid_theta
        parameters_arr(2) = asteroid_omega
        parameters_arr(3:4) = asteroid_pos
        parameters_arr(5:6) = asteroid_vel
        !$OMP PARALLEL IF(my_threads > 1 .AND. Nparticles > 10) DEFAULT(SHARED) &
        !$OMP PRIVATE(i)
        !$OMP DO SCHEDULE (STATIC)
        do i = 1, Nparticles
            parameters_arr(3 + 4 * i) = particles_coord(i,1) ! xP
            parameters_arr(4 + 4 * i) = particles_coord(i,2) ! yP
            parameters_arr(5 + 4 * i) = particles_coord(i,3) ! vPx
            parameters_arr(6 + 4 * i) = particles_coord(i,4) ! vPy
        end do
        !$OMP END DO
        !$OMP END PARALLEL
        first_particle = 7
    end if


    ! CHECKS
    if (use_screen) then
        write (*,*) ACHAR(5)
        write (*,*) "---------- CHECKS ----------"
        write (*,*) ACHAR(5) 
    end if
    
    !!!! Checkeo rápido
    !!!! Check. Do I have to work?
    if (all(particles_MMR <= tini)) then
        if (use_screen) then
            write (*,*) ACHAR(10)
            write (*,*) "Initial particles condition a=0. Nothing to do here."
            write (*,*) "Saliendo."
        end if
        stop 1
    end if


    !!!! Mensaje de inicio
    if (use_screen) then
        if (error_tolerance <= 1.d-16) write (*,*) " WARNING: e_tol might be too low (<= 10⁻¹⁶)"
    end if


    !!!! Ahora si, si no hay que hacer más nada entonces terminamos.
    if (.not. any((/use_datascreen, use_datafile, use_chaosfile, use_potential_map, use_multiple_outputs/))) then ! No tiene sentido hacer nada
        write (*,*) ACHAR(10)
        write (*,*) "Como no se guardará nada, no se integrará"
        write (*,*) "Saliendo."
        stop 1
    end if

    !!!! Store initial conditions
    call save_initial(particles_initial_conditions, m0_and_boulders_initial_conditions, asteroid_initial_conditions)

    !!! Check if not too much multiple files
    if (use_multiple_outputs) then
        if (Nparticles > 1000) then
            write (*,*) ACHAR(5)
            write (*,*) "WARNING: Se crearán más de 1000 archivos de salida."
            write (*,*) "         Esto puede ser un problema."
            write (*,*) "Desear continuar? (s/n)"
            read (*,*) aux_character30
            aux_character1 = trim(adjustl(aux_character30))
            if ((aux_character1 == "y") .or. (aux_character1 == "s")) then
                if (use_screen) then
                    write (*,*) "Continuando..."
                    write (*,*) ACHAR(5)
                end if
            else
                write (*,*) "Saliendo."
                stop 1
            end if
        end if
    end if

    !!! Mensaje
    if (use_screen) then
        if (use_datafile) write (*,*) "Archivo general de salida: ", trim(datafile)
        if (use_multiple_outputs) write (*,*) "Archivos individuales de salida: ", trim(multfile) // "_*"
        if (use_datascreen) then
            write (*,*) "Habrán salidas en pantalla."
        else 
            write (*,*) "No habrán salidas en pantalla."
        end if
        write (*,*) ACHAR(5)
        write (*,*) "CKECHS OK"
        write (*,*) ACHAR(5)
    end if

    !!!! Mensaje de inicio
    if (use_screen) then
        write (*,*) ACHAR(5)
        write (*,*) "---------- INTEGRANDO ----------"
        write (*,*) ACHAR(5) 
    end if

    
    ! ABRIMOS ARCHIVOS
    !! Archivo de salida general
    if (use_datafile) open (unit=20, file=trim(datafile), status='unknown', action='write', access="append")
    !! Archivos individuales
    if (use_multiple_outputs) then
        do i = 0, Nboulders + Nparticles
            write (aux_character20, *) i
            open (unit=200+i, file=trim(multfile) // "_" // trim(adjustl(aux_character20)), &
                & status='unknown', action='write', access="append")
        end do
    end if
    if (use_update_chaos) open (unit=40, file=trim(chaosfile), status='unknown', action='readwrite', access="append")
    

     !!!!!!!!!!!!!!!!!!! Definimos vectores !!!!!!!!!!!!!!!!!!!!!!!!

    !! Definimos vector derivada
    if (use_explicit_method) then
        if (use_version_1) then
            dydt => dydt_explicit_v1
        else
            dydt => dydt_explicit_v2            
        end if
    else
        if (use_torque) then
            domegadt => dydt_single_null ! domegadt se calcula en la subrutina misma
        else if (omega_exp_damping_time < infinity) then
            domegadt => domega_dt_exponential
        else if (omega_linear_damping_time < infinity) then
            domegadt => domega_dt_linear
        else 
            domegadt => dydt_single_null
        end if
        if (use_version_1) then
            dydt => dydt_implicit_v1
        else
            dydt => dydt_implicit_v2
        end if
    end if
    !! Definimos vector merge
    if (use_merge) then
        resolve_merge => accumulate_mass_and_angmom
    else
        resolve_merge => do_not_accumulate_mass_and_angmom
    end if
    !! Definimos vectores de salidas
    if (use_datascreen) then
        if (use_elements_output) then
            write_b_to_screen => do_not_write
            write_i_to_screen => write_elements
        else
            write_b_to_screen => write_coordinates_boulders
            write_i_to_screen => write_coordinates_particle
        end if
    else 
        write_b_to_screen => do_not_write
        write_i_to_screen => do_not_write
    end if
    if (use_datafile) then
        if (use_elements_output) then
            write_b_to_general => do_not_write
            write_i_to_general => write_elements
        else
            write_b_to_general => write_coordinates_boulders
            write_i_to_general => write_coordinates_particle
        end if
    else 
        write_b_to_general => do_not_write
        write_i_to_general => do_not_write
    end if
    if (use_multiple_outputs) then
        if (use_elements_output) then
            write_b_to_individual => do_not_write
            write_i_to_individual => write_elements
        else
            write_b_to_individual => write_coordinates_boulders
            write_i_to_individual => write_coordinates_particle
        end if
    else 
        write_b_to_individual => do_not_write
        write_i_to_individual => do_not_write
    end if
    !! Definimos vector para elementos
    if (use_elements) then
        get_elements_i => calculate_elements_i
    else 
        get_elements_i => do_nothing_i
    end if
    !! Definimos vector para chaos
    if (use_chaos) then
        get_chaos_i => calculate_chaos_i
    else 
        get_chaos_i => do_nothing_i
    end if
    !!! Escribir chaos en cada output?
    if (use_update_chaos) then
        flush_chaos => write_chaos
    else
        flush_chaos => do_nothing_i
    end if


    ! Set parameters new to the derivative
    parameters_arr_new = dydt(initial_time, parameters_arr)

    !! Initial conditions Output
    !$OMP PARALLEL IF(my_threads > 1) DEFAULT(SHARED) &
    !$OMP PRIVATE(i)
    !$OMP DO
    do i = 0, Nboulders
        call write_b_to_individual(i, 200+i)
    end do 
    !$OMP END DO NOWAIT
    !$OMP DO
    do i = 1, Nparticles
        call write_i_to_individual(i, 200+i+Nboulders)
    end do 
    !$OMP END DO 
    !$OMP SECTIONS
    !$OMP SECTION
    do i = 0, Nboulders
        call write_b_to_screen(i, 6)
    end do
    do i = 1, Nparticles
        call write_i_to_screen(i, 6)
    end do
    !$OMP SECTION
    do i = 0, Nboulders
        call write_b_to_general(i, 20)
    end do
    do i = 1, Nparticles
        call write_i_to_general(i, 20)
    end do
    !$OMP END SECTIONS
    !$OMP END PARALLEL

    !!!!!! MAIN LOOP INTEGRATION !!!!!!!
    ! MAIN LOOP
    tom_index_number = 2 !!!! Inicializamos en 2 porque el primer checkpoint es el IC (t0)
    j = 2! From 2 because 1 is the IC (t0) !! +1 por si hay HardExit en el último
    main_loop: do while (.True.)

        ! Check for colissions/escapes
        staying_particles = 0
        discarded_particles = 0
        !!! Pre-set merges
        use_merge = use_merge .and. any(particles_mass(1:Nactive) > cero)
        mass_to_merge = cero
        angular_momentum_to_merge = cero
        !! Calculate distances to the asteroid
        !$OMP PARALLEL IF(my_threads > 1 .AND. Nactive > 10) DEFAULT(SHARED) &
        !$OMP PRIVATE(i) &
        !$OMP SHARED(discarded_particles)
        !$OMP DO REDUCTION(+:mass_to_merge,angular_momentum_to_merge) SCHEDULE (STATIC)
        do i = 1, Nactive
            if ((particles_dist(i) < min_distance) .or. (particles_hexit(i) .eq. 1)) then
                if (use_screen) then
                    write (*,*) "Colisión de la partícula ", i, "(", particles_index(i), ")", &
                    & " en t = ", time / unit_time, "[días]"
                    write (*,*) ACHAR(5)
                end if
                particles_outcome(i) = 1
                !$OMP CRITICAL (discard)
                discarded_particles = discarded_particles + 1
                ij_to_swap(discarded_particles,1) = i
                !$OMP END CRITICAL (discard)
                call resolve_merge(i, mass_to_merge, angular_momentum_to_merge)
            else if ((particles_dist(i) > max_distance) .or. (particles_hexit(i) .eq. 2)) then
                if (use_screen) then
                    write (*,*) ACHAR(5)
                    write (*,*) "Escape de la partícula ", i, "(", particles_index(i), ")", &
                    & " en t = ", time / unit_time, "[días]"
                    write (*,*) ACHAR(5)
                end if
                particles_outcome(i) = 2
                !$OMP CRITICAL (discard)
                discarded_particles = discarded_particles + 1
                ij_to_swap(discarded_particles,1) = i
                !$OMP END CRITICAL (discard)
            end if
            if (use_chaos) particles_times(i) = time ! Update particle times
        end do
        !$OMP END DO
        !$OMP END PARALLEL
        
        staying_particles = Nactive - discarded_particles
        particles_hexit = 0 ! Resetear HardExit
        !! Discard particles
        if (discarded_particles > 0) then
            if (use_merge) call merge_into_asteroid(mass_to_merge, angular_momentum_to_merge) ! Merge particles into asteroid
            call quicksort_int(ij_to_swap(1:discarded_particles,1), 1, discarded_particles) ! Sort particles to discard
            if (staying_particles == 0) then ! All particles are out
                if (use_screen) then
                    do i = Nactive, 1, -1
                        write (*,*) "  - Eliminando partícula ", i, "(", particles_index(i), ")"
                        write (*,*) ACHAR(5)
                    end do
                end if
                Nactive = 0 ! All particles are out
            else
                aux_integer = 1
                !$OMP PARALLEL IF((my_threads > 1) .AND. (discarded_particles > 20)) DEFAULT(SHARED) &
                !$OMP PRIVATE(i)
                !$OMP DO SCHEDULE (STATIC)
                do i = 1, discarded_particles
                    if (particles_outcome(Nactive + 1 - i) .ne. 0) then
                        ij_to_swap(discarded_particles + 1 - i, 2) = Nactive + 1 - i
                    else
                        !$OMP CRITICAL
                        ij_to_swap(aux_integer, 2) = Nactive + 1 - i
                        aux_integer = aux_integer + 1
                        !$OMP END CRITICAL
                    end if
                end do
                !$OMP END DO
                !$OMP BARRIER
                !$OMP DO SCHEDULE (STATIC)
                do i = 1, discarded_particles
                    if (use_screen) then
                        write (*,*) "  - Eliminando partícula ", ij_to_swap(i,1), "(", particles_index(ij_to_swap(i,1)), ")"
                        write (*,*) ACHAR(5)
                    end if
                    if (ij_to_swap(i,1) < (Nactive - discarded_particles)) then
                        call swap_particles(ij_to_swap(i,1), ij_to_swap(i,2), use_chaos)
                    end if
                end do
                !$OMP END DO
                !$OMP END PARALLEL
                Nactive = Nactive - discarded_particles
            end if
            if (use_screen) then
                write (*,*) "Quedan ", Nactive, " partículas activas."
                write (*,*) ACHAR(5)
            end if
            call argsort_int(particles_index(1:Nactive), sorted_particles_index(1:Nactive)) ! Get the sorted index
        end if
        
        ! Check if all done
        !! Time end
        if (j == checkpoint_number + 1) then
            if (use_screen) then
                write (*,*) ACHAR(5) 
                write (*,*) "Finalizó la integración en t = ", time / unit_time, "[días]"
            end if
            exit main_loop
        end if
        !! Particles left
        if (Nactive == 0) then
            if (use_screen) then
                write (*,*) ACHAR(5) 
                write (*,*) "No quedan partículas activas."
                write (*,*) "Finalizó la integración en t = ", time / unit_time, "[días]"
            end if
            exit main_loop
        end if

        
        ! Update dt
        timestep = checkpoint_times(j) - time
        !! Check TOM
        if (checkpoint_is_tom(j) .and. (tom_index_number <= tom_total_number)) then
            print*, "Dentro de TOM. Tiempo:", time
            print*, "Old Omega: ", asteroid_omega
            print*, "Old Mass: ", asteroid_mass
            if (allocated(tom_deltaomega)) then
                print*, "Delta Omega: ", tom_deltaomega(tom_index_number)
                asteroid_omega = asteroid_omega + tom_deltaomega(tom_index_number)
                asteroid_omega2 = asteroid_omega * asteroid_omega
                if (use_explicit_method) asteroid_theta_correction = asteroid_theta - asteroid_omega * (time - initial_time)
                if (use_version_1) then
                    do i = 0, Nboulders
                        aux_integer = i * equation_size
                        parameters_arr(aux_integer+3) = - pos_ast_arr(i,2) * asteroid_omega
                        parameters_arr(aux_integer+4) = pos_ast_arr(i,1) * asteroid_omega
                    end do
                else
                    parameters_arr(2) = asteroid_omega
                end if
            end if
            if (allocated(tom_deltamass)) then
                print*, "Delta Mass: ", tom_deltamass(tom_index_number)
                tom_mass_growth_param = uno + (tom_deltamass(tom_index_number) / asteroid_mass)
                mass_primary = mass_primary * tom_mass_growth_param
                mass_ast_arr = mass_ast_arr * tom_mass_growth_param
                asteroid_mass = asteroid_mass * tom_mass_growth_param
                asteroid_inertia = asteroid_inertia * tom_mass_growth_param
                Gasteroid_mass = G * asteroid_mass
            end if
            tom_index_number = tom_index_number + 1
            print*, "New Omega: ", asteroid_omega
            print*, "New Mass: ", asteroid_mass
            pause
        end if

        !!! Execute an integration method (uncomment/edit one of these)
        ! call integ_caller (time, parameters_arr, adaptive_timestep, dydt, &
        !     & Runge_Kutta4, timestep, parameters_arr_new, particles_hexitptr)
        ! call rk_half_step_caller (time, parameters_arr, adaptive_timestep, dydt, &
        !     & Runge_Kutta5, 5, error_tolerance, learning_rate, min_timestep, timestep, parameters_arr_new, particles_hexitptr)
        ! call embedded_caller (time, parameters_arr, adaptive_timestep, dydt, Dormand_Prince8_7, &
        !    & error_tolerance, learning_rate, min_timestep, timestep, parameters_arr_new, particles_hexitptr)
        call BStoer_caller (time, parameters_arr, adaptive_timestep, dydt, &
            & error_tolerance, timestep, parameters_arr_new, particles_hexitptr)

        !! If so, the dt used is in dt_adap
        if (particles_hexit(0) .ne. 0) timestep = adaptive_timestep

        ! Update parameters
        time = time + timestep
        parameters_arr = parameters_arr_new
        parameters_arr_new = dydt(time, parameters_arr)

        ! Asteroid and boulders
        if (use_explicit_method) then
            ! Constantes de Asteroid: pos, vel, omega, inertia, masa
            asteroid_theta = asteroid_omega * (time - initial_time) + asteroid_theta_correction!! No lo cambio ahora porque está en explicit_v2
            do i = 0, Nboulders
                pos_ast_arr(i,1) = cos(asteroid_theta + theta_ast_arr(i)) * dist_ast_arr(i)
                pos_ast_arr(i,2) = sin(asteroid_theta + theta_ast_arr(i)) * dist_ast_arr(i)
            end do
            vel_ast_arr(0:,1) = - asteroid_omega * pos_ast_arr(0:,2)
            vel_ast_arr(0:,2) =   asteroid_omega * pos_ast_arr(0:,1)
            acc_ast_arr = - asteroid_omega2 * pos_ast_arr
        else
            if (use_version_1) then
                ! Constantes de Asteroid (hasta ahora): masa
                !! Tendremos que obtener Asteroid (center of mass) properties (pos, vel, omega, ...)
                !!! Implicit V1
                do i = 0, Nboulders
                    aux_integer = i * equation_size
                    pos_ast_arr(i,:) = parameters_arr(aux_integer+1 : aux_integer+2)
                    vel_ast_arr(i,:) = parameters_arr(aux_integer+3 : aux_integer+4)
                    acc_ast_arr(i,:) = parameters_arr_new(aux_integer+3 : aux_integer+4)
                end do
                call get_asteroid_from_boulders(mass_ast_arr, pos_ast_arr, vel_ast_arr, &
                                                & asteroid_pos, asteroid_vel, asteroid_omega, asteroid_inertia)
                do i = 0, Nboulders
                    dist_ast_arr(i) = sqrt(sum((pos_ast_arr(i,:) - asteroid_pos)**2))
                end do
                ! Definimos theta como la variación del ángulo de m0 respecto del asteroide, respecto a la condición inicial
                asteroid_theta = mod(atan2(pos_ast_arr(0,2) - asteroid_pos(2), &
                                        & pos_ast_arr(0,1) - asteroid_pos(1)) - &
                                        & theta_ast_arr(0), twopi)
            else 
                ! Constantes de Asteroid (hasta ahora): masa, inertia
                !! Las Asteroid (center of mass) properties están servidas
                !!! Implicit V2
                parameters_arr(1) = mod(parameters_arr(1), twopi)
                asteroid_theta = parameters_arr(1)
                asteroid_omega = parameters_arr(2)
                asteroid_pos = parameters_arr(3:4)
                asteroid_vel = parameters_arr(5:6)
                pos_ast_arr(0:,1) = cos(asteroid_theta + theta_ast_arr(0:)) * dist_ast_arr(0:)
                pos_ast_arr(0:,2) = sin(asteroid_theta + theta_ast_arr(0:)) * dist_ast_arr(0:)
                vel_ast_arr(0:,1) = - asteroid_omega * pos_ast_arr(0:,2)
                vel_ast_arr(0:,2) =   asteroid_omega * pos_ast_arr(0:,1)
                do i = 0, Nboulders
                    pos_ast_arr(i,:) = asteroid_pos + pos_ast_arr(i,:)
                    vel_ast_arr(i,:) = asteroid_vel + vel_ast_arr(i,:)
                end do
                acc_ast_arr = - asteroid_omega2 * pos_ast_arr
            end if
        end if

        ! Update asteroid properties
        asteroid_omega2 = asteroid_omega * asteroid_omega
        asteroid_a_corot = (G * asteroid_mass / asteroid_omega2)**(1/3.)
        
        !! Particles
        !$OMP PARALLEL IF(my_threads > 1 .AND. Nactive > 1) DEFAULT(SHARED) &
        !$OMP PRIVATE(i,aux_integer)
        !$OMP DO SCHEDULE (STATIC)
        do i = 1, Nactive
            aux_integer =  first_particle + 4 * (i - 1)
            !!! Coordinates
            particles_coord(i,1:2) = parameters_arr(aux_integer   : aux_integer+1)
            particles_coord(i,3:4) = parameters_arr(aux_integer+2 : aux_integer+3)
            call get_elements_i(i) !!! Elements
            call get_chaos_i(i) !!! Chaos
            particles_dist(i) = sqrt(sum((particles_coord(i,1:2) - asteroid_pos)**2)) ! sqrt(x^2 + y^2)
        end do
        !$OMP END DO
        !$OMP END PARALLEL
        
        
        ! Output
        if ((checkpoint_is_output(j)) .and. (particles_hexit(0) .eq. 0)) then
            
            !$OMP PARALLEL IF(my_threads > 1) DEFAULT(SHARED) &
            !$OMP PRIVATE(i,aux_integer)
            !$OMP DO
            do i = 0, Nboulders
                call write_b_to_individual(i, 200+i)
            end do 
            !$OMP END DO NOWAIT
            !$OMP DO
            do i = 1, Nparticles
                call write_i_to_individual(i, 200+i+Nboulders)
            end do 
            !$OMP END DO 
            !$OMP SECTIONS
            !$OMP SECTION
            do i = 0, Nboulders
                call write_b_to_screen(i, 6)
            end do
            do i = 1, Nparticles
                call write_i_to_screen(i, 6)
            end do
            !$OMP SECTION
            do i = 0, Nboulders
                call write_b_to_general(i, 20)
            end do
            do i = 1, Nparticles
                call write_i_to_general(i, 20)
            end do
            !$OMP SECTION
            call flush_chaos(40) ! Update chaos
            !$OMP END SECTIONS
            !$OMP END PARALLEL
        end if

        if (use_percentage) call percentage(time, final_time)
        

        ! Update j
        j = j + 1

    end do main_loop

    !! Update surviving particles times
    if (use_chaos) particles_times(1:Nactive) = final_time ! Update particle times

    !! Porcentaje final
    if (use_percentage) call percentage(final_time+1., final_time)

    !! Cerrar archivo de salida
    if (use_datafile) then
        close (2)
        if (use_screen) then
            write (*,*) ACHAR(10)
            write (*,*) "Se guardó la salida en el archivo: ", trim(datafile)
        end if
    end if

    !! Cerrar archivos individuales
    if (use_multiple_outputs) then
        do i = 1, Nparticles
            close (200+i)
        end do
        if (use_screen) then
            write (*,*) ACHAR(10)
            write (*,*) "Se guardaron las salidas individuales en: ", trim(multfile) // "_*"
        end if
    end if

    !! Caos
    if (use_chaos) then
        if (use_chaosfile) then
            !! Guardar archivo de caos
            call argsort_int(particles_index, sorted_particles_index)
            open (unit=40, file=trim(chaosfile), status='unknown', action='write')!, access="append")
            do i = 1, Nparticles
                aux_integer = sorted_particles_index(i)
                write (40,*) aux_integer, & ! i
                & particles_outcome(aux_integer), & ! bad
                & final_time / unit_time, & ! total time to integrate
                & asteroid_initial_conditions(10) / (unit_mass * unit_dist * unit_vel), & ! initial (Asteroid): angular momentum
                & particles_initial_conditions(aux_integer,1) / unit_mass, & ! initial: mass
                & particles_initial_conditions(aux_integer,2) / unit_dist, & ! initial: a
                & particles_initial_conditions(aux_integer,3), & ! initial: e
                & particles_initial_conditions(aux_integer,4) / radian, & ! initial: M
                & particles_initial_conditions(aux_integer,5) / radian, & ! initial: omega
                & particles_initial_conditions(aux_integer,6), & ! initial: MMR
                & sqrt(particles_initial_conditions(aux_integer,2) * &
                  & (uno - particles_initial_conditions(aux_integer,3)**2) &
                  & * Gasteroid_mass) / (unit_dist * unit_vel), & ! initial: angular momentum per unit mass
                & particles_times(aux_integer) / unit_time, & ! surviving time
                & particles_elem(aux_integer,1) / unit_dist, particles_elem(aux_integer,2), & ! final: a, e
                & particles_elem(aux_integer,3) / radian, particles_elem(aux_integer,4) / radian, & ! final: M, omega
                & particles_MMR(aux_integer), & ! final: MMR
                & sqrt(particles_elem(aux_integer,1) * &
                  & (uno - particles_elem(aux_integer,2)**2) &
                  & * Gasteroid_mass) / (unit_dist * unit_vel), & ! final: angular momentum per unit mass
                & particles_min_a(aux_integer) / unit_dist, particles_max_a(aux_integer) / unit_dist, & ! a_min, a_max
                & particles_min_e(aux_integer), particles_max_e(aux_integer), & ! e_min, e_max
                & (particles_max_a(aux_integer) - particles_min_a(aux_integer)) / unit_dist, & ! Delta a
                & (particles_max_e(aux_integer) - particles_min_e(aux_integer)) ! Delta e
            end do
            close (40)
            !! Mensaje
            if (use_screen) then 
                write (*,*) ACHAR(10)
                write (*,*) "Se guardó el archivo de caos en: ", trim(chaosfile)
            end if
        end if
        if (use_datascreen) then
            write (*,*) ACHAR(10)
            if (use_single_particle) then            
                write (*,*) "Caos de la partícula simple [i, bad, tmax, MMR_ini, MMR_fin, a_min, a_max, Delta a, &
                & e_min, e_max, Delta_e]:"
            else 
                write (*,*) "Caos de las partículas [i, bad, tmax, MMR_ini, MMR_fin, a_min, a_max, Delta a, &
                & e_min, e_max, Delta_e]:"
            end if
            do i = 1, Nparticles
                aux_integer = sorted_particles_index(i)
                write (*,*) aux_integer, & ! i
                & particles_outcome(aux_integer), & ! bad
                & particles_times(aux_integer), & ! surviving time
                & particles_initial_conditions(aux_integer,6), & ! initial: MMR
                & particles_MMR(aux_integer), & ! final: MMR
                & particles_min_a(aux_integer) / unit_dist, particles_max_a(aux_integer) / unit_dist, & ! a_min, a_max
                & particles_min_e(aux_integer), particles_max_e(aux_integer), & ! e_min, e_max
                & (particles_max_a(aux_integer) - particles_min_a(aux_integer)) / unit_dist, & ! Delta a
                & (particles_max_e(aux_integer) - particles_min_e(aux_integer)) ! Delta e
            end do
        end if
    end if

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!! LIBERACION MEMORIA!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (use_screen) then
        write (*,*) ACHAR(10)
        write (*,*) "Liberando memoria"
        write (*,*) ACHAR(5)
    end if
    !!!!!!!!! TIEMPOS !!!!!!!!
    call free_times_arrays()
    if (use_tomfile) then
        deallocate(tom_times)
        if (allocated(tom_deltaomega)) deallocate(tom_deltaomega)
        if (allocated(tom_deltamass)) deallocate(tom_deltamass)
    end if
    !!!!!!!! Asteroide  !!!!!!!
    call free_asteroid_arrays()
    !!!!!!!! PARTICULAS !!!!!!!
    call free_particles()
    !!!!!!!! PARAMETROS !!!!!!!
    deallocate(parameters_arr)
    deallocate(parameters_arr_new)
    !!!!!!!!! INITIALS !!!!!!!!
    call free_initial_conditions()
    !!!!!!!!! POINTERS !!!!!!!!
    call nullify_pointers()
    nullify(dydt)
    nullify(domegadt)
    nullify(dmassdt)

    if (use_screen) then 
        write (*,*) ACHAR(5)
        write (*,*) "Fin del programa"
        write (*,*) ACHAR(5)
    end if
    
end program main

subroutine percentage (tout,tstop)   ! version para gfortran
    implicit none
    real(kind=8), intent(in) :: tout, tstop
    integer(kind=4) :: iper, i
    character(len=1)  :: cret
    character(len=99) :: guiones
    character(len=3)  :: cestado
    save :: iper,guiones
    
    cret = achar(13)          ! generate carriage return
    
    iper = int(100.0*tout/tstop)
    guiones = ''
    do i = 1,iper
        guiones = trim(guiones)//'.'
    end do
    
    open (66,status='scratch')
    if (iper < 10) then
        write (66,'(i2)') iper
    else
        write (66,'(i3)') iper
    end if
    rewind (66)
    read (66,'(a)') cestado
    close (66)
    
    110 format (2a)
    if (iper < 100) then
        write (*,110,advance='no') cret,trim(guiones)//trim(cestado)//'%'
    else
        write (*,110,advance='no') cret,trim(guiones)//'. FIN'
        write (*,*)
    end if
end subroutine percentage

subroutine create_map(Nboul,ngx,ngy,xmin,xmax,ymin,ymax,mib,rib,vib,map_file)
    use parameters, only: G, cero, radius_primary, use_screen
    use forces, only: apply_force
    implicit none
    integer(kind=4), intent(in) :: Nboul, ngx, ngy
    real(kind=8), intent(in)    :: xmin, xmax, ymin, ymax
    real(kind=8), intent(in)    :: rib(0:Nboul,2), vib(0:Nboul,2), mib(0:Nboul)
    character(len=*), intent(in) :: map_file
    real(kind=8) :: xyb(ngx,ngy), acb(ngx,ngy,2)
    real(kind=8) :: rb(2), dist
    real(kind=8) :: dummy, dummy2(2)
    real(kind=8) :: pot_at_R, asteroid_mass
    integer(kind=4) :: i, j, k, dumint
    real(kind=8), dimension(2) :: cero2=cero

    !!! Check
    if ((ngx < 2) .or. (ngy < 2)) then
        write (*,*) "ERROR: Verificar que ngx > 2 y ngy > 2"
        stop 1
    end if
    if ((xmin >= xmax) .or. (ymin >= ymax)) then
        write (*,*) "ERROR: Verificar que xmin < xmax y ymin < ymax"
        stop 1
    end if
    xyb = cero
    acb = cero
    pot_at_R = - mib(0) / radius_primary
    asteroid_mass = sum(mib)
    do i = 1, ngx
        do j = 1, ngy
            rb(1) = xmin + i * (xmax - xmin) / ngx
            rb(2) = ymin + j * (ymax - ymin) / ngy
            dist = sqrt((rb(1)-rib(0,1))**2 + (rb(2)-rib(0,2))**2)
            if (dist .le. radius_primary) then
                xyb(i,j) = pot_at_R
                acb(i,j,:) = cero
                cycle
            end if
            xyb(i,j) = - mib(0) / dist
            do k = 1, Nboul
                dist = sqrt((rb(1)-rib(k,1))**2 + (rb(2)-rib(k,2))**2)
                xyb(i,j) = xyb(i,j) - mib(k) / dist
            end do
            call apply_force(cero, mib, cero, dumint, rb, cero2, rib, vib, asteroid_mass, &
                            &  cero2, cero2, dummy, dummy, dummy2, acb(i,j,:))
        end do
    end do
    xyb = xyb * G
    if (use_screen) write (*,*) "   Potencial calculado. Escribiendo..."
    open (unit=50, file=trim(map_file), status='replace', action='write')
    do i = 1, ngx
        do j = 1, ngy
            write (50,*) xmin + i * (xmax - xmin) / ngx, ymin + j * (ymax - ymin) / ngy, &
            xyb(i,j), acb(i,j,:)
        end do
    end do
    close (50)
end subroutine create_map

! module iso_fortran_env

!     ! Nonintrinsic version for Lahey/Fujitsu Fortran for Linux. 
!     ! See Subclause 13.8.2 of the Fortran 2003 standard. 
  
!     implicit none 
!     public 
  
!     integer, parameter :: Character_Storage_Size = 8 
!     integer, parameter :: Error_Unit = 0 
!     integer, parameter :: File_Storage_Size = 8 
!     integer, parameter :: Input_Unit = 5 
!     integer, parameter :: IOSTAT_END = -1 
!     integer, parameter :: IOSTAT_EOR = -2 
!     integer, parameter :: Numeric_Storage_Size = 32 
!     integer, parameter :: Output_Unit = 6 
  
!   end module iso_fortran_env

