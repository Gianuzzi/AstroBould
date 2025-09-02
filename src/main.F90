program main
    use parameters
    use integrators
    use times
    ! use forces
    
    implicit none
    ! Parameters
    type(input_params_st) :: input_params  ! This is the structure with the default -> input parameters
    type(input_params_st) :: sim  ! This is the structure with the input -> derived parameters
    ! Map
    external :: create_map   ! Declare the external function
    logical :: only_potential_map = .False.
    ! Auxiliar
    integer(kind=4) :: aux_int
    character(1) :: aux_character1
    character(20) :: aux_character20
    character(30) :: aux_character30
    real(kind=8) :: aux_real
    real(kind=8), dimension(4) :: aux_real4
    real(kind=8), dimension(:,:), allocatable :: boulders_coords !! mass, radius, theta  | (Nb+1, 4)
    real(kind=8), dimension(:,:), allocatable :: aux_2D  ! To read files
    ! Loops
    integer(kind=4) :: i

    call init_default_parameters()  ! Inicializamos parámetros globales básicos

    use_configfile = .True. ! Usar archivo de configuración
    arguments_number = command_argument_count()
    do i = 1, arguments_number
        call get_command_argument(i, aux_character30)
        if (trim(aux_character30) .eq. "--noconfig") then
            use_configfile = .False.
            exit
        end if
    end do

    if (use_configfile) call read_config_file(input_params, "config.ini", configfile_exists) ! Leemos archivo de configuración

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!! Editar aquí abajo de ser necesario !!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (.not. use_configfile) then ! Usaremos parámetros por defecto y 
        
        ! Asteroide central
        !! Primary
        input_params%mass_primary = 6.3d18 ! Masa del cuerpo 0 [kg] ! -x =>  mAst = x
        input_params%radius_primary = 129.d0 ! Radio del cuerpo 0 [km]

        !! ROTACIÓN
        input_params%asteroid_rotational_period = 7.004d0/24.d0  ! Periodo de rotación [day]
        !lambda_kep = 0.471d0      ! Cociente omega/wk

        !! Boulders        
        input_params%Nboulders = 1 ! Número de boulders

        if (input_params%Nboulders > 0) then
            call allocate_params_asteroid(input_params%Nboulders)!! Alocatamos (No tocar)

            boulders_in(1,1) = 1.d-1   ! Cociente de masas entre boulder 1 y primary
            boulders_in(1,2) = cero    ! Ángulo de fase del boulder 1 [deg]
            boulders_in(1,3) = 2.5d0   ! Radio del boulder 1 [km]       

        end if

        !! Moons
        input_params%Nmoons = 0 ! Número de boulders

        if (input_params%Nmoons > 0) then
            call allocate_params_moons(input_params%Nmoons)!! Alocatamos (No tocar)

            moons_in(1,1) = 1.d-6   ! Cociente de masas entre luna 1 y asteroide
            moons_in(1,2) = cero    ! Semieje [km]
            moons_in(1,3) = cero    ! Eccentricidad
            moons_in(1,4) = cero    ! M [deg]
            moons_in(1,5) = cero    ! w [deg]
            moons_in(1,6) = 11.1d0  ! MMR
            moons_in(1,7) = cero    ! radius [km]

        end if

        !! Particles
        input_params%Nparticles = 0 ! Número de boulders

        if (input_params%Nparticles > 0) then
            call allocate_params_particles(input_params%Nparticles)!! Alocatamos (No tocar)

            particles_in(1,1) = cero   ! Semieje [km]
            particles_in(1,2) = cero   ! Eccentricidad
            particles_in(1,3) = cero   ! M [deg]
            particles_in(1,4) = cero   ! w [deg]
            particles_in(1,5) = 8.1d0  ! MMR
            
        end if

        !!!! Merges (colliding particles into asteroid)
        input_params%use_merge = .False. ! Merge particles into asteroid
        
        !!!! Stokes
        input_params%use_stokes = .False.
        input_params%stokes_a_damping_time = infinity                                  ! [day]
        input_params%stokes_e_damping_time = input_params%stokes_a_damping_time / 1.d2 ! [day]
        input_params%stokes_charac_time = cero                                         ! [day] Tiempo que actúa stokes

        !!!! Naive-Stokes (drag)
        input_params%use_naive_stokes = .False.
        input_params%drag_coefficient = cero  ! Eta
        input_params%drag_charac_time = cero  ! [day] Tiempo que actúa drag

        !!!! Geo-Potential (J2)
        input_params%J2_coefficient = cero

        !!! Parámetros corrida

        !!!! Tiempos
        input_params%initial_time = cero      ! Initial time [day]
        input_params%final_time = 2.d3        ! Final time [day]
        input_params%case_output_type = 0     ! 0: Linear ; 1: Logarithmic ; 2: Combination
        input_params%output_number = 10000    ! Number of outputs (if dt_out=0)
        input_params%output_timestep = cero   ! Output timestep [day] (used if case_output_type != 1)

        !!!! Error
        input_params%learning_rate = 0.85d0   ! [For adaptive step integrators] Learning rate
        input_params%error_digits = 12        ! [For adaptive step integrators] Digits for relative error

        !!!! Colision y escape
        input_params%min_distance = -1.5d0    ! Min distance before impact [km] ! 0 => R0 + max(Rboul)
        input_params%max_distance = -1.d2     ! Max distance before escape [km] ! -x => R0 * x

        !!! Output: "" or "no", if not used
        input_params%datafile = ""
        input_params%chaosfile = ""
        input_params%mapfile = ""
        input_params%multfile = ""

        !!!!! Screeen
        input_params%use_screen = .True. ! Print info in screen
        input_params%use_datascreen = .True. ! Print data in screen

        !!! Input: "" or "no", if not used
        input_params%tomfile = ""
        input_params%moonsfile = ""
        input_params%particlesfile = ""

        !!! Parallel
        input_params%use_parallel = .False.
        input_params%requested_threads = 1 ! Number of threads to use !! -1 => all available
    end if


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!! No tocar de aquí a abajo !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    ! Load command line arguments
    call load_command_line_arguments(input_params, use_configfile)
    
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
    sim = input_params
    call set_derived_parameters(sim) ! Inicializamos parámetros derivados

    
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
    if (sim%use_moonsfile) then
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
            write (*,f12) "  Read ", sim%Nmoons, " rows."
            write (*,f12) "  Read ", aux_int, " columns."
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
    if (sim%use_particlesfile .and. use_command_particle) then
        write (*,*) ACHAR(10)
        write (*,*) "WARNING: Command line particle used. Ignoring particles file: ", trim(sim%particlesfile)
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
            write (*,f12) "  Read ", sim%Nparticles, " rows."
            write (*,f12) "  Read ", aux_int, " columns."
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
    sim%Nactive = sim%Ntotal

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

    if (.not. sim%use_boulders) then
        write (*,*) "WARNING: No boulders to integrate. Rotation effects dissabled."
    end if

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  Parallel  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (sim%use_screen) then
        write (*,*) ACHAR(5)
        if (sim%use_parallel) then
            write (*,*) "---------- Parallel integration----------"
            write (*,f12) "  Available threads:", available_threads
            write (*,f12) "  Requested threads:", sim%requested_threads
            write (*,f12) "  Used threads     :", my_threads
        else
            write (*,*) " Serial integration (No parallel)"
        end if
        write (*,*) ACHAR(5)
    end if
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!! COMIENZO DE CÁLCULOS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (sim%use_screen) then
        write (*,f12) "Starting simulation number: ", simulation_number
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
                    & sim%radius_primary * unit_dist)  ! radius
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
                         & particles_in(i,5))                ! MMR
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
        write (*,f13) "  a_corot               :", system%asteroid%a_corotation / unit_dist, "[km]"
        write (*,f13) "  omega                 :", system%asteroid%omega * unit_time, "[rad/day]"
        write (*,f13) "  omega_kep             :", system%asteroid%omega_kep * unit_time, "[rad/day]"
        write (*,f13) "  lambda_kep            :", system%asteroid%lambda_kep
        write (*,f13) "  Period                :", system%asteroid%rotational_period * 24 * unit_time, "[hs]"
        write (*,f13) "  Inertial momentum     :", system%asteroid%inertia / (unit_mass * unit_dist**2), "[kg km^2]"
        write (*,f13) "  Angular momentum (rot):", system%asteroid%ang_mom_rot / unit_mass / &
                                    & (unit_dist**2) * unit_time, "[kg km^2 / day]"
        write (*,*) ACHAR(5)
        if (.not. sim%use_boulders) then
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
            write (*,f233) i, &
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
            write (*,f233) i, &
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
            write (*,f233) i, &
                     & system%particles(i)%coordinates(1:2) / unit_dist, &
                     & system%particles(i)%coordinates(3:4) / unit_vel, &
                     & system%particles(i)%dist_to_cm / unit_dist
        end do
        write (*,*) ACHAR(5)
    end if

    ! Energy and ang_mom
    if (sim%use_screen) then
        write (*,f13) "   Energy          : ", system%energy
        write (*,f13) "   Angular Momentum: ", system%ang_mom
        write (*,*) ACHAR(5)
    end if

    ! <<<< Save initial data >>>>
    initial_system = system
    if (sim%use_boulders) then
        allocate(boulders_theta_ini(0:sim%Nboulders))
        do i = 0, sim%Nboulders
            boulders_theta_ini(i) = system%asteroid%boulders(i)%initial_theta
        end do 
    end if


    ! <> Parallel revisted
    ! Redefine OMP threads if necessary. Lower equal to particle number
    if ((sim%use_parallel) .and. ((sim%Ntotal - 1) < my_threads)) then
        my_threads = sim%Ntotal - 1
        if (sim%use_screen) then
            write (*,f12) "WARNING: Threads reduced to ", my_threads
            write (*,*) ACHAR(5)
        end if
        !$ call omp_set_num_threads(my_threads)
    end if


    !Test MERGE
    ! call calculate_coordinates(system, aux_real, aux_real4)
    print*, "Merge Ast"
    call merge_moon_i_into_ast(system, 1)
    print*, "Merge Moon"
    call merge_2_moons(system, 1, 2)
    ! call calculate_coordinates(system, aux_real, aux_real4)


!     !!!!!!!!!!!!!!!!!!!!!!!!!!!! EFECTOS EXTERNOS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (sim%use_screen) then
        write (*,*) ACHAR(5)
        write (*,*) ("---------- External/Internal forces ----------")
        write (*,*) ACHAR(5)
    end if

!     !! Variación en omega (funciones del tiempo)
!     if (use_omega_damping) then
!         ! Linear omega damping
!         if ((abs(omega_linear_damping_time) > tini) .and. (omega_linear_damping_time < infinity)) then
!             if (use_explicit_method) then
!                 write (*,*) ACHAR(10)
!                 write (*,*) "ERROR: No se puede usar el método explícito con tau_o finito."
!                 write (*,f13) "tau_o [Prot]:", omega_linear_damping_time * unit_time / asteroid_rotational_period
!                 stop 1
!             end if
!             if (.not. use_boulders) then
!                 write (*,*) "WARNING: Omega damping has no effect without boulders. It won't be used"
!                 omega_linear_damping_time = infinity
!                 omega_linear_damping_slope = cero
!             else
!                 omega_linear_damping_time = omega_linear_damping_time * unit_time
!                 omega_linear_damping_slope = - asteroid_omega / (omega_linear_damping_time - initial_time)
!                 if (use_screen) then
!                     write (*,*) "Explicit omega linear damping"
!                     write (*,f13)  "    tau_o :", omega_linear_damping_time / asteroid_rotational_period, "[Prot]"
!                 end if
!             end if
!         else
!             omega_linear_damping_time = infinity
!             omega_linear_damping_slope = cero
!         end if
!         ! Exponential omega damping
!         if ((abs(omega_exp_damping_time) > tini) .and. (omega_exp_damping_time < infinity)) then
!             if (use_explicit_method) then
!                 write (*,*) ACHAR(10)
!                 write (*,*) "ERROR: No se puede usar el método explícito con tau_o finito."
!                 write (*,f13) "tau_o [Prot]:", omega_exp_damping_time * unit_time / asteroid_rotational_period
!                 stop 1
!             end if
!             if (.not. use_boulders) then
!                 write (*,*) "WARNING: Omega damping has no effect without boulders. It won't be used"
!                 omega_exp_damping_time = infinity
!             else
!                 omega_exp_damping_time = omega_exp_damping_time * unit_time
!                 if (use_screen) then
!                     write (*,*) "Explicit omega exponential damping"
!                     write (*,f13)  "    tau_o :", omega_exp_damping_time / asteroid_rotational_period, "[Prot]"
!                 end if
!             end if
!         else
!             omega_exp_damping_time = infinity
!         end if
!         ! Poly_exponential omega damping
!         if ((abs(omega_exp_poly_A) > tini) .and. (abs(omega_exp_poly_B) > tini)) then
!             if (use_explicit_method) then
!                 write (*,*) ACHAR(10)
!                 write (*,*) "ERROR: No se puede usar el método explícito con tau_o finito."
!                 write (*,f13) "poly exp coeff A:", omega_exp_poly_A
!                 write (*,f13) "poly exp coeff B:", omega_exp_poly_B
!                 stop 1
!             end if
!             if (.not. use_boulders) then
!                 write (*,*) "WARNING: Omega damping has no effect without boulders. It won't be used"
!                 omega_exp_poly_A = cero
!                 omega_exp_poly_B = cero
!             else
!                 if (use_screen) then
!                     write (*,*) "Explicit omega poly-exponential damping"
!                     write (*,f13)  "    A :", omega_exp_poly_A
!                     write (*,f13)  "    B :", omega_exp_poly_B
!                 end if
!             end if
!             omega_exp_poly_AB = omega_exp_poly_A * omega_exp_poly_B ! Define AB
!         else
!             omega_exp_poly_A = cero
!             omega_exp_poly_B = cero
!         end if
!     else
!         omega_linear_damping_time = infinity
!         omega_linear_damping_slope = cero
!         omega_exp_damping_time = infinity
!         omega_exp_poly_A = cero
!         omega_exp_poly_B = cero
!     end if

!     !! Mass exponential damping
!     mass_exp_damping_time = cero !!! Deprecado
!     if ((abs(mass_exp_damping_time) > tini) .and. (mass_exp_damping_time < infinity)) then
!         if (use_explicit_method) then
!             write (*,*) ACHAR(10)
!             write (*,*) "ERROR: No se puede usar el método explícito con tau_m finito."
!             write (*,f13) "tau_m [Prot]:", mass_exp_damping_time * unit_time / asteroid_rotational_period
!             stop 1
!         end if
!         mass_exp_damping_time = mass_exp_damping_time * unit_time
!         if (use_screen) then
!             write (*,*) "Explicit mass exponential damping"
!             write (*,f13)  "    tau_m :", mass_exp_damping_time / asteroid_rotational_period, "[Prot]"
!         end if
!     else
!         mass_exp_damping_time = infinity
!     end if

!     !!! Check not both omega dampings
!     if ((omega_exp_damping_time < infinity .and. omega_linear_damping_time < infinity) .or. &
!       & (omega_exp_damping_time < infinity .and. abs(omega_exp_poly_A) > tini) .or. &
!       & (omega_linear_damping_time < infinity .and. abs(omega_exp_poly_A) > tini)) then
!         write (*,*) ACHAR(10)
!         write (*,*) "ERROR: No se puede usar ambos decaimientos de omega (linear and exp) finitos."
!         stop 1
!     end if

!     !!! Check que no esté torque también
!     if (use_torque .and. use_omega_damping) then
!         write (*,*) ACHAR(10)
!         write (*,*) "ERROR: No se puede usar torque y tau_o finito."
!         stop 1
!     end if
!     !!!

!     !! Stokes
!     if (use_stokes) then
!         call set_stokes_C_and_alpha(stokes_a_damping_time, stokes_e_damping_time, stokes_C, stokes_alpha)
!         stokes_a_damping_time = stokes_a_damping_time * unit_time
!         stokes_e_damping_time = stokes_e_damping_time * unit_time
!         stokes_charac_time = stokes_charac_time * unit_time
!         !!! Mensaje
!         if (use_screen) then
!             write (*,*) "Stokes"
!             write (*,f13)  "    t_stokes: ", stokes_charac_time / unit_time, "[days]"
!             write (*,f13)  "    C       : ", stokes_C
!             write (*,f13)  "    alpha   : ", stokes_alpha
!         end if
!     else ! Just to be sure
!         stokes_a_damping_time = infinity
!         stokes_e_damping_time = infinity
!         stokes_charac_time = cero
!     end if

!     !! J2
!     J2_effective = 1.5d0 * J2_coefficient ! Define effective J2

!     !! Torque
!     if ((.not. use_boulders) .and. use_torque) then
!         if (use_screen) write (*,*) "WARNING: Torque has no effect without boulders. It won't be used"
!         use_torque = .False.
!     end if

!     !! Mensaje !
!     if (use_screen) then
!         if (use_naive_stokes) then !!!! Naive-Stokes
!             write (*,*) "Naive-Stokes (drag radial)"
!             write (*,f13) "    Eta :", drag_coefficient
!             write (*,f13) " t_naive:", stokes_charac_time / unit_time, "[days]"
!         end if
!         if (use_J2) then
!             write (*,*) "Geo-Potential (J2)"
!             write (*,f13) "    J2 :", J2_coefficient !!!! Geo-Potential (J2)
!         end if
!         if (use_torque) write (*,*) "Torque from particles to asteroid ACTIVATED"

!         if (.not. any((/use_stokes, use_naive_stokes, use_J2, use_torque, &
!             & use_omega_damping, &
!             & mass_exp_damping_time < infinity/))) then
!             write (*,*) "No se aplicarán efectos ex/internos."
!         end if
!     end if

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
    if (sim%max_distance <= sim%min_distance) then
        write (*,*) ACHAR(10)
        write (*,*) "ERROR: rmax <= rmin"
        stop 1
    end if
    if (sim%use_screen) then
        write (*,*) "Conditions for escape/collision"
        write (*,f13) "    rmin : ", sim%min_distance / unit_dist, "[km] =", sim%min_distance / system%asteroid%radius, "[Rast]"
        write (*,f13) "    rmax : ", sim%max_distance / unit_dist, "[km] =", sim%max_distance / system%asteroid%radius, "[Rast]"
        if (sim%use_merge) then
            write (*,*) "  Collisions will be solved as mergers to the asteroid."
        else
            write (*,*) "  Collisions will not be sovled."
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
        call create_map(system, &
                      & sim%map_grid_size_x, sim%map_grid_size_y, &
                      & sim%map_min_x, sim%map_max_x, &
                      & sim%map_min_y, sim%map_max_y, &
                      & sim%mapfile)
        if (sim%use_screen) then
            write (*,*) "Maps stored into file: ", trim(sim%mapfile)
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

!     !! Chaos and Outcome
!     particles_outcome = 0
!     particles_hexit = 0



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
        if (allocated(tom_deltaomega) .and. (.not. sim%use_boulders)) then
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
        adaptive_timestep = 1.d0 ! Heuristic
    end if

    !! Mensaje
    if (sim%use_screen) then
        write (*,*) "Times:"
        if (system%asteroid%rotational_period > tini) then
            write (*,f13) "    t0    : ", sim%initial_time / system%asteroid%rotational_period, "[Prot] = ", &
            & sim%initial_time / unit_time, "[day]"
            write (*,f13) "    tf    : ", sim%final_time / system%asteroid%rotational_period, "[Prot] = ", &
            & sim%final_time / unit_time, "[day]"
            write (*,f13) "    dt_out: ", sim%output_timestep / system%asteroid%rotational_period, "[Prot] = ", &
            & sim%output_timestep / unit_time, "[day]"
            write (*,f13) "    dt_min: ", min_timestep / system%asteroid%rotational_period, "[Prot] = ", &
            & min_timestep / unit_time, "[day]"
        else 
            write (*,f13) "    t0    : ", sim%initial_time / unit_time, "[day]"
            write (*,f13) "    tf    : ", sim%final_time / unit_time, "[day]"
            write (*,f13) "    dt_out: ", sim%output_timestep / unit_time, "[day]"
            write (*,f13) "    dt_min: ", min_timestep / unit_time, "[day]"
        end if
        write (*,f12) "    n_out : ", sim%output_number
        write (*,*) ACHAR(5)
    end if




    !!!!!!!!!!!!!!!!!!!!!!!!!!! Integration Arrays !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Allocate
    allocate(m_arr(sim%Ntotal))
    allocate(y_arr(2 + sim%Ntotal * 4))
    allocate(y_arr_new(2 + sim%Ntotal * 4))
    allocate(y_der(2 + sim%Ntotal * 4))

    ! Get values
    call generate_mass_and_coordinates_arrays(system, m_arr, y_arr)

    ! first_particle = 2 + 4 * (system%Nmoons_active) + 1
    ! last_particle = 2 + 4 * (system%Nmoons_active + system%Nparticles_active) + 1

    
    
    ! !!!!!!!! VECTOR A INTEGRAR !!!!!!!
    ! !!!! Inicializamos vector
    ! if (use_version_1 .or. (.not. use_boulders)) then
    !     !!!!! Version 1: [x0, y0, vx0, vy0, x1, y1, vx1, vy1, ...]
    !     allocate(parameters_arr(4 * Ntotal))
    !     allocate(parameters_arr_new(4 * Ntotal))
    !     allocate(parameters_der(4 * Ntotal))
    !     do i = 0, Nboulders
    !         parameters_arr(1 + 4 * i) = pos_ast_arr(i,1)
    !         parameters_arr(2 + 4 * i) = pos_ast_arr(i,2)
    !         parameters_arr(3 + 4 * i) = vel_ast_arr(i,1)
    !         parameters_arr(4 + 4 * i) = vel_ast_arr(i,2)
    !     end do
    !     !$OMP PARALLEL DEFAULT(SHARED) &
    !     !$OMP PRIVATE(i)
    !     !$OMP DO SCHEDULE (STATIC)
    !     do i = 1, Nparticles
    !         parameters_arr(1 + 4 * (i + Nboulders)) = particles_coord(i,1) ! xP
    !         parameters_arr(2 + 4 * (i + Nboulders)) = particles_coord(i,2) ! yP
    !         parameters_arr(3 + 4 * (i + Nboulders)) = particles_coord(i,3) ! vPx
    !         parameters_arr(4 + 4 * (i + Nboulders)) = particles_coord(i,4) ! vPy
    !     end do
    !     !$OMP END DO
    !     !$OMP END PARALLEL
    ! else 
    !     !!!!! Version 2: [theta, omega, xA, yA, vxA, vyA, Part...]
    !     allocate(parameters_arr(6 + 4 * Nparticles))
    !     allocate(parameters_arr_new(6 + 4 * Nparticles))
    !     parameters_arr(1) = asteroid_theta
    !     parameters_arr(2) = asteroid_omega
    !     parameters_arr(3:4) = asteroid_pos
    !     parameters_arr(5:6) = asteroid_vel
    !     !$OMP PARALLEL DEFAULT(SHARED) &
    !     !$OMP PRIVATE(i)
    !     !$OMP DO SCHEDULE (STATIC)
    !     do i = 1, Nparticles
    !         parameters_arr(3 + 4 * i) = particles_coord(i,1) ! xP
    !         parameters_arr(4 + 4 * i) = particles_coord(i,2) ! yP
    !         parameters_arr(5 + 4 * i) = particles_coord(i,3) ! vPx
    !         parameters_arr(6 + 4 * i) = particles_coord(i,4) ! vPy
    !     end do
    !     !$OMP END DO
    !     !$OMP END PARALLEL
    !     first_particle = 7
    ! end if
    ! !! Define last particle ( y -> parameters(:last_particle) )


    ! CHECKS
    if (sim%use_screen) then
        write (*,*) ACHAR(5)
        write (*,*) "---------- CHECKS ----------"
        write (*,*) ACHAR(5) 
    end if
    
!     !!!! Checkeo rápido
!     !!!! Check. Do I have to work?
!     if (all(particles_MMR <= tini)) then
!         if (use_screen) then
!             write (*,*) ACHAR(10)
!             write (*,*) "Initial particles condition a=0. Nothing to do here."
!             write (*,*) "Saliendo."
!         end if
!         stop 1
!     end if

    !!! Check workers
    if (sim%use_parallel) then
        aux_int = MIN(MAX(sim%Nactive, 3), my_threads)
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

!     !!!! Store initial conditions
!     call save_initial(particles_initial_conditions, m0_and_boulders_initial_conditions, asteroid_initial_conditions)

    !!! Check if not too much multiple files
    if (sim%use_multiple_outputs) then
        if (sim%Nactive > 1000) then
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
        write (*,*) "CKECHS OK"
        write (*,*) ACHAR(5)
    end if

    !!!! Mensaje de inicio
    if (sim%use_screen) then
        write (*,*) ACHAR(5)
        write (*,*) "---------- INTEGRATING ----------"
        write (*,*) ACHAR(5) 
    end if

    
    ! ABRIMOS ARCHIVOS
    !! Archivo de salida general
    if (sim%use_datafile) open (unit=20, file=trim(sim%datafile), status='unknown', action='write', access="append")
    !! Archivos individuales
    if (sim%use_multiple_outputs) then
        do i = 0, sim%Nactive
            write (aux_character20, *) i
            open (unit=200+i, file=trim(sim%multfile) // "_" // trim(adjustl(aux_character20)), &
                & status='unknown', action='write', access="append")
        end do
    end if
    if (use_update_chaos) open (unit=40, file=trim(sim%chaosfile), status='unknown', action='readwrite', access="append")
    

!      !!!!!!!!!!!!!!!!!!! Definimos vectores !!!!!!!!!!!!!!!!!!!!!!!!

!     !! Definimos vector derivada
!     if (use_explicit_method) then
!         if (use_version_1) then
!             dydt => dydt_explicit_v1
!             check_continue_ptr => check_continue_v1
!         else
!             dydt => dydt_explicit_v2
!             check_continue_ptr => check_continue_v2
!         end if
!     else
!         if (use_torque) then
!             domegadt => dydt_single_null ! domegadt se calcula en la subrutina misma
!         else if (omega_exp_damping_time < infinity) then
!             domegadt => domega_dt_exponential
!         else if (omega_linear_damping_time < infinity) then
!             domegadt => domega_dt_linear
!         else if (abs(omega_exp_poly_A) > tini) then
!             domegadt => domega_dt_expoly
!         else 
!             domegadt => dydt_single_null
!         end if
!         if (use_version_1) then
!             dydt => dydt_implicit_v1
!             check_continue_ptr => check_continue_v1
!         else
!             dydt => dydt_implicit_v2
!             check_continue_ptr => check_continue_v2
!         end if
!     end if
!     !! Redeinimos vector si hay boulders
!     if (.not. use_boulders) dydt => dydt_no_boulders
!     !! Definimos vector merge
!     if (use_merge) then
!         resolve_merge => accumulate_mass_and_angmom
!     else
!         resolve_merge => do_not_accumulate_mass_and_angmom
!     end if
!     !! Definimos vectores de salidas
!     if (use_datascreen) then
!         if (use_elements_output) then
!             write_b_to_screen => do_not_write
!             write_i_to_screen => write_elements
!         else
!             write_b_to_screen => write_coordinates_boulders
!             write_i_to_screen => write_coordinates_particle
!         end if
!     else 
!         write_b_to_screen => do_not_write
!         write_i_to_screen => do_not_write
!     end if
!     if (use_datafile) then
!         if (use_elements_output) then
!             write_b_to_general => do_not_write
!             write_i_to_general => write_elements
!         else
!             write_b_to_general => write_coordinates_boulders
!             write_i_to_general => write_coordinates_particle
!         end if
!     else 
!         write_b_to_general => do_not_write
!         write_i_to_general => do_not_write
!     end if
!     if (use_multiple_outputs) then
!         write_b_to_individual => write_coordinates_boulders
!         if (use_elements_output) then
!             write_i_to_individual => write_elements
!         else
!             write_i_to_individual => write_coordinates_particle
!         end if
!     else 
!         write_b_to_individual => do_not_write
!         write_i_to_individual => do_not_write
!     end if
!     !! Definimos vector para elementos
!     if (use_elements) then
!         get_elements_i => calculate_elements_i
!     else 
!         get_elements_i => do_nothing_i
!     end if
!     !! Definimos vector para chaos
!     if (use_chaos) then
!         get_chaos_i => calculate_chaos_i
!     else 
!         get_chaos_i => do_nothing_i
!     end if
!     !!! Escribir chaos en cada output?
!     if (use_update_chaos) then
!         flush_chaos => write_chaos
!     else
!         flush_chaos => do_nothing_i
!     end if
!     !! Escribir salida en cada output?
!     if (use_flush_output) then
!         flush_output => flush_to_file
!     else
!         flush_output => do_nothing_i
!     end if


!     ! Set parameters new to the derivative
!     parameters_der = dydt(initial_time, parameters_arr)
!     parameters_arr_new = parameters_der

!     !!! Get particles accelerations (just for the output)
!     do i = 1, Nparticles
!         aux_integer = first_particle + 4 * (i - 1)
!         particles_acc(i,1:2) = parameters_der(aux_integer+2 : aux_integer+3)
!     end do


!     !! Initial conditions Output
!     !$OMP PARALLEL DEFAULT(SHARED) &
!     !$OMP PRIVATE(i)
!     !$OMP DO
!     do i = 0, Nboulders
!         call write_b_to_individual(i, 200+i)
!     end do 
!     !$OMP END DO NOWAIT
!     !$OMP DO
!     do i = 1, Nparticles
!         call write_i_to_individual(i, 200+i+Nboulders)
!     end do 
!     !$OMP END DO NOWAIT
!     !$OMP SECTIONS
!     !$OMP SECTION
!     do i = 0, Nboulders
!         call write_b_to_screen(i, 6)
!     end do
!     do i = 1, Nparticles
!         call write_i_to_screen(i, 6)
!     end do
!     !$OMP SECTION
!     do i = 0, Nboulders
!         call write_b_to_general(i, 20)
!     end do
!     do i = 1, Nparticles
!         call write_i_to_general(i, 20)
!     end do
!     !$OMP END SECTIONS
!     !$OMP END PARALLEL
    
!     !!!!!! MAIN LOOP INTEGRATION !!!!!!!
!     ! MAIN LOOP
!     tom_index_number = 2 !!!! Inicializamos en 2 porque el primer checkpoint es el IC (t0)
!     j = 2! From 2 because 1 is the IC (t0) !! +1 por si hay HardExit en el último
!     main_loop: do while (.True.)

!         ! Check for colissions/escapes
!         staying_particles = 0
!         discarded_particles = 0
!         !!! Pre-set merges
!         use_merge = use_merge .and. any(particles_mass(1:Nactive) > cero)
!         merged_particles = 0
!         mass_to_merge = cero
!         angular_momentum_to_merge = cero
!         !! Calculate distances to the asteroid
!         !$OMP PARALLEL DEFAULT(SHARED) &
!         !$OMP PRIVATE(i)
!         !$OMP DO REDUCTION(+:mass_to_merge,angular_momentum_to_merge) SCHEDULE (STATIC)
!         do i = 1, Nactive
!             if ((particles_dist(i) < min_distance) .or. (particles_hexit(i) .eq. 1)) then
!                 if (use_screen) then
!                     write (*,f125131) " Colisión de la partícula ", i, "(", particles_index(i), ") en t = ", &
!                     & time / unit_time, "[días], y r = ", particles_dist(i) / unit_dist, "[km]"
!                     write (*,*) ACHAR(5)
!                 end if
!                 particles_outcome(i) = 1
!                 !$OMP CRITICAL (discard)
!                 discarded_particles = discarded_particles + 1
!                 ij_to_swap(discarded_particles,1) = i
!                 !$OMP END CRITICAL (discard)
!                 call resolve_merge(i, mass_to_merge, angular_momentum_to_merge)
!             else if ((particles_dist(i) > max_distance) .or. (particles_hexit(i) .eq. 2)) then
!                 if (use_screen) then
!                     write (*,f125131) " Escape de la partícula ", i, "(", particles_index(i), ") en t = ", &
!                     & time / unit_time, "[días], y r = ", particles_dist(i) / unit_dist, "[km]"
!                     write (*,*) ACHAR(5)
!                 end if
!                 particles_outcome(i) = 2
!                 !$OMP CRITICAL (discard)
!                 discarded_particles = discarded_particles + 1
!                 ij_to_swap(discarded_particles,1) = i
!                 !$OMP END CRITICAL (discard)
!             end if
!             if (use_chaos) particles_times(i) = time ! Update particle times
!         end do
!         !$OMP END DO
!         !$OMP END PARALLEL

!         !! Staying particles
!         staying_particles = Nactive - discarded_particles
!         !! Discard particles
!         if (discarded_particles > 0) then
!             merged_particles = count(particles_outcome(:Nactive) .eq. 1)
!             if (use_merge .and. (merged_particles > 0)) then
!                 call merge_into_asteroid(mass_to_merge, angular_momentum_to_merge) ! Merge particles into asteroid
!                 if (use_screen) then
!                     write(*,f12) " Merged", merged_particles, "particles into the asteroid"
!                     write(*,f13) "  New mass :", asteroid_mass / unit_mass, "[kg]"
!                     write(*,f13) "  New omega:", asteroid_omega * unit_time, "[rad/día]"
!                     write(*,*) ACHAR(5)
!                 end if                
!             end if
!             call quicksort_int(ij_to_swap(1:discarded_particles,1), 1, discarded_particles) ! Sort particles to discard
!             if (staying_particles == 0) then ! All particles are out
!                 if (use_screen) then
!                     do i = Nactive, 1, -1
!                         write (*,f12) "  - Eliminando partícula ", i, "(", particles_index(i), ")"
!                         write (*,*) ACHAR(5)
!                     end do
!                 end if
!                 Nactive = 0 ! All particles are out
!             else
!                 aux_integer = 1
!                 !$OMP PARALLEL DEFAULT(SHARED) &
!                 !$OMP PRIVATE(i)
!                 !$OMP DO SCHEDULE (STATIC)
!                 do i = 1, discarded_particles
!                     if (particles_outcome(Nactive + 1 - i) .ne. 0) then
!                         ij_to_swap(discarded_particles + 1 - i, 2) = Nactive + 1 - i
!                     else
!                         !$OMP CRITICAL (swap)
!                         ij_to_swap(aux_integer, 2) = Nactive + 1 - i
!                         aux_integer = aux_integer + 1
!                         !$OMP END CRITICAL (swap)
!                     end if
!                 end do
!                 !$OMP END DO
!                 !$OMP BARRIER
!                 !$OMP DO SCHEDULE (STATIC)
!                 do i = 1, discarded_particles
!                     if (use_screen) then
!                         write (*,f12) "  - Eliminando partícula ", ij_to_swap(i,1), "(", particles_index(ij_to_swap(i,1)), ")"
!                         write (*,*) ACHAR(5)
!                     end if
!                     if (ij_to_swap(i,1) .le. (Nactive - discarded_particles)) then
!                         call swap_particles(ij_to_swap(i,1), ij_to_swap(i,2), use_chaos)
!                     end if
!                 end do
!                 !$OMP END DO
!                 !$OMP END PARALLEL
!                 Nactive = Nactive - discarded_particles
                
!                 if (use_screen) then
!                     write (*,f12) " Quedan ", Nactive, " partículas activas."
!                     write (*,*) ACHAR(5)
!                 end if
                
!                 ! Reset sorted index
!                 do i = 1, Nactive
!                     sorted_particles_index(i) = i
!                 end do
!                 call quickargsort_int(particles_index(1:Nactive), sorted_particles_index(1:Nactive), 1, Nactive) ! Get the sorted index
!                 last_particle = first_particle + 4 * Nactive - 1 ! Redefine last_particle

!                 if (use_parallel) then
!                     aux_integer = MIN(MAX(Nactive, 3), my_threads)
!                     if (aux_integer .ne. my_threads) then
!                         my_threads = aux_integer
!                         !$ call OMP_SET_NUM_THREADS(my_threads)
!                         if (use_screen) then
!                             write (*,f12) " WARNING: Se reducen a ", my_threads, " hilos para la integración."
!                             write (*,*) ACHAR(5)
!                         end if
!                     end if
!                 end if
!             end if
!         end if
!         particles_hexit = 0 ! Resetear HardExit
        
        
!         ! Check if all done
!         !! Time end
!         if (j == checkpoint_number + 1) then
!             if (use_screen) then
!                 write (*,*) ACHAR(5) 
!                 write (*,f13) "Finalizó la integración en t = ", time / unit_time, "[días]"
!             end if
!             exit main_loop
!         end if
!         !! Particles left
!         if (Nactive == 0) then
!             if (use_screen) then
!                 write (*,*) ACHAR(5) 
!                 write (*,*) " No quedan partículas activas."
!                 write (*,*) ACHAR(5) 
!                 write (*,f13) "Finalizó la integración en t = ", time / unit_time, "[días]"
!             end if
!             exit main_loop
!         end if

        
!         ! Update dt
!         timestep = checkpoint_times(j) - time
!         !! Check TOM
!         if (checkpoint_is_tom(j) .and. (tom_index_number <= tom_total_number)) then
!             if (allocated(tom_deltaomega)) then
!                 asteroid_omega = asteroid_omega + tom_deltaomega(tom_index_number)
!                 asteroid_omega2 = asteroid_omega * asteroid_omega
!                 asteroid_rotational_period = twopi / asteroid_omega
!                 if (use_explicit_method) asteroid_theta_correction = asteroid_theta - asteroid_omega * (time - initial_time)
!                 if (use_version_1) then
!                     do i = 0, Nboulders
!                         aux_integer = i * equation_size
!                         parameters_arr(aux_integer+3) = - pos_ast_arr(i,2) * asteroid_omega
!                         parameters_arr(aux_integer+4) = pos_ast_arr(i,1) * asteroid_omega
!                     end do
!                 else
!                     parameters_arr(2) = asteroid_omega
!                 end if
!             end if
!             if (allocated(tom_deltamass)) then
!                 tom_mass_growth_param = uno + (tom_deltamass(tom_index_number) / asteroid_mass)
!                 mass_primary = mass_primary * tom_mass_growth_param
!                 mass_ast_arr = mass_ast_arr * tom_mass_growth_param
!                 asteroid_mass = asteroid_mass * tom_mass_growth_param
!                 asteroid_inertia = asteroid_inertia * tom_mass_growth_param
!                 Gasteroid_mass = G * asteroid_mass
!             end if
!             tom_index_number = tom_index_number + 1
!             if (use_screen .and. (allocated(tom_deltaomega) .or. allocated(tom_deltamass))) then
!                 write (*,*) ACHAR(5)
!                 write (*,f13) "Se actualizó Omega | Masa según archivo TOM, en t = ", time / unit_time, "[días]"
!                 write (*,*) ACHAR(5)
!             end if
!         end if

!         !!! Execute an integration method (uncomment/edit one of these)
!         ! call integ_caller (time, parameters_arr(:last_particle), adaptive_timestep, dydt, &
!         !     & Ralston4, timestep, parameters_arr_new(:last_particle), check_continue_ptr)
!         ! call rk_half_step_caller (time, parameters_arr(:last_particle), adaptive_timestep, dydt, &
!         !     & Runge_Kutta5, 5, error_tolerance, learning_rate, min_timestep, timestep, parameters_arr_new(:last_particle), check_continue_ptr)
!         ! call embedded_caller (time, parameters_arr(:last_particle), adaptive_timestep, dydt, Dormand_Prince8_7, &
!         !    & error_tolerance, learning_rate, min_timestep, timestep, parameters_arr_new(:last_particle), check_continue_ptr)
!         call BStoer_caller (time, parameters_arr(:last_particle), adaptive_timestep, dydt, &
!             & error_tolerance, timestep, parameters_arr_new(:last_particle), check_continue_ptr)
!         ! call BStoer_caller2 (time, parameters_arr(:last_particle), adaptive_timestep, dydt, &
!         !     & error_tolerance, timestep, parameters_arr_new(:last_particle), check_continue_ptr)
!         ! call leapfrog_caller (time, parameters_arr(:last_particle), adaptive_timestep, dydt, &
!         !     & leapfrof_KDK, error_tolerance, learning_rate, min_timestep, timestep, parameters_arr_new(:last_particle), check_continue_ptr)
    
!         ! Check if it might be hard_exit
!         if (particles_hexit(0) .ne. 0) then
!             !! If so, the dt used is in dt_adap
!             if (is_premature_exit) timestep = adaptive_timestep
!             !! Check if premature_exit: If exit at less than 1% of finishing timestep
!             if (timestep > cero) then
!                 is_premature_exit = (checkpoint_times(j) - (time+timestep))/timestep > 0.01d0
!             else
!                 is_premature_exit = checkpoint_times(j) - time < tini
!             end if
!         end if

!         ! Update parameters
!         time = time + timestep
!         parameters_der = dydt(time, parameters_arr_new)
!         parameters_arr = parameters_arr_new

!         ! Asteroid and boulders
!         if (use_boulders) then
!             ! Con boulders
!             if (use_explicit_method) then
!                 ! Constantes de Asteroid: pos, vel, omega, inertia, masa
!                 asteroid_theta = asteroid_omega * (time - initial_time) + asteroid_theta_correction!! No lo cambio ahora porque está en explicit_v2
!                 do i = 0, Nboulders
!                     pos_ast_arr(i,1) = cos(asteroid_theta + theta_ast_arr(i)) * dist_ast_arr(i)
!                     pos_ast_arr(i,2) = sin(asteroid_theta + theta_ast_arr(i)) * dist_ast_arr(i)
!                 end do
!                 vel_ast_arr(0:,1) = - asteroid_omega * pos_ast_arr(0:,2)
!                 vel_ast_arr(0:,2) =   asteroid_omega * pos_ast_arr(0:,1)
!                 acc_ast_arr = - asteroid_omega2 * pos_ast_arr
!             else
!                 if (use_version_1) then
!                     ! Constantes de Asteroid (hasta ahora): masa
!                     !! Tendremos que obtener Asteroid (center of mass) properties (pos, vel, omega, ...)
!                     !!! Implicit V1
!                     do i = 0, Nboulders
!                         aux_integer = i * equation_size
!                         pos_ast_arr(i,:) = parameters_arr(aux_integer+1 : aux_integer+2)
!                         vel_ast_arr(i,:) = parameters_arr(aux_integer+3 : aux_integer+4)
!                         acc_ast_arr(i,:) = parameters_der(aux_integer+3 : aux_integer+4)
!                     end do
!                     call get_asteroid_from_boulders(mass_ast_arr, pos_ast_arr, vel_ast_arr, &
!                                                     & asteroid_pos, asteroid_vel, asteroid_omega, asteroid_inertia)
!                     do i = 0, Nboulders
!                         dist_ast_arr(i) = sqrt(sum((pos_ast_arr(i,:) - asteroid_pos)**2))
!                         ! acc_ast_arr(i,:) = - (pos_ast_arr(i,:) - asteroid_pos) * asteroid_omega**2
!                     end do
!                     ! Definimos theta como la variación del ángulo de m0 respecto del asteroide, respecto a la condición inicial
!                     asteroid_theta = mod(atan2(pos_ast_arr(0,2) - asteroid_pos(2), &
!                                             & pos_ast_arr(0,1) - asteroid_pos(1)) - &
!                                             & theta_ast_arr(0), twopi)
!                 else 
!                     ! Constantes de Asteroid (hasta ahora): masa, inertia
!                     !! Las Asteroid (center of mass) properties están servidas
!                     !!! Implicit V2
!                     parameters_arr(1) = mod(parameters_arr(1), twopi)
!                     asteroid_theta = parameters_arr(1)
!                     asteroid_omega = parameters_arr(2)
!                     asteroid_pos = parameters_arr(3:4)
!                     asteroid_vel = parameters_arr(5:6)
!                     pos_ast_arr(0:,1) = cos(asteroid_theta + theta_ast_arr(0:)) * dist_ast_arr(0:)
!                     pos_ast_arr(0:,2) = sin(asteroid_theta + theta_ast_arr(0:)) * dist_ast_arr(0:)
!                     vel_ast_arr(0:,1) = - asteroid_omega * pos_ast_arr(0:,2)
!                     vel_ast_arr(0:,2) =   asteroid_omega * pos_ast_arr(0:,1)
!                     do i = 0, Nboulders
!                         pos_ast_arr(i,:) = asteroid_pos + pos_ast_arr(i,:)
!                         vel_ast_arr(i,:) = asteroid_vel + vel_ast_arr(i,:)
!                     end do
!                     acc_ast_arr = - asteroid_omega2 * pos_ast_arr
!                 end if
!             end if

!             ! Update asteroid properties
!             asteroid_omega2 = asteroid_omega * asteroid_omega
!             asteroid_a_corot = (G * asteroid_mass / asteroid_omega2)**(1/3.)
!             asteroid_rotational_period = twopi / asteroid_omega
        
!         else 
!             ! Sin boulders, casi todo constante
!             asteroid_theta = mod(asteroid_omega * (time - initial_time), twopi)
!             asteroid_pos = parameters_arr(1:2)
!             asteroid_vel = parameters_arr(3:4)
!             pos_ast_arr(0,:) = asteroid_pos
!             vel_ast_arr(0,:) = asteroid_vel
!         end if

        
!         !! Particles
!         !$OMP PARALLEL DEFAULT(SHARED) &
!         !$OMP PRIVATE(i,aux_integer)
!         !$OMP DO SCHEDULE (STATIC)
!         do i = 1, Nactive
!             aux_integer = first_particle + 4 * (i - 1)
!             !!! Coordinates
!             particles_coord(i,1:2) = parameters_arr(aux_integer   : aux_integer+1)
!             particles_coord(i,3:4) = parameters_arr(aux_integer+2 : aux_integer+3)
!             particles_acc(i,1:2) = parameters_der(aux_integer+2 : aux_integer+3)
!             call get_elements_i(i) !!! Elements
!             call get_chaos_i(i) !!! Chaos
!             particles_dist(i) = sqrt(sum((particles_coord(i,1:2) - asteroid_pos)**2)) ! sqrt(x^2 + y^2)
!         end do
!         !$OMP END DO
!         !$OMP END PARALLEL
        
!         ! Output
!         if ((checkpoint_is_output(j)) .and. (.not. is_premature_exit)) then
!             !$OMP PARALLEL DEFAULT(SHARED) &
!             !$OMP PRIVATE(i,aux_integer)
!             !$OMP DO
!             do i = 0, Nboulders
!                 call write_b_to_individual(i, 200+i)
!             end do 
!             !$OMP END DO NOWAIT
!             !$OMP DO
!             do i = 1, Nactive
!                 call write_i_to_individual(i, 200+i+Nboulders)
!             end do 
!             !$OMP END DO NOWAIT
!             !$OMP SECTIONS
!             !$OMP SECTION
!             do i = 0, Nboulders
!                 call write_b_to_screen(i, 6)
!             end do
!             do i = 1, Nactive
!                 call write_i_to_screen(i, 6)
!             end do
!             !$OMP SECTION
!             do i = 0, Nboulders
!                 call write_b_to_general(i, 20)
!             end do
!             do i = 1, Nactive
!                 call write_i_to_general(i, 20)
!             end do
!             call flush_output(20)
!             !$OMP SECTION
!             call flush_chaos(40) ! Update chaos
!             !$OMP END SECTIONS
!             !$OMP END PARALLEL
!         end if

!         if (use_percentage) call percentage(time, final_time)

!         ! Update j; only if not HardExit
!         if (.not. is_premature_exit) j = j + 1
        
!     end do main_loop
    
!     !! Update surviving particles times
!     if (use_chaos) particles_times(1:Nactive) = final_time ! Update particle times

    !! Porcentaje final
    if (sim%use_percentage) call percentage(sim%final_time + uno, sim%final_time)

    !! Cerrar archivo de salida
    if (sim%use_datafile) then
        close (2)
        if (sim%use_screen) then
            write (*,*) ACHAR(10)
            write (*,*) "Se guardó la salida en el archivo: ", trim(sim%datafile)
        end if
    end if

    !! Cerrar archivos individuales
    if (sim%use_multiple_outputs) then
        do i = 0, sim%Nactive
            close (200+i)
        end do
        if (sim%use_screen) then
            write (*,*) ACHAR(10)
            write (*,*) "Se guardaron las salidas individuales en: ", trim(sim%multfile) // "_*"
        end if
    end if

!     !! Caos
!     if (use_chaos) then
!         if (use_chaosfile) then
!             !! Guardar archivo de caos
!             do i = 1, Nactive
!                 sorted_particles_index(i) = i
!             end do
!             call quickargsort_int(particles_index, sorted_particles_index, 1, Nparticles)
!             inquire (unit=40, opened=aux_logical)
!             if (.not. aux_logical) then
!                 open (unit=40, file=trim(chaosfile), status='unknown', action='write')
!             else
!                 call fseek(40, 0, 0)       ! move to beginning
!             end if
!             do i = 1, Nparticles
!                 aux_integer = sorted_particles_index(i)
!                 write (40,f2233) particles_index(aux_integer), & ! i
!                 & particles_outcome(aux_integer), & ! bad
!                 & final_time / unit_time, & ! total time to integrate
!                 & asteroid_initial_conditions(10) / (unit_mass * unit_dist * unit_vel), & ! initial (Asteroid): angular momentum
!                 & particles_initial_conditions(i,1) / unit_mass, & ! initial: mass
!                 & particles_initial_conditions(i,2) / unit_dist, & ! initial: a
!                 & particles_initial_conditions(i,3), & ! initial: e
!                 & particles_initial_conditions(i,4) / radian, & ! initial: M
!                 & particles_initial_conditions(i,5) / radian, & ! initial: omega
!                 & particles_initial_conditions(i,6), & ! initial: MMR
!                 & sqrt(particles_initial_conditions(i,2) * &
!                   & (uno - particles_initial_conditions(i,3)**2) * &
!                   & G * (particles_initial_conditions(i,1) + asteroid_initial_conditions(1))) / &
!                   & (unit_dist * unit_vel), & ! initial: angular momentum per unit mass
!                 & particles_times(aux_integer) / unit_time, & ! surviving time
!                 & particles_elem(aux_integer,1) / unit_dist, particles_elem(aux_integer,2), & ! final: a, e
!                 & particles_elem(aux_integer,3) / radian, particles_elem(aux_integer,4) / radian, & ! final: M, omega
!                 & particles_MMR(aux_integer), & ! final: MMR
!                 & sqrt(particles_elem(aux_integer,1) * &
!                   & (uno - particles_elem(aux_integer,2)**2) * &
!                   & G * (particles_mass(aux_integer) + asteroid_mass)) / &
!                   & (unit_dist * unit_vel), & ! final: angular momentum per unit mass
!                 & particles_min_a(aux_integer) / unit_dist, particles_max_a(aux_integer) / unit_dist, & ! a_min, a_max
!                 & particles_min_e(aux_integer), particles_max_e(aux_integer), & ! e_min, e_max
!                 & (particles_max_a(aux_integer) - particles_min_a(aux_integer)) / unit_dist, & ! Delta a
!                 & (particles_max_e(aux_integer) - particles_min_e(aux_integer)) ! Delta e
!             end do
!             close (40)
!             !! Mensaje
!             if (use_screen) then 
!                 write (*,*) ACHAR(10)
!                 write (*,*) "Se guardó el archivo de caos en: ", trim(chaosfile)
!             end if
!         end if
!         if (use_datascreen) then
!             write (*,*) ACHAR(10)
!             if (use_single_particle) then            
!                 write (*,*) "Caos de la partícula simple [i, bad, tmax, MMR_ini, MMR_fin, a_min, a_max, Delta a, &
!                 & e_min, e_max, Delta_e]:"
!             else 
!                 write (*,*) "Caos de las partículas [i, bad, tmax, MMR_ini, MMR_fin, a_min, a_max, Delta a, &
!                 & e_min, e_max, Delta_e]:"
!             end if
!             do i = 1, Nparticles
!                 aux_integer = sorted_particles_index(i)
!                 write (*,f2233) particles_index(aux_integer), & ! i
!                 & particles_outcome(aux_integer), & ! bad
!                 & particles_times(aux_integer), & ! surviving time
!                 & particles_initial_conditions(aux_integer,6), & ! initial: MMR
!                 & particles_MMR(aux_integer), & ! final: MMR
!                 & particles_min_a(aux_integer) / unit_dist, particles_max_a(aux_integer) / unit_dist, & ! a_min, a_max
!                 & particles_min_e(aux_integer), particles_max_e(aux_integer), & ! e_min, e_max
!                 & (particles_max_a(aux_integer) - particles_min_a(aux_integer)) / unit_dist, & ! Delta a
!                 & (particles_max_e(aux_integer) - particles_min_e(aux_integer)) ! Delta e
!             end do
!         end if
!     end if

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!! LIBERACION MEMORIA!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (sim%use_screen) then
        write (*,*) ACHAR(10)
        write (*,*) "Liberando memoria"
        write (*,*) ACHAR(5)
    end if
    !!!!!!!!! TIMES !!!!!!!!
    call free_times_arrays()

    !!!!!!!!! PARAMS !!!!!!!!
    call free_parameters_arays()

    !!!!!!!!! BODIES !!!!!!!!
    call free_asteroid(asteroid)
    call free_system(system)
    call free_system(initial_system)

!     !!!!!!!! Asteroide  !!!!!!!
!     call free_asteroid_arrays()
!     !!!!!!!! PARTICULAS !!!!!!!
!     call free_particles()
!     !!!!!!!! PARAMETROS !!!!!!!
!     deallocate(parameters_arr)
!     deallocate(parameters_arr_new)
!     deallocate(parameters_der)
!     !!!!!!!!! INITIALS !!!!!!!!
!     call free_initial_conditions()
!     !!!!!!!!! POINTERS !!!!!!!!
!     call nullify_pointers()
!     nullify(dydt)
!     nullify(domegadt)
!     nullify(dmassdt)
    

    if (sim%use_screen) then 
        write (*,*) ACHAR(5)
        write (*,*) "Fin del programa"
        write (*,*) ACHAR(5)
    end if
    
end program main

subroutine create_map(system,ngx,ngy,xmin,xmax,ymin,ymax,map_file)
    use constants, only: G, cero
    use bodies, only: system_st, get_acc_and_pot_xy
    implicit none
    type(system_st), intent(in) :: system
    integer(kind=4), intent(in) :: ngx, ngy
    real(kind=8), intent(in)    :: xmin, xmax, ymin, ymax
    character(len=*), intent(in) :: map_file
    real(kind=8) :: pot(ngx,ngy), acc(ngx,ngy,2)
    real(kind=8) :: rb(2)
    real(kind=8) :: pot_at_R
    integer(kind=4) :: i, j

    !!! Check
    if ((ngx < 2) .or. (ngy < 2)) then
        write (*,*) "ERROR: Check that (ngx > 2) and (ngy > 2)"
        stop 1
    end if
    if ((xmin >= xmax) .or. (ymin >= ymax)) then
        write (*,*) "ERROR: Check that (xmin < xmax) and (ymin < ymax)"
        stop 1
    end if
    pot = cero
    acc = cero

    pot_at_R = - system%asteroid%boulders(0)%mass / system%asteroid%boulders(0)%radius  ! G at bottom

    !$OMP PARALLEL DEFAULT(SHARED) &
    !$OMP PRIVATE(i,j,k,rb,dist,dummy,dummy2,dumint,cero2)
    !$OMP DO SCHEDULE (STATIC)
    do i = 1, ngx
        do j = 1, ngy
            rb(1) = xmin + i * (xmax - xmin) / ngx
            rb(2) = ymin + j * (ymax - ymin) / ngy
            call get_acc_and_pot_xy(system, rb, acc(i,j,:), pot(i,j), pot_at_R)
        end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

    pot = pot * G  ! Convertimos a unidades correctas
    acc = acc * G  ! Convertimos a unidades correctas
    
    write (*,*) "   Potential calculated. Writing file..."
    open (unit=50, file=trim(map_file), status='replace', action='write')
    do i = 1, ngx
        do j = 1, ngy
            write (50,*) xmin + i * (xmax - xmin) / ngx, ymin + j * (ymax - ymin) / ngy, &
            pot(i,j), acc(i,j,:)
        end do
    end do
    close (50)
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

