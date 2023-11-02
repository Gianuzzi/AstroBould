!!! Versión 5.0
!!! Compilación:
!!! gfortran -O2 const.F90 run.F90 parameters.F90 integrators.F90 stokes.F90 gravity.F90 derivates.F90 main.F90 -o main

module init
    use run
    use integrators
    use derivates ! Incluye const, gravity, extern

    implicit none

    integer(kind=4) :: i, j, nin, nsim, bad, ineqs
    integer(kind=4) :: ngx, ngy 
    real(kind=8)    :: xmin, xmax, ymin, ymax 
    character(30)   :: chn, cha, che, chM, chw, chR, datafile, chaosfile, infofile, map_file
    integer(kind=4) :: auxin
    character(30)   :: auxch
    character(1)    :: auxch1
    character(2)    :: auxch2
    logical         :: auxlo, is_number
    logical         :: screen, perc, datas, eleout
    logical         :: map_pot, infoo, datao, chaos
    logical         :: exact
    real(kind=8)    :: ang_mom

    contains

        subroutine init_vars()
            !!! Inicializamos (Variables default)
            nsim = 0    ; bad = 0
            m    = cero ; mu  = cero ; GM = cero ; mcm = cero ; Gm = cero
            mucm = cero
            radius = cero
            rib = cero ; ria = cero ; rcm = cero
            vib = cero ; via = cero ; vcm = cero
            aib = cero ; aia = cero ; acm = cero
            rb = cero  ; vb = cero  ; ab = cero
            ra = cero  ; va = cero  ; aa = cero
            Prot = uno ; R0 = uno   ; theta_a = cero ; lambda = cero
            r_b = cero ; theta_b = cero
            omega = cero ; wk = cero   ; lambda2 = cero ; omega2 = cero 
            yb = cero ; ybnew = cero ; ya = cero ; yanew = cero
            xc = cero ; xe = cero  ; xc0 = cero  ; xe0 = cero
            Res0 = cero
            ea = cero ; ee = cero  ; ei = cero   ; eM = cero 
            ew = cero ; eO = cero  ; eR = cero
            da = cero ; de = cero  ; amax = cero ; amin = inf
            emax = cero ; emin = inf
            tau_a = inf  ; tau_e = inf ; t_stokes = cero
            t0 = cero    ; tf = cero   ; dt_min = cero
            dt_out = cero ; n_points = 1
            logsp = .FALSE.
            beta = 0.9d0 ; e_tol = 1.d-12
            ngx = 500   ; ngy = 500
            xmin = -300 ; xmax = 300 ; ymin = -300 ; ymax = 300
            screen  = .False. ; perc = .False. !Pantalla, porcentaje
            map_pot = .False. ! Mapa de potencial
            infoo   = .False. ! Archivo de info
            datao   = .False. ! Archivo de datos
            chaos   = .True.  ! Archivo de caos
            exact   = .True.  ! Método exacto
            datas   = .False. ! Salida de datos en pantalla
            eleout  = .False. ! Salida en elementos orbitales 
        end subroutine init_vars

end module init

program main
    use init
    use run
    use integrators
    use derivates ! Incluye const, gravity, extern

    implicit none

    call init_vars() ! Inicializamos variables
 
    !!! Definir parámetros:
    m(0)       = 6.3d18       ! Masa del cuerpo 0 [kg]
    mu(1)      = 1.d-7        ! Cociente de masas
    mu(2)      = 1.d-6        ! Cociente de masas
    mu(3)      = 1.d-4        ! Cociente de masas
    ! mu(4)      = 1.d-7        ! Cociente de masas
    theta_a(1) = cero         ! Ángulo de fase respecto del origen angular [rad]
    theta_a(2) = pi           ! Ángulo de fase [rad]
    theta_a(3) = pi * uno3    ! Ángulo de fase [rad]
    ! theta_a(4) = 250 * rad    ! Ángulo de fase [rad]
    radius(0)  = 129.d0       ! Radio del cuerpo 0 [km]
    radius(1)  = 2.5d0        ! Radio del boulder 1 [km]
    radius(2)  = 1.5d0        ! Radio del boulder 2 [km]
    radius(3)  = 3.5d0        ! Radio del boulder 3 [km]
    ! radius(4)  = 2.5d0        ! Radio del boulder 3 [km]
    ! Prot       = 7.004d0/24.d0  ! Periodo de rotación de los cuerpos [days]
    lambda     = 0.471d0      ! Cociente spin/wk

    !!!! Stokes
    tau_a = inf           ! [Prot]
    tau_e = tau_a / 1.d2  ! [Prot]
    t_stokes = cero       ! [Prot] Tiempo que actua stokes

    !!! Parámetros corrida
    t0       = cero            ! Initial time [Prot]
    tf       = 5.d3            ! Final time [Prot]
    dt_min   = 1.d-3           ! Min timestep [Prot] ! Almost unused
    logsp    = .FALSE.         ! LogSpaced outputs
    n_points = 100             ! Number of outputs (if logsp=.TRUE. or dt_out=0)
    dt_out   = 1.d-1            ! Output timestep [Prot] (if logsp = .False.)
    beta     = 0.85d0          ! [For adaptive step integrators] Learning rate
    e_tol    = 1.d-11          ! [For adaptive step integrators] Relative error
    rmax     = 1.d2*radius(0)  ! Max distance before escape [km]
    

    !! Default
    !!! Output ()"" or "no", if not used)
    infofile  = ""
    datafile  = ""
    chaosfile = "chaos.dat"
    map_file  = ""

    !!!! Particle elements
    ea = cero                  ! Element a of the particle
    ee = cero !0.1d0           ! ecc
    eM = cero                  ! Mean anomaly
    ew = cero                  ! Pericenter argument
    eR = cero !1.001d0 !3.3d0  ! resonancia nominal correspondiente

    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !Leemos de la línea de comandos
    nin = command_argument_count()
    !!!!! PARTÍCULA (en caso de entrada por terminal)
    auxlo = .False.
    auxin = 0
    if (nin > 0) then
        do i = 1, nin
            if (auxin /= 0) then
                auxin = auxin - 1
                cycle
            end if
            call get_command_argument(i, auxch)
            select case (trim(auxch))
            case ("--nodata")
                datao    = .False.
                datafile = ""
            case ("-datafile")
                datao    = .True.
                call get_command_argument(i+1, datafile)
                auxin = 1
            case ("--noinfo")
                infoo    = .False.
                infofile = ""
            case ("-infofile")
                infoo    = .True.
                call get_command_argument(i+1, infofile)
                auxin = 1
            case ("--nochaos")
                chaos    = .False.
                chaosfile = ""
            case ("-chaosfile")
                chaos    = .True.
                call get_command_argument(i+1, chaosfile)
                auxin = 1
            case ("--screen")
                screen = .True.
            case ("--noscreen")
                screen = .False.
            case ("--perc")
                perc = .True.
            case ( "--noperc")
                perc = .False.
            case ("--datascr")
                datas = .True.
            case ( "--nodatascr")
                datas = .False.
            case ("--nomap")
                map_pot  = .False.
                map_file = ""
            case ("-mapfile")
                map_pot  = .True.
                call get_command_argument(i+1, map_file)
                auxin = 1
            case ("--noexact")
                exact = .False.
            case ("--exact")
                exact = .True.
            case ("--elem")
                eleout = .True.
            case ("--noelem")
                eleout = .True.
            case ("--help")
                call get_command_argument(0, auxch)
                write(*,*) "Uso: " // trim(auxch) // " [nsim] [ea] [ee] [eM] [ew] [eR]"
                write(*,*) "    nsim: Número de simulación"
                write(*,*) "    ea  : Elemento a de la partícula"
                write(*,*) "    ee  : Elemento e de la partícula"
                write(*,*) "    eM  : Elemento M de la partícula"
                write(*,*) "    ew  : Elemento w de la partícula"
                write(*,*) "    eR  : Elemento R de la partícula"
                write(*,*) "    --nodata   : No guardar datos"
                write(*,*) "    -datafile  : Guardar datos en el archivo que sigue"
                write(*,*) "    --noinfo   : No guardar información"
                write(*,*) "    -infofile  : Guardar información en el archivo que sigue"
                write(*,*) "    --nochaos  : No guardar caos"
                write(*,*) "    -chaosfile : Guardar caos en el archivo que sigue"
                write(*,*) "    --screen   : Imprimir información en pantalla"
                write(*,*) "    --noscreen : No imprimir en pantalla"
                write(*,*) "    --perc     : Imprimir porcentaje de integración"
                write(*,*) "    --noperc   : No imprimir porcentaje de integración"
                write(*,*) "    --datascr  : Imprimir datos en pantalla"
                write(*,*) "    --nodatascr: No imprimir datos en pantalla"
                write(*,*) "    --nomap    : No guardar mapas de potencial"
                write(*,*) "    -mapfile   : Guardar mapas de potencial en el archivo que sigue"
                write(*,*) "    --noexact  : No usar método exacto"
                write(*,*) "    --exact    : Usar método exacto"
                write(*,*) "    --elem     : Imprimir elementos orbitales (solo partícula)"
                write(*,*) "    --noelem   : Imprimir coordenadas baricéntricas"
                stop 0
            case default
                is_number = .False.
                call get_command_argument(i, auxch)
                do j = 0, 9
                    if (auxch(1:1) == char(48 + j)) then !check if it's a number
                        is_number = .True.
                        exit
                    end if
                end do
                if (.not. is_number) then ! No es un número
                    write(*, '(A, I0, A, A)') "Error: No se reconoce el argumento ", i, ": ", trim(auxch)
                    call get_command_argument(0, auxch)
                    write(*,*) "Para ayuda, ejecute: ", trim(auxch), " --help"
                    write(*,*) "Saliendo."
                    stop 1
                end if
                if (.not. auxlo)  then! No leí los parámetros aún
                    if (nin < i+4) then
                        write(*,*) "Error: Error al intentar leer 'nsim'."
                        write(*,*) "Se requieren al menos 5 argumentos más de línea de comandos"
                        write(*,*) "Saliendo."
                        stop 1
                    else ! Leo los argumentos numéricos
                        call get_command_argument(i, chn)
                        read(chn,*) nsim
                        call get_command_argument(i+1, cha)
                        read(cha,*) ea
                        call get_command_argument(i+2, che)
                        read(che,*) ee
                        call get_command_argument(i+3, chM)
                        read(chM,*) eM
                        call get_command_argument(i+4, chw)
                        read(chw,*) ew
                        auxin = 4
                        auxlo = .True.
                        eR = cero
                    end if
                else ! Ya leí los numéricos. Falta leer eR
                    call get_command_argument(i, chR)
                    read(chR,*) eR
                end if
            end select
        end do
    else
        print*, "WARNING: Se utilizan los parámetros del código."
    end if


!!!!!!!!!!!!!!!!!!!!!!!!!!    OUTPUT     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if ((trim(infofile) /= "") .and. (trim(infofile) /= "no")) then
        infoo = .True.
    else 
        infoo = .False.
    end if
    if ((trim(datafile) /= "") .and. (trim(datafile) /= "no")) then
        datao = .True.
    else 
        datao = .False.
    end if
    if ((trim(chaosfile) /= "") .and. (trim(chaosfile) /= "no")) then
        chaos = .True.
    else 
        chaos = .False.
    end if
    if ((trim(map_file) /= "") .and. (trim(map_file) /= "no")) then
        map_pot = .True.
    else 
        map_pot = .False.
    end if

    if (.not. any((/screen, datao, chaos, map_pot/))) then ! No tiene sentido hacer nada
        write(*,*)  "EXITING: No se guardará ni imprimirá ninguna salida."
        stop 1
    end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    !!! Calcular variables:
    !!!! Radio
    radius = radius * unit_r ! Unidades según G
    R0 = radius(0)

    !!!! Masas
    do i = 1, Nboul
        m(i) = mu(i) * m(0) ! Masa del boulder i
    end do
    m = m * unit_m ! Unidades según G
    mcm  = sum(m)  ! Masa del sistema
    mucm = m / mcm ! Mues respecto a mcm
    Gmi  = G * m
    GM   = G * mcm

    !!!! Rotaciones
    wk = sqrt(GM / R0**3) / unit_t     ! Movimiento medio (kepleriano) que tendrían los boulder
    if (lambda > tini) then    
        omega = wk * lambda            ! Velocidad angular del cuerpo 0 [rad/days]
        Prot  = twopi / omega          ! Periodo de rotación del cuerpo 0 [days]
    else
        omega  = twopi / Prot / unit_t ! Velocidad angular del cuerpo 0 (rad/days)    
        lambda = omega / wk            ! Cociente de velocidades
    end if
    omega2  = omega * omega
    lambda2 = omega2 * R0**3 / GM      ! Coeficiente de Sicardy y Madeira
    a_corot = (GM / omega2)**(1/3.)
    if (screen) then
        write(*,*) "a_corot [km]:", a_corot/unit_r, "R0 [km]:", R0/unit_r
        write(*,*) "omega [rad/days]:", omega*unit_t, "omega_kep [rad/days]:", wk*unit_t, "lambda:", lambda, &
        & "Prot [hs]:", Prot*24*unit_t
    end if

    !!!! Coordenadas
    if (screen) write(*,*) "Masa cuerpo central [kg]:", m(0) / unit_m

    !!!!! Astrocentricas
    do i = 1, Nboul
        ria(i,1) = R0 * cos(theta_a(i)) ! Posición x del boulder i
        ria(i,2) = R0 * sin(theta_a(i)) ! Posición y del boulder i
        via(i,1) = - omega * ria(i,2)   ! Velocidad x del boulder i
        via(i,2) =   omega * ria(i,1)   ! Velocidad y del boulder i
        aia(i,1) = - omega2 * ria(i,1)  ! Aceleración x del boulder i
        aia(i,2) = - omega2 * ria(i,2)  ! Aceleración y del boulder i
    end do
    if (screen) then
        write(*,*) "Astrocentricas:"
        do i = 1, Nboul
            write(*,*) i, ria(i,:)/unit_r, via(i,:)/unit_v, aia(i,:)/unit_a, mu(i), theta_a(i)
        end do
    end if

    !!!!! Centro de masas
    call rcmfromast(ria(1:Nboul,:), via(1:Nboul,:), aia(1:Nboul,:), m(0:Nboul), rcm, vcm, acm)
    theta_acm = atan2(rcm(2), rcm(1))
    if (screen) then
        write(*,*) "Centro de masas:"
        write(*,*) rcm/unit_r, vcm/unit_v, acm/unit_a, mcm/unit_m
    end if

    !!!!! Baricentricas
    rib(0,:) = -rcm
    vib(0,:) = -vcm
    aib(0,:) = -acm
    r_b(0)     = sqrt(rib(0,1)*rib(0,1) + rib(0,2)*rib(0,2))
    theta_b(0) = atan2(rib(0,2), rib(0,1))
    do i = 1, Nboul
        rib(i,:)   = ria(i,:) - rcm
        vib(i,:)   = via(i,:) - vcm
        aib(i,:)   = aia(i,:) - acm
        r_b(i)     = sqrt(rib(i,1)*rib(i,1) + rib(i,2)*rib(i,2))
        theta_b(i) = atan2(rib(i,2), rib(i,1))
    end do
    if (screen) then
        write(*,*) "Baricentricas:"
        do i = 0, Nboul
            write(*,*) i, rib(i,:)/unit_r, vib(i,:)/unit_v, aib(i,:)/unit_a, m(i)/unit_m
        end do
    end if

    !!!!! Momento angular
    call ang_mom_bar(rib, vib, m, radius, omega, ang_mom)
    if (screen) write(*,*) "Momento angular:", ang_mom, "kg km^2 / s"
    

    !!!! MAPAS Potenciales y Aceleraciones
    if (map_pot) then
        if (screen) write(*,*) "Calculando MAPAS potencial..."
        call mapas_pot(Nboul,ngx,ngy,xmin,xmax,ymin,ymax,rib,m,omega,map_file)
        if (screen) write(*,*) "Guardado en el archivo: ", trim(map_file)
    end if
 
    !! Configuramos partícula
    ea = ea * unit_r
    if (eR > tini) then
        ea = eR**(2/3.) * a_corot
    else
        eR = (ea/a_corot)**(1.5)
    end if
    !! Set initial vectors
    xe0  = (/ea, ee, ei, eM, ew, eO/)
    xe   = xe0
    call coord(mcm, ea, ee, ei, eM, ew, eO, xc0)
    call chaosvals(ea, ee, amax, amin, emax, emin)
    xc   = xc0
    Res0 = eR

    rb = xc(1:2)
    vb = xc(4:5)

    !!! INTEGRACIÓN

    !!!! Set Stokes
    tau_a = tau_a*Prot ; tau_e = tau_e*Prot ; t_stokes = t_stokes*Prot
    call Candalpha(tau_a, tau_e, C_stk, a_stk)

    if (screen) then
        write(*,*) "Stokes:"
        write(*,*) "    t_stokes [Prot]:", t_stokes, "C    : ", C_stk,  "alpha: ", a_stk
    end if
    
    if (exact) then
        if ((tau_m < inf) .or. (tau_o < inf)) then
            write(*,*) "Error: No se puede usar el método exacto con tau_m o tau_o finitos."
            if (screen) then
                if (tau_m < inf) write(*,*) "tau_m [units by G]:", tau_m
                if (tau_o < inf) write(*,*) "tau_o [units by G]:", tau_o
            end if
            stop 1
        end if
    end if
    if (screen) write(*,*) "Preparando integración..."
    
    !!!! Yb[aricentric]
    ! call ast2bar(ra, va, aa, rcm, vcm, acm, rb, vb, ab)
    yb(1) = omega
    do i = 0, Nboul
        ineqs = i * neqs
        yb(ineqs+2)           = m(i)
        yb(ineqs+3 : ineqs+4) = rib(i,:)
        yb(ineqs+5 : ineqs+6) = vib(i,:)
    end do
    yb(NP+2)        = m(FP)
    yb(NP+3 : NP+4) = rb
    yb(NP+5 : NP+6) = vb
    
    !!!! Ya[strocentric]
    ya(1)   = omega
    ya(2)   = m(0)
    ya(3)   = radius(0)
    ya(4:6) = cero
    do i = 1, Nboul
        ineqs = i * neqs
        ya(ineqs+2)           = m(i)
        ya(ineqs+3 : ineqs+4) = ria(i,:)
        ya(ineqs+5 : ineqs+6) = via(i,:)
    end do
    call bar2ast(rb, vb, ab, rcm, vcm, acm, ra, va, aa)
    ya(NP+2)        = m(FP)
    ya(NP+3 : NP+4) = ra
    ya(NP+5 : NP+6) = va
    

    !!!! Set times
    t0 = t0 * Prot ; tf = tf * Prot ; dt_out = dt_out * Prot
    call get_t_outs (t0, tf, n_points, dt_out, logsp, t_out) ! Get LOOP checkpoints
    t       = t0                                             ! Init time
    dt_adap = dt_min                                         ! For adaptive step
    dt      = t_out(1) - t0                                  ! This should be == 0
    dt_min  = min(dt_min * Prot, dt_out)

    if (screen) write(*,*) "    t0 [Prot]: ", t0/Prot, " tf [Prot]: ", tf/Prot, " dt_out [Prot]: ", dt_out/Prot, &
        &" dt_min [Prot]: ", dt_min/Prot, "n_points: ", n_points


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  INFO FILE  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (infoo) then
        i = 0
        inquire(file=trim(infofile), exist=auxlo)
        do while (auxlo)
            i = i + 1
            if (i > 10) then
                if (screen) write(*,*) "Error: No se pudo crear el archivo de resumen."
                i = -1
                exit
            end if
            write (auxch2,'(I2.2)') i
            inquire(file=trim(infofile)//"_"//trim(auxch2), exist=auxlo)
        end do
        if (i >= 0) then
            if (i == 0) then
                if (screen) write(*,*) "Se guardará el resumen en el archivo: ", trim(infofile)
                open(unit=1, file=trim(infofile), status="new", action="write")
            else
                if (screen) write(*,*) "Se guardará el resumen en el archivo: ", trim(infofile)//"_"//trim(auxch2)
                open(unit=1, file=trim(infofile)//"_"//trim(auxch2), status="new", action="write")
            end if
            write(1,*) "Nboul,m0,mu,theta,R0,Prot,tau_a,tau_e,t_stokes,tf,dt_out"
            write(1,*) Nboul,m(0),mu,theta_a,R0,Prot,tau_a,tau_e,t_stokes,tf,dt_out
            close(1)
        end if
    end if
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    if (screen) write(*,*) "Integrando..."

    if (datao) open(unit=2, file=trim(datafile), status='unknown', action='write', access="append")
    !!! t, x, y, vx, vy, ax, ay, a, e, w
    
    if (exact) then
        dydt => dydt_bar_ex
    else
        dydt => dydt_bar_im
    end if
    !!! BARYCENTRIC
    if (datas .and. .not. perc) write(*,*) t, ea, ee, ew, eR
    if (datao) then
        if (eleout) then
            write(2,*) t, ea, ee, eM, ew, eR
        else
            do i = 0, Nboul
                write(2,*) i, t, rib(i,:), vib(i,:), aib(i,:), m(i), radius(i)
            end do
            call accbar(m(0:Nboul), rb, rib, ab)   ! Este lo setea a 0 al principio
            call accsto(t, rb, vb, ab) ! Este suma aceleración
            write(2,*) FP, t, rb, vb, ab, cero, cero
        end if
    end if
    ! MAIN LOOP
    do j = 2, n_points ! From 2 because 1 is the IC (t0)
        ra  = rb + rcm
        da0 = sqrt(ra(1)*ra(1) + ra(2)*ra(2))
        if (da0 < R0) then
            if (screen) write(*,*) ACHAR(10) // "Impacto en t = ", t
            bad = 1
            exit
        else if (da0 > rmax) then
            if (screen) write(*,*) ACHAR(10) // "Escape en t = ", t
            bad = 2
            exit
        end if

        ! Update dt
        dt = t_out(j) - t

        !!! Execute an integration method (uncomment/edit one of theese)
        ! call integ_caller (t, yb, dt_adap, dydt, Runge_Kutta4, dt, ybnew)
        ! call rk_half_step_caller (t, yb, dt_adap, dydt, Runge_Kutta5, 5, e_tol, beta, dt_min, dt, ybnew)
        ! call embedded_caller (t, yb, dt_adap, dydt, Dormand_Prince8_7, e_tol, beta, dt_min, dt, ybnew)
        call BStoer_caller (t, yb, dt_adap, dydt, e_tol, dt_min, dt, ybnew)
        ! call BStoer_caller2 (t, yb, dt_adap, dydt, e_tol, dt_min, dt, ybnew)
        
        ! Update parameters
        t  = t + dt
        yb = ybnew

        ybnew = dydt(t, yb)

        ! Asteroid and boulders
        if (exact) then
            do i = 0, Nboul
                ineqs = i * neqs
                rib(i,1) = cos(omega * t + theta_b(i)) * r_b(i)
                rib(i,2) = sin(omega * t + theta_b(i)) * r_b(i)
                vib(i,1)   = - omega * rib(i,2)
                vib(i,2)   =   omega * rib(i,1)
            end do
        else
            omega  = yb(1)
            omega2 = omega * omega
            do i = 0, Nboul
                ineqs = i * neqs
                m(i)      = yb(ineqs+2)
                rib(i,:)  = yb(ineqs+3 : ineqs+4)
                vib(i,:)  = yb(ineqs+5 : ineqs+6)
            end do
            a_corot = (GM / omega2)**(1/3.)
        end if
        aib = -omega2 * rib

        !! Center of mass
        rcm = - rib(0,:)
        vcm = - vib(0,:)
        acm = - aib(0,:)
        
        !! Particle
        rb = yb(NP+3 : NP+4)
        vb = yb(NP+5 : NP+6)
        call accbar(m(0:Nboul), rb, rib(0:Nboul,:2), ab) ! Este lo setea a 0 al principio
        call accsto(t, rb, vb, ab)  ! Este suma aceleración
        rib(FP,:) = rb
        vib(FP,:) = vb
        aib(FP,:) = ab

        !! Elements
        xc = (/rb(1),rb(2),cero,vb(1),vb(2),cero/)
        call elem(mcm, xc, ea, ee, ei, eM, ew, eO)
        eR = (ea/a_corot)**(1.5)

        !! Chaos
        call chaosvals(ea, ee, amax, amin, emax, emin)

        ! Output
        if (datas .and. .not. perc) write(*,*) t, ea, ee, ew, eR
        if (datao) then
            if (eleout) then
                write(2,*) t, ea, ee, eM, ew, eR
            else
                do i = 0, Nboul
                    write(2,*) i, t, rib(i,:), vib(i,:), aib(i,:), m(i), radius(i)
                end do
                write(2,*) FP, t, rb, vb, ab, cero, cero
            end if
        end if
        if (perc) call percentage(t, tf)
    end do
    if (perc) call percentage(tf+1., tf)

    if (datao) close(2)
    if (screen .and. datao) write(*,*) "Se guardó la salida en el archivo: ", trim(datafile)

    ! Output chaos
    da = amax - amin
    de = emax - emin
    if (chaos) then
        open (3, file=trim(chaosfile), status='unknown', position='append')
        ! numero simu, mala?, t total, a ini, e ini, M ini, w ini, Res ini, t integrado, da, de, a final, e final, M final, w final, R final 
        write (3,*) nsim, bad, tf, xe0(1), xe0(2), xe0(4), xe0(5), Res0, t, da, de, ea, ee, eM, ew, eR        
        close(3)
        if (screen) then
            write(*,*) "amin [km]:", amin, "amax [km]:", amax, "da [km]:", da
            write(*,*) "emin     :", emin, "emax     :", emax, "de     :", de
            write(*,*) "Se guardó el caos en el archivo: ", trim(chaosfile)
        end if
    end if

end program main

subroutine mapas_pot(N,ngx,ngy,xmin,xmax,ymin,ymax,rib,m,omega,map_file)
    use gravity, only: potast, potbar, potrot, accbar, accast, accrot
    implicit none
    integer(kind=4), intent(in) :: N,ngx, ngy
    real(kind=8), intent(in)    :: xmin, xmax, ymin, ymax, rib(0:N,2), m(0:N), omega
    real(kind=8) :: xyb(ngx,ngy), xya(ngx,ngy), xyr(ngx,ngy)
    real(kind=8) :: acb(ngx,ngy,2), aca(ngx,ngy,2), acr(ngx,ngy,2)
    real(kind=8) :: rcm(2), rb(2), ra(2), va(2), ab(2), aa(2), ar(2)
    real(kind=8) :: ria(N,2)
    integer(kind=4) :: i,j
    character(len=*) :: map_file

    !! Mapa de potencial
    rcm(1) = -rib(0,1)
    rcm(2) = -rib(0,2)
    do i = 1, N
        ria(i,1) = rib(i,1) + rcm(1)
        ria(i,2) = rib(i,2) + rcm(2)
        write(*,*) rcm, rib(i,:), ria(i,:)
    end do
    va = 0.d0
    do i = 1, ngx
        do j = 1, ngy
            rb(1) = xmin + i * (xmax-xmin)/ngx
            rb(2) = ymin + j * (ymax-ymin)/ngy
            xyb(i,j) = potbar(rb, rib)
            ra = rb + rcm
            xya(i,j) = potast(ra, ria)
            xyr(i,j) = potrot(ra, ria)
            call accbar(m, rb, rib, ab)
            acb(i,j,:) = ab
            call accast(omega, m, ra, ria, aa)
            aca(i,j,:) = aa
            call accrot(omega, m, ra, va, ria, ar)
            acr(i,j,:) = ar
        end do
    end do
    ! write(*,*) "Potencial calculado. Escribiendo..."
    open(unit=7, file=trim(map_file), status='replace', action='write')
    do i = 1, ngx
        do j = 1, ngy
            write(7,*) xmin + i * (xmax-xmin)/ngx, ymin + j * (ymax-ymin)/ngy, &
            xyb(i,j), xya(i,j), xyr(i,j), acb(i,j,:), aca(i,j,:), acr(i,j,:)
        end do
    end do
    close(7)
end subroutine mapas_pot

subroutine coord(msum, a, e, inc, capm, omega, capom, xc)
    use const, only: G, uno
    implicit none
    real(kind=8), intent(in)  :: msum, a, e, inc, capm, omega, capom
    real(kind=8), intent(out) :: xc(6)
    real(kind=8) :: sp, cp, so, co, si, ci
    real(kind=8) :: d11, d12, d13, d21, d22, d23
    real(kind=8) :: cape, dummy, scap, ccap, sqe, sqgma
    real(kind=8) :: ri, xfac1, xfac2, vfac1, vfac2
    
    ! Generate rotation matrices (on p. 42 of Fitzpatrick)
    sp = sin(omega)
    cp = cos(omega)
    so = sin(capom)
    co = cos(capom)
    si = sin(inc)
    ci = cos(inc)
    d11 = cp * co - sp * so * ci
    d12 = cp * so + sp * co * ci
    d13 = sp * si
    d21 = -sp * co - cp * so * ci
    d22 = -sp * so + cp * co * ci
    d23 = cp * si
    
    ! Get the other quantities depending on orbit type (i.e. ialpha)
    call aver(capm, e, cape, dummy)
    scap = sin(cape)
    ccap = cos(cape)
    sqe = sqrt(uno - e * e)
    sqgma = sqrt(G * msum * a)
    xfac1 = a * (ccap - e)
    xfac2 = a * sqe * scap
    ri = uno / (a * (uno - e * ccap))
    vfac1 = -ri * sqgma * scap
    vfac2 = ri * sqgma * sqe * ccap
    
    xc(1) = d11 * xfac1 + d21 * xfac2
    xc(2) = d12 * xfac1 + d22 * xfac2
    xc(3) = d13 * xfac1 + d23 * xfac2
    xc(4) = d11 * vfac1 + d21 * vfac2
    xc(5) = d12 * vfac1 + d22 * vfac2
    xc(6) = d13 * vfac1 + d23 * vfac2
end subroutine coord

subroutine aver(dm, e, u, f)
    use const, only: uno, uno2, dos, eps
    implicit none
    real(kind=8), intent(in)  :: dm, e
    real(kind=8), intent(out) :: u, f
    real(kind=8) :: u0, dif, seno, cose
    
    u0 = dm
    dif = uno
    do while (dif > eps)
        u = dm + e * sin(u0)
        dif = abs(u - u0)
        u0 = u
    end do
    seno = sqrt(uno + e) * sin(uno2 * u)
    cose = sqrt(uno - e) * cos(uno2 * u)
    f = dos * atan2(seno, cose)
end subroutine aver

subroutine elem(msum, xc, a, e, inc, capm, omega, capom)
    use const
    implicit none
    real(kind=8), intent(in)  :: msum, xc(6)
    real(kind=8), intent(out) :: a, e, inc, capm, omega, capom
    real(kind=8) :: gmsum
    real(kind=8) :: x, y, z, vx, vy, vz
    real(kind=8) :: hx, hy, hz, h2, h, fac, u
    real(kind=8) :: r, v, v2, vdotr, energy
    real(kind=8) :: cape, cw, sw, w, face, capf, tmpf
    integer(kind=4) :: ialpha = 0
        
    gmsum = G * msum
    x = xc(1)
    y = xc(2)
    z = xc(3)
    vx = xc(4)
    vy = xc(5)
    vz = xc(6)
    
    hx = y * vz - z * vy
    hy = z * vx - x * vz
    hz = x * vy - y * vx
    h2 = hx * hx + hy * hy + hz * hz
    h = sqrt(h2)
    inc = acos(hz / h)
    
    fac = sqrt(hx*hx + hy*hy) / h
    if (fac < eps) then
        capom = cero
        u = atan2(y, x)
        if (abs(inc - pi) < 10.d0 * eps) then
            u = -u
        end if
    else
        capom = atan2(hx, -hy)
        u = atan2(z / sin(inc), x * cos(capom) + y * sin(capom))
    end if
    
    if (capom < cero) then
        capom = capom + twopi
    end if
    if (u < cero) then
        u = u + twopi
    end if
    
    r = sqrt(x * x + y * y + z * z)
    v2 = vx * vx + vy * vy + vz * vz
    v = sqrt(v2)
    vdotr = x * vx + y * vy + z * vz
    energy = uno2 * v2 - gmsum / r
    
    if (abs(energy * r / gmsum) < sqrt(eps)) then
        ialpha = 0
    else
        if (energy < cero) then
            ialpha = -1
        else if (energy > cero) then
            ialpha = 1
        end if
    end if
    
    !! Ellipse
    if (ialpha == -1) then
        a = -uno2 * gmsum / energy
        fac = uno - h2 / (gmsum * a)
        if (fac > eps) then
            e = sqrt(fac)
            face = (a - r) / (a * e)
            if (face > uno) then
                cape = cero
            else
                if (face > -uno) then
                    cape = acos(face)
                else
                    cape = pi
                end if
            end if
            if (vdotr < cero) then
                cape = twopi - cape
            end if
            cw = (cos(cape) - e) / (uno - e * cos(cape))
            sw = sqrt(uno - e * e) * sin(cape) / (uno - e * cos(cape))
            w  = atan2(sw, cw)
            if (w < cero) then
                w = w + twopi
            end if
        else
            e    = cero
            w    = u
            cape = u
        end if
        capm  = cape - e * sin(cape)
        omega = u - w
        if (omega < cero) then
            omega = omega + twopi
        end if
        omega = omega - int(omega / twopi) * twopi
    end if
    
    !! Hypérbola
    if (ialpha == 1) then
        a = uno2 * gmsum / energy
        fac = h2 / (gmsum * a)
        if (fac > eps) then
            e = sqrt(uno + fac)
            tmpf = (a + r) / (a * e)
            if (tmpf < uno) then
                tmpf = uno
            end if
            capf = log(tmpf + sqrt(tmpf * tmpf - uno))
            if (vdotr < cero) then
                capf = -capf
            end if
            cw = (e - cosh(capf)) / (e * cosh(capf) - uno)
            sw = sqrt(e * e - uno) * sinh(capf) / (e * cosh(capf) - uno)
            w  = atan2(sw, cw)
            if (w < cero) then
                w = w + twopi
            end if
        else
            e = uno
            tmpf = uno2 * h2 / gmsum
            w = acos(dos * tmpf / r - uno)
            if (vdotr < cero) then
                w = twopi - w
            end if
            tmpf = (a + r) / (a * e)
            capf = log(tmpf + sqrt(tmpf * tmpf - uno))
        end if
        capm = e * sinh(capf) - capf
        omega = u - w
        if (omega < cero) then
            omega = omega + twopi
        end if
        omega = omega - int(omega / twopi) * twopi
    end if

    !! Parábola    
    if (ialpha == 0) then
        a = uno2 * h2 / gmsum
        e = uno
        w = acos(dos * a / r - uno)
        if (vdotr < cero) then
            w = twopi - w
        end if
        tmpf = tan(uno2 * w)
        capm = tmpf * (uno + tmpf * tmpf / 3.d0)
        omega = u - w
        if (omega < cero) then 
            omega = omega + twopi
        end if
        omega = omega - int(omega / twopi) * twopi
    end if
end subroutine elem

subroutine chaosvals(ea,ee,amax,amin,emax,emin)
    implicit none
    real(kind=8), intent(in)    :: ea, ee
    real(kind=8), intent(inout) :: amax,amin,emax,emin
    amax = max(amax,ea)
    amin = min(amin,ea)
    emax = max(emax,ee)
    emin = min(emin,ee)
end subroutine chaosvals

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
    if (iper.lt.10) then
        write (66,'(i2)') iper
    else
        write (66,'(i3)') iper
    end if
    rewind (66)
    read (66,'(a)') cestado
    close (66)
    
    110 format (2a)
    if (iper.lt.100) then
        write (*,110,advance='no') cret,trim(guiones)//trim(cestado)//'%'
    else
        write (*,110,advance='no') cret,trim(guiones)//'. FIN'
        write (*,*)
    end if
end subroutine percentage
