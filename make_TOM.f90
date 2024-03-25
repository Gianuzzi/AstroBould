!!!  analisis.f90
!!!
!!!  Analiza resultados de una simulacion con AstroBould, y construye el
!!!  archivo de tiempos de encuentro, Omega0 y aumento de masa del asteroide.
!!!
module common
  implicit real*8 (a-h,k-z)
  parameter (imax=100000)
  real*8 twopi,cero,uno,dos,uno2,pi,G,rad,uno3,error
  save
end module common


program main ; use common
  implicit real*8 (a-h,k-z)
  integer iorden(imax),ip(imax),ibad(imax)
  real*8 m(0:10),mp(imax)
  real*8 a0(imax),e0(imax),anom0(imax),w0(imax),nrat0(imax),Lp0(imax)
  real*8 tfin(imax),a(imax),e(imax),anom(imax),w(imax),nrat(imax),Lp(imax)
  real*8 Da(imax),De(imax),temp(imax),omega(0:imax),Lambda(0:imax),Mast(0:imax)
  character(len=200) arch_in,arch_3,arch_4,add
  character(len=10) arq_isim,arq_e,arq_g
  !
  twopi = 8.0d0*datan(1.0d0)  ;  cero = 0.0d0      ;  uno   = 1.0d0
  dos   = 2.0d0               ;  uno2 = 0.5d0      ;  pi    = uno2*twopi
  rad   = pi/180.0d0          ;  uno3 = uno/3.0d0  ;  error = 1.0d-13

!!! masa total del disco (en unidades de m0).
  mdisc = 0.1

!!! lazo sobre diferentes simulaciones.  
  do igam = 1,2  ! potencia del perfil de denisdad superficial de masa
    !
    if (igam == 1) then
      gamma = -1.0  ;  arq_g = '1'
    end if
    if (igam == 2) then
      gamma = -0.5  ;  arq_g = '05'
    end if
    !
    do iexc0 = 0,2 ! excentricidad inicial de las particulas
      if (iexc0 == 0) arq_e = '0'
      if (iexc0 == 1) arq_e = '01'
      if (iexc0 == 2) arq_e = '02'
      !
      do isim = 1,4 ! elige simulacion de la fase 1
        if (isim == 1)  arq_isim = 'm1_1e-2'
        if (isim == 2)  arq_isim = 'm1_3e-2'
        if (isim == 3)  arq_isim = 'm1_1e-1'
        if (isim == 4)  arq_isim = 'm1_m2_5e-2'
  
!!! archivo con resultados de la simulación inicial con AstroBould.
        add = 'e_'//trim(arq_e)//'_'//trim(arq_isim)
        arch_in = 'Astrobould/sump_'//trim(add)//'.out'
        open(1,file=arch_in,status='old')

!!! nombre de los archivos de salida.
        add = trim(arq_g)//'_e_'//trim(arq_e)//'_'//trim(arq_isim)
        arch_3 = 'analisis_gamma_'//trim(add)//'.dat'
        arch_4 = 'tout_omega0_deltm_gamma_'//trim(add)//'.dat'
        
!!! archivo de salida del analisis de la corrida inicial con AstroBould.
        open (3,file=arch_3,status='replace')
        
!!! archivo de (t_out,Omega,delta_m) para la simulación de la Fase 2.
        open (4,file=arch_4,status='replace')
        
!!! parametros segun Emma.
        unit_m = 1.d-15         ! [Kg]
        unit_r = 1.d0           ! [Km]
        unit_t = 1.d0           ! [dia]
        unit_v = unit_r/unit_t  ! [Km/dia]
        unit_a = unit_v/unit_t  ! [Km/dia^2]
        G_aux  = 4.9823394d-10  ! [km^3 kg^(-1) dia^(-2)]
        G = G_aux*(unit_r**3)/unit_m/unit_t ! [unit_d^3 unit_m^(-1) unit_t^(-2)]
  
!!! datos del cuerpo central.
        m(0)   = 6.3d18         ! masa  de m0 [kg]
        R0     = 129.d0         ! radio de m0 [km]
        O0_kep = 0.471d0        ! espin de m0 en unidades de omega_kep (lambda)
        
!!! datos de la(s) anomalía(s) de masa.
        if (isim == 1) m(1) = 0.01*m(0)
        if (isim == 2) m(1) = 0.03*m(0)
        if (isim == 3) m(1) = 0.1*m(0)
        if (isim == 4) m(1) = 0.1*m(0)
  
!!! cambio de unidades y construcción del factor GM.
        m    = m*unit_m         ! reescaleo de las masas segun unit_m
        mcm  = sum(m)           ! masa total del sistema
        mbol = mcm - m(0)       ! masa de todos los boulders
        GM   = G*mcm
        
!!! momento de inercia energia del cuerpo central.
        C0 = (2.0/5.0)*R0*m(0)
        
!!! frecuencia de rotación inicial del cuerpo central.
        omega_kep = sqrt(GM/R0**3)/unit_t   ! frec kepleriana sobre m0
        omega0    = O0_kep*omega_kep
  
!!! lectura de los datos de la simulación.
        i = 1
        rmin = 1.0e33
        rmax = cero
1       read (1,*,end=2) ip(i),ibad(i),L0,tstop,a0(i),e0(i),anom0(i),w0(i), &
             nrat0(i),Lp0(i),tfin(i),a(i),e(i),anom(i),w(i),nrat(i),Lp(i), &
             Da(i),De(i)
        rmin = min(rmin,a0(i)) ; rmax = max(rmax,a0(i))
        i = i + 1
        goto 1
2       itot = i - 1
        close (1)
        
!!! asigna masa a cada particula.
        deno = (rmax/rmin)**(gamma+dos) - uno
        Sigma0 = mdisc*(m(0)/unit_m)*(gamma+dos)/twopi/rmin/rmin/deno
        do i = 1,itot
          delr = a0(2) - a0(1) ; if (i > 1) delr = a0(i) - a0(i-1)
          mp(i) = twopi*a0(i)*Sigma0*((a0(i)/rmin)**gamma)*delr
        end do
        
!!! ordena las partículas segun tiempo de escape/colisión.
        temp(1:itot) = tfin(1:itot)
        do i = 1,itot
          tmin = 1e30
          do j = 1,itot
            if (temp(j) < tmin) then
              jmin = j  ;  tmin = temp(j)
            end if
          end do
          iorden(i) = jmin
          temp(jmin) = 1e30
        end do
        
!!! valor inicial del momento angular del asteroide por unidad de espín.
        Lambda(0) = L0/omega0
        
!!! valor inicial de la masa total del asteroide.
        Mast(0) = mcm/unit_m
        
!!! cambio en el espín del asteroide despues de cada partícula eliminada.
        omega(0) = omega0
        do i = 1,itot
          j = iorden(i)
          deltaL = Lp(j) - Lp0(j)    ! cambio en momento angular de la partícula
          if (ibad(j) == 0) then     ! la partícula sobrevivió
            Mast(i)   = Mast(i-1)
            Lambda(i) = Lambda(i-1)
            omega(i)  = omega(i-1)
            delta_m   = cero
          else if (ibad(j) == 1) then  ! la partícula chocó con el asteroide
            Mast(i)   = Mast(i-1) + mp(j)
            Lambda(i) = Lambda(i-1)*(uno + mp(j)/Mast(i-1))
            omega(i) = omega(i-1)*(Lambda(i-1)/Lambda(i))+mp(j)*Lp0(j)/Lambda(i)
            delta_m   = mp(j)
          else if (ibad(j) == 2) then  ! la partícula escapó
            Mast(i)   = Mast(i-1)
            Lambda(i) = Lambda(i-1)
            omega(i)  = omega(i-1) - mp(j)*deltaL/Lambda(i)
            delta_m   = cero
          end if
          ratL = Lp(j)/Lp0(j)
          write (3,100) i,j,ibad(j),nrat0(j),tfin(j),De(j),ratL,omega(i),Mast(i)
          write (4,*) tfin(j),omega(i),delta_m
        end do
        close (3)
        close (4)
        !
      end do
    end do
  end do
  !
100 format (3i6,1p10e15.5)
  !
end program main
