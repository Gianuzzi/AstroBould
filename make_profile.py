import os
import numpy as np

## Este programa sirve para crear el archivo
##  TOM (Tiempos, Omega(TIempos), Masa agregada(Tiempos)

## Para crearlo, necesitamos:
##  (1) Nombre de archivo TOM a crear
##  (2) El archuvo 'sump' creado con parallel.
##  (3) Parámetros del disco a crear (masa, y perfil).
##  (4) Información de las condiciones iniciales. 
##        Para esto último hay 2 opciones de uso:
##         (a) Leer parámetros iniciales del archivo config.ini [Default].
##         (b) Introducir valores en este código manualmente.

### (1) Archivo TOM
tomfile = "tomfile.dat"

### (2) Archivo sump
sump_file = "sump.dat"

### (3) Parámetros del disco a crear
mu_D = 0.01 # Cociente de masas de disco a asteroide (mDisk = mu_D * m0)
alpha = 0 # Sigma(r) = r**alpha  [alpha=0 == Homogéneo]

### (4) Condiciones iniciales
# (a)
confini = "config.ini" # "" o None, para usar lo explícito en (b)

# (b)
## Primary
m0 = 6.3e18 # [kg]
lamdba_k = 0.471e0 # ratio of spin to keplerian rotation [spin/wk] (0 if not used)
Prot = 0.2916667e0 # [day]
## Boulders
mu_B = [0.1] # Cociente de masas de boulders a asteroide (mBoul = mu_B * m0)
##### Dato: Si son varios boulders, se introducen todos como lista: [mu1, mu2, ...]

## Units [No tocar, a menos que se haya modificado en const.F90]
G = 4.9823394e-10 # [km^3 kg^(-1) day^(-2)]

####################################################################################
####################################################################################
####################################################################################

def dM(r1, r2, alpha=0, Sigma0=1):
    if Sigma0 != 0:
        return 2 * np.pi * Sigma0  / (alpha + 2) * (r2**(alpha + 2) - r1**(alpha + 2))
    return r2**(alpha + 2) - r1**(alpha + 2) # Después se normaliza por fuera
    
def Lines2015(r, alpha=0, Rgap=1, ratio=0.1): # Genera Sigma(r)
    fgap = 1 / (1 + np.exp(- (r - Rgap) / (Rgap * ratio)))
    Sigma = fgap * r**alpha
    return Sigma

def main():
    ##### Set parameters #####
    # Checks
    ## Alpha
    if alpha <= -2:
        print("Alpha debe ser mayor que -2.")
        print("Saliendo.")
        exit(1)
    ## Mu List
    if not hasattr(mu_B, "__iter__"): mu_B = [mu_B]

    # Get inferred parameters
    ## Mass
    mcm = m0 * (1 + sum(mu_B))
    mD  = mu_D * mcm
    ## Spin
    Omega_k = np.sqrt(G * mcm / R0**3)
    if lambda_k != 0:
        Omega0 = lamdba_k * Omega_k
    else:
        Omega0 = 2 * np.pi / Prot

    # Read sump sump_file
    data = np.genfromtxt("sump.out")
    bad  = data[:,1].astype(int)
    L0 = data[0,2]
    ttot = data[:,3]
    aini = data[:,4]
    n_Oini = data[:,8]
    lini = data[:,9]
    tfin = data[:,10]
    lfin = data[:,16]
    de = data[:,18]

    ###### Set each particle mass #####

    ## Check if sorted a
    srtda = all(aini[i] <= aini[i+1] for i in range(len(aini) - 1))
    if not srtda:
        asort = aini
        de_sort = np.argsort(np.argsort(aini))
    else:
        asort = aini
        de_sort = np.arange(len(aini))
        
    ## Get dr for all
    dr = np.diff(asort)
    amin = asort[0]
    amax = asort[-1]

    ## Set bin edges positions
    rmin = max(amin - 0.5 * (asort[1] - amin), 0)
    rmax = amax + 0.5 * (amax - asort[-2])
    redges = 0.5 * (asort[:-1] + asort[1:])
    redges = np.concatenate(([rmin], redges, [rmax]))

    ## Calculate area of each bin
    abins = np.pi * (redges[1:]**2 - redges[:-1]**2)
    acum = np.cumsum(abins)

    ## Calculate mass at bins
    ### Method 1 (Se integra Sigma * dA, desde r1 a r2)
    Sigma0 = mD / (2 * np.pi) * (alpha + 2) / (rmax**(alpha + 2) - rmin**(alpha + 2))
    mpart = dM(redges[:-1], redges[1:], alpha=alpha, Sigma0=Sigma0)
    ### Method 2 (better performance?) [Default]
    mpart2 = asort**alpha * abins
    mpart2 = mbins2 / np.sum(mbins2) * mD
    ### Method 3 (smoother, made for a disk with a gap)
    aux = rmin if rmin > 0 else amin * 0.5
    Sigma3 = Lines2015(asort, alpha=alpha, Rgap=aux, ratio=0.01)
    mpart3 = Sigma3 * abins
    mpart3 = mbins3 / np.sum(mbins3) * mD
    ###

    ## De-Order mass, in case not sorted a
    if not srtda:
        mpart = mpart[de_sort]
        mpart2 = mpart2[de_sort]
        mpart3 = mpart3[de_sort]
        
    ###### Set each particle event Omega, and Delta M #####
    
    ## Get time events order
    argsrt = np.argsort(tfin)

    ## Initial Lambda = L/Omega
    Lambda0 = L0/Omega0

    ## Initial asteroid mass
    Mast0 = mcm

    ## Arrays to fill
    times_tom = np.zeros(tfin+1)
    omega_tom = np.zeros(tfin+1)
    dmass_tom = np.zeros(tfin+1)
    Lambda = np.zeros(tfin+1)
    Mast = np.zeros(tfin+1)

    ## Initial values
    omega_tom[0] = omega0
    Lambda[0] = Lambda0
    Mast[0] = Mast0


    ## Get beauge's file name
    beauge = "beauge" 
    aux = len(beauge)
    i = 1
    while os.path.isfile(beauge+".dat"):
        beauge = beauge[:aux] + str(i)
        i += 1

    ## Loop and write beauges file
    with open(beauge, "w") as beaugef:
        for i,j in enumerate(argsrt, 1):
            times_tom[i] = tfin[j]
            deltal = lfin[j] - lini[j]
            ratio_l = lfin[j] / lini(j)
            if bad[j] == 0: # Survived
                Mast[i] = Mast[i-1]
                Lambda[i] = Lambda[i-1]
                omega_tom[i] = omega_tom[i-1]
            elif bad[j] == 1: # Collision
                Mast[i] = Mast[i-1] + mpart2[j]
                Lambda[i] = Lambda[i-1] * (1 + mpart2[j]/Mast[i-1])
                omega_tom[i] = omega_tom[i-1] * (Lambda[i-1]/Lambda[i]) + (mpart2[j] * lini[j]) / Lambda[i]
                dmass_tom[i] = mpart2[j]
            elif bad[j] == 2: # Escape
                Mast[i] = Mast[i-1]
                Lambda[i] = Lambda[i-1]
                omega_tom[i] = omega_tom[i-1] - mpart2[j] * deltal / Lambda[i]
            else: # ERROR
                raise ValueError("Bad '%d' no reconocido en partícula '%d'."%(bad[j],i-1))
            beaugef.write(i,j,bad[j],n_Oini[j],tfin[j],de[j],ratio_l,omega_tom[i],Mast[i])
    print("Se ha escrito el archivo %s"%beauge)

    ## Create TOM data, and write to file
    tom = np.concatenate((times_tom[1:], omega_tom[1:], mass_tom[1:]), axis=1)
    np.savetxt(tomfile, tom, delimiter=' ')
    print("Se ha creado el archivo TOM: %s"%tomfile)


if name == "__main__":
    main()
