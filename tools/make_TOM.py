import os
import warnings

import numpy as np

# Este programa sirve para crear el archivo
#  TOM: Tiempos, Omega(Tiempos), Masa agregada(Tiempos)

# Para crearlo, necesitamos:
#  (1) Nombre de archivo TOM a crear
#  (2) El archivo 'chaos' creado con parallel.
#  (3) Parámetros del disco a crear (masa, y perfil).
#  (4) Información de las condiciones iniciales.
#        Para esto último hay 2 opciones de uso:
#         (a) Leer parámetros iniciales del archivo config.ini [Default].
#         (b) Introducir valores en este código manualmente.

# (1) Archivo TOM a crear
tomfile = "tomfile.dat"

# (2) Archivo chaos a utilizar para crear TOM
chaos_file = "chaos.out"

# (3) Parámetros del disco a crear
mu_d = 0.01  # Cociente de masas de disco a asteroide (mDisk = mu_D * m0)
alpha = 0  # Sigma(r) = r**alpha  [alpha=0 == Homogéneo]

# (4) Condiciones iniciales
# (a)
confini = "config.ini"  # "" o None, para usar lo explícito en (b)

# (b)
# Primary
m0 = 6.3e18  # [kg]
R0 = 129  # [km]
lambda_k = (
    0.471e0  # ratio of spin to keplerian rotation [spin/wk] (0 if not used)
)
Prot = 0.2916667e0  # [day]
# Boulders
mu_b = [0.1]  # Cociente de masas de boulders a asteroide (mBoul = mu_B * m0)
# Dato: Si son varios boulders, se introducen todos como lista: [mu1, mu2, ...]


# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# ----------------------- No tocar a partir de aquí -----------------------
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------


def read_config(config_file="config.ini"):
    global m0, R0, lambda_k, Prot, mu_b
    m0 = 0.0
    R0 = 0.0
    lambda_k = 0.0
    Prot = 0.0
    mu_b = []
    nboul = 0
    # Read config file
    with open(config_file, "r") as f:
        lines = f.readlines()
    get_val = (
        lambda x: (x.replace("d", "e").split(":")[1]).split("!")[0].strip()
    )
    j = 0
    for i, line in enumerate(lines):
        if j > 0:
            j -= 1
            continue
        this = line.strip()[:15]
        if this == "mass of primary":
            m0 = float(get_val(line))
        elif this == "radius of prima":
            R0 = float(get_val(line))
        elif this == "ratio of spin t":
            lambda_k = float(get_val(line))
        elif this == "rotational peri":
            Prot = float(get_val(line))
        elif this == "number of bould":
            nboul = int(get_val(line))
        elif this[:7] == "mass_m0":
            if nboul == 0:
                msg = "Error en número de boulders en config.ini."
                raise ValueError(msg)
            io = i
            for j in range(nboul):
                while True:
                    nextl = lines[io + 1].strip()
                    if (
                        (len(nextl) == 0)
                        or (nextl[:2] == "c ")
                        or (nextl[:2] == "! ")
                    ):
                        io += 1
                        continue
                    val = (nextl.split()[0]).replace("d", "e")
                    mu_b.append(float(val))
                    break
                io = io + 1
        else:
            continue


def get_dm(r1, r2, alpha=0, sigma0=1):
    if sigma0 != 0:
        return (
            2
            * np.pi
            * sigma0
            / (alpha + 2)
            * (r2 ** (alpha + 2) - r1 ** (alpha + 2))
        )
    return r2 ** (alpha + 2) - r1 ** (
        alpha + 2
    )  # Después se normaliza por fuera


def lines2015(r, alpha=0, rgap=1, ratio=0.1, sigma0=1):  # Genera Sigma(r)
    fgap = 1.0 / (1.0 + np.exp(-(r - rgap) / (rgap * ratio)))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sigma = sigma0 * fgap * np.where(r > 0, np.power(r, alpha), 0.0)
    return sigma


def lynden_bell1974(
    r, dzeta=0.75, r0=1, sigma0=1, t=0, tdiff=1
):  # Genera Sigma(r)
    t_s = 1 + t / tdiff
    ft_s = np.power(t_s, -(2.5 - dzeta) / (2 - dzeta))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        fexp = np.exp(-np.where(r > 0, np.power(r / r0, 2 - dzeta), 0.0) / t_s)
        sigma = (
            sigma0 * ft_s * np.where(r > 0, np.power(r, -dzeta), 0.0) * fexp
        )
    return sigma


if __name__ == "__main__":
    # Set parameters #
    # Checks
    # Alpha
    if alpha <= -2:
        print("Alpha debe ser mayor que -2.")
        print("Saliendo.")
        exit(1)
    # Read config file
    if confini:
        read_config(confini)
    # Mu List
    if not hasattr(mu_b, "__iter__"):
        mu_b = [mu_b]
    # Units [No tocar, a menos que se haya modificado en const.F90]
    unit_m = 1.0e-15  # [kg]
    unit_r = 1.0e0  # [km]
    unit_t = 1.0e0  # [day]
    unit_v = unit_r / unit_t  # [km/day]
    unit_a = unit_v / unit_t  # [km/day^2]
    G_aux = 4.9823394e-10  # [km^3 kg^(-1) day^(-2)]
    G = (
        G_aux * (unit_r**3) / unit_m / unit_t
    )  # [unit_d^3 unit_m^(-1) unit_t^(-2)]
    m0 *= unit_m
    R0 *= unit_r
    Prot *= unit_t

    # Get inferred parameters
    # Mass
    mcm = m0 * (1 + sum(mu_b))
    md = mu_d * mcm
    # Spin
    Omega_k = np.sqrt(G * mcm / R0**3)
    if lambda_k != 0:
        Omega0 = lambda_k * Omega_k
    else:
        Omega0 = 2 * np.pi / Prot

    # Read chaos chaos_file
    # 0    ! i
    # 1    ! bad
    # 2    ! total time to integrate
    # 3    ! initial (Asteroid): angular momentum
    # 4-9  ! initial: mass, a, e, M, omega, MMR
    # 10   ! initial: angular momentum per unit mass
    # 11   ! surviving time
    # 12-16! final: a, e, M, omega, MMR
    # 17   ! final: angular momentum per unit mass
    # 18-19! a_min, a_max
    # 20-21! e_min, e_max
    # 22   ! Delta a
    # 23   ! Delta e
    data = np.genfromtxt("chaos.out")

    idx = data[:, 0].astype(int)
    bad = data[:, 1].astype(int)
    L0 = data[0, 3] * unit_m * unit_r**2 / unit_t
    aini = data[:, 5] * unit_r
    MMR = data[:, 9]
    lini = data[:, 10] * unit_r**2 / unit_t
    tfin = data[:, 11] * unit_t
    lfin = data[:, 17] * unit_r**2 / unit_t
    de = data[:, 23]

    # Set each particle mass #

    # Check if sorted a
    srtda = all(aini[i] <= aini[i + 1] for i in range(len(aini) - 1))
    if not srtda:
        asort = aini
        de_sort = np.argsort(np.argsort(aini))
    else:
        asort = aini
        de_sort = np.arange(len(aini))

    # Get dr for all
    dr = np.diff(asort)
    amin = asort[0]
    amax = asort[-1]

    # Set bin edges positions
    rmin = max(amin - 0.5 * (asort[1] - amin), 1e-10)
    rmax = amax + 0.5 * (amax - asort[-2])
    redges = 0.5 * (asort[:-1] + asort[1:])
    redges = np.concatenate(([rmin], redges, [rmax]))

    # Calculate area of each bin
    abins = np.pi * (redges[1:] ** 2 - redges[:-1] ** 2)

    # Calculate mass at bins
    # Method 1 (Se integra Sigma * dA, desde r1 a r2)
    Sigma0 = (
        md
        / (2 * np.pi)
        * (alpha + 2)
        / (rmax ** (alpha + 2) - rmin ** (alpha + 2))
    )
    mbin1 = get_dm(redges[:-1], redges[1:], alpha=alpha, sigma0=Sigma0)
    # Method 2 (better performance?) [Default]
    mbin2 = asort**alpha * abins
    mbin2 = mbin2 / np.sum(mbin2) * md
    # Method 3 (smoother, made for a disk with a gap)
    aux = rmin if rmin > 0 else amin * 0.5
    Sigma3 = lines2015(asort, alpha=alpha, rgap=aux, ratio=0.01)
    mbin3 = Sigma3 * abins
    mbin3 = mbin3 / np.sum(mbin3) * md
    #

    # De-Order mass, in case not sorted a
    if not srtda:
        mbin1 = mbin1[de_sort]
        mbin2 = mbin2[de_sort]
        mbin3 = mbin3[de_sort]
    mpart = mbin3  # Asignamos el método 3

    # Create mass profile file
    massdata = np.vstack((aini / unit_r, mpart / unit_m, abins / unit_r**2)).T
    np.savetxt("massfile.dat", massdata, delimiter=" ")
    print(
        "Se ha creado el archivo de perfil de masa: "
        + "massfile.dat [a, mpart, dA]"
    )

    # Set each particle event Omega, and Delta M #

    # Get time events order
    argsrt = np.argsort(tfin)

    # Initial Inertia = L/Omega
    Inertia0 = L0 / Omega0

    # Initial asteroid mass
    Mast0 = mcm

    # Arrays to fill
    times_tom = np.zeros(len(tfin) + 1)
    omega_tom = np.zeros(len(tfin) + 1)
    dmass_tom = np.zeros(len(tfin) + 1)
    Inertia = np.zeros(len(tfin) + 1)
    Mast = np.zeros(len(tfin) + 1)

    # Initial values
    omega_tom[0] = Omega0
    Inertia[0] = Inertia0
    Mast[0] = Mast0

    # Get beauge's file name
    beauge = "beauge.dat"

    # Loop and write beauges file
    formato = "%5d %5d %d %10.5e %10.5e %10.5e %10.5e %10.5e %10.5e\n"
    with open(beauge, "w") as beaugef:
        for i, j in enumerate(argsrt, 1):
            times_tom[i] = tfin[j]
            deltal = lfin[j] - lini[j]
            ratio_l = lfin[j] / lini[j]
            if bad[j] == 0:  # Survived
                Mast[i] = Mast[i - 1]
                Inertia[i] = Inertia[i - 1]
                omega_tom[i] = omega_tom[i - 1]
            elif bad[j] == 1:  # Collision
                Mast[i] = Mast[i - 1] + mpart[j]
                Inertia[i] = Inertia[i - 1] * (1 + mpart[j] / Mast[i - 1])
                omega_tom[i] = (
                    omega_tom[i - 1] * (Inertia[i - 1] / Inertia[i])
                    + (mpart[j] * lini[j]) / Inertia[i]
                )
                dmass_tom[i] = mpart[j]
            elif bad[j] == 2:  # Escape
                Mast[i] = Mast[i - 1]
                Inertia[i] = Inertia[i - 1]
                omega_tom[i] = (
                    omega_tom[i - 1] - mpart[j] * deltal / Inertia[i]
                )
            else:  # ERROR
                raise ValueError(
                    "Bad '%d' no reconocido en partícula '%d'."
                    % (bad[j], i - 1)
                )
            beaugef.write(
                formato
                % (
                    i,
                    j,
                    bad[j],
                    MMR[j],
                    tfin[j] / unit_t,
                    de[j],
                    ratio_l,
                    omega_tom[i] * unit_t,
                    Mast[i] / unit_m,
                )
            )
    print("Se ha escrito el archivo %s" % beauge)

    # Create TOM data, and write to file
    # Eliminamos la condición inicial,
    # y los eventos que no cambian nada (survived)
    event = bad[argsrt] > 0
    times_tom = times_tom[1:][event]
    omega_tom = omega_tom[1:][event]
    dmass_tom = dmass_tom[1:][event]

    # Acoplamos los datos que tengan mismo tiempo de evento
    aux = np.unique(times_tom, return_counts=True)
    if np.any(aux[1] - 1):
        times_tom = aux[0]
        j = 0
        for i, n in enumerate(aux[1]):
            omega_tom[i] = omega_tom[j + n - 1]
            dmass_tom[i] = np.sum(dmass_tom[j : j + n])
            j += n
        omega_tom = omega_tom[: len(times_tom)]
        dmass_tom = dmass_tom[: len(times_tom)]
    # Check if anythig to write
    if len(times_tom) == 0:
        print("No hay eventos para escribir en el archivo TOM.")
        print("Saliendo.")
        exit(1)
    # CAMBIAMOS OMEGA POR DELTA OMEGA
    domega_tom = np.diff(np.insert(omega_tom, 0, Omega0))
    tom = np.vstack(
        (times_tom / unit_t, domega_tom * unit_t, dmass_tom / unit_m)
    ).T
    if os.path.isfile(tomfile):
        print("WARNING: Output file {} already exist.".format(tomfile))
        yes_no = input("Do you want to overwrite it? y/[n]\n")
        if yes_no.lower() not in ["y", "yes", "s", "si"]:
            i = 1
            aux = tomfile.split(".")
            suf = aux[-1] if len(aux) > 1 else ""
            while os.path.isfile(tomfile):
                tomfile = (
                    ".".join(aux[:-1]) + str(i) + ("." + suf if suf else "")
                )
                i += 1
    np.savetxt(tomfile, tom, delimiter=" ")
    print("Se ha creado el archivo TOM: %s" % tomfile, " [t, dOmega, dM]")
