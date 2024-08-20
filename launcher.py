# Version: 8.0

# DENTRO DE UN ENTORNO PYTHON
# Ejecución: $ python launcher.py

# Este programa integra cada sistema (condición inicial) en
# un archivo independiente,
# y luego concatena los archivos de caos en un solo archivo.

# El archivo de partículas a introducir debe tener formato:
# mass a e M w
# Pero también se puede introducir una última columna con el valor de R,
# el cual reemplaza a a: [a = R**(2/3.) * a_corot]. Entonces quedaría:
# mass a e M w R

# IMPORTANTE : Todos los archivos deben estar en la misma carpeta

# En caso de dejar puesta la salida en archivo (<datafile>),
# se creará un archivo llamado salida[id].out por cada partícula.

# Este código crea un archivo de caos (chaos[id].out) por cada partícula.
# Luego los concatena en un solo archivo <final_chaos>,con el siguiente formato:
## 0    ! Numero de partícula / simulación
## 1    ! bad? (0: no, 1: collision, 2: ejection)
## 2    ! total time to integrate
## 3    ! initial (Asteroid): angular momentum
## 4-9  ! initial: mass, a, e, M, omega, MMR
## 10   ! initial: angular momentum per unit mass
## 11   ! surviving time
## 12-16! final: a, e, M, omega, MMR
## 17   ! final: angular momentum per unit mass
## 18-19! a_min, a_max
## 20-21! e_min, e_max
## 22   ! Delta a
## 23   ! Delta e

# Excepto <final_chaos>, el resto de los archivos estarán en carpetas creadas
# con nombre 'dpi[pid]' asociado al ID del procesador que ejecutó el sistema.
# El máximo de carpetas creadas será igual a
# MAX (número de procesadores disponibles en el sistema, workers).

# Si no son necesarias, se recomienda BORRAR las carpetas creadas luego de
# terminar la ejecución.
# Esto puede hacerse con: $ rm -rf dpy*

# Observación: Si es integración con tomfile, los directorios creados
# serán tomd[pid] en vez de dpy[pid].

##############################################################################
##############################################################################

# Importamos #

import os
import sys
import subprocess
from concurrent.futures import ProcessPoolExecutor

# Editar estas líneas según corresponda.
# Estos valores tienen privilegio ante los de config.ini.

# Parámetros #
program = "ASTROBOULD"  # Nombre del ejecutable

## Configuración de integración ##
workers = 2  # Número de procesadores a usar (workers)
explicit = False  # Método: True (cos, sin), False (integra boulders y m0)
version1 = False  # Versión del código (1 o 2)

## Torque and merge ##
torque = False  # Si se quiere usar torque
merge = False  # Si se quiere usar merge

## Input ## ("" o False si no se usa)
config = "config.ini"  # Nombre del archivo de configuración
partfile = "particles.in"  # Nombre del archivo de partículas
tomfile = ""  # Nombre del archivo de valores de t_i, delta_omega(t_i), y delta_masa(t_i)

## Output ## ("" o False si no se usa)
new_dir = True  # Puede ser un string con el nombre de la carpeta a crear. Default: "simulation"
datafile = "salida"  # Nombre del archivo de salida de datos (sin extensión)
final_chaos = "sump"  # Final Chaos Summary Output file name (sin extensión)
suffix = ""  # Suffix for the output files
## Screen ##
screen_info = True  # Si se quiere ver la información en pantalla (requires all_in_one=True)
screen_data = False  # Si se quiere ver los datos en pantalla # "%" para porcentaje (requires all_in_one=True)
## Elements ##
elements = True  # Si se quiere devolver elementos orbitales (en datafile)

### If all particles must be run in same integration (requires parallel)
all_in_one = False

# ----------------------------------------------------------------------
# -------------------- No tocar de aquí en adelante --------------------
# ----------------------------------------------------------------------

# Iniciamos #

# Obtener el path actual
cwd = os.getcwd()

# Define work_dir
if isinstance(new_dir, str):
    wrk_dir = os.path.join(cwd, new_dir)
elif isinstance(new_dir, bool):
    wrk_dir = os.path.join(cwd, "simulation") if new_dir else cwd
else:
    raise ValueError("new_dir must be either a string or a boolean")

# Redefine partfile if bool
if isinstance(datafile, bool) and datafile:
    datafile = "salida"
# Redefine partffinal_chaosile if bool
if isinstance(final_chaos, bool) and final_chaos:
    final_chaos = "sump"


# Archivos
ocini = os.path.join(cwd, config)  # Archivo de configuración
oprogr = os.path.join(cwd, program)  # Ejecutable
oparticles = os.path.join(cwd, partfile)  # Archivo de partículas
otom = os.path.join(
    cwd, tomfile
)  # Archivo de valores de t_i, delta_omega(t_i), y delta_masa(t_i)

# Checkeamos los archivos
## Configuración
existe_ocini = os.path.isfile(ocini)
if not existe_ocini:
    print("WARNING: Configuration file {} does not exist.".format(ocini))
    print("         Se utilizarán los parámetros explicitados en el código, ")
    print("          en vez de los de algún archivo de parámetros. ")
    yes_no = input("¿Desea continuar? [y/[n]]\n")
    if yes_no.lower() not in ["y", "yes", "s", "si"]:
        print("Saliendo.")
        sys.exit()
## Ejecutable
if not os.path.isfile(oprogr):
    msg = "Executable file {} does not exist.".format(oprogr)
    raise FileNotFoundError(msg)
## Partículas
if not os.path.isfile(oparticles):
    msg = "Particles file {} does not exist.".format(oparticles)
    raise FileNotFoundError(msg)
## Tomfile
existe_otom = os.path.isfile(otom)
if tomfile and (not existe_otom):
    print("ERROR: Tau-Omega-Mass file {} does not exist.".format(otom))
    print("Saliendo.")
    sys.exit()
## Chaosfile
if final_chaos is None:
    final_chaos = "sump"
if os.path.isfile(os.path.join(wrk_dir, "%s.out" % final_chaos)):
    print("WARNING: Chaos Output file {} already exist.".format(final_chaos))
    if datafile:
        print(
            "         Independently, <%s> will be appended (if existing)"
            % datafile
        )
    yes_no = input("Do you want to overwrite it? y/[n]\n")
    if yes_no.lower() not in ["y", "yes", "s", "si"]:
        i = 1
        unique_file = final_chaos
        while os.path.isfile(os.path.join(wrk_dir, "%s.out" % unique_file)):
            unique_file = f"{final_chaos}_{i}"
            i += 1
        final_chaos = unique_file

# Checks #
## Si hay torque, entonces explicit debe ser false
if torque and explicit:
    print("WARNING: Torque is active. Explicit mode will be deactivated.")
    explicit = False
    yes_no = input("Do you want to continue? y/[n]\n")
    if yes_no.lower() not in ["y", "yes", "s", "si"]:
        print("Saliendo.")
        sys.exit()

if all_in_one and (not torque):
    print(
        "WARNING: All particles will be integrated together, but without torque."
    )
    yes_no = input("Do you want to continue? y/[n]\n")
    if yes_no.lower() not in ["y", "yes", "s", "si"]:
        print("Saliendo.")
        sys.exit()

if torque and (not all_in_one):
    print("WARNING: All particles will be integrated independently,")
    print(" but each of them will induce torque to its asteroid.")
    yes_no = input("Do you want to continue? y/[n]\n")
    if yes_no.lower() not in ["y", "yes", "s", "si"]:
        print("Saliendo.")
        sys.exit()


# Leemos input
## Partículas
with open(oparticles, "r") as f:
    lines = f.readlines()
# Arreglamos por si hay "e" en vez de "d"
for i in range(len(lines)):
    lines[i] = lines[i].replace("e", "d")

# Obtener el número de líneas del archivo de partículas
nsys = len(lines)
if nsys == 0:
    print("No hay partículas para integrar en el archivo %s." % partfile)
    print("Saliendo.")
    sys.exit()
else:
    print("Cantidad total de partículas: {}".format(nsys))

# Ver si hay que hacer todo, o ya hay alguna realizadas
new_simulation = True
missing_lines = range(1, nsys + 1)
if (
    not all_in_one
):  # Si hay all_in_one, cagaste. La única prueba es si ya hay un sump.
    # Definimos prefijo de directorio, de acuerdo a TOMfile o no
    pref = "tomd" if tomfile else "dpy"

    # Obtener los sistemas realizados
    if os.path.isdir(wrk_dir):
        print(
            "Checkeando integraciones ya completadas"
            + " dentro de %s..."
            % (
                os.path.basename(wrk_dir)
                if cwd != wrk_dir
                else "este directorio"
            )
        )

        # Lista de archivos done.txt
        file_list = []

        # Recorre todas las subcarpetas en la carpeta raíz
        for subdir in os.listdir(wrk_dir):
            # Verifica si el nombre de la subcarpeta comienza con 'pref'
            this_list = []
            if subdir.startswith(pref):
                subdir_path = os.path.join(wrk_dir, subdir)
                # Recorre todos los archivos en la subcarpeta
                this_list = [os.path.join(subdir_path, "done.txt")
                             for f in os.listdir(subdir_path)
                             if f == "done.txt"]
            file_list.extend(this_list)
            
        done = []
        for file in file_list:
            with open(file, "r") as f:
                done.extend([int(cint) for cint in f.readlines() if cint != ""])
        missing_lines = [x for x in range(1, nsys + 1) if x not in done]
        print("   Cantidad de sistemas ya integrados: {}".format(len(done)))
        nsys = len(missing_lines)
        print("   Cantidad de sistemas a integrar: {}".format(nsys))
        new_simulation = False


# Hay que hacer?
if len(missing_lines) == 0:
    print("Ya se han integrado todas las partículas.")
else:
    # Obtener el número de workers
    workers = min(max(1, min(int(workers), len(os.sched_getaffinity(0)))), nsys)
    print("Workers: {}\n".format(workers))

# Argumentos. Estos son:
args = " --nomapf"
if all_in_one:
    args += " --screen" if screen_info else " --noscreen"
    if screen_data == "%":
        args += " --perc --nodatascr"
    else:
        args += " --datascr" if screen_data else " --nodatascr"
        args += " --noperc"
    args += " -partfile %s" % partfile
    args += " -datafile %s.out" % datafile if datafile else " --nodataf"
    args += " -parallel %d" % workers
    args += " -chaosfile %s.out" % final_chaos if final_chaos else " --nochaosf"
else:
    args += " --noscreen --nodatascr --noperc --noparallel --nopartfile"
args += " --explicit" if explicit else " --implicit"
args += " -tomfile %s" % tomfile if tomfile else " --notomfile"
args += " --elem" if elements else " --noelem"
args += " --version1" if version1 else " --version2"
args += " --torque" if torque else " --notorque"
args += " --merge" if merge else " --nomerge"


# Función general
def integrate_n(i):
    # Get processor ID
    pid = os.getpid()
    dirp = os.path.join(wrk_dir, "%s%d" % (pref, pid))
    nprogr = os.path.join(dirp, program)
    ncini = os.path.join(dirp, "config.ini")
    ntom = os.path.join(dirp, tomfile)
    if not os.path.exists(dirp):  # Si no existe el directorio
        subprocess.run(["mkdir", dirp], check=True)
        subprocess.run(["cp", oprogr, nprogr], check=True)
        if existe_ocini:
            subprocess.run(["cp", ocini, ncini], check=True)
        if existe_otom:
            subprocess.run(["cp", otom, ntom], check=True)
    this_args = " -nsim %d" % i
    this_args += " -chaosfile chaos%d%s.out" % (i, suffix)
    this_args += "%s" % (
        " --nodataf"
        if not datafile
        else " -datafile %s"
        % (
            "%s%d%s.out"
            % (datafile if isinstance(datafile, str) else "salida", i, suffix)
        )
    )
    print("Running system %d\n" % (i))
    my_line = " ".join(
        lines[i - 1].split()[1:]
    )  # -1 porque la lsita arranca de 0 #Elimina el primer valor (en caso que fuese la masa...)
    # ESTO SE ESTÁ EJECUTANDO EN LA SHELL #
    # print("Running: ./%s %s %s %s"%(program, args, this_args, my_line))
    # (Lines debe ser último porque termina en "\n") #
    try:
        p = subprocess.run(
            ["./%s %s %s %s" % (program, args, this_args, my_line)],
            cwd=dirp,
            check=True,
            shell=True,
        )
    except:
        print("The system %d has failed." % i)
        return
    print("System %d has been integrated." % i)
    with open(os.path.join(dirp, "done.txt"), "a") as f:
        f.write("%d\n" % i)
    return


# Crear archivo sump
def make_sum(final_chaos, suffix=""):
    # Ruta a la carpeta raíz que contiene las subcarpetas con los archivos
    root_dir = wrk_dir

    # Lista para almacenar los nombres de los archivos
    file_list = []

    # Recorre todas las subcarpetas en la carpeta raíz
    for subdir in os.listdir(root_dir):
        # Verifica si el nombre de la subcarpeta comienza con 'pref'
        if subdir.startswith(pref):
            subdir_path = os.path.join(root_dir, subdir)
            # Recorre todos los archivos en la subcarpeta
            for filename in os.listdir(subdir_path):
                # Verifica si el nombre del archivo comienza con
                # "chaos" y termina con ".out"
                if (
                    filename.startswith("chaos")
                    and filename.endswith(".out")
                    and (suffix in filename)
                ):
                    filepath = os.path.join(subdir_path, filename)
                    file_list.append(filepath)

    # Ordena los nombres de los archivos por el valor de i en "chaos%d%s.out"
    if suffix == "":
        file_list = sorted(
            file_list, key=lambda x: int(x.split("chaos")[1].split(".out")[0])
        )
    else:
        file_list = sorted(
            file_list,
            key=lambda x: int(
                x.split("chaos")[1].split(".out")[0].split(suffix)[0]
            ),
        )

    # Concatena los archivos
    outs = final_chaos.split(".")
    outs.insert(-1, suffix) if len(outs) > 1 else outs.insert(1, suffix)
    outs.insert(-1, ".out")
    final_chaos = "".join(outs)
    command = ' '.join([file for file in file_list])
    try:
        p = subprocess.run(
            ["cat %s > %s" % (command, final_chaos)],
            cwd=wrk_dir, check=True, shell=True
        )
    except:
        print("Could not create %s with cat."%final_chaos)
        print(" Trying with pure python...")
        with open(os.path.join(wrk_dir, final_chaos), "w") as f_out:
            for file in file_list:
                with open(file, "r") as f_in:
                    for line in f_in:
                        f_out.write("{}".format(line))

# Crear archivo salida final
def make_sal(salida, suffix=""):
    # Ruta a la carpeta raíz que contiene las subcarpetas con los archivos
    root_dir = wrk_dir

    # Lista para almacenar los nombres de los archivos
    file_list = []

    # Recorre todas las subcarpetas en la carpeta raíz
    for subdir in os.listdir(root_dir):
        # Verifica si el nombre de la subcarpeta comienza con 'pref'
        if subdir.startswith(pref):
            subdir_path = os.path.join(root_dir, subdir)
            # Recorre todos los archivos en la subcarpeta
            for filename in os.listdir(subdir_path):
                # Verifica si el nombre del archivo comienza con
                # <salida> y termina con ".out"
                if (
                    filename.startswith(salida)
                    and filename.endswith(".out")
                    and (suffix in filename)
                ):
                    filepath = os.path.join(subdir_path, filename)
                    file_list.append(filepath)

    # Ordena los nombres de los archivos por el valor de i en "<salida>.out"
    if suffix == "":
        file_list = sorted(
            file_list, key=lambda x: int(x.split(salida)[1].split(".out")[0])
        )
    else:
        file_list = sorted(
            file_list,
            key=lambda x: int(
                x.split(salida)[1].split(".out")[0].split(suffix)[0]
            ),
        )

    # Concatena los archivos
    outs = salida.split(".")
    outs.insert(-1, suffix) if len(outs) > 1 else outs.insert(1, suffix)
    outs.insert(-1, ".out")
    salida = "".join(outs)
    command = ' '.join([file for file in file_list])
    try:
        p = subprocess.run(
            ["cat %s > %s" % (command, salida)],
            cwd=wrk_dir, check=True, shell=True
        )
    except:
        print("Could not create %s with cat."%salida)
        print(" Trying with pure python...")
        with open(os.path.join(wrk_dir, salida), "w") as f_out:
            for file in file_list:
                with open(file, "r") as f_in:
                    for line in f_in:
                        f_out.write("{}".format(line))



# Definir nombre único para wrk_dir
def generate_unique_dir(base_dir):
    i = 1
    unique_dir = base_dir
    while os.path.isdir(unique_dir):
        unique_dir = f"{base_dir}_{i}"
        i += 1
    return unique_dir


if __name__ == "__main__":
    if new_simulation and (not os.path.exists(wrk_dir)):
        os.mkdir(wrk_dir)  # Creamos directorio donde volcaremos todo
        print("Directory  '% s' created\n" % os.path.basename(wrk_dir))
    if not os.path.samefile(wrk_dir, cwd): # Pasamos todo al dir de trabajo (organizado)
        ## Programa
        nprogr = os.path.join(wrk_dir, program)
        subprocess.run(["cp", oprogr, nprogr], check=True)
        oprogr = nprogr
        ## Archivo de configuracion inicial
        if existe_ocini:
            ncini = os.path.join(wrk_dir, config)
            subprocess.run(["cp", ocini, ncini], check=True)
            ocini = ncini
        ## Partículas
        nparticles = os.path.join(wrk_dir, partfile)
        subprocess.run(["cp", oparticles, nparticles], check=True)
        oparticles = nparticles
        ## Archivo TOM
        if existe_otom:
            ntom = os.path.join(wrk_dir, tomfile)
            subprocess.run(["cp", otom, ntom], check=True)
            otom = ntom
    if not all_in_one:
        with ProcessPoolExecutor(max_workers=workers) as executor:
            results = executor.map(integrate_n, missing_lines)
        if datafile:
            print("")
            print("Creando archivo resumen %s.out\n"%datafile)
            if os.path.isfile(os.path.join(wrk_dir, "%s.out"%datafile)):
                print("WARNING: Se ha reemplazando archivo %s.out ya existente."%datafile)
            make_sal(datafile, suffix)
        if final_chaos:
            print("")
            print("Creando archivo resumen %s.out"%final_chaos)
            if os.path.isfile(os.path.join(wrk_dir, "%s.out"%final_chaos)):
                print("WARNING: Se ha reemplazando archivo %s.out ya existente."%final_chaos)
            make_sum(final_chaos, suffix)
    else:
        print("Running all systems in one process.")
        print("./%s %s" % (program, args))
        subprocess.run(
            ["./%s %s" % (program, args)], cwd=wrk_dir, check=True, shell=True
        )
    print("LISTO!")
