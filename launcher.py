# Version: 7.1

# DENTRO DE UN ENTORNO PYTHON
# Ejecución: $ python parallel.py

# Este programa integra cada sistema (condición inicial) en
# un archivo independiente,
# y luego concatena los archivos de caos en un solo archivo.

# El archivo de partículas a introducir debe tener formato:
# a e M w
# Pero también se puede introducir una última columna con el valor de R,
# el cual reemplaza a a: [a = R**(2/3.) * a_corot]. Entonces quedaría:
# a e M w R

# IMPORTANTE : Todos los archivos deben estar en la misma carpeta

# En caso de dejar puesta la salida en archivo (<datafile>),
# se creará un archivo llamado salida[id].dat por cada partícula.

# Este código crea un archivo de caos (chaos[id].dat) por cada partícula.
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
import subprocess
from concurrent.futures import ProcessPoolExecutor

# Editar estas líneas según corresponda.
# Estos valores tienen privilegio ante los de config.ini.

# Parámetros #
program = "ASTROBOULD"  # Nombre del ejecutable

## Configuración de integración ##
workers = 4  # Número de procesadores a usar (workers)
explicit = False  # Método: True (cos, sin), False (integra boulders y m0)
version1 = False  # Versión del código (1 o 2)

## Torque and merge ##
torque = False  # Si se quiere usar torque
merge = False  # Si se quiere usar merge

## Input ## ("" o False si no se usa)
config = "config.ini"  # Nombre del archivo de configuración
partfile = "particles.in"  # Nombre del archivo de partículas
tomfile = "" # Nombre del archivo de valores de t_i, delta_omega(t_i), y delta_masa(t_i)

## Output ## ("" o False si no se usa)
datafile = "salida" # Nombre del archivo de salida de datos (sin extensión)
final_chaos = "sump.out"  # Final Chaos Summary Output file name
suffix = ""  # Suffix for the output files
### Screen ###
screen_info = True  # Si se quiere ver la información en pantalla (solo con torque = True)
screen_data = False  # Si se quiere ver los datos en pantalla # "%" para porcentaje (solo con torque = True)
### Elements ###
elements = True  # Si se quiere devolver elementos orbitales (en datafile)



# ----------------------------------------------------------------------
# -------------------- No tocar de aquí en adelante --------------------
# ----------------------------------------------------------------------

# Iniciamos #

# Obtener el path actual
cwd = os.getcwd()

# Archivos
ocini = os.path.join(cwd, config) # Archivo de configuración
oprogr = os.path.join(cwd, program) # Ejecutable
oparticles = os.path.join(cwd, partfile) # Archivo de partículas
otom = os.path.join(cwd, tomfile) # Archivo de valores de t_i, delta_omega(t_i), y delta_masa(t_i)

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
        exit(1)
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
    exit(1)
## Datafile
if os.path.isfile(final_chaos):
    print("WARNING: Output file {} already exist.".format(final_chaos))
    yes_no = input("Do you want to overwrite it? y/[n]\n")
    if yes_no.lower() not in ["y", "yes", "s", "si"]:
        i = 1
        aux = final_chaos.split(".")
        suf = aux[-1] if len(aux) > 1 else ""
        while os.path.isfile(final_chaos):
            final_chaos = ".".join(aux[:-1]) + str(i) + ("." + suf if suf else "")
            i += 1

# Checks #
## Si hay torque, entonces explicit debe ser false
if torque and explicit:
    print("WARNING: Torque is active. Explicit mode will be deactivated.")
    explicit = False
    yes_no = input("Do you want to continue? y/[n]\n")
    if yes_no.lower() not in ["y", "yes", "s", "si"]:
        print("Saliendo.")
        exit(1)


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
    print("No hay partículas para integrar.")
    exit(1)
print("Cantidad total de partículas: {}".format(nsys))

if not torque:
    # Definimos prefijo de directorio, de acuerdo a TOMfile o no
    pref = "tomd" if tomfile else "dpy"

    # Obtener los sistemas realizados
    if any(
        [os.path.isdir(name) and name.startswith(pref) for name in os.listdir(cwd)]
    ):
        print("Checkeando integraciones ya completadas.")
        command = (
            f"find {pref}* -name 'chaos*{suffix}.dat' "
            f"| sed -e 's/.*chaos\\([0-9]*\\){suffix}\\.dat/\\1/' "
            f"| sort -n"
        )
        result = subprocess.run(
            command, shell=True, stdout=subprocess.PIPE, text=True
        )
        if result.returncode == 0:
            output_lines = result.stdout.splitlines()
        else:
            raise IOError("Error al leer integraciones ya realizadas.")
        done = set([int(cint) for cint in output_lines if cint != ""])
        missing_lines = [x for x in range(nsys) if x not in done]
        print("   Cantidad de sistemas ya integrados: {}".format(len(done)))
        nsys = len(missing_lines)
        print("   Cantidad de sistemas a integrar: {}".format(nsys))
    else:
        missing_lines = range(0, nsys)

    # Hay que hacer?
    if len(missing_lines) == 0:
        print("Ya se han integrado todas las partículas.")
        exit(1)


# Obtener el número de workers
workers = min(max(1, min(int(workers), len(os.sched_getaffinity(0)))), nsys)
print("Workers: {}\n".format(workers))

# Argumentos. Estos son:
if torque:
    args = "--torque --nomapf"
    args += " --screen" if screen_info else " --noscreen"
    if screen_data=="%":
        args += " --perc --nodatascr"
    else:
        args += " --datascr" if screen_data else " --nodatascr"
        args += " --noperc"        
    args += " -partfile %s" % partfile
    args += " -datafile %s.dat" % datafile if datafile else " --nodataf"
    args += " -parallel %d" % workers
else:
    args = " --notorque  --nomapf --noscreen --nodatascr --noperc --noparallel"
args += " --explicit" if explicit else " --implicit"
args += " -tomfile %s" % tomfile if tomfile else " --notomfile"
args += " --elem" if elements else " --noelem"
args += " --version1" if version1 else " --version2"
args += " --merge" if merge else " --nomerge"



# Función general
def integrate_n(i):
    # Get processor ID
    pid = os.getpid()
    dirp = os.path.join(cwd, "%s%d" % (pref, pid))
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
    this_args += " -chaosfile chaos%d%s.dat" % (i, suffix)
    this_args += "%s" % (" --nodataf"
        if not datafile else
        " -datafile %s" % ("%s%d%s.dat" % (
            datafile if isinstance(datafile,str) else "salida", i, suffix)
            )
    )
    print("Running system %d\n" % (i))
    my_line = " ".join(lines[i].split()[1:]) # Elimina el primer valor (en caso que fuese la masa...)
    # ESTO SE ESTÁ EJECUTANDO EN LA SHELL #
    # print("Running: ./%s %s %s %s"%(program, args, this_args, my_line))
    # (Lines debe ser último porque termina en "\n") #
    p = subprocess.run(
        ["./%s %s %s %s" % (program, args, this_args, my_line)],
        cwd=dirp,
        check=True,
        shell=True,
    )
    if p.returncode != 0:
        print("The system %d has failed." % i)
        return
    print("System %d has been integrated." % i)
    return


def make_sum(final_chaos, suffix=""):
    # Ruta a la carpeta raíz que contiene las subcarpetas con los archivos
    root_dir = cwd

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
                # "chaos" y termina con ".dat"
                if (
                    filename.startswith("chaos")
                    and filename.endswith(".dat")
                    and (suffix in filename)
                ):
                    filepath = os.path.join(subdir_path, filename)
                    file_list.append(filepath)

    # Ordena los nombres de los archivos por el valor de i en "chaos%d%s.dat"
    if suffix == "":
        file_list = sorted(
            file_list, key=lambda x: int(x.split("chaos")[1].split(".dat")[0])
        )
    else:
        file_list = sorted(
            file_list,
            key=lambda x: int(
                x.split("chaos")[1].split(".dat")[0].split(suffix)[0]
            ),
        )

    # Concatena los archivos
    outs = final_chaos.split(".")
    outs.insert(-1, suffix) if len(outs) > 1 else outs.insert(1, suffix)
    outs.insert(-1, ".")
    final_chaos = "".join(outs)
    with open(final_chaos, "w") as f_out:
        for file in file_list:
            with open(file, "r") as f_in:
                for line in f_in:
                    f_out.write("{}".format(line))


if __name__ == "__main__":
    if not torque:
        with ProcessPoolExecutor(max_workers=workers) as executor:
            results = executor.map(integrate_n, missing_lines)
        if final_chaos:
            print("Creando archivo resumen {}".format(final_chaos))
            if os.path.isfile(final_chaos):
                print("WARNING: Se ha reemplazando archivo ya existente.")
            make_sum(final_chaos, suffix)
    else:
        print("Running all systems in one process.")
        subprocess.run(["./%s %s" % (program, args)], cwd=cwd, check=True, shell=True)
    print("LISTO!")
