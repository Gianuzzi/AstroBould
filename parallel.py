# Version: 5.0

# Este programa integra cada sistema en un archivo de partículas
# y luego concatena los archivos de salida en un solo archivo.


# El archivo de partículas tiene formato:
# a e M w
# Pero también se puede introducir una última columna con el valor de R, 
#  el cual reemplaza a a: [a = R**(2/3.) * a_corot]. Entonces quedaría:
# a e M w R

## IMPORTANTE!!!!
# 1) Todos los archivos deben estar en la misma carpeta
# 2) El modo CHAOS debe estar ACTIVADO: [chaos = .TRUE.] (Ahora está automático en este script)
# 3) El ejecutable Debe tener DESactivados:
#   - El modo de salida pantalla: [screen = .FALSE.] (Ahora está automático en este script)
#   - El modo de cálculo de mapa de potencial: [map_pot = .FALSE.] (Ahora está automático en este script)

# En caso de dejar puesta la salida en archivo, se creará un archivo
#  llamado salida[id].dat por cada partícula.

# Finalmente se creará un archivo de caos llamado chaos[id].dat por cada partícula, 
#  y se concatenarán todos los archivos en un solo archivo <outfile> con el siguiente formato:
## 1    : numero simu
## 2    : mala? (0 == OK, 1 == Colisión, 2 == Escape)
## 3    : Momento angular total
## 4    : t total
## 5-9  : a ini, e ini, M ini, w ini, Res ini
## 10   : Momento angular inicial de partícula por unidad de masa
## 11   : t integrado
## 12-16: a final, e final, M final, w final, R final
## 17   : Momento angular final de partícula por unidad de masa
## 18-19: da, de 

# Excepto <outfile>, todos los otros archivos se encontrarán en carpetas creadas con el nombre
# dpy[pid], donde pid es el ID del procesador que ejecutó el sistema. El máximo de carpetas
# creadas será igual al MAX(número de procesadores disponibles en el sistema , workers).

# Si no son necesarias, se recomienda BORRAR las carpetas creadas luego de terminar la ejecución.
## Esto puede hacerse con: $ rm -rf dpy*



import os
import subprocess
from concurrent.futures import ProcessPoolExecutor as PPE

particles = "particles.in"   # Nombre del archivo de partículas
program = "main"             # Nombre del ejecutable
chaosfile = "chaos.dat"      # Nombre de cada archivo de caos (chaosfile en el programa)
datafile = ""                # Nombre de cada archivo de salida (datafile en el programa)  ["" si no se usa]
workers = 5                  # Número de procesadores a usar (workers)
suffix = ""                  # Suffix for the output files
outfile = "sump.out"         # Final Summary Output file name
exact = False                # Método: Exacto (cos, sin), o NO exacto (integra boulders y m0)



##### Iniciamos ####

# Obtener el path actual y renombrar
cwd = os.getcwd()
oparticles = os.path.join(cwd, particles)
oprogr = os.path.join(cwd, program)
ocini = os.path.join(cwd, "config.ini")

## Checkeamos los archivos
if not os.path.isfile(oparticles):
    msg = "Particles file {} does not exist.".format(oparticles)
    raise FileNotFoundError(msg)
if not os.path.isfile(oprogr):
    msg = "Executable file {} does not exist.".format(oprogr)
    raise FileNotFoundError(msg)
existe_ocini = os.path.isfile(ocini)
if not existe_ocini:
    print("WARNING: Configuration file {} does not exist.".format(ocini))
    print("         Se utilizarán los parámetros explicitados en el código, ")
    print("          en vez de los de algún archivo de parámetros. ")
    yes_no = input("¿Desea continuar? [y/n]")
    if yes_no.lower() not in ["y", "yes"]:
        print("Saliendo.")
        exit(1)
    

# Leemos el archivo de partículas
with open(oparticles, "r") as f:
    lines = f.readlines()

# Obtener el número de líneas del archivo de partículas
nsys = len(lines)
if nsys == 0:
    print("No hay partículas para integrar.")
    exit(1)
print("Cantidad total de partículas: {}".format(nsys))

# Obtener los sistemas realizados
if any([os.path.isdir(name) and name.startswith("dpy") for name in os.listdir(cwd)]):
    print("Checkeando integraciones ya completadas.")
    command = (
        f"find dpy* -name 'chaos*{suffix}.dat' "
        f"| sed -e 's/.*chaos\\([0-9]*\\){suffix}\\.dat/\\1/' "
        f"| sort -n"
        )
    result = subprocess.run(command, shell=True,
                            stdout=subprocess.PIPE, text=True)
    if result.returncode == 0:
        output_lines = result.stdout.splitlines()
    else:
        raise IOError("Error al leer integraciones ya realizadas.")
    done = set([int(cint) for cint in output_lines if cint!=''])
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

# Arreglamos por si hay "e" en vez de "d"
for i in range(nsys):
    lines[i] = lines[i].replace("e", "d")

## Argumentos. Estos son
args = "--noinfo --noscreen --nomap --nodatascr --noperc"
args += "%s"%(" -chaosfile %s"%chaosfile if chaosfile else "")
args += "%s"%(" -datafile %s"%datafile if datafile else " --nodata")
args += "%s"%(" --exact" if exact else " --noexact")

# Función general
def integrate_n(i):
    # Get processor ID
    PID = os.getpid()
    dirp = os.path.join(cwd, "dpy%d"%PID)
    nprogr = os.path.join(dirp, program)
    ncini = os.path.join(dirp, "config.ini")
    if not os.path.exists(dirp): # Si no existe el directorio
        subprocess.run(["mkdir", dirp], check=True)
        p = subprocess.run(["cp", oprogr, nprogr], check=True)
        if existe_ocini:
            p = subprocess.run(["cp", ocini, ncini], check=True)
    print("Running system %d\n"%(i))
    ### ESTO SE ESTÁ EJECUTANDO EN LA SHELL
    # print("Running: ./%s %s %d %s"%(program, args, i, lines[i]))
    p = subprocess.run(["./%s %s %d %s"%(program, args, i, lines[i])],
                         cwd=dirp, check=True, shell=True)
    subprocess.run(["mv", "-f",
                    os.path.join(dirp, chaosfile),
                    os.path.join(dirp, "chaos%d%s.dat"%(i, suffix))])
    if datafile != "":
        subprocess.run(["mv", "-f",
                        os.path.join(dirp, datafile),
                        os.path.join(dirp, "salida%d%s.dat"%(i, suffix))])
    return

def make_sum(outfile, suffix=""):
    # Ruta a la carpeta raíz que contiene las subcarpetas con los archivos
    root_dir = cwd

    # Lista para almacenar los nombres de los archivos
    file_list = []

    # Recorre todas las subcarpetas en la carpeta raíz
    for subdir in os.listdir(root_dir):
        # Verifica si el nombre de la subcarpeta comienza con "dpy"
        if subdir.startswith("dpy"):
            subdir_path = os.path.join(root_dir, subdir)
            # Recorre todos los archivos en la subcarpeta
            for filename in os.listdir(subdir_path):
                # Verifica si el nombre del archivo comienza con "chaos" y termina con ".dat"
                if filename.startswith("chaos") and filename.endswith(".dat") and (suffix in filename):
                    filepath = os.path.join(subdir_path, filename)
                    file_list.append(filepath)


    # Ordena los nombres de los archivos por el valor de i en "chaos%d%s.dat"
    if suffix == "":
        file_list = sorted(file_list, key=lambda x: int(x.split("chaos")[1].split(".dat")[0]))
    else:
        file_list = sorted(file_list, key=lambda x: int(x.split("chaos")[1].split(".dat")[0].split(suffix)[0]))

    # Concatena los archivos en el archivo de salida utilizando el comando 'cat' de Unix
    outs = outfile.split(".")
    outs.insert(-1, suffix) if len(outs) > 1 else outs.insert(1, suffix)
    outs.insert(-1, ".")
    outfile = "".join(outs)
    with open(outfile, "w") as f_out:
        for file in file_list:
            # # subprocess.run(["cat", file], stdout=f_out)
            # Extract the value of %d from the filename
            if suffix == "":
                i_val = int(file.split("chaos")[1].split(".")[0])
            else:
                i_val = int(file.split("chaos")[1].split(".")[0].split(suffix)[0])
            # Add the i value as the first column in the output file [NOT NOW]
            with open(file, "r") as f_in:
                for line in f_in:
                    # f_out.write("{}\t{}".format(i_val, line))
                    f_out.write("{}".format(line))


if __name__ == "__main__":
    with PPE(max_workers=workers) as executor:
        results = executor.map(integrate_n, missing_lines)
    if outfile:
        print("Creando archivo resumen {}".format(outfile))
        if os.path.isfile(outfile):
            print("WARNING: Se ha reemplazando archivo ya existente.")
        make_sum(outfile, suffix)
    print("LISTO!")
