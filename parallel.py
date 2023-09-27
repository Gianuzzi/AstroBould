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
# 2) El modo CHAOS debe estar ACTIVADO: [chaos = .TRUE.]
# 3) El ejecutable Debe tener DESactivados:
#   - El modo de salida pantalla: [screen = .FALSE.]
#   - El modo de cálculo de mapa de potencial: [map_pot = .FALSE.]

# En caso de dejar puesta la salida en archivo, se creará un archivo
#  llamado salida[id].dat por cada partícula.

# Finalmente se creará un archivo de caos llamado chaos[id].dat por cada partícula, 
#  y se concatenarán todos los archivos en un solo archivo <outfile> con el siguiente formato:
# id, bad?, tmax, x_0, y_0, vx_0, vy_0, R_0, t, da, de, a_0, e_0, M_0, w_0, R

# Excepto <outfile>, todos los otros archivos se encontrarán en carpetas creadas con el nombre
# dpy[pid], donde pid es el ID del procesador que ejecutó el sistema. El máximo de carpetas
# creadas será igual al MAX(número de procesadores disponibles en el sistema , workers).
# Si no son necesarias, se recomienda borrarlas luego de terminar la ejecución.

import os
import subprocess
from concurrent.futures import ProcessPoolExecutor as PPE

particles = "fort.67"     # Nombre del archivo de partículas
program = "main"          # Nombre del ejecutable
chaosfile = "chaos.dat"   # Nombre de cada archivo de caos (chaosfile en el programa)
datafile = "" # Nombre de cada archivo de salida (datafile en el programa)  ["" si no se usa]
workers = 120             # Número de workers
suffix = ""               # Suffix for the output files
outfile = "sump.out"      # Output file name


# Iniciamos
# Obtener el path actual y renombrar
cwd = os.getcwd()
oparticles = os.path.join(cwd, particles)
oprogr = os.path.join(cwd, program)
# Checkeamos los archivos
if not os.path.isfile(oparticles):
    msg = "Particles file {} does not exist".format(oparticles)
    raise FileNotFoundError(msg)
if not os.path.isfile(oprogr):
    msg = "Executable file {} does not exist".format(oprogr)
    raise FileNotFoundError(msg)
# Leemos el archivo de partículas
with open(oparticles, "r") as f:
    lines = f.readlines()
nsys = len(lines)
# p = subprocess.run(["wc", "-l", oparticles], capture_output=True, check=True)
# nsys = int(p.stdout.strip().split()[0])
print("Number of particles: {}".format(nsys))
# Obtener el número de workers
workers = min(max(1, min(int(workers), len(os.sched_getaffinity(0)))), nsys)
print("Workers: {}\n".format(workers))

for i in range(nsys):
    lines[i] = lines[i].replace("e", "d")

# Función general
def integrate_n(i):
    # Get processor ID
    PID = os.getpid()
    dirp = os.path.join(cwd, "dpy%d"%PID)
    nprogr = os.path.join(dirp, program)
    if not os.path.exists(dirp): # Si no existe el directorio
        subprocess.run(["mkdir", dirp], check=True)
        p = subprocess.run(["cp", oprogr, nprogr], check=True)
    print("Running system %d\n"%(i))
    # print("./%s"%program, "%d"%i, "%s"%lines[i]) ## ESTO SE ESTÁ EJECUTANDO EN LA SHELL
    p = subprocess.run(["./%s %d %s"%(program,i,lines[i])],
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
        results = executor.map(integrate_n, range(0, nsys))
    if outfile:
        make_sum(outfile, suffix)
