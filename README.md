# AstroBould

**AstroBould** is a code designed to integrate particles that gravitationally interact with an asteroid featuring mass anomalies.

The initial configuration and integration parameters can be set in the file: [config.ini](./config.ini)

If this file is not used, the default parameters are those defined in [main.F90](./src/main.F90), between lines [24 and 134](./src/main.F90#L24-L134).

---

# Usage

## 🧑‍💻 Compiling

To compile the code, run:

```bash
make
```

> ⚠️ **Note**: Remember to run `make clean` followed by `make` every time you modify any `.F90` source file.

Default copilator is GNU (_gfortran_). If available, the compilation could be made using the IntelFortranCompiler (_ifx_), by running

```bash
make IFORT=1
```

## 🏃🏼 Running

To execute the code with a single particle, use the following command:

```bash
./ASTROBOULD <a> <e> <M> <w> <R> [args]
```

Here:

- `<a>` = semi-major axis
- `<e>` = eccentricity
- `<M>` = mean anomaly
- `<w>` = argument of periapsis
- `<R>` = initial spin-orbit ratio \( n_\text{part} / \Omega_\text{ast} \) (optional)

> If `<R>` is provided, it overrides the value of `<a>`.

## 🆘 Help
For basic help in Spanish just run:

``` console
$ ./ASTROBOULD --help

 Uso: ./ASTROBOULD <ea> <ee> <eM> <ew> <eR> [args]
     ea  : Elemento a de la partícula (km)
     ee  : Elemento e de la partícula
     eM  : Elemento M de la partícula (deg)
     ew  : Elemento w de la partícula (deg)
     eR  : Elemento R de la partícula [Optional]
     -mpart       : Masa de la partícula individual
     -nsim        : Asignar como número de simulación al 'int' que sigue
     -datafile    : Guardar datos en el archivo que sigue
     --nodataf    : No guardar datos
     -chaosfile   : Guardar caos en el archivo que sigue
     --nochaosf   : No guardar caos
     --screen     : Imprimir información en pantalla
     --noscreen   : No imprimir en pantalla
     --perc       : Imprimir porcentaje de integración
     --noperc     : No imprimir porcentaje de integración
     --datascr    : Imprimir datos en pantalla
     --nodatascr  : No imprimir datos en pantalla
     -multifile   : Guardar datos en archivos individuales
     --nomultif   : No guardar datos en archivos individuales
     -mapfile     : Guardar mapas de potencial en el archivo que sigue
     --nomapf     : No guardar mapas de potencial
     --implicit   : Usar método implícito (integra) [default]
     --explicit   : Usar método explícito (cos, sen)
     --elem       : Imprimir elementos orbitales (solo partículas) [default]
     --noelem     : Imprimir coordenadas baricéntricas
     -tomfile     : Utilizar archivo de (t)iempos|omega|masa que sigue
     --notomfile  : No utilizar archivo de (t)iempos|omega|masa
     -partfile    : Utilizar archivo de partículas que sigue
     --nopartfile : No utilizar archivo de partículas
     --noconfig   : No leer archivo de configuración
     --version1   : Usar versión 1: y=(B0, B1, ..., P0, ...)
     --version2   : Usar versión 2: y=(θ, ω, P0, ...) [default]
     --merge      : Incluir asociar colisiones de partículas a asteroide
     --nomerge    : No asociar colisiones de partículas
     --torque     : Incluir torque de partículas hacia el asteroide
     --notorque   : No incluir torque de partículas hacia el asteroide
     -parallel    : Paralelizar usando la cantida de thread que sique
     --parallel   : Paralelizar usando todos los threads disponibles
     --noparallel : No usar paralelización para partículas
     --help       : Mostrar esta ayuda

``` 

## 🛠️ Configuration

To define the parameters of an integration, you can edit the [configuration file](./config.ini). This file is read by the executable at runtime.

> ⚠️ Note: Some settings defined in `config.ini` may be overridden by those specified in [launcher.py](./launcher.py) if using parallel execution via Python (see below).

## ⛓️ Parallel Execution

There are two modes for parallel execution:
- Dependant 
- Independant (default)

⌨️ **Particles Input File**

Both modes require the existence of a particles file (e.g. _particles.in_) containing all particles (initial conditions) to be integrated. You can generate this file using [make_particles.py](./tools/make_particles.py). Configuration is found between lines [87 and 93](./tools/make_particles.py#L87#L93).

### 1. **Dependent Mode**

- A single integration is performed, including the asteroid and many particles simultaneously.
- Tasks related to each particle (e.g., force calculations, orbital elements, etc.) are parallelized internally.
- Requires OpenMP support. Compile using:
```console
$ make parallel
```
To run:
```console
$ ./ASTROBOULD [args] -parallel <number_of_cpus> -partfile <particles_file>
```
or edit [config.ini](./config.ini) and set "use parallel threads" to the desired value (see line [17](./config.ini#L17)), and "particles input file" to the particles file name (see line [58](./config.ini#L58)). Then simply run
```console
$ ./ASTROBOULD [args]
```

### 2. **Independent Mode (default)**

- Multiple independent integrations are performed, one particle per run.
- Each integration is executed in parallel (e.g., across multiple CPU cores).
- Does not require _-fopenmp_.

The file [launcher.py](./launcher.py) provides all available configurations for this parallel execution mode. Configure the run by editing lines [60 to 95](./launcher.py#L60#L95). More information (in spannish) is available at the top of the file(see lines [3 to 48](./launcher.py#L3#L48)).

To run:
```console
$ python launcher.py
```
💡 Make sure to run inside a Python virtual environment. 🐍
> ⚠️ Make sure that ```all_in_one = False``` in [launcher.py](./launcher.py#L92) for this mode.



# Author
Emmanuel Gianuzzi

# Notes
This code is **under development**.


