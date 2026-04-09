# AstroBould

[![License](https://img.shields.io/badge/License-MIT-blue.svg)](https://tldrlegal.com/license/mit-license)

**AstroBould** is a code designed to integrate particles that gravitationally interact with an asteroid featuring mass anomalies or a central tri-axial ellipsoid.

The initial configuration and integration parameters can be set in the file: [config.ini](./config.ini)

If this file is not used, the default parameters are those defined in [main.F90](./src/main.F90), between lines [39 and 227](./src/main.F90#L39#L227).

---

# Usage

## 🧑‍💻 Compiling

To compile the code, run:

```bash
make
```

> ⚠️ **Note**: Remember to run `make clean` followed by `make` every time you modify any `.F90` source file.

Default copilator is GNU (_gfortran_). If available, the compilation could be made using the IntelFortranCompiler (_ifx_) or AMDFortranCompiler (_flang_), by running

```bash
make intel
```

or


```bash
make amd
```

## 🏃🏼 Running

To execute the code with a single particle, use the following command:

```bash
./ASTROBOULD <a> <e> <M> <w> <mmr> [args]
```

Here:

- `<a>`   = semi-major axis of the particle (km)
- `<e>`   = eccentricity of the particle
- `<M>`   = mean anomaly of the particle (degrees)
- `<w>`   = argument of periapsis of the particle (degrees)
- `<mmr>` = initial spin-orbit ratio ($$\Omega_\text{ast}  / n_\text
{part}$$) (optional)

> If `<mmr>` is provided, it overrides the value of `<a>`.

## 🆘 Help
For basic help in Spanish just run:

``` console
$ ./ASTROBOULD --help

 Uso: .\ASTROBOULD <ea> <ee> <eM> <ew> <mmr> [args]"
    ea  : Elemento a de la partícula/luna (km)
    ee  : Elemento e de la partícula/luna
    eM  : Elemento M de la partícula/luna (deg)
    ew  : Elemento w de la partícula/luna (deg)
    mmr : Valor de MMR de la partícula/luna [Opcional. 0 si no se utiliza.]
    --onlyprint   : No integrar; solo imprimir configuraciones
    -nsim         : Número de simulación [int]
    -mumoon       : Cociente de masa entre la luna individual y el asteroide
    -radius       : Radio de la partícula/luna individual (km)
    -datafile     : Nombre de archivo de salida de datos
    --nodataf     : No guardar datos de salida
    -chaosfile    : Nombre de archivo de salida caos
    --nochaosf    : No guardar salida de caos
    -geomfile     : Nombre de archivo de salida de elementos geométricos
    --nogeomf     : No guardar salida de elementos geométricos
    -filtfile     : Prefijo a agregar en archivos con filtro (y activar)
    --nofilter    : No utilizar filtro
    --screen      : Imprimir información en pantalla
    --noscreen    : No imprimir en pantalla
    --perc        : Imprimir porcentaje de integración
    --noperc      : No imprimir porcentaje de integración
    --datascr     : Imprimir datos de salida en pantalla
    --nodatascr   : No imprimir datos de salida en pantalla
    --diagnostic  : Imprimir datos de diagnostico en pantalla
    --nodiagnostic: No imprimir datos de diagnostico en pantalla
    -multifile    : Nombre base de archivo de salida de datos individuales
    --nomultif    : No guardar datos en archivos individuales
    -mapfile      : Nombre de archivo de mapa
    --nomapf      : No guardar mapa de potencial
    --elem        : Imprimir elementos orbitales (lunas/partículas) [default]
    --noelem      : Imprimir coordenadas baricéntricas
    -tomfile      : Nombre de archivo de (t)iempos|omega|masa a utilizar
    --notomfile   : No utilizar archivo de (t)iempos|omega|masa
    -moonfile     : Nombre de archivo de lunas a utilizar
    --nomoonfile  : No utilizar archivo de lunas
    -partfile     : Nombre de archivo de partículas a utilizar
    --nopartfile  : No utilizar archivo de partículas
    --noconfig    : No leer archivo de configuración
    -merge        : Tipo de colisiones (merges) permitidas [int]: 
                    0: Ninguno, 1: Partícula-Masivo, 2: Masivo-Masivo, 3: Todos
    -stopif       : Detener la integración sin más objetos del tipo [int]:
                    0: No detener, 1: Luna, 2: Partícula, 3: Ambos
    --megno       : Calcular MEGNO para partículas
    --nomegno     : No calcular MEGNO
    --sinodic     : Integar en sistema rotante
    --nosinodic   : Integar en sistema NO rotante
    -parallel     : Cantida de thread a utilizar en paralelo [int]
    --parallel    : Paralelizar usando todos los threads disponibles
    --noparallel  : No usar paralelización para lunas/partículas
    --version     : Imprimir versión del código al iniciar la ejecución
    --help        : Mostrar esta ayuda

``` 

## ⚙️ Configuration

To define the parameters of an integration, you can edit the [configuration file](./config.ini). This file is read by the executable at runtime.

> ⚠️ Note: Some settings defined in `config.ini` may be overridden by those specified in `launcher.py` if using parallel execution via Python (see below).

## ⛓️ Parallel Execution

There are two modes for parallel execution:
- Dependant 
- Independant (default)

⌨️ **Particles Input File**

Both modes require the existence of a particles/moons file (e.g. _particles.in_) containing all bodies (initial conditions) to be integrated in parallel. You can generate this file using [make_particles.py](./tools/make_particles.py). Configuration is found between lines [78 and 118](./tools/make_particles.py#L78#L118).

### 1. **Independent Mode**

- Multiple independent integrations are performed, one particle/moon per run.
- Each integration is executed in parallel (e.g., across multiple CPU cores).
- Does not require _-fopenmp_.

The file [launcher.py](./launcher.py) provides all available configurations for this parallel execution mode. Configure the run by editing lines [67 to 116](./launcher.py#L67#L116). More information (in spannish) is available at the top of the file(see lines [1 to 54](./launcher.py#L1#L54)).

To run:
```console
$ python launcher.py
```
💡 Make sure to run inside a Python virtual environment. 🐍

### 2. **~~Dependent Mode~~ (DEPRECATED)**

- A single integration is performed, including the asteroid and many particles/moons simultaneously.
- Tasks related to each particle (e.g., force calculations, orbital elements, etc.) are parallelized internally.
- Requires OpenMP support (_-fopenmp_). Compile using:
```console
$ make parallel
```
To run:
```console
$ ./ASTROBOULD [args] -parallel <number_of_cpus> -partfile <particles_file>
```
or edit [config.ini](./config.ini) and set "use parallel threads" to the desired value (see line [16](./config.ini#L16)), and "particles input file" to the particles file name (see line [71](./config.ini#L71)). Then simply run
```console
$ ./ASTROBOULD [args]
```


## 🔗 Integrators avaiable

The directory [integrators](./src/integrators/) includes the integrators modules. A brief [README](./src/integrators/README.md) can be found with minimal explanation of available integrators and properties.

## 🛠️ Additional tools

The directory [tools](./tools/) includes some extra scrips that could be useful for the user, such as:

- [get_omega_damp.py](./tools/get_omega_damp.py): Calculate $\tau_\Omega$ for $\Omega$ damping.
- [make_particles.py](./tools/make_particles.py)
: Create a particles/moons file usable as input for the code.
- [make_cj.py](./tools/make_cj.py)
: Create a particles file with fixed Jacobi constant usable as input for the code.
- [read_data.py](./tools/read_data.py): Functions to read the output generated by the code.
- [plot_salida.ipynb](./tools/plot_salida.ipynb): Notebook to plot outputs.
- [plot_chaos.ipynb](./tools/plot_chaos.ipynb): Notebook to plot chaos output.
- [plot_map.ipynb](./tools/plot_map.ipynb): Notebook to plot potential and accceleration map output.
- [plot_surface.ipynb](./tools/plot_surface.ipynb): Notebook to plot surface section output.
- ~~[make_TOM.py](./tools/make_TOM.py)~~: Generate a TOM file from output code. [**Not available yet**]

# Author
Emmanuel Gianuzzi [egianuzzi@unc.edu.ar](egianuzzi@unc.edu.ar) ([IATE-OAC-CONICET][], [FaMAF-UNC][])


  [IATE-OAC-CONICET]: http://iate.oac.uncor.edu/
  [OAC-CONICET]: https://oac.unc.edu.ar/
  [FaMAF-UNC]: https://www.famaf.unc.edu.ar/

# Notes
This code is **under development**.


