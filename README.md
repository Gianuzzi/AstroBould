# AstroBould
Repositorio con programa para integrar partículas gravitatoriamente alrededor de un asteroide con boulders

Los parámetros iniciales del sistema y la integración se modifican en el archivo: [config.ini](./config.ini)

En caso de no utilizarlo, se usan los parámetros definidos en [main.F90](./src/main.F90), entre las líneas [22 a 118](./src/main.F90#L22#L118).

# Usage

## Compile
``` console
make
```

Recuerda volver a ejecutar `make` cada vez que se realice un cambio en algún .F90.

Asimismo, recuerda ejecutar `make clean` en caso de realizar algún cambio de compilación.

## Run

La ejecución básica de una particula individual se realiza con

``` console
$ ./main <a> <e> <M> <w> <R> [args]
```

Los argumentos $a, e, M, \omega$ son los elementos orbitales, y $R$ (cociente de movimientos inicial $n_{part}/\Omega_{ast}$) es opcional.
En caso de usarlo, calcula y reemplaza el valor de $a$.

Una ayuda puede obtenerse con 

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

### Cofigurations file

Para definir los parámetros de una integración a realizar, es posible editar el [archivo de configuraciones](./config.ini). Este archivo será leído por el ejecutable al comenzar la ejecución.

Es importante tener en cuenta que algunas de las configuraciones aquí seteadas serán anuladas y sobreescritas por las introducidas en [launcher.py](./launcher.py), en caso de utilizar este segundo método para integración en paralelo (ver abajo).


### Parallel 

Hay 2 posibilidad de ejecución en paralelo:
     
- _Dependiente_: Se realiza una sola integración, incluyendo un asteroide y muchas partículas al mismo tiempo. Aquí se paralelizan las tareas relacionadas a cada partícula por separado (i.e, cálculo de aceleraciones, elementos orbitales, etc.).
- _Independiente_: Se realizan múltiples integraciones independientes, con una partícula cada una. Aquí cada integración se ejecuta en paralelo.

En este caso se debe tener un **archivo con los datos de las partículas a integrar**.
Este archivo se puede crear con el código [make_particles.py](./tools/make_particles.py). Su configuración se realiza en las líneas [70 a 104](./tools/make_particles.py#L70#L104), y si ejecució se realiza con:

``` console
$ python make_particles.py
```

Para la integración _dependiente_ es necesario compilar usando _-fopenmp_, mientras que para la independiente (default) no es necesario. Para el primer caso, se debe compilar ejecutando

``` console
$ make parallel
```

Para realizar la integración _independiente_, se provee el código en Python [launcher.py](./launcher.py). Más información se encuentra [al inicio](./launcher.py#L3#L47) de este archivo. Luego de configurarlo, editando las líneas [61 a 90](./launcher.py#L61#L90), su ejecución se realiza con:

``` console
$ python launcher.py
```

(Recordar estar en algún entorno de python)

Para realizar una integración _dependiente_, hay dos posibilidades: 
1) Editar el número de cpus en el [archivo de configuraciones](./config.ini) ([use parallel threads](./config.ini#L17)), o ejecutar directamente en la terminal:

``` console
$ ./ASTROBOULD [args] -parallel <number of cpus to use>
```

2) Utilzar el código [launcher.py](./launcher.py) para esta integración, seteando

``` python
all_in_one = True
```
en la línea [90](./launcher.py#L90) de este archivo.


### Archivo TOM (Tiempos, Omega, Masa agregada)

Para crear este archivo se puede usar [make_TOM.py](./tools/make_TOM.py), luego de configurarlo (líneas [18 a 42](./tools/make_TOM.py#L18#L42)). Se puede encontrar más información [al inicio](./tools/make_TOM.py#L6#L16) de este archivo. 

Por otro lado, también  se puede usar [make_TOM.f90](./tools/make_TOM.f90), de C. Beaugé. Se debe configurar manualmente, para luego compilar y ejecutar.

# Author
Emmanuel Gianuzzi

# Notes
This code is **under development**.


