# AstroBould
Repositorio con programa para integrar partículas gravitatoriamente alrededor de un asteroide con boulders

Los parámetros iniciales del sistema y la integración se modifican en el archivo: [config.ini](./config.ini)

En caso de no utilizarlo, se usan los parámtros definidos en [main.F90](./main.F90#18#89), entre las líneas 18 y 89.

# Usage

## Compile
``` console
make
```

Recuerda volver a ejecutar `make`  cada vez que se realice un cambio en algún .F90.

## Run

La ejecución básica de una particula individual se realiza con

``` console
$ ./main <a> <e> <M> <w> <R> [args]
```

Los argumentos "$a$, $e$, $M$, $\omega$" son los elementos orbitales, y $R$ (cociente inicial $n_{part}/\Omega_{ast}$) es opcional.
En caso de usarlo, calcula y reemplaza el valor de a.

Una ayuda puede obtenerse con 

``` console
 $ ./main --help

  Uso: ./main <ea> <ee> <eM> <ew> <eR> [args]
  ea   : Elemento a de la partícula (km)
  ee   : Elemento e de la partícula
  eM   : Elemento M de la partícula (deg)
  ew   : Elemento w de la partícula (deg)
  eR   : Elemento R de la partícula [Opcional]
  -nsim       : Asignar como número de simulación al 'int' que sigue
  --nodata    : No guardar datos
  -datafile   : Guardar datos en el archivo que sigue
  --noinfo    : No guardar información
  -infofile   : Guardar información en el archivo que sigue
  --nochaos   : No guardar caos
  -chaosfile  : Guardar caos en el archivo que sigue
  --screen    : Imprimir información en pantalla
  --noscreen  : No imprimir en pantalla
  --perc      : Imprimir porcentaje de integración
  --noperc    : No imprimir porcentaje de integración
  --datascr   : Imprimir datos en pantalla
  --nodatascr : No imprimir datos en pantalla
  --nomap     : No guardar mapas de potencial
  -mapfile    : Guardar mapas de potencial en el archivo que sigue
  --implicit  : Usar método implícito (integra boulders)
  --explicit  : Usar método explícito (cos, sen)
  --elem      : Imprimir elementos orbitales (solo partícula)
  --noelem    : Imprimir coordenadas baricéntricas
  -tomfile    : Utilizar archivo de (t)iempos|omega(t)|dmasa(t) que sigue
  --notomfile : No utilizar archivo de (t)iempos|omega(t)|dmasa(t)
```


### PARALELO (Recordar estar en algún entorno de python)
En este caso se debe tener un **archivo con los datos de las partículas a integrar**.
Este archivo se puede crear con el código [make_particles.py](./make_particles.py). Su configuración se realiza en las líneas [67 a 97](./make_particles.py#L67#L97).

Para ejecutar en paralelo, se puede utilizar el código [parallel.py](./parallel.py). Más información se encuentra [al inicio](./parallel.py#L2#L55) de este archivo. Luego de configurarlo, su ejecución se realiza con

``` console
$ python parallel.py
```

### Archivo TOM (Tiempos, Omega, Masa agregada)

Para crear este archivo se puede usar [make_TOM.py](./make_TOM.py), luego de configurarlo (líneas [17 a 41](./make_TOM.py#L17#L41)). Se puede encontrar más información [al inicio](./make_TOM.py#L5#L15) de este archivo. 

Por otro lado, también  se puede usar [make_TOM.f90](./make_TOM.f90), de C. Beaugé. Se debe configurar manualmente, para luego compilar y ejecutar.

# Author
Emmanuel Gianuzzi

# Notes
This code is **under development**.


