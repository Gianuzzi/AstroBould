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

Los argumentos "a, e, M, w" son los elementos orbitales, y R (cociente n_part/Omega_ast) es opcional.
En caso de usarlo, calcula y reemplaza el valor de a.

Una ayuda puede obtenerse con 

``` console
 $ ./main --help

  ./main <ea> <ee> <eM> <ew> <eR> [args]
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
  -tomfile    : Utilizar archivo de (t)iempos|omega|masa que sigue
  --notomfile : No utilizar archivo de (t)iempos|omega|masa
```


### PARALELO (Recordar estar en algún entorno de python)
En este caso se debe tener un archivo con los datos de las partículas a integrar.
Este archivo se puede crear con el código [make_particles.py](./make_particles.py)

Se puede encontrar más información al inicio del archivo [parallel.py](./parallel.py#L2#L55)

Luego de configurarlo, su ejecución se realiza con

``` console
$ python parallel.py
```

# Author
Emmanuel Gianuzzi

# Notes
This code is **under development**.


