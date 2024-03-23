# AstroBould
Repositorio con programa para integrar partículas gravitatoriamente alrededor de un asteroide con boulders

Los parámetros iniciales del sistema y la integración se modifican en el archivo: [config.ini](./config.ini)

En caso de no utilizarlo, se usan los parámtros definidos en [main.F90](./main.F90#18#89), entre las líneas 18 y 89.

COMPILAR:
``` console
make
```

Uso: 

``` console
  ./main [nsim] [ea] [ee] [eM] [ew] [eR] [args]
  nsim : Número de simulación"
  ea   : Elemento a de la partícula (km)"
  ee   : Elemento e de la partícula"
  eM   : Elemento M de la partícula (deg)"
  ew   : Elemento w de la partícula (deg)"
  eR   : Elemento R de la partícula"
  --nodata    : No guardar datos"
  -datafile   : Guardar datos en el archivo que sigue"
  --noinfo    : No guardar información"
  -infofile   : Guardar información en el archivo que sigue"
  --nochaos   : No guardar caos"
  -chaosfile  : Guardar caos en el archivo que sigue"
  --screen    : Imprimir información en pantalla"
  --noscreen  : No imprimir en pantalla"
  --perc      : Imprimir porcentaje de integración"
  --noperc    : No imprimir porcentaje de integración"
  --datascr   : Imprimir datos en pantalla"
  --nodatascr : No imprimir datos en pantalla"
  --nomap     : No guardar mapas de potencial"
  -mapfile    : Guardar mapas de potencial en el archivo que sigue"
  --implicit  : Usar método implícito (integra boulders)"
  --explicit  : Usar método explícito (cos, sen)"
  --elem      : Imprimir elementos orbitales (solo partícula)"
  --noelem    : Imprimir coordenadas baricéntricas"
  -tomfile    : Utilizar archivo de (t)iempos|omega|masa que sigue"
  --notomfile : No utilizar archivo de (t)iempos|omega|masa"
```

PARTÍCULA INDIVIDUAL: 
``` console
$ ./main 0 <a> <e> <M> <w> <R> [args]
```
El primer 0 (o cualquier otro "int") es obligatorio

"a, e, M, w" son los elementos orbitales.

R es opcional. En caso de usarlo, calcula y reemplaza el valor de a.

PARALELO: (Recordar estar en algún entorno de python)
En este caso se debe tener un archivo con los datos de las partículas. La info está al inicio del archivo [parallel.py](./parallel.py#L2#L55)
``` console
$ python parallel.py
```
