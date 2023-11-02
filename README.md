# AstroBould
Repositorio con programa para integrar partículas gravitatoriamente alrededor de un asteroide con boulders

Los parámetros iniciales del sistema y la integración se modifican en [main.F90](./main.F90#L90#L122), entre las líneas 90 y 122.

COMPILAR:
``` console
make
```

Uso: 

``` console
  ./main [nsim] [ea] [ee] [eM] [ew] [eR] [args]
   nsim: Número de simulación
   ea  : Elemento a de la partícula
   ee  : Elemento e de la partícula
   eM  : Elemento M de la partícula
   ew  : Elemento w de la partícula
   eR  : Elemento R de la partícula
   --nodata   : No guardar datos
   -datafile  : Guardar datos en el archivo que sigue
   --noinfo   : No guardar información
   -infofile  : Guardar información en el archivo que sigue
   --nochaos  : No guardar caos
   -chaosfile : Guardar caos en el archivo que sigue
   --screen   : Imprimir información en pantalla
   --noscreen : No imprimir en pantalla
   --perc     : Imprimir porcentaje de integración
   --noperc   : No imprimir porcentaje de integración
   --datascr  : Imprimir datos en pantalla
   --nodatascr: No imprimir datos en pantalla
   --nomap    : No guardar mapas de potencial
   -mapfile   : Guardar mapas de potencial en el archivo que sigue
   --noexact  : No usar método exacto
   --exact    : Usar método exacto
   --elem     : Imprimir elementos orbitales (solo partícula)
   --noelem   : Imprimir coordenadas baricéntricas
```

PARTÍCULA INDIVIDUAL: 
``` console
$ ./main 0 a e M w R
```
0 (o cualquier otro "int") es obligatorio

"a, e, M, w" son los elementos orbitales.

R es opcional. En caso de usarlo, calcula y reemplaza el valor de a.

PARALELO: (Recordar estar en algún entorno)
En este caso se debe tener un archivo con los datos de las partículas. La info está al inicio del archivo [parallel.py](./parallel.py#L3#L30))
``` console
$ python parallel.py
```
