# AstroBould
Repositorio con programa para integrar partículas gravitatoriamente alrededor de un asteroide con boulders

Los parámetros iniciales del sistema y la integración se modifican en [main.F90](./main.F90#L90#L122), entre las líneas 90 y 122.

COMPILAR:
``` console
$ gfortran -O3 const.F90 run.F90 parameters.F90 integrators.F90 stokes.F90 gravity.F90 derivates.F90 main.F90 -o main
```

PARTÍCULA INDIVIDUAL: 
``` console
$ main 0 a e M w R
```
0 (o cualquier otro "int") es obligatorio

"a, e, M, w" son los elementos orbitales.

R es opcional.

PARALELO: (Recordar estar en algún entorno)
En este caso se debe tener un archivo con los datos de las partículas. La info está al inicio del archivo [parallel.py](./parallel.py#L3#L30))
``` console
$ python parallel.py
```
