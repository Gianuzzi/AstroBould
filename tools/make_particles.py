import itertools
import os

import numpy as np
from scipy.stats import norm, rayleigh

# ---------------------------------------------------------------------
# Constantes
pi = np.pi
rad = pi / 180.0
deg = 1.0

# Semilla para el generador de números aleatorios
rng = np.random.default_rng(seed=42)
# ---------------------------------------------------------------------


# ---------------------------------------------------------------------
# Funciones útiles
def n_steps(x0, xf, npts, log=False, log_base=True):
    """Comienzo, final (incluído), cantidad de puntos"""
    if log:
        if log_base:
            return np.logspace(x0, xf, npts)
        return np.logspace(np.log10(x0), np.log10(xf), npts)
    return np.linspace(x0, xf, npts)


def size_steps(x0, xf, stepsize):
    """Comienzo, final (incluído de ser posible), tamaño del paso"""
    return np.arange(x0, xf + stepsize, stepsize)


def concat(*arrs):
    """Concatena arrays o listas"""
    return np.concatenate([*arrs])


def rndm(x0=0.0, xf=360.0, func=None):
    """Número aleatorio entre x0 y xf"""
    while True:
        if callable(func):
            yield func(rng.random() * (xf - x0) + x0)
        else:
            yield rng.random() * (xf - x0) + x0


def rayleigh_dist(x0=0, var=0.1, xmax=1.0):
    """Rayleigh distribution
    x0: Valor mínimo
    var: Varianza (define el pico de la distribución)
    xmax: Valor máximo
    """
    while True:
        x = rayleigh.rvs(loc=x0, scale=var)
        if x <= xmax:
            yield x


def norm_dist(x0=0.3, scale=0.01, xmin=0, xmax=1.0):
    """Normal distribution
    x0: Media
    scale: Desviación estándar
    xmax: Valor máximo
    """
    while True:
        x = norm.rvs(loc=x0, scale=scale)
        if xmin <= x <= xmax:
            yield x


type_gen = type(rndm())

# ---------------------------------------------------------------------

# -----------------
# INPUT
# -----------------
use_mass = True

# UNIDADES
unit_angle = deg

# Data
data_in = {}
names = ["mass", "a", "e", "M", "w", "radius"]
# PARÁMETROS A VARIAR (Unidades según código original)
data_in["mass"] = [rndm(0.0, 1.0)]
data_in["a"] = [0.0]
data_in["e"] = [0.0]  # Podría ser: rayleigh_dist(0, 0.1, 1)
data_in["M"] = [rndm(0.0, 360.0)]
data_in["w"] = [rndm(0.0, 360.0)]
data_in["mmr"] = [n_steps(1.4, 4, 5)]
data_in["radius"] = [0.0]

# Disk mass (in kg)
disk_mass = 6.3e15

# -----------------
# OUTPUT
# -----------------

# ARCHIVO DE SALIDA
filename = "particles.in"

# OUTPUT FORMAT:
fmt = "%.11e"

# Shuffle
shuffle = False

# -- No tocar después de esta línea --

# --------Función para dropear. Editar rem de ser necesario-------------


def drop(data):
    rem = np.full(len(data), False)  # = (data[:,5] == 90) & (data[:,6] == 90)
    drop = np.full(len(data), False)
    where_rem = np.where(rem)[0]
    rem = np.unique(np.hstack([where_rem]))
    rem = rem[rem >= 0]
    drop[rem] = True
    return data[~drop]


# ---------------------------------------------------------------------
# Funciones auxiliares
# ---------------------------------------------------------------------


def rad2deg(x):
    return np.mod(x / rad, 360)


def deg2rad(x):
    return np.mod(x * rad, 2 * pi)


def product_dict(**kwargs):
    keys = kwargs.keys()
    for key in keys:
        val = kwargs[key]
        if not isinstance(val, (list, set, dict, np.ndarray)):
            kwargs[key] = [val]
    vals = kwargs.values()
    for instance in itertools.product(*vals):
        yield dict(zip(keys, instance))


# Función para evaluar un generador y obtener los valores en una lista
@np.vectorize
def evaluate_generator(generator):
    if isinstance(generator, type_gen):
        return next(generator)
    return float(generator)


def check(data_in):
    for d in data_in.values():
        assert len(d) == 1
        for i, val in enumerate(d):
            if hasattr(val, "__len__"):
                assert len(val) > 0


def make_ic(data_in):
    check(data_in)
    data_aux = {}
    data_full = {}
    data_aux[0] = [d[0] for d in data_in.values()]
    spl = str(0)
    data_full[spl] = list(
        product_dict(**{names[i]: data_aux[0][i] for i in range(len(names))})
    )
    data_full[spl] = [list(d.values()) for d in data_full[spl]]
    arr = [list(d.values()) for d in product_dict(**data_full)]
    return np.asarray(arr).reshape(-1, len(names))


def check_continue(outfile):
    if os.path.isfile(outfile):
        print(f"WARNING: File {outfile} already exist.")
        q = input("Do yo want to overwrite it? [y/[n]]: ")
        ntry = 3
        while q.lower() not in ["y", "yes", "s", "si", "n", "no"]:
            print(f"{q} is not a valid answer.")
            print(f" ({ntry:1d} attempts left)")
            q = input("Do yo want to overwrite it? [y/[n]]: ")
            ntry -= 1
            if ntry == 0:
                raise Exception("Wrong input.")
        if q.lower() in ["y", "yes", "s", "si"]:
            os.remove(outfile)
            return True
        print("File is not overwritten.")
        return False
    return True


def shuffle_matrix(matrix):
    num_rows, num_cols = matrix.shape
    grouped_matrix = matrix.reshape(-1, 1, num_cols)  # group the rows together
    np.random.shuffle(grouped_matrix)  # shuffle the rows within each group
    shuffled_matrix = grouped_matrix.reshape(
        num_rows, num_cols
    )  # flatten the grouped matrix back into a single matrix
    return shuffled_matrix


# ---------------------------------------------------------------------
# ---------------------------------------------------------------------


# MAIN PROGRAM

if __name__ == "__main__":
    # Generamos las condiciones
    data_out = make_ic(data_in)

    # Iterar sobre el arreglo y evaluar los generadores
    data_out = evaluate_generator(data_out)

    # Ponemos unidades
    data_out[:, [2, 3]] /= unit_angle

    # Dropeamos si es necesario
    data_out = drop(data_out)

    # Mezclamos si se quiere
    if shuffle:
        data_out = shuffle_matrix(data_out)

    # Remove mass if not needed
    if not use_mass:
        data_out = data_out[:, 1:]
    elif disk_mass > 0:
        temp = np.sum(data_out[:, 0])
        if temp > 0:
            data_out[:, 0] = disk_mass * data_out[:, 0] / temp

    # Guardamos
    if not check_continue(filename):
        i = 1
        while os.path.isfile(filename):
            filename = filename + f"_{i}"
            i = i + 1
    np.savetxt(filename, data_out, fmt=fmt, delimiter=" ")
    print(f"File {filename} created.")
