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
use_mu_and_radius = False  # True if a moon is expected

# UNIDADES
unit_angle = deg

# Data
data_in = {}
# Nombres ordenados de las columnas posibles
order = ["mu_to_disk", "a", "e", "M", "w", "mmr", "radius"]
# PARÁMETROS A VARIAR
## Only Moon parameters
disk_to_ast_mass_ratio = 0.1  # Unused if use_mu_and_radius is False.
data_in["mu_to_disk"] = [rndm(0.0, 1.0)]  # mass moon / mass disk (total moons)
data_in["radius"] = [0.0]  # [km]
## Particle parameters
data_in["a"] = [0.0]  # [km]
data_in["e"] = [n_steps(0.0, 0.2, 20)]  # Podría ser: rayleigh_dist(0, 0.1, 1)
data_in["M"] = [rndm(0.0, 360.0)]  # [unit_angle]
data_in["w"] = [rndm(0.0, 360.0)]  # [unit_angle]
data_in["mmr"] = [n_steps(1.45, 4.6, 50)]

# Nota:
#     'mu_to_disk' es solo para input; el output será 'mu_to_asteroid'

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
    """Cartesian product of keyword arguments, yielding dicts."""
    normalized = {
        k: v if isinstance(v, (list, set, dict, np.ndarray)) else [v]
        for k, v in kwargs.items()
    }
    keys, values = zip(*normalized.items())
    for instance in itertools.product(*values):
        yield dict(zip(keys, instance))


# Función para evaluar un generador y obtener los valores en una lista
@np.vectorize
def evaluate_generator(generator):
    if isinstance(generator, type_gen):
        return next(generator)
    return float(generator)


def check(data_in: dict):
    """Validate input dict values."""
    for values in data_in.values():
        assert len(values) == 1, "Each entry must contain exactly one element"
        val = values[0]
        if hasattr(val, "__len__"):
            assert len(val) > 0, "Inner value cannot be empty"


def make_ic(data_in: dict):
    """Expand input dict into combinations array."""
    check(data_in)  # Assuming this validates data_in

    # Get first element for each name (in correct order)
    first_values = {
        name: data_in[name][0] for name in order if name in data_in
    }

    # Base combinations
    base = [list(combo.values()) for combo in product_dict(**first_values)]

    return np.asarray(base).reshape(-1, len(order))


def check_continue(outfile: str):
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
        else:
            print("File is not overwritten.")
            return False
    return True


def shuffle_matrix(matrix: np.ndarray):
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
    # Remove mass if not needed
    if not use_mu_and_radius:
        if "mu_to_disk" in order:
            order.remove("mu_to_disk")
        if "radius" in order:
            order.remove("radius")

    # Use only keys in "order"
    sub_data_in = {key: data_in[key] for key in order if key in data_in}

    # Generamos las condiciones
    data_out = make_ic(sub_data_in)

    # Iterar sobre el arreglo y evaluar los generadores
    data_out = evaluate_generator(data_out)

    # Ponemos unidades
    target = ["inc", "M", "w", "O"]
    angles_indices = [idx for idx, elem in enumerate(order) if elem == target]
    data_out[:, angles_indices] /= unit_angle

    # Dropeamos si es necesario
    data_out = drop(data_out)

    # Mezclamos si se quiere
    if shuffle:
        data_out = shuffle_matrix(data_out)

    # Pasamos de kg a mass ratio si es necesario
    if use_mu_and_radius and (disk_to_ast_mass_ratio > 0):
        mass_idx = order.index("mu_to_disk")
        data_out[:, mass_idx] = data_out[:, mass_idx] * disk_to_ast_mass_ratio

    # Guardamos
    if not check_continue(filename):
        i = 1
        while os.path.isfile(filename):
            filename = filename + f"_{i}"
            i = i + 1
    np.savetxt(filename, data_out, fmt=fmt, delimiter=" ")

    print(f"File {filename} created.")
