import os
import numpy as np
import itertools


# ---------------------------------------------------------------------
# Constantes y conversiones -> Ms | AU | yr | deg
pi = np.pi
rad = pi / 180.0
deg = 1.0
G = 4 * pi * pi  # AU³ Ms⁻¹ yr⁻¹
Msol = 1.0
Mjup = 9.547919e-4
Mear = 3.00273e-6
Parsec = 4.84814e-6
AU = 1.0
Rsol = 4.65047e-3
Rjup = 4.77895e-4
Rear = 4.26352e-5
Gyr = 1e9
Myr = 1e6
yr = 1.0
Day = 2.73785078e-3
Sec = 3.17098e-8

rng = np.random.default_rng(seed=42)
# ---------------------------------------------------------------------


# ---------------------------------------------------------------------
# Funciones útiles
def n_steps(x0, xf, npts):
    """Comienzo, final (incluido), cantidad de puntos"""
    return np.linspace(x0, xf, npts)


def size_steps(x0, xf, stepsize):
    """Comienzo, final (incluido de ser posible), tamaño del paso"""
    return np.arange(x0, xf + stepsize, stepsize)


def concat(*arrs):
    """Concatena arrays o listas"""
    return np.concatenate([*arrs])


def rndm(x0=0.0, xf=360.0):
    """Número aleatorio entre x0 y xf"""
    while True:
        yield rng.random() * (xf - x0) + x0


type_gen = type(rndm())

# ---------------------------------------------------------------------

# -----------------
# INPUT
# -----------------

## UNIDADES
unit_dist = AU
unit_angle = deg

## Data
data_in = {}
names = ["a", "e", "M", "w", "R"]
## PARÁMETROS A VARIAR (Poner unidades de ser necesario)
data_in["a"] = [0.0]
data_in["e"] = [0.0]
data_in["M"] = [rndm(0.0, 360.0)]
data_in["w"] = [0.0]
data_in["R"] = [n_steps(0.9, 7.5, 3000)]

# -----------------
# OUTPUT
# -----------------

## ARCHIVO DE SALIDA
filename = "particles.in"

# OUTPUT FORMAT:
fmt = "%.11e"

## Shuffle
shuffle = False


# --------Función para droper. Editar rem de ser necesario-------------


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
        print("WARNING: File {} already exist.".format(outfile))
        print("Do yo want to overwrite it? Y/[N]")
        q = input()
        ntry = 3
        while q.lower() not in ["y", "yes", "s", "si", "n", "no"]:
            print("{} is not a valid answer ({} attempts left)".format(q, ntry))
            print("Do yo want to overwrite it? Y/[N]")
            q = input()
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


def read_IC(path):
    IC = pd.read_csv(
        path, delim_whitespace=True, header=None, names=["a", "e", "M", "w", "R"]
    )
    return IC


# MAIN PROGRAM

if __name__ == "__main__":
    # Generamos las condiciones
    data_out = make_ic(data_in)

    # Iterar sobre el arreglo y evaluar los generadores
    data_out = evaluate_generator(data_out)

    # Ponemos unidades
    data_out[:, 0] /= unit_dist
    data_out[:, [2, 3]] /= unit_angle

    # Dropeamos si es necesario
    data_out = drop(data_out)

    # Mezclamos si se quiere
    if shuffle:
        data_out = shuffle_matrix(data_out)

    # Guardamos
    if not check_continue(filename):
        i = 1
        while os.path.isfile(filename):
            filename = filename + "_{}".format(i)
            i = i + 1
    np.savetxt(filename, data_out, fmt=fmt, delimiter=" ")
    print("File {} created.".format(filename))
