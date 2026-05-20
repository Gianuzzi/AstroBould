import os

import numpy as np
import pandas as pd

chaos_names = [
    "nsim",
    "idx",
    "btype",
    "res",
    "time_f",
    "theta_i",
    "omega_i",
    "a_i",
    "e_i",
    "M_i",
    "w_i",
    "mmr_i",
    "mass_i",
    "radius_i",
    "theta_f",
    "omega_f",
    "a_f",
    "e_f",
    "M_f",
    "w_f",
    "mmr_f",
    "mass_f",
    "radius_f",
    "amin",
    "amax",
    "emin",
    "emax",
    "megno",
]

outnames_coord = {
    "idx": "i4",
    "btype": "i4",
    "t": "f8",
    "theta": "f8",
    "omega": "f8",
    "x": "f8",
    "y": "f8",
    "vx": "f8",
    "vy": "f8",
    "mass": "f8",
    "radius": "f8",
}

outnames_elem = {
    "idx": "i4",
    "btype": "i4",
    "t": "f8",
    "theta": "f8",
    "omega": "f8",
    "a": "f8",
    "e": "f8",
    "M": "f8",
    "w": "f8",
    "mmr": "f8",
    "mass": "f8",
    "radius": "f8",
    "dist": "f8",
    "amin": "f8",
    "amax": "f8",
    "emin": "f8",
    "emax": "f8",
    "megno": "f8",
}


def read_chaos(chaos_name="sump.out", work_dir=None, use_megno=True):
    if work_dir is None:
        full_name = chaos_name
    else:
        full_name = os.path.join(work_dir, chaos_name)
    df = pd.read_csv(full_name, delimiter="\\s+", header=None)
    if len(df.columns) == 28:
        df.columns = chaos_names
    elif len(df.columns == 27):
        if use_megno:
            df.columns = chaos_names[1:]
        else:
            df.columns = chaos_names[:-1]
    return df


def read_outfile(
    out_name="salida.out",
    work_dir=None,
    nboulders=0,
    coord=False,
    binary=False,
):
    if work_dir is None:
        full_name = out_name
    else:
        full_name = os.path.join(work_dir, out_name)
    if coord:
        if binary:
            dt = np.dtype([(k, v) for k, v in outnames_coord.items()])
            data = np.fromfile(full_name, dtype=dt)
            df = pd.DataFrame(data, columns=list(outnames_coord.keys()))
        else:
            df = pd.read_csv(
                full_name,
                delimiter="\\s+",
                header=None,
                names=list(outnames_coord.keys()),
            )
        df.nt = df.idx.nunique()
        df.nb = nboulders
        df.np = df.nt - nboulders
        df["r"] = np.sqrt(df["x"] ** 2 + df["y"] ** 2)
        df["v"] = np.sqrt(df["vx"] ** 2 + df["vy"] ** 2)
    else:
        if binary:
            dt = np.dtype([(k, v) for k, v in outnames_elem.items()])
            data = np.fromfile(full_name, dtype=dt)
            df = pd.DataFrame(data, columns=list(outnames_elem.keys()))
        else:
            df = pd.read_csv(
                full_name,
                delimiter="\\s+",
                header=None,
            )
        df.columns = list(outnames_elem.keys())[: len(df.columns)]
        df["lam"] = np.mod(df.M + df.w, 360.0)
        df.tmax = df.t.iloc[-1]
    return df


def read_map(map_name="map.out", work_dir=None):
    if work_dir is None:
        full_name = map_name
    else:
        full_name = os.path.join(work_dir, map_name)
    df = pd.read_csv(
        full_name,
        delimiter="\\s+",
        header=None,
        names=["x", "y", "pot", "ax", "ay"],
    )
    df["r"] = np.sqrt(df["x"] ** 2 + df["y"] ** 2)
    df["acc"] = np.sqrt(df["ax"] ** 2 + df["ay"] ** 2)
    return df


def read_surface(
    surface_name="surface.out", work_dir=None, omega_rotation=0.0
):
    if work_dir is None:
        full_name = surface_name
    else:
        full_name = os.path.join(work_dir, surface_name)
    df = pd.read_csv(
        full_name,
        header=None,
        delimiter="\\s+",
        names=["idx", "t", "theta", "omega", "x", "y", "vx", "vy", "CJ"],
    )
    df["r"] = np.sqrt(df["x"] ** 2 + df["y"] ** 2)
    df["vr"] = (df["x"] * df["vx"] + df["y"] * df["vy"]) / df["r"]
    df["psi"] = np.mod(np.arctan2(df["y"], df["x"]), 2 * np.pi)
    if omega_rotation != 0.0:
        df_inertial = __rotating_to_inertial_df(df, omega_rotation)
    else:
        print(
            "WARNING: omega_rotation = 0.0, no se hará la transformación a inercial."
        )
        df_inertial = df.copy()
    return df_inertial


def surface_to_elements(
    surface_df, omega_rotation=0.0, msum=0.0, input_inertial=False
):
    from coord_elem import elem

    if not input_inertial and omega_rotation != 0.0:
        surface_df = __rotating_to_inertial_df(surface_df, omega_rotation)
    if msum == 0.0:
        print(
            "ERROR: msum = 0.0, no se pueden calcular los elementos orbitales."
        )
        return None
    elements = pd.DataFrame(
        np.array(
            [
                elem(msum, (r.x, r.y, 0.0, r.vx, r.vy, 0.0))
                for r in surface_df.itertuples(index=False)
            ]
        ),
        columns=["a", "e", "i", "M", "omega", "Omega"],
    )
    return elements


def __rotating_to_inertial_df(df, omega):
    df = df.copy()
    df["vx"] = df["vx"] - omega * df["y"]
    df["vy"] = df["vy"] + omega * df["x"]
    return df
