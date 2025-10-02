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
]


def read_chaos(chaos_name="sump.out", work_dir=None):
    if work_dir is None:
        full_name = chaos_name
    else:
        full_name = os.path.join(work_dir, chaos_name)
    df = pd.read_csv(full_name, delimiter="\\s+", header=None)
    df.columns = chaos_names if len(df.columns == 28) else chaos_names[1:]
    return df


def read_outfile(
    out_name="salida.out", work_dir=None, nboulders=0, coord=False
):
    if work_dir is None:
        full_name = out_name
    else:
        full_name = os.path.join(work_dir, out_name)
    if coord:
        df = pd.read_csv(
            full_name,
            delimiter="\\s+",
            header=None,
            names=[
                "idx",
                "btype",
                "t",
                "theta",
                "omega",
                "x",
                "y",
                "vx",
                "vy",
                "mass",
                "radius",
            ],
        )
        df.nt = df.idx.nunique()
        df.nb = nboulders
        df.np = df.nt - nboulders
        df["r"] = np.sqrt(df["x"] ** 2 + df["y"] ** 2)
        df["v"] = np.sqrt(df["vx"] ** 2 + df["vy"] ** 2)
    else:
        df = pd.read_csv(
            full_name,
            delimiter="\\s+",
            header=None,
            names=[
                "idx",
                "btype",
                "t",
                "theta",
                "omega",
                "a",
                "e",
                "M",
                "w",
                "mmr",
                "mass",
                "radius",
                "dist",
            ],
        )
        df["lam"] = np.mod(df.M + df.w, 360.)
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
