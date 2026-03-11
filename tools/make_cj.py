import numpy as np
import coord_elem as cm

# Masas [x, y, z=0, masa]
masses = np.array(
    [
        [-0.114885114885115, 0.000000000000000e00, 0.0, 6.293706293706294e18],
        [114.885114885115, 0.000000000000000e00, 0.0, 6.293706293706294e15],
    ]
)

# Lambda rotacional
lam = 21.5423471626332

# Lambda keplerian
wkep = 45.42973133885822

# Radio mínimo
R = 115

# A cuánto aspiro?
CJ_target = 2.032 * (wkep * R) ** 2


# ---------------------------------------------------------------
# ---------------------------------------------------------------
# ---------------------------------------------------------------
# ---------------------------------------------------------------

# Masa total
msum = np.sum(masses[:, -1])


def write_elements_to_file(filename, particles, msum, coord=False):

    with open(filename, "w") as f:

        # Header

        for p in particles:

            x, y, z, vx, vy, vz = cm.rotating_to_inertial(*p, lam)

            # Build 6D inertial state vector
            xc = np.array([x, y, z, vx, vy, vz])

            if not coord:

                # Compute elements
                a, e, inc, capm, omega, capom = cm.elem(msum, xc)
                inc = np.rad2deg(inc)
                capm = np.rad2deg(capm)
                omega = np.rad2deg(omega)
                capom = np.rad2deg(capom)

                # Skip if pericenter too low
                if a * (1.0 - e) < R:
                    continue

                # Write to file
                f.write(
                    f"{a:.13e} {e:.13e} "
                    # f"{inc:.16e} "
                    f"{capm:.13e} {omega:.13e} "
                    # f"{capom:.16e}"
                    "\n"
                )

            else:

                # Write to file
                f.write(
                    f"{x:.13e} {y:.13e} "
                    # f"{inc:.16e} "
                    f"{vx:.13e} {vy:.13e} "
                    # f"{capom:.16e}"
                    "\n"
                )


if __name__ == "__main__":

    # Genero partículas

    ## Method 1
    particles = cm.generate_particles_fixed_CJ(
        N=100,
        CJ=CJ_target,
        masses=masses,
        lam=lam,
        r_range=(R * 4.5, R * 5),
        phi_range=(0, 2 * np.pi),
        z_range=(0, 0),
        clockwise=True,
        verbose=True,
    )

    ## Method 2
    # particles = cm.generate_grid_x_vx(
    #     CJ=CJ_target,
    #     masses=masses,
    #     lam=lam,
    #     x_range=(R * 2.7, R * 5.1),
    #     vx_range=(R * wkep * (-2), R * wkep * (2)),
    #     Nx=20,
    #     Nvx=20,
    #     y0=0.0,
    #     clockwise=True,
    #     verbose=True,
    # )

    # Write to file
    write_elements_to_file("particles_cj.dat", particles, msum, coord=False)

    # Report
    for particle in particles:
        x, y, z, vx, vy, vz = cm.rotating_to_inertial(*particle, lam)
        print(
            cm.jacobi_constant_inertial(
                x,
                y,
                z,
                vx,
                vy,
                vz,
                masses=masses,
                lam=lam,
            ),
            x,
            y,
            z,
            vx,
            vy,
            vz,
        )
