import numpy as np

# ---- Constants ----
G = 4.9823394e-10  # km^3 / day^2 / kg
pi = np.pi
twopi = 2.0 * np.pi
myepsilon = np.finfo(float).eps
sqepsilon = np.sqrt(myepsilon)


def aver(dm, e, max_iter=1000):
    """
    Solve Kepler equation:
        u = dm + e sin(u)
    Returns:
        u  : eccentric anomaly (ellipse)
        f  : true anomaly
    """
    u0 = dm
    dif = 1.0
    i = 0
    u = 0.0

    while dif > myepsilon:
        u = dm + e * np.sin(u0)
        dif = abs(u - u0)
        u0 = u
        i += 1
        if i > max_iter:
            print("aver: too many iterations. Error:", dif)
            break

    seno = np.sqrt(1.0 + e) * np.sin(0.5 * u)
    cose = np.sqrt(1.0 - e) * np.cos(0.5 * u)
    f = 2.0 * np.arctan2(seno, cose)

    return u, f


def coord(msum, a, e, inc, capm, omega, capom):
    """
    Convert orbital elements to Cartesian state vector.
    Returns:
        xc : array(6) -> [x,y,z,vx,vy,vz]
    """

    # Rotation matrix components
    sp = np.sin(omega)
    cp = np.cos(omega)
    so = np.sin(capom)
    co = np.cos(capom)
    si = np.sin(inc)
    ci = np.cos(inc)

    d11 = cp * co - sp * so * ci
    d12 = cp * so + sp * co * ci
    d13 = sp * si

    d21 = -sp * co - cp * so * ci
    d22 = -sp * so + cp * co * ci
    d23 = cp * si

    # Solve Kepler
    cape, _ = aver(capm, e)

    scap = np.sin(cape)
    ccap = np.cos(cape)

    sqe = np.sqrt(1.0 - e * e)
    sqgma = np.sqrt(G * msum * a)

    xfac1 = a * (ccap - e)
    xfac2 = a * sqe * scap

    ri = 1.0 / (a * (1.0 - e * ccap))

    vfac1 = -ri * sqgma * scap
    vfac2 = ri * sqgma * sqe * ccap

    xc = np.zeros(6)

    xc[0] = d11 * xfac1 + d21 * xfac2
    xc[1] = d12 * xfac1 + d22 * xfac2
    xc[2] = d13 * xfac1 + d23 * xfac2

    xc[3] = d11 * vfac1 + d21 * vfac2
    xc[4] = d12 * vfac1 + d22 * vfac2
    xc[5] = d13 * vfac1 + d23 * vfac2

    return xc


def elem(msum, xc):
    """
    Convert Cartesian state vector to orbital elements.
    Returns:
        a, e, inc, capm, omega, capom
    """

    gmsum = G * msum

    x, y, z, vx, vy, vz = xc

    # Angular momentum
    hx = y * vz - z * vy
    hy = z * vx - x * vz
    hz = x * vy - y * vx

    h2 = hx * hx + hy * hy + hz * hz
    h = np.sqrt(h2)

    inc = np.arccos(hz / h)

    fac = np.sqrt(hx * hx + hy * hy) / h

    if fac < myepsilon:
        capom = 0.0
        u = np.arctan2(y, x)
        if abs(inc - pi) < myepsilon:
            u = -u
    else:
        capom = np.arctan2(hx, -hy)
        u = np.arctan2(z / np.sin(inc), x * np.cos(capom) + y * np.sin(capom))

    if capom < 0:
        capom += twopi
    if u < 0:
        u += twopi

    r = np.sqrt(x * x + y * y + z * z)
    v2 = vx * vx + vy * vy + vz * vz
    v = np.sqrt(v2)
    vdotr = x * vx + y * vy + z * vz

    energy = 0.5 * v2 - gmsum / r

    # Determine orbit type
    if abs(energy * r / gmsum) < myepsilon:
        ialpha = 0
    elif energy < 0:
        ialpha = -1
    else:
        ialpha = 1

    # =========================
    # ELLIPSE
    # =========================
    if ialpha == -1:
        a = -0.5 * gmsum / energy
        fac = 1.0 - h2 / (gmsum * a)

        if fac > myepsilon:
            e = np.sqrt(fac)
            face = (a - r) / (a * e)

            face = np.clip(face, -1.0, 1.0)
            cape = np.arccos(face)

            if vdotr < 0:
                cape = twopi - cape

            cw = (np.cos(cape) - e) / (1 - e * np.cos(cape))
            sw = np.sqrt(1 - e * e) * np.sin(cape) / (1 - e * np.cos(cape))

            w = np.arctan2(sw, cw)
            if w < 0:
                w += twopi
        else:
            e = 0.0
            w = u
            cape = u

        capm = cape - e * np.sin(cape)
        omega = u - w

    # =========================
    # HYPERBOLA
    # =========================
    elif ialpha == 1:
        a = 0.5 * gmsum / energy
        fac = h2 / (gmsum * a)

        if fac > myepsilon:
            e = np.sqrt(1.0 + fac)
            tmpf = (a + r) / (a * e)
            tmpf = max(tmpf, 1.0)

            capf = np.log(tmpf + np.sqrt(tmpf * tmpf - 1.0))

            if vdotr < 0:
                capf = -capf

            cw = (e - np.cosh(capf)) / (e * np.cosh(capf) - 1.0)
            sw = (
                np.sqrt(e * e - 1.0)
                * np.sinh(capf)
                / (e * np.cosh(capf) - 1.0)
            )

            w = np.arctan2(sw, cw)
            if w < 0:
                w += twopi
        else:
            e = 1.0
            w = np.arccos(2 * h2 / gmsum / r - 1.0)
            if vdotr < 0:
                w = twopi - w

            tmpf = (a + r) / (a * e)
            capf = np.log(tmpf + np.sqrt(tmpf * tmpf - 1.0))

        capm = e * np.sinh(capf) - capf
        omega = u - w

    # =========================
    # PARABOLA
    # =========================
    else:
        a = 0.5 * h2 / gmsum
        e = 1.0
        w = np.arccos(2 * a / r - 1.0)

        if vdotr < 0:
            w = twopi - w

        tmpf = np.tan(0.5 * w)
        capm = tmpf * (1 + tmpf * tmpf / 3.0)
        omega = u - w

    # Normalize omega
    if omega < 0:
        omega += twopi
    omega = omega % twopi

    return a, e, inc, capm, omega, capom


# --------------------------------------------------
# Inertial <--> Rotational
# --------------------------------------------------


def inertial_to_rotating(x, y, z, vx, vy, vz, lam):
    """
    Convert inertial state vector to rotating frame.

    lam : angular velocity (rad/day)
    """
    # --- Rotate velocity ---
    vx_rot = vx + lam * y
    vy_rot = vy - lam * x

    return x, y, z, vx_rot, vy_rot, vz


def rotating_to_inertial(x, y, z, vx, vy, vz, lam):
    """
    Convert rotating frame state vector to inertial frame.

    lam : angular velocity (rad/day)
    """

    # --- Add frame motion ---
    vx_in = vx - lam * y
    vy_in = vy + lam * x

    return x, y, z, vx_in, vy_in, vz


# --------------------------------------------------
# Gravitational Potential
# --------------------------------------------------


def gravitational_potential(x, y, z, masses):
    """
    masses: array-like of shape (N, 3)
            columns = [x_i, y_i, m_i]
    """
    U = 0.0
    for xi, yi, zi, mi in masses:
        r = np.sqrt((x - xi) ** 2 + (y - yi) ** 2 + (z - zi) ** 2)
        U += G * mi / r
    return U


# --------------------------------------------------
# Jacobi Constant
# --------------------------------------------------


def jacobi_constant(x, y, z, vx, vy, vz, masses, lam):
    U = gravitational_potential(x, y, z, masses)
    return lam**2 * (x**2 + y**2) + 2 * U - (vx**2 + vy**2 + vz**2)


def jacobi_constant_inertial(x, y, z, vx, vy, vz, masses, lam):
    U = gravitational_potential(x, y, z, masses)
    x, y, z, vx_rot, vy_rot, vz = inertial_to_rotating(
        x, y, z, vx, vy, vz, lam
    )
    return lam**2 * (x**2 + y**2) + 2 * U - (vx_rot**2 + vy_rot**2 + vz**2)


# --------------------------------------------------
# Generate particles with fixed CJ
# --------------------------------------------------


def generate_particles_fixed_CJ(
    N,
    CJ,
    masses,
    lam,
    r_range,
    phi_range,
    z_range,
    max_iter=100000,
    clockwise=False,
    verbose=False,
):

    particles = []
    count = 0

    while len(particles) < N and count < max_iter:
        count += 1

        # Sample position
        r = np.random.uniform(*r_range)
        phi_pos = np.random.uniform(*phi_range)
        z = np.random.uniform(*z_range)

        x = r * np.cos(phi_pos)
        y = r * np.sin(phi_pos)
        # Compute gravitational potential
        U = gravitational_potential(x, y, z, masses)

        # Allowed velocity magnitude squared
        V2 = 2 * U + lam**2 * r**2 - CJ

        if V2 >= 0:
            V = np.sqrt(V2)

            # --- Random 3D velocity direction (uniform on sphere) ---
            phi_vel = np.random.uniform(0, 2 * np.pi)
            cos_theta = 0  # np.random.uniform(-1, 1)
            sin_theta = np.sqrt(1 - cos_theta**2)

            vx = V * sin_theta * np.cos(phi_vel)
            vy = V * sin_theta * np.sin(phi_vel)
            vz = V * cos_theta

            # This ensures that the velocity is in the plane (z=0) and that it is tangential (vr=0)
            vr = 0
            vt = np.sqrt(V**2 - vr**2)

            vx = vr * np.cos(phi_pos) - vt * np.sin(phi_pos)
            vy = vr * np.sin(phi_pos) + vt * np.cos(phi_pos)
            vz = 0

            # Force clockwise (so it is direct motion in inertial frame)
            if (x * vy - y * vx > 0) and clockwise:
                vx = -vx
                vy = -vy

            particles.append([x, y, z, vx, vy, vz])

    # Report
    if verbose:
        print(
            f"Generated {len(particles)} particles with CJ={CJ:.3f} in {count} attempts."
        )

    return np.array(particles)


def generate_grid_x_vx(
    CJ,
    masses,
    lam,
    x_range,
    vx_range,
    Nx,
    Nvx,
    y0=0.0,
    clockwise=False,
    verbose=False,
):
    """
    Generate regular grid in (x, vx) plane (z=0, vz=0, y=y0)

    direction: "clockwise" or "anticlockwise"
    """

    particles = []

    x_vals = np.linspace(x_range[0], x_range[1], Nx)
    vx_vals = np.linspace(vx_range[0], vx_range[1], Nvx)

    z = 0.0
    vz = 0.0
    y = y0
    for x in x_vals:
        for vx in vx_vals:

            r2 = x * x + y * y
            U = gravitational_potential(x, y, z, masses)

            V2 = 2 * U + lam**2 * r2 - CJ

            vy2 = V2 - vx**2

            if vy2 >= 0:

                vy = np.sqrt(vy2)

                # --- enforce direction ---
                cross = x * vy - y * vx

                if clockwise and (cross > 0):
                    vy = -vy
                elif not clockwise and (cross < 0):
                    vy = -vy

                particles.append([x, y, z, vx, vy, vz])

    # Report
    if verbose:
        print(
            f"Generated {len(particles)} particles on grid with CJ={CJ:.3f}."
        )

    return np.array(particles)
