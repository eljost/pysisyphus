from collections import namedtuple
import logging

import numpy as np

from pysisyphus.constants import BOHR2M, AMU2KG, AU2J, AU2MPERSEC, FORCE2ACC
from pysisyphus.dynamics.helpers import kinetic_energy_from_velocities, \
                                        temperature_for_kinetic_energy


MDResult = namedtuple("MDResult",
                      "coords t",
)


def md(geom, v0, t, dt, term_funcs=None, rm_vcom=False):
    """TODO: dump coords, velocities; check energy conservation.
    Align geometry to avoid drift."""
    """Velocity verlet integrator.

    Parameters
    ----------
    geom : Geometry
        The system to be integrated.
    v0 : np.array, floats
        Initial velocities.
    t : float
        Total simulation time i n fs.
    dt : float
        Timestep in fs.
    term_funcs : optional
        Add description.
    rm_vcom : bool, optional, default=False
        Remove center of mass velocity to avoid drift.
    """

    assert geom.coord_type == "cart"

    steps = int(t/dt)
    print(f"Doing {steps} steps of {dt:.1f} fs for a total of {steps*dt:.1f} fs.")

    # Convert to seconds
    t = t  * 1e-15
    dt = dt * 1e-15

    if term_funcs is None:
        term_funcs = list()

    m = geom.masses_rep * AMU2KG
    M = geom.masses.sum() * AMU2KG
    # x = geom.cart_coords * BOHR2M
    x = geom.coords * BOHR2M
    v = v0 * AU2MPERSEC
    a_prev = 0
    xs = list()

    t_cur = 0
    while t_cur < t:
        print(t_cur)
        t_cur += dt
        xs.append(x.copy())
        try:
            f = geom.forces * AU2J / BOHR2M
        except:
            logging.exception("Force calculation in MD failed.")
            break
        a = f / m
        v += .5 * (a + a_prev) * dt
        # Remove center of mass velocity
        if rm_vcom:
            v -= np.sum(m*v / M)
        x += v*dt + .5*a*dt**2
        geom.coords = x / BOHR2M
        a_prev = a
        if any([tf(x) for tf in term_funcs]):
            break
    md_result = MDResult(
                    coords=np.array(xs)/BOHR2M,
                    t=t,
    )

    return md_result


def md(geom, v0, steps, dt):
    """Velocity verlet integrator.

    Parameters
    ----------
    geom : Geometry
        The system to be integrated.
    v0 : np.array, floats
        Initial velocities in Bohr/fs.
    steps : float
        Number of simulation steps.
    dt : float
        Timestep in fs.
    """

    assert geom.coord_type == "cart"

    # Total simulation time
    t_ps = steps * dt * 1e-3
    print(f"Doing {steps} steps of {dt:.1f} fs for a total of {t_ps:.2f} ps.")

    # In amu
    masses = geom.masses
    masses_rep = geom.masses_rep

    x = geom.cart_coords
    # v is given in Bohr/fs
    v = v0
    a_prev = np.zeros_like(x)
    xs = list()

    t_cur = 0
    for i in range(steps):
        xs.append(x.copy())
        t_cur += dt
        try:
            forces = geom.forces
            E_pot = geom.energy
        except:
            logging.exception("Force calculation in MD failed.")
            break

        if (i % 25) == 0:
            E_kin = kinetic_energy_from_velocities(masses, v.reshape(-1, 3))
            T = temperature_for_kinetic_energy(len(masses), E_kin)
            E_tot = E_pot + E_kin
            print(f"Step {i:05d}  {t_cur*1e-3: >6.2f} ps  E={E_tot: >8.6f} E_h  "
                  f"T={T: >8.2f} K"
            )

        # Acceleration, convert from Hartree / (Bohr * amu) to Bohr/fs²
        a = forces / masses_rep * FORCE2ACC
        v += .5 * (a + a_prev) * dt
        # v*dt = Bohr/fs * fs -> Bohr
        # a*dt**2 = Bohr/fs² * fs² -> Bohr
        x += v*dt + .5*a*dt**2
        geom.coords = x
        a_prev = a

    md_result = MDResult(
                    coords=np.array(xs),
                    t=t_ps,
    )

    return md_result
