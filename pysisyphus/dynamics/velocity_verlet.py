from collections import namedtuple
import logging

import numpy as np

from pysisyphus.constants import FORCE2ACC
from pysisyphus.dynamics.helpers import kinetic_energy_from_velocities, \
                                        temperature_for_kinetic_energy

logger = logging.getLogger("dynamics")

MDResult = namedtuple("MDResult",
                      "coords t",
)


def md(geom, v0, steps, dt, term_funcs=None, verbose=True):
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
    term_funcs : dict, optional
        Iterable of functions that are called with the atomic
        coordinates in every MD cycle and result in termination
        of the MD integration when they evaluate to true.
    """

    assert geom.coord_type == "cart"

    if term_funcs is None:
        term_funcs = dict()

    if verbose:
        t_ps = steps * dt * 1e-3  # Total simulation time
        print(f"Doing {steps} steps of {dt:.4f} fs for a total of {t_ps:.2f} ps.")

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

        if verbose and (i % 25) == 0:
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

        terminate = False
        for name, func in term_funcs.items():
            if func(x):
                logger.debug(f"Termination function '{name}' evaluted to True. Breaking.")
                terminate = True
                break
        if terminate:
            break

    md_result = MDResult(
                    coords=np.array(xs),
                    t=t_cur*1e-3,
    )

    return md_result
