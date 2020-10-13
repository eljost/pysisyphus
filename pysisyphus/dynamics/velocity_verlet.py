from collections import namedtuple
import logging

import numpy as np

from pysisyphus.constants import FORCE2ACC
from pysisyphus.dynamics.helpers import kinetic_energy_from_velocities, \
                                        temperature_for_kinetic_energy, \
                                        energy_forces_getter_closure, \
                                        kinetic_energy_for_temperature
from pysisyphus.dynamics.csvr import resample_kin
from pysisyphus.dynamics.rattle import rattle_closure
from pysisyphus.helpers import check_for_stop_sign

logger = logging.getLogger("dynamics")

MDResult = namedtuple("MDResult",
                      "coords t terminated T E_tot",
)
THERMOSTATS = {
    "csvr": resample_kin,
}


def md(geom, v0, steps, dt, remove_com_v=True, thermostat=None, T=298.15,
       timecon=100, term_funcs=None, constraints=None, constraint_kwargs=None,
       verbose=True):
    """Velocity verlet integrator.

    Parameters
    ----------
    geom : Geometry
        The system for which the dynamics are to be run.
    v0 : np.array, floats
        Initial velocities in Bohr/fs.
    steps : float
        Number of simulation steps.
    dt : float
        Timestep in fs.
    remove_com_v : bool, default=True
        Remove center-of-mass velocity.
    thermostat : str, optional, default None
        Which and whether to use a thermostat.
    T : float, optional, default=None
        Desired temperature in thermostated runs.
    timecon : float
        Timeconsanst of the thermostat in fs.
    term_funcs : dict, optional
        Iterable of functions that are called with the atomic
        coordinates in every MD cycle and result in termination
    constraints : 2d iterable, optional
        2D iterable containing atom indices describing constrained
        bond lengths.
        of the MD integration when they evaluate to true.
    constraint_kwargs : dict, optional
        Keyword arguments for the constraint algorithm.
    verbose : bool, default=True
        Do additional printing when True.
    """

    assert geom.coord_type == "cart"

    if term_funcs is None:
        term_funcs = dict()

    if verbose:
        t_ps = steps * dt * 1e-3  # Total simulation time
        print(f"Doing {steps} steps of {dt:.4f} fs for a total of {t_ps:.2f} ps.")

    energy_forces_getter = energy_forces_getter_closure(geom)

    if constraint_kwargs is None:
        constraint_kwargs = dict()

    # Fixed degrees of freedom
    fixed_dof = 0

    if remove_com_v:
        fixed_dof += 3

    constrained_md = constraints is not None
    # Get RATTLE function from closure for constrained MD
    if constrained_md:
        fixed_dof += len(constraints)
        rattle = rattle_closure(geom, constraints, dt,
                                energy_forces_getter=energy_forces_getter,
                                **constraint_kwargs)

    if thermostat is not None:
        thermo_func = THERMOSTATS[thermostat]
        tau_t = dt / timecon
        sigma = kinetic_energy_for_temperature(len(geom.atoms), T,
                                               fixed_dof=fixed_dof)

    # In amu
    masses = geom.masses
    masses_rep = geom.masses_rep
    total_mass = masses.sum()

    x = geom.cart_coords
    # v is given in Bohr/fs
    v = v0
    a_prev = np.zeros_like(x)
    xs = list()
    Ts = list()
    E_tots = list()

    E_pot, forces = energy_forces_getter(geom.coords)

    t_cur = 0
    terminate = False
    terminate_key = None
    T_avg = 0
    for i in range(steps):
        xs.append(x.copy())

        E_kin = kinetic_energy_from_velocities(masses, v.reshape(-1, 3))
        T = temperature_for_kinetic_energy(len(masses), E_kin, fixed_dof=fixed_dof)
        T_avg += T
        Ts.append(T)
        E_tot = E_pot + E_kin
        E_tots.append(E_tot)

        if verbose and (i % 25) == 0:
            print(f"Step {i:05d}  {t_cur*1e-3: >6.2f} ps  E={E_tot: >8.6f} E_h  "
                  f"T={T: >8.2f} K <T>={T_avg/(i+1): >8.2f}"
            )

        if thermostat:
            E_kin_new = thermo_func(E_kin, sigma, v.size-fixed_dof, tau_t)
            scale = (E_kin_new / E_kin)**0.5
            v *= scale

        # RATTLE algorithm
        if constrained_md:
            x, v, E_pot, forces = rattle(x, v, forces)
        # Simple Velocity-Verlet integration
        else:
            E_pot, forces = energy_forces_getter(geom.coords)
            # Acceleration, convert from Hartree / (Bohr * amu) to Bohr/fs²
            a = forces / masses_rep * FORCE2ACC
            v += .5 * (a + a_prev) * dt
            if remove_com_v:
                v -= v * masses_rep / total_mass
            # v*dt = Bohr/fs * fs -> Bohr
            # a*dt**2 = Bohr/fs² * fs² -> Bohr
            x += v*dt + .5*a*dt**2
            a_prev = a

        # Update coordinates
        geom.coords = x

        for name, func in term_funcs.items():
            if func(x):
                terminate = True
                terminate_key = name
                break
        if terminate:
            logger.debug(f"Termination function '{name}' evaluted to True. Breaking.")
            break

        if check_for_stop_sign():
            break

        # Advance time
        t_cur += dt

    md_result = MDResult(
                    coords=np.array(xs),
                    t=t_cur*1e-3,
                    terminated=terminate_key,
                    T=np.array(Ts),
                    E_tot=np.array(E_tots),
    )

    return md_result
