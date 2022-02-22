from collections import namedtuple
import logging

import numpy as np

from pysisyphus.config import T_DEFAULT
from pysisyphus.constants import FORCE2ACC
from pysisyphus.dynamics.helpers import (
    kinetic_energy_from_velocities,
    temperature_for_kinetic_energy,
    energy_forces_getter_closure,
    kinetic_energy_for_temperature,
    remove_com_velocity,
)
from pysisyphus.dynamics.thermostats import csvr_closure, csvr_closure_2, berendsen_closure
from pysisyphus.dynamics.rattle import rattle_closure
from pysisyphus.helpers import check_for_end_sign
from pysisyphus.helpers_pure import log
from pysisyphus.io.hdf5 import get_h5_group

logger = logging.getLogger("dynamics")

MDResult = namedtuple(
    "MDResult",
    "coords t_ps step terminated T E_tot",
)

THERMOSTATS = {
    "csvr": csvr_closure,
    "csvr_2": csvr_closure_2,
    "berendsen": berendsen_closure,
}


def get_data_model(atoms, dump_steps):
    coord_size = 3 * len(atoms)
    _1d = (dump_steps,)
    _2d = (dump_steps, coord_size)

    data_model = {
        "cart_coords": _2d,
        "step": _1d,
        "energy_tot": _1d,
        "energy_pot" : _1d,
        "energy_kin": _1d,
        "energy_conserved": _1d,
        "T": _1d,
        "T_avg": _1d,
        "velocity": _2d,
    }

    return data_model


def md(
    geom,
    v0,
    steps,
    dt,
    remove_com_v=True,
    thermostat=None,
    T=T_DEFAULT,
    timecon=100,
    term_funcs=None,
    constraints=None,
    constraint_kwargs=None,
    gaussians=None,
    verbose=True,
    print_stride=50,
    dump_stride=None,
    h5_group_name="run",
):
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
    gaussians: list, optional, default=None
        List of Gaussians to be used in a metadynamics run.
    verbose : bool, default=True
        Do additional printing when True.
    print_stride : int, default=50
        Report every n-th step.
    dump_stride : int, default=None
        If given, MD data will be dumped to a HDF5 file every n-th step.
    h5_group_name = str, optional
        Name of the HDF5 group, used for dumping.
    """

    assert geom.coord_type == "cart"

    if term_funcs is None:
        term_funcs = dict()

    if dump_stride:
        dump_steps = steps // dump_stride
        data_model = get_data_model(geom.atoms, dump_steps)
        h5_group = get_h5_group("md.h5", h5_group_name, data_model, reset=True)
        h5_group.attrs["atoms"] = geom.atoms
        h5_group.attrs["dt"] = dt
        h5_group.attrs["steps"] = steps
        h5_group.attrs["masses"] = geom.masses
        h5_group.attrs["T_target"] = T

    if verbose:
        t_ps = steps * dt * 1e-3  # Total simulation time
        print(f"Doing {steps} steps of {dt:.4f} fs for a total of {t_ps:.2f} ps.\n")

    energy_forces_getter = energy_forces_getter_closure(geom)

    if gaussians is None:
        gaussians = tuple()
    gau_centers = list()
    gau_center_num = np.zeros(len(gaussians), dtype=int)
    for _, gaussian, gau_stride in gaussians:
        num_centers = steps // gau_stride
        gau_centers.append(np.zeros(num_centers))

    def update_gaussians(step, coords):
        for i, (gau_key, gaussian, gau_stride) in enumerate(gaussians):
            if (step + 1) % gau_stride == 0:
                new_center = gaussian.colvar.value(coords.reshape(-1, 3))
                gau_centers[i][gau_center_num[i]] = new_center
                gau_center_num[i] += 1
                log(
                    logger,
                    f"Added {gau_center_num[i]: >6d}. '{gau_key}' Gaussian in step {step}.",
                )
                gaussian.dump(step, gaussian.s, gaussian.w, new_center)

    def gaussian_wrapper(coords):
        E_pot, forces = energy_forces_getter(geom.coords)
        for i, (_, gaussian, _) in enumerate(gaussians):
            center_num = gau_center_num[i]
            gau_pot, gau_grad = gaussian.eval(coords, gau_centers[i][:center_num])
            # print("\tgau_pot", gau_pot)
            gau_forces = -gau_grad.flatten()
            E_pot += gau_pot
            forces += gau_forces
        return E_pot, forces

    if constraint_kwargs is None:
        constraint_kwargs = dict()

    # Fixed degrees of freedom
    fixed_dof = 0
    if remove_com_v:
        print("Removing center-of-mass velocity.")
        fixed_dof += 3
    constrained_md = constraints is not None
    # Get RATTLE function from closure for constrained MD
    if constrained_md:
        fixed_dof += len(constraints)
        rattle = rattle_closure(
            geom,
            constraints,
            dt,
            energy_forces_getter=energy_forces_getter,
            **constraint_kwargs,
        )
    print(f"Fixed degrees of freedom: {fixed_dof}")
    dof = len(geom.coords) - fixed_dof

    if thermostat is not None:
        sigma = kinetic_energy_for_temperature(len(geom.atoms), T, fixed_dof=fixed_dof)
        thermo_func = THERMOSTATS[thermostat](sigma, dof, dt=dt, tau=timecon)
    # In amu
    masses = geom.masses
    masses_rep = geom.masses_rep

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
    log(logger, f"Running MD with Δt={dt:.2f} fs for {steps} steps.")
    print()
    thermo_corr = 0.0
    for step in range(steps):
        xs.append(x.copy())

        E_kin = kinetic_energy_from_velocities(masses, v.reshape(-1, 3))
        T = temperature_for_kinetic_energy(len(masses), E_kin, fixed_dof=fixed_dof)
        T_avg += T
        Ts.append(T)
        E_tot = E_pot + E_kin
        E_tots.append(E_tot)
        E_conserved = E_pot + E_kin - thermo_corr

        status_msg = (
            f"Step {step:06d}  {t_cur*1e-3: >6.2f} ps  E={E_tot: >8.6f} E_h  "
            f"T={T: >8.2f} K <T>={T_avg/(step+1): >8.2f} K"
        )
        if step % print_stride == 0:
            log(logger, status_msg)
            if verbose:
                print(status_msg)

        if dump_stride and (step % dump_stride == 0):
            ind = step // dump_stride
            h5_group["step"][ind] = step
            h5_group["cart_coords"][ind] = x
            h5_group["velocity"][ind] = v
            h5_group["energy_tot"][ind] = E_tot
            h5_group["energy_pot"][ind] = E_pot
            h5_group["energy_kin"][ind] = E_kin
            h5_group["energy_conserved"][ind] = E_conserved
            h5_group["T"][ind] = T
            h5_group["T_avg"][ind] = T_avg / (step + 1)

        if thermostat:
            alpha = thermo_func(E_kin)
            thermo_corr += (alpha**2 - 1) * E_kin
            v *= alpha

        update_gaussians(step, x)

        # RATTLE algorithm
        if constrained_md:
            x, v, E_pot, forces = rattle(x, v, forces)
        # Simple Velocity-Verlet integration
        else:
            # E_pot, forces = energy_forces_getter(geom.coords)
            E_pot, forces = gaussian_wrapper(geom.coords)
            # Acceleration, convert from Hartree / (Bohr * amu) to Bohr/fs²
            a = forces / masses_rep * FORCE2ACC
            v += 0.5 * (a + a_prev) * dt

            if remove_com_v:
                v = remove_com_velocity(v.reshape(-1, 3), masses).flatten()
            # v*dt = Bohr/fs * fs -> Bohr
            # a*dt**2 = Bohr/fs² * fs² -> Bohr
            x += v * dt + 0.5 * a * dt ** 2
            a_prev = a

        # Update coordinates
        geom.coords = x

        for name, func in term_funcs.items():
            if func(x.reshape(-1, 3)):
                terminate = True
                terminate_key = name
                break
        if terminate:
            log(logger, f"Termination function '{name}' evaluted to True. Breaking.")
            break

        if check_for_end_sign():
            break

        # Advance time
        t_cur += dt
    log(logger, "")

    md_result = MDResult(
        coords=np.array(xs),
        t_ps=t_cur * 1e-3,
        step=step,
        terminated=terminate_key,
        T=np.array(Ts),
        E_tot=np.array(E_tots),
    )

    return md_result
