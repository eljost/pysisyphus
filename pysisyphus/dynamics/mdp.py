# [1] https://doi.org/10.1063/1.5082885
#     Minimum dynamic path
#     Unke, 2019

from collections import namedtuple
import operator
import re

import numpy as np

from pysisyphus.constants import AU2KJPERMOL
from pysisyphus.dynamics.driver import md, MDResult
from pysisyphus.dynamics.helpers import (
    dump_coords,
    get_mb_velocities_for_geom,
    temperature_for_kinetic_energy,
)
from pysisyphus.helpers_pure import highlight_text


def parse_raw_term_func(raw_term_func):
    funcs = {
        "<": operator.lt,
        ">": operator.gt,
        "<=": operator.le,
        ">=": operator.ge,
        "==": operator.eq,
    }

    def comp_closure(indices, op, ref_value):
        a_ind, b_ind = indices
        func = funcs[op]

        def comp_func(coords3d):
            a = coords3d[a_ind]
            b = coords3d[b_ind]
            dist = np.linalg.norm(a - b)
            return func(dist, ref_value)

        return comp_func

    operator_re = re.compile("([<>=]+)")
    mobj = operator_re.split(raw_term_func)
    if mobj is None:
        print(f"Could not parse term_func '{raw_term_func}!'")
        return None
    indices, op, ref_value = mobj
    ref_value = float(ref_value)
    indices = [int(ind) for ind in indices.split(",")]
    return comp_closure(indices, op, ref_value)


def parse_raw_term_funcs(raw_term_funcs):
    comp_funcs = {}
    for k, v in raw_term_funcs.items():
        comp_func = parse_raw_term_func(v)
        if comp_func:
            comp_funcs[k] = comp_func
    return comp_funcs


def run_md(geom, dt, steps, v0=None, term_funcs=None, external=False):
    if external and hasattr(geom.calculator, "run_md"):
        t = dt * steps
        t_ps = t * 1e-3
        md_kwargs = {
            "atoms": geom.atoms,
            "coords": geom.coords,
            "t": t,
            "dt": dt,
            "velocities": v0,
            "dump": dt,
        }
        print("Running MD with external calculator implementation.")
        if term_funcs is not None:
            print("Termination functions are not supported in external MD!")
        geoms = geom.calculator.run_md(**md_kwargs)
        md_result = MDResult(
            coords=[geom.coords for geom in geoms],
            t_ps=t_ps,
            step=int(t / dt - 1),
            terminated=None,
            T=None,
            E_tot=None,
        )
    else:
        md_kwargs = {
            "v0": v0,
            "steps": steps,
            "dt": dt,
            "term_funcs": term_funcs,
            "verbose": False,
            "remove_com_v": False,
        }
        print("Running MD with internal implementation.")
        md_result = md(geom, **md_kwargs)

    return md_result


MDPResult = namedtuple(
    "MDResult",
    "ascent_xs md_init_plus md_init_minus "
    "md_fin_plus md_fin_minus "
    "md_fin_plus_term md_fin_minus_term",
)


def mdp(
    geom,
    steps,
    dt,
    term_funcs=None,
    steps_init=None,
    E_excess=0.0,
    displ_length=0.1,
    epsilon=5e-4,
    ascent_alpha=0.05,
    max_ascent_steps=25,
    max_init_trajs=10,
    dump=True,
    seed=None,
    external_md=False,
):
    # Sanity checks and forcing some types
    dt = float(dt)
    assert dt > 0.0
    steps = int(steps)
    t = dt * steps
    # assert t > dt
    if steps_init is None:
        steps_init = steps // 10
        print(f"No 'steps_init' provided! Using {steps_init}")
    E_excess = float(E_excess)
    assert E_excess >= 0.0
    displ_length = float(displ_length)
    assert displ_length >= 0.0
    if term_funcs is None:
        term_funcs = {}
    for k, v in term_funcs.items():
        if callable(v):
            continue
        elif isinstance(v, str):
            term_funcs[k] = parse_raw_term_func(v)
        else:
            raise Exception(f"Invalid term function '{k}: {v}' encountered!")

    print(highlight_text("Minimum dynamic path calculation"))

    if seed is None:
        # 2**32 - 1
        seed = np.random.randint(4294967295)
    np.random.seed(seed)
    print(f"Using seed {seed} to initialize the random number generator.\n")

    E_TS = geom.energy
    E_tot = E_TS + E_excess
    # Distribute E_excess evenly on E_pot and E_kin
    E_pot_diff = 0.5 * E_excess
    E_pot_desired = E_TS + E_pot_diff

    print(f"E_TS={E_TS:.6f} au")

    # Determine transition vector
    w, v = np.linalg.eigh(geom.hessian)
    assert w[0] < -1e-8
    trans_vec = v[:, 0]

    # Disable removal of translation/rotation for analytical potentials
    remove_com_v = remove_rot_v = geom.cart_coords.size > 3

    if E_excess == 0.0:
        print("MDP without excess energy.")
        # Without excess energy we have to do an initial displacement along
        # the transition vector to get a non-vanishing gradient.
        initial_displacement = displ_length * trans_vec
        x0_plus = geom.coords + initial_displacement
        x0_minus = geom.coords - initial_displacement

        v0_zero = np.zeros_like(geom.coords)
        md_kwargs = {
            "v0": v0_zero.copy(),
            "t": t,
            "dt": dt,
            "term_funcs": term_funcs,
            "external": external_md,
        }

        geom.coords = x0_plus
        md_fin_plus = run_md(geom, **md_kwargs)

        geom.coords = x0_minus
        md_fin_minus = run_md(geom, **md_kwargs)

        if dump:
            dump_coords(geom.atoms, md_fin_plus.coords, "mdp_plus.trj")
            dump_coords(geom.atoms, md_fin_minus.coords, "mdp_minus.trj")

        mdp_result = MDPResult(
            ascent_xs=None,
            md_init_plus=None,
            md_init_minus=None,
            md_fin_plus=md_fin_plus,
            md_fin_minus=md_fin_minus,
        )
        return mdp_result

    print(f"E_excess={E_excess:.6f} au, ({E_excess*AU2KJPERMOL:.1f} kJ/mol)")
    print(f"E_pot,desired=E_TS + {E_pot_diff*AU2KJPERMOL:.1f} kJ/mol")
    print()

    # Generate random vector perpendicular to transition vector
    perp_vec = np.random.rand(*trans_vec.shape)
    # Zero last element if we have an analytical surface
    if perp_vec.size == 3:
        perp_vec[2] = 0
    # Orthogonalize vector
    perp_vec = perp_vec - (perp_vec @ trans_vec) * trans_vec
    perp_vec /= np.linalg.norm(perp_vec)

    # Initial displacement from x_TS to x, generating a point with
    # non-vanishing gradient.
    x = geom.coords + epsilon * perp_vec
    geom.coords = x

    # Do steepest ascent until E_tot is reached
    E_pot = geom.energy
    ascent_xs = list()
    for i in range(max_ascent_steps):
        ascent_xs.append(geom.coords.copy())
        ascent_converged = E_pot >= E_pot_desired
        if ascent_converged:
            break
        gradient = geom.gradient
        E_pot = geom.energy

        direction = gradient / np.linalg.norm(gradient)
        step = ascent_alpha * direction
        new_coords = geom.coords + step
        geom.coords = new_coords

    # calc = geom.calculator
    # class Opt:
    # pass
    # _opt = Opt()
    # _opt.coords = np.array(ascent_xs)
    # calc.plot_opt(_opt, show=True)

    assert ascent_converged, "Steepest ascent didn't converge!"
    assert (E_tot - E_pot) > 0.0, (
        "Potential energy after steepst ascent is greater than the desired "
        f"total energy ({E_pot:.6f} > {E_tot:.6f}). Maybe try a smaller epsilon? "
        f"The current value Ɛ={epsilon:.6f} may be too big!"
    )

    ascent_xs = np.array(ascent_xs)
    if dump:
        dump_coords(geom.atoms, ascent_xs, "mdp_ee_ascent.trj")
    x0 = geom.coords.copy()

    print(highlight_text("Runninig initialization trajectories", level=1))
    for i in range(max_init_trajs):
        # Determine random momentum vector for the given kinetic energy
        E_kin = E_tot - E_pot
        T = temperature_for_kinetic_energy(len(geom.atoms), E_kin)
        v0 = get_mb_velocities_for_geom(
            geom, T, remove_com_v=remove_com_v, remove_rot_v=remove_rot_v
        ).flatten()

        # Zero last element if we have an analytical surface
        if v0.size == 3:
            v0[2] = 0

        # Run initial MD to check if both trajectories run towards different
        # basins of attraction.

        # First MD with positive v0
        md_init_kwargs = {
            "v0": v0.copy(),
            "steps": steps_init,
            "dt": dt,
            "external": external_md,
        }
        geom.coords = x0.copy()
        md_init_plus = run_md(geom, **md_init_kwargs)
        # Second MD with negative v0
        geom.coords = x0.copy()
        md_init_kwargs["v0"] = -v0.copy()
        md_init_minus = run_md(geom, **md_init_kwargs)

        dump_coords(geom.atoms, md_init_plus.coords, f"mdp_ee_init_plus_{i:02d}.trj")
        dump_coords(geom.atoms, md_init_minus.coords, f"mdp_ee_init_minus_{i:02d}.trj")

        # Check if both MDs run into different basins of attraction.
        # We (try to) do this by calculating the overlap between the
        # transition vector and the normalized vector defined by the
        # difference between x0 and the endpoint of the respective
        # test trajectory. Both overlaps should have different sings.
        end_plus = md_init_plus.coords[-1]
        pls = end_plus - x0
        pls /= np.linalg.norm(pls)
        end_minus = md_init_minus.coords[-1]
        minus = end_minus - x0
        minus /= np.linalg.norm(minus)
        p = trans_vec @ pls
        m = trans_vec @ minus
        init_trajs_converged = np.sign(p) != np.sign(m)

        if init_trajs_converged:
            print("Trajectories ran into different basins. Breaking.")
            break
    if dump:
        dump_coords(geom.atoms, md_init_plus.coords, "mdp_ee_init_plus.trj")
        dump_coords(geom.atoms, md_init_minus.coords, "mdp_ee_init_minus.trj")
    assert init_trajs_converged
    print(f"Ran 2*{i+1} initialization trajectories.")
    print()

    # Run actual trajectories, using the supplied termination functions if possible.
    print(highlight_text("Running actual full trajectories.", level=1))

    def print_status(terminated, step):
        if terminated:
            msg = f"\tTerminated by '{terminated}' in step {step}."
        else:
            msg = "\tMax time steps reached!"
        print(msg)

    # "Production"/Final MDs
    md_fin_kwargs = {
        "v0": v0.copy(),
        "steps": steps,
        "dt": dt,
        "term_funcs": term_funcs,
        "external": external_md,
    }
    # MD with positive v0.
    geom.coords = x0.copy()
    md_fin_plus = run_md(geom, **md_fin_kwargs)
    print_status(md_fin_plus.terminated, md_fin_plus.step)

    # MD with negative v0.
    geom.coords = x0.copy()
    md_fin_kwargs["v0"] = -v0
    md_fin_minus = run_md(geom, **md_fin_kwargs)
    print_status(md_fin_minus.terminated, md_fin_minus.step)

    md_fin_plus_term = md_fin_plus.terminated
    md_fin_minus_term = md_fin_minus.terminated

    if dump:
        dump_coords(geom.atoms, md_fin_plus.coords, "mdp_ee_fin_plus.trj")
        dump_coords(geom.atoms, md_fin_minus.coords, "mdp_ee_fin_minus.trj")

    mdp_result = MDPResult(
        ascent_xs=ascent_xs,
        md_init_plus=md_init_plus,
        md_init_minus=md_init_minus,
        md_fin_plus=md_fin_plus,
        md_fin_minus=md_fin_minus,
        md_fin_plus_term=md_fin_plus_term,
        md_fin_minus_term=md_fin_minus_term,
    )
    return mdp_result
