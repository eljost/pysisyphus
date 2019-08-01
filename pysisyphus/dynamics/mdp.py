#!/usr/bin/env python3

# [1] https://doi.org/10.1063/1.5082885
#     Unke, 2019

from collections import namedtuple

import numpy as np

from pysisyphus.constants import AU2KJPERMOL
from pysisyphus.dynamics.velocity_verlet import md, MDResult
from pysisyphus.dynamics.helpers import (dump_coords, get_velocities,
                                         temperature_for_kinetic_energy_au)
from pysisyphus.helpers import highlight_text


def run_md(geom, t, dt, v0=None, term_funcs=None, external=False):
    if external and hasattr(geom.calculator, "run_md"):
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
        md_result = MDResult(coords=[geom.coords for geom in geoms], t=t)
    else:
        md_kwargs = {
            "v0": v0,
            "t": t,
            "dt": dt,
            "term_funcs": term_funcs,
        }
        print("Running MD with internal implementation.")
        md_result = md(geom, **md_kwargs)

    return md_result


MDPResult = namedtuple("MDResult",
                       "ascent_xs md_init_plus md_init_minus "
                       "md_fin_plus md_fin_minus"
)


def mdp(geom, t, dt, term_funcs, t_init=None, E_excess=0.,
        displ_length=.1, epsilon=5e-4, ascent_alpha=0.05,
        max_ascent_steps=25, max_init_trajs=10, dump=True,
        seed=None, external_md=False):
    # Sanity checks and forcing some types
    dt = float(dt)
    assert dt > 0.
    t = float(t)
    assert t > dt
    if t_init is None:
        t_init = t / 10
    E_excess = float(E_excess)
    assert E_excess >= 0.
    displ_length = float(displ_length)
    assert displ_length >= 0.

    print(highlight_text("Minimum dynamic path calculation"))

    if seed is not None:
        seed = int(seed)
        print(f"Using seed {seed} to initialize the random number generator.")
        print()
        np.random.seed(seed)


    E_TS = geom.energy
    E_tot = E_TS + E_excess
    # Distribute E_excess evenly on E_pot and E_kin
    E_pot_diff = 0.5*E_excess
    E_pot_desired = E_TS + E_pot_diff

    print(f"E_TS={E_TS:.6f} au")

    # Determine transition vector
    w, v = np.linalg.eigh(geom.hessian)
    assert w[0] < -1e-8
    trans_vec = v[:,0]

    if E_excess == 0.:
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
    assert ascent_converged, "Steepest ascent didn't converge!"
    assert (E_tot - E_pot) > 0., \
         "Potential energy after steepst ascent is greater than the desired " \
        f"total energy ({E_pot:.6f} > {E_tot:.6f}). Maybe try a smaller epsilon? " \
        f"The current value Ɛ={epsilon:.6f} may be too big!"

    ascent_xs = np.array(ascent_xs)
    if dump:
        dump_coords(geom.atoms, ascent_xs, "mdp_ee_ascent.trj")
    x0 = geom.coords.copy()

    masses = geom.masses_rep

    print("Runninig initialization trajectories")
    for i in range(max_init_trajs):
        # Determine random momentum vector for the given kinetic energy
        E_kin = E_tot - E_pot
        T = temperature_for_kinetic_energy_au(len(geom.atoms), E_kin)
        v0 = get_velocities(geom, T).flatten()

        # Zero last element if we have an analytical surface
        if v0.size == 3:
            v0[2] = 0

        # Run initial MD to check if both trajectories run towards different
        # basins of attraction.

        # First MD with positive v0
        md_init_kwargs = {
            "v0": v0.copy(),
            "t": t_init,
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
        init_trajs_converged = (np.sign(p) != np.sign(m))

        if init_trajs_converged:
            break
    if dump:
        dump_coords(geom.atoms, md_init_plus.coords, "mdp_ee_init_plus.trj")
        dump_coords(geom.atoms, md_init_minus.coords, "mdp_ee_init_minus.trj")
    assert init_trajs_converged
    print(f"Ran 2*{i+1} trajectories.")
    print("Initialization completed.")
    print()

    # Run actual trajectories, using the supplied termination functions if possible.
    print("Running actual full trajectories.")

    # MD with positive v0.
    md_fin_kwargs = {
        "v0": v0.copy(),
        "t": t,
        "dt": dt,
        "term_funcs": term_funcs,
        "external": external_md,
    }
    geom.coords = x0.copy()
    md_fin_plus = run_md(geom, **md_fin_kwargs)

    geom.coords = x0.copy()
    # MD with negative v0.
    md_fin_kwargs["v0"] = -v0
    md_fin_minus = run_md(geom, **md_fin_kwargs)

    if dump:
        dump_coords(geom.atoms, md_fin_plus.coords, "mdp_ee_fin_plus.trj")
        dump_coords(geom.atoms, md_fin_minus.coords, "mdp_ee_fin_minus.trj")

    mdp_result = MDPResult(
                    ascent_xs=ascent_xs,
                    md_init_plus=md_init_plus,
                    md_init_minus=md_init_minus,
                    md_fin_plus=md_fin_plus,
                    md_fin_minus=md_fin_minus,
    )
    return mdp_result
