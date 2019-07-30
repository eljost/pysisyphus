#!/usr/bin/env python3

from collections import namedtuple

import numpy as np

from pysisyphus.dynamics.velocity_verlet import md
from pysisyphus.xyzloader import make_trj_str
from pysisyphus.constants import BOHR2ANG


MDPResult = namedtuple("MDResult",
                       "ascent_xs md_init_plus md_init_minus "
                       "md_fin_plus md_fin_minus"
)


def mdp(geom, E_excess, t_init, t, dt,
        term_funcs, epsilon, ascent_alpha=0.05, max_ascent_steps=10,
        max_init_trajs=10, dump=True):
    def dump_coords(coords, trj_fn):
        coords = np.array(coords)
        coords = coords.reshape(-1, len(geom.atoms), 3) * BOHR2ANG
        trj_str = make_trj_str(geom.atoms, coords)
        if dump:
            with open(trj_fn, "w") as handle:
                handle.write(trj_str)
    E_excess = float(E_excess)
    assert E_excess >= 0.
    E_TS = geom.energy
    E_tot = E_TS + E_excess
    E_pot_desired = E_TS + 0.5*E_excess

    # print("E_TS", E_TS)
    # print("E_excess", E_excess)
    # print("E_tot", E_tot)

    # Determine transition vector
    w, v = np.linalg.eigh(geom.hessian)
    assert w[0] < -1e-8
    trans_vec = v[:,0]

    if E_excess > 0.:
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
            _ = geom.coords + step
            geom.coords = _
        assert ascent_converged, "Steepest ascent didn't converge!"
        ascent_xs = np.array(ascent_xs)
        dump_coords(ascent_xs, "ascent.trj")
        x0 = geom.coords.copy()
    else:
        # Without excess energy we have to do an initial displacement along
        # the transition vector to get a non-vanishing gradient.
        x0 = geom.coords.copy()
        # We need displacements in both directions along the transition
        # vector but right now we only use one variable x0 to store the
        # coordiantes to initiate the MD. So right now this is not possible.
        # I first have to reworkt he function.
        assert E_excess > 0., \
            "E_excess = 0 not yet supported! Please see the comment above " \
            "in the code."

    masses = geom.masses_rep

    def get_E_kin(v):
        return np.sum(masses * v**2 / 2)

    for i in range(max_init_trajs):
        # Determine random momentum vector
        v0 = np.random.rand(*trans_vec.shape)
        # Zero last element if we have an analytical surface
        if v0.size == 3:
            v0[2] = 0
        E_kin = get_E_kin(v0)
        # This determines the scaling of the kinetic energy. As we want to
        # scale the velocities we have to use the square root of the factor.
        factor = (E_tot - E_pot) / E_kin
        factor = factor**0.5
        # Scale initial velocities to yield E_kin + E_pot = E_tot
        v0 *= factor
        E_kin = get_E_kin(v0)
        np.testing.assert_allclose(E_pot + E_kin, E_tot)

        # Run initial MD to check if both trajectories run towards different
        # basins of attraction.
        geom.coords = x0.copy()
        # First MD with positive v0
        md_kwargs = {
            "v0": v0.copy(),
            "t": t_init,
            "dt": dt,
        }
        md_init_plus = md(geom, **md_kwargs)

        # Second MD with negative v0
        geom.coords = x0.copy()
        md_kwargs["v0"] = -v0.copy()
        md_init_minus = md(geom, **md_kwargs)

        # Check if both MDs run into different basins of attraction.
        # We (try to) do this by calculating the overlap between the
        # transition vector and the normalized vector defined by the
        # difference between x0 and the endpoint of the respective 
        # test trajectory. Both overlaps should yield different sings.
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
    dump_coords(md_init_plus.coords, "init_plus.trj")
    dump_coords(md_init_minus.coords, "init_minus.trj")
    assert init_trajs_converged

    geom.coords = x0.copy()
    # Run actual trajectories, using the supplied termination functions.
    # MD with positive v0.
    md_kwargs = {
        "v0": v0.copy(),
        "t": t,
        "dt": dt,
        "term_funcs": term_funcs,
    }
    md_fin_plus = md(geom, **md_kwargs)

    geom.coords = x0.copy()
    # MD with positive v0.
    md_kwargs["v0"] = -v0
    md_fin_minus = md(geom, **md_kwargs)

    dump_coords(md_fin_plus.coords, "fin_plus.trj")
    dump_coords(md_fin_minus.coords, "fin_minus.trj")

    mdp_result = MDPResult(
                    ascent_xs=ascent_xs,
                    md_init_plus=md_init_plus,
                    md_init_minus=md_init_minus,
                    md_fin_plus=md_fin_plus,
                    md_fin_minus=md_fin_minus,
    )
    return mdp_result
