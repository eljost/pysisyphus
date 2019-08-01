#!/usr/bin/env python3

from collections import namedtuple

import numpy as np

from pysisyphus.constants import BOHR2M, AMU2KG, AU2J, AU2MPERSEC


MDResult = namedtuple("MDResult",
                      "coords t",
)


# def md(geom, v0, t, dt, term_funcs=None):
    # """TODO: dump coords, velocities; check energy conservation.
    # Align geometry to avoid drift."""
    # t = t*FS2AU
    # dt = dt*FS2AU

    # if term_funcs is None:
        # term_funcs = list()

    # m = geom.masses_rep
    # x = geom.coords
    # v = v0
    # a_prev = 0
    # xs = list()

    # t_cur = 0
    # while t_cur < t:
        # t_cur += dt
        # xs.append(x.copy())
        # f = geom.forces
        # a = f / m
        # v += .5 * (a + a_prev) * dt
        # x += v*dt + .5*a*dt**2
        # geom.coords = x
        # a_prev = a
        # if any([tf(x) for tf in term_funcs]):
            # break
    # md_result = MDResult(
                    # coords=np.array(xs),
                    # t=t,
    # )
    # return md_result


def md(geom, v0, t, dt, term_funcs=None):
    """TODO: dump coords, velocities; check energy conservation.
    Align geometry to avoid drift."""
    steps = int(t/dt)
    print(f"Doing {steps} steps of {dt:.1f} fs for a total of {steps*dt:.1f} fs.")

    # Convert to seconds
    t = t  * 1e-15
    dt = dt * 1e-15

    if term_funcs is None:
        term_funcs = list()

    m = geom.masses_rep * AMU2KG
    x = geom.coords * BOHR2M
    v = v0 * AU2MPERSEC
    a_prev = 0
    xs = list()

    t_cur = 0
    while t_cur < t:
        t_cur += dt
        xs.append(x.copy())
        f = geom.forces * AU2J / BOHR2M
        a = f / m
        v += .5 * (a + a_prev) * dt
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
