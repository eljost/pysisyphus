#!/usr/bin/env python3

from collections import namedtuple

import numpy as np


MDResult = namedtuple("MDResult",
                      "coords t",
)


def md(geom, v0, t, dt, term_funcs=None):
    if term_funcs is None:
        term_funcs = list()
    m = geom.masses_rep
    x = geom.coords
    v = v0
    a_prev = 0
    t_cur = 0
    xs = list()
    while t_cur < t:
        t_cur += dt
        xs.append(x.copy())
        f = geom.forces
        a = f / m
        v += .5 * (a + a_prev) * dt
        x += v*dt + .5*a*dt**2
        geom.coords = x
        a_prev = a
        if any([tf(x) for tf in term_funcs]):
            break
    md_result = MDResult(
                    coords=np.array(xs),
                    t=t,
    )
    return md_result
