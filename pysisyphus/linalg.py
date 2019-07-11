#!/usr/bin/env python3

import numpy as np


def gram_schmidt(vecs, thresh=1e-8):
    def proj(v1, v2): return v1.dot(v2)/v1.dot(v1)
    ortho = [vecs[0], ]
    for v1 in vecs[1:]:
        tmp = v1.copy()
        for v2 in ortho:
            tmp -= proj(v2, v1)*v2
        norm = np.linalg.norm(tmp)
        # Don't append linear dependent vectors, as their norm will be
        # near zero. Renormalizing them to unity would lead to numerical
        # garbage and to erronous results later on, when we orthgonalize
        # against this 'arbitrary' vector.
        if norm <= thresh:
            continue
        ortho.append(tmp / norm)
    return np.array(ortho)
