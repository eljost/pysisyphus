#!/usr/bin/env python3

import numpy as np


def gram_schmidt(vecs):
    def proj(v1, v2): return v1.dot(v2)/v1.dot(v1)
    ortho = [vecs[0], ]
    for v1 in vecs[1:]:
        tmp = v1.copy()
        for v2 in ortho:
            tmp -= proj(v2, v1)*v2
        ortho.append(tmp / np.linalg.norm(tmp))
    return ortho
