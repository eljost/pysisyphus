#!/usr/bin/env python3

import numpy as np

def abls(cur_forces, prev_forces):
    """Accelerated backtracking line search."""
    epsilon = 1e-3
    alpha0 = 0.1
    gamma = 0.5
    Nback = 0
    N0 = 0

    chk = (cur_forces - prev_forces) / np.linalg.norm(cur_forces+prev_forces)
    skip = False

    if chk > epsilon:
        # Slow alpha
        alpha *= gamma
        skip = True
        Nback = N0
    else:
        Nback -= 1
        if Nback < 0:
            Nback = N0
            if alpha < alpha0:
                # Reset alpha
                alpha = alpha0
                skip = True
            else:
                # Accelerate alpha
                alpha /= gamma
    return alpha, skip
