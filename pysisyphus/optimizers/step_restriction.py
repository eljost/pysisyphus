import numpy as np


def scale_by_max_step(steps, max_step):
    steps_max = np.abs(steps).max()
    if steps_max > max_step:
        steps *= max_step / steps_max
    return steps
