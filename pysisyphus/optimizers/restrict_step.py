import numpy as np


def scale_by_max_step(steps, max_step):
    steps_max = np.abs(steps).max()
    if steps_max > max_step:
        steps *= max_step / steps_max
    return steps


def get_scale_max(max_element):
    def scale_max(step):
        step_max = np.abs(step).max()
        if step_max > max_element:
            step *= max_element / step_max
        return step
    return scale_max
