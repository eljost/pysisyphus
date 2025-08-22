import numpy as np
from scipy.stats import special_ortho_group


def get_random_U(N):
    """Get random rotation matrix."""
    return special_ortho_group.rvs(N)


def fmt_tensor(tensor):
    return np.array2string(tensor, formatter={"float": lambda f: f"{f: >10.4e}"})
