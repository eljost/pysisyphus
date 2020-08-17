from math import exp, sqrt

from numpy.random import default_rng


"""
[1] https://aip.scitation.org/doi/10.1063/1.2408420
[2] https://sci-hub.tw/10.1016/j.cpc.2008.01.006
    Reformulation fo the algorithm. This is implemented e.g., in YAFF.
    https://github.com/molmod/yaff/blob/master/yaff/sampling/nvt.py

The present implementation is based on the code provided by Bussi on his
homepage (https://sites.google.com/site/giovannibussi/downloads/resamplekin.tgz)
"""


RNG = default_rng()


def sum_noises(num, rng=None):
    """
    Parameters
    ----------
    num : int
        Number of independent Gaussian noises to be squared.
    rng : numpy.random.Generator, optional
        Instances of a random number generator (RNG). If it is not provided the module-level
        RNG will be used.
    """

    if rng is None:
        rng = RNG

    if num == 0:
        sum_ = 0.0
    elif num == 1:
        sum_ = rng.normal()**2
    # nn even, dof - 1 odd
    elif (num % 2) == 0:
        sum_ = 2.0 * rng.gamma(shape=num/2)
    # nn odd, dof - 1 even
    else:
        sum_ = 2.0 * rng.gamma(shape=(num-1)/2) + rng.normal()**2

    return sum_


def resample_kin(cur_kinetic_energy, sigma, dof, tau_t, rng=None):
    """
    Parameters
    ----------
    cur_kinetic_energy : float
        Present value of the kinetic energy of the atoms to be thermalized
        in arbitrary units.
    sigma : float
        Target average value of the kinetic energy (1/2 dof k_b T) in the same units as
        cur_kinetic_energy.
    dof : int
        Degrees of freedom.
    tau_t : float
        Relaxation time of the thermostat in units of "how often this routine is called"
        (dt / timeconstant).
    rng : numpy.random.Generator, optional
        Instances of a random number generator (RNG). If it is not provided the module-level
        RNG will be used.
    """
    if tau_t > 0.1:
        factor = exp(-1.0 / tau_t)
    else:
        factor = 0.0

    if rng is None:
        rng = RNG

    rr = rng.normal()
    new_kinetic_energy = (
        cur_kinetic_energy
        + (1.0 - factor)
        * (sigma * (sum_noises(dof-1) + rr**2) / dof - cur_kinetic_energy)
        + 2.0 * rr * sqrt(cur_kinetic_energy * sigma / dof * (1.0 - factor) * factor)
    )
    return new_kinetic_energy
