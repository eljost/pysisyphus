from math import exp, sqrt

import numpy as np
from numpy.random import default_rng

from pysisyphus.constants import KBAU


"""
[1] https://aip.scitation.org/doi/10.1063/1.2408420
[2] https://dx.doi.org/10.1016/j.cpc.2008.01.006
    Reformulation fo the algorithm. This is implemented e.g., in YAFF.
    https://github.com/molmod/yaff/blob/master/yaff/sampling/nvt.py

csvr_closure() is based on the implementation provided on Bussis homepage:
    https://sites.google.com/site/giovannibussi/downloads/resamplekin.tgz
    (At least in my implementation) there seems to be a problem with
    the conserved quantity, which is not conserved at all ...
csvr_closure_2() is based on [2]
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


def csvr_closure(sigma, dof, dt, tau=100, rng=None):
    """
    Parameters
    ----------
    sigma : float
        Target average value of the kinetic energy (1/2 dof k_b T) in the same units as
        cur_kinetic_energy.
    dof : int
        Degrees of freedom.
    tau : float
        Timeconstant of the thermostat.    tau : float
        Timeconstant of the thermostat.
    rng : numpy.random.Generator, optional
        Instances of a random number generator (RNG). If it is not provided the module-level
        RNG will be used.
    """
    # Relaxation time of the thermostat in units of "how often this routine
    # is called" (dt / timeconstant).
    tau_t = dt / tau

    if tau_t > 0.1:
        factor = exp(-1.0 / tau_t)
    else:
        factor = 0.0

    if rng is None:
        rng = RNG

    def resample_kin(cur_kinetic_energy):
        """
        Parameters
        ----------
        cur_kinetic_energy : float
            Present value of the kinetic energy of the atoms to be thermalized
            in arbitrary units.
        """
        rr = rng.normal()
        new_kinetic_energy = (
            cur_kinetic_energy
            + (1.0 - factor)
            * (sigma * (sum_noises(dof-1) + rr**2) / dof - cur_kinetic_energy)
            + 2.0 * rr * sqrt(cur_kinetic_energy * sigma / dof * (1.0 - factor) * factor)
        )
        alpha = sqrt(new_kinetic_energy / cur_kinetic_energy)
        return alpha
    return resample_kin


def csvr_closure_2(sigma, dof, dt, tau=100, rng=None):
    if rng is None:
        rng = RNG
    c = exp(-dt / tau)

    def resample_kin(cur_kinetic_energy):
        """Canonical stocastical velocity rescaling.

        See dx.doi.org/10.1016/j.cpc.2008.01.006
        """
        R = rng.normal()
        S = np.sum(rng.normal(size=dof-1)**2)
        quot = (1 - c) * sigma / (dof * cur_kinetic_energy)
        alpha =  sqrt(c + quot * (S + R**2) + 2 * R * sqrt(c*quot))
        sign = np.sign(R + sqrt(c / quot))
        return sign * alpha
    return resample_kin


def berendsen_closure(sigma, dof, dt, tau=100, rng=None):
    """ https://doi.org/10.1063/1.448118"""
    tau_t = dt / tau

    def resample_kin(cur_kinetic_energy):
      alpha = sqrt(1 + tau_t * (sigma / cur_kinetic_energy - 1))
      return alpha
    return resample_kin
