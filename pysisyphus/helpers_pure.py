import numpy as np

from pysisyphus.constants import AU2J, AMU2KG, BOHR2M


"""Functions defined here don't import anything from pysisyphus, besides
the constants module, but only from the stdlib and from third parties."""


def eigval_to_wavenumber(ev):
    conv = AU2J/(AMU2KG*BOHR2M**2)

    return np.sign(ev) * np.sqrt(np.abs(ev)*conv)/(2*np.pi*3e10)
