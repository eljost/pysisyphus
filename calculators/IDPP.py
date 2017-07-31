#!/usr/bin/env python3

import numpy as np
from scipy.spatial.distance import pdist

from calculators.Calculator import Calculator

def idpp_interpolate(geometries, images_between):
    coords = [geometry.coords.reshape((-1, 3)) for geometry in geometries]
    print("coords")
    for c in coords:
        print(c)
    pdists = [pdist(c) for c in coords]
    print("pdists")
    for pd in pdists:
        print(pd)
    """
    from scipy.spatial.distance squareform
    squareforms = [squareform(pd) for pd in pdists]
    print("squareforms")
    for sq in squareforms:
        print(sq)
    """
    # Iterate over pairs of distance matrices and interpolate between them
    for i in range(len(pdists)-1):
        from_pd = pdists[i]
        to_pd = pdists[i+1]
        pd_diff = (to_pd - from_pd) / (images_between + 1)
        print(i, pd_diff)
        for i, geometry in enumerate(geometries):
            geometry.set_calculator(IDPP(from_pd + i * pd_diff))

class IDPP(Calculator):

    def __init__(self, target): 
        self.target = target
        print("target", target)

        super(Calculator, self).__init__()

    def get_energy(self, atoms, coords):
        raise Exception("Not implemented!")

    def get_forces(self, atoms, coords):
        raise Exception("Not implemented!")

    def __str__(self):
        return "IDPP calculator"
