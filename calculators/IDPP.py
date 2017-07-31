#!/usr/bin/env python3

import numpy as np
from scipy.spatial.distance import pdist, squareform

from cos.ChainOfStates import ChainOfStates
from calculators.Calculator import Calculator
from optimizers.FIRE import FIRE


def idpp_interpolate(geometries, images_between):
    cos = ChainOfStates(geometries)
    cos.interpolate(image_num=images_between)
    images = cos.images

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
    chunk_length = images_between + 1
    for i in range(len(pdists)-1):
        from_pd = pdists[i]
        to_pd = pdists[i+1]
        pd_diff = (to_pd - from_pd) / chunk_length
        print(i, pd_diff)
        slice_start = i * chunk_length
        slice_end = slice_start + chunk_length
        slice_ = slice(slice_start, slice_end)
        images_slice = images[slice_]
        print(images_slice)
        for j, image in enumerate(images_slice):
            image.set_calculator(IDPP(from_pd + j * pd_diff))

        kwargs = {
            "max_cycles":10,
        }
        opt = FIRE(ChainOfStates(images_slice), **kwargs)
        opt.run()

    # Add last images

class IDPP(Calculator):

    def __init__(self, target): 
        self.target = target
        print("target", target)

        super(Calculator, self).__init__()

    def get_energy(self, atoms, coords):
        raise Exception("Not implemented!")

    def get_forces(self, atoms, coords):
        coords_reshaped = coords.reshape((-1, 3))
        current_pdist = pdist(coords_reshaped)
        square_dists = squareform(current_pdist)
        diff = current_pdist - self.target
        print(diff)

        energy = 0.5 * (diff**2 / current_pdist**4).sum()
        #forces = -2 * ((diff * (1 -2 * dd / current_pdist) / cur
        forces = 0.1*coords

        results = {
            "energy" : energy,
            "forces": forces
        }
        return results

    def __str__(self):
        return "IDPP calculator"
