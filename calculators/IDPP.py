#!/usr/bin/env python3

import numpy as np
from scipy.spatial.distance import pdist, squareform

from cos.ChainOfStates import ChainOfStates
from cos.NEB import NEB
from calculators.Calculator import Calculator
from optimizers.FIRE import FIRE

# [1] http://aip.scitation.org/doi/full/10.1063/1.4878664


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
    squareforms = [squareform(pd) for pd in pdists]
    print("squareforms")
    for sq in squareforms:
        print(sq)
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
            "max_cycles":100,
        }
        opt = FIRE(NEB(images_slice), **kwargs)
        opt.run()

    # Add last images

class IDPP(Calculator):

    def __init__(self, target): 
        self.target = squareform(target)

        super(Calculator, self).__init__()

    def get_energy(self, atoms, coords):
        raise Exception("Not implemented!")

    def get_forces(self, atoms, coords):
        coords_reshaped = coords.reshape((-1, 3))

        D = []
        for c in coords_reshaped:
            Di = coords_reshaped - c
            D.append(Di)
        D = np.array(D)

        curr_pdist = pdist(coords_reshaped)
        curr_square = squareform(curr_pdist)
        curr_diff = curr_square - self.target

        curr_square = curr_square + np.eye(curr_square.shape[0])


        # The bigger the differences 'curr_diff', the bigger the energy.
        # The smaller the current distances 'current_pdist', the bigger
        # the energy.
        energy = 0.5 * (curr_diff**2 / curr_square**4).sum()

        forces = -2 * ((curr_diff *
                       (1 - 2 * curr_diff / curr_square) /
                        curr_square**5)[...,np.newaxis] * D).sum(0)
        #print("forces")
        #print(forces)

        results = {
            "energy" : energy,
            "forces": forces.flatten()
        }
        return results

    def __str__(self):
        return "IDPP calculator"
