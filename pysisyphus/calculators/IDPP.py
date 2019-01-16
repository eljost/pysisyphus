#!/usr/bin/env python3

import numpy as np
from scipy.spatial.distance import pdist, squareform

from pysisyphus.calculators.Calculator import Calculator
from pysisyphus.constants import BOHR2ANG, ANG2BOHR
from pysisyphus.cos.ChainOfStates import ChainOfStates
from pysisyphus.cos.NEB import NEB
from pysisyphus.optimizers.BFGS import BFGS
from pysisyphus.optimizers.FIRE import FIRE
from pysisyphus.helpers import procrustes


# [1] http://aip.scitation.org/doi/full/10.1063/1.4878664
# See https://gitlab.com/ase/ase/blob/master/ase/neb.py

def idpp_interpolate(geometries, images_between, keep_cycles=False):
    # Do an initial linear interpolation to generate all geometries/images
    # that will be refined later by IDPP interpolation.
    cos = ChainOfStates(geometries)
    procrustes(cos)
    cos.interpolate(image_num=images_between)
    images = cos.images

    # Interestingly IDPP calculations work much better when done
    # in Angstroem instead of in Bohr.
    for image in cos.images:
        image.coords *= BOHR2ANG

    coords = [geometry.coords.reshape((-1, 3)) for geometry in geometries]
    # Calculate the condensed distances matrices
    pdists = [pdist(c) for c in coords]

    # A chunk consists of an initial image and the interpolated images.
    chunk_length = 1 + images_between
    # We don't consider the last image because it gets appended at the end,
    # or if more than two initial geometries are present it gets added as
    # first image of the next chunk.
    idpp_interpolated_images = list()
    for i in range(len(pdists)-1):
        # We want to interpolate between these two condensed distance matrices
        from_pd = pdists[i]
        to_pd = pdists[i+1]
        pd_diff = (to_pd - from_pd) / chunk_length

        slice_start = i * chunk_length
        slice_end = slice_start + chunk_length
        slice_ = slice(slice_start, slice_end)
        images_slice = images[slice_]
        for j, image in enumerate(images_slice):
            image.set_calculator(IDPP(from_pd + j * pd_diff))

        kwargs = {
            "max_cycles":1000,
            "rms_force": 1e-2,
            "keep_cycles": keep_cycles,
            "align": False,
        }
        # opt = FIRE(NEB(images_slice, parallel=False), **kwargs)
        opt = BFGS(NEB(images_slice, parallel=False), **kwargs)
        opt.run()
        idpp_interpolated_images.extend(images_slice)

    # Add last image
    idpp_interpolated_images.append(images[-1])
    assert(len(idpp_interpolated_images) == len(images))

    for image in idpp_interpolated_images:
        # Delete IDPP calculator, energies and forces
        image.clear()
        image.coords *= ANG2BOHR

    return idpp_interpolated_images


class IDPP(Calculator):

    def __init__(self, target): 
        self.target = squareform(target)

        super(IDPP, self).__init__(base_name="idpp")

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

        # Adapted from ASE IDPP calculator
        # https://gitlab.com/ase/ase/blob/master/ase/neb.py, GPL2
        forces = -2 * ((curr_diff *
                       (1 - 2 * curr_diff / curr_square) /
                        curr_square**5)[...,np.newaxis] * D).sum(0)

        results = {
            "energy" : energy,
            "forces": forces.flatten()
        }
        return results

    def __str__(self):
        return "IDPP calculator"
