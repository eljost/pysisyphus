#!/usr/bin/env python3

from scipy.spatial.distance import pdist

# from pysisyphus.constants import BOHR2ANG, ANG2BOHR
from pysisyphus.calculators.IDPPCalculator import IDPPCalculator
from pysisyphus.constants import BOHR2ANG, ANG2BOHR
from pysisyphus.cos.NEB import NEB
from pysisyphus.helpers import align_geoms
from pysisyphus.optimizers.FIRE import FIRE
from pysisyphus.interpolate.Interpolator import Interpolator


# [1] http://aip.scitation.org/doi/full/10.1063/1.4878664
# See https://gitlab.com/ase/ase/blob/master/ase/neb.py


class IDPP(Interpolator):

    def interpolate(self, initial_geom, final_geom):
        # Do an initial linear interpolation to generate all geometries/images
        # that will be refined later by IDPP interpolation.
        linear_interpol = super().interpolate(initial_geom, final_geom)
        idpp_geoms = [initial_geom] + linear_interpol + [final_geom]
        align_geoms(idpp_geoms)

        # Interestingly IDPP calculations work much better when done
        # in Angstroem instead of in Bohr.
        for geom in idpp_geoms:
            geom.coords *= BOHR2ANG

        # We want to interpolate between these two condensed distance matrices
        initial_pd = pdist(initial_geom.coords3d)
        final_pd = pdist(final_geom.coords3d)
        steps = 1 + self.between
        pd_diff = (final_pd - initial_pd) / steps

        for i, geom in enumerate(idpp_geoms):
            geom.set_calculator(IDPPCalculator(initial_pd + i * pd_diff))

        neb = NEB(idpp_geoms, parallel=False, fix_ends=True)
        opt_kwargs = {
            "max_cycles": 1000,
            "rms_force": 1e-2,
            "keep_cycles": False,
            "align": False,
        }
        opt = FIRE(neb, **opt_kwargs)
        opt.run()

        for geom in idpp_geoms:
            # Delete IDPP calculator, energies and forces
            geom.clear()
            geom.coords *= ANG2BOHR

        interpolated_geoms = idpp_geoms[1:-1]
        return interpolated_geoms
