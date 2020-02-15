import numpy as np

from pysisyphus.helpers import geom_from_library
from pysisyphus.optimizers.PreconLBFGS import PreconLBFGS
from pysisyphus.calculators import XTB

np.set_printoptions(suppress=True, precision=4, linewidth=120)


def run():
    # geom = geom_from_library("h2o.xyz")
    geom = geom_from_library("h2o_shaken.xyz")
    geom.set_calculator(XTB())
    en = geom.energy
    grad = geom.gradient
    import pdb; pdb.set_trace()

    opt = PreconLBFGS(geom, thresh="gau_tight")
    opt.run()


def verify_P():
    geom = geom_from_library("h2o_shaken.xyz")
    from pysisyphus.InternalCoordinates import RedundantCoords
    int_ = RedundantCoords(geom.atoms, geom.cart_coords)
    bonds = int_.bond_indices
    # bends = int_.bending_indices
    bends = list()

    from pysisyphus.optimizers.precon import get_lindh_k, get_precon

    ks = get_lindh_k(geom.atoms, geom.coords3d, bonds, bends)

    from ase import Atoms
    from ase.utils import ff 
    from ase.optimize.precon.precon import FF

    from pysisyphus.constants import BOHR2ANG
    # Convert from hartree/bohr² to eV/ang²
    ks_ase = np.array(ks) / 0.0103

    from pysisyphus.intcoords.findbonds import get_pair_covalent_radii
    from scipy.spatial.distance import squareform

    pair_cov_radii = get_pair_covalent_radii(geom.atoms)
    pcr_square = squareform(pair_cov_radii)
    c3d = geom.coords3d
    reqs = [np.linalg.norm(c3d[i]-c3d[j]) for i, j in bonds]
    reqs = np.array(reqs) * BOHR2ANG
    # bond_lengths = np.array([pcr_square[i, j] for i, j in bonds])
    # bond_lengths *= BOHR2ANG

    # ase_bonds = [ff.Bond(i, j, k, b0) for (i, j), k, b0
                 # in zip(bonds, ks_ase, bond_lengths)]
    ase_bonds = [ff.Bond(i, j, k, b0) for (i, j), k, b0
                 in zip(bonds, ks_ase, reqs)]
    # import pdb; pdb.set_trace()
    ase_precon = FF(apply_positions=True, apply_cell=False, bonds=ase_bonds)

    ang_coords = geom.coords3d*BOHR2ANG
    atoms = Atoms(symbols=geom.atoms, positions=ang_coords, cell=(10.,10.,10.))
    ase_P = ase_precon.make_precon(atoms)

    ase_P_ = ase_P.todense() * 0.0103
    wa, va = np.linalg.eigh(ase_P_)

    P3 = get_precon(geom.atoms, geom.coords, bonds, bends)
    P3_ = P3.todense()
    import pdb; pdb.set_trace()


if __name__ == "__main__":
    # run()
    verify_P()
