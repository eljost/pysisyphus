import matplotlib.pyplot as plt
import numpy as np
from scipy import sparse

from pysisyphus.Geometry import Geometry
from pysisyphus.helpers import geom_loader
from pysisyphus.intcoords.setup_fast import find_bonds


def test_reverse_cuthill_mckee():
    geom = geom_loader("lib:thr75_from_1bl8.xyz")
    inds = np.arange(len(geom.atoms))
    np.random.shuffle(inds)
    ats = [geom.atoms[i] for i in inds]
    c3d = geom.coords3d.copy()
    # import pdb; pdb.set_trace()
    c3d = c3d[inds]
    geom_ = Geometry(ats, c3d)
    # geom = geom_loader("lib:h2o.xyz")
    # geom.jmol()
    bonds = find_bonds(geom_)
    # data = np.ones(bonds.shape[0])
    data = np.ones(bonds.size)
    rows, cols = bonds.T
    dim = len(geom_.atoms)
    rows_cols = np.concatenate((rows, cols), axis=0)
    cols_rows = np.concatenate((cols, rows), axis=0)
    sp = sparse.csc_matrix((data, (rows_cols, cols_rows)), shape=(dim, dim))
    perm = sparse.csgraph.reverse_cuthill_mckee(sp, symmetric_mode=True)

    sp_permd = sp[perm[:,None], perm[None,:]]

    # fig, (ax0, ax1) = plt.subplots(ncols=2)
    # ax0.imshow(sp.A)
    # ax0.set_title("org")
    # _ = sparse.csgraph.reverse_cuthill_mckee(sp, symmetric_mode=True)
    # ax1.imshow(sp_permd.A)
    # ax1.set_title("sorted")
    # plt.show()
