import pytest

from pysisyphus.helpers import geom_loader
from pysisyphus.intcoords.setup import get_fragments
import pysisyphus.linalg_affine3 as aff3
from pysisyphus.wrapper import jmol


def test_reflect_cycle():
    """Reflect phenyl-fragment on plane containing the other atoms."""

    geom = geom_loader("lib:affine3_input.xyz")
    # geom.jmol()
    atoms = geom.atoms
    c3d = geom.coords3d

    # Determine fragments and associated coordinates
    frag1, frag2 = get_fragments(atoms, geom.coords3d)
    # print(f"{len(frag1)=}")
    # print(f"{len(frag2)=}")
    frag1 = list(frag1)
    frag2 = list(frag2)
    c3d1 = c3d[frag1]
    c3d2 = c3d[frag2]

    # Determine normal of best-fit plane containing fragment 1
    n = aff3.svd_plane_normal(c3d1)
    # Jmol arrow
    point = c3d1[0]
    stdin = jmol.arrow(point, point + 4 * n)

    # Determine reflection matrix and reflect second fragment.
    R = aff3.reflection_matrix(n, c3d1[0])
    rc3d2 = aff3.reflect_coords(R, c3d2)

    # from pysisyphus.Geometry import Geometry
    # atoms2 = [atoms[i] for i in frag2]
    # frag_geom = Geometry(atoms2, rc3d2)
    # union = geom + frag_geom
    # union.jmol(stdin=stdin)

    c3d[frag2] = rc3d2
    geom.coords3d = c3d
    # geom.dump_xyz("reflected.xyz")

    geom_ref = geom_loader("lib:affine3_output.xyz")
    assert geom.rmsd(geom_ref) == pytest.approx(0.0, abs=1e-8)
