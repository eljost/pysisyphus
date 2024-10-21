import pytest

from pysisyphus.constants import BOHR2ANG
from pysisyphus.helpers import geom_loader
from pysisyphus.partfuncs import inertmom


# Values from Table VI on p. 6660 of https://doi.org/10.1063/1.473958
REF = {
    (1, 1): (3.1814, 0.8640),
    (1, 2): (0.6794, 0.6794),
    (1, 3): (0.6406, 0.6766),
    #
    (2, 1): (3.1810, 0.8128),
    (2, 2): (0.6474, 0.6474),
    (2, 3): (0.6408, 0.6469),
    #
    (3, 1): (3.1862, 0.7915),
    (3, 2): (0.6340, 0.6340),
    (3, 3): (0.6321, 0.6339),
}


@pytest.mark.parametrize("n", (1, 2, 3))
@pytest.mark.parametrize("m", (1, 2, 3))
def test_inertias(n, m, this_dir):
    # Geometry was obtained w/ ORCA 5.0.4
    geom = geom_loader(this_dir / "meoh_mp2_631gd_opt.xyz")
    indO = 0
    indC = 1
    # C H H H
    top = [1, 2, 3, 4]
    bond_inds = [indO, indC]

    fmt = "> 8.4f"
    print(f"{n=}, {m=}")
    Itop, Irest = inertmom.get_top_moment_of_inertia(
        geom.coords3d, geom.masses, top, bond_inds, m=m, n=n
    )
    Itop *= BOHR2ANG**2
    Irest *= BOHR2ANG**2
    print(f"{top}: ({Itop=:{fmt}}, {Irest:{fmt}}) amu Å²")
    Itop_ref, Irest_ref = REF[(n, m)]
    assert Itop == pytest.approx(Itop_ref, abs=6e-5)
    assert Irest == pytest.approx(Irest_ref, abs=6e-5)
