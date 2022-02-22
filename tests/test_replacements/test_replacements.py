import pytest

from pysisyphus.drivers.replace import replace_atom, replace_atoms
from pysisyphus.Geometry import Geometry
from pysisyphus.helpers import geom_loader
from pysisyphus.testing import using


@pytest.mark.parametrize("opt", [pytest.param(True, marks=using("obabel")), (False)])
def test_merge_methans(opt):
    geom = geom_loader("lib:methane.xyz")
    ind = 3  # H
    # Merge two methanes to H3C-CH3
    union = replace_atom(geom, ind, geom, ind, opt=opt)
    assert len(union.atoms) == 8


@pytest.mark.parametrize("opt", [pytest.param(True, marks=using("obabel")), (False)])
def test_sature_methane(opt):
    geom = geom_loader("lib:methane.xyz")
    inds = (1, 2, 3, 4)
    repl_inds = (1, 2, 3, 4)
    frags = [geom] * len(inds)
    replacements = zip(inds, repl_inds, frags)
    union = replace_atoms(geom, replacements, opt=opt)
    assert len(union.atoms) == 17


def test_replace_atom_with_atom():
    geom = geom_loader("lib:methane.xyz")
    repl_geom = Geometry(("I",), (0.0, 0.0, 0.0))
    union = replace_atom(geom, 3, repl_geom, 0)
    assert union.atoms[3] == "I"
