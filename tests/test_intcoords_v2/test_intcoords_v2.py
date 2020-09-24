import numpy as np
import pytest

from pysisyphus.calculators import XTB
from pysisyphus.calculators.PySCF import PySCF
from pysisyphus.helpers import geom_loader
from pysisyphus.intcoords import Bend, LinearBend, Stretch, Bend, Torsion
from pysisyphus.optimizers.RFOptimizer import RFOptimizer
from pysisyphus.testing import using


def compare_internals(xyz_fn):
    geom = geom_loader(xyz_fn, coord_type="redund")

    cls = {
        2: Stretch,
        3: Bend,
        4: Torsion,
    }
    for pi in geom.internal._prim_internals:
        inds = pi.inds
        obj = cls[len(inds)](inds)
        print(obj)
        v, g = obj.calculate(geom.coords3d, gradient=True)
        assert v == pytest.approx(pi.val)
        np.testing.assert_allclose(g, pi.grad)


def test_compare_h2o2():
    xyz_fn = "lib:h2o2_hf_321g_opt.xyz"
    compare_internals(xyz_fn)


def test_compare_biaryl():
    xyz_fn = "lib:biaryl_bare_pm6_splined_hei.xyz"
    compare_internals(xyz_fn)


def test_bmat_allene():
    xyz_fn = "lib:08_allene.xyz"
    ref_geom = geom_loader(xyz_fn, coord_type="redund")
    print(ref_geom)
    ref_int = ref_geom.internal
    Bref = ref_int.B

    geom = geom_loader(xyz_fn, coord_type="redund",
                       coord_kwargs={
                           "lb_min_deg": None,
                        }
    )
    print(geom)
    int_ = geom.internal
    B = int_.B
    
    np.testing.assert_allclose(B, Bref)


# @pytest.mark.skip
@using("pyscf")
def test_allene_opt():
    geom = geom_loader("lib:08_allene.xyz", coord_type="redund")

    calc = PySCF(basis="321g", pal=1)
    geom.set_calculator(calc)
    opt = RFOptimizer(geom, thresh="gau_tight")
    opt.run()

    assert opt.is_converged
    assert opt.cur_cycle == 5


def test_hydrogen_bonds_fragments():
    geom = geom_loader("lib:hydrogen_bond_fragments_test.xyz",
                       coord_type="redund")
    int_ = geom.internal
    assert len(int_.hydrogen_bond_indices) == 1
    assert len(int_.fragments) == 3


def test_linear_bend_allene():
    np.set_printoptions(precision=5)
    xyz_fn = "lib:08_allene.xyz"
    geom = geom_loader(xyz_fn, coord_type="redund")
    c3d = geom.coords3d

    from pysisyphus.intcoords.LinearBend import LinearBend
    inds = (0, 3, 4)# 1)
    lb = LinearBend(inds)
    _ = lb.calculate(c3d)
    __ = np.rad2deg(_)
    print(_, __)
    g, grad = lb.calculate(c3d, gradient=True)
    grad = grad.reshape(-1, 3)
    print(grad)
    print("huhu")
    from pysisyphus.constants import BOHR2ANG
    g16g = grad[list(inds)].reshape(-1, 1) #/ BOHR2ANG
    print(g16g)

    co = LinearBend(inds, complement=True)
    cv, cg = co.calculate(c3d, gradient=True)
    print("complement")
    cvd = np.rad2deg(cv)
    print(f"\t{cv} {cvd}")
    c16g = cg.reshape(-1, 3)[list(inds)].reshape(-1, 1)
    print(c16g)
    print(c16g.flatten().dot(g16g.flatten()))
    print()
    print(np.concatenate((g16g, c16g), axis=1))


def test_linear_bend_c4():
    np.set_printoptions(precision=5)
    # xyz_fn = "lib:c4_lb_test.xyz"
    xyz_fn = "lib:c4_lb_test_02.xyz"
    geom = geom_loader(xyz_fn, coord_type="redund")
    c3d = geom.coords3d

    from pysisyphus.intcoords.LinearBend import LinearBend
    inds = (0, 1, 2)
    lb = LinearBend(inds)
    _ = lb.calculate(c3d)
    __ = np.rad2deg(_)
    print("val=", _, "grad=", __)
    g, grad = lb.calculate(c3d, gradient=True)
    grad = grad.reshape(-1, 3)
    print(grad)
    print("huhu")
    from pysisyphus.constants import BOHR2ANG
    g16g = grad[list(inds)].reshape(-1, 1) #/ BOHR2ANG
    print(g16g)

    # import pdb; pdb.set_trace()
    co = LinearBend(inds, complement=True)
    cv, cg = co.calculate(c3d, gradient=True)
    print("complement")
    cvd = np.rad2deg(cv)
    print(f"\t{cv} {cvd}")
    c16g = cg.reshape(-1, 3)[list(inds)].reshape(-1, 1)
    print(c16g)
    print(c16g.flatten().dot(g16g.flatten()))

    print(np.concatenate((g16g, c16g), axis=1))


def test_molcas_lb():
    np.set_printoptions(precision=5)
    xyz_fn = "lib:08_allene.xyz"
    geom = geom_loader(xyz_fn, coord_type="redund")
    # c3d = geom.coords3d


@using("xtb")
@pytest.mark.parametrize(
    "make_complement", [
        True,
        False
    ]
)
@pytest.mark.parametrize(
    "fn", [
        "lib:co2_bent.xyz",
        "lib:h2o.xyz",
    ]
)
def test_complement(make_complement, fn):
    coord_kwargs = {
        "make_complement": make_complement,
    }
    geom = geom_loader(fn, coord_type="redund",
                       coord_kwargs=coord_kwargs,
    )
    int_ = geom.internal
    prims = int_._primitives
    # List primitives
    for pi in prims:
        print(pi)
    print()

    if make_complement and (len(prims) == 5):
        lb = prims[3]
        comp = prims[4]
        c3d = geom.coords3d
        # Gradients
        _, lbg = lb.calculate(c3d, gradient=True)
        lbg = lbg.reshape(-1, 3)
        _, cg = comp.calculate(c3d, gradient=True)
        cg = cg.reshape(-1, 3)

        for pi in (lb, comp):
            q, grad = pi.calculate(c3d, gradient=True)
            print(pi)
            print(q)
            print(grad)
    geom.set_calculator(XTB())

    opt_kwargs = {
        "thresh": "gau_loose",
        "dump": True,
        "line_search": False,
        "gdiis": False,
        # "max_cycles": 15,
        "max_cycles": 11,
    }
    opt = RFOptimizer(geom, **opt_kwargs)
    opt.run()

    print()


def test_derivs():
    coord_kwargs = {
        "make_complement": True,
    }
    fn = "lib:co2_bent.xyz"
    # fn = "lib:co2.xyz"
    geom = geom_loader(fn, coord_type="redund",
                       coord_kwargs=coord_kwargs,
    )
    int_ = geom.internal
    prims = int_._primitives


    c3d = geom.coords3d
    indices = (1, 0, 2)
    slb = LinearBend(indices)
    lb, glb = slb.calculate(c3d, gradient=True)
    clb = LinearBend(indices, complement=True)
    cb, gcb = clb.calculate(c3d, gradient=True)

    return
    # List primitives
    for pi in prims:
        print(pi)
    print()

    if len(prims) == 5:
        lb = prims[3]
        comp = prims[4]
        c3d = geom.coords3d
        # Gradients
        _, lbg = lb.calculate(c3d, gradient=True)
        lbg = lbg.reshape(-1, 3)
        _, cg = comp.calculate(c3d, gradient=True)
        cg = cg.reshape(-1, 3)

        for pi in (lb, comp):
            q, grad = pi.calculate(c3d, gradient=True)
            print(pi)
            print(q)
            print(grad)


def test_lb():
    from pysisyphus.Geometry import Geometry

    coords = np.array((
        ( 0.,  0., 0.),
        ( 0.,  0., 2.),
        ( 0., 0., -1.),
    ))
    atoms = ("O", "O", "O")
    geom = Geometry(atoms, coords, coord_type="redund")
    indices = (1, 0, 2)
    lb = LinearBend(indices)
    print(lb)
    q, g = lb.calculate(coords, gradient=True)
    g = g.reshape(-1, 3)

    print("val")
    print(q)
    print("grad")
    print(g)


    print()
    print("COMPLEMENT")
    cb = LinearBend(indices, complement=True)
    print(cb)
    qc, gc = cb.calculate(coords, gradient=True)
    gc = gc.reshape(-1, 3)

    print("complement, val")
    print(qc)
    print("complement, grad")
    print(gc)
