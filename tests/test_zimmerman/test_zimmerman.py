import shutil

import pytest

from pysisyphus.benchmarks import Benchmark
from pysisyphus.calculators import XTB
from pysisyphus.helpers import geom_loader, align_geoms
from pysisyphus.optimizers.RFOptimizer import RFOptimizer
from pysisyphus.testing import using
from pysisyphus.tsoptimizers import RSPRFOptimizer
from pysisyphus.xyzloader import write_geoms_to_trj


def calc_getter(charge, mult):
    calc = XTB(charge=charge, mult=mult, pal=6)
    return calc


TS_FAILS = [7, 8, 13, 14, 15, 34, 39]
SET_1 = list(range(72))
SET_2 = list(range(72, 105))


def getZm(coord_type="cart"):
    Zm = Benchmark(
        "zimmerman",
        coord_type=coord_type,
        # calc_getter=calc_getter,
        # only=(0, ),
        exclude=TS_FAILS + SET_2,
    )
    return Zm


ZmTS = getZm("cart")


@pytest.mark.benchmark
@using("xtb")
@pytest.mark.parametrize("fn, geoms, charge, mult, ref_energy", ZmTS.geom_iter)
def test_zimmerman_ts_opt(fn, geoms, charge, mult, ref_energy):
    for geom in geoms:
        geom.set_calculator(calc_getter(charge, mult))
    ts = geoms[1]
    ts_ref = ts.copy()
    id_ = fn[-7:-4]

    tsopt_kwargs = {
        "thresh": "gau",
        "dump": True,
        "max_cycles": 75,
        "hessian_recalc": 1,
        "trust_max": 0.3,
    }
    opt = RSPRFOptimizer(ts, **tsopt_kwargs)
    opt.run()

    rmsd = ts.rmsd(ts_ref)
    print(f"@{id_}: converged={opt.is_converged}, rmsd={rmsd:.6f} au")
    assert opt.is_converged

    shutil.copy("final_geometry.xyz", f"{id_}_ts.xyz")


ZmGS = getZm("cart")


@pytest.mark.benchmark
@using("xtb")
@pytest.mark.parametrize("fn, geoms, charge, mult, ref_energy", ZmGS.geom_iter)
def test_zimmerman_gs_opt(fn, geoms, charge, mult, ref_energy):
    for geom in geoms:
        geom.set_calculator(calc_getter(charge, mult))
    start, _, end = geoms
    id_ = fn[-7:-4]

    for prefix, geom in (("start", start), ("end", end)):
        opt_kwargs = {
            "dump": True,
            "max_cycles": 75,
            "hessian_recalc": 10,
        }
        opt = RFOptimizer(geom, **opt_kwargs)
        opt.run()

        print(f"@{id_}: converged={opt.is_converged}")
        assert opt.is_converged

        shutil.copy("final_geometry.xyz", f"{id_}_{prefix}.xyz")


def join_gs_and_ts_geoms():
    Zm = getZm()
    for i, (fn, geom, _) in enumerate(Zm):
        id_ = fn[-7:-4]
        start_ = f"{id_}_start.xyz"
        ts_ = f"{id_}_ts.xyz"
        end_ = f"{id_}_end.xyz"
        geoms = [geom_loader(fn) for fn in (start_, ts_, end_)]
        align_geoms(geoms)

        new_id = f"{i:02d}"
        trj_ = f"{new_id}_zm_xtb.trj"
        comment = f"org_id={id_}, gfn-xtb reopt, "
        comments = [comment + suffix for suffix in ("start", "ts", "end")]
        write_geoms_to_trj(geoms, trj_, comments=comments)


# Exceute the tests
#
#   test_zimmerman_ts_opt()
#   test_zimmerman_gs_opt()
#
# before calling the function below.
#
# join_gs_and_ts_geoms()
