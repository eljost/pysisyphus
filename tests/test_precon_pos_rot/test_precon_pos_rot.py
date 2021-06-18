from pathlib import Path

import pytest

from pysisyphus.benchmarks import Benchmark
from pysisyphus.drivers.precon_pos_rot import precon_pos_rot
from pysisyphus.helpers import geom_loader
# from pysisyphus.helpers_pure import filter_fixture_store


Bm = Benchmark("precon_pos_rot")
# Bm = Benchmark("precon_pos_rot", only=0)


def test_precon_pos_rot_figure2(this_dir):
    educt, product = geom_loader(this_dir / "figure2_mod.trj")
    rgeom, pgeom = precon_pos_rot(educt, product, prefix="figure2_")


@pytest.mark.parametrize("fn, geoms, ref_energy", Bm)
def test_birkholz_benchmark(fn, geoms, ref_energy, results_bag):
    prefix = Path(fn).stem + "_"
    reactants, _, products = geoms
    rcp = reactants.copy()
    pcp = products.copy()
    rgeom, pgeom = precon_pos_rot(reactants, products, prefix=prefix)
    rr = rgeom.rmsd(rcp)
    pr = pgeom.rmsd(pcp)
    print(f"@@@ {fn} R_RMSD={rr:.4f} P_PRMSD={pr:.4f}")

    r_org_xyz = rcp.as_xyz(comment="R, original")
    p_org_xyz = pcp.as_xyz(comment="P, original")
    r_xyz = rgeom.as_xyz(comment="R, precon")
    p_xyz = pgeom.as_xyz(comment="P, precon")
    trj = "\n".join((r_org_xyz, r_xyz, p_org_xyz, p_xyz))
    with open(prefix + "comp.trj", "w") as handle:
        handle.write(trj)

    results_bag.fn = fn
    results_bag.r_rmsd = rr
    results_bag.p_rmsd = pr
    results_bag.trj = trj


# @filter_fixture_store("test_birkholz_benchmark")
def test_birkholz_benchmark_synthesis(fixture_store):
    for i, fix in enumerate(fixture_store):
        print(i, fix)

    bags = fixture_store["results_bag"]
    trjs = list()
    for k, v in bags.items():
        if not k.startswith("test_birkholz_benchmark"):
            continue
        print(v.fn)
        print(f"\tRMSD(R)={v.r_rmsd:.4f}")
        print(f"\tRMSD(P)={v.p_rmsd:.4f}")
        trjs.append(v.trj)
    with open("birkholz.trj", "w") as handle:
        handle.write("\n".join(trjs))
