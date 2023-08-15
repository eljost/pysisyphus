from concurrent.futures import ThreadPoolExecutor
from os import kill
from signal import SIGTERM
from subprocess import check_output
from time import sleep

import numpy as np
import pytest

from pysisyphus.helpers import geom_loader
from pysisyphus.init_logging import init_logging
from pysisyphus.calculators import ORCA
from pysisyphus.calculators.ORCA import (
    parse_orca_cis,
    parse_orca_densities,
    parse_orca_gbw_new,
    set_mo_coeffs_in_gbw,
)
from pysisyphus.config import WF_LIB_DIR
from pysisyphus.testing import using
from pysisyphus.wavefunction import Wavefunction


@pytest.fixture
def retry_geom():
    geom = geom_loader("lib:h2o.xyz")
    calc_kwargs = {
        "keywords": "hf sto-3g",
        "retry_calc": 1,
    }
    calc = ORCA(**calc_kwargs)
    geom.set_calculator(calc)
    return geom


@using("orca")
def test_orca_check_termination(retry_geom, this_dir):
    calc = retry_geom.calculator
    normal_term_out = this_dir / "normal_term.orca.out"

    assert calc.check_termination(normal_term_out)

    with open(normal_term_out) as handle:
        text = handle.read()
    assert calc.check_termination(text)

    fail_term_out = this_dir / "fail_term.orca.out"
    assert not calc.check_termination(fail_term_out)


@using("orca")
def test_print_capabilities(retry_geom):
    retry_geom.calculator.print_capabilities()


# Disabled for now as it uses 'pidof' which may not be available in Nix
@pytest.mark.skip
@using("orca")
def test_orca_restart():
    """Spawn two threads. One executes ORCA with retry_calc=1, the other thread
    kills the ORCA process. Then retry_calc kicks in and ORCA is executed again."""
    init_logging()

    geom = geom_loader("lib:benzene.xyz")
    calc = ORCA(
        keywords="b3lyp def2-svp tightscf",
        pal=1,
        retry_calc=1,
    )
    geom.set_calculator(calc)

    def run_orca():
        geom.forces  # Execute ORCA by requesting an engrad calculation
        return "forces returned"

    def kill_orca():
        secs = 2
        sleep(secs)
        print(f"Waited {secs} s ")
        pid = check_output(["pidof", "orca"], text=True)
        pid = int(pid)
        kill(pid, SIGTERM)
        print("Killed orca")
        return "kill returned"

    with ThreadPoolExecutor() as executor:
        running_tasks = [executor.submit(task) for task in (run_orca, kill_orca)]
        for running_task in running_tasks:
            print(running_task.result())

    # Run ORCA here and kill it manually with `pkill orca`. This seems quite
    # reliable.
    # run_orca()


@using("orca")
def test_orca_parse_triplet_energies(this_dir):
    out_fn = this_dir / "orca_tddft_triplets.out"
    calc = ORCA("")
    calc.root = 1
    calc.do_tddft = True
    calc.out = out_fn
    ens = calc.parse_all_energies(triplets=True)
    ref_ens = (-394.21337084, -394.086255836, -394.05756284)
    np.testing.assert_allclose(ens, ref_ens)


@pytest.mark.parametrize(
    "method",
    ("rhf", "uhf"),
)
@pytest.mark.parametrize(
    "tda",
    (True, False),
)
@pytest.mark.parametrize(
    "triplets",
    (True, False),
)
def test_parse_orca_cis(method, tda, triplets, this_dir):
    if (method == "uhf") and triplets:
        return
    base_name = "_".join(
        [method] + (["tda"] if tda else []) + (["triplets"] if triplets else [])
    )
    # geom = geom_loader("h2dimer.xyz")
    # nroots = 2
    # calc = ORCA(
    # base_name=base_name,
    # keywords=f"{method} hf sto-3g",
    # blocks=f"%tddft tda {tda} nroots {nroots} triplets {triplets} end",
    # charge=0,
    # mult=1,
    # pal=2,
    # out_dir="parse_cis",
    # )
    # geom.set_calculator(calc)
    # Create before the calculation is carried out the have the correct index
    # geom.energy  # Calc

    cis_fn = this_dir / "ref_cis" / f"{base_name}_000.000.orca.cis"
    Xa, Ya, Xb, Yb = parse_orca_cis(cis_fn)

    ac = Xa**2 - Ya**2
    a_sum = np.sum(ac, axis=(1, 2))
    bc = Xb**2 - Yb**2
    b_sum = np.sum(bc, axis=(1, 2))

    Y0 = np.zeros_like(Ya)
    if method == "rhf":
        np.testing.assert_allclose(a_sum, np.ones_like(a_sum))
        np.testing.assert_allclose(b_sum, np.zeros_like(b_sum))
    elif method == "uhf":
        half = np.ones_like(a_sum) * 0.5
        np.testing.assert_allclose(a_sum, half)
        np.testing.assert_allclose(b_sum, half)
    else:
        raise Exception(f"Unknown method='{method}'")

    # Y-matrix is zero for TDA
    if tda:
        np.testing.assert_allclose(Ya, Y0)
        np.testing.assert_allclose(Yb, Y0)

    # β-spin matrices will be 0, as the excitation is from a closed-shell
    # reference.
    if triplets:
        np.testing.assert_allclose(Xb, Y0)
        np.testing.assert_allclose(Yb, Y0)

    print(f"{method=}, {tda=}, {triplets=}")

    def print_summary(X, Y, aorb):
        print(f"X{aorb}")
        print(X)
        print(f"Y{aorb}")
        print(Y)
        norms = X**2 - Y**2
        norms_sum = np.sum(norms, axis=(1, 2))
        print(f"@{method}, {tda=}, {aorb} coeffs, {norms_sum}")
        print(norms)

    print_summary(Xa, Ya, "a")
    print_summary(Xb, Yb, "b")
    print("@")


def assert_dens_mats(dens_dict, json_fn):
    wf = Wavefunction.from_orca_json(json_fn)
    Pa_ref, Pb_ref = wf.P  # Always unrestricted density matrices
    P = Pa_ref + Pb_ref  # Electronic density
    Pr = Pa_ref - Pb_ref  # Spin density
    np.testing.assert_allclose(P, dens_dict["scfp"], atol=1e-14)
    if wf.unrestricted:
        np.testing.assert_allclose(Pr, dens_dict["scfr"], atol=1e-14)
    return wf


@pytest.mark.parametrize(
    "dens_fn, json_fn",
    (
        ("orca_ch4_sto3g_rhf.densities", "orca_ch4_sto3g_rhf.json"),
        ("orca_ch4_sto3g_uhf.densities", "orca_ch4_sto3g_uhf.json"),
    ),
)
def test_orca_gs_densities(dens_fn, json_fn):
    dens_dict = parse_orca_densities(WF_LIB_DIR / dens_fn)
    _ = assert_dens_mats(dens_dict, WF_LIB_DIR / json_fn)


@pytest.mark.parametrize(
    "dens_fn, json_fn, ref_dpm",
    (
        (
            "orca_ch4_sto3g_rhf_cis.densities",
            "orca_ch4_sto3g_rhf_cis.json",
            (0.00613, 0.00867, -0.00000),
        ),
        (
            "orca_ch4_sto3g_uhf_cis.densities",
            "orca_ch4_sto3g_uhf_cis.json",
            (0.0, 0.0, 0.0),
        ),
    ),
)
def test_orca_es_densities(dens_fn, json_fn, ref_dpm):
    dens_dict = parse_orca_densities(WF_LIB_DIR / dens_fn)
    wf = assert_dens_mats(dens_dict, WF_LIB_DIR / json_fn)
    cisp = dens_dict["cisp"]
    dpm = wf.get_dipole_moment(cisp)
    np.testing.assert_allclose(dpm, ref_dpm, atol=2e-4)


@using("orca")
def test_orca_hf():
    geom = geom_loader("lib:h2o.xyz")
    calc = ORCA("hf sto-3g")
    geom.set_calculator(calc)
    energy = geom.energy
    assert energy == pytest.approx(-74.960702484)


@using("orca")
def test_orca_stored_wavefunction():
    geom = geom_loader("lib:h2o.xyz")
    calc = ORCA("hf sto-3g", wavefunction_dump=True)
    geom.set_calculator(calc)
    geom.energy
    # Wavefunction already does some internal sanity checking
    calc.get_stored_wavefunction()


def test_set_gbw_restricted(this_dir, tmp_path):
    """Do a roundtrip."""
    gbw_in = this_dir / "restricted.gbw"
    gbw_out = tmp_path / "restricted.gbw"
    moc = parse_orca_gbw_new(gbw_in)
    set_mo_coeffs_in_gbw(gbw_in, gbw_out, moc.Ca)
    moc2 = parse_orca_gbw_new(gbw_out)
    np.testing.assert_allclose(moc2.Ca, moc.Ca)


def test_set_gbw_unrestricted(this_dir, tmp_path):
    """Do a roundtrip."""
    gbw_in = this_dir / "unrestricted.gbw"
    gbw_out = tmp_path / "unrestricted.gbw"
    moc = parse_orca_gbw_new(gbw_in)
    set_mo_coeffs_in_gbw(gbw_in, gbw_out, moc.Ca, moc.Cb)
    moc2 = parse_orca_gbw_new(gbw_out)
    np.testing.assert_allclose(moc2.Ca, moc.Ca)
    np.testing.assert_allclose(moc2.Cb, moc.Cb)


def test_set_gbw_empty(this_dir, tmp_path):
    """Do a roundtrip."""
    gbw_in = this_dir / "unrestricted.gbw"
    gbw_out = tmp_path / "unrestricted.gbw"
    moc = parse_orca_gbw_new(gbw_in)
    set_mo_coeffs_in_gbw(gbw_in, gbw_out)
    moc2 = parse_orca_gbw_new(gbw_out)
    np.testing.assert_allclose(moc2.Ca, moc.Ca)
    np.testing.assert_allclose(moc2.Cb, moc.Cb)
