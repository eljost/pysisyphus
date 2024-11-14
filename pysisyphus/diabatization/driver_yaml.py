import argparse
import itertools as it
import sys

import matplotlib.pyplot as plt
import numpy as np
import yaml

from pysisyphus.constants import AU2EV
from pysisyphus.diabatization.multipole import dq_diabatization


def make_array(nstates, components, lines):
    if len(lines) == 0:
        return None

    arr = np.zeros((components, nstates, nstates))
    expect_diag = nstates
    expect_off_diag = sum(range(nstates))

    for line in lines:
        from_, to_, *props = line
        assert (from_ >= 0) and (to_ >= 0)
        assert (
            len(props) == components
        ), f"Expected line of length {components} but got '{line}'!"

        if from_ == to_:
            expect_diag -= 1
        elif from_ != to_:
            expect_off_diag -= 1

        arr[:, from_, to_] = arr[:, to_, from_] = props
    assert expect_diag == 0, "Some diagonal elements are missing!"
    assert expect_off_diag == 0, "Some off-diagonal elements are missing!"
    return arr


def dq_diabatization_from_run_dict(run_dict):
    adia_ens = np.array(run_dict["adiabatic_energies"], dtype=float)
    nstates = adia_ens.size

    # Dipole moments must be present
    dip_moms = make_array(nstates, 3, run_dict["dipoles"])
    # Quadrupole moments and electrostatic potential are optional.
    quad_moms = make_array(nstates, 1, run_dict.get("quadrupoles", list()))
    epots = make_array(nstates, 1, run_dict.get("epots", list()))

    kwargs = {}
    if "alpha" in run_dict:
        kwargs["alpha"] = run_dict["alpha"]
    if "beta" in run_dict:
        kwargs["beta"] = run_dict["beta"]

    return dq_diabatization(
        adia_ens, dip_moms, quad_moms=quad_moms, epots=epots, **kwargs
    )


def diabatize_path(adia_ens, dip_moms, tr_quad_moms=None, epots=None, **kwargs):
    nones = [None for _ in adia_ens]
    if tr_quad_moms is None:
        tr_quad_moms = nones
    if epots is None:
        epots = nones
    assert len(adia_ens) == len(dip_moms) == len(tr_quad_moms) == len(epots)

    for aens, dpm, qpm, epot in zip(adia_ens, dip_moms, tr_quad_moms, epots):
        yield dq_diabatization(aens, dpm, qpm, epot, **kwargs)


def plot_dia_res(dia_res, show=False):
    nstates = dia_res[0].nstates
    adia_ens = np.zeros((len(dia_res), nstates))
    dia_ens = np.zeros((len(dia_res), nstates))
    keys = list(it.combinations(range(nstates), 2))
    couplings = np.zeros((len(dia_res), len(keys)))
    for i, dr in enumerate(dia_res):
        adia_ens[i] = dr.adia_ens
        dia_ens[i] = dr.dia_ens
        for j, key in enumerate(keys):
            couplings[i, j] = dr.couplings[key]

    adia_min = adia_ens.min()
    adia_ens -= adia_min
    adia_ens *= AU2EV
    dia_ens -= adia_min
    dia_ens *= AU2EV
    couplings *= AU2EV

    fig, (ax0, ax1) = plt.subplots(nrows=2, sharex=True)
    for i, state in enumerate(adia_ens.T):
        ax0.plot(state, "o--", label=f"$V_{i}$")
    for i, state in enumerate(dia_ens.T):
        ax0.plot(state, "x-", label=f"$U_{i}$")
    ax0.legend()
    ax0.set_xlabel("Step")
    ax0.set_ylabel(r"$\Delta E$ / eV")

    # Couplings
    for i, cpls in enumerate(couplings.T):
        from_to = "".join([str(_) for _ in keys[i]])
        ax1.plot(cpls, "o-", label=f"$|D_{{{from_to}}}|$")
    ax1.axhline(0.0, ls="--", c="k")
    ax1.legend()

    fig.tight_layout()
    if show:
        plt.show()
    return fig, (ax0, ax1)


def parse_args(args):
    parser = argparse.ArgumentParser()
    parser.add_argument("yaml")
    return parser.parse_args(args)


def run():
    args = parse_args(sys.argv[1:])

    yaml_fn = args.yaml
    with open(yaml_fn) as handle:
        run_dict = yaml.load(handle, Loader=yaml.SafeLoader)
    adia_labels = run_dict.pop("adiabatic_labels", None)
    unit = run_dict.pop("unit", "eV")
    dia_res = dq_diabatization_from_run_dict(run_dict)
    report = dia_res.render_report(adia_labels=adia_labels, unit=unit)
    print(report)


if __name__ == "__main__":
    run()
