import argparse
import matplotlib.pyplot as plt
import numpy as np
import sys
import yaml

from pysisyphus.wavefunction.diabatization import dq_diabatization


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


"""
def diabatize_path(adia_ens, dip_moms, quad_moms_tr=None, epots=None, **kwargs):
    assert len(adia_ens) > 0
    arrs = [arr for arr in (dip_moms, quad_moms_tr, epots) if arr is not None]
    # All present arrays must be of same length
    assert [len(arr) == len(adia_ens) for arr in arrs]

    for aens, dpm, qpm, epot in zip(adia_ens, dip_moms, quad_moms_tr, epots):
        yield dq_diabatization(aens, dpm, qpm, epot, **kwargs)


def plot_adia_dia(adia_energies, dia_energies, show=True):
    fig, ax = plt.subplots()
    for i, state in enumerate(adia_energies.T):
        ax.plot(state, "o--", label=f"$V_{i}$")
    for i, state in enumerate(dia_energies.T):
        ax.plot(state, "x-", label=f"$U_{i}$")
    ax.legend()
    ax.set_xlabel("Step")
    ax.set_ylabel(r"$\Delta E$ / eV")
    if show:
        plt.show()
    return fig, ax
"""


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
