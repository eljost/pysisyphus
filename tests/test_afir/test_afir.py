#!/usr/bin/env python3

from pysisyphus.calculators.AFIR import AFIR
from pysisyphus.calculators.XTB import XTB
from pysisyphus.helpers import geom_from_library
from pysisyphus.optimizers.RFOptimizer import RFOptimizer


def test_afir():
    # geom = geom_from_library("ohch3f_anion_cs.xyz", coord_type="redund")
    geom = geom_from_library("ohch3f_anion_cs.xyz")
    fragment_indices = ([0, 1, 2, 3, 4], [5, 6])
    calc = XTB(charge=-1)
    # gamma = 100 * 3.8088e-4
    gamma = 100
    afir = AFIR(calc, geom.atoms, fragment_indices, gamma)
    geom.set_calculator(afir)

    opt = RFOptimizer(geom, dump=True, trust_max=.3)
    opt.run()
    assert opt.is_converged


def plot_afir():
    from pysisyphus.constants import AU2KJPERMOL
    from pysisyphus.peakdetect import peakdetect
    import matplotlib.pyplot as plt
    import numpy as np
    import yaml

    with open("image_results.yaml") as handle:
        res = yaml.load(handle.read(), Loader=yaml.loader.Loader)

    afir_ens = [_["energy"] for _ in res]
    true_ens = [_["true_energy"] for _ in res]
    afir_ens = np.array(afir_ens) * AU2KJPERMOL
    afir_ens -= afir_ens.min()
    true_ens = np.array(true_ens) * AU2KJPERMOL
    true_ens -= true_ens.min()

    afir_forces = np.linalg.norm([_["forces"] for _ in res], axis=1)
    true_forces = np.linalg.norm([_["true_forces"] for _ in res], axis=1)
    afir_forces = np.array(afir_forces)
    true_forces = np.array(true_forces)

    peak_inds, _ = peakdetect(true_ens, lookahead=2)
    peak_xs, peak_ys = zip(*peak_inds)
    highest = np.argmax(peak_ys)

    fig, (en_ax, forces_ax) = plt.subplots(nrows=2, sharex=True)

    style1 = "ro-"
    style2 = "go-"

    l1 = en_ax.plot(afir_ens, style1, label="AFIR")
    en_ax2 = en_ax.twinx()
    l2 = en_ax2.plot(true_ens, style2, label="True")
    en_ax2.scatter(peak_xs, peak_ys, s=100, marker="X", c="k", zorder=10)
    en_ax2.scatter(peak_xs[highest], peak_ys[highest],
                s=150, marker="X", c="k", zorder=10)
    en_ax.axvline(peak_xs[highest], c="k", ls="--")

    # lines = l1 + l2
    # labels = [l.get_label() for l in lines]
    # en_ax.legend(lines, labels, loc=0)

    en_ax.set_title("Energies")
    en_ax.set_ylabel("$\Delta$E kJ / mol")

    forces_ax.set_title("||Forces||")
    forces_ax.plot(afir_forces, style1)
    forces_ax2 = forces_ax.twinx()
    forces_ax2.plot(true_forces, style2)
    forces_ax.axvline(peak_xs[highest], c="k", ls="--")

    fig.legend(loc="upper center")
    plt.show()


if __name__ == "__main__":
    test_afir()
    plot_afir()
