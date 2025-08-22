import matplotlib.pyplot as plt
import pytest

import pysisyphus.numint.radint as radint


@pytest.mark.skip
def test_euler_maclaurin():
    fig, ax = plt.subplots()
    m_r = 2
    nn = 20
    for n in range(nn, 5 * nn, nn):
        rw = radint.euler_maclaurin_dma(n, m_r=m_r, atomic_radius=1.0)
        r, w = rw.T
        ax.plot(r, w, "o-", label=f"$n={{{n}}}$")
    ax.legend()
    ax.set_xscale("log")
    ax.set_xlabel("log(r)")
    ax.set_yscale("log")
    ax.set_ylabel("log(w)")
    plt.show()
