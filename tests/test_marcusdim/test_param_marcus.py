# [1] https://doi.org/10.1039/D3SC01402A
#     Identifying the Marcus dimension of electron transfer from
#     ab initio calculations
#     Šrut, Lear, Krewald, 2023, actual published version


import math

import numpy as np

from pysisyphus.drivers.marcusdim_param import solve_marcus_wavenums_and_ang


def assert_parametrization(model, lambda_ref, dG_ref, _2Vab_ref, d_ref):
    *energies, _, _, d = model.as_wavenums_and_ang_tuple()
    np.testing.assert_allclose(
        energies, (lambda_ref, dG_ref, _2Vab_ref / 2.0), atol=0.5
    )
    assert math.isclose(d, d_ref, abs_tol=5e-4)


def test_mdnb():
    """Test parametrization of Marcus Model with data from [1].

    Ab initio scan data from Table 1 in [1] for m-DNB·⁻, second to last row.
    """
    dG_ref = 1233
    en_exc_min_ref = 7434
    Vab_ref = 2752 / 2.0
    res = solve_marcus_wavenums_and_ang(
        R=0.06, Vab=Vab_ref, dG=dG_ref, en_exc_min=en_exc_min_ref
    )
    ma = res["a"]
    mb = res["b"]

    print("Parametrization A")
    print("\t", ma.pretty())
    print("Parametrization B")
    print("\t", mb.pretty())

    assert_parametrization(ma, 9651, dG_ref, 2 * Vab_ref, 0.125)
    assert_parametrization(mb, en_exc_min_ref, 737, 2 * Vab_ref, 0.129)

    # import matplotlib.pyplot as plt
    # xs = np.linspace(-0.10, 0.2)
    # figa, _ = ma.plot_diabatic(xs, show=False)
    # figa.suptitle("Model A")
    # figb, _ = mb.plot_diabatic(xs, show=False)
    # figb.suptitle("Model B")
    # plt.show()
