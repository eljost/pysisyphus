# [1] https://doi.org/10.26434/chemrxiv-2022-253hc-v2
#     Identifying the Marcus dimension of electron transfer from
#     ab initio calculations
#     Šrut, Lear, Krewald, 2022 on chemrxiv
# [2] https://doi.org/10.1039/D3SC01402A
#     Identifying the Marcus dimension of electron transfer from
#     ab initio calculations
#     Šrut, Lear, Krewald, 2023, actual published version


from distributed import LocalCluster
import matplotlib.pyplot as plt
import psutil

from pysisyphus.helpers_pure import highlight_text
from pysisyphus.marcus_dim.fit import epos_from_wf, fit_marcus_dim
from pysisyphus.marcus_dim.param import param_marcus
from pysisyphus.marcus_dim.scan import plot_scan, scan


def run_marcus_dim(
    geom, fragments, calc_getter, fit_kwargs=None, cluster=True, scan_kwargs=None
):
    """Fit Marucs dimension, scan along it and try to parametrize Marcus models."""
    if fit_kwargs is None:
        fit_kwargs = {}
    if scan_kwargs is None:
        scan_kwargs = {}
    fit_kwargs = fit_kwargs.copy()
    scan_kwargs = scan_kwargs.copy()

    # 1.) Fit Marcus dimension
    print(highlight_text("Fitting of Marcus Dimension"))

    # Use scheduler/cluster/cluster_kwargs as in ChainOfStates.py?
    if cluster:
        n_workers = psutil.cpu_count(logical=False)
        scheduler = LocalCluster(n_workers=n_workers, threads_per_worker=1)
        fit_kwargs["scheduler"] = scheduler

    marcus_dim = fit_marcus_dim(geom, calc_getter, fragments, **fit_kwargs)

    if cluster:
        scheduler.close()

    # with np.load("marcusdim_results.npz") as handle:
    # marcus_dim = handle["marcus_dim"]

    # 2.) Scan along Marcus dimension
    print()
    print(highlight_text("Scan along Marcus Dimension"))
    pos_calc = calc_getter(base_name="scan_pos")
    neg_calc = calc_getter(base_name="scan_neg")

    def get_properties(factor, coords_i):
        # Use one of two different calculators, depending on the scan direction.
        calc = pos_calc if factor > 0.0 else neg_calc

        results = calc.get_energy(geom.atoms, coords_i)
        energy = results["energy"]
        wf = calc.get_stored_wavefunction()
        tot_epos, alpha_epos = epos_from_wf(wf, fragments)
        return energy, alpha_epos

    scan_factors, scan_coords, scan_energies, scan_properties = scan(
        coords_init=geom.cart_coords.copy(),
        direction=marcus_dim,
        get_properties=get_properties,
        **scan_kwargs,
    )
    scan_fig, scan_axs = plot_scan(scan_factors, scan_energies, scan_properties)
    scan_fig.savefig("scan.png")
    plt.close(scan_fig)

    # 3.) Parametrization of Marcus model
    print()
    print(highlight_text("Parametrization of Marcus Model"))
    models = param_marcus(scan_factors, scan_energies)
    for model in models:
        print(model)
