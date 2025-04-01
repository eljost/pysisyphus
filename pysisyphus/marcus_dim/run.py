# [1] https://doi.org/10.26434/chemrxiv-2022-253hc-v2
#     Identifying the Marcus dimension of electron transfer from
#     ab initio calculations
#     Šrut, Lear, Krewald, 2022 on chemrxiv
# [2] https://doi.org/10.1039/D3SC01402A
#     Identifying the Marcus dimension of electron transfer from
#     ab initio calculations
#     Šrut, Lear, Krewald, 2023, actual published version


from collections.abc import Callable
from pathlib import Path
import os
from typing import Optional, Union
import warnings

import distributed
import matplotlib.pyplot as plt
import numpy as np


from pysisyphus.executors import MaybeScheduler
from pysisyphus.exceptions import (
    CalculationFailedException,
    RunAfterCalculationFailedException,
)
from pysisyphus.Geometry import Geometry
from pysisyphus.helpers_pure import get_state_label, highlight_text
from pysisyphus.marcus_dim.config import DUMMY_2EV, FIT_RESULTS_FN, SCAN_RESULTS_FN
from pysisyphus.marcus_dim.fit import epos_from_wf, fit_marcus_dim
from pysisyphus.marcus_dim.param import param_marcus
from pysisyphus.marcus_dim.scan import plot_scan, scan, scan_parallel
import pysisyphus.marcus_dim.types as mdtypes


def dump_scan_result(geom, scan_result, dummy_scan, cwd="."):
    cwd = Path(cwd)
    # Dump scan geometries into trj file
    xyzs = list()
    for sc, (se_gs, *se_es) in zip(scan_result.coords, scan_result.energies):
        xyz = geom.as_xyz(cart_coords=sc, comment=f"{se_gs:.6f}")
        xyzs.append(xyz)
    scan_trj_path = cwd / "marcus_dim_scan.trj"
    with open(scan_trj_path, "w") as handle:
        handle.write("\n".join(xyzs))

    scan_fig, scan_axs = plot_scan(
        scan_result.factors,
        scan_result.energies,
        scan_result.properties,
        dummy_scan=dummy_scan,
    )
    scan_fig.savefig(cwd / "scan.svg")
    plt.close(scan_fig)


def run_scan(
    geom: Geometry,
    marcus_dim: np.ndarray,
    calc_getter: Callable,
    fragments: list[list[int]],
    dummy_scan: bool = False,
    cwd: Union[os.PathLike, str] = ".",
    **kwargs,
) -> mdtypes.MarcusDimScanResult:
    """Scan along Marcus dimension using serial calculations in two directions.

    Parameters
    ----------
    geom
        Equilibrium geometry with N atoms and 3N Cartesian coordinates.
    marucs_dim
        1d array of shape (3N, ) containing the Marcus dimension in
        unweighted Cartesian coordinates.
    calc_getter
        Callable that returns a new calculator. Takes optional kwargs,
        that are passed to the underlying Calculator class.
    fragments
        List of lists of integers, containing the different fragment indices.
        Order is important, as the prefactors for weighting the alpha-electron
        spin density depend on the order.
    dummy_scan
        Optional; when True the scan can be carried out using a calculator that
        supports only GS calculations, e.g. xTB. In such a case, a dummy value
        for the first ES is used (GS + 2 eV).
    cwd
        String or Path object that controls the directory, where various files
        are written.

    Returns
    -------
    scan_result
        MarcusDimScanResult, carrying the factors, Cartesian coordinates,
        calculated energies and properties along the scan.
    """
    cwd = Path(cwd)
    pos_calc = calc_getter(base_name="scan_pos")
    neg_calc = calc_getter(base_name="scan_neg")

    has_all_energies = hasattr(pos_calc, "get_all_energies")
    if dummy_scan and not has_all_energies:
        warnings.warn(
            f"Calculator '{pos_calc}' does not implement 'get_all_energies()'! "
            f"Using dummy values for first excited state!'",
        )

    def get_properties(factor: int, coords_i: np.ndarray):
        # Use one of two different calculators, depending on the scan direction.
        calc = neg_calc if factor < 0 else pos_calc

        try:
            results = calc.get_all_energies(geom.atoms, coords_i)
            all_energies = results["all_energies"]
        except AttributeError as err:
            if not dummy_scan:
                raise err
            results = calc.get_energy(geom.atoms, coords_i)
            energy = results["energy"]
            all_energies = np.array((energy, energy + DUMMY_2EV))
        energies = all_energies[:2]
        assert len(energies) == 2
        wf = calc.get_stored_wavefunction()
        tot_epos, alpha_epos = epos_from_wf(wf, fragments)
        return energies, alpha_epos

    scan_result = scan(
        coords_init=geom.cart_coords.copy(),
        direction=marcus_dim,
        get_properties=get_properties,
        **kwargs,
    )
    dump_scan_result(geom, scan_result, dummy_scan, cwd=cwd)
    return scan_result


def run_scan_parallel(
    geom: Geometry,
    marcus_dim_data: dict,
    calc_getter: Callable,
    fragments: list[list[int]],
    rd_class: mdtypes.RobinDay,
    dummy_scan: bool = False,
    cwd: Union[os.PathLike, str] = ".",
    scheduler: distributed.LocalCluster = None,
    **scan_kwargs,
) -> mdtypes.MarcusDimScanResult:
    """Scan along Marcus dimension using parallel calculations.

    Parameters
    ----------
    geom
        Equilibrium geometry with N atoms and 3N Cartesian coordinates.
    marucs_dim_data
        Dictionary, as returned by fit_marcus_dim.
    calc_getter
        Callable that returns a new calculator. Takes optional kwargs,
        that are passed to the underlying Calculator class.
    fragments
        List of lists of integers, containing the different fragment indices.
        Order is important, as the prefactors for weighting the alpha-electron
        spin density depend on the order.
    rd_class
        RobinDay enum, denoting the Robin-Day class.
    dummy_scan
        Optional; when True the scan can be carried out using a calculator that
        supports only GS calculations, e.g. xTB. In such a case, a dummy value
        for the first ES is used (GS + 2 eV).
    cwd
        String or Path object that controls the directory, where various files
        are written.
    scheduler
        Dask.distributed LocalCluster object.

    Returns
    -------
    scan_result
        MarcusDimScanResult, carrying the factors, Cartesian coordinates,
        calculated energies and properties along the scan.
    """
    cwd = Path(cwd)
    # 2.) Scan along Marcus dimension
    print(highlight_text("Scan along Marcus Dimension"))

    test_calc = calc_getter()
    has_all_energies = hasattr(calc_getter(), "get_all_energies")
    if dummy_scan and not has_all_energies:
        warnings.warn(
            f"Calculator '{test_calc}' does not implement 'get_all_energies()'! "
            f"Using dummy values for first excited state!'",
        )

    atoms = geom.atoms

    def get_properties(calc_number, coords):
        print(f"Started scan calculation {calc_number}")
        calc = calc_getter(
            with_td=True, pal=1, calc_number=calc_number, base_name="scan"
        )
        all_energies = np.nan
        alpha_epos = np.nan
        try:
            # TODO: move actual calculation into property functions?!
            results = calc.get_all_energies(geom.atoms, coords)
            all_energies = results["all_energies"]
            assert all_energies.size >= 2, (
                "Energies of at least two states must be calculated!"
            )
        except AttributeError as err:
            if not dummy_scan:
                raise err
            result = calc.get_energy(atoms, coords)
            energy = result["energy"]
            all_energies = np.array((energy, energy + DUMMY_2EV))
        except CalculationFailedException as err:
            print(err)
            print(f"Calculation {calc_number:03d} failed!")
        except RunAfterCalculationFailedException as err:
            print(err)
            print(f"Postprocessing of calculation {calc_number:03d} failed!")

        try:
            wf = calc.get_stored_wavefunction()
            _, alpha_epos = epos_from_wf(wf, fragments)
        except AttributeError:
            pass
        return all_energies, alpha_epos

    marcus_dim = marcus_dim_data["marcus_dim"]
    normal_coords = marcus_dim_data["normal_coordinates"]
    properties = marcus_dim_data["properties"]
    eigvecs = marcus_dim_data["eigvecs"]
    masses = marcus_dim_data["masses"]

    int_geom = geom.copy(
        coord_type="redund",
        coord_kwargs={
            "bonds_only": True,
        },
    )
    B = int_geom.internal.B

    scan_result = scan_parallel(
        coords_init=geom.cart_coords.copy(),
        get_properties=get_properties,
        marcus_dim=marcus_dim,
        normal_coords=normal_coords,
        properties=properties,
        eigvecs=eigvecs,
        masses=masses,
        B=B,
        rd_class=rd_class,
        scheduler=scheduler,
        out_dir=cwd,
        **scan_kwargs,
    )
    dump_scan_result(geom, scan_result, dummy_scan, cwd=cwd)
    return scan_result


def run_param(
    scan_factors: np.ndarray,
    scan_energies: np.ndarray,
    mult: int,
    cwd: Union[os.PathLike, str] = ".",
) -> dict[str, mdtypes.MarcusModel]:
    """Try to parametrize quadratic Marcus model from scan data (factors & energies).

    Parameters
    ----------
    scan_factors
        1d-array of shape (n, ) containing the multiplicative factors, describing the displacment
        along the Marcus dimension.
    scan_energies
        2d-array of shape(n, 2) containing the energies of the ground and first excited state.
    mult
        Integer representing the multiplicity of the system under study. Used to create
        meaningful labels in the plots.
    cwd
        String or Path object that controls the directory, where various files
        are written.

    Returns
    -------
    models
        Dictionary with two str keys ("a", "b") containing the (possibly) parametrized
        Marcus models. Parametrization is only possible for class II systems.
    """
    cwd = Path(cwd)
    # 3.) Parametrization of Marcus model
    print(highlight_text("Parametrization of Marcus Model"))
    models = param_marcus(scan_factors, scan_energies)
    shifted_scan_energies = scan_energies - scan_energies.min()
    for para, model in models.items():
        print(f"Parametrization {para}: {model.pretty()}")
        fig, ax = model.plot(scan_factors)
        for i, state in enumerate(shifted_scan_energies.T):
            label = get_state_label(mult, i)
            ax.plot(scan_factors, state, label=f"${label}$ scan")
            ax.legend()
        fig_fn = cwd / f"marcus_model_{para}.svg"
        fig.savefig(fig_fn)
        plt.close(fig)
        print(f"Saved plot of Marcus model to '{fig_fn}'")
    return models


def run_marcus_dim(
    geom: Geometry,
    calc_getter: Callable,
    fragments: list[list[int]],
    rd_class: mdtypes.RobinDay,
    fit_kwargs: Optional[dict] = None,
    scan_kwargs: Optional[dict] = None,
    T: float = 300.0,
    cluster: bool = True,
    force: bool = False,
    dummy_scan: bool = False,
    cwd: Union[os.PathLike, str] = ".",
) -> dict[str, mdtypes.MarcusModel]:
    """Fit Marucs dimension, scan along it and try to parametrize the Marcus model

    Parameters
    ----------
    geom
        Equilibrium geometry with N atoms and 3N Cartesian coordinates.
    calc_getter
        Callable that returns a new calculator. Takes optional kwargs,
        that are passed to the underlying Calculator class.
    fragments
        List of lists of integers, containing the different fragment indices.
        Order is important, as the prefactors for weighting the alpha-electron
        spin density depend on the order.
    rd_class
        RobinDay enum, denoting the Robin-Day class.
    fit_kwargs
        Optional dict, containing additional arguments for fitting the Marcus dimension.
    scan_kwargs
        Optional dict, containing additional arguments for scanning along the Marcus dimension.
    T
        Temperature in Kelvin.
    cluster
        Boolean; controls wheter all calculations are done in serial or in parallel
        using dask.distributed.
    force
        Boolean; when True fit & scan are done, even when the results are already present.
    dummy_scan
        Optional; see docstring of run_scan()/run_scan_parallel().
    cwd
        String or Path object that controls the directory, where various files
        are written.
    """
    rd_class = mdtypes.get_rd_class(rd_class)

    if fit_kwargs is None:
        fit_kwargs = {}
    if scan_kwargs is None:
        scan_kwargs = {}
    fit_kwargs = fit_kwargs.copy()
    scan_kwargs = scan_kwargs.copy()
    assert T > 0.0, f"Temperature {T=} must not be negative!"
    cwd = Path(cwd)

    # When no explicit property function choice is made, we derive the correct
    # property function from the given Robin-Day class.
    try:
        fit_kwargs["property"]
    except KeyError:
        property = {
            mdtypes.RobinDay.CLASS2: mdtypes.Property.EEXC,
            mdtypes.RobinDay.CLASS3: mdtypes.Property.EPOS_IAO,
        }[rd_class]
        fit_kwargs["property"] = property
        print(
            "No property was selected for the fit! Using the default property "
            f"for {rd_class}: {property}"
        )

    #############################
    # 1.) Fit Marcus dimension  #
    #############################

    print(highlight_text("Fitting of Marcus Dimension"))

    fit_results_path = cwd / FIT_RESULTS_FN
    rms_converged = False
    one_batch = False
    if fit_results_path.exists():
        marcus_dim_data = np.load(fit_results_path)
        try:
            rms_converged = bool(marcus_dim_data["rms_converged"])
            one_batch = int(marcus_dim_data["max_batches"]) == 1
        except KeyError:
            print("Could not determine if Marcus dimension already converged!")
    # When rms_converged is True md_results will always be present
    if not force and (rms_converged or one_batch):
        marcus_dim = marcus_dim_data["marcus_dim"]
        print(f"Loaded Marcus dimension from '{fit_results_path}'. Skipping fit.")
        print(
            f"Final Marcus dimension:\n{geom.as_xyz(cart_coords=marcus_dim, comment=None)}"
        )
    else:
        with MaybeScheduler(cluster) as scheduler:
            marcus_dim_data = fit_marcus_dim(
                geom, calc_getter, fragments, T=T, scheduler=scheduler, **fit_kwargs
            )
    print()

    ###################################
    # 2.) Scan along Marcus dimension #
    ###################################

    print(highlight_text("Scan along Marcus Dimension"))

    scan_results_path = cwd / SCAN_RESULTS_FN
    scan_done = False
    if scan_results_path.exists():
        scan_results = np.load(scan_results_path)
        scan_factors = scan_results["factors"]
        scan_coords = scan_results["coords"]
        scan_energies = scan_results["energies"]
        scan_properties = scan_results["properties"]
        scan_result = mdtypes.MarcusDimScanResult(
            factors=scan_factors,
            coords=scan_coords,
            energies=scan_energies,
            properties=scan_properties,
        )
        try:
            scan_done = bool(scan_results["scan_done"])
        except KeyError:
            pass

    if not scan_done:
        with MaybeScheduler(cluster) as scheduler:
            scan_result = run_scan_parallel(
                geom,
                marcus_dim_data,
                calc_getter,
                fragments,
                rd_class=rd_class,
                dummy_scan=dummy_scan,
                scheduler=scheduler,
                cwd=".",
                **scan_kwargs,
            )
        """
        scan_result = run_scan(
            geom,
            marcus_dim_data["marcus_dim"],
            calc_getter,
            fragments,
            dummy_scan=dummy_scan,
            cwd=".",
        )
        """
    else:
        print(f"Loaded scan data from '{scan_results_path}'. Skipping scan.")
    print()

    #######################################
    # 3.) Parametrization of Marcus model #
    #######################################

    # Try to determine the multiplicity, for the creation of state labels as D_0, D_1, etc.
    # for the following plots.
    mult_calc = calc_getter()
    try:
        mult = mult_calc.mult
    except AttributeError:
        mult = -1
    # Parametrize using only the energies of ground and 1st excited state.
    models = run_param(scan_result.factors, scan_result.energies[:, :2], mult, cwd=cwd)
    print()
    return models
