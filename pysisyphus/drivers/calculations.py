import hashlib
import sys
import time
from typing import Callable, Optional, Sequence

import numpy as np

from pysisyphus.exceptions import (
    CalculationFailedException,
    DifferentAtomOrdering,
    RunAfterCalculationFailedException,
)
from pysisyphus.Geometry import Geometry
from pysisyphus.helpers_pure import (
    recursive_extract,
    highlight_text,
    json_to_results,
    results_to_json,
)
from pysisyphus.io import save_hessian
from pysisyphus import timing


def hash_atoms_and_coords(
    atoms: Sequence[str],
    coords: np.ndarray,
    precision: int,
    charge: Optional[int] = None,
    mult: Optional[int] = None,
    hash_func=hashlib.md5,
) -> str:
    """Wrapper function that returns a hash for given atoms and coordinates."""

    # Truncate coordinates to selected precision. No rounding is done.
    coords_trunc = np.trunc(coords.flatten() * 10**precision).astype(int)
    coords_str = "".join(map(str, coords_trunc))

    # Normalize atom symbols to lowercase
    atom_symbols = [symbol.lower() for symbol in atoms]
    atom_str = "".join(atom_symbols)

    hash_inp = atom_str + coords_str
    if charge is not None:
        hash_inp += str(charge)
    if mult is not None:
        hash_inp += str(mult)
    hash_obj = hash_func(hash_inp.encode("utf-8"))
    hash_hex = hash_obj.hexdigest()
    return hash_hex


def calc_wrapper(geom: Geometry, func_name: str, skip_existing: bool = False):
    calc = geom.calculator
    # Unfortunately, this neglects most of the calculator details...
    cur_hash = hash_atoms_and_coords(
        geom.atoms, geom.cart_coords, precision=8, charge=calc.charge, mult=calc.mult
    )
    hash_fn = (calc.out_dir / cur_hash).with_suffix(".json")
    # Check if results are already present; if so, skip the calculation.
    if skip_existing and hash_fn.exists():
        with open(hash_fn) as handle:
            results = json_to_results(handle.read())
        print(
            f"Found matching hash '{cur_hash}'!\n"
            f"Returning results from '{hash_fn}'."
        )
        return results
    # Otherwise, do the calculation.
    try:
        results = getattr(geom.calculator, func_name)(geom.atoms, geom.cart_coords)
        hash_fn.touch()
    except (
        CalculationFailedException,
        RunAfterCalculationFailedException,
    ) as err:
        print(err)
        return None

    # results dict of MultiCalc will contain keys that can't be dumped yet. So
    # we skip the JSON dumping when KeyError is raised.
    try:
        as_json = results_to_json(results)
        calc = geom.calculator
        # Decrease counter, because it will be increased by 1, w.r.t to the
        # calculation.
        with open(hash_fn, "w") as handle:
            handle.write(as_json)
        json_fn = calc.make_fn("results.json", counter=calc.calc_counter - 1)
        json_fn.symlink_to(hash_fn)
        msg = f"Dumped JSON results to '{hash_fn}'."
    except KeyError:
        msg = "Skipped JSON dump of calculation results!"
    print(msg)

    # Check if Hessians were calculated; if so, save them to HDF5 files.
    hessian_results = recursive_extract(results, "hessian")
    energy_results = recursive_extract(results, "energy")
    for (*rest, _), hessian in hessian_results.items():
        energy_key = (*rest, "energy")
        energy = energy_results[energy_key]
        prefix = ("_".join(rest) if rest else calc.name) + "_"
        h5_fn = f"{prefix}hessian.h5"
        save_hessian(
            h5_fn,
            geom,
            cart_hessian=hessian,
            energy=energy,
        )
        print(f"Dumped hessian to '{h5_fn}'.")
    return results


def run_calculations(
    geoms: list[Geometry],
    calc_getter: Callable,
    run_func: Optional[str] = None,
    one_calculator: bool = False,
    rmsd_thresh: float = 0.75,
    skip_existing: bool = True,
) -> tuple[list[Geometry], list[dict]]:
    """Run calculations for all given geometries.

    Sometimes one just wants to run a series of calculations for list of geometries,
    w/o any optimization or IRC integration etc. This function will be executed when
    the YAML inputs lacks any opt/tsopt/irc/... section.

    Parameters
    ----------
    geoms
        List of geometries. Depending on the values of the other arguments all
        geometries must be/don't have to be compatible (same atom order). When
        only one calculator is utilized all geometries must be compatible.
    calc_getter
        Function that returns a new calculator.
    run_func
        By default 'run_calculation()' will be called for the calculator, but with
        run_func another method name can be specified, e.g., 'get_hessian'.
    one_calculator
        Defaults to false. When enabled all calculations will be carried out using
        the same calculator object. All geometries must be compatible (same atoms
        and same atom order) and parallel calculations are not possible. This is a useful
        option for excited state calculations along a path, e.g., from a COS calculation,
        as all ES data will be dumped in one 'overlap_data.h5' file.
    rmsd_thresh
        RMSD threshold in Bohr. Postive floating point number that is used to determine,
        if two geometries are similar enough for propagating chkfiles.
    skip_existing
        Whether calculations for which results already exist should be carried out again,
        or not. If enabled, calculations results will be read from dumped JSON files.

    Returns
    -------
    geoms
        List of geometries that were used for the calculations.
    all_resuls
        List of calculation results from the different calculations.
    """
    print(highlight_text("Running calculations"))

    func_name = "run_calculation" if run_func is None else run_func

    if one_calculator:
        geom0 = geoms[0]
        for other_geom in geoms[1:]:
            geom0.assert_compatibility(other_geom)
        print("Doing all calculations with the same calculator.")
        calc = calc_getter()
        # Overwrite calc_getter to always return the same calculator
        calc_getter = lambda: calc

    # Set calculators on all geometries
    for geom in geoms:
        geom.set_calculator(calc_getter())

    durations = list()
    ncalcs = len(geoms)

    all_results = list()
    i_fmt = "02d"
    for i, geom in enumerate(geoms):
        print(highlight_text(f"Calculation {i:{i_fmt}}", level=1))
        print(geom)

        dur = time.time()
        results = calc_wrapper(geom, func_name, skip_existing=skip_existing)
        dur = time.time() - dur
        # When a calculation failed results will be None
        if results is None:
            # Subtract one from the number of calculations, as we are now
            # lacking a durations because of the failed calculation.
            ncalcs -= 1
            continue
        all_results.append(results)

        # Report timings etc. when calculation succeeded
        print(f"Calculation took {timing.render(dur)}.")
        durations.append(dur)
        estimate = timing.estimate_runtime(durations, ncalcs)
        if estimate.nremain:
            print(
                f"Estimated remaining runtime for {estimate.nremain} calculations: "
                f"from mean={estimate.render_mean()}; "
                f"from median={estimate.render_median()}"
            )
        sys.stdout.flush()

        # Check if we can propagate chkfiles to the next calculator.
        # There are several requirements:
        #   1.) Atom ordering must match
        #   2.) Geometries must be similar enough (small RMSD)
        #   3.) Propagating chkfiles must be implemented.
        if i < (len(geoms) - 1):
            next_geom = geoms[i + 1]
            try:
                rmsd = geom.rmsd(next_geom)
            except DifferentAtomOrdering:
                continue
            if rmsd > rmsd_thresh:
                print(
                    f"Skipping chkfile propagating as geometries {i} and {i+1} are "
                    f"too dissimilar (rmsd={rmsd: >8.4f} au)."
                )
                continue
            try:
                cur_calculator = geom.calculator
                next_calculator = geoms[i + 1].calculator
                next_calculator.set_chkfiles(cur_calculator.get_chkfiles())
                msg = f"Set chkfiles of calculator {i:{i_fmt}} on calculator {i+1:{i_fmt}}"
            except AttributeError:
                msg = "Calculator does not support set/get_chkfiles!"
            print(msg)
        print()
    print()

    for geom, results in zip(geoms, all_results):
        try:
            geom.set_results(results)
        except KeyError:
            pass

    return geoms, all_results
