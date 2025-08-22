import dataclasses
import functools
from typing import Optional

import h5py
import numpy as np

from pysisyphus.constants import AU2KJPERMOL
from pysisyphus.drivers.opt import run_opt
from pysisyphus.Geometry import Geometry
from pysisyphus.helpers_pure import highlight_text
from pysisyphus.optimizers.RFOptimizer import RFOptimizer
from pysisyphus.TablePrinter import TablePrinter


def relaxed_scan(
    geom,
    calc_getter,
    constrain_prims,
    target_values,
    title,
    max_cycles=25,
    trust_radius=0.5,
    thresh=1e-2,
    dump=True,
):
    """Relaxed scan, allowing fixing of multiple primitive internals."""

    # Constrain desired primitives
    copy_kwargs = {
        "coord_type": "redund",
        "coord_kwargs": {"constrain_prims": constrain_prims},
    }
    scan_geom = geom.copy(**copy_kwargs)
    calc = calc_getter()
    scan_geom.set_calculator(calc)

    lengths_ = [len(inds) for prim_type, *inds in constrain_prims]
    if all([True for l in lengths_ if l == 2]):
        unit = "au"
    elif all([True for l in lengths_ if l in (3, 4)]):
        unit = "rad"
    else:
        unit = "au (rad)"

    typed_prims = scan_geom.internal.typed_prims
    constr_inds = [typed_prims.index(cp) for cp in constrain_prims]
    constr_fmt = lambda arr: f"({np.array2string(arr, precision=4)}) {unit}"

    def print_constraints():
        print(f"Desired constraints: {constr_fmt(target_values)}")
        print(f" Actual constraints: {constr_fmt(scan_geom.coords[constr_inds])}")

    trj = ""
    scan_cart_coords = [scan_geom.cart_coords.copy()]
    scan_energies = [scan_geom.energy]
    actual_values = [scan_geom.coords[constr_inds]]

    for i in range(max_cycles):
        print(highlight_text(f"Step {i}", level=1))
        cur_diff = target_values - scan_geom.coords[constr_inds]
        step = cur_diff
        step_norm = np.linalg.norm(step)
        print_constraints()
        print(f"         Difference: {step_norm:.6f} {unit}\n")
        if step_norm <= thresh:
            print(
                f"Relaxed scan converged! Norm of proposed step <= {thresh:.4f} {unit}"
            )
            break
        if step_norm > trust_radius:
            step = trust_radius * step / step_norm
        new_coords = scan_geom.coords.copy()
        new_coords[constr_inds] += step
        scan_geom.set_coords(new_coords, update_constraints=True)

        prefix = f"{title}_step_{i}"
        opt_kwargs = {
            "prefix": prefix,
            "h5_group_name": prefix,
            "dump": dump,
        }
        opt = RFOptimizer(scan_geom, **opt_kwargs)
        opt.run()

        scan_cart_coords.append(scan_geom.cart_coords.copy())
        scan_energies.append(scan_geom.energy)
        actual_values.append(scan_geom.coords[constr_inds])
        trj += scan_geom.as_xyz() + "\n"

    scan_cart_coords = np.array(scan_cart_coords)
    scan_energies = np.array(scan_energies)
    actual_values = np.array(actual_values)

    if dump:
        trj_fn = f"{title}_scan.trj"
        with open(trj_fn, "w") as handle:
            handle.write(trj)
        print(f"Dumped optimized geometries to '{trj_fn}'.")

        group_name = f"scan_{title}"
        with h5py.File("scan.h5", "a") as handle:
            try:
                del handle[group_name]
            except KeyError:
                pass
            group = handle.create_group(group_name)
            group.create_dataset("energies", data=scan_energies)
            group.create_dataset("cart_coords", data=scan_cart_coords)
            group.create_dataset("target_values", data=target_values)
            group.create_dataset("actual_values", data=actual_values)
            dt = h5py.vlen_dtype(np.dtype("int32"))
            cp = np.array([(pt.value, *ind) for pt, *ind in constrain_prims])
            cp_dataset = group.create_dataset(
                "constrain_prims", (len(target_values),), dtype=dt
            )
            cp_dataset[:] = cp
    return scan_geom, scan_cart_coords, scan_energies


def concat(self, other):
    return np.concatenate((self[::-1], other))


@dataclasses.dataclass
class ScanResult:
    vals: np.ndarray
    energies: np.ndarray
    target_vals: np.ndarray
    geoms: list[Geometry]
    unit: str
    name: Optional[str] = None

    def reverse_and_add(self, other, name: Optional[str] = None):
        assert self.unit == other.unit
        return ScanResult(
            vals=concat(self.vals, other.vals),
            energies=concat(self.energies, other.energies),
            target_vals=concat(self.target_vals, other.target_vals),
            geoms=self.geoms[::-1] + other.geoms,
            unit=self.unit,
            name=name,
        )


def relaxed_1d_scan(
    geom,
    calc_getter,
    constrain_prims,
    start,
    step_size,
    steps,
    opt_key,
    opt_kwargs,
    pref=None,
    callback=None,
):
    assert geom.coord_type == "redund", "'coord_type: redund' is required!"
    pref_org = pref
    if pref is None:
        pref = ""
    else:
        pref = f"{pref}_"

    constr_prim = constrain_prims[0]
    copy_kwargs = {
        "coord_type": "redund",
        "coord_kwargs": {
            "constrain_prims": constrain_prims,
            "recalc_B": True,
        },
    }
    constr_geom = geom.copy(**copy_kwargs)
    constr_ind = constr_geom.internal.get_index_of_typed_prim(constrain_prims[0])
    calc = calc_getter()

    # Always return same calculator, so orbitals can be reused etc.
    def _calc_getter():
        return calc

    # Drop PrimType at index 0
    # TODO: just to check for number of indices does not seem very robust ...
    unit = "au" if len(constr_prim[1:]) == 2 else "rad"

    org_val = constr_geom.coords[constr_ind]
    cur_val = start
    end_val = start + step_size * steps
    target_scan_vals = np.linspace(cur_val, end_val, steps + 1)
    print(
        f"    Coordinate: {constr_prim}\n"
        f"Original value: {org_val:.4f} {unit}\n"
        f" Initial value: {cur_val:.4f} {unit}\n"
        f"   Final value: {end_val:.4f} {unit}\n"
        f"         Steps: {steps}+1\n"
        f"     Step size: {step_size:.4f} {unit}\n"
    )
    scan_vals = list()
    scan_geoms = list()
    scan_energies = list()
    scan_xyzs = list()

    # Start of scan loop
    for cycle, cur_val in enumerate(target_scan_vals):
        opt_kwargs_ = opt_kwargs.copy()
        name = f"{pref}relaxed_scan_{cycle:04d}"
        opt_kwargs_["prefix"] = name
        opt_kwargs_["h5_group_name"] = name

        # Update constrained coordinate
        new_coords = constr_geom.coords
        new_coords[constr_ind] = cur_val
        constr_geom.set_coords(new_coords, update_constraints=True)

        title = f"{pref}Step {cycle:02d}, coord={cur_val:.4f} {unit}"
        opt_result = run_opt(
            constr_geom, _calc_getter, opt_key, opt_kwargs_, title=title, level=1
        )
        if callback is not None:
            callback(opt_result)
        # Continue w/ updated geometry
        constr_geom = opt_result.geom
        # Keep (copy of) optimized XYZ coordinates, energy and geometry
        scan_xyzs.append(constr_geom.as_xyz())
        scan_vals.append(constr_geom.coords[constr_ind])
        scan_energies.append(constr_geom.energy)
        scan_geoms.append(constr_geom.copy())

        # This will rewrite the data in every cycle, but this is not a problem
        # as the amount of data is very limited.
        scan_data = np.stack((scan_vals, scan_energies), axis=1)
        np.savetxt(f"{pref}relaxed_scan.dat", scan_data)

        if not opt_result.opt.is_converged:
            print(f"Step {cycle} did not converge. Breaking!")
            break
    # End of scan loop

    with open(f"{pref}relaxed_scan.trj", "w") as handle:
        handle.write("\n".join(scan_xyzs))

    scan_vals, scan_energies = scan_data.T
    scan_res = ScanResult(
        vals=scan_vals,
        energies=scan_energies,
        target_vals=target_scan_vals,
        geoms=scan_geoms,
        unit=unit,
        name=pref_org,
    )
    return scan_res


@functools.singledispatch
def print_summary(
    target_scan_vals: np.ndarray,
    scan_vals,
    scan_energies,
    unit,
    name: Optional[str] = None,
):
    if name is None:
        name = ""
    else:
        name += " "

    scan_energies_rel = scan_energies - scan_energies.min()
    scan_energies_rel *= AU2KJPERMOL

    print(highlight_text(f"{name}scan summary", level=1))
    header = (
        "Step",
        f"Target / {unit}",
        f"Actual / {unit}",
        f"Δ / {unit}",
        "E / au",
        "ΔE / kJ mol⁻¹",
    )
    col_fmts = ("{:3d}", "float", "float", "{: >12.4e}", "float", "float")
    table = TablePrinter(header, col_fmts, width=15)
    table.print_header()

    for step, (target, act, en, en_rel) in enumerate(
        zip(target_scan_vals, scan_vals, scan_energies, scan_energies_rel)
    ):
        delta = target - act
        table.print_row((step, target, act, delta, en, en_rel))


@print_summary.register
def _(scan_result: ScanResult):
    print_summary(
        scan_result.target_vals,
        scan_result.vals,
        scan_result.energies,
        scan_result.unit,
        scan_result.name,
    )
