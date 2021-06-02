import h5py
import numpy as np

from pysisyphus.helpers_pure import highlight_text
from pysisyphus.optimizers.RFOptimizer import RFOptimizer


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
    """Relaxed scan driver."""

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
        scan_geom.coords = new_coords

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
