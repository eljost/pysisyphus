import collections.abc
from difflib import SequenceMatcher
from enum import Enum
import itertools as it
import json
import logging
import math
from pathlib import Path
from typing import List, Tuple
import re
import time
from typing import Optional

import numpy as np
import psutil

from pysisyphus.config import p_DEFAULT, T_DEFAULT
from pysisyphus.constants import AU2J, BOHR2ANG, C, R, AU2KJPERMOL, NA


"""Functions defined here don't import anything from pysisyphus, besides
the constants module, but only from the stdlib and from third parties."""


HASH_PREC = 4


def eigval_to_wavenumber(ev):
    # This approach seems numerically more unstable
    # conv = AU2J / (AMU2KG * BOHR2M ** 2) / (2 * np.pi * 3e10)**2
    # w2nu = np.sign(ev) * np.sqrt(np.abs(ev) * conv)
    # The two lines below are adopted from Psi4 and seem more stable,
    # compared to the approach above.
    conv = np.sqrt(NA * AU2J * 1.0e19) / (2 * np.pi * C * BOHR2ANG)
    w2nu = np.sign(ev) * np.sqrt(np.abs(ev)) * conv
    return w2nu


def hash_arr(arr, precision=HASH_PREC):
    str_ = np.array2string(arr, precision=precision)
    return hash(str_)


def hash_args(*args, precision=HASH_PREC):
    hashes = list()
    for arg in args:
        try:
            hash_ = hash(arg)
        except TypeError:
            hash_ = hash_arr(arg, precision=precision)
        hashes.append(hash_)
    return hash(tuple(hashes))


def log(logger, msg, level=logging.DEBUG):
    if logger is not None:
        logger.log(level, msg)


def sort_by_central(set1, set2):
    """Determines a common index in two sets and returns a length 3
    tuple with the central index at the middle position and the two
    terminal indices as first and last indices."""
    central_set = set1 & set2
    union = set1 | set2
    assert len(central_set) == 1
    terminal1, terminal2 = union - central_set
    (central,) = central_set
    return (terminal1, central, terminal2), central


def merge_sets(fragments):
    """Merge a list of iterables."""
    # Hold the final fragments that can't be merged further, as they
    # contain distinct atoms.
    fragments = [frozenset(frag) for frag in fragments]
    merged = list()
    while len(fragments) > 0:
        popped = fragments.pop(0)
        # Look for an intersection between the popped unmerged fragment
        # and the remaining unmerged fragments.
        for frag in fragments:
            if popped & frag:
                fragments.remove(frag)
                # If a intersecting unmerged fragment is found merge
                # both fragments and append them at the end.
                fragments.append(popped | frag)
                break
        else:
            # Add the unmerged fragment into merged if it doesn't
            # intersect with any other unmerged fragment.
            merged.append(popped)
    return merged


def remove_duplicates(seq):
    tuples = [tuple(itm) for itm in seq]
    seen = set()
    seen_add = seen.add
    return [itm for itm in tuples if not (itm in seen or seen_add(itm))]


class OrderedEnum(Enum):
    def __ge__(self, other):
        if self.__class__ is other.__class__:
            return self.value >= other.value
        return NotImplemented

    def __gt__(self, other):
        if self.__class__ is other.__class__:
            return self.value > other.value
        return NotImplemented

    def __le__(self, other):
        if self.__class__ is other.__class__:
            return self.value <= other.value
        return NotImplemented

    def __lt__(self, other):
        if self.__class__ is other.__class__:
            return self.value < other.value
        return NotImplemented

    def __str__(self):
        return self.name


def timed(logger):
    def decorator(func):
        def wrapper(*args, **kwargs):
            start = time.time()
            result = func(*args, **kwargs)
            end = time.time()
            duration = end - start
            log(logger, f"Execution of '{func.__name__}' took {duration:.2f} s.")
            return result

        return wrapper

    return decorator


def get_input(data, prompt, lbl_func=None):
    if lbl_func is None:
        lbl_func = lambda _: _
    labels = [lbl_func(d) for d in data]
    print(prompt)
    while True:
        for i, l in enumerate(labels):
            print(f"{i: >3d}: {l}")
        try:
            inp = int(input("Selection: "))
            if not (0 <= inp < len(labels)):
                raise ValueError
            break
        except ValueError:
            print("Invalid input!")
            print()
    return data[inp]


def expand(to_expand):
    if any([isinstance(to_expand, cls) for cls in (list, tuple, np.ndarray)]):
        return to_expand
    elif ".." in to_expand:
        start, end = [int(i) for i in to_expand.split("..")]
        return list(range(start, end))
    # Numbers
    else:
        return [int(to_expand)]


def full_expand(to_expand):
    try:
        split = to_expand.strip().split(",")
        expanded = list(it.chain(*[expand(te) for te in split]))
    except AttributeError:
        try:
            expanded = list(to_expand)
        except TypeError:
            expanded = [
                to_expand,
            ]
    return expanded


def file_or_str(*args, method=False):
    exts = args

    def inner_func(func):
        def wrapped(inp, *args, **kwargs):
            if method:
                obj = inp
                inp, *args = args
            p = Path(inp)
            looks_like_file = exts and (p.suffix in exts)
            if looks_like_file and p.is_file():
                with open(p) as handle:
                    inp = handle.read()
            elif looks_like_file and not p.exists():
                raise FileNotFoundError(
                    f"{inp} looks like a file/path, but it does not exist!"
                )
            if method:
                res = func(obj, inp, *args, **kwargs)
            else:
                res = func(inp, *args, **kwargs)
            return res

        return wrapped

    return inner_func


def recursive_update(d, u):
    """Recursive update of d with keys/values from u.

    From https://stackoverflow.com/questions/3232943
    """
    if u is None:
        return d

    # Best I can do is ... secretly transform hierarchical input
    # to old-style flat input.
    #
    # Try to transform new-style hierarchical input
    # to the old-style flat input. If this new-style proves useful
    # and everything is refactored accordingly the try/except could
    # be dropped. But now everything is still in the old-style.
    try:
        keys = list(u["type"].keys())
        if len(keys) == 1:
            key = keys[0]
            kwargs = u["type"][key]
            u["type"] = key
            try:
                u.update(kwargs)
            # Raised when kwargs is None, e.g., the input is hierarchical,
            # but no further keywords are provided for the given type.
            except TypeError:
                pass
    # Raised when u has no "type" key
    except KeyError:
        pass
    except AttributeError:
        pass

    for k, v in u.items():
        if isinstance(v, collections.abc.Mapping):
            d[k] = recursive_update(d.get(k, {}), v)
        else:
            d[k] = v
    return d


def report_isotopes(geom, affect_str):
    if (geom.isotopes is not None) and len(geom.isotopes) > 0:
        print(f"Different isotopes were requested! This will affect {affect_str}.")
        atoms = geom.atoms
        masses = geom.masses
        for atom_ind, _ in geom.isotopes:
            print(f"\tAtom {atom_ind}{atoms[atom_ind]}: {masses[atom_ind]:.6f} au")
        print()


def report_frozen_atoms(geom):
    if hasattr(geom, "freeze_atoms") and len(geom.freeze_atoms) > 0:
        print(f"Coordinates of {len(geom.freeze_atoms)} atoms are frozen:")
        atoms = geom.atoms
        coords3d = geom.coords3d * BOHR2ANG
        fmt = " >12.8f"
        for atom_ind in geom.freeze_atoms:
            atom = atoms[atom_ind]
            x, y, z = coords3d[atom_ind]
            print(f"\t{atom_ind:03d} {atom} {x:{fmt}} {y:{fmt}} {z:{fmt}}")
        print()


def highlight_text(text, width=80, level=0):
    levels = {
        #  horizontal
        #       vertical
        #           corner
        0: ("#", "#", "#"),
        1: ("-", "|", "+"),
    }
    full_length = len(text) + 4
    pad_len = width - full_length
    pad_len = (pad_len - (pad_len % 2)) // 2
    pad = " " * pad_len
    hchar, vchar, cornerchar = levels[level]
    full_row = cornerchar + (hchar * (full_length - 2)) + cornerchar
    highlight = (
        f"""{pad}{full_row}\n{pad}{vchar} {text.upper()} {vchar}\n{pad}{full_row}"""
    )
    return highlight


def interpolate_colors(values, c1, c2, num=32):
    """Expects two RGB colors c1 and c2."""
    c_diff = c2 - c1
    step = c_diff / (num - 1)
    colors = (c1 + np.arange(num)[:, None] * step).astype(int)

    # Map value interval onto interval range(num)
    # y = m*x + n
    val_min = values.min()
    val_max = values.max()
    m = abs((num - 1) / (val_min - val_max))
    n = -m * val_min
    inds = np.around(m * values + n).astype(int)
    rgb_colors = colors[inds]
    hex_colors = [f"#{r:02x}{g:02x}{b:02x}" for r, g, b in rgb_colors]
    return rgb_colors, hex_colors


def get_molecular_radius(coords3d, min_offset=0.9452):
    coords3d = coords3d.copy()
    mean = coords3d.mean(axis=0)
    coords3d -= mean[None, :]
    distances = np.linalg.norm(coords3d, axis=1)
    std = max(min_offset, np.std(distances))  # at least 2 angstrom apart
    radius = distances.mean() + 2 * std
    return radius


def filter_fixture_store(test_name):
    def inner_function(function):
        def wrapper(fixture_store):
            rb = fixture_store["results_bag"]
            filtered = {
                "results_bag": {k: v for k, v in rb.items() if k.startswith(test_name)}
            }
            return function(filtered)

        return wrapper

    return inner_function


def get_clock():
    ref = time.time()

    def clock(msg=""):
        nonlocal ref
        now = time.time()
        dur = now - ref
        ref = now
        print(f"{msg: >32}, {dur:.3f} s since last call!")

    return clock


def chunks(l, n):
    """Yield successive n-sized chunks from l.
    https://stackoverflow.com/a/312464
    """
    for i in range(0, len(l), n):
        yield l[i : i + n]


def describe(arr):
    shape = arr.shape
    min_ = arr.min()
    max_ = arr.max()
    mean = np.mean(arr)
    median = np.median(arr)
    var = np.var(arr)
    fmt = ".4f"
    return (
        f"shape={shape}, min={min_:{fmt}}, mean={mean:{fmt}}, "
        f"median={median:{fmt}}, max={max_:{fmt}}, variance={var:.4e}"
    )


def touch(fn):
    try:
        Path(fn).touch()
    except IsADirectoryError:
        pass


def approx_float(
    num: float, expected: float, abs_tol: float = 1e-6, rel_tol: float = 1e-12
) -> bool:
    def pos_or_none(num, name):
        assert (num > 0.0) or (num is None), f"{name} must be positive or None!"

    pos_or_none(abs_tol, "Absolute tolerance")
    pos_or_none(rel_tol, "Relative tolerance")
    assert abs_tol or rel_tol

    if rel_tol is None:
        rel_tol = 0.0
    rel_tol = rel_tol * expected
    tolerance = max(abs_tol, rel_tol)
    return abs(num - expected) <= tolerance


_CONV_FUNCS = {
    # First item in value converts to JSON dumpable type,
    # second items converts from JSON to original type.
    "energy": (lambda en: float(en), lambda en: float(en)),
    "forces": (
        lambda forces: forces.tolist(),
        lambda forces: np.array(forces, dtype=float).flatten(),
    ),
    "hessian": (
        lambda hessian: hessian.tolist(),
        lambda hessian: np.array(hessian, dtype=float),
    ),
}


def results_to_json(results):
    conv_results = {key: _CONV_FUNCS[key][0](val) for key, val in results.items()}
    return json.dumps(conv_results)


def json_to_results(as_json):
    results = {
        key: _CONV_FUNCS[key][1](val) for key, val in json.loads(as_json).items()
    }
    return results


def standard_state_corr(T=T_DEFAULT, p=p_DEFAULT, n=1):
    """dG for change of standard state from gas to solution of 1 mol/l"""
    Vm = n * R * T / p  # in m³
    Vm_litres = Vm * 1000
    dG = R * T * math.log(Vm_litres) / 1000  # kJ mol⁻¹
    dG_au = dG / AU2KJPERMOL
    return dG_au


def check_mem(mem, pal, avail_frac=0.85, logger=None):
    virt_mem = psutil.virtual_memory()
    mb_available = virt_mem.available * avail_frac / 1024 / 1024
    mb_requested = mem * pal
    msg = f"{mb_available:.2f} MB memory available, {mb_requested:.2f} MB requested."
    if mb_requested > mb_available:
        mb_corr = int(mb_available / pal)
        msg += f" Too much memory requested. Using smaller value of {mb_corr} MB."
    else:
        mb_corr = mem
    log(logger, msg)

    return mb_corr


def get_ratio(str_: str, comp_str: str) -> str:
    """See https://stackoverflow.com/a/17388505"""
    return SequenceMatcher(None, str_, comp_str).ratio()


def find_closest_sequence(str_: str, comp_strs: List[str]) -> Tuple[str, float]:
    ratios = [get_ratio(str_, comp_str) for comp_str in comp_strs]
    argmax = np.argmax(ratios)
    max_ratio = ratios[argmax]
    best_match = comp_strs[argmax]
    return (best_match, max_ratio)


def increment_fn(org_fn: str, suffix: Optional[str]=None) -> str:
    """
    Append, or increase a suffixed counter on a given filename.
    If no counter is present it will be set to zero. Otherwise
    it is incremented by one.

    >>> increment_fn("opt", "rebuilt")
    'opt_rebuilt_000'
    >>> increment_fn("opt")
    'opt_000'
    >>> increment_fn("opt_rebuilt_000", "rebuilt")
    'opt_rebuilt_001'

    Parameters
    ----------
    org_fn
        The original, unaltered filename.
    suffix
        Optional suffix to be append.

    Returns
    -------
    incr_fn
        Modified filename with optional suffix and incremented counter.
    """
    if suffix is not None and len(suffix) > 0:
        suffix = f"_{suffix}"
    else:
        suffix = ""

    # Check if prefix is present
    regex = re.compile(rf"(.+)({suffix}_(\d+))")
    if mobj := regex.match(org_fn):
        org_fn, _, counter = mobj.groups()
        counter = int(counter)
        counter += 1
    else:
        counter = 0

    incr_fn = f"{org_fn}{suffix}_{counter:03d}"
    return incr_fn


if __name__ == "__main__":
    import doctest

    doctest.testmod()
