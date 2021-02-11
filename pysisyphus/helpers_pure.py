import collections.abc
from enum import Enum
import logging
from pathlib import Path
import time

import numpy as np

from pysisyphus.constants import AU2J, AMU2KG, BOHR2M


"""Functions defined here don't import anything from pysisyphus, besides
the constants module, but only from the stdlib and from third parties."""


def eigval_to_wavenumber(ev):
    conv = AU2J/(AMU2KG*BOHR2M**2)

    return np.sign(ev) * np.sqrt(np.abs(ev)*conv)/(2*np.pi*3e10)


def hash_arr(arr, precision=4):
    str_ = np.array2string(arr, precision=4)
    return hash(str_)


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
    else:
        start, end = [int(i) for i in to_expand.split("..")]
        return list(range(start, end))


def file_or_str(*args):
    exts = args
    def inner_func(func):
        def wrapped(inp, *args, **kwargs):
            p = Path(inp)
            looks_like_file = exts and (p.suffix in exts)
            if looks_like_file and p.is_file():
                with open(p) as handle:
                    inp = handle.read()
            elif looks_like_file and not p.exists():
                raise Exception(f"{inp} looks like a file/path, but it does not exist!")
            return func(inp, *args, **kwargs)

        return wrapped

    return inner_func


def recursive_update(d, u):
    """From https://stackoverflow.com/questions/3232943"""
    if u is None:
        return d

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


