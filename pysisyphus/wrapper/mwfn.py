import logging
import os
from pathlib import Path
import shutil
from subprocess import PIPE, Popen

import numpy as np

from pysisyphus.config import get_cmd
from pysisyphus.constants import AU2EV
from pysisyphus.wrapper.exceptions import SegfaultException


logger = logging.getLogger("mwfn")


def log(msg):
    logger.debug(msg)


def wrap_stdin(stdin):
    return f"<< EOF\n{stdin}\nEOF"


def call_mwfn(inp_fn, stdin, cwd=None):
    if cwd is None:
        cwd = Path(".")
    mwfn_cmd = get_cmd("mwfn")
    cmd = [mwfn_cmd, inp_fn]
    log(f"\n{mwfn_cmd} {inp_fn} {wrap_stdin(stdin)}")
    proc = Popen(
        cmd, universal_newlines=True, stdin=PIPE, stdout=PIPE, stderr=PIPE, cwd=cwd
    )
    stdout, stderr = proc.communicate(stdin)
    if "segmentation fault occurred" in stderr:
        raise SegfaultException(
            "Multiwfn segfaulted! Multiwfn seems to have problems "
            "with systems >= 1000 basis functions. Maybe your system is too big."
        )
    proc.terminate()
    return stdout, stderr


def make_cdd(inp_fn, state, log_fn, cwd=None, keep=False, quality=2, prefix="S"):
    """Create CDD cube in cwd.

    Parameters
    ----------
    inp_fn : str
        Filename of a .molden/.fchk file.
    state : int
        CDD cubes will be generated up to this state.
    log_fn : str
        Filename of the .log file.
    cwd : str or Path, optional
        If a different cwd should be used.
    keep : bool
        Wether to keep electron.cub and hole.cub, default is False.
    quality : int
        Quality of the cube. (1=low, 2=medium, 3=high).
    """

    assert quality in (1, 2, 3)

    msg = (
        f"Requested CDD calculation from Multiwfn for state {state} using "
        f"{inp_fn} and {log_fn}"
    )
    log(msg)

    stdin = f"""18
    1
    {log_fn}
    {state}
    1
    {quality}
    10
    1
    11
    1
    15
    0
    0
    0
    q
    """
    stdout, stderr = call_mwfn(inp_fn, stdin, cwd=cwd)

    if cwd is None:
        cwd = "."
    cwd = Path(cwd)

    cube_fns = ("electron.cub", "hole.cub", "CDD.cub")
    if not keep:
        # always keep CDD.cub
        for fn in cube_fns[:2]:
            full_path = cwd / fn
            os.remove(full_path)
    # Rename cubes according to the current state
    new_paths = list()
    for fn in cube_fns:
        old_path = cwd / fn
        root, ext = os.path.splitext(fn)
        new_path = cwd / f"{prefix}_{state:03d}_{root}{ext}"
        try:
            shutil.copy(old_path, new_path)
            os.remove(old_path)
            new_paths.append(new_path)
        except FileNotFoundError:
            pass
    return new_paths


def get_mwfn_exc_str(energies, Xa, Ya=None, Xb=None, Yb=None, thresh=1e-3):
    """Write plain text input for MWFN according to 3.21.

    As of version 3.8 MWFN does not seem to handle this file when spin
    labels are present, even though it is given as an example in the manual.
    """

    # We only use the spin label for unrestricted calculations, where Xb is present.
    spin_a = "A" if (Xb is not None) else ""
    assert len(energies) == (
        len(Xa) + 1
    ), "Found too few energies. Is the GS energy missing?"
    exc_energies = (energies[1:] - energies[0]) * AU2EV
    # states, occ, virt
    nstates, occ_mos, _ = Xa.shape

    exc_str = ""
    mult = 1
    log(f"Using dummy multiplicity={mult} in get_mwfn_exc_str")

    def set_default(mat):
        if mat is None:
            mat = [None] * nstates
        return mat

    Ya = set_default(Ya)
    Xb = set_default(Xb)
    Yb = set_default(Yb)

    def get_exc_lines(ci_coeffs, arrow, spin=""):
        exc_lines = list()
        for (occ, virt), coeff in np.ndenumerate(ci_coeffs):
            if abs(coeff) < thresh:
                continue
            occ_mo = occ + 1
            virt_mo = occ_mos + 1 + virt
            exc_line = f"{occ_mo:>8d}{spin} {arrow} {virt_mo}{spin}       {coeff: .5f}"
            exc_lines.append(exc_line)
        return exc_lines

    for root_, (xa, ya, xb, yb, exc_en) in enumerate(
        zip(Xa, Ya, Xb, Yb, exc_energies), 1
    ):
        exc_str += f"Excited State {root_} {mult} {exc_en:.4f}\n"
        # Excitations (X vector)
        exc_lines = get_exc_lines(xa, "->", spin_a)
        # De-Excitations (Y vector), if present
        if ya is not None:
            exc_lines += get_exc_lines(ya, "<-", spin_a)
        if xb is not None:
            exc_lines += get_exc_lines(xb, "->", "B")
        if yb is not None:
            exc_lines += get_exc_lines(yb, "<-", "B")
        exc_str += "\n".join(exc_lines)
        exc_str += "\n\n"
    return exc_str
