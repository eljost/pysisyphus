#!/usr/bin/env python3

import os
from pathlib import Path
import shutil
from subprocess import PIPE, Popen

from pysisyphus.config import get_cmd

def call_mwfn(inp_fn, stdin, cwd=None):
    if cwd is None:
        cwd = Path(".")
    mwfn_cmd = [get_cmd("mwfn"), inp_fn]
    proc = Popen(mwfn_cmd, universal_newlines=True,
                 stdin=PIPE, stdout=PIPE, stderr=PIPE, cwd=cwd)
    stdout, stderr = proc.communicate(stdin)
    proc.terminate()
    return stdout, stderr


def make_cdd(inp_fn, state, log, cwd, keep=False, quality=2):
    """Create CDD cube in cwd.

    Parameters
    ----------
    inp_fn : str
        Filename of a .molden/.fchk file.
    state : int
        CDD cubes will be generated up to this state.
    log : str
        Filename of the .log file.
    cwd : str or Path
        If a different cwd should be used.
    keep : bool
        Wether to keep electron.cub and hole.cub, default is False.
    quality : int
        Quality of the cube. (1=low, 2=medium, 3=high).
    """

    assert quality in (1,2,3)

    stdin = f"""18
    1
    {log}
    {state}

    1
    {quality}
    10
    1
    11
    1
    15
    """
    stdout, stderr = call_mwfn(inp_fn, stdin, cwd=cwd)

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
        new_path = cwd / f"S_{state:03d}_{root}{ext}"
        try:
            shutil.copy(old_path, new_path)
            os.remove(old_path)
            new_paths.append(new_path)
        except FileNotFoundError:
            pass
    return new_paths
