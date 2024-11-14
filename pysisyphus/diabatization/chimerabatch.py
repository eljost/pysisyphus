import os
from pathlib import Path
import stat

import numpy as np

from pysisyphus.xyzloader import make_xyz_str_au


def set_executable_flag(fn):
    current_permissions = os.stat(fn).st_mode
    # new_permissions = current_permissions | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH
    new_permissions = current_permissions | stat.S_IXUSR
    os.chmod(fn, new_permissions)


# Python 3.14 template strings would be useful here, as I only want to set the
# color and not cubs/xyz yet.
# DENS_CUB_TPL = """
# chimerabatch.py chimera.py --iso 0.001 --cubes {cubs} --xyzs {xyz} --color2 {color}
# """.strip()
# Just supplying color does not work, as cubs and xyz are missing ...
# DET_TPL = DENS_CUB_TPL.format(color="1.0 0.0 0.0")
# ATT_TPL = DENS_CUB_TPL.format(color="0.0 0.0 1.0")

DENS_CUB_TPL = """
chimerabatch.py chimera.py --iso 0.001 --cubes {cubs} --xyzs {xyz} --color2 
""".strip()
DET_TPL = DENS_CUB_TPL + " 1.0 0.0 0.0"
ATT_TPL = DENS_CUB_TPL + " 0.0 0.0 1.0"


def dens_cub_call(cub_fns, xyz_fn, tpl):
    cub_names = " ".join(map(lambda path: path.name, cub_fns))
    cmd = tpl.format(cubs=cub_names, xyz=xyz_fn)
    return cmd


def det_att_call(fns, xyz_fn):
    nfns = len(fns)
    assert nfns % 2 == 0
    npairs = nfns // 2
    det_fns = fns[:npairs]
    att_fns = fns[npairs:]

    det_cmd = dens_cub_call(det_fns, xyz_fn, DET_TPL)
    att_cmd = dens_cub_call(att_fns, xyz_fn, ATT_TPL)
    cmds = [det_cmd, att_cmd]
    return cmds


def write_inp(
    cube_fns: dict[str, list[Path]],
    atoms: tuple[str, ...],
    coords: np.ndarray,
    base_name: str,
    out_dir: Path,
):
    # Dump geometry if it does not exist
    xyz_fn = out_dir / Path(base_name).with_suffix(".xyz")
    if not xyz_fn.exists():
        xyz = make_xyz_str_au(atoms, coords.reshape(-1, 3))
        with open(xyz_fn, "w") as handle:
            handle.write(xyz)
    else:
        print(f"'{str(xyz_fn)}' already exists. Skipping creation.")

    cmds = list()
    # Detachment/attachment densities
    keys = cube_fns.keys()
    da_keys = [key for key in keys if key.endswith("_DA")]
    # Only pass the actual name of the xyz file, as it will be written
    # out_dir, next to the run_chimerabatch.sh.
    xyz_name = xyz_fn.name
    for da_key in da_keys:
        fns = cube_fns[da_key]
        da_cmds = det_att_call(fns, xyz_name)
        cmds.extend(da_cmds)
    # TODO: also handle spin densities

    cb_fn = out_dir / "run_chimerabatch.sh"
    with open(cb_fn, "w") as handle:
        text = "\n".join(["#!/usr/bin/env bash"] + cmds)
        handle.write(text)
    set_executable_flag(cb_fn)
    return cb_fn
