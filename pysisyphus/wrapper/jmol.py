import subprocess
import tempfile
from pathlib import Path

import numpy as np

from pysisyphus.config import get_cmd
from pysisyphus.constants import BOHR2ANG
from pysisyphus.helpers_pure import interpolate_colors


TPL_BASE = """
{orient}

function _setModelState() {{

    select;
    Spacefill 0.0;

    frank off;
    font frank 16.0 SansSerif Plain;
    select *;
    set fontScaling false;
    background white
    frank off
    set showhydrogens True;

}}

_setModelState;

"""
CUBE_TPL = (
    "load {cube_fn}"
    + TPL_BASE
    + """
isosurface cutoff {isoval} sign {colors} "{cube_fn}"
write image pngt "{png_fn}"
"""
)


def call_jmol(spt_str, show=False):
    with tempfile.NamedTemporaryFile("w", suffix=".spt") as spt_handle:
        spt_handle.write(spt_str)
        spt_handle.flush()
        jmol_cmd = [get_cmd("jmol"), "-n", spt_handle.name]
        if show:
            del jmol_cmd[1]
        proc = subprocess.Popen(
            jmol_cmd,
            universal_newlines=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        proc.wait()
    stdout = proc.stdout.read()
    stderr = proc.stderr.read()
    return stdout, stderr


def render_cdd_cube(cdd_cube, isoval=0.002, orient=""):
    png_fn = Path(cdd_cube).with_suffix(".png")
    spt = CUBE_TPL.format(
        orient=orient,
        cube_fn=cdd_cube,
        isoval=isoval,
        colors="red blue",
        png_fn=png_fn,
    )
    with open("jmol.spt", "w") as handle:
        handle.write(spt)
    stdout, stderr = call_jmol(spt)
    return png_fn


def render_geom_and_charges(geom, point_charges):
    point_charges = point_charges.copy()
    point_charges[:, :3] *= BOHR2ANG
    charges = point_charges[:, -1]

    cr = np.array((255, 0, 0))  # red
    cw = np.array((255, 255, 255))  # white
    cb = np.array((0, 0, 255))  # blue
    chrg_min = charges.min()
    chrg_max = charges.max()
    print(
        f"charges:\n{np.array2string(point_charges, precision=4)}\n"
        f"min(charges): {chrg_min: .4f}\nmax(charges): {chrg_max: .4f}\n"
    )
    c1 = cb if chrg_min < 0.0 else cw
    c2 = cr if chrg_max > 0.0 else cw
    rgb_colors, _ = interpolate_colors(charges, c1, c2)

    # Dump geometry to temporary file
    with tempfile.NamedTemporaryFile("w", suffix=".xyz") as tmp_xyz:
        tmp_xyz.write(geom.as_xyz())
        tmp_xyz.flush()
        spt = f"load {tmp_xyz.name};\n"
        for i, ((x, y, z, pc), (r, g, b)) in enumerate(zip(point_charges, rgb_colors)):
            id_ = f"chrg{i}"
            spt += (
                f"isosurface {id_} center {{{x} {y} {z}}} sphere 0.25;\n"
                f"color ${id_} [{r} {g} {b}];\n"
            )
        print(f"SPT:\n\n{spt}")
        call_jmol(spt, show=True)


if __name__ == "__main__":
    cdd = "/scratch/turbontos/11_bz_pure/image_000.005.S_004_CDD.cub"
    render_cdd_cube(cdd)
