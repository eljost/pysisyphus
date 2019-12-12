#!/usr/bin/env python3

import subprocess
import tempfile
from pathlib import Path

from pysisyphus.config import get_cmd


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
CUBE_TPL= ("load {cube_fn}" + TPL_BASE + """
isosurface cutoff {isoval} sign {colors} "{cube_fn}"
write image pngt "{png_fn}"
""")


def call_jmol(spt_str):
    spt_fn = "jmol.spt"
    with tempfile.NamedTemporaryFile("w", suffix=".spt") as spt_handle:
        spt_handle.write(spt_str)
        spt_handle.flush()
        jmol_cmd = [get_cmd("jmol"), "-n", spt_handle.name]
        proc = subprocess.Popen(jmol_cmd, universal_newlines=True,
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE)
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


if __name__ == "__main__":
    cdd = "/scratch/turbontos/11_bz_pure/image_000.005.S_004_CDD.cub"
    render_cdd_cube(cdd)
