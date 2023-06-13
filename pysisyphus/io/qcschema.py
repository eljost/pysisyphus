import functools
import json
from pathlib import Path
from typing import Dict

import numpy as np

from pysisyphus.Geometry import Geometry


@functools.singledispatch
def geom_from_qcschema(qcschema: Dict, **geom_kwargs):
    mol = qcschema["molecule"]
    coords = np.array((mol["geometry"]))
    atoms = mol["symbols"]
    return Geometry(atoms, coords, **geom_kwargs)


# These definitions may seem a bit wild, but I was not able to make
# singledispatch work with 'file_or_str'.


@geom_from_qcschema.register
def _(text: str, **geom_kwargs):
    if (p := Path(text).exists()):
        with open(p, "r") as handle:
            text = handle.read()
    qcschema = json.loads(text)
    return geom_from_qcschema(qcschema, **geom_kwargs)


@geom_from_qcschema.register
def _(path: Path, **geom_kwargs):
    return geom_from_qcschema(str(path), **geom_kwargs)
