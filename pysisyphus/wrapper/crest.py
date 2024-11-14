import dataclasses
from pathlib import Path
import subprocess
import tempfile
from typing import Literal, Optional

import numpy as np

from pysisyphus.config import get_cmd
from pysisyphus.constants import AU2KCALPERMOL
from pysisyphus.Geometry import Geometry
from pysisyphus.io import xyz, zip
from pysisyphus.typing import PathLike


# Far from complete ...
Solvent = Literal[
    "ch2cl2",
    "chcl3",
    "acetone",
    "acetonitrile",
    "dmf",
    "dmso",
    "water",
    "thf",
    "methanol",
]


class CrestException(Exception):
    pass


@dataclasses.dataclass
class CrestResult:
    best: Geometry
    conformers: list[Geometry]
    energies: np.ndarray
    clustered: Optional[list[Geometry]] = None
    clustered_energies: Optional[np.ndarray] = None

    def __post_init__(self):
        if self.clustered is None:
            self.clustered = list()
            self.clustered_energies = np.array(())


def run_crest(
    fn: PathLike,
    charge: int,
    alpb: Optional[Solvent] = None,
    pal=1,
    cluster: Optional[int] = None,
    zip_fn: Optional[PathLike] = None,
    gfn: Literal["ff", "0", "1", "2"] = "ff",
) -> CrestResult:
    inp = Path(fn)
    crest_cmd = get_cmd("crest")
    args = [crest_cmd]

    def add_arg(arg):
        args.extend(arg.split())

    add_arg(str(inp.absolute()))  # Input filename
    add_arg(f"-chrg {charge}")  # System's charge
    add_arg(f"-T {pal}")  # Number of threads

    # Solvent
    if alpb is not None:
        add_arg(f"-alpb {alpb}")
    # Clustering of conformers
    if cluster is not None:
        add_arg(f"-cluster {cluster}")
    # Level of theory
    gfn_arg = {
        "ff": "-gfnff",
        "0": "-gfn0",
        "1": "-gfn1",
        "2": "-gfn2",
    }[gfn]
    add_arg(gfn_arg)

    cmd = " ".join(args)
    # Run command like
    #   crest inp.xyz -gfnff -chrg 0 -alpb ch2cl2 -T 4 -cluster 5
    # in a temporary directory. Zip the whole directory if desired.
    with tempfile.TemporaryDirectory() as tmp:
        result = subprocess.run(
            cmd, shell=True, cwd=tmp, capture_output=True, text=True
        )
        if not result.returncode == 0:
            raise CrestException(result.stderr)
        if zip_fn is not None:
            zip_fn = Path(zip_fn)
            zip.zip_directory(tmp, zip_fn)
        crest_result = process_crest_run(tmp)
    return crest_result


def process_crest_run(crest_dir: PathLike) -> CrestResult:
    crest_path = Path(crest_dir)
    best_fn = crest_path / "crest_best.xyz"
    best = xyz.geom_from_xyz(best_fn)
    conformers_fn = crest_path / "crest_conformers.xyz"
    conformers = xyz.geoms_from_xyz(conformers_fn)
    energies_fn = crest_path / "crest.energies"
    # array is 1d when only 1 conformer is found; we force at least 2d
    energies_kcal = np.atleast_2d(np.loadtxt(energies_fn))
    energies = energies_kcal[:, 1] / AU2KCALPERMOL
    clustered_fn = crest_path / "crest_clustered.xyz"
    if clustered_fn.exists():
        clustered, clustered_comments = xyz.geoms_and_comments_from_xyz(clustered_fn)
        clustered_energies = np.array(clustered_comments, dtype=float)
        # To be consistent with the CrestResults.energies attribute we only save relative
        # energies. Energies are already in Hartree.
        clustered_energies -= clustered_energies.min()
    else:
        clustered = None
        clustered_energies = None

    result = CrestResult(
        best=best,
        conformers=conformers,
        energies=energies,
        clustered=clustered,
        clustered_energies=clustered_energies,
    )
    return result
