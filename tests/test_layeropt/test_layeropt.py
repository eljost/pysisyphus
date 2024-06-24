from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from queue import Empty, Queue
import subprocess
import tempfile
import time

import pytest

from pysisyphus.calculators import ONIOM
from pysisyphus.calculators.IPIClient import calc_ipi_client
from pysisyphus.helpers import geom_loader
from pysisyphus.init_logging import init_logging
from pysisyphus.testing import using


init_logging()


SERVER_YAML = """
geom:
 type: cart
 fn: {fn}
calc:
 type: dummy
opt:
 type: layers
 thresh: gau
 layers:
  - address: {address0}
    geom:
     type: redund
  - indices: {indices}
    address: {address1}
    geom:
     type: redund
assert:
 opt_geom.energy: -115.535533
"""


@pytest.mark.skip_ci
@using("xtb")
@using("orca")
def test_ethanal_oniom_layeropt(this_dir):
    fn = "lib:acetaldehyd_oniom.xyz"
    indices = [4, 5, 6]  # 0 is link atom host
    link_hosts = [
        0,
    ]

    geom = geom_loader(fn)
    oniom_kwargs = {
        "calcs": {
            "real": {
                "type": "xtb",
                "quiet": True,
            },
            "high": {
                "type": "orca",
                "keywords": "hf sto-3g norijcosx",
            },
        },
        "models": {
            "high": {
                "inds": indices,
                "calc": "high",
            },
        },
    }
    calc = ONIOM(**oniom_kwargs, geom=geom)
    geom.set_calculator(calc)
    (real_model,), (high_model,) = calc.layers
    real_geom = real_model.as_geom(geom.atoms, geom.coords3d.copy())

    addrs = [f"./layer{i}sock" for i, _ in enumerate(oniom_kwargs["calcs"])]
    address0, address1 = addrs

    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_path = Path(tmp_dir)

        def run_pysis():
            yaml = SERVER_YAML.format(
                fn=fn,
                address0=address0,
                address1=address1,
                indices=link_hosts + indices,
            )
            yaml_fn = "inp.yaml"
            with open(tmp_path / yaml_fn, "w") as handle:
                handle.write(yaml)
            args = f"pysis {yaml_fn}".split()
            subprocess.run(args, cwd=tmp_path)

        def layer0sock():
            time.sleep(3)
            addr = str(tmp_path / address0)
            print(f"Starting {addr}")
            calc_ipi_client(addr, real_geom.atoms, real_geom.calculator)

        coords_queue = Queue(maxsize=0)

        def layer1sock():
            time.sleep(3)
            addr = str(tmp_path / address1)
            print(f"Starting {addr}")
            calc_ipi_client(addr, geom.atoms, geom.calculator, queue=coords_queue)

        with ThreadPoolExecutor() as executor:
            _ = [executor.submit(task) for task in (run_pysis, layer0sock, layer1sock)]
        all_cart_coords = list()
        while True:
            try:
                _, cart_coords = coords_queue.get(False)
                all_cart_coords.append(cart_coords)
            except Empty:
                break
        xyzs = [geom.as_xyz(cart_coords=cc) for cc in all_cart_coords]
        with open(this_dir / "optimization.trj", "w") as handle:
            handle.write("\n".join(xyzs))


if __name__ == "__main__":
    test_ethanal_oniom_layeropt()
