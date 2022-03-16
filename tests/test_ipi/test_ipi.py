from pathlib import Path
import subprocess
import tempfile
import threading
import time

import numpy as np
import pytest

from pysisyphus.calculators import XTB
from pysisyphus.calculators.IPIClient import calc_ipi_client
from pysisyphus.helpers import geom_loader
from pysisyphus.testing import using

SERVER_YAML = """
geom:
 type: redund
 fn: {fn}
calc:
 type: ipiserver
 address: {address}
 #verbose: True
opt:
 thresh: gau
 # Comment this out or set to False, if your client does not support GETHESSIAN
 do_hess: True
assert:
 opt_geom.energy: -5.07054442
"""


@pytest.mark.skip_ci
@using("xtb")
def test_ipi_server_client_opt():
    fn = "lib:h2o.xyz"
    address = "./opt_socket"

    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_path = Path(tmp_dir)

        def run_pysis():
            yaml = SERVER_YAML.format(fn=fn, address=tmp_path / address)
            print(yaml)
            yaml_fn = "inp.yaml"
            with open(tmp_path / yaml_fn, "w") as handle:
                handle.write(yaml)
            args = f"pysis {yaml_fn}".split()
            subprocess.run(args, cwd=tmp_path)

        pysis_thread = threading.Thread(target=run_pysis)
        pysis_thread.start()
        time.sleep(10)

        geom = geom_loader(fn)
        calc = XTB()
        calc_ipi_client(tmp_path / address, geom.atoms, calc)
        ref_coords = np.array(
            (
                -0.0,
                -0.2094041,
                0.0,
                1.48149972,
                0.8376164,
                0.0,
                -1.48149972,
                0.8376164,
                0.0,
            )
        )
        np.testing.assert_allclose(geom.coords, ref_coords)
