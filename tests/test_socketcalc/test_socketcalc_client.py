import json
import socket

import numpy as np
import pytest

from pysisyphus.calculators.XTB import XTB


@pytest.mark.skip
def test_socketcalc():
    host = "localhost"
    port = 8080

    calc = XTB(pal=2)
    funcs = {
        "energy": calc.get_energy,
        "forces": calc.get_forces,
        "hessian": calc.get_hessian,
    }

    with socket.socket() as s:
        s.connect((host, port))
        data = s.recv(1024)

        try:
            js = json.loads(data.decode("utf-8"))
            print("Decoded json", js)
        except json.JSONDecodeError:
            js = None
            print("JSONDecode error")

        if js:
            print("Got valid JSON")
            atoms = js["atoms"]
            print("atoms", atoms)
            coords = np.array(js["coords"], dtype=float)
            print("coords", coords)

            result = funcs[js["request"]](atoms, coords)
            if "forces" in result:
                result["forces"] = result["forces"].tolist()
            if "hessian" in result:
                result["hessian"] = result["hessian"].tolist()
            print("result", result)

            response = (json.dumps(result) + "\n").encode("utf-8")
        s.sendall(response)


if __name__ == "__main__":
    test_socketcalc()
