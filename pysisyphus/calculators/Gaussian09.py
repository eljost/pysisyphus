#!/usr/bin/env python3

from pysisyphus.calculators.Gaussian16 import Gaussian16


class Gaussian09(Gaussian16):

    conf_key = "gaussian09"

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.fn_base = "gaussian09"
        self.inp_fn = f"{self.fn_base}.com"
        self.out_fn = f"{self.fn_base}.log"
        self.chk_fn = f"{self.fn_base}.chk"
        self.base_cmd = self.get_cmd("cmd")
        self.formchk_cmd = self.get_cmd("formchk_cmd")

    def __str__(self):
        return "Gaussian09 calculator"


if __name__ == "__main__":
    from pysisyphus.helpers import geom_from_library
    from pysisyphus.init_logging import init_logging
    init_logging()
    geom = geom_from_library("h2.xyz")
    kwargs = {
        "method": "HF",
        "basis": "STO-3G",
    }
    g09 = Gaussian09(**kwargs)
    geom.set_calculator(g09)
    #import pdb; pdb.set_trace()
    forces = geom.forces
    print(forces)
