import re

import numpy as np

from pysisyphus.calculators.Calculator import Calculator


class DFTD3(Calculator):
    
    conf_key = "dftd3"

    def __init__(self, geom, functional, bjdamping=False, **kwargs):
        super().__init__(**kwargs)
        self.atoms = geom.atoms
        self.functional = functional.lower()
        self.bjdamping = bjdamping

        self.parser_funcs = {
            "energy": self.parse_energy,
            "grad": self.parse_gradient,
        }

    def parse_energy(self, path):
        with open(path / self.out_fn) as handle:
            text = handle.read()
        mobj = re.search(r"Edisp /kcal,au:\s+([\d\-\.]+)\s+([\d\-\.]+)", text)
        results={
            "energy": float(mobj.group(2)),
        }
        return results

    def parse_gradient(self, path):
        grad = np.loadtxt(path / "dftd3_gradient")
        energy = self.parse_energy(path)["energy"]
        results={
            "energy": energy,
            "grad": grad.flatten(),
        }
        return results
  
    def calc(self, coords3d, gradient=False):
        inp = self.prepare_turbo_coords(self.atoms, coords3d.flatten())
        
        args = ["-func", self.functional]
        if self.bjdamping:
            args.append("-bj")
        if gradient:
            args.append("-grad")
        
        if gradient:
            calc_type = "grad"
        else:
            calc_type = "energy"

        results = self.run(inp, calc=calc_type, cmd="dftd3", add_args=args)
        
        if not gradient:
            return results["energy"]
        else:
            return results["energy"], results["grad"] 
    
    def __str__(self):
        return f"DFT-D3({self.name})"



