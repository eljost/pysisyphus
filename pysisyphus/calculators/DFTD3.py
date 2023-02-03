from pysisyphus.calculators.Calculator import Calculator
import re
import numpy as np
import os.path

class DFTD3(Calculator):
    
    conf_key = "dftd3"

    def __init__(self, geom, functional, bjdamping=False, **kwargs):
        super(DFTD3, self).__init__(**kwargs)
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
            "energy": float(mobj[-1]),
        }
        return results

    def parse_gradient(self, path):
        grad = np.loadtxt(os.path.join(path,"dftd3_gradient"))
        energy = self.parse_energy(path)
        results={
            "energy": energy,
            "grad": grad.flatten(),
        }
        return results
  
    def calc(self, coords3d, gradient=False):
        inp = self.prepare_turbo_coords(self.atoms, coords3d.flatten())
        
        args = [f"-func {self.functional}"]
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



