#!/usr/bin/env python3

import re

import numpy as np

from pysisyphus.calculators.Calculator import Calculator


class Psi4(Calculator):

    conf_key = "psi4"

    def __init__(self, method, basis, to_set=None, mem=2000, **kwargs):
        super().__init__(**kwargs)

        self.method = method
        self.basis = basis
        self.to_set = {} if to_set is None else dict(to_set)
        self.mem = mem

        self.inp_fn = "psi4.inp"
        self.out_fn = "psi4.out"
        self.to_keep = ("inp", "psi4.out", "grad.npy", "hessian.npy")

        self.parser_funcs = {
            "energy": self.parse_energy,
            "grad": self.parse_grad,
            "hessian": self.parse_hessian,
        }

        self.base_cmd = self.get_cmd("cmd")

        self.inp = """
        molecule mol{{
          {xyz}
          {charge} {mult}
        symmetry c1
        }}

        set_num_threads({pal})
        memory {mem} MB


        set basis {basis}
        {to_set}

        {method}
        """

    def prepare_input(self, atoms, coords, calc_type):
        xyz = self.prepare_coords(atoms, coords)

        calc_types = {
            "energy": "energy('{}')",
            # Right now we don't need the wavefunction
            # "grad": "G, wfn = gradient('{}', return_wfn=True)\n" \
                    # "G_arr = np.array(G)\n" \
                    # "np.save('grad', G_arr)",
            # "hessian": "H, wfn = hessian('{}', return_wfn=True)\n" \
                       # "H_arr = np.array(H)\n" \
                       # "np.save('hessian', H_arr)",
            "grad": "G = gradient('{}')\n" \
                    "G_arr = np.array(G)\n" \
                    "np.save('grad', G_arr)",
            "hessian": "H = hessian('{}')\n" \
                       "H_arr = np.array(H)\n" \
                       "np.save('hessian', H_arr)",
        }
        method = calc_types[calc_type].format(self.method)
        set_strs = [f"set {key} {value}" for key, value in self.to_set.items()]
        set_strs = "\n".join(set_strs)

        inp = self.inp.format(
                xyz=xyz,
                charge=self.charge,
                mult=self.mult,
                basis=self.basis,
                to_set=set_strs,
                method=method,
                pal=self.pal,
                mem=self.mem,
        )
        inp = "\n".join([line.strip() for line in inp.split("\n")])
        return inp

    def get_energy(self, atoms, coords):
        calc_type = "energy"
        inp = self.prepare_input(atoms, coords, calc_type)
        results = self.run(inp, calc="energy")
        return results

    def get_forces(self, atoms, coords):
        calc_type = "grad"
        inp = self.prepare_input(atoms, coords, calc_type)
        results = self.run(inp, calc="grad")
        return results

    def get_hessian(self, atoms, coords):
        calc_type = "hessian"
        inp = self.prepare_input(atoms, coords, calc_type)
        results = self.run(inp, calc="hessian")
        return results

    def parse_energy(self, path):
        with open(path / "psi4.out") as handle:
            text = handle.read()
        en_regex = re.compile("Total Energy =\s*([\d\-\.]+)")
        mobj = en_regex.search(text)
        result = {
            "energy": float(mobj[1])
        }
        return result


    def parse_grad(self, path):
        gradient = np.load(path / "grad.npy")
        forces = -gradient.flatten()
        result = {
            "forces": forces,
        }
        result.update(self.parse_energy(path))
        return result

    def parse_hessian(self, path):
        hessian = np.load(path / "hessian.npy")
        result = {
            "hessian": hessian,
        }
        result.update(self.parse_energy(path))
        return result

    def keep(self, path):
        kept_fns = super().keep(path)

    def __str__(self):
        return f"Psi4({self.name})"
