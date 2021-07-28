import re
import textwrap

import numpy as np

from pysisyphus.calculators.Calculator import Calculator


class Psi4(Calculator):

    conf_key = "psi4"

    def __init__(self, method, basis, to_set=None, pcm="iefpcm",
                 solvent=None, mem=2000, write_fchk=False, **kwargs):
        super().__init__(**kwargs)

        self.method = method
        self.basis = basis
        self.to_set = {} if to_set is None else dict(to_set)
        self.pcm = pcm
        self.solvent = solvent
        self.write_fchk = write_fchk
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

        self.inp = textwrap.dedent("""
        molecule mol{{
          {xyz}
          {charge} {mult}
        symmetry c1
        }}

        set_num_threads({pal})
        memory {mem} MB

        {basis}
        {to_set}
        {pcm}

        {method}

        {fchk}
        """)

    def get_fchk_str(self):
        fchk_str = ""
        if self.write_fchk:
            fchk_fn = self.make_fn("wfn.fchk")
            fchk_str = (
                "fchk_writer = psi4.FCHKWriter(wfn)\n"
                f"fchk_writer.write('{fchk_fn}')"
            )
        return fchk_str

    def prepare_input(self, atoms, coords, calc_type):
        xyz = self.prepare_coords(atoms, coords)

        calc_types = {
            "energy": "E, wfn = energy('{}', return_wfn=True)",
            # Right now we don't need the wavefunction
            # "grad": "G, wfn = gradient('{}', return_wfn=True)\n" \
                    # "G_arr = np.array(G)\n" \
                    # "np.save('grad', G_arr)",
            # "hessian": "H, wfn = hessian('{}', return_wfn=True)\n" \
                       # "H_arr = np.array(H)\n" \
                       # "np.save('hessian', H_arr)",
            "grad": "G, wfn = gradient('{}', return_wfn=True)\n" \
                    "G_arr = np.array(G)\n" \
                    "np.save('grad', G_arr)",
            "hessian": "H, wfn = hessian('{}', return_wfn=True)\n" \
                       "H_arr = np.array(H)\n" \
                       "np.save('hessian', H_arr)",
        }
        method = calc_types[calc_type].format(self.method)
        wfn_path = self.make_fn("wfn.npy")
        method += f"\nWavefunction.to_file(wfn, '{wfn_path}')"
        method += "\nprint('PARSE ENERGY:', wfn.energy())"
        set_strs = [f"set {key} {value}" for key, value in self.to_set.items()]
        set_strs = "\n".join(set_strs)

        # Basis section
        basis = self.basis
        # Construct more complex basis input
        if isinstance(basis, dict):
            # Check if a global basis is given for all atoms. This must come
            # first, otherwise Psi4 throws an error.
            basis_lines = ["basis {", ]
            try:
                basis_lines.append(f"assign {basis['assign']}")
            except KeyError:
                pass
            # Add remaining lines
            basis_lines.extend(
                [f"assign {atms} {bas}" for atms, bas in basis.items()
                 if atms != "assign"]
            )
            basis_lines.append("}")

            basis = "\n".join(basis_lines)
        # Use set when self.basis is a string
        else:
            basis =  f"set basis {basis}"

        # PCM section
        pcm = ""
        if self.solvent:
            pcm = textwrap.dedent(f"""
            set pcm true

            pcm = {{
                Medium {{
                    SolverType = {self.pcm}
                    Solvent = {self.solvent}
                }}

                Cavity {{
                    Type = GePol
                }}
            }}
            """)

        inp = self.inp.format(
                xyz=xyz,
                charge=self.charge,
                mult=self.mult,
                basis=basis,
                to_set=set_strs,
                pcm=pcm,
                method=method,
                pal=self.pal,
                mem=self.mem,
                fchk=self.get_fchk_str(),
        )
        # inp = "\n".join([line.strip() for line in inp.split("\n")])
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

    def run_calculation(self, atoms, coords):
        return self.get_energy(atoms, coords)

    def parse_energy(self, path):
        with open(path / "psi4.out") as handle:
            text = handle.read()
        en_regex = re.compile("PARSE ENERGY: ([\d\-\.]+)")
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

    def __str__(self):
        return f"Psi4({self.name})"
