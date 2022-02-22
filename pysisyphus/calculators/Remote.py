from pathlib import Path
import tarfile

from fabric import Connection
import yaml

from pysisyphus.calculators.Calculator import Calculator
from pysisyphus.helpers_pure import json_to_results


class Remote(Calculator):
    def __init__(self, remote_calc, host, prefix="", **kwargs):
        super().__init__(**kwargs)

        self.remote_calc = remote_calc
        self.host = host
        self.prefix = prefix

        self.calc_inp = yaml.dump(self.remote_calc)
        self.tar_fn = "pysis.tar.gz"
        self.yaml_fn = "inp.yaml"
        self.run_dict = {
            "geom": {
                "fn": None,
            },
            "calc": {
                "run_func": None,
            },
        }
        self.run_dict["calc"].update(remote_calc)

    def run_calculation(self, atoms, coords, run_func="get_energy"):
        con = Connection(self.host)
        res = con.run("mktemp -d", hide=True)
        tmp_dir = res.stdout.strip()

        xyz_str = self.prepare_xyz_string(atoms, coords)
        run_dict = self.run_dict.copy()
        run_dict["geom"]["fn"] = xyz_str
        # Inline xyz coordinates
        run_dict["calc"]["run_func"] = run_func

        # Create YAML input for remote calculation
        with open(self.yaml_fn, "w") as handle:
            yaml.dump(run_dict, handle)

        with con.cd(tmp_dir):
            tmp_parent = Path(tmp_dir).parent
            tar_target = tmp_parent / self.tar_fn
            con.put(self.yaml_fn, tmp_dir)
            # yaml_str = yaml.dump(run_dict)
            # con.run(f"cat > {self.yaml_fn} <<EOL\n{yaml_str}\nEOL")
            # Execute pysis with prefix, e.g. activation of conda env or venv
            with con.prefix(self.prefix):
                con.run(f"pysis {self.yaml_fn}", hide=True)
            # Pack the whole directory
            con.run(f"tar -czf {tar_target} .")
            # and download it.
            con.get(tar_target, self.tar_fn)
        con.run(f"rm -r {tmp_dir}")
        return self.parse_results()

    def parse_results(self):
        with tarfile.open(self.tar_fn, "r:gz") as tfile:
            as_json = tfile.extractfile("./calculator_000.000.results").read()
            results = json_to_results(as_json)
        return results

    def get_energy(self, atoms, coords, **prepare_kwargs):
        return self.run_calculation(atoms, coords, run_func="get_energy")

    def get_forces(self, atoms, coords, **prepare_kwargs):
        return self.run_calculation(atoms, coords, run_func="get_forces")

    def get_hessian(self, atoms, coords, **prepare_kwargs):
        return self.run_calculation(atoms, coords, run_func="get_hessian")
