#!/usr/bin/env python3

import logging
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

logging.basicConfig()
logging.getLogger().setLevel(logging.INFO)

class Calculator:

    def __init__(self):

        self._energy = None
        self._forces = None
        self._hessian = None

        self.inp_fn = "calc.inp"
        self.out_fn = "calc.out"

    def get_energy(self, coords):
        raise Exception("Not implemented!")

    def get_hessian(self):
        raise Exception("Not implemented!")

    def prepare(self, inp, path=None):
        if not path:
            path = tempfile.mkdtemp()
        inp_path = os.path.join(path, self.inp_fn)
        with open(inp_path, "w") as handle:
            handle.write(inp)

        return path

    def run(self, inp, calc, add_args=None):
        path = self.prepare(inp)
        args = [self.base_cmd, self.inp_fn]
        if add_args:
            args.extend(add_args)
        logging.debug("Running calculation in {}".format(path))
        with open(os.path.join(path, self.out_fn), "w") as handle:
            result = subprocess.Popen(args, cwd=path, stdout=handle)
            result.wait()
        #logging.info("Calculation finished".format(path))
        try:
            results = self.parser_funcs[calc](path)
        except Exception as err:
            print(err)
            print()
            print("Crashed input:")
            print(inp)
            backup_dir = Path(os.getcwd()) / "crashed"
            if backup_dir.exists():
                shutil.rmtree(backup_dir)
            shutil.copytree(path, backup_dir)
            sys.exit()
        finally:
            self.clean(path)

        return results

    def clean(self, path):
        shutil.rmtree(path)
        #logging.info("Removed {}".format(path))
