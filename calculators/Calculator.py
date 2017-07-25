#!/usr/bin/env python3

import os
import shutil
import subprocess
import tempfile

class Calculator:

    def __init__(self):

        self._energy = None
        self._forces = None
        self._hessian = None

        self.inp_fn = "calc.inp"
        self.out_fn = "calc.out"

    def get_energy(self, coords):
        raise Exception("Not implemented!")

    def get_grad(self, coords):
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

    def run(self, inp):
        path = self.prepare(inp)
        args = [self.base_cmd, self.inp_fn]
        with open(os.path.join(path, self.out_fn), "w") as handle:
            result = subprocess.Popen(args, cwd=path, stdout=handle)
            result.wait()
        results = self.parse(path)
        self.clean(path)
        return results

    def clean(self, path):
        shutil.rmtree(path)
