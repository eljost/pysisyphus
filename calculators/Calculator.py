#!/usr/bin/env python3

class Calculator:

    def __init__(self):

        self._energy = None
        self._forces = None
        self._hessian = None

    def get_energy(self, coords):
        raise Exception("Not implemented!")

    def get_grad(self, coords):
        raise Exception("Not implemented!")

    def get_hessian(self):
        raise Exception("Not implemented!")
