#!/usr/bin/env python3

from calculators.Calculator import Calculator

inp="""BP86 def2-SV(P) def2-SV(P)/J R

class ORCA(Calculator):

    def __init__(self): 
        super(ORCA, self).__init__()

        self.cmd = "/home/carpx/programme/orca_3_0_3_linux_x86-64/orca"

    def get_energy(self, coords):
        raise Exception("Not implemented!")

    def get_grad(self, coords):
        raise Exception("Not implemented!")

    def get_hessian(self):
        raise Exception("Not implemented!")
