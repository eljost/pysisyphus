#!/usr/bin/env/python3

# [1] http://www.cims.nyu.edu/~eve2/string_jcp_simplified07.pdf

import numpy as np

from pysisyphus.calculators.AnaPotBase import AnaPotBase


class MullerBrownPot(AnaPotBase):

    def __init__(self): 
        A  = (-200, -100, -170, 15)
        x0 = (1.0, 0.0, -0.5, -1.0)
        y0 = (0.0, 0.5, 1.5, 1.0)
        a  = (-1.0, -1.0, -6.5, 0.7)
        b  = (0.0, 0.0, 11.0, 0.6)
        c  = (-10.0, -10.0, -6.5, 0.7)

        #V_str_base = """{Ai} * exp(
        #                        {ai}*(x-{xi})**2
        #                        + {bi}*(x-{xi})*(y-{yi})
        #                        + {ci}*(y-{yi})**2
        #)"""
        V_str_base = "{Ai}*exp({ai}*(x-{xi})**2 + {bi}*(x-{xi})*(y-{yi}) + {ci}*(y-{yi})**2)"
        V_str = ""
        V_strs = [V_str_base.format(
                                Ai=A[i],
                                ai=a[i],
                                xi=x0[i],
                                bi=b[i],
                                yi=y0[i],
                                ci=c[i])
                  for i in range(4)
        ]
        V_str = " + ".join(V_strs)
        xlim = (-1.75, 1.25)
        ylim = (-0.5, 2.25)
        levels = np.linspace(-200, 400, 100)

        super(MullerBrownPot, self).__init__(V_str=V_str, xlim=xlim, ylim=ylim,
                                             levels=levels)

    def __str__(self):
        return "MullerBrownPot calculator"
