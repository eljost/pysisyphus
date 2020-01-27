from pysisyphus.calculators.AnaPotBase2D import AnaPotBase2D

# [1] http://www.cims.nyu.edu/~eve2/string_jcp_simplified07.pdf

class MullerBrownSympyPot2D(AnaPotBase2D):

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
                                        
        super(MullerBrownSympyPot2D, self).__init__(V_str=V_str)

    def __str__(self):
        return "MullerBrownSympyPot2D calculator"
