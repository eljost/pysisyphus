from pysisyphus.calculators.AnaPotBase2D import AnaPotBase2D

class AnaPot2D(AnaPotBase2D):

    def __init__(self): 
        V_str = "4 + 4.5*x - 4*y + x**2 + 2*y**2-2*x*y + x**4 - 2*x**2*y"
        super(AnaPot2D, self).__init__(V_str=V_str)

    def __str__(self):
        return "AnaPot2D calculator"
