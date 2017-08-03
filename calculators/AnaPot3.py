from pysisyphus.calculators.AnaPotBase import AnaPotBase

# [1] http://www.cims.nyu.edu/~eve2/string_jcp_simplified07.pdf

class AnaPot3(AnaPotBase):

    def __init__(self): 
        V_str = "(1 - x**2 - y**2)**2 + (y**2) / (x**2 + y**2)"
        super(AnaPot3, self).__init__(V_str=V_str)

    def __str__(self):
        return "AnaPot3 calculator"
