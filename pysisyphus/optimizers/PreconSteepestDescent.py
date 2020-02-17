from pysisyphus.optimizers.PreconLBFGS import PreconLBFGS


class PreconSteepestDescent(PreconLBFGS):

    def __init__(self, geometry, alpha_init=0.5, **kwargs):
        super().__init__(geometry, alpha_init=alpha_init, history=0, **kwargs)
