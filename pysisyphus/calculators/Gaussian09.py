from pysisyphus.calculators.Gaussian16 import Gaussian16


class Gaussian09(Gaussian16):

    conf_key = "gaussian09"

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.fn_base = "gaussian09"
        self.inp_fn = f"{self.fn_base}.com"
        self.out_fn = f"{self.fn_base}.log"
        self.chk_fn = f"{self.fn_base}.chk"
        self.base_cmd = self.get_cmd()
        self.formchk_cmd = self.get_cmd("formchk")

    def __str__(self):
        return "Gaussian09 calculator"
