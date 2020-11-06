from pysisyphus.calculators.OverlapCalculator import OverlapCalculator


class DFTBp(OverlapCalculator):

    conf_key = "dftpb"

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.parser_funcs = {
            "grad": self.parse_forces,
        }

        self.base_cmd = self.get_cmd("cmd")

    def get_forces(self, atoms, coords, prepare_kwargs=None):
        pass

    def parse_forces(self, path):
        pass
