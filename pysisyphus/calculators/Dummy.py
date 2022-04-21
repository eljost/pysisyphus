from pysisyphus.calculators.Calculator import Calculator


class Dummy(Calculator):
    def raise_exception(self):
        raise Exception(
            "Dummy calculator does not implement any calculation capabilities."
        )

    def run_calculation(self, *args, **kwargs):
        self.raise_exception()

    def get_energy(self, *args, **kwargs):
        self.raise_exception()

    def get_forces(self, *args, **kwargs):
        self.raise_exception()

    def get_hessian(self, *args, **kwargs):
        self.raise_exception()
