from pysisyphus.calculators.Calculator import Calculator


class EnergyMin(Calculator):

    def __init__(self, calculator1, calculator2, **kwargs):
        super().__init__(**kwargs)
        self.calc1 = calculator1
        self.calc2 = calculator2

    def do_calculations(self, name, atoms, coords, **prepare_kwargs):
        results1 = getattr(self.calc1, name)(atoms, coords, **prepare_kwargs)
        results2 = getattr(self.calc2, name)(atoms, coords, **prepare_kwargs)

        energy1 = results1["energy"]
        energy2 = results2["energy"]
        en1_or_en2 = "calc1" if energy1 < energy2 else "calc2"
        self.log(
            f"energy_calc1={energy1:.6f} au, energy_calc2={energy2:.6f} au, returning "
            f"results for {en1_or_en2}."
        )

        return results1 if (energy1 < energy2) else results2

    def get_energy(self, atoms, coords, **prepare_kwargs):
        return self.do_calculations("get_energy", atoms, coords, **prepare_kwargs)

    def get_forces(self, atoms, coords, **prepare_kwargs):
        return self.do_calculations("get_forces", atoms, coords, **prepare_kwargs)

    def get_hessian(self, atoms, coords, **prepare_kwargs):
        return self.do_calculations("get_hessian", atoms, coords, **prepare_kwargs)
