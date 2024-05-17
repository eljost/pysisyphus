import numpy as np
import sympy as sym

# from sympy import sympify, lambdify,

from pysisyphus.calculators.Calculator import Calculator

from pysisyphus.calculators import ORCA, Turbomole, DFTD4


CALC_CLASSES = {
    "dftd4": DFTD4.DFTD4,
    "orca": ORCA.ORCA,
    "turbomole": Turbomole.Turbomole,
}
try:
    from pysisyphus.calculators import PySCF

    CALC_CLASSES["pyscf"] = PySCF.PySCF
except (ModuleNotFoundError, OSError):
    pass


class Composite(Calculator):
    def __init__(
        self, final, keys_calcs=None, calcs=None, remove_translation=False, **kwargs
    ):
        # Either directly supply a dictionary with Calculator objects (key_calcs)
        # or a dictionary containing information to set up calculators.
        assert keys_calcs or calcs

        super().__init__(**kwargs)

        if calcs:
            keys_calcs = {}
            calc_kwargs_ = {
                "charge": self.charge,
                "mult": self.mult,
                "pal": self.pal,
                "mem": self.mem,
            }
            for i, (key, kwargs) in enumerate(calcs.items()):
                calc_kwargs = calc_kwargs_.copy()
                calc_kwargs["calc_number"] = i
                calc_kwargs["base_name"] = f"composite_{key}"
                # Don't modify original dict
                kwargs = kwargs.copy()
                type_ = kwargs.pop("type")
                calc_kwargs.update(**kwargs)
                calc_cls = CALC_CLASSES[type_]
                keys_calcs[key] = calc_cls(**calc_kwargs)

        self.keys_calcs = keys_calcs
        assert all([key in final for key in self.keys_calcs.keys()])
        self.final = final
        self.remove_translation = remove_translation

        # The energies are just numbers that we can easily substitute in
        expr = sym.sympify(self.final)
        free_symbs = expr.free_symbols
        symb_names = {symb.name for symb in free_symbs}
        keys = list(self.keys_calcs.keys())
        if not symb_names == set(keys):
            unknown_symbs = free_symbs - keys
            raise Exception(f"Found unknown symbol(s): {unknown_symbs} in 'final'!")

        x = sym.symbols("x", real=True)
        key_funcs = {symb: sym.Function(symb.name)(x) for symb in free_symbs}

        # Map between key string and sympy symbol
        key_symbs = {}
        for key in keys:
            for fsymb in free_symbs:
                if fsymb.name == key:
                    key_symbs[key] = fsymb

        expr_x = expr.subs(key_funcs)
        # Take derivatives w.r.t. coordinates
        dexpr_x = sym.diff(expr_x, x)  # Gradient
        d2expr_x = sym.diff(expr_x, x, x)  # Hessian

        """
        Energies will always be present and gradients/derivatives never appear in the energy
        expression. The gradient expression may depend on the energy, but it will always be
        available when a gradient is calculated.
        Special care is be needed when dealing with the Hessian, as it may also depend on the
        gradient, that is not necessarily calculated when calculating a Hessian. So we create
        a list of first derivatives that appear in the Hessian (derivs_in_d2).
        """
        derivs = {}
        derivs2 = {}
        derivs_in_d2 = {}
        # Inverted dict 'key_funcs'
        func_subs = {value: key for key, value in key_funcs.items()}
        deriv_subs = {}
        deriv2_subs = {}
        deriv_args = list()
        deriv2_deriv_args = list()
        deriv2_args = list()
        for key in keys:
            symb = key_symbs[key]
            name = symb.name
            func = key_funcs[symb]
            deriv = sym.Derivative(func, x)
            deriv2 = sym.Derivative(func, (x, 2))
            # Squared derivatives appearing in the expression actually correspond to Hessians
            d2expr_x = d2expr_x.subs(deriv**2, deriv2)
            derivs[symb] = deriv
            derivs2[symb] = deriv2
            deriv_symb = sym.symbols(f"{name}_deriv", real=True)
            deriv_args.append(deriv_symb.name)
            deriv2_symb = sym.symbols(f"{name}_deriv2", real=True)
            deriv2_args.append(deriv2_symb.name)

            # Check if gradient appears in Hessian expressions
            if has_deriv := d2expr_x.has(deriv):
                derivs_in_d2[symb] = has_deriv
                deriv2_deriv_args.append(deriv_symb.name)
            deriv_subs[deriv] = deriv_symb
            deriv2_subs[deriv2] = deriv2_symb
        dexpr = dexpr_x.subs(deriv_subs).subs(func_subs)
        d2expr = d2expr_x.subs(deriv2_subs).subs(deriv_subs).subs(func_subs)
        assert (
            len(deriv2_deriv_args) == 0
        ), "Hessian expressions that depend on the first derivative are not yet supported!"

        deriv_args = keys + deriv_args
        deriv2_args = keys + deriv2_deriv_args + deriv2_args

        # Setup function that will be used to evaluate the energy and its derivatives
        # from the different calculators.
        self.energy_func = sym.lambdify(keys, expr)
        self.grad_func = sym.lambdify(deriv_args, dexpr)
        self.hessian_func = sym.lambdify(deriv2_args, d2expr)

    def get_energy(self, atoms, coords, **prepare_kwargs):
        energy_kwargs = {}
        for key, calc in self.keys_calcs.items():
            energy = calc.get_energy(atoms, coords, **prepare_kwargs)["energy"]
            energy_kwargs[key] = energy

        final_energy = self.energy_func(**energy_kwargs)
        results = {
            "energy": final_energy,
        }
        return results

    def get_forces(self, atoms, coords, **prepare_kwargs):
        energy_kwargs = {}
        deriv_kwargs = dict()
        for key, calc in self.keys_calcs.items():
            results = calc.get_forces(atoms, coords, **prepare_kwargs)
            energy_kwargs[key] = results["energy"]
            deriv_kwargs[key] = results["energy"]
            deriv_kwargs[key + "_deriv"] = -results["forces"]
        keys = self.keys_calcs.keys()
        for key in keys:
            self.log(
                f"|forces_{key}|={np.linalg.norm(deriv_kwargs[key + "_deriv"]):.6f}"
            )
        self.log("")

        final_energy = self.energy_func(**energy_kwargs)
        final_forces = -self.grad_func(**deriv_kwargs)

        # Remove overall translation
        if self.remove_translation:
            f3d = final_forces.reshape(-1, 3)
            f3d -= f3d.mean(axis=0)[None, :]

        results = {
            "energy": final_energy,
            "forces": final_forces,
        }
        return results

    def get_hessian(self, atoms, coords, **prepare_kwargs):
        energy_kwargs = {}
        deriv2_kwargs = {}
        for key, calc in self.keys_calcs.items():
            results = calc.get_hessian(atoms, coords, **prepare_kwargs)
            energy_kwargs[key] = results["energy"]
            deriv2_kwargs[key] = results["energy"]
            deriv2_kwargs[key + "_deriv2"] = results["hessian"]

        final_energy = self.energy_func(**energy_kwargs)
        final_hessian = self.hessian_func(**deriv2_kwargs)

        results = {
            "energy": final_energy,
            "hessian": final_hessian,
        }
        return results

    def run_calculation(self, atoms, coords, **prepare_kwargs):
        return self.get_energy(atoms, coords, **prepare_kwargs)
