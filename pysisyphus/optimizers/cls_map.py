from pysisyphus.optimizers import (
    BFGS,
    ConjugateGradient,
    FIRE,
    LayerOpt,
    LBFGS,
    MicroOptimizer,
    NCOptimizer,
    PreconLBFGS,
    PreconSteepestDescent,
    QuickMin,
    RFOptimizer,
    SteepestDescent,
    StabilizedQNMethod,
    StringOptimizer,
)
from pysisyphus.tsoptimizers import (
    RSPRFOptimizer,
    TRIM,
    RSIRFOptimizer,
)


OPT_DICT = {
    "bfgs": BFGS.BFGS,
    "cg": ConjugateGradient.ConjugateGradient,
    "fire": FIRE.FIRE,
    "layers": LayerOpt.LayerOpt,
    "lbfgs": LBFGS.LBFGS,
    "micro": MicroOptimizer,
    "nc": NCOptimizer.NCOptimizer,
    "plbfgs": PreconLBFGS.PreconLBFGS,
    "psd": PreconSteepestDescent.PreconSteepestDescent,
    "qm": QuickMin.QuickMin,
    "rfo": RFOptimizer.RFOptimizer,
    "sd": SteepestDescent.SteepestDescent,
    "sqnm": StabilizedQNMethod.StabilizedQNMethod,
    "string": StringOptimizer.StringOptimizer,
}

TSOPT_DICT = {
    "rsprfo": RSPRFOptimizer,
    "prfo": RSPRFOptimizer,  # Shortcut
    "trim": TRIM,
    "rsirfo": RSIRFOptimizer,
    "irfo": RSIRFOptimizer,  # Shortcut
}
OPT_DICT.update(TSOPT_DICT)


def get_opt_cls(opt_key):
    return OPT_DICT[opt_key]


def key_is_tsopt(opt_key):
    return opt_key in TSOPT_DICT
