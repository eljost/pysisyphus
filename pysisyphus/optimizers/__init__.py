import logging

__all__ = [
    "BFGS",
    "ConjugateGradient",
    "CubicNewton",
    "FIRE",
    "LayerOpt",
    "LBFGS",
    "MicroOptimizer",
    "NCOptimizer",
    "PreconLBFGS",
    "PreconSteepestDescent",
    "QuickMin",
    "RFOptimizer",
    "SteepestDescent",
    "StringOptimizer",
    "StabilizedQNMethod",
]

from pysisyphus.optimizers.CubicNewton import CubicNewton
from pysisyphus.optimizers.MicroOptimizer import MicroOptimizer

logger = logging.getLogger("optimizer")
logger.setLevel(logging.DEBUG)
# delay = True prevents creation of empty logfiles
handler = logging.FileHandler("optimizer.log", mode="w", delay=True)
# fmt_str = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
fmt_str = "%(message)s"
formatter = logging.Formatter(fmt_str)
handler.setFormatter(formatter)
logger.addHandler(handler)

# logger = logging.getLogger("gdiis")
# logger.setLevel(logging.DEBUG)
# # delay = True prevents creation of empty logfiles
# handler = logging.FileHandler("gdiis.log", mode="w", delay=True)
# fmt_str = "%(message)s"
# formatter = logging.Formatter(fmt_str)
# handler.setFormatter(formatter)
# logger.addHandler(handler)
