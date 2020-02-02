import logging

__all__ = [
    "ConjugateGradient",
    "FIRE",
    "LBFGS",
    #"LBFGS_mod",
    "QuickMin",
    "RFOptimizer",
    "SteepestDescent",
    "ANCOptimizer",
    "StringOptimizer",
    "StabilizedQNMethod",
]

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
