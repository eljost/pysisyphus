import logging

__all__ = [
    "BFGS",
    "ConjugateGradient",
    "FIRE",
    "QuickMin",
    "SteepestDescent",
    "SciPyOptimizer",
    "RFOptimizer",
]

logger = logging.getLogger("optimizer")
logger.setLevel(logging.DEBUG)
# delay = True prevents creation of empty logfiles
handler = logging.FileHandler("optimizer.log", mode="w", delay=True)
fmt_str = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
formatter = logging.Formatter(fmt_str)
handler.setFormatter(formatter)
logger.addHandler(handler)
