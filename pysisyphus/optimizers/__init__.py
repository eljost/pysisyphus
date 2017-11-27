import logging

__all__ = [
    "BFGS",
    "ConjugateGradient",
    "FIRE",
    "QuickMin",
    "SteepestDescent",
    "SciPyOptimizer",
]

logger = logging.getLogger("optimizer")
logger.setLevel(logging.DEBUG)
handler = logging.FileHandler("optimizer.log", mode="w")
fmt_str = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
formatter = logging.Formatter(fmt_str)
handler.setFormatter(formatter)
logger.addHandler(handler)
