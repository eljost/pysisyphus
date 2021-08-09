import logging

from pysisyphus.tsoptimizers.RSPRFOptimizer import RSPRFOptimizer
from pysisyphus.tsoptimizers.TRIM import TRIM
from pysisyphus.tsoptimizers.RSIRFOptimizer import RSIRFOptimizer
from pysisyphus.tsoptimizers.dimer import dimer_method

__all__ = [
    "RSPRFOptimizer",
    "TRIM",
    "RSIRFOptimizer",
    "TSHessianOptimizer",
    "dimer_method",
]


logger = logging.getLogger("tsoptimizer")
logger.setLevel(logging.DEBUG)
# delay = True prevents creation of empty logfiles
handler = logging.FileHandler("tsoptimizer.log", mode="w", delay=True)
# fmt_str = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
fmt_str = "%(message)s"
formatter = logging.Formatter(fmt_str)
handler.setFormatter(formatter)
logger.addHandler(handler)
