import logging

__all__ = [
    "ORCA",
    "XTB",
    "IDPP",
    "OpenMolcas",
    "Gaussian09",
    "Turbomole",
]

logger = logging.getLogger("calculator")
logger.setLevel(logging.DEBUG)
# delay = True prevents creation of empty logfiles
handler = logging.FileHandler("calculator.log", mode="w")
fmt_str = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
formatter = logging.Formatter(fmt_str)
handler.setFormatter(formatter)
logger.addHandler(handler)
