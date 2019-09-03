import logging

__all__ = [
    "DampedVelocityVerlet",
    "Euler",
    "GonzalesSchlegel",
    "IMKMod",
    "RK4",
    "LQA",
    "ModeKill",
]

from pysisyphus.irc.RK4 import RK4
from pysisyphus.irc.LQA import LQA
from pysisyphus.irc.ModeKill import ModeKill

logger = logging.getLogger("irc")
logger.setLevel(logging.DEBUG)
# delay = True prevents creation of empty logfiles
handler = logging.FileHandler("irc.log", mode="w", delay=True)
fmt_str = "%(levelname)s - %(message)s"
formatter = logging.Formatter(fmt_str)
handler.setFormatter(formatter)
logger.addHandler(handler)
