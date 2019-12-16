import logging

__all__ = [
    "DampedVelocityVerlet",
    "Euler",
    "EulerPC",
    "LQA",
    "GonzalesSchlegel",
    "IMKMod",
    "ModeKill",
    "RK4",
]

from pysisyphus.irc.DampedVelocityVerlet import DampedVelocityVerlet
from pysisyphus.irc.Euler import Euler
from pysisyphus.irc.EulerPC import EulerPC
from pysisyphus.irc.GonzalesSchlegel import GonzalesSchlegel
from pysisyphus.irc.IMKMod import IMKMod
from pysisyphus.irc.LQA import LQA
from pysisyphus.irc.ModeKill import ModeKill
from pysisyphus.irc.RK4 import RK4

logger = logging.getLogger("irc")
logger.setLevel(logging.DEBUG)
# delay = True prevents creation of empty logfiles
handler = logging.FileHandler("irc.log", mode="w", delay=True)
# fmt_str = "%(levelname)s - %(message)s"
fmt_str = "%(message)s"
formatter = logging.Formatter(fmt_str)
handler.setFormatter(formatter)
logger.addHandler(handler)
