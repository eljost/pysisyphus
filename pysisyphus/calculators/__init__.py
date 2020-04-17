import logging

__all__ = [
    "AFIR",
    "Dimer",
    "ExternalPotential",
    "FakeASE",
    "Gaussian09",
    "Gaussian16",
    "LennardJones",
    "MOPAC",
    "OpenMolcas",
    "ORCA",
    "Psi4",
    "PyPsi4",
    "PyXTB",
    "Turbomole",
    "XTB",
]


from pysisyphus.calculators.AFIR import AFIR
from pysisyphus.calculators.Dimer import Dimer
from pysisyphus.calculators.ExternalPotential import ExternalPotential
from pysisyphus.calculators.FakeASE import FakeASE
from pysisyphus.calculators.Gaussian16 import Gaussian16
from pysisyphus.calculators.LennardJones import LennardJones
from pysisyphus.calculators.MOPAC import MOPAC
from pysisyphus.calculators.Psi4 import Psi4
from pysisyphus.calculators.ORCA import ORCA
from pysisyphus.calculators.PyPsi4 import PyPsi4
from pysisyphus.calculators.PyXTB import PyXTB
from pysisyphus.calculators.Turbomole import Turbomole
from pysisyphus.calculators.XTB import XTB


logger = logging.getLogger("dimer")
logger.setLevel(logging.DEBUG)
handler = logging.FileHandler("dimer.log", mode="w", delay=True)
fmt_str = "%(message)s"
formatter = logging.Formatter(fmt_str)
handler.setFormatter(formatter)
logger.addHandler(handler)
