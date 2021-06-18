import logging

__all__ = [
    "AFIR",
    "AtomAtomTransTorque",
    "Composite",
    "DFTBp",
    "Dimer",
    "ExternalPotential",
    "FakeASE",
    "Gaussian09",
    "Gaussian16",
    "HardSphere",
    "IPIServer",
    "LennardJones",
    "MOPAC",
    "MultiCalc",
    "ONIOM",
    "OpenMolcas",
    "ORCA",
    "Psi4",
    "PyPsi4",
    "PyXTB",
    "TIP3P",
    "Turbomole",
    "TransTorque",
    "XTB",
]


from pysisyphus.calculators.AFIR import AFIR
from pysisyphus.calculators.AtomAtomTransTorque import AtomAtomTransTorque
from pysisyphus.calculators.Composite import Composite
from pysisyphus.calculators.DFTBp import DFTBp
from pysisyphus.calculators.Dimer import Dimer
from pysisyphus.calculators.ExternalPotential import ExternalPotential
from pysisyphus.calculators.FakeASE import FakeASE
from pysisyphus.calculators.Gaussian16 import Gaussian16
from pysisyphus.calculators.IPIServer import IPIServer
from pysisyphus.calculators.HardSphere import HardSphere
from pysisyphus.calculators.LennardJones import LennardJones
from pysisyphus.calculators.MultiCalc import MultiCalc
from pysisyphus.calculators.MOPAC import MOPAC
from pysisyphus.calculators.Psi4 import Psi4
from pysisyphus.calculators.ONIOMv2 import ONIOM
from pysisyphus.calculators.OpenMolcas import OpenMolcas
from pysisyphus.calculators.ORCA import ORCA
from pysisyphus.calculators.PyPsi4 import PyPsi4
from pysisyphus.calculators.PyXTB import PyXTB
from pysisyphus.calculators.TIP3P import TIP3P
from pysisyphus.calculators.TransTorque import TransTorque
from pysisyphus.calculators.Turbomole import Turbomole
from pysisyphus.calculators.XTB import XTB


logger = logging.getLogger("dimer")
logger.setLevel(logging.DEBUG)
handler = logging.FileHandler("dimer.log", mode="w", delay=True)
fmt_str = "%(message)s"
formatter = logging.Formatter(fmt_str)
handler.setFormatter(formatter)
logger.addHandler(handler)
