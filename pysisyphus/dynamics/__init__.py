import logging

from pysisyphus.dynamics.colvars import get_colvar as get_colvar
from pysisyphus.dynamics.Gaussian import Gaussian as Gaussian
from pysisyphus.dynamics.helpers import get_mb_velocities_for_geom as get_mb_velocities_for_geom
from pysisyphus.dynamics.mdp import mdp as mdp
from pysisyphus.dynamics.rattle import rattle_closure as rattle_closure
from pysisyphus.dynamics.driver import md as md
from pysisyphus.dynamics.wigner import get_wigner_sampler as get_wigner_sampler


logger = logging.getLogger("dynamics")
logger.setLevel(logging.DEBUG)
# delay = True prevents creation of empty logfiles
handler = logging.FileHandler("dynamics.log", mode="w", delay=True)
# fmt_str = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
fmt_str = "%(asctime)s - %(message)s"
formatter = logging.Formatter(fmt_str)
handler.setFormatter(formatter)
logger.addHandler(handler)
