import logging

from pysisyphus.dynamics.mdp import mdp
from pysisyphus.dynamics.rattle import rattle_closure


logger = logging.getLogger("dynamics")
logger.setLevel(logging.DEBUG)
# delay = True prevents creation of empty logfiles
handler = logging.FileHandler("dynamics.log", mode="w", delay=True)
# fmt_str = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
fmt_str = "%(asctime)s - %(message)s"
formatter = logging.Formatter(fmt_str)
handler.setFormatter(formatter)
logger.addHandler(handler)
