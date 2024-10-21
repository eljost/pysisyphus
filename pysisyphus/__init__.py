import logging
import sys

from .version import version as __version__


logger = logging.getLogger("pysis")
logger.setLevel(logging.INFO)

file_handler = logging.FileHandler("pysisyphus.log", mode="w", delay=True)
file_handler.setLevel(logging.INFO)
logger.addHandler(file_handler)

stdout_handler = logging.StreamHandler(sys.stdout)
stdout_handler.setLevel(logging.INFO)
logger.addHandler(stdout_handler)
