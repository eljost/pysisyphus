import logging
import sys

from .version import version as __version__

logger = logging.getLogger("pysisyphus")
logger.setLevel(logging.DEBUG)

file_handler = logging.FileHandler("pysisyphus.log", mode="w", delay=True)
logger.addHandler(file_handler)

stdout_handler = logging.StreamHandler(sys.stdout)
stdout_handler.setLevel(logging.INFO)
logger.addHandler(stdout_handler)
