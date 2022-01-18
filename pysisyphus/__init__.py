from importlib.metadata import version, PackageNotFoundError
import logging
import sys

try:
    __version__ = version("pysisyphus")
except PackageNotFoundError:
    # package is not installed
    pass

logger = logging.getLogger("pysisyphus")
logger.setLevel(logging.DEBUG)

file_handler = logging.FileHandler("pysisyphus.log", mode="w", delay=True)
logger.addHandler(file_handler)

stdout_handler = logging.StreamHandler(sys.stdout)
stdout_handler.setLevel(logging.INFO)
logger.addHandler(stdout_handler)
