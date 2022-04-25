import logging

from pysisyphus import logger as pysis_logger


logger = pysis_logger.getChild("internal_coords")
logger.setLevel(logging.DEBUG)
handler = logging.FileHandler("internal_coords.log", mode="w", delay=True)
handler.setLevel(logging.DEBUG)
logger.addHandler(handler)
