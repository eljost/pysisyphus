import logging


logger = logging.getLogger("pysis.internal_coords")
logger.setLevel(logging.DEBUG)
handler = logging.FileHandler("internal_coords.log", mode="w", delay=True)
handler.setLevel(logging.DEBUG)
logger.addHandler(handler)
