import logging

logger = logging.getLogger("pipeline")
logger.setLevel(logging.DEBUG)
fh = logging.FileHandler("pipeline.log", mode="w", delay=True)
fh.setLevel(logging.DEBUG)
# fmt_str = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
# formatter = logging.Formatter(fmt_str)
# fh.setFormatter(formatter)
logger.addHandler(fh)
