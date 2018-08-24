import logging
import sys

logger = logging.getLogger("stocastic")
logger.setLevel(logging.DEBUG)
fh = logging.FileHandler("stocastic.log", mode="w", delay=True)
fh.setLevel(logging.DEBUG)
logger.addHandler(logging.StreamHandler(sys.stdout))
# fmt_str = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
# formatter = logging.Formatter(fmt_str)
# fh.setFormatter(formatter)
logger.addHandler(fh)

__all__ = [
    "FragmentKick",
    "Kick",
]
