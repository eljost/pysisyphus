import logging
import sys


logger = logging.getLogger("overlapper")
logger.setLevel(logging.DEBUG)
# delay = True prevents creation of empty logfiles
handler = logging.FileHandler("overlapper.log", mode="w", delay=True)
fmt_str = "%(message)s"
formatter = logging.Formatter(fmt_str)
handler.setFormatter(formatter)
logger.addHandler(handler)
# Prints to stdout
logger.addHandler(logging.StreamHandler(sys.stdout))
