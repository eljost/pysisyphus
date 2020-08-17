import logging


logger = logging.getLogger("mwfn")
logger.setLevel(logging.DEBUG)
handler = logging.FileHandler("mwfn.log", mode="w", delay=True)
# fmt_str = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
fmt_str = "%(asctime)s %(message)s"
formatter = logging.Formatter(fmt_str, datefmt='%Y-%m-%d %H:%M:%S')
handler.setFormatter(formatter)
logger.addHandler(handler)
