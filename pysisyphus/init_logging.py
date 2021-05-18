import logging
from pathlib import Path

from distributed import Client

from pysisyphus.helpers import slugify_worker

LOGGERS = {
    "calculator": "calculator.log",
    "wfoverlap": "wfoverlap.log",
}


def get_fh_logger(name, log_fn):
    """Initialize a logger with 'name', level DEBUG and a FileHandler."""
    logger = logging.getLogger(name)
    logger.setLevel(logging.DEBUG)
    if len(logger.handlers) == 0:
        fh = logging.FileHandler(log_fn, mode="w", delay=True)
        fh.setLevel(logging.DEBUG)
        # fmt_str = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
        # fmt_str = "%(asctime)s - %(message)s"
        fmt_str = "%(asctime)s - %(message)s"
        datefmt = "%y-%m-%d %H:%M:%S"
        formatter = logging.Formatter(fmt_str, datefmt=datefmt)
        fh.setFormatter(formatter)
        logger.addHandler(fh)
        # Uncommented this for now as the host is already in the filename
        # logger.debug(f"Initialized logging on {platform.node()}")


def init_logging_base(dask_worker, log_path):
    """Prepare individual loggers for one dask_worker."""
    slug = slugify_worker(dask_worker.worker_address)
    for name, log_fn_base in LOGGERS.items():
        log_fn = log_path / f"{slug}_{log_fn_base}"
        get_fh_logger(name, log_fn)


def init_logging(log_dir="./", scheduler=None):
    """Prepare the logger in log_path. When called with scheduler
    loggers for every worker are prepared."""
    log_path = Path(log_dir)
    if scheduler:
        client = Client(scheduler)
        client.run(init_logging_base, log_path=log_path)
    else:
        for name, log_fn_base in LOGGERS.items():
            log_fn = log_path / log_fn_base
            get_fh_logger(name, log_fn)
