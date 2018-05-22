#!/usr/bin/env python3

import logging
import os
from pathlib import Path
import pickle
import platform
import socketserver
import struct
import threading

from distributed import Client

from pysisyphus.helpers import slugify_worker


def get_fh_logger(name, log_fn):
    logger = logging.getLogger(name)
    logger.setLevel(logging.DEBUG)
    if len(logger.handlers) == 0:
        fh = logging.FileHandler(log_fn, mode="w", delay=True)
        fh.setLevel(logging.DEBUG)
        fmt_str = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
        formatter = logging.Formatter(fmt_str)
        fh.setFormatter(formatter)
        logger.addHandler(fh)
        logger.debug(f"Initialized logging on {platform.node()}")


def init_logging_base(dask_worker, log_path):
    slug = slugify_worker(dask_worker.worker_address)
    log_fn = log_path / f"{slug}_calculator.log"
    get_fh_logger("calculator", log_fn)
    log_fn = log_path / f"{slug}_wfoverlap.log"
    get_fh_logger("wfoverlap", log_fn)


def init_logging(log_dir="./", scheduler=None):
    log_path = Path(log_dir)
    if scheduler:
        client = Client(scheduler)
        client.restart()
        client.run(init_logging_base, log_path=log_path)
    else:
        log_fn = log_path / "calculator.log"
        get_fh_logger("calculator", log_fn)
        log_fn = log_path / "wfoverlap.log"
        get_fh_logger("wfoverlap", log_fn)
