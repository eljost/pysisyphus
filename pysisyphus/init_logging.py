#!/usr/bin/env python3

import logging
import os
from pathlib import Path
import platform

from distributed import Client


def init_logging(log_dir="./", scheduler=None):
    log_dir = Path(log_dir)
    init_logging_base("calculator", log_dir / "calculator.log",
                      scheduler=scheduler)
    init_logging_base("wfoverlap", log_dir / "wfoverlap.log",
                      scheduler=scheduler)


def init_logging_base(name, log_path, mode="w", scheduler=None):
    logger = logging.getLogger(name)
    logger.setLevel(logging.DEBUG)
    """If we reuse an existing scheduler we have to reset the handlers
    everytime. This shouldn't be an issue in production calculations,
    as we always start a new scheduler."""
    logger.handlers = []
    if not len(logger.handlers):
        fh = logging.FileHandler(log_path, mode=mode)
        fh.setLevel(logging.DEBUG)
        fmt_str = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
        formatter = logging.Formatter(fmt_str)
        fh.setFormatter(formatter)
        logger.addHandler(fh)
        logger.debug(f"Initialized logging on {platform.node()}")
    # Also initalize logging on all connected workers, but now we append
    # the the previously created logfile.
    if scheduler:
        client = Client(scheduler)
        client.run(init_logging_base, name, log_path, "a")
