#!/usr/bin/env python3

import configparser
import pathlib
import sys

DEFAULTS = {
    "mwfn": "Multiwfn",
    "jmol": "jmol",
}

home = pathlib.Path.home()
Config = configparser.ConfigParser()
config_fn = home / ".pysisyphusrc"
if not config_fn.is_file():
    print(f"Couldn't find configuration file. Expected it at '{config_fn}'.")
read_fns = Config.read(config_fn)


def get_cmd(key, use_defaults=True):
    try:
        return Config[key]["cmd"]
    except KeyError:
        if (not use_defaults) or (key not in DEFAULTS):
            print(f"Failed to load key 'cmd' from section '{key}' "
                   "in ~/.pysisyphusrc and no default was specified. Exiting!"
            )
            sys.exit()
        return DEFAULTS[key]
