import configparser
import os
from pathlib import Path
import sys



CONFIG_DIR = Path(os.path.abspath(os.path.dirname(__file__)))
LIB_DIR = CONFIG_DIR / "geom_library"
T_DEFAULT = 298.15  # Kelvin
p_DEFAULT = 101325  # Pascal
OUT_DIR_DEFAULT = "qm_calcs"

DEFAULTS = {
    "mwfn": "Multiwfn",
    "jmol": "jmol",
    "packmol": "packmol",
    "formchk": "formchk",
    "unfchk": "unfchk",
    "rwfdump": "rwfdump",
}


# First try to read pysisyphusrc path an environment variable
try:
    pysisrc_env = os.getenv("PYSISRC", default=None)
    config_fn = Path(pysisrc_env).resolve()
    print(f"Read pysisyphus configuration from '{config_fn}'")
# Fallback to $HOME/.pysisyphusrc
except TypeError:
    home = Path.home()
    config_fn = Path.home() / ".pysisyphusrc"

if not config_fn.is_file():
    print(f"Couldn't find configuration file. Expected it at '{config_fn}'.")

Config = configparser.ConfigParser()
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
