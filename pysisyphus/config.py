import argparse
import configparser
from io import StringIO
import os
from pathlib import Path
import shutil
import sys

CONFIG_DIR = Path(os.path.abspath(os.path.dirname(__file__)))
LIB_DIR = CONFIG_DIR / "geom_library"
T_DEFAULT = 298.15  # Kelvin
p_DEFAULT = 101325  # Pascal
OUT_DIR_DEFAULT = "qm_calcs"

DEFAULTS = {
    # .pysisyphusrc key: command
    "mwfn": "Multiwfn",
    "jmol": "jmol",
    "packmol": "packmol",
    "formchk": "formchk",
    "unfchk": "unfchk",
    "rwfdump": "rwfdump",
    # QC codes
    "orca": "orca",
    "orca5": "orca",
    "gaussian16": "g16",
    "wfoverlap": "wfoverlap.x",
    "openmolcas": "pymolcas",
    "gamess": "rungms",
    "xtb": "xtb",
    "mopac": "mopac",
}


# First try to read path to .pysisyphusrc from environment variable
try:
    pysisrc_env = os.getenv("PYSISRC", default=None)
    config_fn = Path(pysisrc_env).resolve()
    print(f"Read pysisyphus configuration from '{config_fn}'")
# Fallback to $HOME/.pysisyphusrc
except TypeError:
    config_fn = Path.home() / ".pysisyphusrc"

if not config_fn.is_file():
    print(f"Couldn't find configuration file. Expected it at '{config_fn}'.")

Config = configparser.ConfigParser()
read_fns = Config.read(config_fn)


def get_cmd(key, use_defaults=True):
    try:
        cmd = Config[key]["cmd"]
    except KeyError:
        if use_defaults:
            try:
                cmd = shutil.which(DEFAULTS[key])
            except KeyError:
                cmd = None
    if not cmd:
        print(
            f"Failed to load key 'cmd' from section '{key}' "
            "in ~/.pysisyphusrc and no default was specified."
        )
        sys.exit()
    return cmd


def detect_paths():
    config = configparser.ConfigParser()

    for k, v in DEFAULTS.items():
        print(f"{k: >16}: ... ", end="")
        if not (path := shutil.which(v)):
            print("not found.")
        else:
            path = Path(path).resolve()
            config[k] = {"cmd": str(path)}
            print(path)

    fp = StringIO()
    config.write(fp, space_around_delimiters=False)
    fp.seek(0)
    config_text = fp.read()
    return config_text


def parse_args(args):
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--out", "-o", help="Filename to which the detected cmds are dumped."
    )

    return parser.parse_args(args)


def run_detect_paths():
    args = parse_args(sys.argv[1:])
    out = args.out

    config_text = detect_paths()
    if out:
        print("Dumped detected cmds to '{out}'")
        with open(out, "w") as handle:
            handle.write(config_text)
    else:
        print(config_text)
