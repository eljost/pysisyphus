import argparse
import configparser
from io import StringIO
import os
from pathlib import Path
import shutil
import sys


from pysisyphus import logger


CONFIG_DIR = Path(os.path.abspath(os.path.dirname(__file__)))
LIB_DIR = CONFIG_DIR / "geom_library"
T_DEFAULT = 298.15  # Kelvin
p_DEFAULT = 101325  # Pascal
OUT_DIR_DEFAULT = "qm_calcs"
BEND_MIN_DEG = 15
LB_MIN_DEG = 175
DIHED_MAX_DEG = 175
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


def get_cmd(section, key="cmd", use_defaults=True):
    cmd = None
    msg = f" and no default was specified for '{section}'."

    # First we try to load the command from .pysisyphusrc
    try:
        cmd = Config[section][key]
        msg = "."
    # When the command is not available on .pysisyphusrc we check the defaults
    except KeyError:
        if use_defaults:
            default_key = section if key == "cmd" else key
            try:
                default_cmd = DEFAULTS[default_key]
                cmd = shutil.which(default_cmd)
                # 'msg' will only be printed when 'cmd' is None, so the msg is
                # always negative.
                msg = f" and default cmd='{default_cmd}' was not found in $PATH."
            except KeyError:
                pass

    if cmd is None:
        logger.warning(
            f"Failed to load '{key}' from [{section}] "
            f"in ~/.pysisyphusrc{msg}"
        )
        cmd = None
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
        print("\nExample .pysisyphusrc:\n")
        print(config_text)
