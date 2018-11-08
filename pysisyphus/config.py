import configparser
import pathlib
import sys

home = pathlib.Path.home()
Config = configparser.ConfigParser()
config_fn = home / ".pysisyphusrc"
if not config_fn.is_file():
    print(f"Couldn't find configuration file. Expected it at '{config_fn}'.")
    sys.exit()
read_fns = Config.read(config_fn)
