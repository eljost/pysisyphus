import configparser
import pathlib

home = pathlib.Path.home()
Config = configparser.ConfigParser()
config_fn = home / ".pysisyphusrc"
read_fns = Config.read(config_fn)
if read_fns == []:
    print(f"Couldn't read configuration file. Expected it at '{config_fn}'.")
