import configparser
import pathlib

home = pathlib.Path.home()
Config = configparser.ConfigParser()
Config.read(home / ".pysisyphusrc")
