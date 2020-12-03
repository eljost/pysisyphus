import os
from pathlib import Path

from pysisyphus.db.level import LEVELS
from pysisyphus.db.molecules import MOLECULES


THIS_DIR = Path(os.path.dirname(os.path.realpath(__file__)))
GUESS_DIR = THIS_DIR / "guess"
LEVEL_DIR = THIS_DIR / "levels"
OPT_FLAG = "##OPT##"
