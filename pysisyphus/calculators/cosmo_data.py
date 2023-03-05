from dataclasses import dataclass
from enum import IntEnum


_JSON_DATA = {
    # Generated from 'radii.cosmo'
    #
    # Version 1.0 080606
    # entries per line:
    # 1) symbol                  (two lower case characters) 
    # 2) radius in Angstroem     (double)  
    # 3) quality of rad. (0|1|2) (integer)  0=some guess; 1=reasonable guess; 2=optimized
    # 4) isodens in e/bohr**3    (double)
    # For entry 2 and 4, -1.0 denotes an undefined value
    "h": [1.3, 2, 0.0018],
    "he": [1.638, 1, -1.0],
    "li": [1.404, 1, -1.0],
    "be": [1.053, 1, -1.0],
    "b": [2.0475, 1, 0.001],
    "c": [2.0, 2, 0.0016],
    "n": [1.83, 2, 0.0017],
    "o": [1.72, 2, 0.0013],
    "f": [1.72, 2, 0.0011],
    "ne": [1.8018, 1, -1.0],
    "na": [1.755, 0, -1.0],
    "mg": [1.638, 0, -1.0],
    "al": [2.153, 1, 0.0035],
    "si": [2.2, 1, 0.0013],
    "p": [2.106, 1, 0.0016],
    "s": [2.16, 2, 0.0009],
    "cl": [2.05, 2, 0.0014],
    "ar": [2.223, 0, -1.0],
    "k": [2.223, 0, -1.0],
    "ca": [2.223, 0, -1.0],
    "sc": [2.223, 0, -1.0],
    "ti": [2.223, 0, -1.0],
    "v": [2.223, 0, -1.0],
    "cr": [2.223, 0, -1.0],
    "mn": [2.223, 0, -1.0],
    "fe": [2.223, 0, -1.0],
    "co": [2.223, 0, -1.0],
    "ni": [2.223, 0, -1.0],
    "cu": [2.223, 0, -1.0],
    "zn": [1.626, 1, -1.0],
    "ga": [2.223, 0, -1.0],
    "ge": [2.7, 1, 0.001],
    "as": [2.35, 1, 0.0011],
    "se": [2.2, 1, 0.0013],
    "br": [2.16, 2, 0.0012],
    "kr": [2.223, 0, -1.0],
    "rb": [2.223, 0, -1.0],
    "sr": [2.223, 0, -1.0],
    "y": [2.223, 0, -1.0],
    "zr": [2.223, 0, -1.0],
    "nb": [2.223, 0, -1.0],
    "mo": [2.223, 0, -1.0],
    "tc": [2.223, 0, -1.0],
    "ru": [2.223, 0, -1.0],
    "rh": [2.223, 0, -1.0],
    "pd": [2.223, 0, -1.0],
    "ag": [2.223, 0, -1.0],
    "cd": [2.223, 0, -1.0],
    "in": [2.223, 0, -1.0],
    "sn": [2.55, 0, -1.0],
    "sb": [2.223, 0, -1.0],
    "te": [2.223, 0, -1.0],
    "i": [2.32, 2, 0.0016],
    "xe": [2.223, 0, -1.0],
    "cs": [2.223, 0, -1.0],
    "ba": [2.223, 0, -1.0],
    "la": [2.223, 0, -1.0],
    "ce": [2.223, 0, -1.0],
    "pr": [2.223, 0, -1.0],
    "nd": [2.223, 0, -1.0],
    "pm": [2.223, 0, -1.0],
    "sm": [2.223, 0, -1.0],
    "eu": [2.223, 0, -1.0],
    "gd": [2.223, 0, -1.0],
    "tb": [2.223, 0, -1.0],
    "dy": [2.223, 0, -1.0],
    "ho": [2.223, 0, -1.0],
    "er": [2.223, 0, -1.0],
    "tm": [2.223, 0, -1.0],
    "yb": [2.223, 0, -1.0],
    "lu": [2.223, 0, -1.0],
    "hf": [2.223, 0, -1.0],
    "ta": [2.223, 0, -1.0],
    "w": [2.223, 0, -1.0],
    "re": [2.223, 0, -1.0],
    "os": [2.223, 0, -1.0],
    "ir": [2.223, 0, -1.0],
    "pt": [2.223, 0, -1.0],
    "au": [2.223, 0, -1.0],
    "hg": [2.223, 0, -1.0],
    "tl": [2.223, 0, -1.0],
    "pb": [2.36, 1, -1.0],
    "bi": [2.223, 0, -1.0],
    "po": [2.223, 0, -1.0],
    "at": [2.223, 0, -1.0],
    "rn": [2.223, 0, -1.0],
    "fr": [2.223, 0, -1.0],
    "ra": [2.223, 0, -1.0],
    "ac": [2.223, 0, -1.0],
    "th": [2.223, 0, -1.0],
    "pa": [2.223, 0, -1.0],
    "u": [2.223, 0, -1.0],
    "np": [2.223, 0, -1.0],
    "pu": [2.223, 0, -1.0],
    "q": [0.0, 0, -1.0],
}


class COSMORadiusQuality(IntEnum):
    GUESS = 0
    REASONABLE_GUESS = 1
    OPTIMIZED = 2


cosmo_qualities ={
    0: COSMORadiusQuality["GUESS"],
    1: COSMORadiusQuality["REASONABLE_GUESS"],
    2: COSMORadiusQuality["OPTIMIZED"],
}


@dataclass
class COSMORadius:
    atom: str
    radius: float  # in Angstrom
    quality: COSMORadiusQuality
    isodens: float  # e/bohr³


COSMO_RADII = {
    atom: COSMORadius(
        atom=atom, radius=radius, quality=cosmo_qualities[quality], isodens=isodens
    )
    for atom, (radius, quality, isodens) in _JSON_DATA.items()
}
