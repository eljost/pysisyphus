from dataclasses import dataclass


@dataclass
class IRCDummy:
    all_coords: list
    atoms: tuple
    forward: bool = True
    backward: bool = True
    downhill: bool = False
