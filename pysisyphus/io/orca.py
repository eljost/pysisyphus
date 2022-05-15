import json

import numpy as np

from pysisyphus.integrals import get_l, Shell, Shells


def parse_json_shells(text):
    data = json.loads(text)
    atoms = data["Molecule"]["Atoms"]

    _shells = list()
    for atom in atoms:
        center = np.array(atom["Coords"])
        center_ind = atom["Idx"]
        bfs = atom["BasisFunctions"]
        for bf in bfs:
            exps = np.array(bf["Exponents"])
            coeffs = np.array(bf["Coefficients"])
            L = get_l(bf["Shell"])
            shell = Shell(
                L=L,
                center=center,
                coeffs=coeffs,
                exps=exps,
                center_ind=center_ind,
            )
            _shells.append(shell)
    shells = Shells(_shells)
    return shells
