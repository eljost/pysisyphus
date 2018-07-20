#!/usr/bin/env python3

from math import sin, cos

import numpy as np
import rmsd

from pysisyphus.Geometry import Geometry
from pysisyphus.stocastic.Kick import Kick



np.set_printoptions(suppress=True, precision=2)


class FragmentKick(Kick):

    def __init__(self, geom, fragments, rmsd_thresh=0.25, **kwargs):
        self.fragments = [np.array(frag) for frag in fragments]
        # atoms_arr = np.array(geom.atoms)
        # self.fragment_atoms = [atoms_arr[frag].tolist()
                               # for frag in self.fragments]
        super().__init__(geom, rmsd_thresh=rmsd_thresh, **kwargs)

    def get_rot_mat(self):
        # Euler angles
        a, b, c = np.random.rand(3)*np.pi*2
        R = np.array((
            (cos(a)*cos(b)*cos(c)-sin(a)*sin(c),
            -cos(a)*cos(b)*sin(c)-sin(a)*cos(c),
             cos(a)*sin(b)),
            (sin(a)*cos(b)*cos(c)+cos(a)*sin(c),
            -sin(a)*cos(b)*sin(c)+cos(a)*cos(c),
             sin(a)*sin(b)),
            (-sin(b)*cos(c),
              sin(b)*sin(c),
              cos(b)))
        )
        return R

    def kick_fragment(self, frag_coords):
        R = self.get_rot_mat()
        # Fragment rotation
        rot_coords = R.dot(frag_coords.T).T
        kick = self.get_kick()[:3]
        # Fragment translation
        rot_kicked_coords = rot_coords + kick
        return rot_kicked_coords

    def get_frag_coords(self, geom):
        return [geom.coords3d[frag] for frag in self.fragments]

    def get_kicked_geom(self, geom):
        frag_coords = self.get_frag_coords(geom)
        kicked_frags = [self.kick_fragment(fc) for fc in frag_coords]
        new_coords3d = np.concatenate(kicked_frags)
        new_coords = rmsd.kabsch_rotate(new_coords3d,
                                        self.initial_coords3d
        ).flatten()
        new_geom = Geometry(geom.atoms, new_coords)
        return new_geom
