#!/usr/bin/env python3

# See [1] 10.1002/jcc.21026

import bisect
from math import cos, sin

import numpy as np
import rmsd

from pysisyphus.calculators.XTB import XTB
from pysisyphus.Geometry import Geometry
from pysisyphus.stocastic.Kick import Kick
from pysisyphus.xyzloader import make_trj_str_from_geoms


np.set_printoptions(suppress=True, precision=2)


class FragmentKick(Kick):

    def __init__(self, geom, fragments, fix_fragments=list(), **kwargs):
        super().__init__(geom, **kwargs)

        self.fragments = [np.array(frag) for frag in fragments]
        self.fix_fragments = fix_fragments

        # Shift fragment coordinates into their centroid
        frag_coords = self.get_frag_coords(self.initial_geom)
        self.frag_coords = [fc - fc.mean(axis=0) for fc in frag_coords]

        new_coords3d = np.concatenate(self.frag_coords)
        geom = Geometry(self.atoms, new_coords3d.flatten())
        with open("fragments_in_origin.xyz", "w") as handle:
            handle.write(geom.as_xyz())

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
        kick = self.get_kick()[:3]
        # Fragment rotation
        rot_coords = R.dot(frag_coords.T).T
        # Fragment translation
        rot_kicked_coords = rot_coords + kick
        return rot_kicked_coords

    def get_frag_coords(self, geom):
        return [geom.coords3d[frag] for frag in self.fragments]

    # def get_kick(self):
        # # Interval [0, 1)
        # kick = np.random.random(self.coords_size)
        # # Stretch [0, 1) to  [-r ... r)
        # kick = self.radius * (2*kick - 1)
        # kick = kick.reshape(-1, 3)
        # return kick.flatten()

    def get_input_geom(self, geom):
        # frag_coords = self.get_frag_coords()
        frag_coords = self.frag_coords
        # kicked_frags = [self.kick_fragment(fc) for fc in frag_coords]
        kicked_frags = list()
        for i, fc in enumerate(frag_coords):
            if i in self.fix_fragments:
                kicked_frag = fc
            else:
                kicked_frag = self.kick_fragment(fc)
            kicked_frags.append(kicked_frag)
        new_coords3d = np.concatenate(kicked_frags)
        new_coords = rmsd.kabsch_rotate(new_coords3d,
                                        self.initial_coords3d
        ).flatten()
        new_geom = Geometry(self.atoms, new_coords)
        return new_geom
