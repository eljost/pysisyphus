#!/usr/bin/env python3

# See [1] 10.1002/jcc.21026

from math import cos, sin

import numpy as np
import rmsd

from pysisyphus.Geometry import Geometry
from pysisyphus.stocastic.Kick import Kick


np.set_printoptions(suppress=True, precision=2)


class FragmentKick(Kick):

    def __init__(self, geom, fragments, fix_fragments=list(),
                 displace_from=list(), random_displacement=False,
                 **kwargs):
        super().__init__(geom, **kwargs)

        self.fragments = self.get_fragments(fragments)
        self.fix_fragments = fix_fragments
        self.displace_from = displace_from
        # True when displace_from is given
        self.random_displacement = self.displace_from or random_displacement

        # Shift fragment coordinates into their centroid
        frag_coords = [self.initial_geom.coords3d[frag] for frag in self.fragments]
        self.frag_coords = [fc - fc.mean(axis=0) for fc in frag_coords]

        new_coords3d = np.concatenate(self.frag_coords)
        geom = Geometry(self.atoms, new_coords3d.flatten())
        with open("fragments_in_origin.xyz", "w") as handle:
            handle.write(geom.as_xyz())

    def get_fragments(self, fragments):
        if len(fragments) == 1:
            fragment = fragments[0]
            self.log("Fragments were not fully specified. "
                     "Determining remaining fragment.")
            all_indices = set(range(len(self.initial_geom.atoms)))
            second_frag = tuple(all_indices - set(fragment))
            fragments = (fragment, second_frag)
        return [np.array(frag) for frag in fragments]

    def get_origin(self):
        if not self.random_displacement:
            return (0, 0, 0)

        # Coordinates of fixed fragment
        frag_coords = self.frag_coords[self.fix_fragments[0]]
        """Here we have two options:
        1.) displace_from is not given, so we use all possible
            atom indices in the fragment.
        2.) displace_from is given, so we restrict the valid
            atom indices. This way some degree of symmetry can
            be considered and the generated structures should
            be more similar. This should reduce the effort in
            matching the geometries later on with the hungarian
            method.

        Consider a stocastic search for minima between benzene and
        a chlorine molecule. As all carbons/hydrogens in the benzene
        are equivalent one could greatly reduce the effort of matching
        atoms later on, when displacing only from one carbon and one
        hydrogen."""
        if not self.displace_from:
            frag_indices = np.arange(frag_coords.shape[0])
        else:
            frag_indices = self.displace_from

        # Select random atom of fixed fragment
        random_ind = np.random.choice(frag_indices, size=1)[0]
        return frag_coords[random_ind]

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
        rot_kicked_coords = rot_coords + self.get_origin() + kick
        return rot_kicked_coords

    def get_input_geom(self, geom):
        kicked_frags = list()
        for i, fc in enumerate(self.frag_coords):
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
