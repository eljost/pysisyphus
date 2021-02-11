# See [1] 10.1002/jcc.21026

import itertools as it

import numpy as np
import rmsd

from pysisyphus.Geometry import Geometry
from pysisyphus.linalg import get_rot_mat
from pysisyphus.stocastic.Kick import Kick


class FragmentKick(Kick):
    def __init__(
        self,
        geom,
        fragments,
        fix_fragments=None,
        displace_from=None,
        random_displacement=False,
        **kwargs,
    ):
        super().__init__(geom, **kwargs)

        self.fragments = self.get_fragments(fragments)
        if fix_fragments is None:
            fix_fragments = (0,)
        self.fix_fragments = fix_fragments

        if displace_from is None:
            displace_from = list()
        self.displace_from = displace_from
        # True when displace_from is given
        self.random_displacement = self.displace_from or random_displacement

        frag_coords = [self.initial_geom.coords3d[frag] for frag in self.fragments]
        # Shift centroids of the fragment into the cartesian origin
        self.frag_coords = [fc - fc.mean(axis=0) for fc in frag_coords]
        atoms_arr = np.array(self.initial_geom.atoms)
        self.frag_atoms = [atoms_arr[frag] for frag in self.fragments]

        new_coords3d = np.concatenate(self.frag_coords)
        geom = Geometry(list(it.chain(*self.frag_atoms)), new_coords3d.flatten())
        with open("fragments_in_origin.xyz", "w") as handle:
            handle.write(geom.as_xyz())

        self.frag_inds_flat = list(it.chain(*self.fragments))
        # Given an array with fragment coordinates these indices will sort the
        # array into the original atom order.
        self.sort_inds = np.argsort(self.frag_inds_flat)
        self.print_fragments()

    def print_fragments(self):
        lines = self.initial_geom.as_xyz().split("\n")[2:]
        for i, fragment in enumerate(self.fragments):
            self.log(f"Fragment {i}:")
            frag_lines = [lines[j] for j in fragment]
            self.log("\n".join(frag_lines) + "\n")

    def get_fragments(self, fragments):
        fragments = list(fragments)
        # Compare number of atoms defined in fragments with the total
        # number of atoms in the molecule.
        fragment_atoms = list(it.chain(*fragments))
        if len(fragment_atoms) != len(self.atoms):
            self.log(
                "Fragments were not fully specified. " "Determining remaining fragment."
            )
            all_indices = set(range(len(self.atoms)))
            new_fragment = list(all_indices - set(fragment_atoms))
            fragments.append(new_fragment)

        fragments = [np.array(frag) for frag in fragments]
        return fragments

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

    def kick_fragment(self, frag_coords):
        R = get_rot_mat()
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
        # As the fragments are not necessarily defined in order, especially
        # when we generate the second fragment automatically, we have to
        # sort the coordinates, so they are in the original order. E.g. the
        # fragments ((4, 3), (2, 1, 0)) would result 'new_coords3d' with the
        # coordinates of atom 4 at index 0, coordinates of atom 3 at index 1, etc.
        new_coords3d = new_coords3d[self.sort_inds]
        new_coords = rmsd.kabsch_rotate(new_coords3d, self.initial_coords3d).flatten()
        new_geom = Geometry(self.atoms, new_coords)
        return new_geom
