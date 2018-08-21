#!/usr/bin/env python3

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

    def __init__(self, geom, fragments, rmsd_thresh=0.25,
                 fix_fragments=list(), **kwargs):
        self.fragments = [np.array(frag) for frag in fragments]
        self.fix_fragments = fix_fragments
        # atoms_arr = np.array(geom.atoms)
        # self.fragment_atoms = [atoms_arr[frag].tolist()
                               # for frag in self.fragments]
        super().__init__(geom, rmsd_thresh=rmsd_thresh, **kwargs)

        # Shift fragment coordinates into their centroid
        frag_coords = self.get_frag_coords(self.initial_geom)
        self.frag_coords = [fc - fc.mean(axis=0) for fc in frag_coords]

        new_coords3d = np.concatenate(self.frag_coords)
        geom = Geometry(self.atoms, new_coords3d.flatten())
        with open("origin.xyz", "w") as handle:
            handle.write(geom.as_xyz())

        self.energies = list()
        self.opt_geoms = list()


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

    def get_kick(self):
        # Interval [0, 1)
        kick = np.random.random(self.coords_size)
        # Stretch [0, 1) to  [-r ... r)
        kick = self.radius * (2*kick - 1)
        kick = kick.reshape(-1, 3)
        return kick.flatten()

    def get_kicked_geom(self):
        # frag_coords = self.get_frag_coords()
        frag_coords = self.frag_coords
        # kicked_frags = [self.kick_fragment(fc) for fc in frag_coords]
        kicked_frags = list()
        for i, fc in enumerate(frag_coords):
            if i in self.fix_fragments:
                kicked_frag = fc
                print(f"keep fragment {i} fixed.")
            else:
                kicked_frag = self.kick_fragment(fc)
            kicked_frags.append(kicked_frag)
        new_coords3d = np.concatenate(kicked_frags)
        print(new_coords3d)
        new_coords = rmsd.kabsch_rotate(new_coords3d,
                                        self.initial_coords3d
        ).flatten()
        new_geom = Geometry(self.atoms, new_coords)
        return new_geom

    def run_kicked_geom(self, geom):
        calc = XTB(calc_number=self.calc_counter)
        self.calc_counter += 1
        opt_geom = calc.run_opt(geom.atoms, geom.coords, keep=False)
        return opt_geom

    def geoms_to_trj(self, geoms, fn):
        with open(fn, "w") as handle:
            handle.write(make_trj_str_from_geoms(geoms))

    def run(self):
        def are_close(en1, en2, thresh=1e-4):
            return abs(en1 - en2) < thresh
        while self.cur_cycle < self.cycles:
            print(f"Starting cycle {self.cur_cycle} with " \
                  f"{len(self.geoms_to_kick)} geometries.")
            kicked_geoms = [self.get_kicked_geom() for _ in range(self.cycle_size)]
            # Write input geometries to disk
            self.geoms_to_trj(kicked_geoms, f"cycle_{self.cur_cycle:03d}_input.trj")
            opt_geoms = [self.run_kicked_geom(geom) for geom in kicked_geoms]
            # Filter out None and reject geoms where atoms are too close
            opt_geoms = [geom for geom in opt_geoms
                         if geom and not self.reject_by_distance(geom)]
            # opt_geoms = opt_geoms[not to_reject]
            opt_num = len(opt_geoms)
            print(f"Kicks in cycle {self.cur_cycle} produced "
                  f"{opt_num} new geometries.")

            # sorted_geoms = sorted(opt_geoms, key=lambda g: g.energy)
            # ens = [geom.energy for geom in sorted_geoms]
            # en_diffs = np.abs(np.diff(ens))
            # import pdb; pdb.set_trace()

            inds = list()
            for geom in opt_geoms:
                energy = geom.energy
                ind = bisect.bisect_left(self.energies, energy)
                try:
                    en_before= self.energies[ind-1]
                    if are_close(energy, en_before):
                        print("close energy before")
                        print(ind, energy, en_before)
                        continue
                except IndexError:
                    pass
                try:
                    en_after = self.energies[ind+1]
                    if are_close(energy, en_after):
                        print("close energy after")
                        print(ind, energy, en_after)
                        continue
                except IndexError:
                    pass
                bisect.insort(self.energies, energy)
                self.opt_geoms.insert(ind, geom)
                # if ind == 0 and self.energies and are_close(energy, self.energies[0]):
                    # print("low")
                    # continue
                # elif ind == len(self.energies) and are_close(energy, self.energies[-1]):
                    # print("high")
                    # continue
                # elif are_close(energy, self.energies[ind-1]) or are_close(energy, self.energies[ind+1]):
                    # print("middle")
                    # continue
                # else:
                    # bisect.insort(self.energies, energy)
                    # self.opt_geoms.insert(ind, geom)
            gens = [geom.energy for geom in self.opt_geoms]
            print("gens", gens)
            print("self.energies", self.energies)

            # import pdb; pdb.set_trace()
            trj_filtered_fn = f"cycle_{self.cur_cycle:03d}.trj"
            # Sort by energy
            opt_geoms = sorted(opt_geoms, key=lambda g: g.energy)
            self.geoms_to_trj(opt_geoms, trj_filtered_fn)

            self.cur_cycle += 1
        fn = "final.trj"
        self.geoms_to_trj(self.opt_geoms, fn)
        en_diffs = np.abs(np.diff(self.energies))
        print(en_diffs)
        assert en_diffs > 1e-4
