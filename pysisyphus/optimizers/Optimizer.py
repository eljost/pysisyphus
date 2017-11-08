#!/usr/bin/env python3

import os
from pathlib import Path
import time

import numpy as np
import pandas as pd
import scipy as sp

from pysisyphus.cos.ChainOfStates import ChainOfStates


gauss_loose = {
    "max_force_thresh": 2.5e-3,
    "rms_force_thresh": 1.7e-3,
    "max_step_thresh": 1.0e-2,
    "rms_step_thresh": 6.7e-3
}


class Optimizer:

    def __init__(self, geometry, convergence=gauss_loose,
                 align=False, **kwargs):
        self.geometry = geometry

        self.convergence = convergence
        self.align = align
        for key, value in convergence.items():
            setattr(self, key, value)

        # Setting some default values
        self.max_cycles = 50
        self.max_step = 0.04
        self.rel_step_thresh = 1e-3
        self.keep_cycles = True
        self.out_dir = os.getcwd()

        assert(self.max_step > self.rel_step_thresh)

        self.is_cos = issubclass(type(self.geometry), ChainOfStates)
        self.is_zts = getattr(self.geometry, "reparametrize", None)

        # Overwrite default values if they are supplied as kwargs
        for key, value in kwargs.items():
            setattr(self, key, value)

        self.out_dir = Path(self.out_dir)
        if not self.out_dir.exists():
            os.mkdir(self.out_dir)
        self.cur_cycle = 0

        self.coords = list()
        self.energies = list()
        self.forces = list()
        self.steps = list()

        self.max_forces = list()
        self.rms_forces = list()
        self.max_steps = list()
        self.rms_steps = list()
        self.cycle_times = list()

        self.tangents = list()

    def procrustes(self):
        # http://nghiaho.com/?page_id=671#comment-559906
        image = self.geometry.images[0]
        coords3d = image.coords.reshape(-1, 3)
        centroid = coords3d.mean(axis=0)
        last_centered = coords3d - centroid
        image.coords = last_centered.flatten()

        rot_mats = [np.eye(3), ]
        for image in self.geometry.images[1:]:
            coords3d = image.coords.reshape(-1, 3)
            centroid = coords3d.mean(axis=0)
            # Center next image
            centered = coords3d - centroid
            mat = centered.T.dot(last_centered)
            U, W, Vt = np.linalg.svd(mat)
            rot_mat = U.dot(Vt)
            # Avoid reflections
            if np.linalg.det(rot_mat) < 0:
                U[:, -1] *= -1
                rot_mat = U.dot(Vt)
            # Rotate the coords
            rotated3d = centered.dot(rot_mat)
            image.coords = rotated3d.flatten()
            last_centered = rotated3d
            rot_mats.append(rot_mat)
        return rot_mats

    """
    def fit_rigid(self, vectors=(), hessian=None):
        rot_mats = self.procrustes()
        image_num = len(self.geometry.images)
        coord_len = len(self.geometry.images[0].coords)
        # Iterate over all vectors to be rotated
        results = list()
        for vec in vectors:
            # Reshape into a per image vectors
            rvec = vec.reshape(image_num, -1)
            rot = list()
            for rv, rm in zip(rvec, rot_mats):
                a = rv.dot(sp.linalg.block_diag(*([rm]*42)))
                #print(a.shape)
                rot.append(a)
            results.append(np.array(rot).flatten())
        return results

        if hessian:
            hessian = mat.dot(hessian).dot(mat.transpose())
            results = (velocities, hessian)
        return results
    """

    def fit_rigid_vector(self, vector, rot_mats):
        image_num = len(self.geometry.images)
        coord_len = len(self.geometry.images[0].coords)
        # Reshape into a per image vectors
        per_image = vector.reshape(image_num, -1, 3)
        rotated_vector = np.array(
            [vec.dot(mat) for vec, mat in zip(per_image, rot_mats)]
        ).flatten()
        """
        for rvec, mat in zip(vector.reshape(image_num, -1, 3), rot_mats):
            print(rvec.dot(mat).shape)
        for rv, rm in zip(rvec, rot_mats):
            a = rv.dot(sp.linalg.block_diag(*([rm]*42)))
            #print(a.shape)
            rot.append(a)
        results.append(np.array(rot).flatten())
        """
        return rotated_vector

    def fit_rigid(self, vectors=(), hessian=None):
        rot_mats = self.procrustes()
        rotated_vectors = [self.fit_rigid_vector(vec, rot_mats)
                           for vec in vectors]
        return rotated_vectors

    def check_convergence(self):
        # Only use forces perpendicular to the mep
        if self.is_cos:
            forces = self.geometry.perpendicular_forces
        else:
            forces = self.forces[-1]
        step = self.steps[-1]

        max_force = forces.max()
        rms_force = np.sqrt(np.mean(np.square(forces)))
        self.max_forces.append(max_force)
        self.rms_forces.append(rms_force)

        max_step = step.max()
        rms_step = np.sqrt(np.mean(np.square(step)))
        self.max_steps.append(max_step)
        self.rms_steps.append(rms_step)

        this_cycle = {
            "max_force_thresh": max_force,
            "rms_force_thresh": rms_force,
            "max_step_thresh": max_step,
            "rms_step_thresh": rms_step
        }

        self.is_converged = all(
            [this_cycle[key] <= getattr(self, key)
             for key in self.convergence.keys()
            ]
        )

    def print_header(self):
        hs = "max(force) rms(force) max(step) rms(step) s/cycle".split()
        header = "cycle" + " ".join([h.rjust(13) for h in hs])
        print(header)

    def print_convergence(self):
        int_fmt = "{:>5d}"
        float_fmt = "{:>12.6f}"
        conv_str = int_fmt + " " + (float_fmt + " ") * 4 + "{:>12.1f}"
        print(conv_str.format(
            self.cur_cycle, self.max_forces[-1], self.rms_forces[-1],
            self.max_steps[-1], self.rms_steps[-1], self.cycle_times[-1])
        )

    def scale_by_max_step(self, steps):
        steps_max = steps.max()
        if steps_max > self.max_step:
            steps *= self.max_step / steps_max
        return steps

    def optimize(self):
        raise Exception("Not implemented!")

    def write_to_out_dir(self, out_fn, content, mode="w"):
        out_path = self.out_dir / out_fn
        with open(out_path, mode) as handle:
            handle.write(content)

    def write_cycle_to_file(self):
        as_xyz_str = self.geometry.as_xyz()

        if self.is_cos:
            out_fn = "cycle_{:03d}.trj".format(self.cur_cycle)
            self.write_to_out_dir(out_fn, as_xyz_str)
        else:
            # Append to .trj file
            out_fn = "opt.trj"
            self.write_to_out_dir(out_fn, as_xyz_str+"\n", mode="a")

    def write_opt_data_to_file(self):
        energy_df = pd.DataFrame(self.energies)
        energy_df.to_csv(self.out_dir / "energies.csv", index=False)

    def run(self):
        self.print_header()
        while True:
            start_time = time.time()
            if self.cur_cycle == self.max_cycles:
                print("Number of cycles exceeded!")
                break

            self.coords.append(self.geometry.coords)

            steps = self.optimize()

            if steps is None:
                continue

            if self.is_cos:
                self.tangents.append(self.geometry.get_tangents())

            self.steps.append(steps)
            self.energies.append(self.geometry.energy)

            self.check_convergence()

            new_coords = self.geometry.coords + steps
            self.geometry.coords = new_coords

            self.cur_cycle += 1
            end_time = time.time()
            elapsed_seconds = end_time - start_time
            self.cycle_times.append(elapsed_seconds)

            if self.keep_cycles:
                self.write_cycle_to_file()
                self.write_opt_data_to_file()
            self.print_convergence()
            if self.is_converged:
                print("Converged!")
                print()
                break

            if self.is_zts:
                self.geometry.reparametrize()
