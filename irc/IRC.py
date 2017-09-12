#!/usr/bin/env python3

import copy
import logging

import matplotlib.pyplot as plt
import numpy as np

class IRC:

    def __init__(self, geometry, max_step=0.1, max_cycles=10,
                 only_forward=False, only_backward=False):
        assert(max_step > 0), "max_step has to be > 0"

        self.geometry = geometry
        self.max_step = max_step
        self.max_cycles = max_cycles
        self.forward = not only_backward
        self.backward = not only_forward

        self.hessian = self.geometry.hessian
        # Backup TS data
        self.ts_coords = copy.copy(self.geometry.coords)
        self.ts_energy = copy.copy(self.geometry.energy)
        self.ts_hessian = copy.copy(self.hessian)


        self.coords_list = list()


    def initial_displacement(self, energy_lowering=2.5e-4):
        mm_sqr_inv = self.geometry.mass_mat_sqr_inv
        eigvals, eigvecs = np.linalg.eig(self.geometry.mw_hessian)

        # Find smallest eigenvalue to get the imaginary mode
        logging.warning("Only considering smallest eigenvalue for now!")
        eigval_min = np.min(eigvals)
        img_index = np.where(eigvals == eigval_min)[0][0]
        logging.info(f"Smallest eigenvalue: {eigval_min}, index {img_index}")
        """
        # Zero small eigenvalues
        eigvals = np.abs(eigvals)
        all_indices = np.arange(eigvals.size)
        keep_indices = eigvals > 1e-4
        freqs = np.sqrt(eigvals[keep_indices])
        # We determine img_index before we cut the arrays. Imagine the inital
        # img_index would be 4, but after cutting the small eigenvalues there
        # are only 3 eigenvalues/eigenvectors left. Than the imaginary mode
        # couldn be accessed anymore.
        assert(img_index < freqs.size)
        # Flip sign on the imaginary eigenvalue
        freqs[img_index] *= -1
        print("Frequencies:")
        for i, f in enumerate(freqs):
            print(f"{i}: {f}")
        """

        # Calculate cartesian displacement vectors. Right now the eigenvectors
        # are mass-weighted.
        cart_displs = [mm_sqr_inv.dot(nm) for nm in eigvecs.transpose()]
        cart_displs = [cd/np.linalg.norm(cd) for cd in cart_displs]
        transition_vector = cart_displs[img_index]
        logging.info(f"Transition vector: {transition_vector}")
        """
        print("Cartesian displacement vectors, normalized:")
        for i, cd in enumerate(cart_displs):
            print(f"{i}: {cd}")
        """
        # Calculate the length of the initial step away from the TS to initiate
        # the IRC/MEP. We assume a quadratic potential and calculate the
        # displacement for a given energy lowering.
        # dE = (k*dq**2)/2 (dE = energy lowering, k = eigenvalue corresponding to
        # the transition vector/imaginary mode, dq = step length)
        # dq = sqrt(dE*2/k)
        # See 10.1021/ja00295a002 and 10.1063/1.462674
        step = np.sqrt(energy_lowering*2/np.abs(eigvals[img_index]))
        print(f"Inital step of {step:.4f} away from the TS.")
        # Forward and backward step from the TS
        """
        ts_coords = geom.coords
        forward_step = step*cart_displs[img_index]
        backward_step = -forward_step
        """
        return step*transition_vector

    def irc(self):
        last_energy = None
        irc_coords = [self.geometry.coords, ]
        irc_energies = []
        i = 0
        while True:
            if i == self.max_cycles:
                print("Cycles exceeded!")
                break

            logging.info(f"Step {i}")
            # Do macroiteration/IRC step
            self.step()
            irc_coords.append(self.geometry.coords)
            this_energy = self.geometry.energy
            irc_energies.append(this_energy)
            if last_energy and (this_energy > last_energy):
                print("Energy increased!")
                break
            elif last_energy and abs(last_energy - this_energy) <= 1e-4:
                print("Energy converged!")
                break
            last_energy = this_energy
            i += 1

        return irc_coords, irc_energies

    def run(self):
        self.coords = list()
        self.energies = list()
        init_displ = self.initial_displacement()
        logging.warning("Conversion angstroem/bohr!")

        if self.forward:
            logging.info("IRC forward")
            # Do inital displacement from the TS, forward
            self.geometry.coords = self.ts_coords + init_displ
            # Update hessian after initial displacement?
            self.geometry.hessian = self.ts_hessian
            forward_coords, forward_energies = self.irc()
            self.coords.extend(forward_coords[::-1])
            self.energies.extend(forward_energies[::-1])

        # Add TS data
        self.coords.append(self.ts_coords)
        self.energies.append(self.ts_energy)

        if self.backward:
            logging.info("IRC backward")
            # Do inital displacement from the TS, backward
            self.geometry.coords = self.ts_coords - init_displ
            # Update hessian after initial displacement?
            self.geometry.hessian = self.ts_hessian
            backward_coords, backward_energies = self.irc()
            self.coords.extend(backward_coords)
            self.energies.extend(backward_energies)

        self.coords = np.array(self.coords)
        self.energies = np.array(self.energies)
        self.postprocess()

    def postprocess(self):
        pass
