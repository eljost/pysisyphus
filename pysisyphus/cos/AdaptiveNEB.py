#!/usr/bin/env python3

import numpy as np

from pysisyphus.cos.NEB import NEB


# [1] https://aip.scitation.org/doi/pdf/10.1063/1.1495401

# See /scratch/projekte/biaryl/me_cn/cycloadd22/guess01/neb2

class AdaptiveNEB(NEB):

    def __init__(self, images, adapt_fact=.25, adapt_between=1, **kwargs):
        """Adaptive Nudged Elastic Band.

        Parameters
        ----------
        images : list of Geometry objects
            Images of the band.
        adapt_fact : positive float
            Factor that is used to decide wether to adapt. The inital
            threshold is calculated by multiplying the RMS force of the
            band with this factor. When the RMS of the force falls
            below this threshold adaption takes place.
        adapat_between : positive integer
            Number of images to interpolate between the highest energy
            image and its neighbours. The number of starting images
            must be higher or equal then 2*adapt_between+3, as we
            reuse/transfer the calculators from the starting images onto
            the new ones.
        """
        kwargs["fix_ends"] = True
        super().__init__(images, **kwargs)

        self.adapt_fact = adapt_fact
        self.adapt_between = adapt_between

        self.adapt_thresh = None
        self.level = 1


    def rms(self, arr):
        """Root mean square

        Returns the root mean square of the given array.

        Parameters
        ----------
        arr : iterable of numbers

        Returns
        -------
        rms : float
            Root mean square of the given array.
        """
        return np.sqrt(np.mean(np.square(arr)))

    def update_adapt_thresh(self, forces):
        """Update the adaption threshold.

        Parameters
        ----------
        forces : np.array
            Forces of the previous optimization cycle.
        """
        old_thresh = self.adapt_thresh

        self.adapt_thresh = self.rms(forces) * self.adapt_fact
        self.log(f"Updating adapt_thres. Old thresh was {old_thresh:}. "
                 f"New threshold is {self.adapt_thresh:.03f}")

    def adapt_this_cycle(self, forces):
        """Decide wether to adapt.

        Parameters
        ----------
        forces : np.array
            Forces of the previous optimization cycle.

        Returns
        -------
        adapt : bool
            Flag that indicates if adaption should take place in this cycle.
        """
        cur_rms = self.rms(forces)
        adapt = cur_rms <= self.adapt_thresh
        self.log(f"Current RMS of forces is {cur_rms:03f}. Current thresh "
                 f"{self.adapt_thresh:03f}. Adapt = {adapt}")
        return adapt

    def prepare_opt_cycle(self, *args, **kwargs):
        """Check for adaption and adapt if needed.

        See ChainOfStates.prepare_opt_cycle for a complete docstring.
        """
        super().prepare_opt_cycle(*args, **kwargs)

        # Transferring Calculators including WFOWrapper objects
        # in excited state calculations may be problematic.
        # A problem would be the interpolation between two images with
        # different roots. What root should be taken for the interpolated
        # image? Additionally we probably would have to shift around the
        # iterations stored in the WFOWrapper.
        img0 = self.images[0]
        if hasattr(img0, "track") and (img0.track == True):
            raise Exception("track = True and interpolating new images "
                "may give problems with excited state tracking, so this "
                "is disabled for now.")

        # Initialize adapt_thresh
        if not self.adapt_thresh:
            self.update_adapt_thresh(self.forces_list[-1])

        if not self.adapt_this_cycle(self.forces_list[-1]):
            return

        # Adapation from here on
        # Determine highest energy index and image (HEI)
        hei_index = self.get_hei_index(self.energies_list[-1])
        self.log(f"Index of highest energy image is {hei_index}")
        if (hei_index == 0) or (hei_index == len(self.images)-1):
            self.log("Cant adapt, HEI is first or last!")
            return False
        prev_index = hei_index - 1
        next_index = hei_index + 1
        prev_image = self.images[prev_index]
        hei_image = self.images[hei_index]
        next_image = self.images[next_index]

        # Interpolation of new images between previous neighbour
        # and the HEI.
        new_images_1 = self.interpolate_between(prev_index, hei_index,
                                                self.adapt_between)
        # Between next neighbour and the HEI.
        new_images_2 = self.interpolate_between(hei_index, next_index,
                                                self.adapt_between)
        all_new_images = ([prev_image] + new_images_1
                          + [hei_image]
                          + new_images_2 + [next_image])
        assert len(all_new_images) <= len(self.images), "The number of new " \
            f"images ({len(all_new_images)}) is smaller than the number of " \
            f"current images ({len(self.images)}). Increase the number of " \
             "starting images or decrease 'adapt_between'."

        # Backup old calculators
        calcs = [img.calculator for img in self.images]
        # Transfer calculators to the new images
        for new_image, calc in zip(all_new_images, calcs):
            new_image.set_calculator(calc)
        self.images = all_new_images

        # Reset adapt_thresh so it will be set again in the beginning
        # of the next iteration.
        self.adapt_thresh = None
        self.level += 1
        return True
