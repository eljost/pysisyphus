#!/usr/bin/env python3

import numpy as np

from pysisyphus.cos.NEB import NEB


# [1] https://aip.scitation.org/doi/pdf/10.1063/1.1495401

# I need a way to add calculators to the new images ...
# or just reuse the old ones?

# See /scratch/projekte/biaryl/me_cn/cycloadd22/guess01/neb2

class AdaptiveNEB(NEB):

    def __init__(self, images, adapt_fact=.25, adapt_moving=3, **kwargs):
        kwargs["fix_ends"] = True
        super().__init__(images, **kwargs)

        self.adapt_fact = adapt_fact
        self.adapt_moving = adapt_moving

        self.adapt_thresh = None
        self.level = 1


    def rms(self, arr):
        return np.sqrt(np.mean(np.square(arr)))

    def update_adapt_thresh(self, forces):
        """Base on RMS forces."""
        old_thresh = self.adapt_thresh

        self.adapt_thresh = self.rms(forces) * self.adapt_fact
        self.log(f"Updating adapt_thres. Old thresh was {old_thresh:}. "
                 f"New threshold is {self.adapt_thresh:.03f}")
        self.level += 1

    def adapt_this_cycle(self, forces):
        cur_rms = self.rms(forces)
        adapt = cur_rms <= self.adapt_thresh
        self.log(f"Current RMS of forces is {cur_rms:03f}. Current thresh "
                 f"{self.adapt_thresh:03f}. Adapt = {adapt}")
        return adapt

    def prepare_opt_cycle(self, *args, **kwargs):
        super().prepare_opt_cycle(*args, **kwargs)

        # Check if track = True
        img = self.images[0]
        if hasattr(img, "track") and (img.track == True):
            raise Exception("track = True and interpolating new images "
                "may give problems with excited state tracking, so this "
                "is disabled for now.")

        if not self.adapt_thresh:
            self.update_adapt_thresh(self.forces_list[-1])

        if not self.adapt_this_cycle(self.forces_list[-1]):
            return

        hei_index = self.get_hei_index(self.energies_list[-1])
        if (hei_index == 0) or (hei_index == len(self.images)-1):
            raise Exception("Cant adapt, HEI is first or last!")

        prev_index = hei_index - 1
        next_index = hei_index + 1
        prev_image = self.images[prev_index]
        hei_image = self.images[hei_index]
        next_image = self.images[next_index]

        self.log(f"Index of highest energy image is {hei_index}")
        # Reset adapt_thresh so it will be set again in the beginning
        # of the next iteration.
        self.adapt_thresh = None

        new_images_1 = self.interpolate_between(prev_index, hei_index, 1)
        new_images_2 = self.interpolate_between(hei_index, next_index, 1)

        all_new_images = ([prev_image] + new_images_1
                          + [hei_image]
                          + new_images_2 + [next_image])
        assert len(all_new_images) == len(self.images), "This needs some " \
            "additional thought"

        # Backup old calculators
        calcs = [img.calculator for img in self.images]
        # Transfer calculators to the new images
        for new_image, calc in zip(all_new_images, calcs):
            new_image.set_calculator(calc)
        self.images = all_new_images
        return True
