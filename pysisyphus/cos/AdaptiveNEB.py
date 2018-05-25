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
        self.log(f"Updating adapt_thres. Old thresh was {old_thresh}. "
                 f"New threshold is {self.adapt_thresh:03f}")
        self.level += 1

    def adapt_this_cycle(self, forces):
        cur_rms = self.rms(forces)
        self.log(f"Current RMS of forces is {cur_rms:03f}.")
        return cur_rms <= self.adapt_thresh

    def prepare_opt_cycle(self, *args, **kwargs):
        super().prepare_opt_cycle(*args, **kwargs)

        if not self.adapt_thresh:
            self.update_adapt_thresh(self.forces_list[-1])

        if not self.adapt_this_cycle(self.forces_list[-1]):
            return

        hei_index = self.get_hei_index(self.energies_list[-1])
        if (hei_index == 0) or (hei_index == len(self.images)-1):
            raise Exception("Cant adapt, HEI is first or last!")

        prev_neighbour = hei_index - 1
        next_neighbour = hei_index + 1
        self.log(f"hei_index is {hei_index}")
        self.log("Adapting!")
        self.update_adapt_thresh(self.forces_list[-1])

        new_images = self.interpolate_between(prev_neighbour, next_neighbour,
                                              self.adapt_moving)
        assert len(new_images)+2 == len(self.images), "Need a way to add new calcs" \
                "to the images..."
        all_new_images = ([self.images[prev_neighbour]] + new_images
                          + [self.images[next_neighbour]])
        # Transfer coords to the old images.
        for i, image in enumerate(all_new_images):
            self.images[i].coords = image.coords
        print(new_images) 

        # This means the optimizer should reset itself
        return True
