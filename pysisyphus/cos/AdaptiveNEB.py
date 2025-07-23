# [1] https://aip.scitation.org/doi/pdf/10.1063/1.1495401
# [2] http://dx.doi.org/10.1063/1.4962019
#     Zhang, 2016
#     FreeEnd Adaptive NEB

# See /scratch/projekte/biaryl/me_cn/cycloadd22/guess01/neb2

import numpy as np

from pysisyphus.cos.NEB import NEB
from pysisyphus.interpolate import interpolate


class AdaptiveNEB(NEB):
    def __init__(
        self,
        images,
        adapt=True,
        adapt_fact=0.25,
        adapt_between=1,
        scale_fact=False,
        keep_hei=True,
        free_ends=True,
        **kwargs,
    ):
        """(Free-End) Adaptive Nudged Elastic Band.

        Parameters
        ----------
        images : list of Geometry objects
            Images of the band.
        adapt : bool, default True
            Whether to adapt the image number or not. This switch is included
            to support the FreeEndNEB class, that is just a thin wrapper around
            this class.
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
        scale_fact : bool, default False
            Whether to increase adapt_fact in deeper levels. This may lead
            to earlier adapation.
        keep_hei : bool, optional
            Whether to keep the highest energy image (usually a very good
            idea) or to interpolate only between the neighbouring images.
        free_ends : bool, default True
            Whether to use modified forces on the end images.
        """
        super().__init__(images, **kwargs)

        self.adapt = adapt
        self.adapt_fact = adapt_fact
        self.adapt_between = adapt_between
        self.scale_fact = scale_fact
        self.keep_hei = keep_hei
        self.free_ends = free_ends

        self.adapt_thresh = None
        self.level = 1
        self.coords_backup = list()

    @NEB.forces.getter
    def forces(self):
        """See Eq. (7) in [2]."""

        forces = super().forces
        forces_size = self.images[-1].forces.size

        if self.free_ends and (not self.fix_first):
            mod_forces = self.get_perpendicular_forces(0)
            forces[:forces_size] = mod_forces

        if self.free_ends and (not self.fix_last):
            mod_forces = self.get_perpendicular_forces(self.last_index)
            forces[-forces_size:] = mod_forces

        self._forces = forces
        return self._forces

    def update_adapt_thresh(self, forces):
        """Update the adaption threshold.

        Parameters
        ----------
        forces : np.array
            Forces of the previous optimization cycle.
        """
        old_thresh = self.adapt_thresh

        # Dividing by (1 / level)**1/2 scales the adapt_fact as
        # level 1: 1. (No scaling)
        # level 2: 1.414
        # level 3: 1.732
        # level 4: 2.
        # level 5: 2.236
        # ...
        # This ensures that the adapt_fact increases as we recurse
        # deeper.
        if self.scale_fact:
            self.adapt_thresh = (
                self.rms(forces) * self.adapt_fact / np.sqrt(1 / self.level)
            )
        else:
            self.adapt_thresh = self.rms(forces) * self.adapt_fact
        # arr / np.sqrt(1/np.arange(1, 6))
        self.log(
            f"Updating adapt_thres. Old thresh was {old_thresh:}. "
            f"New threshold is {self.adapt_thresh:.03f}"
        )

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
        self.log(
            f"Current RMS of forces is {cur_rms:03f}. Current thresh "
            f"{self.adapt_thresh:03f}. Adapt = {adapt}"
        )
        return adapt

    def prepare_opt_cycle(self, *args, **kwargs):
        """Check for adaption and adapt if needed.

        See ChainOfStates.prepare_opt_cycle for a complete docstring.
        """
        base_reset, base_messages = super().prepare_opt_cycle(*args, **kwargs)
        if not self.adapt:
            return base_reset, base_messages

        # Transferring Calculators including WFOWrapper objects
        # in excited state calculations may be problematic.
        # A problem would be the interpolation between two images with
        # different roots. What root should be taken for the interpolated
        # image? Additionally we probably would have to shift around the
        # iterations stored in the WFOWrapper.
        img0 = self.images[0]
        if hasattr(img0, "track") and (img0.track == True):
            raise Exception(
                "track = True and interpolating new images "
                "may give problems with excited state tracking, so this "
                "is disabled for now."
            )

        # Initialize adapt_thresh
        if not self.adapt_thresh:
            self.update_adapt_thresh(self.forces_list[-1])

        if not self.adapt_this_cycle(self.forces_list[-1]):
            return base_reset, base_messages

        #
        # Adapation from here on
        #
        # Backup coords if we have to step back
        self.coords_backup.append(self.coords)
        messages = base_messages
        # Determine highest energy index and image (HEI)
        hei_index = self.get_hei_index(self.all_energies[-1])
        self.log(f"Index of highest energy image is {hei_index}")
        if (hei_index == 0) or (hei_index == len(self.images) - 1):
            self.log("Cant adapt, HEI is first or last!")
            return base_reset
        else:
            self.fix_first = False
            self.fix_last = False
            self.log("First and last image are now free to move.")
        prev_index = hei_index - 1
        next_index = hei_index + 1
        prev_image = self.images[prev_index]
        hei_image = self.images[hei_index]
        next_image = self.images[next_index]

        # Two interpolations
        #   prev. neighbour - HEI
        #   HEI - next neighbour
        # Usually the better idea
        kwargs = {
            "kind": "lst",
            "only_between": True,
            "interpol_kwargs": {"silent": True},
        }
        if self.keep_hei:
            # Interpolation of new images between previous neighbour
            # and the HEI.
            kwargs["between"] = self.adapt_between
            new_images_1 = interpolate(prev_image, hei_image, **kwargs)
            new_images_2 = interpolate(hei_image, next_image, **kwargs)
            # Between next neighbour and the HEI.
            all_new_images = (
                [prev_image] + new_images_1 + [hei_image] + new_images_2 + [next_image]
            )
        # One interpolation
        #   prev. neighbour - next neighbour
        else:
            kwargs["between"] = 2 * self.adapt_between + 1
            new_images = interpolate(prev_image, next_image, **kwargs)
            all_new_images = [prev_image] + new_images + [next_image]

        assert len(all_new_images) <= len(self.images), (
            f"The number of new images ({len(all_new_images)}) is smaller than "
            f"the number of current images ({len(self.images)}). Increase the number "
            "of starting images or decrease 'adapt_between'."
        )
        self.level += 1

        messages.append(
            f"Adapted images! New number of images is {len(all_new_images)}. "
            f"Current level is {self.level}."
        )

        # Backup old calculators
        calcs = [img.calculator for img in self.images]
        # Get numbered indices of the current calcs, so we can continue with
        # higher numbers for the new images.
        calc_numbers = [calc.calc_number for calc in calcs]
        new_calc_number_start = max(calc_numbers) + 1
        # Transfer calculators to the new images
        for i, (new_image, calc) in enumerate(zip(all_new_images, calcs)):
            calc.calc_number = new_calc_number_start + i
            new_image.set_calculator(calc)
        self.images = all_new_images
        self.set_zero_forces_for_fixed_images()

        # Reset adapt_thresh so it will be set again in the beginning
        # of the next iteration.
        self.adapt_thresh = None
        return True, messages
