import numpy as np

from pysisyphus.cos.ChainOfStates import ChainOfStates
from pysisyphus.Geometry import Geometry


class GrowingChainOfStates(ChainOfStates):
    def __init__(self, images, calc_getter, max_nodes=10, **kwargs):
        super().__init__(images, **kwargs)

        self.max_nodes = max_nodes
        self.calc_getter = calc_getter
        self.zero_step = np.zeros_like(self.images[0].coords)

    def get_new_image_from_coords(self, coords, index):
        new_image = Geometry(
            self.image_atoms,
            coords,
            coord_type=self.coord_type,
            coord_kwargs={"typed_prims": self.typed_prims},
        )
        new_image.set_calculator(self.calc_getter())
        self.images.insert(index, new_image)
        self.log(f"Create new image; insert it before index {index}.")
        return new_image

    @property
    def arc_dims(self):
        cds = [
            0,
        ]
        for i, image in enumerate(self.images[:-1]):
            next_image = self.images[i + 1]
            diff = np.linalg.norm(next_image - image)
            cds.append(diff)
        cds = np.cumsum(cds)
        tot_length = cds[-1]
        norm_cds = cds / cds.max()
        return tot_length, norm_cds

    @property
    def max_image_num(self):
        return self.max_nodes + 2

    def new_node_coords(self, k):
        l = (self.max_nodes - k) / (self.max_nodes + 1 - k)
        kth_coords = self.images[k].coords
        last_coords = self.images[-1].coords
        new_coords = l * kth_coords + (1 - l) * last_coords
        return new_coords

    def set_new_node(self, k):
        new_coords = self.new_node_coords(k)
        new_node = Geometry(self.image_atoms, new_coords)
        new_node.set_calculator(self.calc_getter())
        self.images.insert(k + 1, new_node)
        return new_node

    def prepare_opt_cycle(self, *args, **kwargs):
        parent_result, parent_messages = super().prepare_opt_cycle(*args, **kwargs)

        # Compare size of coords arrays to determine if new nodes
        # were added in the last reparametrization.
        last_size = self.coords_list[-1].size
        length_changed = last_size != self.coords.size
        reset_flag = parent_result or length_changed
        return reset_flag, parent_messages

    def propagate(self, force_unique_overlap_data_fns=True):
        # As the COS is not fully grown at the beginning calculating overlaps
        # between both sides of the COS may be difficult/impossible.
        # A way in supporting this for double-sided COS calucaltions may be to ask
        # the user for a second root that is used for the final image in the COS,
        # which is then propagated "inwards".
        raise Exception(
            "propagate() is not yet implemented for growing Chain-Of-States!"
        )
