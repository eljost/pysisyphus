import numpy as np

from pysisyphus.Geometry import Geometry


class FreezingString:

    def __init__(self, images, calc_getter, max_nodes=10, opt_steps=3):
        self.images = images
        self.calc_getter = calc_getter
        self.opt_steps = opt_steps
        self.max_nodes = int(max_nodes)
        assert self.max_nodes % 2 == 0, "max_nodes must be a multiple of 2!"

        left_frontier, right_frontier = self.images
        self.left_string = [left_frontier, ]
        self.right_string = [right_frontier, ]
        self.opt_steps_left = self.opt_steps

        self._forces = None
        self.atoms = left_frontier.atoms
        coord_diff = np.linalg.norm(right_frontier.coords - left_frontier.coords)
        self.step_length = coord_diff / (self.max_nodes+1)
        self.set_new_frontier_nodes()
        self.coord_type = "cart"

    @property
    def left_frontier(self):
        return self.left_string[-1]

    @property
    def right_frontier(self):
        return self.right_string[0]

    def get_tangent(self):
        tangent = self.right_frontier.coords - self.left_frontier.coords
        tangent /= np.linalg.norm(tangent)
        return tangent

    @property
    def forces(self):
        left_forces = self.left_frontier.forces
        right_forces = self.right_frontier.forces
        forces = (left_forces, right_forces)
        tangent = self.get_tangent()
        perp_forces = np.array([f - f.dot(tangent)*tangent for f in forces]).flatten()
        self._forces = perp_forces
        self.energies = [self.left_frontier.energy, self.right_frontier.energy]
        
        return self._forces

    def as_xyz(self):
        return ""

    @property
    def fully_grown(self):
        return (len(self.left_string) + len(self.right_string)
                == self.max_nodes + 2)

    @property
    def energy(self):
        return 10
        # return self.energies

    def set_new_frontier_nodes(self):
        tangent = self.get_tangent()
        new_left_coords = self.left_frontier.coords + tangent*self.step_length
        new_left_frontier = Geometry(self.atoms, new_left_coords)
        new_left_frontier.set_calculator(self.calc_getter())
        self.left_string.append(new_left_frontier)

        new_right_coords = self.right_frontier.coords - tangent*self.step_length
        new_right_frontier = Geometry(self.atoms, new_right_coords)
        new_right_frontier.set_calculator(self.calc_getter())
        self.right_string.insert(0, new_right_frontier)
 
    def get_new_image(self, coords, index):
        new_image = Geometry(self.left_frontier.atoms, coords)
        new_image.set_calculator(self.calc_getter())
        self.images.insert(index, new_image)
        self.log(f"Create new image; insert it before index {index}.")
        return new_image

    def reparametrize(self):
        self.opt_steps_left -= 1

        if not self.fully_grown and self.opt_steps_left == 0:
            self.set_new_frontier_nodes()
            self.opt_steps_left = self.opt_steps

    @property
    def allcoords(self):
        return np.array([img.coords for img in self.left_string + self.right_string]).flatten()
    
    @property
    def coords(self):
        return np.array((self.left_frontier.coords, self.right_frontier.coords)).flatten()

    @coords.setter
    def coords(self, coords):
        left, right = coords.reshape(2, -1)
        self.left_frontier.coords = left
        self.right_frontier.coords = right
