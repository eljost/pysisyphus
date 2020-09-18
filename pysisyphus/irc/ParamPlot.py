import matplotlib.pyplot as plt
import numpy as np


class ParamPlot:

    def __init__(self, coords_list, param1_inds, param2_inds):
        self.coords_list = coords_list
        self.coords_list_3d = [c.reshape(-1, 3) for c in coords_list]
        self.p1inds = param1_inds
        self.p2inds = param2_inds

        self.p1 = [self.get_param(c, param1_inds) for c in self.coords_list_3d]
        self.p2 = [self.get_param(c, param2_inds) for c in self.coords_list_3d]

        self.fig, self.ax = plt.subplots()

    def calc_bond(self, coords, ind1, ind2):
        return np.linalg.norm(coords[ind1] - coords[ind2])

    def calc_angle(self, coords, ind1, ind2, ind3):
        vec1 = coords[ind1] - coords[ind2]
        vec2 = coords[ind3] - coords[ind2]
        vec1n = np.linalg.norm(vec1)
        vec2n = np.linalg.norm(vec2)
        dotp = np.dot(vec1, vec2)
        radians = np.arccos(dotp / (vec1n * vec2n))
        return radians * 180 / np.pi

    def get_param(self, coords, param_inds):
        if len(param_inds) == 2:
            return self.calc_bond(coords, *param_inds)
        elif len(param_inds) == 3:
            return self.calc_angle(coords, *param_inds)
        else:
            raise Exception

    def plot(self):
        self.ax.plot(self.p1, self.p2, "ro", ls="--")
        self.ax.set_xlabel(str(self.p1inds))
        self.ax.set_ylabel(str(self.p2inds))

    def show(self):
        plt.tight_layout()
        plt.show()

    def save(self, path, prefix):
        out_fn = str(path / (prefix + ".pdf"))
        self.fig.tight_layout()
        self.fig.savefig(out_fn)

    def save_coords(self, path, prefix):
        np.savetxt(path / (prefix + ".coords"), self.coords_list)
