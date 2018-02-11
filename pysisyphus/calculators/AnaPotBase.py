import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
from sympy import symbols, diff, lambdify, sympify

from pysisyphus.calculators.Calculator import Calculator

class AnaPotBase(Calculator):

    def __init__(self, V_str, xlim=(-1,1), ylim=(-1,1)): 
        super(AnaPotBase, self).__init__()
        self.xlim = xlim
        self.ylim = ylim
        x, y = symbols("x y")
        V = sympify(V_str)
        dVdx = diff(V, x)
        dVdy = diff(V, y)
        self.V = lambdify((x, y), V, "numpy")
        self.dVdx = lambdify((x, y), dVdx, "numpy")
        self.dVdy = lambdify((x, y), dVdy, "numpy")

        dVdxdx = diff(V, x, x)
        dVdxdy = diff(V, x, y)
        dVdydy = diff(V, y, y)

        self.dVdxdx = lambdify((x, y), dVdxdx, "numpy")
        self.dVdxdy = lambdify((x, y), dVdxdy, "numpy")
        self.dVdydy = lambdify((x, y), dVdydy, "numpy")

    def get_energy(self, atoms, coords):
        x, y, z = coords
        energy = self.V(x, y)
        return {"energy": energy}

    def get_forces(self, atoms, coords):
        x, y, z = coords
        dVdx = self.dVdx(x, y)
        dVdy = self.dVdy(x, y)
        dVdz = np.zeros_like(dVdx)
        forces = -np.array((dVdx, dVdy, dVdz))
        results = self.get_energy(atoms, coords)
        results["forces"] = forces
        return results

    def get_hessian(self, atoms, coords):
        x, y, z = coords
        dVdxdx = self.dVdxdx(x, y)
        dVdxdy = self.dVdxdy(x, y)
        dVdydy = self.dVdydy(x, y)
        hessian = np.array(((dVdxdx, dVdxdy, 0),
                            (dVdxdy, dVdydy, 0),
                            (0, 0, 0))
        )
        results = self.get_forces(atoms, coords)
        results["hessian"] = hessian
        return results

    def plot(self):
        self.fig, self.ax = plt.subplots()
        x = np.linspace(*self.xlim, 100)
        y = np.linspace(*self.ylim, 100)
        X, Y = np.meshgrid(x, y)
        Z = np.full_like(X, 0)
        fake_atoms = ("H", )
        pot_coords = np.stack((X, Y, Z))
        pot = self.get_energy(fake_atoms, pot_coords)["energy"]

        # Draw the contourlines of the potential
        levels = np.linspace(pot.min(), pot.max(), 25)
        contours = self.ax.contour(X, Y, pot, levels)

    def plot_opt(self, opt):
        self.plot()
        # We only do 2d contour plots
        coords = np.array(opt.coords)[:,:2]
        self.ax.plot(*coords.T, "X-")
        # Highlight start and end
        self.ax.plot(*coords[0], "X-", ms=10, c="r")
        self.ax.plot(*coords[-1], "X-", ms=10, c="r")
        for i in range(len(coords)):
            self.ax.annotate(f"{i}", xy=coords[i])

    def plot_rf_opt(self, opt):
        self.plot()
        # We only do 2d contour plots
        self.opt_steps = np.array(opt.coords)[:,:2]
        self.opt_step_lines, = self.ax.plot(*self.opt_steps[0], "X-")
        import pdb; pdb.set_trace()



        #self.images, = self.ax.plot(images_x, images_y, "ro", ls="-")
        #images_x = self.coords[frame][:,0]
        #images_y = self.coords[frame][:,1]
        #self.images.set_xdata(images_x)
        #self.images.set_ydata(images_y)

        ## Highlight start and end
        #self.ax.plot(*coords[0], "X-", ms=10, c="r")
        #self.ax.plot(*coords[-1], "X-", ms=10, c="r")
        #for i in range(len(coords)):
        #    self.ax.annotate(f"{i}", xy=coords[i])

    def anim_rfo(self, frame):
        opt_steps_x = self.opt_steps[:frame, 0]
        opt_steps_y = self.opt_steps[:frame, 1]
        self.opt_step_lines.set_xdata(opt_steps_x)
        self.opt_step_lines.set_ydata(opt_steps_y)
        print(f"frame {frame}")

    def animate(self, func, cycles):
        self.interval = 500
        frames = range(cycles)
        self.animation = animation.FuncAnimation(self.fig,
                                                 func,
                                                 frames=frames,
                                                 interval=self.interval)
