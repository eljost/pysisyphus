import itertools

import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.patches import Circle
import numpy as np

class RFOPlotter():
    def __init__(self, calc, opt):
        self.opt = opt
        self.calc = calc

        self.xlim = calc.xlim
        self.ylim = calc.ylim
        self.coords = np.array(self.opt.coords)[:,:2]
        self.rfo_steps2 = np.array(opt.rfo_steps)[:,:2]
        self.rfo_steps = self.rfo_steps2 + self.coords

        self.fig, (self.ax, self.ax2) = plt.subplots(ncols=2)
        self.pause = True
        self.fig.canvas.mpl_connect('key_press_event', self.on_keypress)
        self.get_frame = itertools.cycle(range(self.opt.cur_cycle))

    def plot(self):
        x = np.linspace(*self.xlim, 100)
        y = np.linspace(*self.ylim, 100)
        X, Y = np.meshgrid(x, y)
        Z = np.full_like(X, 0)
        fake_atoms = ("X", )
        pot_coords = np.stack((X, Y, Z))
        pot = self.calc.get_energy(fake_atoms, pot_coords)["energy"]

        # Draw the contourlines of the potential
        levels = np.linspace(pot.min(), pot.max(), 25)
        contours = self.ax.contour(X, Y, pot, levels)

        # We only do 2d contour plots
        self.coord_lines, = self.ax.plot(*self.coords[0], "X-", label="Geometry")
        self.rfo_lines, = self.ax.plot(*self.rfo_steps[0], "x--", c="r",
                                       label="Pure RFO step")

        def quadratic_approx(cycle, step):
            E0 = self.opt.energies[cycle]
            g = -self.opt.forces[cycle]
            H = self.opt.hessians[cycle]
            return E0 + np.inner(g, step) + 0.5*step.dot(H.dot(step))

        self.lqa_pots = list()
        self.xs = np.linspace(*self.xlim, 100)
        self.ys = np.linspace(*self.ylim, 100)
        for cycle in range(self.opt.cur_cycle):
            cycle_pot = list()
            for x in self.xs:
                for y in self.ys:
                    step = np.array((x, y, 0))
                    E = quadratic_approx(cycle, step)
                    cycle_pot.append(E)
            cycle_pot = np.array(cycle_pot).reshape(-1, self.ys.size)
            self.lqa_pots.append(cycle_pot)
        self.lqa_pots = np.array(self.lqa_pots)
        self.lqa_contour = self.ax2.contour(self.xs, self.ys, self.lqa_pots[0])
        #import pdb; pdb.set_trace()

        circle_kwargs = {
            "fill": False,
            "color": "k",
        }
        #self.trust_region = Circle(self.coords[0], radius=self.opt.trust_radii[0],
        #                           **circle_kwargs)
        self.trust_region2 = Circle((0, 0), radius=self.opt.trust_radii[0],
                                   **circle_kwargs)
        #self.ax.add_patch(self.trust_region)
        self.ax2.add_patch(self.trust_region2)
        # Actual steps
        self.steps = np.array(self.opt.steps)[:,:2]
        act_steps1 = self.steps[0] + self.coords[0]
        self.actual_step_lines, = self.ax.plot(*act_steps1, "x", c="k",
                                               label="Actual step")
        self.actual_step_lines2, = self.ax2.plot(*self.steps[0], "x", c="k")

        self.rfo_lines2, = self.ax2.plot(*self.rfo_steps2[0], "x--", c="r",
                                       label="Pure RFO step")
        self.ax.legend()

    def update_plot(self, frame):
        self.fig.suptitle(f"Cycle {frame}")
        coords_x = self.coords[frame, 0]
        coords_y = self.coords[frame, 1]
        self.coord_lines.set_xdata(coords_x)
        self.coord_lines.set_ydata(coords_y)
        #plt.pause(1)

        # Draw alls steps
        rfo_x = self.rfo_steps[frame, 0]
        rfo_y = self.rfo_steps[frame, 1]
        # Draw only latest step
        #rfo_x = self.rfo_steps[frame, 0]
        #rfo_y = self.rfo_steps[frame, 1]
        self.rfo_lines.set_xdata(rfo_x)
        self.rfo_lines.set_ydata(rfo_y)

        self.lqa_contour.set_array(self.lqa_pots[frame])
        for tp in self.lqa_contour.collections:
            tp.remove()
        cycle_x, cycle_y = self.coords[frame]
        self.lqa_contour = self.ax2.contour(self.xs, self.ys, self.lqa_pots[frame])
        #self.lqa_contour = self.ax2.contour(self.xs+cycle_x,
        #                                    self.ys+cycle_y, self.lqa_pots[frame])
        #self.ax2.plot(0, 0, "X", c="r") # Highlight center of LQA potential
        #self.trust_region.center = (cycle_x, cycle_y)
        #self.trust_region.set_radius(self.opt.trust_radii[frame])
        self.trust_region2.set_radius(self.opt.trust_radii[frame])

        actual_step_x = self.steps[frame, 0]
        actual_step_y = self.steps[frame, 1]
        self.actual_step_lines.set_xdata(actual_step_x+cycle_x)
        self.actual_step_lines.set_ydata(actual_step_y+cycle_y)
        self.actual_step_lines2.set_xdata(actual_step_x)
        self.actual_step_lines2.set_ydata(actual_step_y)

        rfo2_x = self.rfo_steps2[frame, 0]
        rfo2_y = self.rfo_steps2[frame, 1]
        self.rfo_lines2.set_xdata(rfo2_x)
        self.rfo_lines2.set_ydata(rfo2_y)

    def animate(self):
        self.interval = 2000
        frames = range(self.opt.cur_cycle)
        self.animation = animation.FuncAnimation(self.fig,
                                                 self.update_plot,
                                                 frames=frames,
                                                 interval=self.interval)

    def on_keypress(self, event):
        """Advance the plot by one cycle (frame)."""
        if event.key == " ":
            frame = next(self.get_frame)
            self.update_plot(frame)
            plt.draw()

    #def on_keypress(self, event):
    #    """Pause on SPACE press."""
    #    #https://stackoverflow.com/questions/41557578
    #    if event.key == " ":
    #        if self.pause:
    #            self.animation.event_source.stop()
    #        else:
    #            self.animation.event_source.start()
    #        self.pause = not self.pause
