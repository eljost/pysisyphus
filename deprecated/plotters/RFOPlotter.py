import itertools

import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.patches import Circle
import numpy as np

class RFOPlotter():
    def __init__(self, calc, opt, figsize=(8, 6),
                 save=False, title=True):
        self.opt = opt
        self.calc = calc
        self.figsize = figsize
        self.save = save
        self.title = title

        self.xlim = calc.xlim
        self.ylim = calc.ylim
        self.coords = np.array(self.opt.coords)[:,:2]
        self.rfo_steps2 = np.array(opt.rfo_steps)[:,:2]
        self.rfo_steps = self.rfo_steps2 + self.coords

        self.fig, (self.ax, self.ax2) = plt.subplots(ncols=2,
                                                     figsize=figsize)
        self.ax.set_title("True potential")
        self.ax2.set_title("Local quadratic approximation")
        self.pause = True
        self.fig.canvas.mpl_connect('key_press_event', self.on_keypress)
        self.get_frame = itertools.cycle(range(self.opt.cur_cycle))

    def plot(self):
        # Draw potential as contour lines
        self.xs = np.linspace(*self.xlim, 100)
        self.ys = np.linspace(*self.ylim, 100)
        X, Y = np.meshgrid(self.xs, self.ys)
        Z = np.full_like(X, 0)
        fake_atoms = ("X", )
        pot_coords = np.stack((X, Y, Z))
        pot = self.calc.get_energy(fake_atoms, pot_coords)["energy"]
        levels = np.linspace(pot.min(), pot.max(), 25)
        contours = self.ax.contour(X, Y, pot, levels)

        # Calculate the LQA potentials for every cycle (frame)
        def quadratic_approx(cycle, step):
            E0 = self.opt.energies[cycle]
            g = -self.opt.forces[cycle]
            H = self.opt.hessians[cycle]
            return E0 + np.inner(g, step) + 0.5*step.dot(H.dot(step))

        self.lqa_pots = list()
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
        # Draw LQA potential as contour lines
        self.lqa_contour = self.ax2.contour(self.xs, self.ys, self.lqa_pots[0])

        # Draw the actual geometries
        self.coord_lines, = self.ax.plot(*self.coords[0], "X",
                                         label="Geometry")
        # Draw the pure RFO steps
        rfo_kwargs = {
            "marker": "x",
            "c": "r",
            "label": "Pure RFO step",
        }
        self.rfo_lines, = self.ax.plot(*self.rfo_steps[0], **rfo_kwargs)
        self.rfo_lines2, = self.ax2.plot(*self.rfo_steps2[0], **rfo_kwargs)

        # Draw the actual steps taken
        actual_kwargs = {
            "marker": "x",
            "c": "k",
            "label": "Actual step",
        }
        self.steps = np.array(self.opt.steps)[:,:2]
        act_steps1 = self.steps[0] + self.coords[0]
        self.actual_step_lines, = self.ax.plot(*act_steps1, **actual_kwargs)
        self.actual_step_lines2, = self.ax2.plot(*self.steps[0], **actual_kwargs)


        circle_kwargs = {
            "fill": False,
            "color": "k",
        }
        #self.trust_region = Circle(self.coords[0], radius=self.opt.trust_radii[0],
        #                           **circle_kwargs)
        #self.ax.add_patch(self.trust_region)
        self.trust_region2 = Circle((0, 0), radius=self.opt.trust_radii[0],
                                   **circle_kwargs)
        self.ax2.add_patch(self.trust_region2)
        self.ax2.plot(0, 0, "X")

        self.ax.legend()

    def update_plot(self, frame):
        if self.title:
            self.fig.suptitle(f"Cycle {frame}")

        # Update the geometry in ax
        coords_x = self.coords[frame, 0]
        coords_y = self.coords[frame, 1]
        self.coord_lines.set_xdata(coords_x)
        self.coord_lines.set_ydata(coords_y)

        # Update LQA contours in ax2
        self.lqa_contour.set_array(self.lqa_pots[frame])
        for tp in self.lqa_contour.collections:
            tp.remove()
        cycle_x, cycle_y = self.coords[frame]
        self.lqa_contour = self.ax2.contour(self.xs, self.ys, self.lqa_pots[frame])
        #self.trust_region.center = (cycle_x, cycle_y)
        #self.trust_region.set_radius(self.opt.trust_radii[frame])

        # Update trust region in ax2
        self.trust_region2.set_radius(self.opt.trust_radii[frame])

        # Update the actual steps taken
        actual_step_x = self.steps[frame, 0]
        actual_step_y = self.steps[frame, 1]
        self.actual_step_lines.set_xdata(actual_step_x+cycle_x)
        self.actual_step_lines.set_ydata(actual_step_y+cycle_y)
        self.actual_step_lines2.set_xdata(actual_step_x)
        self.actual_step_lines2.set_ydata(actual_step_y)

        # Update the (potentially bigger than allowed) RFO steps
        rfo_x = self.rfo_steps[frame, 0]
        rfo_y = self.rfo_steps[frame, 1]
        self.rfo_lines.set_xdata(rfo_x)
        self.rfo_lines.set_ydata(rfo_y)
        rfo2_x = self.rfo_steps2[frame, 0]
        rfo2_y = self.rfo_steps2[frame, 1]
        self.rfo_lines2.set_xdata(rfo2_x)
        self.rfo_lines2.set_ydata(rfo2_y)
        plt.tight_layout()

        if self.save:
            frame_fn = f"step{frame}.png"
            self.fig.savefig(frame_fn)
            # if not os.path.exists(frame_fn):
                # self.fig.savefig(frame_fn)

    # def animate(self):
        # self.interval = 2000
        # frames = range(self.opt.cur_cycle)
        # self.animation = animation.FuncAnimation(self.fig,
                                                 # self.update_plot,
                                                 # frames=frames,
                                                 # interval=self.interval)

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
