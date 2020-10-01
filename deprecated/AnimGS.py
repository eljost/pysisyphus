#!/usr/bin/env python3

import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np

class AnimGS:

    def __init__(self, gs, calc_getter):
        self.gs = gs
        self.calc = calc_getter()

        self.calc.plot()
        self.fig = self.calc.fig
        self.ax = self.calc.ax

        self.coords = [c.reshape(-1, 3) for c in self.gs.coords_list]
        self.tangents = self.gs.tangent_list
        self.perp_forces = self.gs.perp_force_list

        c0x = self.coords[0][:,0]
        c0y = self.coords[0][:,1]
        self.coord_lines, = self.ax.plot(c0x, c0y, "X-", ms=8, c="k")
        t0x = self.tangents[0][:,0]
        t0y = self.tangents[0][:,1]
        self.tang_quiv = self.ax.quiver(c0x, c0y, t0x, t0y)

    def update_func(self, frame):
        self.fig.suptitle(f"Cycle {frame}")

        coords = self.coords[frame]
        cx = self.coords[frame][:,0]
        cy = self.coords[frame][:,1]
        self.coord_lines.set_xdata(cx)
        self.coord_lines.set_ydata(cy)

        offsets = np.stack((cx, cy), axis=-1).flatten()

        tx = self.tangents[frame][:,0]
        ty = self.tangents[frame][:,1]
        self.tang_quiv.set_offsets(offsets)
        self.tang_quiv.set_UVC(tx, ty)

    def animate(self):
        frames = range(self.gs.cur_cycle)
        self.animation = animation.FuncAnimation(
                                    self.fig,
                                    self.update_func,
                                    frames=frames,
                                    interval=250,
        )
