#!/usr/bin/env python3

#from linesearch import abls

def steepest_descent(geometry, alpha=-0.5, max_cycles=20):
    cycle = 0
    while cycle < max_cycles:
        forces = geometry.forces
        print("max(forces)", forces.max())
        if forces.max() < 0.001:
            break
        step = alpha*forces
        new_coords = geometry.coords + step
        geometry.coords = new_coords
        cycle += 1
