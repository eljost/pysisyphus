#!/usr/bin/env python3

#from linesearch import abls

def steepest_descent(image, force, alpha=-0.05):
    # Step against the force
    return image + alpha*force

"""
def sd(calculator):
    max_iter = 150
    hessian = np.identity
    for i in range(max_iter):
        forces = 
"""
