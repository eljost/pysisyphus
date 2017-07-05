#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

# [1] http://aip.scitation.org/doi/pdf/10.1063/1.1323224
# 
# https://github.com/cstein/neb/blob/master/neb/neb.py

def V(x, y):
    return 4 + 4.5*x - 4*y + x**2 + 2*y**2-2*x*y + x**4 - 2*x**2*y


def get_grad(x, y):
    return (
        4.5 + 2*x -2*y + 4*x**3 - 4*x*y,
        -4 + 4*y - 2*x - 2*x**2
    )


def interpolate_images(initial, final, images=10):
    step = (final-initial) / images
    x_initial, y_initial = initial
    x_final, y_final = final
    # Don't return initial and final image
    x_interpol = np.linspace(x_initial, x_final, num=images+2)[1:-1]
    y_interpol = np.linspace(y_initial, y_final, num=images+2)[1:-1]
    return np.column_stack((x_interpol, y_interpol))


def get_tangent(i, images):
    # [1], Eq. (2)
    prev_image = images[i-1]
    next_image = images[i+1]
    return (next_image-prev_image) / np.linalg.norm(next_image-prev_image)


def get_true_force(image, tangent):
    grad = np.array(get_grad(*image))
    return -grad + np.vdot(grad, tangent)


def get_spring_force(i, images, tangent, k=0.1):
    prev_image = images[i-1]
    image = images[i]
    next_image = images[i+1]
    #return k * np.vdot(((next_image-image)-(image-prev_image)), tangent) * tangent
    return k * (np.linalg.norm(next_image-image) - np.linalg.norm(image-prev_image)) * tangent


def steepest_descent(image, force, alpha=0.05):
    # Step along the force
    return image + alpha*force


def neb_step(all_images):
    grad_xs, grad_ys = get_grad(all_images[:, 0], all_images[:, 1])
    
    inner_indices = list(range(1, len(all_images)-1))
    tangents = np.array([get_tangent(i, all_images) for i in inner_indices])
    true_forces = np.array([get_true_force(all_images[i], tangent)
                   for i, tangent in zip(inner_indices, tangents)]
    )
    spring_forces = np.array([get_spring_force(i, all_images, tangent)
                    for i, tangent in zip(inner_indices, tangents)]
    )
    total_forces = spring_forces + true_forces
    new_images = ([steepest_descent(all_images[i], tf)
                           for i, tf in zip(inner_indices, total_forces)]
    )
    new_images = np.vstack((all_images[0].copy(), new_images, all_images[-1].copy()))

    return grad_xs, grad_ys, tangents, true_forces, spring_forces, total_forces, new_images


def plot_neb_step(old_images, new_images):
    fig, ax = plt.subplots(figsize=(width, height))

    # Potential
    levels = np.linspace(-4, 8, 20)
    contours = ax.contour(X, Y, Z, levels)
    ax.clabel(contours, inline=1, fontsize=10)

    img_xs = old_images[:,0]
    img_ys = old_images[:,1]

    # Image positions
    ax.plot(img_xs, img_ys, "ro", ls="-")

    # Gradient
    ax.quiver(img_xs, img_ys, grad_xs, grad_ys)

    # True force component
    #C = np.hypot(true_forces[:,0], true_forces[:,1])
    #ax.quiver(img_xs[1:-1], img_ys[1:-1], true_forces[:,0], true_forces[:,1], C)

    # Total force component
    C = np.hypot(total_forces[:,0], total_forces[:,1])
    ax.quiver(img_xs[1:-1], img_ys[1:-1], total_forces[:,0], total_forces[:,1], C)

    # New images
    ax.plot(new_images[:,0], new_images[:,1], "bx", ls="-")
    return fig


if __name__ == "__main__":
    # In[3]:

    width = 8
    height = width


    # In[4]:

    x = np.linspace(-2, 2.5, 100)
    y = np.linspace(0, 5, 100)
    X, Y = np.meshgrid(x, y)
    Z = V(X, Y)


    # In[5]:

    initial = np.array((-1.05274, 1.02776))
    final = np.array((1.94101, 3.85427))
    images = interpolate_images(initial, final, images=7)
    all_images = np.vstack((initial, images, final))


    old_images = all_images
    for i in range(1):
        grad_xs, grad_ys, tangents, true_forces, spring_forces, total_forces, new_images = neb_step(old_images)
        fig = plot_neb_step(old_images, new_images)
        #fig.savefig("iter{:02d}.png".format(i))
        old_images = new_images
        plt.show()

