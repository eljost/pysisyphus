# [1] https://doi.org/10.1002/jcc.27169
#     Nudged elastic stiffness band method: A method to solve kinks problems
#     of reaction paths
#     Mitsuta, Asada

# Adapted from the original implementation provided under MIT
# license by Yuki Mitsuta. Modifications by Johannes Steinmetzer.
#
# Main modifications were made in 'get_stiff_stress'; the other two
# functions are mostly original, disregarding some renaming.
#
# Original code is found at https://github.com/YukiMitsuta/Fstiffness


from collections.abc import Sequence
from typing import Union

import numpy as np


def get_stiff_stress(
    bandwidth: float,
    kappa: Union[float, Sequence],
    image_coords: np.ndarray,
    tangents: np.ndarray,
    straight_thresh: float = 1e-8,
) -> np.ndarray:
    """Function to calculate the stress of stiffness

    See [1] for a discussion.

    Parameters
    ----------
    bandwidth
        The band width of stiffness term in Bohr. The larger this parameter is,
        the higher the stress of stiffness.
    kappa
        The NEB spring constant (float or sequence of float). When a float is given
        the same spring constants is used for all images.
    image_coords
        2d array of image coordinates of shape (nimages, ncoords) with nimages
        being the number of images and ncoords the number of coordinates of
        each image the instance of each coordinates are numpy.array
    tangents
        2d array with the same shape as 'image_coords', containing the tangents
        defined at each image.

    Returns
    -------
    Fstlist
        Array containing the stress of stiffness with the same shape as 'image_coords'.
    """
    assert bandwidth > 0.0
    assert image_coords.shape == tangents.shape

    nimages, ncoords = image_coords.shape
    Fstlist = np.zeros((nimages, ncoords))

    # Determine if all images are in a straight line by calculating the dot products
    # of all tangents with the normalized coordinate difference vector between the first
    # and last image.
    start_end_tangent = image_coords[-1] - image_coords[0]
    start_end_tangent /= np.linalg.norm(start_end_tangent)
    ovlps = np.einsum("u,vu->v", start_end_tangent, tangents)
    if (np.abs(ovlps - 1.0) <= straight_thresh).all():
        return Fstlist

    tangents_perp, _ = calcEholo_vert(image_coords, tangents)

    offset = 0.5 * tangents_perp * bandwidth
    # Virtual image coordinates. Eq. (14) and (15) in [1].
    plus_coords = image_coords.copy() + offset
    minus_coords = image_coords.copy() - offset

    # Construct list of kappas, when only one kappa is given.
    # For N images, pysisyphus uses N-1 spring constants
    try:
        len(kappa)
    except TypeError:
        kappa = [kappa] * (nimages - 1)

    assert len(kappa) == (nimages - 1)

    # Terminal images don't feel stiffness force.
    for k in range(1, nimages - 1):
        Fstlist[k] = Fstiffness_k(
            k,
            kappa[k - 1],
            kappa[k],
            minus_coords,
            plus_coords,
            tangents_perp,
        )
    return Fstlist


def Fstiffness_k(
    k: int,
    kappa0: float,
    kappa1: float,
    innerimagelist: np.ndarray,
    outerimagelist: np.ndarray,
    Eholo_vertical: np.ndarray,
):
    """Stress of stiffness for image k."""

    v0inner = innerimagelist[k] - innerimagelist[k + 1]
    tdelta0inner = np.linalg.norm(v0inner)
    v0outer = outerimagelist[k] - outerimagelist[k + 1]
    tdelta0outer = np.linalg.norm(v0outer)
    tdeltadelta = tdelta0inner - tdelta0outer
    Vvert = -Eholo_vertical[k] * tdeltadelta * kappa0 * 0.5

    v1inner = innerimagelist[k] - innerimagelist[k - 1]
    tdelta1inner = np.linalg.norm(v1inner)
    v1outer = outerimagelist[k] - outerimagelist[k - 1]
    tdelta1outer = np.linalg.norm(v1outer)
    tdeltadelta = tdelta1inner - tdelta1outer
    Vvert -= Eholo_vertical[k] * tdeltadelta * kappa1 * 0.5
    return Vvert


def calcEholo_vert(
    image_coords: np.ndarray, tangents: np.ndarray, thresh: float = 1e-3
):
    """Calculate perpendicular tangents for given images and tangents.

    Implements eq. (8) to (13) in [1]. This function seems to crash when all images
    (tangents) lie (point) in a straight line. Currently, this check is done in
    'get_stiff_stress'.
    """
    nimages = len(image_coords)
    tangents_perp = [None for _ in range(nimages)]
    findvertlist = [False for _ in range(nimages)]
    whileN = 0
    while whileN < 1000:
        whileN += 1
        if all(findvertlist):
            break
        for k in range(nimages):
            if findvertlist[k]:
                continue
            tau = tangents[k]
            if k == 0 and findvertlist[1]:
                if findvertlist[1]:
                    tau_nei = tangents_perp[1]
                    a = -np.dot(tau_nei, tau) / np.dot(tau, tau)
                    tangents_perp[k] = a * tau + tau_nei
                    findvertlist[k] = True
            elif k == nimages - 1:
                if findvertlist[-2]:
                    tau_nei = tangents_perp[-2]
                    a = -np.dot(tau_nei, tau) / np.dot(tau, tau)
                    tangents_perp[k] = a * tau + tau_nei
                    findvertlist[k] = True
            else:
                v1 = image_coords[k - 1] - image_coords[k]
                v1 = v1 / np.linalg.norm(v1)
                v1taudot = np.abs(np.dot(v1, tau))
                v2 = image_coords[k + 1] - image_coords[k]
                v2 = v2 / np.linalg.norm(v2)
                v2taudot = np.abs(np.dot(v2, tau))
                if 1.0 - thresh < v1taudot and 1.0 - thresh < v2taudot:
                    if findvertlist[k - 1]:
                        tau_nei = tangents_perp[k - 1]
                        a = -np.dot(tau_nei, tau) / np.dot(tau, tau)
                        tangents_perp[k] = a * tau + tau_nei
                        findvertlist[k] = True
                    elif findvertlist[k + 1]:
                        tau_nei = tangents_perp[k + 1]
                        a = -np.dot(tau_nei, tau) / np.dot(tau, tau)
                        tangents_perp[k] = a * tau + tau_nei
                        findvertlist[k] = True
                elif thresh <= v1taudot:
                    a = -np.dot(v2, tau) / np.dot(v1, tau)
                    tangents_perp[k] = a * v1 + v2
                    tangents_perp[k] /= np.linalg.norm(tangents_perp[k])
                    findvertlist[k] = True
                elif thresh <= v2taudot:
                    a = -np.dot(v1, tau) / np.dot(v2, tau)
                    tangents_perp[k] = a * v2 + v1
                    tangents_perp[k] /= np.linalg.norm(tangents_perp[k])
                    findvertlist[k] = True
                # The elif is from the original code, but it should always evaluted to True.
                # If one of the two conditions wouldn't evalute to True, one of the two
                # conditionals above would run.
                # elif v1taudot < thresh and v2taudot < thresh:
                else:
                    tangents_perp[k] = v1
                    findvertlist[k] = True
                # else:
                # raise Exception("Error in calcEholovertTh")
    for k in range(1, nimages):
        vbefore = tangents_perp[k - 1]
        v = tangents_perp[k]
        if vbefore.dot(v) < 0.0:
            tangents_perp[k] *= -1
    tangents_perp = np.array(tangents_perp)
    return tangents_perp, findvertlist
