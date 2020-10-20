import numpy as np 
import h5py

from pysisyphus.helpers_pure import eigval_to_wavenumber


def save_hessian(h5_fn, geom, cart_hessian=None, energy=None, mult=None):
    if cart_hessian is None:
        cart_hessian = geom.cart_hessian

    if energy is None:
        energy = geom.energy

    if mult is None:
        mult = geom.calculator.mult

    if len(geom.atoms) > 1:
        hessian = geom.eckart_projection(geom.mass_weigh_hessian(cart_hessian))
    else:
        hessian = cart_hessian
    eigvals, eigvecs = np.linalg.eigh(hessian)
    vibfreqs = eigval_to_wavenumber(eigvals)

    masses = geom.masses
    atoms = geom.atoms
    coords3d = geom.coords3d

    with h5py.File(h5_fn, "w") as handle:
        handle.create_dataset("hessian", data=cart_hessian)
        handle.create_dataset("vibfreqs", data=vibfreqs)
        handle.create_dataset("masses", data=masses)
        handle.create_dataset("coords3d", data=coords3d)

        handle.attrs["atoms"] = [atom.lower() for atom in atoms]
        handle.attrs["energy"] = energy
        handle.attrs["mult"] = mult
