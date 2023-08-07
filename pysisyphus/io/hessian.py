import numpy as np
import h5py

from pysisyphus.Geometry import Geometry
from pysisyphus.helpers_pure import eigval_to_wavenumber


def save_hessian(h5_fn, geom, cart_hessian=None, energy=None, mult=None, charge=None):
    if cart_hessian is None:
        cart_hessian = geom.cart_hessian

    if energy is None:
        energy = geom.energy

    if mult is None:
        mult = geom.calculator.mult

    if charge is None:
        charge = geom.calculator.charge

    if len(geom.atoms) > 1:
        proj_hessian = geom.eckart_projection(geom.mass_weigh_hessian(cart_hessian))
    else:
        proj_hessian = cart_hessian
    eigvals, _ = np.linalg.eigh(proj_hessian)
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
        handle.attrs["charge"] = charge


def save_third_deriv(h5_fn, geom, third_deriv_result, H_mw, H_proj):
    with h5py.File(h5_fn, "w") as handle:
        for key, value in third_deriv_result._asdict().items():
            handle.create_dataset(key, data=value)

        handle.create_dataset("coords3d", data=geom.coords3d)
        handle.create_dataset("masses", data=geom.masses)
        handle.create_dataset("H_mw", data=H_mw)
        handle.create_dataset("H_proj", data=H_proj)
        handle.attrs["atoms"] = [atom.lower() for atom in geom.atoms]


def geom_from_hessian(h5_fn: str, with_attrs: bool = False, **geom_kwargs):
    """Construct geometry from pysisyphus Hessian in HDF5 format.

    Parameters
    ----------
    h5_fn
        Filename of HDF5 Hessian.
    with_attrs
        Whether to also return an attributes dictionary. Attributes
        contain charge and multiplicity, as well as atoms and the electronic
        energy.

    Returns
    -------
    geom
        Geometry object with Hessian and electronic energy set.
    attrs
        Dictinoary containing the attributes set in the HDF5 file. Only returned
        when with_attrs is True.
    """
    with h5py.File(h5_fn, "r") as handle:
        coords3d = handle["coords3d"][:]
        cart_hessian = handle["hessian"][:]

        attrs = dict(handle.attrs.items())
        atoms = [atom.capitalize() for atom in attrs["atoms"]]
        energy = attrs["energy"]

    geom = Geometry(atoms=atoms, coords=coords3d, **geom_kwargs)
    geom.cart_hessian = cart_hessian
    geom.energy = energy

    if with_attrs:
        return geom, attrs
    else:
        return geom
