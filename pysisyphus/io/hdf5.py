import h5py


def init_h5_group(f, group_name, data_model):
    """Create group with given name and data model."""
    group = f.create_group(group_name)
    # Create (resizable) datasets by using None in maxshape
    for key, shape in data_model.items():
        assert len(shape) <= 2, "3D not yet supported"
        maxshape = (None, ) if (len(shape) == 1) else (None, shape[-1])
        group.create_dataset(key, shape, maxshape=maxshape)


def get_h5_group(fn, group_name, data_model=None, reset=False):
    """Return (and create if neccesary) group with given name and
    data model."""

    f = h5py.File(fn, mode="a")
    # Shortcut
    if data_model is None:
        return f[group_name]

    if group_name not in f:
        init_h5_group(f, group_name, data_model)
    group = f[group_name]

    # Check compatibility of data_model and group. If they aren't compatible
    # recreate the group with the proper shapes.
    try:
        compatible = [group[key].shape == shape for key, shape in data_model.items()]
    except KeyError:
        compatible = (False, )
    compatible = all(compatible)
    if (not compatible) or reset:
        del f[group_name]
        init_h5_group(f, group_name, data_model)
    group = f[group_name]

    return group


def resize_h5_group(group, max_cycles):
    """Increase size of first dimension of datasets in the given group."""
    for key, dataset in group.items():
        # No need to resize scalar datasets
        if dataset.shape == ():
            continue
        new_shape = list(dataset.shape).copy()
        new_shape[0] = max_cycles
        dataset.resize(new_shape)
