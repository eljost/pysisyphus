from pysisyphus.helpers import Geometry
from pysisyphus.helpers import geom_loader
from pysisyphus.benchmarks.data import *


class Benchmark:
    data_funcs = {
        "baker": get_baker_data,
        "baker_ts": get_baker_ts_data,
        "s22": get_s22_data,
        "zimmerman": get_zimmerman_data,
        "zimmerman_xtb": get_zimmerman_xtb_data,
        "birkholz_rx": get_birkholz_rx_data,
        "xtb_rx": get_xtb_rx_data,
        "precon_pos_rot": get_precon_pos_rot_data,
    }

    def __init__(
        self,
        name,
        exclude=None,
        inv_exclude=False,
        only=None,
        coord_type="cart",
        calc_getter=None,
    ):

        self.name = name
        if exclude is None:
            exclude = tuple()
        if only is None:
            only = tuple()
        self.exclude = exclude
        if isinstance(only, int):
            only = (only, )
        self.only = only
        self.coord_type = coord_type
        self.calc_getter = calc_getter

        self.data_func = self.data_funcs[self.name]
        self.prefix, self.data = self.data_func()
        self.fns = [fn for fn, *_ in self.data]

        def exclude_from(iterable):
            self.exclude = [
                id_ for id_, _ in enumerate(self.data) if id_ not in iterable
            ]

        # Only takes precedence
        if self.only:
            exclude_from(self.only)
        elif inv_exclude:
            exclude_from(self.exclude)

    def get_geoms(self, id_, set_calculator=True):
        fn, charge, mult, ref_energy = self.data[id_]

        geoms = geom_loader(self.prefix + fn, coord_type=self.coord_type)
        # Make atleast 1d
        if isinstance(geoms, Geometry):
                geoms = (geoms, )

        for geom in geoms:
            if set_calculator and self.calc_getter:
                geom.set_calculator(self.calc_getter(charge=charge, mult=mult))

        # Single geomtries are returned directly
        if len(geoms) == 1:
            geoms = geoms[0]
        return geoms

    def __iter__(self):
        self._id = 0
        return self

    def __next__(self):
        """Iterating directly over a Benchmark object, e.g., in
        pytest.mark.parametrize will give problems when calc_getter is set,
        but the calculator is not available/installed.

        Pytest will eagerly collect all geoms, and in this process the
        calculators will be created. But when the calculator command is not
        available, a SystemExit occur.

        In theses cases it may be better to leave the calc_getter and iterate
        over Benchmark.geom_iter, and to construct the calculators manually.
        """
        while self._id < len(self.data):
            if self._id in self.exclude:
                self._id += 1
                continue
            fn, *_, ref_energy = self.data[self._id]
            geoms = self.get_geoms(self._id)
            self._id += 1
            return fn, geoms, ref_energy
        raise StopIteration

    @property
    def geom_iter(self):
        for i, fn in enumerate(self.fns):
            if i in self.exclude:
                continue
            fn, charge, mult, ref_energy = self.data[i]
            geoms = self.get_geoms(i, set_calculator=False)
            yield fn, geoms, charge, mult, ref_energy

    @property
    def geoms(self):
        for i, fn in enumerate(self.fns):
            if i in self.exclude:
                continue
            yield self.get_geoms(i)
