from pysisyphus.helpers import geom_loader
from pysisyphus.benchmarks.data import get_baker_data, get_baker_ts_data, get_s22_data


class Benchmark:
    data_funcs = {
        "baker": get_baker_data,
        "baker_ts": get_baker_ts_data,
        "s22": get_s22_data,
    }

    def __init__(self, name, exclude=None, coord_type="cart", calc_getter=None):

        self.name = name
        if exclude is None:
            exclude = tuple()
        self.exclude = exclude
        self.coord_type = coord_type
        self.calc_getter = calc_getter

        self.data_func = self.data_funcs[self.name]
        self.prefix, self.data = self.data_func()
        self.fns = [fn for fn, *_ in self.data]

    def get_geom(self, id_):
        fn, charge, mult, ref_energy = self.data[id_]

        geom = geom_loader(self.prefix + fn, coord_type=self.coord_type)
        if self.calc_getter:
            geom.set_calculator(self.calc_getter(charge=charge, mult=mult))
        return geom

    def __iter__(self):
        self._id = 0
        return self

    def __next__(self):
        i = self._id
        while i < len(self.data):
            if i in self.exclude:
                continue
            fn, *_, ref_energy = self.data[i]
            geom = self.get_geom(i)
            self._id += 1
            return fn, geom, ref_energy
        raise StopIteration

    @property
    def geoms(self):
        for i, fn in enumerate(self.fns):
            if i in self.exclude:
                continue
            geom = self.get_geom(i)
            yield geom

    @property
    def names_geoms(self):
        for i, fn in enumerate(self.fns):
            if i in self.exclude:
                continue
            geom = self.get_geom(i)
            yield fn, geom
