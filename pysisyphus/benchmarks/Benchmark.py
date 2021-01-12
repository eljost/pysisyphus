from pysisyphus.helpers import geom_loader
from pysisyphus.benchmarks.data import get_baker_data, get_baker_ts_data, get_s22_data


class Benchmark:
    data_funcs = {
        "baker": get_baker_data,
        "baker_ts": get_baker_ts_data,
        "s22": get_s22_data,
    }

    def __init__(
        self, name, exclude=None, inv_exclude=False, coord_type="cart", calc_getter=None
    ):

        self.name = name
        if exclude is None:
            exclude = tuple()
        self.exclude = exclude
        self.coord_type = coord_type
        self.calc_getter = calc_getter

        self.data_func = self.data_funcs[self.name]
        self.prefix, self.data = self.data_func()
        self.fns = [fn for fn, *_ in self.data]

        if inv_exclude:
            self.exclude = [
                id_ for id_, _ in enumerate(self.data) if id_ not in self.exclude
            ]

    def get_geom(self, id_, set_calculator=True):
        fn, charge, mult, ref_energy = self.data[id_]

        geom = geom_loader(self.prefix + fn, coord_type=self.coord_type)
        if set_calculator and self.calc_getter:
            geom.set_calculator(self.calc_getter(charge=charge, mult=mult))
        return geom

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
            geom = self.get_geom(self._id)
            self._id += 1
            return fn, geom, ref_energy
        raise StopIteration

    @property
    def geom_iter(self):
        for i, fn in enumerate(self.fns):
            if i in self.exclude:
                continue
            fn, charge, mult, ref_energy = self.data[i]
            geom = self.get_geom(i, set_calculator=False)
            yield fn, geom, charge, mult, ref_energy
