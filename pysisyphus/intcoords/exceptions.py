class NeedNewInternalsException(Exception):

    def __init__(
        self, coords3d, *args, invalid_inds=None, invalid_prims=None, **kwargs
    ):
        super().__init__(*args, **kwargs)

        self.coords3d = coords3d
        if invalid_inds is None:
            invalid_inds = ()
        self.invalid_inds = invalid_inds
        if invalid_prims is None:
            invalid_prims = ()
        self.invalid_prims = invalid_prims


class RebuiltInternalsException(Exception):

    def __init__(self, *args, typed_prims=None, **kwargs):
        super().__init__(*args, **kwargs)

        self.typed_prims = typed_prims


class DifferentPrimitivesException(Exception):
    pass


class DifferentCoordLengthsException(Exception):
    pass


class PrimitiveNotDefinedException(Exception):

    def __init__(self, typed_prim, *args, **kwargs):
        self.typed_prim = typed_prim
        super().__init__(*args, **kwargs)
