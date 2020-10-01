class NeedNewInternalsException(Exception):
    def __init__(self, coords3d, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.coords3d = coords3d


class RebuiltInternalsException(Exception):
    pass


class DifferentPrimitivesException(Exception):
    pass
