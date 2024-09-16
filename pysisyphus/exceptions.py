class DifferentAtomOrdering(Exception):
    pass


class CalculationFailedException(Exception):
    def __init__(self, msg="", path=None):
        super().__init__(msg)
        self.path = path


class RunAfterCalculationFailedException(Exception):
    def __init__(self, msg="", err=None):
        super().__init__(msg)
        self.err = err


class HEIIsFirstOrLastException(Exception):
    pass
