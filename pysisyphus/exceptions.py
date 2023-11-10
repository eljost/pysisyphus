class CalculationFailedException(Exception):
    def __init__(self, msg, path):
        super().__init__(msg)
        self.path = path


class HEIIsFirstOrLastException(Exception):
    pass
