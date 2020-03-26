import abc


class Primitive(metaclass=abc.ABCMeta):

    def __init__(self, indices):
        self.indices = indices

    @abc.abstractmethod
    def _calculate(*, coords3d, indices, gradient):
        pass

    def calculate(self, coords3d, indices=None, gradient=False):
        if indices is None:
            indices = self.indices

        return self._calculate(
                    coords3d=coords3d,
                    indices=indices,
                    gradient=gradient
        )
