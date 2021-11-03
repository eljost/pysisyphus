from pysisyphus.intcoords.Primitive import Primitive
from pysisyphus.intcoords.Stretch import Stretch


class DistanceFunction(Primitive):
    def __init__(self, indices, coeff=-1, **kwargs):
        self.coeff = float(coeff)
        kwargs["calc_kwargs"] = ("coeff",)
        super().__init__(indices, **kwargs)

    @staticmethod
    def _weight(atoms, coords3d, indices, f_damping):
        return 1

    @staticmethod
    def _calculate(coords3d, indices, gradient=False, coeff=-1):
        """
        m--n  o--p
        """
        m, n, o, p = indices
        val1, grad1 = Stretch._calculate(coords3d, (m, n), gradient=True)
        val2, grad2 = Stretch._calculate(coords3d, (o, p), gradient=True)
        val = val1 + coeff * val2
        grad = grad1 + coeff * grad2
        if gradient:
            return val, grad
        return val

    @staticmethod
    def _jacobian(coords3d, indices, coeff):
        m, n, o, p = indices
        jac1 = Stretch._jacobian(coords3d, [m, n])
        jac2 = Stretch._jacobian(coords3d, [o, p])
        return jac1 + coeff * jac2
