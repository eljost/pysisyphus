from typing import Protocol


class CoordSys(Protocol):

    def transform_forces(self, cart_forces):
        pass
    
    def transform_hessian(self, cart_hessian, cart_gradient):
        pass

    def transform_int_step(self, step, **kwargs):
        pass

    def project_hessian(self, hessian):
        pass

    @property
    def coords(self):
        """Potentially modified coordinates."""
        pass

    @property
    def coords3d(self):
        """Original Cartesians as 2-dimensional array. """
        pass

    @property
    def typed_prims(self):
        pass
