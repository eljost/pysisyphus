from abc import abstractmethod
from typing import List, Optional, Protocol

from numpy.typing import NDArray


class CoordSys(Protocol):
    @abstractmethod
    def transform_forces(self, cart_forces: NDArray) -> NDArray:
        """Transform Cartesian forces to this coordinate system."""
        pass

    @abstractmethod
    def transform_hessian(
        self, cart_hessian: NDArray, int_gradient: Optional[NDArray]
    ) -> NDArray:
        """Transform Cartesian Hessian to this coordinate system."""
        pass

    @abstractmethod
    def transform_int_step(
        self, step: NDArray, update_constraints: bool = False, pure: bool = False
    ) -> NDArray:
        """Transform step in this coordinate system to Cartesians."""
        pass

    @abstractmethod
    def project_hessian(self, hessian: NDArray) -> NDArray:
        """Project Hessian in the current coordinate system."""
        pass

    @property
    @abstractmethod
    def coords(self) -> NDArray:
        """Getter for coordinates in this coordinate system."""
        pass

    @property
    @abstractmethod
    def coords3d(self) -> NDArray:
        """Getter for 3d Cartesian coordinates."""
        pass

    @coords3d.setter
    @abstractmethod
    def coords3d(self, coords3d: NDArray):
        """Setter for 3d Cartesian coordinates."""
        pass

    @property
    @abstractmethod
    def typed_prims(self) -> List:
        """List of (primitive) internal coordinates.

        May be empty, e.g., when the coordinate system is Cartesian."""
        pass
