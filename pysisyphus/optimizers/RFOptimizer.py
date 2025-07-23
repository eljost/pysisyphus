# [1] http://aip.scitation.org/doi/10.1063/1.1515483 Optimization review
# [2] https://doi.org/10.1063/1.450914 Trust region method
# [3] 10.1007/978-0-387-40065-5 Numerical optimization
# [4] 10.1007/s00214-016-1847-3 Explorations of some refinements


import numpy as np

from pysisyphus.Geometry import Geometry
from pysisyphus.helpers_pure import rms
from pysisyphus.optimizers.HessianOptimizer import HessianOptimizer
from pysisyphus.optimizers.poly_fit import poly_line_search
from pysisyphus.optimizers.gdiis import gdiis, gediis


class RFOptimizer(HessianOptimizer):
    def __init__(
        self,
        geometry: Geometry,
        line_search: bool = True,
        gediis: bool = False,
        gdiis: bool = True,
        gdiis_thresh: float = 2.5e-3,
        gediis_thresh: float = 1e-2,
        gdiis_test_direction: bool = True,
        max_micro_cycles: int = 25,
        adapt_step_func: bool = False,
        **kwargs,
    ) -> None:
        """
        Rational function Optimizer.

        Parameters
        ----------
        geometry
            Geometry to be optimized.
        line_search
            Whether to carry out implicit line searches.
        gediis
            Whether to enable GEDIIS.
        gdiis
            Whether to enable GDIIS.
        gdiis_thresh
            Threshold for rms(forces) to enable GDIIS.
        gediis_thresh
            Threshold for rms(step) to enable GEDIIS.
        gdiis_test_direction
            Whether to the overlap of the RFO step and the GDIIS step.
        max_micro_cycles
            Number of restricted-step microcycles. Disabled by default.
        adapt_step_func
            Whether to switch between shifted Newton and RFO-steps.

        Other Parameters
        ----------------
        **kwargs
            Keyword arguments passed to the Optimizer/HessianOptimizer baseclass.
        """
        super().__init__(geometry, max_micro_cycles=max_micro_cycles, **kwargs)

        self.line_search = line_search
        self.gediis = gediis
        self.gdiis = gdiis
        self.gdiis_thresh = gdiis_thresh  # Will be compared to rms(step)
        self.gediis_thresh = gediis_thresh  # Will be compared to rms(forces)
        self.gdiis_test_direction = gdiis_test_direction
        self.adapt_step_func = adapt_step_func

        self.successful_gediis = 0
        self.successful_gdiis = 0
        self.successful_line_search = 0

    def get_step(self, energy, forces, hessian, eigvals, eigvecs, resetted):
        gradient = -forces
        step_func, pred_func = self.get_step_func(eigvals, gradient)

        ref_gradient = gradient.copy()
        # Reference step, used for judging the proposed GDIIS step
        ref_step = step_func(eigvals, eigvecs, gradient)

        # Right everything is in place to check for convergence.  If all values are below
        # the thresholds, there is no need to do additional inter/extrapolations.
        if self.check_convergence(ref_step)[0]:  # Drop conv_info
            self.log("Convergence achieved! Skipping inter/extrapolation.")
            return ref_step

        # Try to interpolate an intermediate geometry, either from GDIIS or line search.
        #
        # Set some defaults
        ip_gradient = None
        ip_step = None
        diis_result = None

        # Check if we can do GDIIS or GEDIIS. If we (can) do a line search is decided
        # after trying GDIIS.
        rms_forces = rms(gradient)
        rms_step = rms(ref_step)
        can_diis = (rms_step <= self.gdiis_thresh) and (not resetted)
        can_gediis = (rms_forces <= self.gediis_thresh) and (not resetted)

        # GDIIS / GEDIIS, prefer GDIIS over GEDIIS
        if self.gdiis and can_diis:
            # Gradients as error vectors
            err_vecs = -np.array(self.forces)
            diis_result = gdiis(
                err_vecs,
                self.coords,
                self.forces,
                ref_step,
                test_direction=self.gdiis_test_direction,
                logger=self.logger,
            )
            self.successful_gdiis += 1 if diis_result else 0
        # Don't try GEDIIS if GDIIS failed. If GEDIIS should be tried after GDIIS failed
        # comment the line below and uncomment the line following it.
        elif self.gediis and can_gediis:
            # if self.gediis and can_gediis and (diis_result == None):
            diis_result = gediis(
                self.coords,
                self.energies,
                self.forces,
                hessian=hessian,
                logger=self.logger,
            )
            self.successful_gediis += 1 if diis_result else 0

        try:
            ip_coords = diis_result.coords
            ip_step = ip_coords - self.geometry.coords
            ip_gradient = -diis_result.forces
        # When diis_result is None
        except AttributeError:
            self.log("GDIIS didn't succeed.")

        # Try line search if GDIIS failed or not requested
        if self.line_search and (diis_result is None) and (not resetted):
            ip_energy, ip_gradient, ip_step = poly_line_search(
                energy,
                self.energies[-2],
                gradient,
                -self.forces[-2],
                self.steps[-1],
                cubic_max_x=-1,
                quartic_max_x=2,
                logger=self.logger,
            )
            self.successful_line_search += 1 if ip_gradient is not None else 0

        # Use the interpolated gradient for the RFO step if interpolation succeeded
        if (ip_gradient is not None) and (ip_step is not None):
            gradient = ip_gradient
        # Keep the original gradient when the interpolation failed, but recreate
        # ip_step, as it will be returned as None from poly_line_search().
        else:
            ip_step = np.zeros_like(gradient)

        step = step_func(eigvals, eigvecs, gradient)
        # Form full step. If we did not interpolate or interpolation failed,
        # ip_step will be zero.
        step = step + ip_step

        # Use the original, actually calculated, gradient
        prediction = pred_func(ref_gradient, hessian, step)
        self.predicted_energy_changes.append(prediction)

        return step

    def postprocess_opt(self):
        msg = (
            f"Successful invocations:\n"
            f"\t     GEDIIS: {self.successful_gediis}\n"
            f"\t      GDIIS: {self.successful_gdiis}\n"
            f"\tLine Search: {self.successful_line_search}\n"
        )
        self.log(msg)
