# [1] http://aip.scitation.org/doi/10.1063/1.1515483 Optimization review
# [2] https://doi.org/10.1063/1.450914 Trust region method
# [3] 10.1007/978-0-387-40065-5 Numerical optimization
# [4] 10.1007/s00214-016-1847-3 Explorations of some refinements


import numpy as np

from pysisyphus.helpers import rms
from pysisyphus.optimizers.HessianOptimizer import HessianOptimizer
from pysisyphus.optimizers.gdiis import gdiis, gediis


class RFOptimizer(HessianOptimizer):

    def __init__(self, geom, line_search=True, gediis=False, gdiis=True,
                 gdiis_thresh=2.5e-3, gediis_thresh=1e-2, max_micro_cycles=1,
                 *args, **kwargs):
        super().__init__(geom, max_micro_cycles=max_micro_cycles, *args, **kwargs)

        self.line_search = line_search
        self.gediis = gediis
        self.gdiis = gdiis
        self.gdiis_thresh = gdiis_thresh  # Will be compared to rms(step)
        self.gediis_thresh = gediis_thresh  # Will be compared to rms(forces)

    def optimize(self):
        energy, gradient, H, big_eigvals, big_eigvecs, resetted = self.housekeeping()

        # Reference RFO step, used for judging the proposed GDIIS step
        ref_gradient = gradient.copy()
        ref_rfo_step = self.get_rs_step(big_eigvals, big_eigvecs, gradient, name="RS-RFO")

        # Right everything is in place to check for convergence.  If all values are below
        # the thresholds, there is no need to do additional inter/extrapolations.
        if self.check_convergence(ref_rfo_step):
            self.log("Convergence achieved! Skipping inter/extrapolation.")
            return  ref_rfo_step


        # Try to interpolate an intermediate geometry, either from GDIIS or line search.
        #
        # Set some defaults
        ip_gradient = None
        ip_step = None
        diis_result = None

        # Check if we can do GDIIS or GEDIIS. If we (can) do a line search is decided
        # after trying GDIIS.
        rms_forces = rms(gradient)
        rms_step = rms(ref_rfo_step)
        can_diis = (rms_step <= self.gdiis_thresh) and (not resetted)
        can_gediis = (rms_forces <= self.gediis_thresh) and (not resetted)

        # GDIIS / GEDIIS, prefer GDIIS over GEDIIS
        if self.gdiis and can_diis:
            # Gradients as error vectors
            err_vecs = -np.array(self.forces)
            diis_result = gdiis(err_vecs, self.coords, self.forces, ref_rfo_step)
        # Don't try GEDIIS if GDIIS failed. If GEDIIS should be tried after GDIIS failed
        # comment the line below and uncomment the line following it.
        elif self.gediis and can_gediis:
        # if self.gediis and can_gediis and (diis_result == None):
            diis_result = gediis(self.coords, self.energies, self.forces, hessian=H)

        try:
            ip_coords = diis_result.coords
            ip_step = ip_coords - self.geometry.coords
            ip_gradient = -diis_result.forces
        # When diis_result is None
        except AttributeError:
            self.log("GDIIS didn't succeed.")

        # Try line search if GDIIS failed or not requested
        if self.line_search and (diis_result is None) and (not resetted):
            ip_energy, ip_gradient, ip_step = self.poly_line_search()

        # Use the interpolated gradient for the RFO step if interpolation succeeded
        if (ip_gradient is not None) and (ip_step is not None):
            gradient = ip_gradient
        # Keep the original gradient when the interpolation failed, but recreate
        # ip_step, as it will be returned as None from poly_line_search().
        else:
            ip_step = np.zeros_like(gradient)

        # RFO step (from intermediate geometry) with (interpolated) gradient
        rfo_step = self.get_rs_step(big_eigvals, big_eigvecs, gradient, name="RS-RFO")
        # Form full step. If we did not interpolate or it failed ip_step will be zero.
        step = rfo_step + ip_step

        # Use the original, actually calculated, gradient
        quadratic_prediction = step @ ref_gradient + 0.5 * step @ H @ step
        rfo_prediction = quadratic_prediction / (1 + step @ step)
        self.predicted_energy_changes.append(rfo_prediction)

        return step
