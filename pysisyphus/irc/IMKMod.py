import numpy as np

from pysisyphus.irc.IRC import IRC


# [1] http://pubs.acs.org/doi/pdf/10.1021/ja00295a002
# [2] https://math.stackexchange.com/questions/1882772/curve-fitting-with-derivatives


class IMKMod(IRC):

    def __init__(self, geometry, corr_first=True, corr_first_size=0.5,
                 corr_first_energy=True, corr_bisec_size=0.0025,
                 corr_bisec_energy=True, **kwargs):
        super().__init__(geometry, **kwargs)

        # Correction of pivot point
        self.corr_first = bool(corr_first)
        self.corr_first_size = float(corr_first_size)
        self.corr_first_energy = bool(corr_first_energy)
        # Correction along bisector
        self.corr_bisec_size = corr_bisec_size
        self.corr_bisec_energy = corr_bisec_energy

        self.corr_bisec_thresh = 2*self.corr_bisec_size

    def fit_parabola(self, x, y):
        y_str = np.array2string(np.array(y), precision=6)
        self.log(f"Fitting parabola with x={x} and y={y_str}")
        fit = np.polyfit(x, y, deg=2)
        parabola = np.poly1d(fit)
        minimum = parabola.deriv().r
        # We are only looking for a real minimum
        real_minimum = (minimum[minimum.imag==0].real)[0]
        self.log(f"Found minimum of fitted parabola at x={real_minimum:.4f}")

        return real_minimum

    def fit_grad_and_energies(self, energy0, grad0, energy1, direction):
        grad0_proj = grad0.dot(direction/np.linalg.norm(direction))
        # f (x) =  ax² + bx + c
        # f'(x) = 2ax  + b
        #
        # Assuming energy0 and grad0 are calculated at x=0
        # and
        # energy is calculated at x=1
        c = energy0
        b = grad0_proj
        a = energy1 - b - c
        coeffs = (a, b, c)
        poly = np.poly1d(coeffs)
        minimum = poly.deriv().r
        # We are only looking for a real minimum
        real_minimum = (minimum[minimum.imag==0].real)[0]
        self.log(f"Found minimum of fitted parabola at x={real_minimum:.4f}")

        return real_minimum

    def step(self):
        # Get gradient at Q0
        mw_coords_0 = self.mw_coords.copy()
        grad_0 = self.mw_gradient
        energy_0 = self.energy
        self.log(f"Current point:\n\tenergy_0={energy_0:.6f} au")


        ##############
        # PIVOT STEP #
        ##############

        grad_0_norm = np.linalg.norm(grad_0)
        dir_0 = -grad_0 / grad_0_norm
        # Step direction is against the gradient to the pivot point Q1.
        step = self.step_length * dir_0
        mw_coords_1 = self.mw_coords + step
        self.log(f"Took initial step of length {self.step_length:.6f} to pivot point")
        self.mw_coords = mw_coords_1
        energy_1 = self.energy

        #########################
        # PIVOT STEP CORRECTION #
        #########################

        self.log(f"Pivot point\n\tenergy_1={energy_1:.6f} au")
        energy_diff = energy_1 - energy_0
        if energy_diff > 0:
            self.log(f"Energy increased at pivot point ΔE={energy_diff:.6f}")
        correct_first_step = self.corr_first and (energy_diff > 0)
        # Fit using energy at q0, energy at q1, and energy at half step between
        # q0 and q1
        corr_min = None
        if correct_first_step and self.corr_first_energy:
            self.log("Computing correction for pivot point from energy_0, "
                     "energy_1, and additional energy calculation.")
            corr_step = self.corr_first_size*step
            mw_coords_1_corr = mw_coords_0 + corr_step
            self.mw_coords = mw_coords_1_corr
            energy_1_corr = self.energy
            x_corr = np.array((0, self.corr_first_size, 1))
            y_corr = (energy_0, energy_1_corr, energy_1)
            corr_min = self.fit_parabola(x_corr, y_corr)
        # Fit using energy & gradient at q0, and energy at q1
        elif correct_first_step:
            self.log("Computing correction for pivot point from energy_0, "
                     "grad_0 and energy_1, along dir_0")

            corr_min = self.fit_grad_and_energies(energy_0, grad_0, energy_1, dir_0)
        else:
            self.log("Skipping correction of pivot step")

        if corr_min:
            step_corr = corr_min * self.step_length * dir_0
            mw_coords_1 = mw_coords_0 + step_corr
            self.mw_coords = mw_coords_1

        # Get gradient at Q1
        self.log("Calculation of gradient at pivot point")
        grad_1 = self.mw_gradient
        energy_1 = self.energy
        grad_1_norm = np.linalg.norm(grad_1)
        dir_1 = -grad_1 / grad_1_norm

        ##############
        # BISECTOR D #
        ##############

        # Determine bisector d and normalized bisector vector D.
        # d = grad_0/grad_0_norm - grad_1/grad_1_norm
        d = -dir_0 + dir_1
        D = d / np.linalg.norm(d)
        self.log("Determined bisector D")

        ##########################
        # MINIMUM SEARCH ALONG D #
        ##########################

        # Take a step along D
        step_2 = self.corr_bisec_size * D
        mw_coords_2 = mw_coords_1 + step_2
        self.mw_coords = mw_coords_2
        energy_2 = self.energy
        self.log(f"1st step along D\n\tenergy_2={energy_2:.6f} au")

        if self.corr_bisec_energy:
            if energy_2 > energy_1:
                factor = 0.5
            else:
                factor = 2
            step_3 = self.corr_bisec_size * factor * D
            mw_coords_3 = mw_coords_1 + step_3
            self.mw_coords = mw_coords_3
            energy_3 = self.energy
            self.log(f"2nd step along D\n\tenergy_3={energy_3:.6f} au")
            x = np.array((0, 1, factor))
            y = (energy_1, energy_2, energy_3)
            corr_min = self.fit_parabola(x, y)
        else:
            self.log("Minimum search along D, by fitting energy_1, "
                     "grad_1 and energy_2.")
            corr_min = self.fit_grad_and_energies(energy_1, grad_1, energy_2, D)

        step_corr = self.corr_bisec_size * corr_min * D
        mw_coords_4 = mw_coords_1 + step_corr
        self.mw_coords = mw_coords_4
