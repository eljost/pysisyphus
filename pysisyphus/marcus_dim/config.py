from pysisyphus.constants import AU2EV


# Dummy energy offset to derive the energy of the first excited state
# for calculators, that only support GS calculations.
DUMMY_2EV = 2 / AU2EV
FIT_RESULTS_FN = "marcus_dim_fit.npz"
SCAN_RESULTS_FN = "marcus_dim_scan.npz"
# Default RMS threshold for the Marcus dimension
RMS_THRESH = 0.005
# Maximum bond length change in scan along the Marcus dimension
MAX_BOND_CHANGE = 0.4
# Default step length for scans along the Marcus dimension
SCAN_STEP_LENGTH = 0.01
