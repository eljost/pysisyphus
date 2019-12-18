import numpy as np

def numhess(geom, step_size=0.0001):
    coords = geom.coords
    cnum = len(coords)
    H = list()
    for i in range(cnum):
        print(f"Step {i+1}/{cnum}")
        # Plus
        step = np.zeros_like(coords)
        step[i] += step_size
        pl_res = geom.get_energy_and_forces_at(coords+step)
        # print(pl_res)
        min_res = geom.get_energy_and_forces_at(coords-step)
        # print(min_res)

        gp = -pl_res["forces"]
        gm = -min_res["forces"]
        fd = (gp - gm) / (2*step_size)
        # print(i, fd)
        H.append(fd)
    H = np.array(H)
    # Symmetrize
    H = (H+H.T)/2
    return H


from pysisyphus.helpers import geom_from_library, eigval_to_wavenumber
from pysisyphus.calculators.XTB import XTB
from pysisyphus.calculators.ORCA import ORCA
from pysisyphus.calculators.AnaPot import AnaPot

# geom = geom_from_library("h2o2_hf_321g_opt.xyz")
# geom = geom_from_library("h2o.xyz")
geom = geom_from_library("hcn_bent.xyz")
# calc = XTB(pal=2)
# calc = ORCA(keywords="BP86 def2-SVP", pal=4)
calc = ORCA(keywords="hf sto-3g tightscf", pal=4)
geom.set_calculator(calc)
# geom = AnaPot.get_geom((-0.428, 0.981, 0.))
H = geom.hessian
print("ref hessian")
print(H)
nH = numhess(geom)
print("num hessian")
print(nH)
diff = np.abs(nH-H).sum()
print(diff)
rms = np.sqrt(np.mean((nH-H)**2))
print("rms", rms)


mw_nH = geom.mass_weigh_hessian(nH)
mwp_nH = geom.eckart_projection(mw_nH)
w, v = np.linalg.eigh(mwp_nH)
nus = eigval_to_wavenumber(w)
print(nus)


assert np.allclose(nH, H) 
# geom.calculator.plot(show=True)
# import matplotlib.pyplot as plt
# plt.show()
# numhes
# np.testing.assert_allclose(nH, H) 
# import pdb; pdb.set_trace()

# for i in range(1, 6):
    # step_size = 10**(-i)
    # nH = numhess(geom, step_size)
    # diff = np.abs(nH-H).sum()
    # try:
        # assert np.allclose(nH, H) 
        # print(f"{i:02d} ss={step_size:.4e} diff={diff:.4e}")
    # except AssertionError:
        # print(f"step_size={step_size:.4e} is not enough!")
