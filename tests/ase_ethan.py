from ase import io
from ase.calculators.orca import ORCA
from ase.optimize.fire import FIRE

#Optimise molecule
ethan = io.read('ethan.xyz')
orca_calc= ORCA(
    label="orca",
    maxiter=200,
    task="gradient",
    orcasimpleinput="HF-3c"
)
ethan.set_calculator(orca_calc)

opt = FIRE(ethan, trajectory='ethan.traj')
opt.run(fmax=0.00241)
