from matplotlib.patches import Circle
import matplotlib.pyplot as plt
import numpy as np
import pytest

from pysisyphus.calculators.AnaPot import AnaPot
from pysisyphus.dynamics.velocity_verlet import md


def test_velocity_verlet():
    geom = AnaPot.get_geom((0.52, 1.80, 0))
    x0 = geom.coords.copy()
    v0 = .1 * np.random.rand(*geom.coords.shape)
    t = 3
    dts = (.005, .01, .02, .04, .08)
    all_xs = list()
    for dt in dts:
        geom.coords = x0.copy()
        md_kwargs = {
            "v0": v0.copy(),
            "t": t,
            "dt": dt,
        }
        md_result = md(geom, **md_kwargs)
        all_xs.append(md_result.coords)
    calc = geom.calculator
    calc.plot()
    ax = calc.ax
    for dt, xs in zip(dts, all_xs):
        ax.plot(*xs.T[:2], "o-", label=f"dt={dt:.3f}")
        # ax.plot(*xs.T[:2], "-", label=f"dt={dt:.3f}")
    ax.legend()
    plt.show()


def ase_md_playground():
    geom = AnaPot.get_geom((0.52, 1.80, 0), atoms=("H", ))
    atoms = geom.as_ase_atoms()
    # ase_calc = FakeASE(geom.calculator)
    # from ase.optimize import BFGS
    # dyn = BFGS(atoms)
    # dyn.run(fmax=0.05)

    import ase
    from ase import units
    from ase.io.trajectory import Trajectory
    from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
    from ase.md.verlet import VelocityVerlet

    MaxwellBoltzmannDistribution(atoms, 300 * units.kB)
    momenta = atoms.get_momenta()
    momenta[0, 2] = 0.
    # Zero 3rd dimension
    atoms.set_momenta(momenta)

    dyn = VelocityVerlet(atoms, .005 * units.fs)  # 5 fs time step.


    def printenergy(a):
        """Function to print the potential, kinetic and total energy"""
        epot = a.get_potential_energy() / len(a)
        ekin = a.get_kinetic_energy() / len(a)
        print('Energy per atom: Epot = %.3feV  Ekin = %.3feV (T=%3.0fK)  '
              'Etot = %.3feV' % (epot, ekin, ekin / (1.5 * units.kB), epot + ekin))

    # Now run the dynamics
    printenergy(atoms)
    traj_fn = 'asemd.traj'
    traj = Trajectory(traj_fn, 'w', atoms)
    dyn.attach(traj.write, interval=5)
    # dyn.attach(bumms().bimms, interval=1)

    dyn.run(10000)
    printenergy(atoms)
    traj.close()

    traj = ase.io.read(traj_fn+"@:")#, "r")
    pos = [a.get_positions() for a in traj]
    from pysisyphus.constants import BOHR2ANG
    pos = np.array(pos) / BOHR2ANG

    calc = geom.calculator
    calc.plot()

    ax = calc.ax
    ax.plot(*pos[:,0,:2].T)

    plt.show()


if __name__ == "__main__":
    ase_md_playground()
