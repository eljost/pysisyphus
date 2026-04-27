"""Test the optional Calculator.get_forces_batch hook for ChainOfStates.

Verifies that when a calculator subclass overrides get_forces_batch, the COS
dispatch uses it once per cycle instead of calling get_forces per image.
The default sequential fallback path is also covered to confirm existing
calculators are unaffected.
"""

import numpy as np

from pysisyphus.calculators.AnaPot import AnaPot
from pysisyphus.cos.NEB import NEB
from pysisyphus.interpolate.helpers import interpolate


class CountingBatchAnaPot(AnaPot):
    """AnaPot with an overridden get_forces_batch hook that records call count."""

    def __init__(self):
        super().__init__()
        self.batch_calls = 0
        self.scalar_calls = 0

    def get_forces(self, atoms, coords, **kw):
        self.scalar_calls += 1
        return super().get_forces(atoms, coords, **kw)

    def get_forces_batch(self, atoms_list, coords_list, **kw):
        self.batch_calls += 1
        # Defer to the per-image implementation but increment scalar_calls so
        # we can assert the batch path is exclusive (it bypasses calc_energy_and_forces
        # and therefore should not increment scalar_calls when used).
        return [super(AnaPot, self).get_forces(a, c, **kw)
                for a, c in zip(atoms_list, coords_list)]


def _build_neb(calculator, n_between=3):
    initial = AnaPot.get_geom((-1.05274, 1.02776, 0))
    final = AnaPot.get_geom((1.94101, 3.85427, 0))
    images = interpolate(initial, final, between=n_between, align=False)
    for img in images:
        img.set_calculator(calculator)
    return NEB(images)


def test_cos_batch_hook_invoked():
    """When the calculator overrides get_forces_batch, COS dispatches once per cycle."""
    calc = CountingBatchAnaPot()
    cos = _build_neb(calc, n_between=3)
    cos.calculate_image_forces()
    assert calc.batch_calls == 1, (
        f"batched path not taken; batch_calls={calc.batch_calls}, "
        f"scalar_calls={calc.scalar_calls}"
    )
    assert calc.scalar_calls == 0, (
        f"per-image fallback ran alongside batch; scalar_calls={calc.scalar_calls}"
    )
    # Sanity: each image's energy and forces are populated.
    for img in cos.images:
        assert img._energy is not None
        assert img._forces is not None


def test_cos_default_calculator_unchanged():
    """A calculator that doesn't override get_forces_batch uses the existing serial path."""
    calc = AnaPot()
    cos = _build_neb(calc, n_between=3)
    cos.calculate_image_forces()
    energies = np.array([img.energy for img in cos.images])
    assert energies.shape == (5,)
    assert np.all(np.isfinite(energies))
