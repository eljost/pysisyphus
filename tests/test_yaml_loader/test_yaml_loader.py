import pytest
import yaml

from pysisyphus.yaml_mods import get_loader
from pysisyphus.constants import AU2KJPERMOL, AU2KCALPERMOL, AU2EV, ANG2BOHR, BOHR2M


@pytest.mark.parametrize(
    "unit, value",
    (
        # Energies
        ("Eh", 1.0),
        ("kJpermol", AU2KJPERMOL),
        ("kcalpermol", AU2KCALPERMOL),
        ("eV", AU2EV),
        # Times
        ("fs", 1.0),
        ("ps", 1e3),
        ("ns", 1e6),
        # Lengths
        ("a0", 1.0),
        ("nm", BOHR2M * 1e-9),
        ("angstrom", ANG2BOHR),
    ),
)
def test_units(unit, value):
    yaml_str = f"value: !{unit} {value}"
    loader = yaml.SafeLoader
    loader = get_loader()
    run_dict = yaml.load(yaml_str, Loader=loader)
    assert run_dict["value"] == pytest.approx(1.0)
