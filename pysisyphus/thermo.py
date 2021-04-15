import h5py
import jinja2

try:
    from thermoanalysis.QCData import QCData
    from thermoanalysis.thermo import thermochemistry
    can_thermoanalysis = True
except ModuleNotFoundError:
    print("Could not import 'thermoanalysis'")
    can_thermoanalysis = False

from pysisyphus.constants import AU2KJPERMOL
from pysisyphus.helpers_pure import highlight_text


def get_thermoanalysis_from_hess_h5(h5_fn, T=298.15, point_group="c1"):
    with h5py.File(h5_fn, "r") as handle:
        masses = handle["masses"][:]
        vibfreqs = handle["vibfreqs"][:]
        coords3d = handle["coords3d"][:]
        energy = handle.attrs["energy"]
        mult = handle.attrs["mult"]

    thermo_dict = {
        "masses": masses,
        "vibfreqs": vibfreqs,
        "coords3d": coords3d,
        "energy": energy,
        "mult": mult,
    }

    qcd = QCData(thermo_dict, point_group=point_group)
    thermo = thermochemistry(qcd, temperature=T)
    return thermo


def get_thermoanalysis(geom, T=298.15, point_group="c1"):
    hessian = geom.cart_hessian
    energy = geom.energy
    vibfreqs, *_ = geom.get_frequencies(hessian)
    try:
        mult = geom.calculator.mult
    except AttributeError:
        mult = 1
        print(f"Multiplicity could not be determined! Using 2S+1 = {mult}.")

    thermo_dict = {
        "masses": geom.masses,
        "vibfreqs": vibfreqs,
        "coords3d": geom.coords3d,
        "energy": energy,
        "mult": mult,
    }

    qcd = QCData(thermo_dict, point_group=point_group)
    thermo = thermochemistry(qcd, temperature=T)

    return thermo


THERMO_TPL = jinja2.Template(
    """
{% if geom -%}
Geometry          : {{ geom }}, {{ geom.atoms|length }} atoms
{%- endif %}
Temperature       : {{ "%0.2f" % thermo.T }} K
Pressure          : {{ thermo.p }} Pa
Total Mass        : {{ "%0.4f" % thermo.M }} amu

! Symmetry is currently not supported in pysisyphus. !
! If not specified c1 and σ = 1 are assumed.         !
Point Group       : {{ thermo.point_group }}
Symmetry Number σ : {{ thermo.sym_num }}  

+
| Inner energy U = U_el + U_vib + U_rot + U_trans
{{ sep }}
{{ fmt("U_el", thermo.U_el) }}
{{ fmt("ZPE", thermo.ZPE) }}
{{ fmt("U_vib (incl. ZPE)", thermo.U_vib) }}
{{ fmt("U_rot", thermo.U_rot) }}
{{ fmt("U_trans", thermo.U_trans) }}
{{ sep }}
{{ fmt("U", thermo.U_tot) }}

+
| Enthalpy H = U + kB*T
{{ sep }}
{{ fmt("U", thermo.U_tot) }}
{{ fmt("kB*T", thermo.kBT) }}
{{ sep }}
{{ fmt("H", thermo.H) }}

+
| Entropy correction T*S = T*(S_el + S_vib + S_rot + S_trans)
{{ sep }}
{{ fmt("T*S_el", thermo.TS_el) }}
{{ fmt("T*S_vib", thermo.TS_vib) }}
{{ fmt("T*S_rot", thermo.TS_rot) }}
{{ fmt("T*S_trans", thermo.TS_trans) }}
{{ sep }}
{{ fmt("T*S", thermo.TS_tot) }}

+
| Gibbs free energy G = H - T*S
{{ sep }}
{{ fmt("H", thermo.H) }}
{{ fmt("T*S", thermo.TS_tot) }}
{{ sep }}
{{ fmt("G", thermo.G) }}
{{ fmt("dG", thermo.dG) }}
""".strip()
)


def print_thermoanalysis(thermo, geom=None, level=0, title=None):
    """Print thermochemical analysis."""

    def fmt(key, hartree):
        """Output key & energy in Hartree and kJ/mol."""
        kjm = hartree * AU2KJPERMOL
        return f"{key:<18}: {hartree: >16.6f} Eh ({kjm: >18.2f} kJ/mol)"

    # Separator
    sep = "+-----------------------------------------------------------------+"

    rendered = THERMO_TPL.render(geom=geom, thermo=thermo, sep=sep, fmt=fmt)

    if title is None:
        title = ""
    else:
        title = f", {title}"
    print(highlight_text(f"Thermochemical corrections{title}", level=level))
    print()
    print(rendered)
