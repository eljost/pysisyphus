import jinja2

try:
    from thermoanalysis.QCData import QCData
    from thermoanalysis.thermo import thermochemistry
except ModuleNotFoundError:
    print("Could not import 'thermoanalysis'")

from pysisyphus.constants import AU2KJPERMOL
from pysisyphus.helpers import highlight_text


def get_thermoanalysis(geom, T=298.15):
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

    qcd = QCData(thermo_dict)
    thermo = thermochemistry(qcd, temperature=298.15)

    return thermo


THERMO_TPL = jinja2.Template(
    """
Temperature       : {{ thermo.T }} K
Pressure          : {{ thermo.p }} Pa
Total Mass        : {{ thermo.M }} amu

! Symmetry is currently not supported in pysisyphus, !
! so point group will always be c1 and σ = 1.        !
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


def print_thermoanalysis(thermo):
    """Print thermochemical analysis."""

    def fmt(key, hartree):
        """Output key & energy in Hartree and kJ/mol."""
        kjm = hartree * AU2KJPERMOL
        return f"{key:<18}: {hartree: >16.6f} Eh ({kjm: >18.2f} kJ/mol)"

    # Separator
    sep = "+-----------------------------------------------------------------+"

    rendered = THERMO_TPL.render(thermo=thermo, sep=sep, fmt=fmt)

    print(highlight_text("Thermochemical corrections"))
    print(rendered)
