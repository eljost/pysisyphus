import h5py
import jinja2

try:
    from thermoanalysis.QCData import QCData
    from thermoanalysis.thermo import thermochemistry

    can_thermoanalysis = True
except ModuleNotFoundError:
    can_thermoanalysis = False

from pysisyphus.config import p_DEFAULT, T_DEFAULT
from pysisyphus.constants import AU2KJPERMOL, AU2KCALPERMOL
from pysisyphus.helpers_pure import highlight_text
from pysisyphus.Geometry import Geometry


def get_thermoanalysis_from_hess_h5(
    h5_fn, T=T_DEFAULT, p=p_DEFAULT, point_group="c1", return_geom=False
):
    with h5py.File(h5_fn, "r") as handle:
        masses = handle["masses"][:]
        vibfreqs = handle["vibfreqs"][:]
        coords3d = handle["coords3d"][:]
        energy = handle.attrs["energy"]
        mult = handle.attrs["mult"]
        atoms = handle.attrs["atoms"]

    thermo_dict = {
        "masses": masses,
        "wavenumbers": vibfreqs,
        "coords3d": coords3d,
        "scf_energy": energy,
        "mult": mult,
    }

    qcd = QCData(thermo_dict, point_group=point_group)
    thermo = thermochemistry(qcd, temperature=T, pressure=p)
    if return_geom:
        geom = Geometry(atoms=atoms, coords=coords3d)
        return thermo, geom
    else:
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
! If not given otherwise, c1 and σ = 1 are assumed.  !
Point Group       : {{ thermo.point_group }}
Symmetry Number σ : {{ thermo.sym_num }}  
Linear            : {{ thermo.linear }}

+
| Normal Mode Wavenumbers
{{ sep }}
{% for nu in used_nus -%}
 {{ "\t%04d" % loop.index }}: {{ nu }}
{% endfor -%}
{{ sep }}
{% if is_ts %}This should be a TS.{% endif %}
Expected {{ expected }} normal modes, got {{ used_nus|length}}.

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


def print_thermoanalysis(
    thermo, geom=None, is_ts=False, unit="joule", level=0, title=None
):
    """Print thermochemical analysis."""

    units = {
        "calorie": ("kcal mol⁻¹", AU2KCALPERMOL),
        "joule": ("kJ mol⁻¹", AU2KJPERMOL),
    }
    unit_key, unit_conv = units[unit]

    def fmt(key, hartree):
        """Output key & energy in Hartree and the chosen unit."""
        kjm = hartree * unit_conv
        return f"{key:<18}: {hartree: >16.8f} Eh ({kjm: >18.2f} {unit_key})"

    # Separator
    sep = "+-----------------------------------------------------------------+"

    sub_modes = 5 if thermo.linear else 6
    # Expect one real mode less if it is a TS
    sub_modes += 1 if is_ts else 0
    expected = 3 * thermo.atom_num - sub_modes

    def fmt_nus(nus):
        return [f"{nu: >8.2f} cm⁻¹" + (", excluded" if nu < 0.0 else "") for nu in nus]

    rendered = THERMO_TPL.render(
        geom=geom,
        thermo=thermo,
        org_nus=thermo.org_wavenumbers,
        used_nus=fmt_nus(thermo.wavenumbers),
        expected=expected,
        is_ts=is_ts,
        sep=sep,
        fmt=fmt,
    )

    if title is None:
        title = ""
    else:
        title = f", {title}"
    print(highlight_text(f"Thermochemical corrections{title}", level=level))
    print()
    print(rendered)
