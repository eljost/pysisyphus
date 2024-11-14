# [1] https://doi.org/10.1063/1.4894472
#     Diabatization based on the dipole and quadrupole: The DQ method
#     Hoyer, Xu, Ma, Gagliardi, Truhlar, 2014
# [2] https://doi.org/10.1063/1.4948728
#     The DQ and DQΦ electronic structure diabatization methods:
#     Validation for general applications
#     Hoyer, Parker, Galgiardi, Truhlar, 2016
# [3] https://doi.org/10.1063/1.3042233
#     Constructing diabatic states from adiabatic states: Extending
#     generalized Mulliken–Hush to multiple charge centers with Boys localization
#     Subotnik, Yeganeh, Cave, Ratner, 2008
# [4]  https://doi.org/10.1007/BF01113521
#      Efficient use of Jacobi rotations for orbital optimization and localization
#      Raffenetti, Ruedenberg, Janssen, Schaefer, 1993
# [5]  https://doi.org/10.1103/RevModPhys.35.457
#      Localized Atomic and Molecular Orbitals
#      Edmiston, Ruedenberg, 1963
# [6] https://doi.org/10.1063/1.3148777
#      The initial and final states of electron and energy transfer processes:
#      Diabatization as motivated by system-solvent interactions
#      Subotnik, Cave, Steele, Shenvi, 2009

import dataclasses
from typing import Optional

from jinja2 import Template
import numpy as np

from pysisyphus.constants import EV2NU
from pysisyphus.helpers_pure import to_subscript_num
from pysisyphus.wavefunction.localization import JacobiSweepResult


DiaResultTemplate = Template(
    """
########################
# DIABATIZATION REPORT #
########################

Kind: {{ kind }}

All energies are given in {{ unit }}.

Every column of the rotation matrix U describes the composition of
a diabatic state in terms of (possibly various) adiabatic states.

Adiabatic energy matrix V
-------------------------
{{ adia_mat }}

Rotation matrix U
-----------------
{{ U }}
det(U)={{ "%.4f"|format(det) }}

Diabatic energy matrix D = UᵀVU
-------------------------------
{{ dia_mat }}

Diabatic states Ξᵢ sorted by energy
-----------------------------------
{%- for ind, dia_state, dia_en in dia_states_sorted %}
{{ ind }}: {{ dia_state }}, {{ "%.4f"|format(dia_en) }} {{ unit }}
{%- endfor %}

Composition of diabatic states Ξᵢ
---------------------------------
{%- for dia_comp in dia_compositions %}
{{ dia_comp }}
{%- endfor %}

Weights U²
----------
{{ U2 }}

Unique absolute diabatic couplings
----------------------------------
{%- for key, coupling in couplings %}
|{{ key }}| = {{ "%.5f"|format(coupling) }} {{ unit }}
{%- if unit == "eV" %}, ({{ "%8.2f"|format(coupling * EV2NU) }} cm⁻¹){% endif %}
{%- endfor %}
"""
)


@dataclasses.dataclass
class DiabatizationResult:
    kind: str
    U: np.ndarray
    adia_ens: np.ndarray
    is_converged: bool
    cur_cycle: int
    P: float
    # Adiabatic properties
    #
    # Dipole moments
    dip_moms: Optional[np.ndarray] = None
    # Quadrupole moments
    quad_moms: Optional[np.ndarray] = None
    # Electrostatic potentials
    epots: Optional[np.ndarray] = None
    # 4d Coulomb tensor  of shape (nstates, nstates, nstates, nstates)
    R_tensor: Optional[np.ndarray] = None  # Coulomb tensor
    # 3d density fitting tensor of shape (naux, nstates, nstates)
    L_tensor: Optional[np.ndarray] = None
    # TODO: add adiabatic labels and use them in the report

    def __post_init__(self):
        self.adia_mat = np.diag(self.adia_ens)
        self.dia_mat = self.U.T @ self.adia_mat @ self.U
        self.dia_ens = np.diag(self.dia_mat)

    def sort(self):
        inds = self.dia_ens.argsort()
        Usort = self.U[:, inds]
        kwargs = dataclasses.asdict(self)
        kwargs["U"] = Usort
        return DiabatizationResult(**kwargs)

    def savez(self, fn, **add_kwargs):
        kwargs = dataclasses.asdict(self)
        kwargs.update(add_kwargs)
        np.savez(fn, **kwargs)

    @property
    def nstates(self):
        return len(self.adia_ens)

    @property
    def couplings(self):
        nstates = self.nstates
        abs_dia_mat = np.abs(self.dia_mat)
        couplings = dict()
        for from_ in range(nstates):
            for to_ in range(from_ + 1, nstates):
                key = (from_, to_)
                couplings[key] = couplings[key[::-1]] = abs_dia_mat[from_, to_]
        return couplings

    def render_report(self, adia_labels=None, unit="eV"):
        U = self.U
        nstates = self.nstates
        U2 = U**2
        det = np.linalg.det(U)

        subscripts = [to_subscript_num(i) for i in range(nstates)]
        adia_wfs = [f"Φ{subscripts[i]}" for i in range(nstates)]
        dia_wfs = [f"Ξ{subscripts[i]}" for i in range(nstates)]

        # Diabatic states, sorted by energy
        dia_ens = self.dia_ens
        sorted_inds = np.argsort(dia_ens)
        dia_states_sorted = list()
        for i, sort_ind in enumerate(sorted_inds):
            dia_states_sorted.append((i, dia_wfs[sort_ind], dia_ens[sort_ind]))

        # Composition of diabatic states
        dia_compositions = list()
        if adia_labels:
            adia_labels = [f"({al})" for al in adia_labels]
        else:
            adia_labels = ["" for _ in range(nstates)]
        for i in range(nstates):
            coeffs = U[:, i]
            signs = ["+" if c > 0 else "-" for c in coeffs]
            sum_ = " ".join(
                [
                    f"{sign} {c:>6.4f}·{awf}{al}"
                    for sign, c, awf, al in zip(
                        signs, np.abs(coeffs), adia_wfs, adia_labels
                    )
                ]
            )
            dia_compositions.append(f"{dia_wfs[i]} = {sum_}")

        # Couplings
        couplings = self.couplings
        flat_couplings = list()
        for from_ in range(nstates):
            for to_ in range(from_ + 1, nstates):
                key = f"D{subscripts[from_]}{subscripts[to_]}"
                flat_couplings.append((key, couplings[(from_, to_)]))

        def fmta(arr, fmt):
            return np.array2string(
                arr,
                formatter={"float": lambda f: f"{f:{fmt}}"},
                max_line_width=200,
                # suppress_small=True,
            )

        def fmtu(arr):
            """For U/U2"""
            return fmta(arr, "8.4f")

        def fmten(arr):
            """For adia_mat/dia_mat"""
            return fmta(arr, "10.6f")

        rendered = DiaResultTemplate.render(
            kind=self.kind,
            unit=unit,
            adia_mat=fmten(self.adia_mat),
            U=fmtu(U),
            det=det,
            dia_mat=fmten(self.dia_mat),
            dia_states_sorted=dia_states_sorted,
            dia_compositions=dia_compositions,
            U2=fmtu(U2),
            couplings=flat_couplings,
            EV2NU=EV2NU,
        )
        return rendered


def dia_result_from_jac_result(
    kind: str,
    adia_ens: np.ndarray,
    jac_res: JacobiSweepResult,
    sort: bool = True,
    **property_tensors,
) -> DiabatizationResult:
    """DiabatizationResult construction wrapper.

    Parameter
    ---------
    kind
        String label containing the name of the diabatization algorithm,
        e.g., 'er' or 'boys'.
    adia_ens
        1d array of shape (nstates, ) holding the electronic energies of
        the adiabatic states.
    jac_res
        JacobiSweepResult from the diabatiatzion function.
    sort
        Boolean flag that indicates whether the diabatic states will be
        sorted by their energies. If true, the columns of U will be reorderd.
    property_tensors
        Various property tensors of varying shape containing e.g., the Coulomb
        tensor or multipole moments.

    Returns
    -------
    dia_result
        DiabatizationResult.
    """
    U = jac_res.C
    dia_result = DiabatizationResult(
        kind=kind,
        U=U,
        adia_ens=adia_ens,
        is_converged=jac_res.is_converged,
        cur_cycle=jac_res.cur_cycle,
        P=jac_res.P,
        **property_tensors,
    )
    if sort:
        dia_result = dia_result.sort()
    return dia_result
