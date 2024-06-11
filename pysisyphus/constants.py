import scipy.constants as spc

PI = spc.pi

# Speed of light in m/s
C = spc.c
# Speed of light in atomic units
CAU = spc.value("inverse fine-structure constant")
# ħ = h / 2π in J/s
HBAR = spc.hbar
# h
PLANCK = spc.Planck
PLANCKAU = 2 * spc.pi  # With ħ_au = 1 and ħ_au=h_au/2π -> h_au = 2π
# Conversion factor for momentum in SI to atomic units
P_AU = spc.value("atomic unit of momentum")

# Bohr radius in m
BOHR2M = spc.value("Bohr radius")
# Bohr -> Å conversion factor
BOHR2ANG = BOHR2M * 1e10
# Å -> Bohr conversion factor
ANG2BOHR = 1 / BOHR2ANG
# Hartree to J
AU2J = spc.value("Hartree energy")
# Hartree to kJ / mol
AU2KJPERMOL = AU2J / 1000 * spc.Avogadro
# Hartree to kcal mol⁻¹
AU2KCALPERMOL = AU2KJPERMOL / spc.calorie
# Hartree to eV
AU2EV = spc.value("Hartree energy in eV")
# Hartree to wavelength in nm
# ! This is not a multiplicative constant but it must be divided
# ! by the energy in Hartree to yield an excitation energy in nm
_1OVER_AU2NM = PLANCK * C * 1e9 / AU2J
# !
# Joule to eV
JOULE2EV = AU2EV / AU2J
# Wavenumber to Hartree
NU2AU = spc.h * C * 1e2 / AU2J
# Hartree to Wavenumber
AU2NU = 1 / NU2AU
# eV/Å -> Hartree/Bohr
EVANG2AUBOHR = 1 / AU2EV / ANG2BOHR
# Boltzman constant, (m² kg s⁻² K⁻¹) or just (J / K)
KB = spc.Boltzmann
KBAU = KB / AU2J
# Ideal gas constant
R = spc.gas_constant
# Atomic unit of time to second
AU2SEC = spc.value("atomic unit of time")
AU2FS = AU2SEC * 1e15
# Atomic mass unit to kg
AMU2KG = spc.value("atomic mass constant")
M_E = spc.electron_mass
AMU2AU = AMU2KG / M_E
# Avogadro constant
NA = spc.Avogadro
# Force constants from au to mdyne/Å (force / Ångstrom)
# mdyne/Å = 10⁻³ * 10⁻⁵ N/Å
#         = 10² kg/s²
#
# 1 au = 1556.89 kg/s²
# 1 au = 15.5689 * 10² kg/s²
AU2MDYNEPERANG = M_E / AU2SEC**2 / 100

##############
# MD related #
##############

# Force/amu to acceleration
#   Hartree/(Bohr*amu) -> Bohr / fs²
FORCE2ACC = AU2J / (AMU2KG * BOHR2M**2 * 1e30)
# Velocity*amu -> Energy
#   Bohr²/fs²*amu -> Hartree
# VELO2E = AMU2KG * BOHR2M**2 / (1e-30 * AU2J)
VELO2E = 1 / FORCE2ACC
# Velocity from Bohr/fs to Bohr/t_atomic_unit
BOHRPERFS2AU = AU2SEC * 1e15
KCALPERMOLPERANG2 = 1 / AU2KCALPERMOL / ANG2BOHR**2
