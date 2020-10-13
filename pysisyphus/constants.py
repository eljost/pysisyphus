BOHR2M = 5.291_772_109e-11
BOHR2ANG = BOHR2M * 1e10
ANG2BOHR = 1.889726125457828
AU2J = 4.359_744_722_207_1e-18
AU2KJPERMOL = 2625.499
AU2EV = 27.211386
# eV/Å -> Hartree/Bohr
EVANG2AUBOHR = 1/AU2EV/ANG2BOHR
# fs -> Bohr * sqrt(amu/Hartree)
FS2AU = 0.9682885864793366
# Boltzman constant
KB = 1.38064852E-23  # (m² kg s⁻² K⁻¹) or just (J / K)
KBAU = KB / AU2J
# Atomic mass unit to kg
AMU2KG = 1.660_539_066_60e-27
AU2SEC = 2.4188843265857e-17
# Hartree to kcal mol⁻¹
AU2KCALMOL = 627.509474

# MD related

# Force/amu to acceleration
#   Hartree/(Bohr*amu) -> Bohr / fs²
FORCE2ACC = AU2J / (AMU2KG * BOHR2M**2 * 1e30)
# Velocity*amu -> Energy
#   Bohr²/fs²*amu -> Hartree
# VELO2E = AMU2KG * BOHR2M**2 / (1e-30 * AU2J)
VELO2E = 1 / FORCE2ACC
# Velocity from Bohr/fs to Bohr/t_atomic_unit
BOHRPERFS2AU = AU2SEC*1e15
