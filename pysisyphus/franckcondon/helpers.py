from pysisyphus.constants import AU2SEC, JOULE2EV, C as speed_of_light, PI, HBAR


def nu2angfreq_au(wavenum):
    """Angular frequency in atomic units from wavenumber in cm⁻¹."""
    angfreq_si = speed_of_light * wavenum * 1e2 * 2 * PI
    return angfreq_si * AU2SEC


def nu2eV(wavenum):
    """Wavenumber to electronvolt."""
    joule = HBAR * (2 * PI) * speed_of_light * wavenum * 1e2
    return joule * JOULE2EV
