from pysisyphus.io.psf import parse_psf


def test_parse_3nir(this_dir):
    fn = this_dir / "step1_pdbreader.psf"
    d = parse_psf(fn)
    assert d["natom"] == 618

    def check_section(key, per_term):
        section = d[key]
        assert section["nterm"] == len(section["inds"]) // per_term

    check_section("nbond", 2)
    check_section("ntheta", 3)
    check_section("nphi", 4)
