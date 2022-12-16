import pyparsing as pp

from pysisyphus.helpers_pure import file_or_str


@file_or_str(".psf")
def parse_psf(text):
    int_ = pp.common.integer
    real = pp.common.sci_real
    psf_style_comment = pp.Regex(r"\*.*")

    psf = pp.CaselessLiteral("PSF")
    token = pp.ZeroOrMore(
        pp.CaselessLiteral("EXT")
        | pp.CaselessLiteral("CMAP")
        | pp.CaselessLiteral("XPLOR")
        | pp.CaselessLiteral("AUTOG")
    ).set_results_name("token")
    ntitle = int_ + pp.CaselessLiteral("!NTITLE")
    natom = int_.set_results_name("natom") + pp.CaselessLiteral("!NATOM")

    # Atom records
    atom_name = pp.Word(pp.alphanums)
    zero = pp.Literal("0")
    atom_record = pp.Group(
        int_.set_results_name("id")
        + pp.Word(pp.alphas).set_results_name("segment")
        + int_.set_results_name("resid")
        + pp.Word(pp.alphanums).set_results_name("resname")  # VAL, ALA, TIP3, etc.
        + atom_name.set_results_name("atom_name")
        + atom_name.set_results_name("atom_type")
        + real.set_results_name("charge")
        + real.set_results_name("mass")
        + zero
    )

    new_section = int_ + pp.Literal("!")

    def get_section(lhs, rhs=None):
        parser = int_.set_results_name("nterm") + pp.CaselessLiteral(lhs)

        if rhs is not None:
            parser += pp.Literal(":") + pp.CaselessLiteral(rhs)
        # Negative lookahead prevents matching 'nterm' of the next section
        parser += pp.Group(pp.ZeroOrMore(~new_section + int_)).set_results_name("inds")
        parser = pp.Group(parser)
        return parser

    nbond = get_section("!NBOND", "bonds").set_results_name("nbond")
    ntheta = get_section("!NTHETA", "angles").set_results_name("ntheta")
    nphi = get_section("!NPHI", "dihedrals").set_results_name("nphi")
    nimphi = get_section("!NIMPHI", "impropers").set_results_name("nimphi")
    ndon = get_section("!NDON", "donors").set_results_name("ndon")
    nacc = get_section("!NACC", "acceptors").set_results_name("nacc")
    nnb = get_section("!NNB").set_results_name("nnb")
    ncrterm = get_section("!NCRTERM", "cross-terms").set_results_name("ncrterm")

    parser = (
        psf
        + token
        + ntitle
        + natom
        + pp.OneOrMore(atom_record).set_results_name("atoms")
        + nbond  # bonds
        + ntheta  # bends
        + nphi  # dihedrals
        + nimphi  # impropers
        + ndon  # donors
        + nacc  # acceptors
        + nnb  # Non bonding something?
        # 0 0 !NGRP NST2
        # 0 0 !NUMLPH NUMLPH
        # + ncrterm
    )
    parser.ignore(psf_style_comment)

    result = parser.parseString(text)
    as_dict = result.asDict()
    return as_dict
