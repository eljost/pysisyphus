import glob
import os

import numpy as np
import pyparsing as pp


def to_float(s, loc, toks):
    match = toks[0].replace("D", "E")
    return float(match)

def make_float_class(**kwargs):
    return pp.Word(pp.nums + ".-DE+", **kwargs).setParseAction(to_float)


def parse_turbo_gradient(path):
    results = {}
    gradient_fn = glob.glob(os.path.join(path, "gradient"))
    if not gradient_fn:
        raise Exception("gradient file not found!")
    assert(len(gradient_fn) == 1)
    gradient_fn = gradient_fn[0]
    with open(gradient_fn) as handle:
        text = handle.read()

    float_ = make_float_class()
    cycle = pp.Word(pp.nums).setResultsName("cycle")
    scf_energy = float_.setResultsName("scf_energy")
    grad_norm = float_.setResultsName("grad_norm")
    float_line = float_ + float_ + float_
    coord_line = pp.Group(float_line + pp.Word(pp.alphas))
    grad_line = pp.Group(float_line)
    cart_grads = pp.Literal("cartesian gradients")
    energy_type = pp.Or((pp.Literal("SCF energy"),
                        pp.Literal("ex. state energy"),
                        pp.Literal("CC2 energy"),
                        pp.Literal("ADC(2) energy"),
                        pp.Literal("MP2 energy"),
    ))

    parser = (
        pp.Or((pp.Literal("$grad"), pp.Literal("$gradient"))) + pp.Optional(cart_grads) +
        pp.Literal("cycle =") + cycle +
        energy_type + pp.Literal("=") + scf_energy +
        pp.Literal("|dE/dxyz| =") + grad_norm +
        pp.OneOrMore(coord_line).setResultsName("coords") +
        pp.OneOrMore(grad_line).setResultsName("grad") +
        pp.Literal("$end")
    )
    parsed = parser.parseString(text)
    gradient = np.array(parsed["grad"].asList()).flatten()

    results["energy"] = parsed["scf_energy"]
    results["forces"] = -gradient
    return results


def parse_turbo_ccre0_ascii(text):
    float_ = make_float_class()
    float_20 = make_float_class(exact=20)
    int_ = pp.Word(pp.nums).setParseAction(lambda s, loc, toks: int(toks[0]))

    title = pp.Literal("$CCRE0-") + pp.Word(pp.nums + "-")
    method = pp.Word(pp.alphanums + "()") + float_ + float_ + int_ + int_ + int_
    data = pp.Group(pp.OneOrMore(float_20))
    end = pp.Literal("$end")

    parser = (title
              + method
              + data.setResultsName("data")
              + end

    )
    result = parser.parseString(text)
    data = np.array(result["data"].asList())
    return data


def parse_turbo_mos(text):
    float_20 = make_float_class(exact=20)
    int_ = pp.Word(pp.nums)
    comment = pp.Literal("#") + pp.restOfLine

    mo_num = int_
    sym = pp.Word(pp.alphanums)
    eigenvalue = pp.Literal("eigenvalue=") + float_20
    nsaos = pp.Literal("nsaos=") + int_
    mo_coeffs = pp.OneOrMore(float_20)

    mo = pp.Group(
        mo_num + sym + eigenvalue + nsaos +
        mo_coeffs.setResultsName("mo_coeffs")
    )

    parser = (
        pp.Literal("$scfmo") + pp.Literal("scfconv=") + pp.Word(pp.nums)
        + pp.Literal("format(4d20.14) ")
        + pp.ZeroOrMore(comment)
        + pp.OneOrMore(mo).setResultsName("mos")
        + pp.Literal("$end")
    )
    parsed = parser.parseString(text)
    mo_coeffs = np.array(
        [mo.mo_coeffs.asList() for mo in parsed.mos]
    )

    return mo_coeffs


def parse_turbo_exstates(text):
    """Parse excitation energies (first blocks) from an exstates file."""
    float_ = make_float_class()
    exc_energies_line = (pp.Literal("$excitation_energies_")
                         + pp.Word(pp.alphanums + "()").setResultsName("model")
                         + pp.Word("_")
                         + pp.restOfLine
    )
    exc_energy = pp.Suppress(pp.Word(pp.nums)) + float_
    exc_energies_block = pp.Group(exc_energies_line
                                  + pp.Group(pp.OneOrMore(exc_energy)).setResultsName("exc_ens"))

    parser = pp.OneOrMore(exc_energies_block).setResultsName("exc_blocks")
    result = parser.parseString(text)
    exc_energies_by_model = [(b.model, b.exc_ens.asList()) for b in result.exc_blocks]
    return exc_energies_by_model
