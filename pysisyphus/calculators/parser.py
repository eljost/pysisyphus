import glob
import os

import numpy as np
import pyparsing as pp
import re

from pysisyphus.helpers_pure import file_or_str


def to_float(s, loc, toks):
    match = toks[0].replace("D", "E")
    return float(match)


def make_float_class(**kwargs):
    return pp.Word(pp.nums + ".-DE+", **kwargs).setParseAction(to_float)


@file_or_str("gradient", exact=True)
def parse_turbo_gradient_str(text):
    float_ = make_float_class()

    title = pp.Or((pp.Literal("$grad"), pp.Literal("$gradient")))
    cycle = pp.Word(pp.nums).setResultsName("cycle")
    scf_energy = float_.setResultsName("scf_energy")
    grad_norm = float_.setResultsName("grad_norm")
    float_line = float_ + float_ + float_
    coord_line = pp.Group(float_line + pp.Word(pp.alphas))
    grad_line = pp.Group(float_line)
    cart_grads = pp.Literal("cartesian gradients")
    energy_type = pp.Or(
        (
            pp.Literal("SCF energy"),
            pp.Literal("ex. state energy"),
            pp.Literal("CC2 energy"),
            pp.Literal("ADC(2) energy"),
            pp.Literal("MP2 energy"),
        )
    )
    cycle_grads = pp.Group(
        pp.Literal("cycle =")
        + cycle
        + energy_type
        + pp.Literal("=")
        + scf_energy
        + pp.Literal("|dE/dxyz| =")
        + grad_norm
        + pp.OneOrMore(coord_line).setResultsName("coords")
        + pp.OneOrMore(grad_line).setResultsName("grad")
    )

    parser = (
        title
        + pp.Optional(cart_grads)
        + pp.OneOrMore(cycle_grads).set_results_name("cycle_grads")
        + pp.Literal("$end")
    )
    parsed = parser.parseString(text)
    return parsed


def parse_turbo_gradient(path):
    gradient_fn = glob.glob(os.path.join(path, "gradient"))
    if not gradient_fn:
        raise Exception("gradient file not found!")
    assert len(gradient_fn) == 1
    gradient_fn = gradient_fn[0]
    parsed = parse_turbo_gradient_str(gradient_fn)
    cycle_grads = parsed["cycle_grads"]

    # There may be multiple gradients of different optimization cycles
    # in the file. Only use latest/last gradient.
    latest_grad = cycle_grads[-1]
    gradient = np.array(latest_grad["grad"].asList()).flatten()

    results = {
        "energy": latest_grad["scf_energy"],
        "forces": -gradient,
    }
    return results


def parse_turbo_ccre0_ascii(text):
    float_ = make_float_class()
    float_20 = make_float_class(exact=20)
    int_ = pp.Word(pp.nums).setParseAction(lambda s, loc, toks: int(toks[0]))

    title = pp.Literal("$CCRE0-") + pp.Word(pp.nums + "-")
    method = pp.Word(pp.alphanums + "()") + float_ + float_ + int_ + int_ + int_
    data = pp.Group(pp.OneOrMore(float_20))
    end = pp.Literal("$end")

    parser = title + method + data.setResultsName("data") + end
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
        mo_num + sym + eigenvalue + nsaos + mo_coeffs.setResultsName("mo_coeffs")
    )

    parser = (
        (pp.Literal("$scfmo") | pp.Literal("$uhfmo_alpha") | pp.Literal("$uhfmo_beta"))
        + pp.Literal("scfconv=")
        + pp.Word(pp.nums)
        + pp.Literal("format(4d20.14)")
        + pp.ZeroOrMore(comment)
        + pp.OneOrMore(mo).setResultsName("mos")
        + pp.Literal("$end")
    )
    parsed = parser.parseString(text)
    # MOs are in columns
    mo_coeffs = np.array([mo.mo_coeffs.asList() for mo in parsed.mos]).T
    return mo_coeffs


def parse_turbo_exstates(text):
    """Parse excitation energies (first blocks) from an exstates file."""
    float_ = make_float_class()
    exc_energies_line = (
        pp.Literal("$excitation_energies_")
        + pp.Word(pp.alphanums + "()").setResultsName("model")
        + pp.Word("_")
        + pp.restOfLine
    )
    exc_energy = pp.Suppress(pp.Word(pp.nums)) + float_
    exc_energies_block = pp.Group(
        exc_energies_line + pp.Group(pp.OneOrMore(exc_energy)).setResultsName("exc_ens")
    )

    parser = pp.OneOrMore(exc_energies_block).setResultsName("exc_blocks")
    result = parser.parseString(text)
    exc_energies_by_model = [(b.model, b.exc_ens.asList()) for b in result.exc_blocks]
    return exc_energies_by_model


def parse_turbo_exstates_re(text):
    results = dict()
    exc_ens_str = "excitation_energies_"
    model_re = re.compile(exc_ens_str + r"(\w+?)_")
    blocks = text.strip().split("$")
    for block in blocks:
        if not block.startswith(exc_ens_str):
            continue
        model_line, *exc_lines = block.strip().split("\n")
        mobj = model_re.search(model_line)
        model = mobj.group(1)
        exc_lines = [line.strip().split() for line in exc_lines]
        assert all([len(line) == 2 for line in exc_lines])
        exc_ens = np.array([exc_en for _, exc_en in exc_lines], dtype=float)
        results[model] = exc_ens
    return results
