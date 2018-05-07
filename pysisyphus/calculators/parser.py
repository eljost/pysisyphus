#!/usr/bin/env python3

import glob
import os

import numpy as np
import pyparsing as pp


def to_float(s, loc, toks):
    match = toks[0].replace("D", "E")
    return float(match)

def make_float_class(**kwargs):
    return pp.Word(pp.nums + ".-D+", **kwargs).setParseAction(to_float)


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
                        pp.Literal("MP2 energy"),
    ))

    parser = (
        pp.Literal("$grad") + pp.Optional(cart_grads) +
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


if __name__ == "__main__":
    from pathlib import Path
    fn = "/tmp/calculator_0_000_ab7q7o9y"
    path = Path(fn)
    parse_turbo_gradient(path)
