#!/usr/bin/env python3

import glob
import os

import numpy as np
import pyparsing as pp


def parse_turbo_gradient(path):
    results = {}
    gradient_fn = glob.glob(os.path.join(path, "gradient"))
    if not gradient_fn:
        raise Exception("gradient file not found!")
    assert(len(gradient_fn) == 1)
    gradient_fn = gradient_fn[0]
    with open(gradient_fn) as handle:
        text = handle.read()

    def to_float(s, loc, toks):
        match = toks[0].replace("D", "E")
        return float(match)

    float_ = pp.Word(pp.nums + ".-D+").setParseAction(to_float)
    cycle = pp.Word(pp.nums).setResultsName("cycle")
    scf_energy = float_.setResultsName("scf_energy")
    grad_norm = float_.setResultsName("grad_norm")
    float_line = float_ + float_ + float_
    coord_line = pp.Group(float_line + pp.Word(pp.alphas))
    grad_line = pp.Group(float_line)
    cart_grads = pp.Literal("cartesian gradients")
    energy_type = pp.Or((pp.Literal("SCF energy"),
                        pp.Literal("ex. state energy")))

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


if __name__ == "__main__":
    from pathlib import Path
    fn = "/tmp/calculator_0_000_ab7q7o9y"
    path = Path(fn)
    parse_turbo_gradient(path)
