#!/usr/bin/env python

import argparse
import json
import sys

import numpy as np
import pyparsing as pp

from pysisyphus.elem_data import ATOMIC_NUMBERS
from pysisyphus.wavefunction.helpers import get_l


def parse_libmol(text):
    integer = pp.common.integer
    sci_real = pp.common.sci_real | integer
    word = pp.Word(pp.printables)
    libmol_style_comment = pp.Regex(r"!.*").set_name("libmol style comment")

    atom = pp.Word(pp.alphas)("atom")
    angmom = pp.one_of("s p d f g h i")("angmom")
    colon = pp.Literal(":")
    basis_label = pp.OneOrMore(~colon + word)("basis_label")
    coeff_num = integer("coeff_num")
    shell_num = integer("shell_num")
    shell_len = ~pp.line_end + pp.Word(pp.nums + ".")
    shell_lens = pp.OneOrMore(shell_len)("shell_lens")
    comment_line = pp.Optional(pp.line_start) + pp.SkipTo(pp.line_end)
    exps_coeffs = pp.OneOrMore(sci_real)("exps_coeffs")
    shell = pp.Group(
        atom
        + angmom
        + basis_label
        + colon
        + coeff_num
        + shell_num
        + shell_lens
        + pp.line_end
        #
        + comment_line
        #
        + exps_coeffs
    )
    parser = pp.OneOrMore(shell)("shells")
    parser.ignore(libmol_style_comment)

    result = parser.parseString(text, parse_all=True)
    as_dict = result.asDict()
    return as_dict


def convert(as_dict):
    shells = as_dict["shells"]
    elements = {}

    for shell in shells:
        atom = shell["atom"]
        atomic_num = ATOMIC_NUMBERS[atom.lower()]
        l = get_l(shell["angmom"])
        coeff_num = shell["coeff_num"]
        shell_num = shell["shell_num"]
        exps_coeffs = np.array(shell["exps_coeffs"])
        try:
            assert exps_coeffs.size == coeff_num * (shell_num + 1)
        except AssertionError as err:
            raise err
        exps, *coeffs = exps_coeffs.reshape(shell_num + 1, coeff_num)

        elements.setdefault(str(atomic_num), dict()).setdefault(
            "electron_shells", list()
        )
        el_shells = elements[str(atomic_num)]["electron_shells"]
        el_shell = {
            "function_type": "gto",
            "region": "",
            "angular_momentum": [
                l,
            ],
            "exponents": list(exps),
            "coefficients": [cs.tolist() for cs in coeffs],
        }
        el_shells.append(el_shell)
    final = {
        "elements": elements,
        "revision_description": (
            "Parsed from 'minao.libmol' as found in ir/wmee 2020-02-28 "
            "Gerald Knizias homepage. "
            "61c0e27f96a98bebf26a98343204901acf0e88dc  minao.libmol"
        ),
        "revision_date": "2022-23-05",
        "description": "minao",
        "name": "minao",
        "version": 1,
    }
    return final


def parse_args(args):
    parser = argparse.ArgumentParser()
    parser.add_argument("--in-fn", default="minao.libmol")
    parser.add_argument("--out-fn", default="minao.json")
    return parser.parse_args(args)


def run():
    args = parse_args(sys.argv[1:])

    in_fn = args.in_fn
    out_fn = args.out_fn

    with open(in_fn) as handle:
        text = handle.read()
    print(f"Read '{in_fn}'.")
    as_dict = parse_libmol(text)
    final = convert(as_dict)

    with open(out_fn, "w") as handle:
        json.dump(final, handle, indent=4, sort_keys=True)
    print(f"Dumped json to '{out_fn}'.")


if __name__ == "__main__":
    run()
