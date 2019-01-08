#!/usr/bin/env python3

import textwrap

class TablePrinter:

    def __init__(self, header, col_fmts, width=12):
        self.header = header
        self.col_fmts = col_fmts
        self.width = width
        w = str(self.width)

        self.fmts = {
                "int": "{:>" + w + "d}",
                "float": "{:>" + w + ".6f}",
                "str": "{:>" + w + "s}",
        }

        self.header_str = " ".join([h.rjust(self.width)
                                    for h in self.header])
        self.conv_str = " ".join([self.fmts[fmt] for fmt in self.col_fmts])
        h0_len = len(self.header[0])
        self.offset = self.width - h0_len
        self.prefix = " " * self.offset
        self.sep = self.prefix + "-" * (len(self.header_str) - self.width + h0_len)

    def print_header(self):
        print(self.header_str)
        print(self.sep)

    def print_row(self, args):
        print(self.conv_str.format(*args))

    def print(self, text):
        print(textwrap.indent(text, self.prefix))


if __name__ == "__main__":
    h = "int int float int float str".split()
    args = [1, 1, 1.0, 3, 1.1453125, "halloo"]
    T = TablePrinter(header=h, col_fmts=h)
    T.print_header()
    T.print_row(args)

    text = "aragawasf garbel\nyololo anaonoasdlo\noi!"
    T.print(text)
    T.print_row(args)
