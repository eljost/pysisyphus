#!/usr/bin/env python3

import random

import numpy as np

class TableFormatter:

    def __init__(self, header, fmts, min_width=7, space=3):
        self.min_width = min_width + (space-1)
        self.space = space

        # Get lengths of header strings
        widths = np.array([len(h) for h in header])
        # Expand entries smaller than min_widht to min_width
        smaller_indices = widths < min_width
        widths[smaller_indices] = min_width
        self.widths = widths

        # Construct header
        #header_fmts = ["{:" + "{}".format(width) + "s}"
        #               for width in self.widths]
        header_fmts = self.min_width_fmts()
        self._header = self.join_format(header_fmts, header)
        self._header += "\n" + (self.space*" ").join(
                ["-"*width for width in self.widths]
        )

        # Modify fmts to consinder min_widths 
        self.fmts = self.min_width_fmts(fmts)

    def min_width_fmts(self, raw_fmts=None):
        if not raw_fmts:
            raw_fmts = ["s" for _ in self.widths]
        return ["{:>" + "{}".format(width) + fmt + "}"
                for width, fmt in zip(self.widths, raw_fmts)]

    def join_format(self, fmts, lst):
        """Format a given iterable lst with formats given in the iterable
        lst and return the joined items of the formatted list."""
        return (self.space*" ").join(
            [fmt.format(item) for fmt, item in zip(fmts, lst)]
        )

    @property
    def header(self):
        return self._header

    def line(self, *args):
        #formatted = " "["
        return self.join_format(self.fmts, args)

def run():
    header = "# |dx| |tangent|".split()
    fmts = ["d", ".2E", ".3E"]
    min_width = 10
    tp = TableFormatter(header, fmts, min_width)
    print(tp.header)
    for i in range(10):
        dx = random.random()
        tangent = random.random()
        print(tp.line(i, dx, tangent))

if __name__ == "__main__":
    run()
