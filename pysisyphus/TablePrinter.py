import logging
import textwrap


class TablePrinter:
    def __init__(
        self,
        header,
        col_fmts,
        width=12,
        sub_underline=True,
        shift_left=0,
        fmts_update=None,
        mark="*",
        logger=None,
        level=logging.INFO,
    ):
        self.header = header
        self.col_fmts = col_fmts
        self.width = width
        self.sub_underline = sub_underline
        self.shift_left = shift_left
        if fmts_update is None:
            fmts_update = {}
        self.fmts_update = fmts_update
        self.mark = mark
        self.logger = logger
        self.level = level

        w = str(self.width - 1)

        self.fmts = {
            "int": "{:>" + w + "d}",
            "float": "{: >" + w + ".6f}",
            "float_short": "{: >" + w + ".2f}",
            "str": "{:>" + w + "s}",
            "mark": "{:>1s}",
        }
        self.fmts.update(fmts_update)
        if self.sub_underline:
            self.header = [h.replace("_", " ") for h in self.header]

        self.header_str = " ".join([h.rjust(self.width) for h in self.header])
        mark_fmt = self.fmts["mark"]
        self.conv_str = " ".join([self.fmts[fmt] + mark_fmt for fmt in self.col_fmts])
        h0_len = len(self.header[0])
        self.offset = self.width - h0_len
        self.prefix = " " * (self.offset - self.shift_left)
        self.sep = (
            self.prefix
            + "-" * (len(self.header_str) - self.width + h0_len)
            + "-" * abs(self.shift_left)
        )

    def _print(self, msg, level=None):
        if level is None:
            level = self.level
        if self.logger:
            self.logger.log(level, msg)
        else:
            print(msg)

    def print_sep(self):
        self._print(self.sep)

    def print_header(self, with_sep=True):
        self._print(self.header_str)
        if with_sep:
            self.print_sep()

    def print_row(self, args, marks=None):
        if marks is None:
            marks = ["" for _ in args]
        marked_args = list()
        for arg, to_mark in zip(args, marks):
            marked_args.append(arg)
            marked_args.append(self.mark if to_mark else " ")
        self._print(self.conv_str.format(*marked_args))

    def print(self, *args, **kwargs):
        text = " ".join([str(a) for a in args])
        try:
            level = kwargs["level"]
        except KeyError:
            level = 0
        level_prefix = "    " * level
        self._print(textwrap.indent(text, self.prefix + level_prefix))
