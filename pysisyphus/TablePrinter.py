import logging
import textwrap

import pyparsing as pp


lit = pp.Literal
start = lit("{:")
end = lit("}")
alignment = lit(">") | lit("<") | lit("^")
sign = lit("+") | lit("-") | lit(" ")
int_num = pp.Word(pp.nums).set_parse_action(lambda s, l, t: int(t[0]))
width = int_num("width")
precision = pp.Suppress(lit(".")) + int_num("precision")
type_ = (
    lit("s")
    | lit("b")
    | lit("c")
    | lit("d")
    | lit("o")
    | lit("x")
    | lit("X")
    | lit("e")
    | lit("E")
    | lit("f")
    | lit("F")
    | lit("g")
    | lit("G")
    | lit("n")
    | lit("%")
)

fparser = (
    start
    + pp.Optional(alignment)("alignment")
    + pp.Optional(sign)("sign")
    + pp.Optional(width)
    + pp.Optional(precision)
    + pp.Optional(type_)("type")
    + end
)


def parse_fstring(fstr) -> dict:
    result = fparser.parse_string(fstr)
    return result.as_dict()


def res_to_fstring(res: dict):
    alignment = res.get("alignment", "")
    sign = res.get("sign", "")
    width = res.get("width", "")
    precision = res.get("precision", "")
    if precision:
        precision = f".{precision}"
    type_ = res["type"]
    return f"{{:{alignment}{sign}{width}{precision}{type_}}}"


def center(string, width):
    length = len(string)
    whitespace = width - length
    if whitespace >= 2:
        before = " " * (whitespace // 2)
        after = before
        centered = f"{before}{string}{after}"
    else:
        centered = string
    return centered


class TablePrinter:
    def __init__(
        self,
        header,
        col_fmts,
        width=12,
        sub_underline=True,
        mark="*",
        offset=8,
        logger=None,
        level=logging.INFO,
    ):
        self.header = header
        self.col_fmts = col_fmts
        self.width = width
        self.sub_underline = sub_underline
        self.mark = mark
        self.logger = logger
        self.level = level

        assert len(header) == len(col_fmts), (
            f"Number of {len(header)} fields does not match the number of "
            f"column formats {len(col_fmts)}"
        )

        w = str(self.width)  # Shortcut
        whalf = str(self.width // 2)
        self.fmts = {
            "str": "{:>" + w + "s}",
            "mark": "{:>1s}",
            #
            "int": "{:>" + w + "d}",
            "int_short": "{:>" + whalf + "d}",
            "int3": "{: >3d}",
            #
            "float": "{: >" + w + ".6f}",
            "float_short": "{: >" + whalf + ".3f}",
            #
            "complex_tdm": "{: >" + w + ".4f}",
            "complex_short": "{: >" + w + ".2f}",
        }
        if self.sub_underline:
            self.header = [h.replace("_", " ") for h in self.header]

        # Determine alignments and widths from given formats
        fmts = list()
        alignments = list()
        widths = list()
        for i, col_fmt in enumerate(self.col_fmts):
            # First, try to look up the format in our prepopulated
            # dictionary.
            try:
                fmt = self.fmts[col_fmt]
            # Second, if not found in the dict, we assume that the
            # string is a valid f-string.
            except KeyError:
                fmt = col_fmt
            # In any case, we parse the given string to determine
            # alignment and width.
            res = parse_fstring(fmt)
            alignment = res.get("alignment", "<")
            alignments.append(alignment)
            # Check if header is longer than the desired width.
            header_width = len(self.header[i])
            width = max(res.get("width", self.width), header_width)
            # Update the fstring with the (possibly) new width and reconstruct it
            res["width"] = width
            fmt = res_to_fstring(res)
            widths.append(width)
            fmts.append(fmt)

        # TODO: Check if header fields are longer than the given width.
        # If so, update the widths in the formats.

        mark_fmt = self.fmts["mark"]
        # self.conv_str will be used to render the given fields in a row
        self.conv_str = " ".join([fmt + mark_fmt for fmt in fmts])

        # Distance from the line start
        self.offset = offset
        # Whitespace that is prepended on every row
        self.prefix = " " * self.offset
        # Length of a given line w/o offset/prefix
        self.line_length = sum(widths) + len(widths) + len(widths) - 1
        # Separator string
        self.sep = self.prefix + "-" * self.line_length

        header_fmts = list()
        for alignment, width in zip(alignments, widths):
            hfmt = "{:" + alignment + str(width) + "}"
            header_fmts.append(hfmt)
        # Join with 2 spaces as we also have the mark field
        header_fmt = "  ".join(header_fmts)
        self.header_str = self.prefix + header_fmt.format(*self.header)

    @property
    def nfields(self) -> int:
        return len(self.header)

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
        row = self.prefix + self.conv_str.format(*marked_args)
        self._print(row)

    def print(self, *args, **kwargs):
        text = " ".join([str(a) for a in args])
        try:
            level = kwargs["level"]
        except KeyError:
            level = 0
        level_prefix = "    " * level
        self._print(textwrap.indent(text, self.prefix + level_prefix))

    def print_rows(self, all_args, marks=None, first_n=None):
        assert len(all_args) == self.nfields
        nrows = len(all_args[0])
        # Verify that provided args all have the same length
        assert all([len(arg) == nrows for arg in all_args])
        # Update the number of rows that will be printed
        if first_n is not None:
            nrows = min(nrows, first_n)

        for i in range(nrows):
            args_i = [arg[i] for arg in all_args]
            self.print_row(args_i, marks=marks)
