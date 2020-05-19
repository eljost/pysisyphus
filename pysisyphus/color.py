def red(str_):
    return f"\033[91m{str_}\033[0m"


def green(str_):
    return f"\033[92m{str_}\033[0m"


BOOL_COLORS = {
    True: green,
    False: red,
}


def bool_color(bool_, str_=None):
    if str_ is None:
        str_ = str(bool_)
    color = BOOL_COLORS[bool_]
    return color(str_)
