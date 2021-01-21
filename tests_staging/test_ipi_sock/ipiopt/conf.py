addr = "./h2o_sock"


def get_fmts(cartesians):
    fmts = {
        "int": "i",
        "float": "d",
        "nine_floats": "d" * 9,
        "floats": "d" * cartesians,
    }
    return fmts
