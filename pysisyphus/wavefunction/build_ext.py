from cffi import FFI


def load_h(h_fn):
    with open(h_fn) as handle:
        text = handle.read().strip()
    return text


def run():
    ffibuilder = FFI()
    dir_ = "devel_ints/"
    names = ("ovlp3d", "dipole3d", "diag_quadrupole3d", "quadrupole3d", "kinetic3d")

    all_sources = list()
    all_includes = list()
    all_cdefs = list()
    for name in names:
        c_name = dir_ + name + ".c"
        all_sources.append(c_name)
        h_name = dir_ + name + ".h"
        all_includes.append(h_name)
        cdefs = load_h(h_name)
        all_cdefs.append(cdefs)

    cdefs = "\n\n".join(all_cdefs)
    includes = "\n".join([f'#include "{h_name}"' for h_name in all_includes])

    ffibuilder.cdef(cdefs)
    ffibuilder.set_source(
        "_ints1el",
        includes,
        sources=all_sources,
        libraries=[
            "m",
        ],
    )
    ffibuilder.compile(verbose=True)


if __name__ == "__main__":
    run()
