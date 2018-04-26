#!/usr/bin/env python3

from pathlib import Path
from pprint import pprint
import struct

from natsort import natsorted
import numpy as np


np.set_printoptions(suppress=True, precision=4)


def find_ccre(path):
    #return [ccre for ccre in path.glob("CCRE*")]
    return natsorted(path.glob("CCRE*"))


def parse_ccre(ccre, mos_num, occ_num):
    """Adapted from TheoDORE 1.6 lib/file_parser.py"""
    handle = open(ccre, "rb")

    handle.read(8)
    method = struct.unpack('8s', handle.read(8))[0]
    handle.read(8)
    nentry = struct.unpack('i', handle.read(4))[0]
    handle.read(20)

    virt_num = mos_num - occ_num

    assert(nentry % virt_num == 0)
    active_num = nentry // virt_num
    frozen_num = occ_num - active_num

    shape = (occ_num-frozen_num, virt_num)
    tden = np.empty(shape)
    assert tden.size == nentry

    for occ_ind in range(occ_num-frozen_num):
        for virt_ind in range(mos_num-occ_num):
            coeff = struct.unpack('d', handle.read(8))[0]
            if abs(coeff) >= .1:
                print(f"[{occ_ind}, {virt_ind}] = {coeff}")
            tden[occ_ind, virt_ind] = coeff
    lbytes = struct.unpack('4s', handle.read(4))


    handle.close()
    pprint(locals())
    #print(tden[-1])
    #print(ccre)

def run():
    p = Path(".")
    ccres = find_ccre(p)
    mos_num = 86
    occ_num = 15
    parse_ccre(ccres[0], mos_num, occ_num)
    pass


if __name__ == "__main__":
    run()
