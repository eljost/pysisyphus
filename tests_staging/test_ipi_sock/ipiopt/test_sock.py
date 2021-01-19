import struct
import socket

import numpy as np

from pysisyphus.calculators import XTB
from pysisyphus.helpers import geom_loader

"""
In Pane 1:
    source ~/Code/i-pi/env.sh
    i-pi ipiopt.xml
In Pane 2
    module load xtb
    pytest -s -v
"""


def ipi_client(addr, atoms, forces_getter, hdrlen=12):
    atom_num = len(atoms)
    # Number of entries in a Caretsian forces/coords vector
    cartesians = 3 * atom_num
    # Formats needed for struct.pack, to cast variables to bytes.
    fmts = {
        "int": "i",
        "float": "d",
        "floats": "d" * cartesians,
    }
    # Bytes needed, to store a forces/coords vector
    floats_bytes = 8 * cartesians
    # Virial is hardcoded to the zero vector.
    VIRIAL = struct.pack("d" * cartesians, *np.zeros(cartesians))
    ZERO = struct.pack("i", 0)

    # Unix socket is hardcoded right now, but may also be inet-socket
    sock = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
    sock.connect(addr)

    def send_msg(msg, fmt=None, packed=False):
        """Send message through socket.

        If fmt is None and the message is not packed (converted to bytes),
        it is either convert to bytes via struct with the given fmt, or it
        is assumed to be a str and encoded as ascii. Strings are padded
        to a length of hdrlen (header length, default=12).

        If it is already packed the message is sent as is, e.g. for the
        virial, which is already in bytes (packed).
        """
        if not packed:
            if fmt == "floats":
                msg = struct.pack(fmts[fmt], *msg)
            elif (fmt is not None):
                msg = struct.pack(fmts[fmt], msg)
            else:
                msg = f"{msg: <{hdrlen}}".encode("ascii")
                print(f"SENDING: {msg}")
        sock.sendall(msg)

    def recv_msg(nbytes=None, fmt=None):
        """Receive a message of given length from the socket.

        If nbytes is not given we expect hdrlen (header length) bytes.
        If fmt is given the message is also unpacked using struct, otherwise
        ascii bytes are assumed and a string is returned.
        """
        if nbytes is None:
            nbytes = hdrlen

        msg = sock.recv(nbytes)
        if fmt:
            msg = struct.unpack(fmts[fmt], msg)
        else:
            msg = msg.decode("ascii").strip()
            print(f"RECEIVED: {msg}")
        return msg

    while True:
        try:
            # The server initiates everything with a STATUS
            status = recv_msg()
            send_msg("READY")
            status = recv_msg()
            # When the optimization converged EXIT will be returned .. not documented!
            if status == "EXIT":
                print("EXITING!")
                break

            # It seems we have to send READY two times ... not documented!
            send_msg("READY")
            # The server then returns POSDATA.
            posdata = recv_msg()
            # Receive cell vectors, inverse cell vectors and number of atoms ...
            cell = sock.recv(72)
            icell = sock.recv(72)
            ipi_atom_num = recv_msg(4, fmt="int")[0]
            assert ipi_atom_num == atom_num
            # ... and the current coordinates.
            coords = np.array(recv_msg(floats_bytes, "floats"))

            # Acutal QC calculation
            forces, energy = forces_getter(coords)

            status = recv_msg()
            # Indicate, that energy and forces are available
            send_msg("HAVEDATA")
            getforces = recv_msg()
            send_msg("FORCEREADY")
            # Send everything to the server
            send_msg(energy, "float")
            send_msg(atom_num, "int")
            send_msg(forces, "floats")
            send_msg(VIRIAL, packed=True)
            # We don't want to send additional information, so just send 0.
            send_msg(ZERO, packed=True)
        except Exception as err:
            raise err


def calc_ipi_client(addr, atoms, calc, **kwargs):

    def forces_getter(coords):
        results = calc.get_forces(atoms, coords)
        forces = results["forces"]
        energy = results["energy"]
        return forces, energy

    ipi_client(addr, atoms, forces_getter, **kwargs)


def test_ipi_client():
    addr = "/tmp/ipi_h2o"
    geom = geom_loader("h2o.xyz")
    calc = XTB()

    calc_ipi_client(addr, geom.atoms, calc) 
