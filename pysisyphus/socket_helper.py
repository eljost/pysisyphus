import struct

import numpy as np


NINE_ZEROS = struct.pack("d" * 9, *[0.0] * 9)
EYE3 = struct.pack("d" * 9, *np.eye(3).flatten())


def send_closure(sock, hdrlen, fmts, verbose=False):
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
            elif fmt is not None:
                msg = struct.pack(fmts[fmt], msg)
            else:
                msg = f"{msg: <{hdrlen}}".encode("ascii")
        if verbose:
            print(f"SENDING: {msg}")
        sock.sendall(msg)

    return send_msg


"""Vllt. recv_msg so verändern, dass es, wenn ein fmt gegeben ist auch
automatisch die korrekte anzahl an bytes erwartet"""


def recv_closure(sock, hdrlen, fmts, verbose=False):
    def recv_msg(nbytes=None, fmt=None, expect=""):
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
        if verbose:
            if expect:
                expect = f", expected '{expect}'"

            print(f"RECEIVED: {msg}{expect}")
        return msg

    return recv_msg


def get_fmts(cartesians):
    fmts = {
        "int": "i",  # Atom number
        "float": "d",  # Energy
        "nine_floats": "d" * 9,  # (Inverse) cell vectors, virial, zeros
        "floats": "d" * cartesians,  # Forces
        "floats_sq": "d" * cartesians**2,  # Hessian
    }
    return fmts
