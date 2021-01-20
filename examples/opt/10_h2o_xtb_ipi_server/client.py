import struct
import socket

import numpy as np
import yaml

from pysisyphus.calculators import XTB
from pysisyphus.helpers import geom_loader

from pysisyphus.socket_helper import send_closure, recv_closure, get_fmts


def ipi_client(addr, atoms, forces_getter, hdrlen=12):
    atom_num = len(atoms)
    # Number of entries in a Caretsian forces/coords vector
    cartesians = 3 * atom_num
    # Formats needed for struct.pack, to cast variables to bytes.
    # Bytes needed, to store a forces/coords vector
    floats_bytes = 8 * cartesians
    # Virial is hardcoded to the zero vector.
    VIRIAL = struct.pack("d" * 9, *np.zeros(cartesians))
    ZERO = struct.pack("i", 0)

    fmts = get_fmts(cartesians)

    # Unix socket is hardcoded right now, but may also be inet-socket
    sock = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
    sock.connect(addr)

    send_msg = send_closure(sock, hdrlen, fmts)
    recv_msg = recv_closure(sock, hdrlen, fmts)

    while True:
        try:
            # The server initiates everything with a STATUS
            status = recv_msg()
            send_msg("READY")
            status = recv_msg()
            # When the optimization converged EXIT will be returned .. not documented!
            if status == "EXIT":
                print("Exited!")
                break

            # It seems we have to send READY two times ... not documented!
            send_msg("READY")
            # The server then returns POSDATA.
            posdata = recv_msg()
            # Receive cell vectors, inverse cell vectors and number of atoms ...
            # but we don't use them here, so we don't even try to convert them to something.
            cell = sock.recv(72)
            icell = sock.recv(72)
            ipi_atom_num = recv_msg(4, fmt="int")[0]
            assert ipi_atom_num == atom_num
            # ... and the current coordinates.
            coords = np.array(recv_msg(floats_bytes, "floats"))

            # Acutal QC calculation
            forces, energy = forces_getter(coords)
            print(f"Calculated energy: {energy:.6f}")

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


yaml_fn = "10_h2o_xtb_ipi_server.yaml"
with open(yaml_fn) as handle:
    run_dict = yaml.load(handle.read(), Loader=yaml.SafeLoader)
geom = geom_loader(run_dict["geom"]["fn"])
address = run_dict["calc"]["address"]
calc = XTB()
calc_ipi_client(address, geom.atoms, calc)
