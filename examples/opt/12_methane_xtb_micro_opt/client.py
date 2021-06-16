import struct
import socket

import numpy as np
import yaml

from pysisyphus.calculators import XTB
from pysisyphus.helpers import geom_loader

from pysisyphus.socket_helper import send_closure, recv_closure, get_fmts


def ipi_client(addr, atoms, forces_getter, hessian_getter=None, hdrlen=12):
    atom_num = len(atoms)
    # Number of entries in a Caretsian forces/coords vector
    cartesians = 3 * atom_num
    # Formats needed for struct.pack, to cast variables to bytes.
    # Bytes needed, to store a forces/coords vector
    floats_bytes = 8 * cartesians
    # Virial is hardcoded to the zero vector.
    VIRIAL = struct.pack("d" * 9, *np.zeros(9))
    ZERO = struct.pack("i", 0)

    fmts = get_fmts(cartesians)

    # Unix socket is hardcoded right now, but may also be inet-socket
    sock = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
    sock.connect(addr)

    send_msg = send_closure(sock, hdrlen, fmts)
    recv_msg = recv_closure(sock, hdrlen, fmts)

    counter = 0
    while True:
        try:
            # Lets start talking
            recv_msg(expect="STATUS")  # The server initiates with a STATUS
            send_msg("READY")

            status = recv_msg(expect="STATUS")
            if status == "NEEDPOS":
                send_msg("HAVEPOS")
                _ = recv_msg(4, fmt="int")[0]  # Recive atom num from IPI
                coords = np.array(recv_msg(floats_bytes, "floats"))  # Receive current coords
                assert coords.size % 3 == 0  # Assert Cartesian coordinates
                send_msg(atom_num, "int")
                # Just send back the current coordinates or translate all atoms in +X
                coords.reshape(-1, 3)[:, 0] += 1
                send_msg(coords, "floats")
                continue

            # When the optimization converged EXIT will be returned .. not documented!
            if status == "EXIT":
                print("Exited!")
                break

            # It seems we have to send READY two times ... not documented!
            send_msg("READY")
            # The server then returns POSDATA.
            recv_msg(expect="POSDATA")
            # Receive cell vectors, inverse cell vectors and number of atoms ...
            # but we don't use them here, so we don't even try to convert them to something.
            sock.recv(72)  # cell
            sock.recv(72)  # icell
            ipi_atom_num = recv_msg(4, fmt="int")[0]
            assert ipi_atom_num == atom_num
            # ... and the current coordinates.
            coords = np.array(recv_msg(floats_bytes, "floats"))

            recv_msg(expect="STATUS")
            # Indicate, that calculation is possible
            send_msg("HAVEDATA")
            get_what = recv_msg()
            # Acutal QC calculations
            if get_what == "GETFORCE":
                forces, energy = forces_getter(coords)
                print(f"Calculated energy & forces: {energy:.6f}, counter={counter}")
                send_msg("FORCEREADY")
                # Send everything to the server
                send_msg(energy, "float")
                send_msg(atom_num, "int")
                send_msg(forces, "floats")
                send_msg(VIRIAL, packed=True)
                # We don't want to send additional information, so just send 0.
                send_msg(ZERO, packed=True)
            elif get_what == "GETHESSIAN":
                hessian, energy = hessian_getter(coords)
                hessian_packed = struct.pack("d" * cartesians ** 2, *hessian.flatten())
                print(f"Calculated energy & Hessian: {energy:.6f}, counter={counter}")
                send_msg("HESSIANREADY")
                # Send everything to the server
                send_msg(energy, "float")
                send_msg(atom_num, "int")
                send_msg(hessian_packed, packed=True)

            counter += 1
        except Exception as err:
            raise err


def calc_ipi_client(addr, atoms, calc, **kwargs):
    def forces_getter(coords):
        results = calc.get_forces(atoms, coords)
        forces = results["forces"]
        energy = results["energy"]
        return forces, energy

    def hessian_getter(coords):
        results = calc.get_hessian(atoms, coords)
        hessian = results["hessian"]
        energy = results["energy"]
        return hessian, energy

    ipi_client(addr, atoms, forces_getter, hessian_getter, **kwargs)


yaml_fn = "12_methane_xtb_micro_opt.yaml"
with open(yaml_fn) as handle:
    run_dict = yaml.load(handle.read(), Loader=yaml.SafeLoader)
geom = geom_loader(run_dict["geom"]["fn"])
address = run_dict["calc"]["address"]
# address = "/tmp/ipi_h2o"  # For local testing with original I-PI server
calc = XTB()
calc_ipi_client(address, geom.atoms, calc)
