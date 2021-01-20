import os
import socket

import numpy as np

from pysisyphus.socket_helper import send_closure, recv_closure, get_fmts

from pysisyphus.calculators.Calculator import Calculator

CELL = "CELL    " * 9
ICELL = "ICELL   " * 9


class IPIServer(Calculator):
    def __init__(self, address, *args, family=socket.AF_UNIX, unlink=True, hdrlen=12, **kwargs):
        super().__init__(*args, **kwargs)
        self.address = address
        self.hdrlen = hdrlen

        if unlink:
            try:
                os.unlink(self.address)
            except OSError as err:
                if os.path.exists(self.address):
                    raise err

        # Create socket
        self.sock = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
        self.sock.bind(self.address)
        self.sock.listen(1)

        self._client_conn = None
        self._client_address = None
        self.send_msg = None
        self.recv_msg = None

    def listen_for(self, atoms, coords):
        atom_num = len(atoms)
        coords_num = len(coords)
        fmts = get_fmts(coords_num)

        if (self._client_conn is None) or (self._client_address is None):
            self.log("Waiting for a connection.")
            self._client_conn, self._client_address = self.sock.accept()
            self.log("Got connection from {self._client_address}")
            # Create send/receive functions for this connection
            self.send_msg = send_closure(self._client_conn, self.hdrlen, fmts)
            self.recv_msg = recv_closure(self._client_conn, self.hdrlen, fmts)

        # Reuse existing connection
        send_msg = self.send_msg
        recv_msg = self.recv_msg
        send_msg("STATUS")
        ready = recv_msg()
        send_msg("STATUS")
        ready = recv_msg()
        send_msg("POSDATA")
        # Send cell vectors, inverse cell vectors, number of atoms and coordinates
        send_msg(CELL)
        send_msg(ICELL)
        send_msg(atom_num, fmt="int")
        send_msg(coords, fmt="floats")
        send_msg("STATUS")
        have_data = recv_msg()
        send_msg("GOTFORCES")
        force_ready = recv_msg()

        energy = recv_msg(8, fmt="float")[0]
        atom_num_ = recv_msg(4, fmt="int")[0]
        forces = recv_msg(coords_num * 8, fmt="floats")
        virial = recv_msg(72, fmt="nine_floats")
        zero = recv_msg(4, fmt="int")
        results = {
            "energy": energy,
            "forces": np.array(forces),
        }
        return results
        # send_msg("STATUS")
        # ready_ = recv_msg()
        # send_msg("EXIT")
        # conn.close()

    def get_energy(self, atoms, coords):
        return self.get_forces(atoms, coords)

    def get_forces(self, atoms, coords):
        result = self.listen_for(atoms, coords)
        return result
