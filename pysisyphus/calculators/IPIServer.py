import os
import socket

import numpy as np

from pysisyphus.calculators.Calculator import Calculator
from pysisyphus.socket_helper import (
    send_closure,
    recv_closure,
    get_fmts,
    EYE3,
)


class IPIServer(Calculator):
    def __init__(
        self,
        *args,
        address=None,
        host=None,
        port=None,
        unlink=True,
        hdrlen=12,
        **kwargs
    ):
        super().__init__(*args, **kwargs)
        self.address = address
        self.host = host
        self.port = port
        if self.host:
            assert self.port is not None
        self.hdrlen = hdrlen

        if self.address and unlink:
            self.unlink(self.address)

        if self.address:
            family = socket.AF_UNIX
            bind_args = (self.address, )
        else:
            family = socket.AF_INET
            bind_args = (self.host, self.port)

        # Create socket
        self.sock = socket.socket(family, socket.SOCK_STREAM)
        self.sock.bind(*bind_args)
        self.sock.listen(1)

        self._client_conn = None
        self._client_address = None
        self.fmts = None
        self.send_msg = None
        self.recv_msg = None

    def unlink(self, address):
        try:
            os.unlink(address)
        except OSError as err:
            if os.path.exists(address):
                raise err

    def listen_for(self, atoms, coords):
        atom_num = len(atoms)
        coords_num = len(coords)

        # Setup connection
        if (self._client_conn is None) or (self._client_address is None):
            self.log("Waiting for a connection.")
            self._client_conn, self._client_address = self.sock.accept()
            self.log("Got connection from {self._client_address}")
            # Create send/receive functions for this connection
            self.fmts = get_fmts(coords_num)
            self.send_msg = send_closure(self._client_conn, self.hdrlen, self.fmts)
            self.recv_msg = recv_closure(self._client_conn, self.hdrlen, self.fmts)

        # Reuse existing connection self._client_conn, wrapped in the
        # functions below.
        send_msg = self.send_msg
        recv_msg = self.recv_msg

        # Lets start talking
        send_msg("STATUS")
        ready = recv_msg()  # ready
        send_msg("STATUS")
        ready = recv_msg()  # ready
        send_msg("POSDATA")
        # Send cell vectors, inverse cell vectors, number of atoms and coordinates
        send_msg(EYE3, packed=True)  # cell vectors
        send_msg(EYE3, packed=True)  # inverse cell vectors
        send_msg(atom_num, fmt="int")
        send_msg(coords, fmt="floats")
        send_msg("STATUS")
        have_data = recv_msg()
        send_msg("GETFORCE")
        force_ready = recv_msg()

        energy = recv_msg(8, fmt="float")[0]
        client_atom_num = recv_msg(4, fmt="int")[0]
        assert atom_num == client_atom_num
        forces = recv_msg(coords_num * 8, fmt="floats")
        virial = recv_msg(72, fmt="nine_floats")
        zero = recv_msg(4, fmt="int")
        results = {
            "energy": energy,
            "forces": np.array(forces),
        }
        return results

    def cleanup(self):
        self.send_msg("STATUS")
        _ = self.recv_msg()
        self.send_msg("EXIT")
        self.conn.close()
        self.unlink(self.address)

    def get_energy(self, atoms, coords):
        return self.get_forces(atoms, coords)

    def get_forces(self, atoms, coords):
        result = self.listen_for(atoms, coords)
        return result
