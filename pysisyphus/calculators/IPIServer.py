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
        max_retries=0,
        **kwargs
    ):
        super().__init__(*args, **kwargs)
        self.address = address
        self.host = host
        self.port = port
        if self.host:
            assert self.port is not None
        self.hdrlen = hdrlen
        self.max_retries = max_retries

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

        self.fmts = None
        self.send_msg = None
        self.recv_msg = None
        self.reset_client_connection()

    def unlink(self, address):
        try:
            os.unlink(address)
        except OSError as err:
            if os.path.exists(address):
                raise err

    def reset_client_connection(self):
        self.log("Resetting client connection info.")
        try:
            self.conn.close()
        except AttributeError:
            self.log("No client connection present.")
            pass
        self._client_conn = None
        self._client_address = None
        self.cur_retries = 0

    def listen_for(self, atoms, coords):
        atom_num = len(atoms)
        coords_num = len(coords)

        # Setup connection
        if (self._client_conn is None) or (self._client_address is None):
            self.log("Waiting for a connection.")
            self._client_conn, self._client_address = self.sock.accept()
            if self._client_address != "":
                conn_msg = f"Got new connection from {self._client_address}."
            else:
                conn_msg = "Got new connection."
            self.log(conn_msg)
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

    def retried_listen_for(self, atoms, coords):
        while self.cur_retries < self.max_retries:
            try:
                result = self.listen_for(atoms, coords)
                break
            except Exception as err:
                self.log(f"Caught exception: {err}.")
                self.cur_retries += 1
                self.reset_client_connection()
                result = self.listen_for(atoms, coords)
        return result

    def cleanup(self):
        self.send_msg("STATUS")
        _ = self.recv_msg()
        self.send_msg("EXIT")
        self.log("Sent EXIT to client.")
        self.reset_client_connection()
        # self.unlink(self.address)

    def get_energy(self, atoms, coords):
        return self.get_forces(atoms, coords)

    def get_forces(self, atoms, coords):
        if self.max_retries:
            result = self.retried_listen_for(atoms, coords)
        else:
            result = self.listen_for(atoms, coords)
        return result
