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
    listen_kinds = ("coords", "energy", "forces", "hessian")

    def __init__(
        self,
        *args,
        address=None,
        host=None,
        port=None,
        unlink=True,
        hdrlen=12,
        max_retries=0,
        verbose=False,
        **kwargs,
    ):
        super().__init__(*args, **kwargs)
        self.address = address
        self.host = host
        self.port = port
        if self.host:
            assert self.port is not None
        self.hdrlen = hdrlen
        self.max_retries = max_retries
        self.verbose = verbose

        if self.address and unlink:
            self.unlink(self.address)

        if self.address:
            family = socket.AF_UNIX
            bind_args = (self.address,)
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
        self._client_conn = None
        self._client_address = None
        self.cur_retries = 0

    def listen_for_client_atom_num(self, atom_num):
        client_atom_num = self.recv_msg(4, fmt="int")[0]
        self.log(
            f"Client sent number of atoms: {client_atom_num}, expecting {atom_num}."
        )
        assert atom_num == client_atom_num
        return atom_num

    def listen_for_energy(self):
        send_msg = self.send_msg
        recv_msg = self.recv_msg

        send_msg("GETENERGY")
        recv_msg(expect="ENERGYREADY")
        energy = recv_msg(8, fmt="float")[0]
        results = {
            "energy": energy,
        }
        return results

    def listen_for_forces(self, atom_num):
        send_msg = self.send_msg
        recv_msg = self.recv_msg

        send_msg("GETFORCE")
        recv_msg(expect="FORCEREADY")
        energy = recv_msg(8, fmt="float")[0]
        self.listen_for_client_atom_num(atom_num)
        coord_num = 3 * atom_num
        forces = recv_msg(coord_num * 8, fmt="floats", expect="forces")
        recv_msg(72, fmt="nine_floats", expect="virial")
        recv_msg(4, fmt="int", expect="zero")
        results = {
            "energy": energy,
            "forces": np.array(forces),
        }
        return results

    def listen_for_hessian(self, atom_num):
        send_msg = self.send_msg
        recv_msg = self.recv_msg

        send_msg("GETHESSIAN")
        recv_msg(expect="HESSIANREADY")
        energy = recv_msg(8, fmt="float")[0]
        self.listen_for_client_atom_num(atom_num)
        coord_num = 3 * atom_num
        hessian = recv_msg(coord_num ** 2 * 8, fmt="floats_sq", expect="Hessian")
        hessian = np.array(hessian).reshape(-1, coord_num)
        results = {
            "energy": energy,
            "hessian": hessian,
        }
        return results

    def listen_for(self, atoms, coords, kind="forces"):
        assert kind in self.listen_kinds

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
            self.send_msg = send_closure(
                self._client_conn, self.hdrlen, self.fmts, verbose=self.verbose
            )
            self.recv_msg = recv_closure(
                self._client_conn, self.hdrlen, self.fmts, verbose=self.verbose
            )

        # Reuse existing connection self._client_conn, wrapped in the
        # functions below.
        send_msg = self.send_msg
        recv_msg = self.recv_msg

        # Lets start talking
        send_msg("STATUS")
        recv_msg(expect="READY")

        # This path handles a coordinate update through the client.
        if kind == "coords":
            send_msg("NEEDPOS")
            # TODO: allow skipping the update
            recv_msg(expect="HAVEPOS")
            # Send current atom number and coordinates
            send_msg(atom_num, fmt="int")
            send_msg(coords, fmt="floats")
            # Receive atom number and potentially modified coordinates from the client.
            self.listen_for_client_atom_num(atom_num)
            new_coords = recv_msg(atom_num * 3 * 8, fmt="floats")
            results = {"coords": np.array(new_coords)}
        # The path below leads to sending of coordinates and calculation of
        # energy and maybe its derivatives by the client.
        else:
            send_msg("STATUS")
            recv_msg(expect="READY")
            send_msg("POSDATA")
            # Send cell vectors, inverse cell vectors, number of atoms and coordinates
            send_msg(EYE3, packed=True)  # cell vectors
            send_msg(EYE3, packed=True)  # inverse cell vectors
            send_msg(atom_num, fmt="int")
            send_msg(coords, fmt="floats")
            send_msg("STATUS")
            recv_msg(expect="HAVEDATA")
            if kind == "energy":
                results = self.listen_for_energy()
            elif kind == "forces":
                results = self.listen_for_forces(atom_num)
            elif kind == "hessian":
                results = self.listen_for_hessian(atom_num)
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

    def __del__(self):
        self.send_msg("STATUS")
        _ = self.recv_msg()
        self.send_msg("EXIT")
        self.log("Sent EXIT to client.")
        self.reset_client_connection()
        # self.unlink(self.address)

    def get_energy(self, atoms, coords):
        return self.listen_for(atoms, coords, kind="energy")

    # def get_forces(self, atoms, coords):
    # if self.max_retries:
    # result = self.retried_listen_for(atoms, coords)
    # else:
    # result = self.listen_for(atoms, coords)

    # def get_coords(self, atoms, coords):
        # return self.listen_for(atoms, coords, kind="coords")

    def get_forces(self, atoms, coords):
        return self.listen_for(atoms, coords, kind="forces")

    def get_hessian(self, atoms, coords):
        return self.listen_for(atoms, coords, kind="hessian")
