import json
import socket

import numpy as np

from pysisyphus.calculators.Calculator import Calculator


class SocketCalc(Calculator):

    valid_requests = ("energy", "forces", "hessian")

    def __init__(self, *args, host="localhost", port=8080, **kwargs):
        super().__init__(*args, **kwargs)

        self.port = port
        self.host = host

    def listen_for(self, atoms, coords, request):
        request = request.lower()
        assert request.lower() in self.valid_requests, \
            f"Invalid request '{request}'! Valid requests are '{self.valid_requests}'."

        request_for = {
            "atoms": atoms,
            "coords": coords.tolist(),
            "request": request,
        }
        request_for = json.dumps(request_for).encode("utf-8")

        # Create socket
        sock = socket.socket()
        # Allow reuse
        sock.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
        sock.bind((self.host, self.port))

        sock.listen(0)
        while True:
            client, address = sock.accept()

            with client:
                self.log(f"Got connection from {address}")

                # Send atom, coordinates and the request type
                client.sendall(request_for)

                to_json = b""
                while True:
                    data = client.recv(1024)
                    to_json += data

                    # Read until linebreak is found
                    if b"\n" in data:
                        self.log("Found linebreaking. Stop listening.")
                        break

                # Try to parse received data as JSON
                results = dict()
                try:
                    to_json.decode("utf-8")
                    results = json.loads(to_json)
                except json.JSONDecodeError:
                    self.log("JSON decode error")

                # Check if all required fields are present in the results. The
                # energy has to be always present.
                if ("energy" in results) and (request in results):
                    self.log("All required fields are present in the received JSON. "
                             "Breaking.")
                    break
                else:
                    self.log("Could not parse received data as JSON!")

        sock.close()

        # energy has to be always present
        results = {key: results[key] for key in ("energy", request)}

        if "forces" in results:
            results["forces"] = np.array(results["forces"], dtype=float)

        if "hessian" in results:
            results["hessian"] = np.array(results["hessian"],
                                          dtype=float).reshape(-1, 3*len(atoms))
        return results

    def get_energy(self, atoms, coords):
        result = self.listen_for(atoms, coords, "energy")
        return result

    def get_forces(self, atoms, coords):
        result = self.listen_for(atoms, coords, "forces")
        return result

    def get_hessian(self, atoms, coords):
        result = self.listen_for(atoms, coords, "hessian")
        return result
