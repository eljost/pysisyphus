import matplotlib.pyplot as plt
import numpy as np
import requests
import simplejson as json

from pysisyphus.calculators.AnaPot import AnaPot
from pysisyphus.server.optimizer import OptState


def run():
    geom = AnaPot.get_geom((-1, 3, 0))
    print(geom)

    def prep_post(geom):
        return {
            "coords": geom.coords.tolist(),
            "gradient": geom.gradient.tolist(),
            "energy": float(geom.energy),
        }

    req_kwargs = {
        "headers": {"Content-type": "application/json", },
    }
    base_url = "http://localhost:8080/"
    def send_json(path, data):
        url = base_url + path
        response = requests.post(url=url, data=json.dumps(data), **req_kwargs)
        return response


    init = {
        "key": "cg",
        "alpha": 0.1,
        "gdiis": False,
        "line_search": False,
    }
    init_resp = send_json("init", init)
    print("Init response", init_resp)

    coords = list()
    opt_state = OptState(alpha=0.3, key="cg")

    for i in range(50):
        coords.append(geom.coords)
        print("\tcur. coords", geom.coords)
        gradient = geom.gradient
        norm = np.linalg.norm(gradient)

        print(f"{i:02d}: norm(grad)={norm:.6f}")
        if norm <= 0.1:
            print("Converged")
            break

        pp = prep_post(geom)
        # # resp = requests.post(**req_kwargs, data=pp)
        # resp = send_json("step", pp)
        # resp_data = resp.json()

        # step = resp_data["step"]
        # status = resp_data["status"]
        step, status = opt_state.step(geom.coords, geom.energy, geom.gradient)
        new_coords = geom.coords + step
        # print("\t", status)
        geom.coords = new_coords
        print()
        print()

        # import pdb; pdb.set_trace()
        # pass

    coords = np.array(coords)
    calc = geom.calculator
    calc.plot()
    ax = calc.ax
    ax.plot(*coords.T[:2], "ro-")
    for i, xy in enumerate(coords[:,:2]):
        ax.annotate(str(i), xy)
    plt.show()


if __name__ == "__main__":
    run()
