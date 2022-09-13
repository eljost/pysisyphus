import numpy as np

from pysisyphus.io.cube import Cube
from pysisyphus.helpers import geom_loader


def test_cube(this_dir):
    geom = geom_loader("lib:h2o.xyz")
    npoints = np.full(3, 10)
    cube = Cube(
        atoms=geom.atoms,
        coords3d=geom.coords3d,
        origin=np.zeros(3),
        axes=np.eye(3),
        npoints=npoints,
        vol_data=np.ones(np.product(npoints)),
    )
    print(cube)
    out_fn = this_dir / "test.cub"
    cube.write(out_fn)

    _ = Cube.from_file(out_fn)
