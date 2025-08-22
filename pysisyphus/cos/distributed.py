from dataclasses import dataclass
from collections.abc import Callable
import math

import distributed
import numpy as np

from pysisyphus.Geometry import Geometry
from pysisyphus.helpers_pure import log


@dataclass
class ParallelData:
    ncores: int
    ncalcs: int
    nbatches: int
    batch_sizes: np.ndarray
    batch_pals: np.ndarray
    calc_pals: np.ndarray


def pal_values_for_parallel_calcs(ncores: int, ncalcs: int) -> list[int]:
    """Determine sensible pal values for parallel calculations.

    Given 'ncores' available CPU cores and 'ncalcs' calculation, determine suitable
    'pal' values for each calculation, so the calculations can run in parallel. If there
    are more calculations than cores, most calculations will run with pal=1, with a few
    remaining calculations with pal > 1.

    Parameters
    ----------
    ncores
        Positive integer, number of available cores.
    ncalcs
        Integer, number of calculations to carry out.

    Returns
    -------
    calc_pals
        List of positive integers.
    """
    assert ncores > 0
    assert ncalcs > 0

    # Most calculations will run with pal=1, i.e., one calculation per core.
    nbatches = math.ceil(ncalcs / ncores)

    # Run one calculation per core in a batch and the remaining calculations in the
    # last batch with potentially more cores.
    batch_sizes = np.full(nbatches, ncores)
    if (rest := ncalcs % ncores) != 0:
        batch_sizes[-1] = rest

    # The number of cores used per calculation in a batch are derived from the number
    # of available cores.
    batch_pals = ncores // batch_sizes
    # print(f"{ncalcs=}, {ncores=}, {nbatches=}, {batch_sizes=}, {batch_pals=}")

    # Distribute pal values of all batches and images
    calc_pals = np.repeat(batch_pals, batch_sizes).tolist()
    assert len(calc_pals) == ncalcs
    pal_data = ParallelData(
        ncores=ncores,
        ncalcs=ncalcs,
        nbatches=nbatches,
        batch_sizes=batch_sizes,
        batch_pals=batch_pals,
        calc_pals=calc_pals,
    )
    return pal_data


def distributed_calculations(
    client: distributed.Client,
    images: list[Geometry],
    func: Callable,
    logger=None,
) -> list[Geometry]:
    """Carray out distributed calculations via dask.

    func should return the modified image."""
    nimages = len(images)

    # Determine number of available CPU resources.
    ncores = 0
    scheduler_info = client.scheduler_info()
    for worker_data in scheduler_info["workers"].values():
        try:
            cpu = worker_data["resources"]["CPU"]
        except KeyError:
            cpu = 0
        ncores += cpu
    assert ncores > 0, "No 'CPU' resources available. Did you forget to specify them?"

    # Backup original pal values for later restoration.
    pals_backup = [image.calculator.pal for image in images]
    # Assert that all pal values are initially the same. In the end this is probably
    # not necessary, but if there are varying pal values a more elaborate strategy
    # may be warranted.
    pal0 = pals_backup[0]
    assert all([pal == pal0 for pal in pals_backup]), (
        "Image calculators have different pal values! This function only supports "
        "images calculators with the same pal value throughout!"
    )
    pal_data = pal_values_for_parallel_calcs(ncores, nimages)
    log(logger, pal_data)
    pals_parallel = pal_data.calc_pals

    futures = list()
    for image, pal_parallel in zip(images, pals_parallel):
        # Set potentially modified pal value
        image.calculator.pal = pal_parallel
        futures.append(
            client.submit(
                func,
                image,
                resources={
                    "CPU": pal_parallel,
                },
            )
        )
        # print(f"submitted task with {pal_parallel=}")
    calculated_images = client.gather(futures)

    # Set original pal values on all calculators
    for image, pal_org in zip(calculated_images, pals_backup):
        image.calculator.pal = pal_org

    return calculated_images
