from concurrent.futures import ProcessPoolExecutor
from typing import Optional

import cloudpickle
from distributed import LocalCluster
import psutil


def apply_cloudpickle(fn, /, *args, **kwargs):
    fn = cloudpickle.loads(fn)
    return fn(*args, **kwargs)


class CloudpickleProcessPoolExecutor(ProcessPoolExecutor):
    """ProcessPoolExecutor using cloudpickle.

    From: https://stackoverflow.com/a/76008866"""

    def submit(self, fn, /, *args, **kwargs):
        return super().submit(apply_cloudpickle, cloudpickle.dumps(fn), *args, **kwargs)


class MaybeScheduler:
    def __init__(
        self,
        use_cluster: bool,
        n_workers: Optional[int] = None,
        threads_per_worker: int = 1,
    ):
        # TODO: support adress of existing scheduler
        if n_workers is None:
            n_workers = psutil.cpu_count(logical=False)

        if use_cluster:
            scheduler = LocalCluster(
                n_workers=n_workers, threads_per_worker=threads_per_worker
            )
        else:
            scheduler = None
        self.scheduler = scheduler

    def __enter__(self):
        return self.scheduler

    def __exit__(self, type, value, traceback):
        try:
            self.scheduler.close()
        except AttributeError:
            pass
