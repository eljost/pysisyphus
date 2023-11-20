from concurrent.futures import ProcessPoolExecutor
import cloudpickle


def apply_cloudpickle(fn, /, *args, **kwargs):
    fn = cloudpickle.loads(fn)
    return fn(*args, **kwargs)


class CloudpickleProcessPoolExecutor(ProcessPoolExecutor):
    """ProcessPoolExecutor using cloudpickle.

    From: https://stackoverflow.com/a/76008866"""

    def submit(self, fn, /, *args, **kwargs):
        return super().submit(apply_cloudpickle, cloudpickle.dumps(fn), *args, **kwargs)
