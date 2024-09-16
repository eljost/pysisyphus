from dataclasses import dataclass
from statistics import median
import time
from typing import Sequence


def render(dur: float) -> str:
    if dur < 60.0:
        fmt = "%Ss"
    elif dur < 3600.0:
        fmt = "%Mm %Ss"
    else:
        fmt = "%Hh %Mm %Ss"
    return time.strftime(fmt, time.gmtime(dur))


@dataclass
class RuntimeEstimation:
    mean_est: float
    median_est: float
    nremain: int

    def render_mean(self):
        return render(self.mean_est)

    def render_median(self):
        return render(self.median_est)


def estimate_runtime(durations: Sequence[float], ncalcs: int) -> RuntimeEstimation:
    ndurations = len(durations)
    if not durations:
        return RuntimeEstimation(0, 0, 0)
    elif ndurations > ncalcs:
        return RuntimeEstimation(float("nan"), float("nan"), 0)

    average_duration = sum(durations) / ndurations
    median_duration = median(durations)
    nremain = ncalcs - ndurations
    mean_estimation = average_duration * nremain
    median_estimation = median_duration * nremain
    return RuntimeEstimation(
        mean_est=mean_estimation,
        median_est=median_estimation,
        nremain=nremain,
    )
