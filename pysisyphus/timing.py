from dataclasses import dataclass
from statistics import median
from typing import Sequence


def render(dur: float) -> str:
    dur = int(dur)
    day = ("d", 60 * 60 * 24)
    hour = ("h", 60 * 60)
    min = ("m", 60)
    sec = ("s", 1)
    tokens = list()
    append = False
    for key, divisor in (day, hour, min, sec):
        quot, dur = divmod(dur, divisor)
        if append or quot > 0:
            token = f"{quot}{key}"
            tokens.append(token)
            append = True
    return " ".join(tokens)


@dataclass
class RuntimeEstimation:
    mean_est: float
    median_est: float
    nremain: int

    def render_mean(self):
        return render(self.mean_est)

    def render_median(self):
        return render(self.median_est)


def estimate_runtime(
    durations: Sequence[float], ncalcs: int, pal: int = 1
) -> RuntimeEstimation:
    assert pal >= 1, f"pal must be a positive integer, but got {pal=}!"
    ndurations = len(durations)
    if not durations:
        return RuntimeEstimation(0, 0, 0)
    elif ndurations > ncalcs:
        return RuntimeEstimation(float("nan"), float("nan"), 0)

    average_duration = sum(durations) / ndurations
    median_duration = median(durations)
    nremain = ncalcs - ndurations
    mean_estimation = average_duration * nremain / pal
    median_estimation = median_duration * nremain / pal
    return RuntimeEstimation(
        mean_est=mean_estimation,
        median_est=median_estimation,
        nremain=nremain,
    )
