from datetime import timedelta
import sys
import time

from pysisyphus.helpers_pure import highlight_text
from pysisyphus.TablePrinter import TablePrinter


def run_perf(geom, calc_getter, pal_range, mem_range, repeat=1):
    assert repeat > 0

    try:
        pal_iter = range(*pal_range)
    except TypeError:
        pal_iter = (pal_range,)
    try:
        mem_iter = range(*mem_range)
    except TypeError:
        mem_iter = (mem_range,)

    results = {}
    for pal in pal_iter:
        print(f"{pal=} and ", end="")
        for mem in mem_iter:
            print(f"{mem=}:")
            for i in range(repeat):
                print(f"\tRunning cycle {i} ... ", end="")
                calc = calc_getter()
                calc.pal = pal
                calc.mem = mem
                start = time.time()
                _ = calc.get_forces(geom.atoms, geom.cart_coords)
                dur = time.time() - start
                td = timedelta(seconds=dur)
                print(f"finished in {td} h")
                sys.stdout.flush()
                key = (pal, mem)
                results.setdefault(key, list()).append(td)
            print()
    return results


def print_perf_results(results):
    pal_avgs = dict()
    pal_min = 1_000_000_000  # Dummy value
    for i, (key, tds) in enumerate(sorted(results.items())):
        pal, mem = key
        pal_min = min(pal, pal_min)
        print(highlight_text(f"{i:03d}: pal={pal}, mem={mem} MB"))
        for j, td in enumerate(tds):
            print(f"Run {j:03d}: {td} h")
        avg = sum(tds, timedelta(0)) / len(tds)
        # pal_avgs.setdefault(pal, list()).append(avg)
        pal_avgs[pal] = avg
        print(f"    Avg: {avg} h")
        print()

    avg_min = pal_avgs[pal_min]
    print(f"   pal_min: {pal_min} core(s)")
    print("efficiency: speedup / (pal / pal_min)\n")
    header = ("pal", "avg. / h", "speedup", "efficieny")
    col_fmts = ("int", "str", "float", "float")
    table = TablePrinter(header, col_fmts, width=20)
    table.print_header()
    for pal, avg in pal_avgs.items():
        speedup = avg_min / avg
        speedup_per_pal = speedup / (pal / pal_min)
        table.print_row((pal, str(avg), speedup, speedup_per_pal))
