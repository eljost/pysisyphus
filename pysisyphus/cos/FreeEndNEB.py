from pysisyphus.cos.AdaptiveNEB import AdaptiveNEB

# [1] https://www.pnas.org/content/pnas/104/9/3031.full.pdf
#     Zhu, 2006
#     Original method
# [2] http://dx.doi.org/10.1063/1.4962019
#     Zhang, 2016
#     FreeEnd Adaptive NEB


class FreeEndNEB(AdaptiveNEB):
    def __init__(self, *args, fix_first=False, fix_last=False, **kwargs):
        """Simple Free-End-NEB method.

        Only implements Eq. (7) from [2]. For other implementations
        please see the commit 01bc8812ca6f1cd3645d43e0337d9e3c5fb0ba55.
        """
        kwargs["adapt"] = False
        super().__init__(*args, fix_first=fix_first, fix_last=fix_last, **kwargs)
