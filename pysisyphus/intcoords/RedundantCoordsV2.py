from pysisyphus.intcoords import RedundantCoords as RC


class RedundantCoords(RC):

    def __init__(self, *args, lb_min_deg=170, **kwargs):
        super().__init__(*args, lb_min_deg=lb_min_deg, **kwargs)
