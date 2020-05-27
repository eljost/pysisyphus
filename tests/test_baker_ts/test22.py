from pysisyphus.helpers import geom_loader

@pytest.mark.skip
def test_22():
    geom = geom_loader("lib:baker_ts/22_hconhoh.xyz", coord_type="redund")
    import pdb; pdb.set_trace()


if __name__ == "__main__":
    test_22()
