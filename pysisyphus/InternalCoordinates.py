#!/usr/bin/env python3

# [1] https://doi.org/10.1063/1.1515483

from scipy.spatial.distance import pdist, squareform

from pysisyphus.helpers import geom_from_library

def dist_mat(geom, factor=1.3):
    coords = geom.coords.reshape(-1, 3)
    # Condensed distance matrix
    cdm = pdist(coords)
    print(cdm)

if __name__ == "__main__":
    geom = geom_from_library("h2o.xyz") 
    dist_mat(geom)
