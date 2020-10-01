#!/usr/bin/env/python3

from pprint import pprint

import matplotlib.pyplot as plt
import networkx as nx
from networkx.algorithms import isomorphism

from pysisyphus.intcoords.setup import get_bond_mat
from pysisyphus.helpers import geom_from_xyz_file


def get_labels(geom):
    return {i: f"{i:02d}_{atom}" for i, atom in enumerate(geom.atoms)}


def run():
    fn = "07_ref_rx_phosphine_def2tzvp_reopt.xyz"
    geom = geom_from_xyz_file(fn)
    bm = get_bond_mat(geom)
    print(bm)
    node_attrs = {i: {"atom": atom} for i, atom in enumerate(geom.atoms)}
    g = nx.from_numpy_array(bm)
    nx.set_node_attributes(g, node_attrs)
    # import pdb; pdb.set_trace()

    # fig, ax = plt.subplots()
    # draw_kwargs = {
        # "ax": ax,
        # "with_labels": True,
        # "node_size": 250,
    # }
    # nx.draw(g, labels=get_labels(geom), **draw_kwargs)
    # plt.show()

    prod_fn = "01_ref_rx_product_opt.xyz"
    prod = geom_from_xyz_file(prod_fn)
    pbm = get_bond_mat(prod)
    gp = nx.from_numpy_array(pbm)
    pnode_attrs = {i: {"atom": atom} for i, atom in enumerate(prod.atoms)}
    nx.set_node_attributes(gp, pnode_attrs)

    # fig, ax = plt.subplots()
    # draw_kwargs["ax"] = ax
    # nx.draw(gp, labels=get_labels(prod), **draw_kwargs)
    # plt.show()

    gm = isomorphism.GraphMatcher(gp, g)
    si = gm.subgraph_is_isomorphic()
    sims = list(gm.subgraph_isomorphisms_iter())
    llens = [len(_) for _ in sims]
    pprint(sims)
    print(llens)
    ms = [i for i, d in enumerate(sims)
          if all([i == j for i, j in d.items()])
    ]
    mapping = sims[ms[0]]
    pprint(mapping)
    # import pdb; pdb.set_trace()
    # for mapping in sims:
    # import pdb; pdb.set_trace()

    pass


if __name__ == "__main__":
    run()
