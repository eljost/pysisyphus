import matplotlib.pyplot as plt
import networkx as nx
from networkx.algorithms import isomorphism
import numpy as np

from pysisyphus.helpers import geom_loader
from pysisyphus.intcoords.findbonds import get_bond_sets
from zmat import zmat_from_fn, geom_from_zmat


def geom_to_graph(geom):
    G = nx.Graph()
    G.add_nodes_from(
        [(i, {"atom": atom}) for i, atom in enumerate(geom.atoms)]
    )
    bonds = get_bond_sets(geom.atoms, geom.coords3d)
    c3d = geom.coords3d
    lengths = [np.linalg.norm(c3d[i]-c3d[j]) for i, j in bonds]
    bonds_with_lengths = [(i, j, {"length": length, })
                          for (i, j), length in zip(bonds, lengths)]
    G.add_edges_from(bonds_with_lengths)

    return G


def draw_graph(G):
    # labels = nx.get_node_attributes(G, "atom") 
    # nx.draw(G, labels=labels, font_weight="bold")
    nx.draw(G, with_labels=True, font_weight="bold")


def node_match(n1, n2):
    return n1["atom"] == n2["atom"]


def bond_dict(G):
    return {
        (i, j): G.edges()[i, j]["length"]
        for (i, j) in G.edges()
    }


def find_subgraph(G1, G2):
    GM = isomorphism.GraphMatcher(G1, G2, node_match=node_match)
    subgraphs = list(GM.subgraph_isomorphisms_iter())
    assert len(subgraphs) > 0

    if len(subgraphs) == 1:
        return subgraphs[0]

    ref_sg = G1.subgraph(subgraphs[0].keys())
    ref_bonds = bond_dict(ref_sg)
    g2_bonds = bond_dict(G2)
    diffs = np.zeros(len(subgraphs))
    for i, map_ in enumerate(subgraphs):
        def translate_edge(i, j):
            return (map_[i], map_[j])
        for edge, length in ref_bonds.items():
            edge_ = translate_edge(*edge)
            diffs[i] += abs(g2_bonds[edge_] - length)
    return subgraphs[diffs.argmin()]


def hydrogen_bond(G, edge):
    i, j = edge
    return (G.nodes[i]["atom"] == "H") or (G.nodes[j]["atom"] == "H")


def missing_hydrogens(G1, G2, sg):
    h_missing = dict()
    total = 0
    for g1, g2 in sg.items():
        e1 = G1.edges(g1)
        e2 = G2.edges(g2)
        e1H = [e for e in e1 if hydrogen_bond(G1, e)]
        a1 = G1.nodes[g1]["atom"]
        a2 = G2.nodes[g2]["atom"]
        h_missing[g2] = len(e1H)
        total += len(e1H)
        print(f"G1: {a1}, G2: {a2}, missing {len(e1H)} hydrogens")
    return h_missing, total


def add_hydrogens(geom, sg, htot, zmat):
    coords3d = np.zeros((len(geom.atoms) + htot, 3))
    print(coords3d.shape)
    # Copy atomic coordinates of heavy atoms at the correct positions
    present = list()
    for g1, g2 in sg.items():
        print(g1, g2)
        coords3d[g1] = geom.coords3d[g2]
        present.append(g1)
    # Assert that there are no "holes".
    # We assume that the heave atoms come before the hydrogens in the Z-Matrix.
    assert set(present) == set(range(len(present)))
    start_at = len(present)
    geom = geom_from_zmat(zmat, coords3d=coords3d, start_at=start_at)
    return geom


# Glycine, reference, with hydrogens 
zmat = zmat_from_fn("glycine_resorted.zmat")
gly = geom_loader("glycine_resorted.xyz")
G = geom_to_graph(gly)
# draw_graph(G)
# plt.show()

# Glycine, without hydrogens
glynoh = geom_loader("glycine_noh.xyz")
Gnh = geom_to_graph(glynoh)
# draw_graph(Gnh)
# plt.show()

sg = find_subgraph(G, Gnh)
print(sg)
# sg0 = sgs[0]
# sg0 = sgs[1]
hms0, htot0 = missing_hydrogens(G, Gnh, sg)
print(hms0, htot0)

geom = add_hydrogens(glynoh, sg, htot0, zmat)

# from pysisyphus.calculators.XTB import XTB
# from pysisyphus.optimizers.PreconLBFGS import PreconLBFGS
# geom.set_calculator(XTB())
# opt = PreconLBFGS(geom)
# opt.run()
