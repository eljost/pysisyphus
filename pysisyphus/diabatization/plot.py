import argparse
from collections.abc import Sequence
import sys
import warnings

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np

from pysisyphus.constants import AU2NU, AU2EV
from pysisyphus.diabatization import logger


EV2NU = AU2NU / AU2EV


def state_graph_from_en_mat(
    en_mat: np.ndarray, state_inds: Sequence[int], thresh_eV: float = 0.0
) -> nx.Graph:
    """Graph representation of energy matrix with state couplings as edges.

    Parameters
    ----------
    en_mat
        Quadratic energy matrix containing electronic energies in eV.
    state_inds
        Sequence of positive integers containing state labels.
    thresh_eV
        Positive floating point number that can be used to filter the
        couplings. If a coupling is below 'thresh_eV', no edge is created.
        Defaults to 0.0 eV, so by default all couplings are included.

    Returns
    -------
    G
        networkx.Graph representation of the energy matrix with states as
        nodes and couplings as edges. The node have an 'energy' attribute,
        containing the state energy in eV and the edges have an "weight"
        attribute, containing the electronic coupling in eV.
    """

    nrows, ncols = en_mat.shape
    assert nrows == ncols
    nstates = len(state_inds)
    assert nstates == nrows
    # Matrix must be symmetric
    np.testing.assert_allclose(en_mat, en_mat.T, atol=1e-12)

    # Create nodes from states
    state_ens = np.diag(en_mat)
    G = nx.Graph()
    for i, state_ind in enumerate(state_inds):
        energy = state_ens[i]
        G.add_node(state_ind, energy=energy)

    # Add couplingss as edges
    for i in range(nstates):
        state_i = state_inds[i]
        for j in range(i + 1, nstates):
            state_j = state_inds[j]
            coupling = abs(en_mat[i, j])
            if coupling > thresh_eV:
                coupling_meV = coupling * 1e3
                G.add_edge(state_i, state_j, weight=coupling_meV)
    return G


def map_array_to_interval(
    arr: np.ndarray, min_new: float, max_new: float, thresh=1e-12
) -> np.ndarray:
    """Map array from interval [array.min(), array.max()] to [min_new, max_new].

    Parameter
    ---------
    arr
        1d array containing floating point numbers.

    min_new
        Lower bound of the new interval.
    max_new
        Upper bound of the new interval.
    thresh
        When 'max_new - min_new' falls between this threshold an exception is
        raised.

    Returns
    -------
    mapped
        1d array w/ original shapped mapped onto the new interval.
    """

    spread_new = abs(max_new - min_new)
    if len(arr) == 1:
        mean = max(spread_new / 2.0, min_new)
        mapped = np.array((mean,))
        return mapped

    min_ = arr.min()
    max_ = arr.max()
    spread_org = max_ - min_

    if (thresh is not None) and (spread_new < thresh):
        raise Exception(
            f"'max_new - min_new = {spread_new: >10.4e}' < {thresh=: >10.4e}! Either set "
            f"thresh to None or decrease thresh below {spread_new: >10.4e}"
        )
    scale = spread_org / spread_new
    if scale <= 1e-14:
        warnings.warn(
            f"Obtained very small scaling factor {scale=: >10.4e}! "
            "Results may be unreliable/wrong!"
        )
    mapped = min_new + (arr - min_) / scale
    return mapped


def draw_state_graph(G: nx.Graph) -> plt.Figure:
    """Draw state graph."""
    fig, ax = plt.subplots(figsize=(10, 10))

    # Spring layout takes "weight" atributes of edges into account.
    pos = nx.spring_layout(G)
    nx.draw(
        G,
        pos=pos,
        ax=ax,
        node_color="white",
        edgecolors="purple",
        node_size=3000,
        linewidths=2,
    )

    # Choose edge widths according to coupling strength.
    #
    # Loop over all edges to determine the minimum and the maximum edge weight.
    min_weight = float("inf")
    max_weight = -1
    widths = list()
    for edge in G.edges:
        weight = G.get_edge_data(*edge)["weight"]
        min_weight = min(min_weight, weight)
        max_weight = min(max_weight, weight)
        widths.append(weight)
    # Map edge weights onto new interval [1, 12], so the highest coupling has
    # width 12 and the lowest coupling has width 1.
    # TODO: what to do when couplings are degenerate?
    widths = map_array_to_interval(np.array(widths), 1.0, 13.0)
    nx.draw_networkx_edges(G, pos, width=widths)

    # Draw node labels. Render one node per state w/ its energy.
    labels = nx.get_node_attributes(G, "energy")
    labels = {k: f"State {k}\n{v: >6.3} eV" for k, v in labels.items()}
    nx.draw_networkx_labels(G, pos, labels=labels, font_size=10)

    # Draw couplings as edge labels in meV
    edge_labels = nx.get_edge_attributes(G, "weight")
    edge_labels = {k: f"{v:.1f}" for k, v in edge_labels.items()}
    label_pos = 0.6
    nx.draw_networkx_edge_labels(G, pos, edge_labels, label_pos=label_pos)

    ax.set_title(f"Couplings in meV; labels at position {label_pos}.")
    return fig


def parse_args(args):
    parser = argparse.ArgumentParser()
    parser.add_argument("fn", help="Path to pysisyphus diabatiatzion result file.")
    parser.add_argument("--out-fn", default=None)
    parser.add_argument("--state-inds", nargs="+", type=int)
    return parser.parse_args(args)


def run():
    args = parse_args(sys.argv[1:])
    fn = args.fn
    out_fn = args.out_fn
    # TODO: read state_inds from npz file
    state_inds = args.state_inds

    logger.info(f"Rendering {fn}")

    data = np.load(fn)
    adia_ens = data["adia_ens"]
    try:
        state_inds = data["states"]
    except KeyError:
        # If state_inds was not set via --state-inds we enumerate them by ourselves,
        # starting from 0.
        if not state_inds:
            state_inds = list(range(len(adia_ens)))
    logger.info(f"Using {state_inds=}")

    with np.printoptions(precision=4, formatter={"float": lambda f: f"{f: >8.4f}"}):
        logger.info(f"Adiabatic energies: {adia_ens} eV")

    adia_mat = np.diag(adia_ens)
    U = data["U"]
    dia_mat = U.T @ adia_mat @ U
    G = state_graph_from_en_mat(dia_mat, state_inds=state_inds)
    fig = draw_state_graph(G)
    fig.suptitle(fn)
    fig.tight_layout()
    if out_fn is not None:
        fig.savefig(out_fn)
        logger.info(f"Saved figure to '{out_fn}'.")
    plt.show()


if __name__ == "__main__":
    run()
