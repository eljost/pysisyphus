import warnings

import matplotlib.pyplot as plt
import matplotlib.colors
import matplotlib.patches
import matplotlib.patheffects
import networkx as nx
import numpy as np

from pysisyphus.elem_data import COVALENT_RADII, CPK_RGB, get_tm_indices
from pysisyphus.drivers.spectrum import spectrum_from_ens_fosc
from pysisyphus.intcoords.setup import get_fragments, get_bond_sets


def get_frags_noh(atoms, coords, bond_inds):
    tm_inds = get_tm_indices(atoms)
    h_inds = [i for i, atom in enumerate(atoms) if atom.lower() == "h"]
    frags = get_fragments(
        atoms, coords, bond_inds=bond_inds, ignore_atom_inds=tm_inds + h_inds
    )
    frags = [[tmi] for tmi in tm_inds] + list(map(list, frags))
    warnings.warn(
        "This still misses unconnected single atoms, e.g., a chlorine ligand!"
    )
    return frags


def plot_geometry_graph(ax, atoms, coords3d, frags, bond_inds, scale=12.5):
    node_color = list()
    node_size = list()
    pos0 = dict()
    G = nx.Graph()
    for i, atom in enumerate(atoms):
        node_color.append(CPK_RGB[atom])
        node_size.append(COVALENT_RADII[atom])
        if atom == "h":
            continue
        pos0[i] = coords3d[i, :2]  # Drop dimension w/ lowest variance?
        G.add_node(i, atom=atom)
    node_color = np.array(node_color)
    node_size = np.array(node_size)
    node_size /= node_size.min()
    node_size *= 350

    for from_, to_ in bond_inds:
        from_atom = atoms[from_]
        to_atom = atoms[to_]
        if from_atom == "h" or to_atom == "h":
            continue
        bond_length = np.linalg.norm(coords3d[from_] - coords3d[to_])
        G.add_edge(from_, to_, weight=bond_length)

    labels = {
        i: f"{atom.capitalize()}{i}" for i, atom in enumerate(atoms) if atom != "h"
    }
    pos = nx.kamada_kawai_layout(G, scale=scale, pos=pos0)
    """
    nx.draw(
        G, pos=pos, ax=ax, labels=labels, node_color=node_color, node_size=node_size
    )
    """
    frag_nodes = list()
    frag_pos = list()
    # Each fragment is drawn separately as the initial idea was to color the fragment
    # nodes different, depending on whether they loose or receive electron density.
    # Now arrows are used and drawing each fragment separately is actually not required.
    # But it doesn't hurt either to keep it like this.
    for fnodes in frags:
        fnode_color = node_color[fnodes]
        fnode_size = node_size[fnodes]
        fpos = sum([pos[n] for n in fnodes]) / len(fnodes)
        fnodes = nx.draw_networkx_nodes(
            G,
            pos=pos,
            ax=ax,
            nodelist=fnodes,
            node_color=fnode_color,
            node_size=fnode_size,
        )
        frag_nodes.append(fnodes)
        frag_pos.append(fpos)
    nx.draw_networkx_labels(G, pos=pos, labels=labels, ax=ax)
    nx.draw_networkx_edges(G, pos=pos, ax=ax)
    plt.axis("off")
    return G, frag_nodes, frag_pos


def circle_arrow(ax, radius, cent_x, cent_y, start_deg, fill_deg, color="black"):
    """Circular arrrow.

    Based on: https://stackoverflow.com/a/38208040
    """
    arc = matplotlib.patches.Arc(
        (cent_x, cent_y),
        radius,
        radius,
        angle=start_deg,
        theta1=0,
        theta2=fill_deg,
        capstyle="round",
        linestyle="-",
        lw=3.5,
        color=color,
        zorder=10,
    )
    ax.add_patch(arc)

    end_x = cent_x + (radius / 2) * np.cos(np.radians(fill_deg + start_deg))
    end_y = cent_y + (radius / 2) * np.sin(np.radians(fill_deg + start_deg))

    # Triangular arrow head
    head = matplotlib.patches.RegularPolygon(
        (end_x, end_y),
        3,
        radius=radius / 9,
        orientation=np.radians(start_deg + fill_deg),
        color=color,
        zorder=10,
    )
    ax.add_patch(head)
    return [arc, head]


def get_render_exc(ax, exc_ens_nm, foscs, ct_numbers, frag_pos):
    nfrags = len(frag_pos)

    def render_exc(ind, thresh=0.1):
        cts = ct_numbers[ind]

        patches = list()
        # Loop over all fragments
        for from_ in range(nfrags):
            # Fragment center coordinates
            from_pos = frag_pos[from_]
            fx, fy = from_pos
            # Loop over all fragments
            for to_ in range(from_, nfrags):
                to_pos = frag_pos[to_]
                ct = cts[from_, to_]
                # Don't draw anything for small contributions
                if ct < thresh:
                    continue
                tx, ty = to_pos
                # Arrow length ... will be 0.0 for local excitations (LE)
                dx = tx - fx
                dy = ty - fy
                width = 1.0 * ct
                # Draw an arrow for charge-transfer (CT) excitations
                if from_ != to_:
                    arrow = ax.arrow(
                        fx,
                        fy,
                        dx,
                        dy,
                        color="k",
                        width=width,
                        length_includes_head=True,
                        zorder=10,
                    )
                    patches.append(arrow)
                # Draw a circular arrow for LE
                else:
                    radius = 5.0 * ct
                    patches.extend(circle_arrow(ax, radius, fx, fy, 0, 300))

                # Annotate previously drawn patch with contribution
                patches.append(
                    ax.annotate(
                        f"{ct:.2%}",
                        (fx + dx / 2.0, fy + dy / 2.0),
                        color="white",
                        path_effects=[
                            matplotlib.patheffects.withStroke(
                                linewidth=4, foreground="black"
                            )
                        ],
                        fontsize=12,
                        ha="center",
                        zorder=15,
                    )
                )
        exc_en_nm = exc_ens_nm[ind]
        fosc = foscs[ind]
        ax.set_title(f"State {ind+1} at {exc_en_nm:.2f} nm with f={fosc:.5f}")
        return patches

    return render_exc


def get_annotate_state(ax, exc_ens_nm, foscs):
    def annotate_state(ind):
        exc_en_nm = exc_ens_nm[ind]
        fosc = foscs[ind]
        text = f"{exc_en_nm:.2f} nm\nf={fosc:.5f}"
        xy = (exc_en_nm, fosc)
        xytext = (exc_en_nm + 15.0, fosc + 0.125)
        annot = ax.annotate(
            text,
            xy,
            xytext=xytext,
            ha="center",
            arrowprops=dict(arrowstyle="->"),
            zorder=10,
            fontsize=12,
            weight="bold",
        )
        return [
            annot,
        ]

    return annotate_state


def ct_number_plot(atoms, coords, exc_ens, foscs, ct_numbers, bond_inds=None):
    if bond_inds is None:
        # Bond indices without interfragment bonds and/or hydrogen bonds
        bond_inds = get_bond_sets(atoms, coords.reshape(-1, 3))

    frags = get_frags_noh(atoms, coords, bond_inds)

    nfrags = len(frags)
    nexcs = len(exc_ens)
    assert ct_numbers.shape == (nexcs, nfrags, nfrags)
    spectrum = spectrum_from_ens_fosc(exc_ens, foscs)

    exc_ens_nm = spectrum.exc_ens_nm
    foscs = spectrum.fosc
    fosc_max = max(foscs.max() * 1.1, 0.1)
    epsilon = spectrum.epsilon
    # The initial idea was to also support picking excited states.
    #
    # Scale down broadened spectrum so the y-values fit the oscillator strengths.
    # Usually, we would plot oscillator strengths and the broadened spectrum on different
    # axis, but this breaks picking.
    epsilon_fosc = epsilon / epsilon.max() * fosc_max

    fig, (ax0, ax1) = plt.subplots(ncols=2, figsize=(12, 6))

    _, _, frag_pos = plot_geometry_graph(
        ax1, atoms, coords.reshape(-1, 3), frags, bond_inds
    )

    # Plot broadened spectrum
    ax0.plot(spectrum.nm, epsilon_fosc)
    # Use bar-plot for oscillator strengths, as stems from stem-plot are not pickable
    ax0.bar(exc_ens_nm, foscs, width=1.2, picker=True)

    exc_ind_max = nexcs - 1

    # Keep track of applied patches in this list so we can later delete them when
    # we cycle through the states.
    cur_patches = list()
    annotate_state = get_annotate_state(ax0, exc_ens_nm, foscs)
    render_exc = get_render_exc(ax1, exc_ens_nm, foscs, ct_numbers, frag_pos)

    def update_plot(ind):
        # Delete patches (arrows, annotations, ...) that are already present
        for patch in cur_patches:
            patch.remove()
        for _ in range(len(cur_patches)):
            cur_patches.pop(0)
        cur_patches.extend(annotate_state(ind))
        cur_patches.extend(render_exc(ind))
        plt.draw()

    # Pointer/index of the currently selected excited state
    exc_ind = 0

    def on_key_release(event):
        nonlocal exc_ind
        key = event.key
        if key == "pagedown":
            exc_ind = max(0, exc_ind - 1)
            update_plot(exc_ind)
        if key == "pageup":
            exc_ind = min(exc_ind_max, exc_ind + 1)
            update_plot(exc_ind)

    okr = fig.canvas.mpl_connect("key_release_event", on_key_release)
    update_plot(exc_ind)
    return fig
