"""
examples.py
===========

Graph constructions referenced in the paper, for use with resistance.py,
sprawling.py, and toughness.py.
"""

from __future__ import annotations

from itertools import combinations

import networkx as nx


def petersen_graph() -> nx.Graph:
    """The Petersen graph: RP but not Hamiltonian (Section 1)."""
    return nx.petersen_graph()


def grid_graph(m: int, n: int) -> nx.Graph:
    """
    The grid graph P_m x P_n (Cartesian product of two paths). Sprawling
    for all m, n >= 2 (Theorem 17), hence RN (Theorem 8).
    """
    return nx.grid_2d_graph(m, n)


def toughness_family(s: list[int]) -> nx.Graph:
    """
    The graph G_t(s_1, ..., s_t) from the proof of Theorem 4 / Lemma 15:
    take K_{t+1} on {v0, v1, ..., vt} and subdivide edge v0-vi exactly
    s_i times. For t = 5 with all s_i >= 1, this is 1-tough but not RN
    (Theorem 4), disproving Fiedler's conjecture that every 1-tough
    graph is RP.

    `s` should have length t (any t >= 3 is meaningful; the paper uses
    t = 5).
    """
    t = len(s)
    G = nx.complete_graph(t + 1)
    G = nx.relabel_nodes(G, {i: f"v{i}" for i in range(t + 1)})

    for i in range(1, t + 1):
        G.remove_edge("v0", f"v{i}")
        prev = "v0"
        for k in range(s[i - 1]):
            node = f"x{i}_{k}"
            G.add_edge(prev, node)
            prev = node
        G.add_edge(prev, f"v{i}")

    return G


def build_minimal_tough_graph(n: int, l: int) -> nx.Graph:
    """
    Advisor's original constructor for the Theorem 4 / Lemma 15 family:
    a hub vertex 0 connected via n branches (each a path of l edges) to
    n "end" vertices, which are then pairwise connected into a clique.

    This is the SAME construction as `toughness_family` above (a hub
    connected to a clique via subdivided spokes, with the clique itself
    left intact) -- just built with integer-labeled branch vertices
    instead of string-labeled ones, and with all branches forced to
    equal length l instead of allowing per-branch lengths s_1,...,s_t.
    Kept here (rather than only using toughness_family) so any code
    referencing this name by the advisor's original convention still
    works; the paper uses the t=5, s_i>=1 case in the proof of Theorem 4.
    """
    interval = 50
    edges = []
    end = []
    for i in range(1, n + 1):
        start = interval * i
        edges += [(0, start)]
        for j in range(l):
            edges += [(start + j, start + j + 1)]
        end += [start + l]
    edges += list(combinations(end, 2))
    G = nx.Graph()
    G.add_edges_from(edges)
    return G


def bowtie() -> nx.Graph:
    """Two triangles sharing a bridge edge -- traceable but not sprawling."""
    G = nx.Graph()
    G.add_edges_from([
        (0, 1), (1, 2), (2, 0),
        (2, 3),
        (3, 4), (4, 5), (5, 3),
    ])
    return G


if __name__ == "__main__":
    import resistance
    import sprawling
    import toughness

    print("=== Petersen graph ===")
    G = petersen_graph()
    rp, rn, cert, _ = resistance.resistance_positive_decision(G, verbose=False)
    print(f"RN={rn}, RP={rp}  [{cert['method']}]  (expect RP=True)\n")

    print("=== Grid graph P_3 x P_3 ===")
    G = grid_graph(3, 3)
    sprawl, info = sprawling.is_sprawling(sprawling.from_networkx(G))
    print(f"sprawling={sprawl}  (expect True), |S|={len(info) if sprawl else '-'}\n")

    print("=== G_5(1,1,1,1,1) [Theorem 4 construction] ===")
    G = toughness_family([1, 1, 1, 1, 1])
    tough, S = toughness.is_one_tough(G)
    rp, rn, cert, _ = resistance.resistance_positive_decision(nx.convert_node_labels_to_integers(G), verbose=False)
    print(f"1-tough={tough} (expect True), RN={rn} (expect False)\n")

    print("=== Bowtie (two triangles sharing a bridge edge) ===")
    G = bowtie()
    sprawl, reason = sprawling.is_sprawling(sprawling.from_networkx(G))
    print(f"sprawling={sprawl}  (expect False), reason={reason}")
