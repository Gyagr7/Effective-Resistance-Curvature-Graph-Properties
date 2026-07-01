"""
examples.py
===========

Graph constructions referenced in the paper, for use with resistance.py,
sprawling.py, and toughness.py.
"""

from __future__ import annotations

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
    rp, rn, t_star, _ = resistance.resistance_positive_decision(G)
    print(f"RN={rn}, RP={rp}, t*={t_star:.4f}  (expect RP=True)\n")

    print("=== Grid graph P_3 x P_3 ===")
    G = grid_graph(3, 3)
    sprawl, witness = sprawling.is_sprawling(sprawling.from_networkx(G))
    print(f"sprawling={sprawl}  (expect True)\n")

    print("=== G_5(1,1,1,1,1) [Theorem 4 construction] ===")
    G = toughness_family([1, 1, 1, 1, 1])
    tough, S = toughness.is_one_tough(G)
    rp, rn, t_star, _ = resistance.resistance_positive_decision(nx.convert_node_labels_to_integers(G))
    print(f"1-tough={tough} (expect True), RN={rn} (expect False)\n")

    print("=== Bowtie (two triangles) ===")
    G = bowtie()
    sprawl, witness = sprawling.is_sprawling(sprawling.from_networkx(G))
    print(f"sprawling={sprawl}  (expect False), witness={witness}")
