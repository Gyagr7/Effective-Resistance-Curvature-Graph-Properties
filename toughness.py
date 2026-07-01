"""
toughness.py
============

Compute the (vertex) toughness of a graph, and check 1-toughness.

Definition. G is t-tough if for every vertex set S whose removal
disconnects G, the number of components of G - S is at most |S| / t.
Equivalently, the toughness of G is

    tau(G) = min_{S} |S| / c(G - S)

over all S subseteq V(G) such that G - S is disconnected, and G is
1-tough iff tau(G) >= 1, i.e. c(G - S) <= |S| for every disconnecting S.

This is brute-force (checks all vertex subsets), so it is only intended
for small graphs, matching the examples used in the paper (e.g. the
Dawes-Rodrigues-style construction G_5(s_1, ..., s_5) used in the proof
of Theorem 4).
"""

from __future__ import annotations

from itertools import combinations

import networkx as nx


def toughness(G: nx.Graph, verbose: bool = False):
    """
    Compute tau(G) exactly by brute force over all vertex subsets S.

    Returns (tau, witness_S) where witness_S is a minimizing cut set
    (None if G has no disconnecting set at all, e.g. G is complete --
    in which case tau(G) is conventionally infinite).
    """
    nodes = list(G.nodes())
    n = len(nodes)

    best_tau = float("inf")
    best_S = None

    # S can range over all proper subsets; only sets whose removal
    # disconnects G are relevant.
    for r in range(1, n):
        for S in combinations(nodes, r):
            S_set = set(S)
            H = G.copy()
            H.remove_nodes_from(S_set)
            if H.number_of_nodes() == 0:
                continue
            c = nx.number_connected_components(H)
            if c <= 1:
                continue  # S does not disconnect G
            ratio = len(S_set) / c
            if ratio < best_tau:
                best_tau = ratio
                best_S = S_set
                if verbose:
                    print(f"new min: |S|={len(S_set)}, components={c}, ratio={ratio:.4f}, S={S_set}")

    return best_tau, best_S


def is_one_tough(G: nx.Graph, verbose: bool = False):
    """
    Decide whether G is 1-tough, i.e. c(G - S) <= |S| for every S whose
    removal disconnects G.

    Returns (result, witness) where witness is the first offending cut
    set found if G is not 1-tough, else None.
    """
    nodes = list(G.nodes())
    n = len(nodes)

    for r in range(1, n):
        for S in combinations(nodes, r):
            S_set = set(S)
            H = G.copy()
            H.remove_nodes_from(S_set)
            if H.number_of_nodes() == 0:
                continue
            c = nx.number_connected_components(H)
            if c > len(S_set):
                if verbose:
                    print(f"Not 1-tough: |S|={len(S_set)}, components(G-S)={c}, S={S_set}")
                return False, S_set

    if verbose:
        print("Graph is 1-tough")
    return True, None


if __name__ == "__main__":
    # Sanity check: K_4 is complete, hence trivially 1-tough (no
    # disconnecting set exists at all).
    K4 = nx.complete_graph(4)
    result, witness = is_one_tough(K4, verbose=True)
    print(f"K4: 1-tough={result}  (expected True)\n")

    # Sanity check: the star graph K_{1,3} is not 1-tough -- removing
    # the center disconnects it into 3 components with |S| = 1.
    star = nx.star_graph(3)  # center 0, leaves 1,2,3
    result, witness = is_one_tough(star, verbose=True)
    print(f"Star K_1,3: 1-tough={result}  (expected False)\n")

    # Sanity check against Lemma 15 in the paper: G_5(1,1,1,1,1), the
    # graph obtained from K_6 by subdividing each spoke v0-vi once,
    # should be 1-tough.
    G = nx.complete_graph(6)  # vertices 0..5, v0 = 0, spokes to 1..5
    G = nx.relabel_nodes(G, {i: f"v{i}" for i in range(6)})
    G.remove_edges_from([("v0", f"v{i}") for i in range(1, 6)])
    for i in range(1, 6):
        G.add_node(f"x{i}")
        G.add_edge("v0", f"x{i}")
        G.add_edge(f"x{i}", f"v{i}")
    tau, S = toughness(G, verbose=False)
    result, witness = is_one_tough(G)
    print(f"G_5(1,1,1,1,1): tau={tau:.4f}, 1-tough={result}  (expected True, tau>=1)")
