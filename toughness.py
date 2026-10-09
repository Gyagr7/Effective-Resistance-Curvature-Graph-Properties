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

This enumerates every vertex subset, so it is only usable on small
graphs, roughly n <= 12. That covers the G_t(s_1, ..., s_t) family of
the appendix, which `examples.py` builds.
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


def _g5() -> nx.Graph:
    """G_5(1,1,1,1,1): K_6 with each spoke v0-vi subdivided once."""
    G = nx.complete_graph(6)
    G = nx.relabel_nodes(G, {i: f"v{i}" for i in range(6)})
    G.remove_edges_from([("v0", f"v{i}") for i in range(1, 6)])
    for i in range(1, 6):
        G.add_edge("v0", f"x{i}")
        G.add_edge(f"x{i}", f"v{i}")
    return G


# Each check is a graph, its toughness, and whether it is 1-tough.
# K_4 is complete, so no set disconnects it and its toughness is
# infinite by convention. Removing the centre of K_{1,3} leaves three
# components, so its toughness is 1/3.
CHECKS = (
    ("K_4", nx.complete_graph(4), float("inf"), True),
    ("K_{1,3}", nx.star_graph(3), 1 / 3, False),
    ("C_5", nx.cycle_graph(5), 1.0, True),
    ("K_{2,3}", nx.complete_bipartite_graph(2, 3), 2 / 3, False),
    ("G_5(1,1,1,1,1)", _g5(), 1.0, True),
)


def main() -> int:
    """Run the checks, and return 1 if any graph disagrees with its row."""
    problems = []
    print(f"{'graph':<16} {'tau':>8}  1-tough")
    for name, G, want_tau, want_tough in CHECKS:
        tau, _ = toughness(G)
        tough, _ = is_one_tough(G)
        print(f"{name:<16} {tau:>8.4f}  {'y' if tough else 'n'}")
        if tau != want_tau and abs(tau - want_tau) > 1e-9:
            problems.append(f"{name}: tau is {tau}, expected {want_tau}")
        if tough is not want_tough:
            problems.append(f"{name}: 1-tough is {_flag(tough)}, "
                            f"expected {_flag(want_tough)}")

    print()
    if problems:
        print(f"{len(problems)} disagreement(s):")
        for line in problems:
            print(f"  {line}")
        return 1
    print(f"all {len(CHECKS)} graphs agree with their known toughness")
    return 0


def _flag(value: bool) -> str:
    return "y" if value else "n"


if __name__ == "__main__":
    raise SystemExit(main())
