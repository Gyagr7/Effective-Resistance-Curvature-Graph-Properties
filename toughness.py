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

from fractions import Fraction
from itertools import combinations
from typing import Iterator, Optional, Set, Tuple, Union

import networkx as nx

# tau(G) is an exact Fraction, except on a graph no vertex set
# disconnects, where it is infinite by convention.
Rational = Union[Fraction, float]


def _disconnecting_sets(G: nx.Graph) -> Iterator[Tuple[Set, int]]:
    """
    Yield every proper vertex set S whose removal disconnects G, paired
    with the number of components it leaves.

    G - S is taken as a subgraph view rather than a copy, so the only
    cost per set is the component count itself.
    """
    nodes = list(G.nodes())
    V = set(nodes)
    for r in range(1, len(nodes)):
        for S in combinations(nodes, r):
            c = nx.number_connected_components(G.subgraph(V.difference(S)))
            if c > 1:
                yield set(S), c


def toughness(G: nx.Graph, verbose: bool = False) -> Tuple[Rational, Optional[Set]]:
    """
    Compute tau(G) exactly by brute force over all vertex subsets S.

    Returns (tau, witness_S), with tau an exact Fraction and witness_S a
    minimizing cut set. A graph with no disconnecting set at all, such as
    a complete graph, has no witness and the conventional tau of
    infinity, which comes back as the float rather than as a Fraction.
    """
    best_tau: Rational = float("inf")
    best_S = None

    for S, c in _disconnecting_sets(G):
        ratio = Fraction(len(S), c)
        if ratio < best_tau:
            best_tau, best_S = ratio, S
            if verbose:
                print(f"new min: |S|={len(S)}, components={c}, "
                      f"ratio={ratio}, S={S}")

    return best_tau, best_S


def is_one_tough(G: nx.Graph, verbose: bool = False):
    """
    Decide whether G is 1-tough, i.e. c(G - S) <= |S| for every S whose
    removal disconnects G.

    Returns (result, witness) where witness is the first offending cut
    set found if G is not 1-tough, else None.
    """
    for S, c in _disconnecting_sets(G):
        if c > len(S):
            if verbose:
                print(f"Not 1-tough: |S|={len(S)}, components(G-S)={c}, S={S}")
            return False, S

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
    ("K_{1,3}", nx.star_graph(3), Fraction(1, 3), False),
    ("C_5", nx.cycle_graph(5), Fraction(1), True),
    ("K_{2,3}", nx.complete_bipartite_graph(2, 3), Fraction(2, 3), False),
    ("G_5(1,1,1,1,1)", _g5(), Fraction(1), True),
)


def main() -> int:
    """Run the checks, and return 1 if any graph disagrees with its row."""
    problems = []
    print(f"{'graph':<16} {'tau':>8}  1-tough")
    for name, G, want_tau, want_tough in CHECKS:
        tau, _ = toughness(G)
        tough, _ = is_one_tough(G)
        print(f"{name:<16} {str(tau):>8}  {'y' if tough else 'n'}")
        if tau != want_tau:
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
