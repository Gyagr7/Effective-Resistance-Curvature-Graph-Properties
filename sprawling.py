"""
sprawling.py
============

Check whether a graph G is "sprawling" in the sense of the paper.

Definition. G = (V, E) is sprawling if there is a collection
S = {H_1, ..., H_q} of Hamiltonian paths of G such that:

    (1) every edge of G lies in at least one H_i;
    (2) for every U with 2 <= |U| <= |V|-1 such that G[U] is connected,
        some H_i restricts to a subgraph of G[U] that is NOT a spanning
        tree of G[U] (equivalently: NOT every vertex of U appears as a
        contiguous interval of H_i, or the restriction has too few edges).

Every sprawling graph is RN (Theorem 8). This module brute-forces over
all Hamiltonian paths of G and all relevant connected vertex subsets U,
so it is only intended for small graphs (roughly n <= 11-12 depending on
density) -- exactly the regime used for the examples in the paper.

Graph representation: a plain dict {vertex: set(neighbors)}, matching
the format used in the exploratory notebooks. Use `from_networkx` to
convert from a networkx.Graph.
"""

from __future__ import annotations

import itertools
from typing import Dict, FrozenSet, List, Set


def from_networkx(G) -> Dict:
    """Convert a networkx.Graph into the adjacency-dict format used here."""
    return {v: set(G.neighbors(v)) for v in G.nodes()}


def _is_connected_subset(G: Dict, vertices: Set) -> bool:
    vertices = set(vertices)
    if not vertices:
        return True
    start = next(iter(vertices))
    visited = set()
    stack = [start]
    while stack:
        v = stack.pop()
        if v not in visited:
            visited.add(v)
            stack.extend((G[v] & vertices) - visited)
    return visited == vertices


def all_hamiltonian_paths(G: Dict) -> List[List]:
    """
    Enumerate all Hamiltonian paths of G by backtracking, up to reversal
    (i.e. a path and its reverse are counted once).
    """
    vertices = list(G.keys())
    n = len(vertices)
    paths = []
    seen = set()

    def extend(path, used):
        if len(path) == n:
            tup, rev = tuple(path), tuple(reversed(path))
            canon = min(tup, rev)
            if canon not in seen:
                seen.add(canon)
                paths.append(list(canon))
            return
        last = path[-1]
        for nbr in G[last]:
            if nbr not in used:
                used.add(nbr)
                path.append(nbr)
                extend(path, used)
                path.pop()
                used.remove(nbr)

    for start in vertices:
        extend([start], {start})

    return paths


def _could_be_spanned_by_a_path(G: Dict, F: Set) -> bool:
    """Necessary condition for G[F] to admit a Hamiltonian path: max degree <= 2."""
    degs = [len(G[v] & F) for v in F]
    return max(degs, default=0) <= 2


def _relevant_connected_subsets(G: Dict) -> List[FrozenSet]:
    """
    All proper, nonempty, connected vertex sets U (2 <= |U| <= n-1) whose
    induced subgraph could plausibly be spanned by a path -- this is the
    set of U's that condition (2) needs to be checked against.
    """
    V = set(G.keys())
    n = len(V)
    out = []
    for r in range(2, n):
        for subset in itertools.combinations(V, r):
            U = set(subset)
            if _is_connected_subset(G, U) and _could_be_spanned_by_a_path(G, U):
                out.append(frozenset(U))
    return out


def _interval_in_path(pos: Dict, U: FrozenSet) -> bool:
    """True iff U occupies a contiguous block of positions in the path."""
    vals = [pos[v] for v in U]
    return max(vals) - min(vals) + 1 == len(vals)


def is_sprawling(G: Dict, verbose: bool = False):
    """
    Decide whether G is sprawling.

    Returns (result, witness) where:
      * result is True/False
      * witness is None if sprawling, or the first uncovered subset U
        (a set of vertices) found if not sprawling / has no Hamiltonian
        path at all.
    """
    H_paths = all_hamiltonian_paths(G)
    if verbose:
        print(f"Number of Hamiltonian paths (up to reversal): {len(H_paths)}")

    if not H_paths:
        return False, "no Hamiltonian path exists"

    # Condition (1): every edge covered by some H_i.
    covered_edges = set()
    for H in H_paths:
        for i in range(len(H) - 1):
            covered_edges.add(frozenset((H[i], H[i + 1])))
    all_edges = {frozenset((u, v)) for u in G for v in G[u]}
    missing_edges = all_edges - covered_edges
    if missing_edges:
        return False, f"edge(s) not covered by any Hamiltonian path: {missing_edges}"

    # Condition (2): every relevant U is "broken" by some H_i.
    Us = _relevant_connected_subsets(G)
    if verbose:
        print(f"Number of relevant connected subsets U: {len(Us)}")

    positions = [{v: i for i, v in enumerate(H)} for H in H_paths]
    uncovered = set(Us)

    for pos in positions:
        if not uncovered:
            break
        broken_now = {U for U in uncovered if not _interval_in_path(pos, U)}
        uncovered -= broken_now

    if uncovered:
        witness = next(iter(uncovered))
        return False, set(witness)

    return True, None


if __name__ == "__main__":
    # Sanity check: the 3-cycle is sprawling (Figure 5a in the paper).
    triangle = {0: {1, 2}, 1: {0, 2}, 2: {0, 1}}
    result, witness = is_sprawling(triangle, verbose=True)
    print(f"Triangle: sprawling={result}  (expected True)")

    # Sanity check: a graph with an edge in every Hamiltonian path is
    # not sprawling (Figure 5c-style counterexample idea: a "bridge"
    # edge shared by every Hamiltonian path).
    bowtie = {0: {1, 2}, 1: {0, 2}, 2: {0, 1, 3}, 3: {2, 4, 5}, 4: {3, 5}, 5: {3, 4}}
    result, witness = is_sprawling(bowtie, verbose=True)
    print(f"Bowtie of two triangles: sprawling={result}  (expected False), witness={witness}")
