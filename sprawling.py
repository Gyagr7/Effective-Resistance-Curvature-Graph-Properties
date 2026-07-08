"""
sprawling.py
============

Check whether a graph G is "sprawling" in the sense of the paper.

Definition. G = (V, E) is sprawling if there is a collection
S = {H_1, ..., H_q} of Hamiltonian paths of G such that:

    (1) every edge of G lies in at least one H_i;
    (2) for every U with 2 <= |U| <= |V|-1 such that G[U] is connected,
        there exists H_i in S such that H_i[U] is NOT a spanning tree
        of G[U].

`verify_sprawling_set(G, S)` checks these two conditions literally
against a concrete, explicit S -- it does not reason about "all
Hamiltonian paths" or rely on any equivalence. `is_sprawling(G)` then
constructs an explicit witness S (a genuine list of Hamiltonian paths)
and calls verify_sprawling_set on it before returning True, so the
result always comes with a checkable certificate rather than an
implicit argument.

Internally, is_sprawling first checks the FULL pool of all Hamiltonian
paths, since conditions (1) and (2) are both monotone under adding
paths to S (more paths can only cover more edges / break more U's).
So if the full pool already fails, no smaller S can succeed either --
this only reduces work, and every claimed witness S is independently
re-verified via verify_sprawling_set before being returned.

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


def _relevant_connected_subsets(G: Dict) -> List[FrozenSet]:
    """
    All proper, nonempty, connected vertex sets U with 2 <= |U| <= n-1.

    Note: condition (2) of the sprawling definition places NO restriction
    on the degree sequence of G[U] -- it only requires that G[U] be
    connected. (An earlier version of this function incorrectly also
    required max degree <= 2 in G[U], reasoning that only "path-like"
    induced subgraphs could ever have H_i[U] be a spanning tree. That's
    wrong: G[U] can have high-degree vertices while some Hamiltonian
    path H_i still happens to traverse U contiguously, or conversely
    every H_i might fail to do so even for low-degree U. The degree
    sequence of G[U] and the interval structure of U within each H_i are
    logically independent, so every connected U must be checked.)
    """
    V = set(G.keys())
    n = len(V)
    out = []
    for r in range(2, n):
        for subset in itertools.combinations(V, r):
            U = set(subset)
            if _is_connected_subset(G, U):
                out.append(frozenset(U))
    return out


def _interval_in_path(pos: Dict, U: FrozenSet) -> bool:
    """True iff U occupies a contiguous block of positions in the path."""
    vals = [pos[v] for v in U]
    return max(vals) - min(vals) + 1 == len(vals)


def verify_sprawling_set(G: Dict, S: List[List]):
    """
    Directly check conditions (1) and (2) of the sprawling definition,
    verbatim, against a concrete, explicit collection S of Hamiltonian
    paths of G. This function does NOT reason about "all Hamiltonian
    paths" or any equivalence -- it takes S exactly as given and checks
    exactly what the two bullet points say.

    Returns (ok, reason):
      * ok = True, reason = None                     if S satisfies both conditions
      * ok = False, reason = "condition (1) fails..." if some edge of G
        is in no H_i in S
      * ok = False, reason = "condition (2) fails..." if some qualifying
        U is a spanning tree of G[U] under every single H_i in S
    """
    all_edges = {frozenset((u, v)) for u in G for v in G[u]}

    # Condition (1), verbatim: every edge of G is contained in at least
    # one path in S.
    covered_edges = set()
    for H in S:
        for i in range(len(H) - 1):
            covered_edges.add(frozenset((H[i], H[i + 1])))
    missing_edges = all_edges - covered_edges
    if missing_edges:
        return False, f"condition (1) fails: edge(s) {set(missing_edges)} in no H_i in S"

    # Condition (2), verbatim: for every U with 2 <= |U| <= |V|-1 and
    # G[U] connected, there exists H_i in S with H_i[U] not a spanning
    # tree of G[U]. H_i[U] is always a linear forest (a subgraph of a
    # path), so it is a spanning tree of G[U] iff (a) it has exactly
    # |U|-1 edges AND (b) it is connected -- equivalently, iff U occupies
    # a contiguous block of positions in H_i. We check this directly,
    # per U and per H_i, without assuming the equivalence in advance.
    positions = [{v: i for i, v in enumerate(H)} for H in S]
    Us = _relevant_connected_subsets(G)

    for U in Us:
        broken = False
        for H, pos in zip(S, positions):
            edges_in_U = sum(
                1 for i in range(len(H) - 1)
                if H[i] in U and H[i + 1] in U
            )
            vals = sorted(pos[v] for v in U)
            contiguous = (vals[-1] - vals[0] + 1 == len(vals))
            is_spanning_tree = contiguous and edges_in_U == len(U) - 1
            if not is_spanning_tree:
                broken = True
                break
        if not broken:
            return False, f"condition (2) fails: U={set(U)} is a spanning tree of G[U] under every H_i in S"

    return True, None


def is_sprawling(G: Dict, verbose: bool = False):
    """
    Decide whether G is sprawling, and if so, exhibit an explicit
    witness collection S of Hamiltonian paths satisfying conditions (1)
    and (2) verbatim (checked via verify_sprawling_set, not inferred).

    Strategy: the candidate pool is all Hamiltonian paths of G. Both
    conditions are monotone under adding paths to S (more paths can only
    cover more edges / break more U's, never fewer), so:
      * if the FULL pool already fails either condition, no subset can
        succeed, and G is not sprawling;
      * otherwise, we greedily select a small explicit subset S from the
        pool that still satisfies both conditions, and verify S directly
        and literally with verify_sprawling_set before returning it.

    Returns (result, info) where:
      * result is True/False
      * if True,  info is the explicit witness S (list of Hamiltonian paths)
      * if False, info is the reason (from verify_sprawling_set, or "no
        Hamiltonian path exists")
    """
    H_paths = all_hamiltonian_paths(G)
    if verbose:
        print(f"Number of Hamiltonian paths (up to reversal): {len(H_paths)}")

    if not H_paths:
        return False, "no Hamiltonian path exists"

    # Feasibility check against the full pool -- if even every
    # Hamiltonian path together can't satisfy (1) and (2), nothing can.
    full_ok, full_reason = verify_sprawling_set(G, H_paths)
    if not full_ok:
        return False, full_reason

    # Greedily build a smaller explicit S from H_paths.
    all_edges = {frozenset((u, v)) for u in G for v in G[u]}
    Us = _relevant_connected_subsets(G)
    positions_all = [{v: i for i, v in enumerate(H)} for H in H_paths]

    def edges_of(H):
        return {frozenset((H[i], H[i + 1])) for i in range(len(H) - 1)}

    def breaks(pos, U):
        vals = sorted(pos[v] for v in U)
        return vals[-1] - vals[0] + 1 != len(vals)

    remaining_edges = set(all_edges)
    remaining_Us = set(Us)
    S: List[List] = []

    while remaining_edges or remaining_Us:
        best_gain, best_H, best_pos = -1, None, None
        for H, pos in zip(H_paths, positions_all):
            gain = len(edges_of(H) & remaining_edges)
            gain += sum(1 for U in remaining_Us if breaks(pos, U))
            if gain > best_gain:
                best_gain, best_H, best_pos = gain, H, pos
        if best_gain <= 0:
            break  # should not happen since the full pool is feasible
        S.append(best_H)
        remaining_edges -= edges_of(best_H)
        remaining_Us = {U for U in remaining_Us if not breaks(best_pos, U)}

    # Final, literal verification of the constructed S before returning
    # it as a witness. Never trust the greedy construction on its own.
    ok, reason = verify_sprawling_set(G, S)
    if not ok:
        # Fall back to the full pool, which we already verified above.
        S, ok, reason = H_paths, True, None

    if verbose:
        print(f"Witness |S| = {len(S)} (out of {len(H_paths)} total Hamiltonian paths)")

    return True, S


if __name__ == "__main__":
    # Sanity check: the 3-cycle is sprawling (Figure 5a in the paper).
    triangle = {0: {1, 2}, 1: {0, 2}, 2: {0, 1}}
    result, info = is_sprawling(triangle, verbose=True)
    print(f"Triangle: sprawling={result}  (expected True), |S|={len(info) if result else '-'}")
    if result:
        ok, reason = verify_sprawling_set(triangle, info)
        print(f"  re-verified witness S directly: ok={ok}")

    # Sanity check: a bowtie of two triangles sharing a bridge edge is
    # not sprawling -- the bridge edge is in every Hamiltonian path, so
    # some U (the two vertices spanning either triangle side) can never
    # be broken.
    bowtie = {0: {1, 2}, 1: {0, 2}, 2: {0, 1, 3}, 3: {2, 4, 5}, 4: {3, 5}, 5: {3, 4}}
    result, info = is_sprawling(bowtie, verbose=True)
    print(f"Bowtie of two triangles: sprawling={result}  (expected False), reason={info}")
