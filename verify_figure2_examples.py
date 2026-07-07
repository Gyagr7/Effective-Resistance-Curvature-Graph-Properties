"""
verify_figure2_examples.py
===========================

Runs resistance.py, sprawling.py, and toughness.py against every
non-trivial example graph in Figure 2 of the paper, and prints the
result labeled by panel letter, so a reviewer can check the paper's
claims against the code's output directly -- either by reading the
captured output block at the bottom of this file, or by re-running the
script themselves.

Graphs (a)-(f) are the constructed examples specific to this paper,
reconstructed exactly from the TikZ source of Figure 2. Panels (g), (h),
(i) -- the Petersen graph, the Thomassen 34-graph, and the path family
-- are handled as follows:

  (g) Petersen graph: included below. Fully checked (all four modules);
      matches its known RP status.
  (i) Path family: included below (representative instance P_7). Fully
      checked; matches its SRN status (RN but not RP) via the "tree"
      special case in resistance.py.
  (h) Thomassen 34-graph: DELIBERATELY NOT reconstructed or checked
      here. Its RP status is already established rigorously in the
      paper by Theorem 5, via a hand-constructed rational edge
      weighting (a=33/68, b=99/136, c=165/272, Figure 4) verified
      directly against the closed-form spanning-tree-polytope
      inequalities in the proof -- this is a complete, independent proof
      and does not need (or benefit from) a computational re-check.
      Separately, at 34 vertices it would also be outside the range
      where the brute-force sprawling/toughness checkers here are
      tractable (they enumerate Hamiltonian paths / vertex subsets),
      though that is not the reason for omitting it -- the reason is
      that Theorem 5 already settles the question directly.

STATUS NOTE on panel (e): reconstructing the K4-hub-plus-4-legs graph
exactly from the TikZ source and checking it two independent ways (the
general cutting-plane LP in resistance.py, and an exact brute-force
solve over all 2197 of its actual spanning trees, cross-checked with
both cvxpy and scipy's HiGHS solver) gives RN = False, contradicting the
panel's "SRN" caption. This is flagged below rather than silently
reported as a pass/fail -- see the printed note for panel (e). All other
panels' code output matches their captions exactly.

NOTE on reproducibility: panel (D)'s "not sprawling" witness set is one
of three structurally symmetric leg-vertex-pairs (any of the 3 legs
would serve equally as a witness); which one gets printed can vary
between runs due to Python's set iteration order, but the underlying
result (sprawling=False) is exact and deterministic.
"""

from __future__ import annotations

import networkx as nx

import resistance
import sprawling
import toughness


# ---------------------------------------------------------------------
# Graph constructions, reconstructed from Figure 2's TikZ source
# ---------------------------------------------------------------------

def fig2_a_bowtie() -> nx.Graph:
    """
    Panel (a): two triangles sharing a single hub vertex.
    TikZ: node (0) at origin; for x in {-1,1}, nodes at (-1,x) and (1,x);
    edges (0)--(-1,x)--(1,x)--(0) for each x.
    Caption: "traceable, not 2-connected, not RN"
    """
    G = nx.Graph()
    G.add_edges_from([("0", "A"), ("A", "B"), ("B", "0"),
                       ("0", "C"), ("C", "D"), ("D", "0")])
    return G


def fig2_b_k23() -> nx.Graph:
    """
    Panel (b): K_{2,3}.
    TikZ: nodes a, b; three nodes at x=1; edges a--n_i--b for each.
    Caption: "sprawling, not 1-tough"
    """
    return nx.complete_bipartite_graph(2, 3)


def fig2_c_banana(legs: int = 4) -> nx.Graph:
    """
    Panel (c): two hub vertices a, b joined by 4 internally-disjoint
    paths of length 3 (a generalized theta graph).
    TikZ: nodes a, b; for x in {-1.5,-.5,.5,1.5}, nodes at (1,x),(2,x);
    edges a--(1,x)--(2,x)--b for each x.
    Caption: "2-connected, not 1-tough, not traceable, not RN"
    """
    G = nx.Graph()
    G.add_node("a")
    G.add_node("b")
    for i in range(legs):
        x, y = f"x{i}", f"y{i}"
        G.add_edges_from([("a", x), (x, y), (y, "b")])
    return G


def fig2_clique_plus_legs(k: int, leg_len: int = 2) -> nx.Graph:
    """
    Panels (d)/(e)/(f): a K_k hub clique, with each of the k hub
    vertices connected via a path of `leg_len` intermediate vertices to
    a single common external vertex v.
    TikZ pattern (shown for K4, panel (e)): node (0) [= v]; for x in
    hub-index-set, nodes at (0,x) [hub], (1,..) [a_x], (2,..) [b_x];
    edges hub_x--a_x--b_x--(0); plus explicit \\draw commands forming
    K_k among the hub vertices.
      k=3: panel (d), caption "traceable, RP, not sprawling"
      k=4: panel (e), caption "traceable, 1-tough, SRN, not sprawling"
           -- SEE STATUS NOTE ABOVE; code currently finds RN=False,
           contradicting this caption. Flagged, not silently passed.
      k=5: panel (f), caption "1-tough, not RN, not traceable"
    """
    G = nx.Graph()
    hubs = [f"h{i}" for i in range(k)]
    G.add_edges_from([(hubs[i], hubs[j]) for i in range(k) for j in range(i + 1, k)])
    G.add_node("v")
    for i, h in enumerate(hubs):
        prev = h
        for s in range(leg_len):
            node = f"p{i}_{s}"
            G.add_edge(prev, node)
            prev = node
        G.add_edge(prev, "v")
    return G


def fig2_g_petersen() -> nx.Graph:
    """Panel (g): the Petersen graph. Caption: "Petersen graph" (RP but
    not Hamiltonian; classical, see Section 1 of the paper)."""
    return nx.petersen_graph()


def fig2_i_path(n: int = 7) -> nx.Graph:
    """
    Panel (i): a path on n >= 3 vertices (n=7 used here as a
    representative instance -- the paper's claim holds for every path
    with n >= 3, and resistance.py's "tree" special case decides any
    instance in O(1) via a direct degree check, so n=7 exercises the
    same code path as any other n).
    Caption: "Path with >= 3 vertices: SRN, traceable, not 2-connected"
    """
    return nx.path_graph(n)


# ---------------------------------------------------------------------
# Verification runner
# ---------------------------------------------------------------------

def _traceable(G: nx.Graph) -> bool:
    return len(sprawling.all_hamiltonian_paths(sprawling.from_networkx(G))) > 0


def verify(label: str, G: nx.Graph, caption: str, check_sprawling: bool = False):
    print(f"{label}  (n={G.number_of_nodes()}, m={G.number_of_edges()})")
    print(f"  Figure 2 caption: \"{caption}\"")

    trace = _traceable(G)
    two_conn = nx.is_biconnected(G)
    is_rp, is_rn, cert, _ = resistance.resistance_positive_decision(G, verbose=False)
    tough, _ = toughness.is_one_tough(G)

    print(f"  traceable    = {trace}")
    print(f"  2-connected  = {two_conn}")
    print(f"  RN           = {is_rn}   [method: {cert['method']}]")
    print(f"  RP           = {is_rp}   [method: {cert['method']}]")
    print(f"  1-tough      = {tough}")

    if check_sprawling:
        sprawl, info = sprawling.is_sprawling(sprawling.from_networkx(G))
        if sprawl:
            print(f"  sprawling    = True   (explicit witness S found, |S|={len(info)}, "
                  f"independently re-verified via verify_sprawling_set)")
        else:
            print(f"  sprawling    = False  (reason: {info})")

    print()


if __name__ == "__main__":
    print("=" * 78)
    print("Figure 2 example verification")
    print("=" * 78)
    print()

    verify(
        "(A) Bowtie (two triangles sharing a hub vertex)",
        fig2_a_bowtie(),
        "traceable, not 2-connected, not RN",
    )

    verify(
        "(B) K_{2,3}",
        fig2_b_k23(),
        "sprawling, not 1-tough",
        check_sprawling=True,
    )

    verify(
        "(C) Banana (4 parallel length-3 paths between two hubs)",
        fig2_c_banana(4),
        "2-connected, not 1-tough, not traceable, not RN",
    )

    verify(
        "(D) K3 hub + 3 legs to a common point",
        fig2_clique_plus_legs(3, 2),
        "traceable, RP, not sprawling",
        check_sprawling=True,
    )

    print("(E) K4 hub + 4 legs to a common point")
    print("  Figure 2 caption: \"traceable, 1-tough, SRN, not sprawling\"")
    print("  *** STATUS: UNRESOLVED -- see module docstring above. ***")
    print("  Reconstructing this graph from the TikZ source and checking it two")
    print("  independent ways (the general LP in resistance.py, and an EXACT")
    print("  brute-force solve over all 2197 of its actual spanning trees, cross-")
    print("  checked with cvxpy and scipy's HiGHS solver) both give RN = False,")
    print("  contradicting the \"SRN\" (RN=True) caption. Not reporting a pass/fail")
    print("  here pending confirmation of the intended graph structure.")
    print()

    verify(
        "(F) K5 hub + 5 legs to a common point",
        fig2_clique_plus_legs(5, 2),
        "1-tough, not RN, not traceable",
    )

    verify(
        "(G) Petersen graph",
        fig2_g_petersen(),
        "Petersen graph (RP but not Hamiltonian)",
        check_sprawling=True,
    )

    print("(H) Thomassen 34-graph")
    print("  Not reconstructed or checked computationally here. Its RP status")
    print("  is proven directly in the paper (Theorem 5) via a hand-constructed")
    print("  rational edge weighting (a=33/68, b=99/136, c=165/272, Figure 4),")
    print("  independent of this codebase. See module docstring for details.")
    print()

    verify(
        "(I) Path P_7 (representative of the path family, n >= 3)",
        fig2_i_path(7),
        "Path with >= 3 vertices: SRN, traceable, not 2-connected",
    )

    print("=" * 78)
    print("Summary: panels (A), (B), (C), (D), (F), (G), (I) all match their")
    print("captions exactly. Panel (E) is flagged as unresolved -- see notes")
    print("above. Panel (H) is intentionally out of scope -- see notes above.")
    print("=" * 78)





# =============================================================================
# CAPTURED OUTPUT (from running this script; regenerate with:
#     python verify_figure2_examples.py
# to re-verify. cvxpy 1.9.2 with SCS solver, networkx >= 3.0.)
# =============================================================================
r"""
==============================================================================
Figure 2 example verification
==============================================================================

(A) Bowtie (two triangles sharing a hub vertex)  (n=5, m=6)
  Figure 2 caption: "traceable, not 2-connected, not RN"
  traceable    = True
  2-connected  = False
  RN           = False   [method: not-2-connected-not-path]
  RP           = False   [method: not-2-connected-not-path]
  1-tough      = False

(B) K_{2,3}  (n=5, m=6)
  Figure 2 caption: "sprawling, not 1-tough"
  traceable    = True
  2-connected  = True
  RN           = True   [method: margin-lp]
  RP           = False   [method: margin-lp]
  1-tough      = False
  sprawling    = True   (explicit witness S found, |S|=3, independently re-verified via verify_sprawling_set)

(C) Banana (4 parallel length-3 paths between two hubs)  (n=10, m=12)
  Figure 2 caption: "2-connected, not 1-tough, not traceable, not RN"
  traceable    = False
  2-connected  = True
  RN           = False   [method: margin-lp]
  RP           = False   [method: margin-lp]
  1-tough      = False

(D) K3 hub + 3 legs to a common point  (n=10, m=12)
  Figure 2 caption: "traceable, RP, not sprawling"
  traceable    = True
  2-connected  = True
  RN           = True   [method: margin-lp]
  RP           = True   [method: margin-lp]
  1-tough      = True
  sprawling    = False  (reason: condition (2) fails: U={'p2_1', 'p2_0'} is a spanning tree of G[U] under every H_i in S)

(E) K4 hub + 4 legs to a common point
  Figure 2 caption: "traceable, 1-tough, SRN, not sprawling"
  *** STATUS: UNRESOLVED -- see module docstring above. ***
  Reconstructing this graph from the TikZ source and checking it two
  independent ways (the general LP in resistance.py, and an EXACT
  brute-force solve over all 2197 of its actual spanning trees, cross-
  checked with cvxpy and scipy's HiGHS solver) both give RN = False,
  contradicting the "SRN" (RN=True) caption. Not reporting a pass/fail
  here pending confirmation of the intended graph structure.

(F) K5 hub + 5 legs to a common point  (n=16, m=25)
  Figure 2 caption: "1-tough, not RN, not traceable"
  traceable    = False
  2-connected  = True
  RN           = False   [method: margin-lp]
  RP           = False   [method: margin-lp]
  1-tough      = True

(G) Petersen graph  (n=10, m=15)
  Figure 2 caption: "Petersen graph (RP but not Hamiltonian)"
  traceable    = True
  2-connected  = True
  RN           = True   [method: margin-lp]
  RP           = True   [method: margin-lp]
  1-tough      = True
  sprawling    = True   (explicit witness S found, |S|=3, independently re-verified via verify_sprawling_set)

(H) Thomassen 34-graph
  Not reconstructed or checked computationally here. Its RP status
  is proven directly in the paper (Theorem 5) via a hand-constructed
  rational edge weighting (a=33/68, b=99/136, c=165/272, Figure 4),
  independent of this codebase. See module docstring for details.

(I) Path P_7 (representative of the path family, n >= 3)  (n=7, m=6)
  Figure 2 caption: "Path with >= 3 vertices: SRN, traceable, not 2-connected"
  traceable    = True
  2-connected  = False
  RN           = True   [method: tree]
  RP           = False   [method: tree]
  1-tough      = False

==============================================================================
Summary: panels (A), (B), (C), (D), (F), (G), (I) all match their
captions exactly. Panel (E) is flagged as unresolved -- see notes
above. Panel (H) is intentionally out of scope -- see notes above.
==============================================================================
"""
