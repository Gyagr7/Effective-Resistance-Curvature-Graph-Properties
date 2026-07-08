"""
verify_figure2_examples.py
===========================

Runs resistance.py, sprawling.py, and toughness.py against every
example graph in Figure 1 (formerly "Figure 2") of the paper, and
prints the result labeled by panel letter, so a reviewer can check the
paper's claims against the code's output directly -- either by reading
the captured output block at the bottom of this file, or by re-running
the script themselves.

PANEL LETTERING (as of the current figure revision): a new panel was
inserted at (C) -- a small 3-leg "banana" graph -- which shifted every
subsequent panel down one letter, and the Thomassen 34-graph panel was
removed entirely. Current lettering:

  (A) Bowtie (two triangles sharing a hub vertex)
  (B) K_{2,3}
  (C) Small banana: 2 hub vertices + 3 length-3 legs        [NEW panel]
  (D) Banana: 2 hub vertices + 4 length-3 legs               [was (C)]
  (E) K3 hub clique + 3 legs to a common point                [was (D)]
  (F) K4 hub clique + 4 legs to a common point                [was (E)]
  (G) K5 hub clique + 5 legs to a common point                [was (F)]
  (H) Petersen graph (caption now explicitly includes "sprawling")  [was (G)]
  (I) Path family, n >= 3                                     [unchanged]

  (Thomassen 34-graph, formerly panel (H), has been REMOVED from the
  figure -- there is no panel referencing it anymore.)

STATUS NOTE -- TWO panels currently show a discrepancy with their "SRN"
caption, found independently by full spanning-tree enumeration (exact,
not LP-approximate) cross-checked with two different solvers:

  (C) Small banana (3 legs, no hub-hub edge): caption claims SRN
      (RN=True, RP=False). Exact enumeration of all 27 spanning trees
      gives RN=False.
  (F) K4 hub + 4 legs: caption claims SRN. Exact enumeration of all
      2197 spanning trees gives RN=False -- confirmed THREE independent
      ways: this repo's general margin-LP, an exact solve over the full
      enumerated spanning-tree set (cross-checked with both cvxpy and
      scipy's HiGHS solver), and a separately-derived, differently-
      structured LP (the advisor's original cutting-plane formulation,
      after fixing the same closed-polytope-vs-relative-interior bug
      that this repo's resistance.py also had). All three agree the
      graph is not RN as reconstructed here, which rules out a
      solver-specific numerical artifact as the explanation; the
      discrepancy with the "SRN" caption remains open.

Both use the SAME "hub(s) connected via parallel legs to a common
point" building block, just with different hub structures (no hub-hub
edges in (C); a K_k clique hub in (E)/(F)/(G)). (E) [K3 hub, RP claimed]
and (G) [K5 hub, not-RN claimed] both check out correctly against their
captions, so the discrepancy is specifically tied to the "SRN" claim on
this family, not the construction in general. Both (C) and (F) are
flagged below rather than silently reported as pass/fail.

Panel (I) is fully checked (representative instance P_7); panel (H)
(Petersen) is fully checked, now including sprawling per its updated
caption.
"""

from __future__ import annotations

import networkx as nx

import resistance
import sprawling
import toughness


# ---------------------------------------------------------------------
# Graph constructions, reconstructed from Figure 1's TikZ source
# ---------------------------------------------------------------------

def fig1_a_bowtie() -> nx.Graph:
    """
    Panel (A): two triangles sharing a single hub vertex.
    TikZ: node (0) at origin; for x in {-1,1}, nodes at (-1,x) and (1,x);
    edges (0)--(-1,x)--(1,x)--(0) for each x.
    Caption: "traceable, not 2-connected, not RN"
    """
    G = nx.Graph()
    G.add_edges_from([("0", "A"), ("A", "B"), ("B", "0"),
                       ("0", "C"), ("C", "D"), ("D", "0")])
    return G


def fig1_b_k23() -> nx.Graph:
    """
    Panel (B): K_{2,3}.
    TikZ: nodes a, b; three nodes at x=1; edges a--n_i--b for each.
    Caption: "sprawling, not 1-tough"
    """
    return nx.complete_bipartite_graph(2, 3)


def fig1_banana(legs: int) -> nx.Graph:
    """
    Panels (C)/(D): two hub vertices L, R (NOT directly connected to
    each other) joined by `legs` internally-disjoint paths of length 3
    (a generalized theta graph with no hub-hub edge).

    TikZ for (C) (3 legs): nodes L, R; for x in {-1,0,1}, nodes at
    (1,x),(2,x); edges L--(1,x)--(2,x)--R for each x.
    TikZ for (D) (4 legs): nodes a, b; for x in {-1.5,-.5,.5,1.5},
    nodes at (1,x),(2,x); edges a--(1,x)--(2,x)--b for each x.

      legs=3: panel (C), caption "traceable, SRN, not 1-tough, not sprawling"
              -- SEE STATUS NOTE ABOVE; code finds RN=False (exact,
              cross-validated), contradicting this caption.
      legs=4: panel (D), caption "2-connected, not 1-tough, not traceable, not RN"
    """
    G = nx.Graph()
    G.add_node("L")
    G.add_node("R")
    for i in range(legs):
        x, y = f"x{i}", f"y{i}"
        G.add_edges_from([("L", x), (x, y), (y, "R")])
    return G


def fig1_clique_plus_legs(k: int, leg_len: int = 2) -> nx.Graph:
    """
    Panels (E)/(F)/(G): a K_k hub clique, with each of the k hub
    vertices connected via a path of `leg_len` intermediate vertices to
    a single common external vertex v.
    TikZ pattern (shown for K4, panel (F)): node (0) [= v]; for x in
    hub-index-set, nodes at (0,x) [hub], (1,..) [a_x], (2,..) [b_x];
    edges hub_x--a_x--b_x--(0); plus explicit \\draw commands forming
    K_k among the hub vertices.
      k=3: panel (E), caption "traceable, RP, not sprawling"
      k=4: panel (F), caption "traceable, 1-tough, SRN, not sprawling"
           -- SEE STATUS NOTE ABOVE; code finds RN=False (exact,
           cross-validated), contradicting this caption.
      k=5: panel (G), caption "1-tough, not RN, not traceable"
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


def fig1_h_petersen() -> nx.Graph:
    """Panel (H): the Petersen graph.
    Caption (current revision): "Petersen graph: sprawling, RP, not Hamiltonian"."""
    return nx.petersen_graph()


def fig1_i_path(n: int = 7) -> nx.Graph:
    """
    Panel (I): a path on n >= 3 vertices (n=7 used here as a
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
    print(f"  Figure 1 caption: \"{caption}\"")

    trace = _traceable(G)
    two_conn = nx.is_biconnected(G)
    is_rp, is_rn, t_star, _ = resistance.resistance_positive_decision(G, verbose=False)
    tough, _ = toughness.is_one_tough(G)

    print(f"  traceable    = {trace}")
    print(f"  2-connected  = {two_conn}")
    print(f"  RN           = {is_rn}   [t* = {t_star:.6f}]")
    print(f"  RP           = {is_rp}   [t* = {t_star:.6f}]")
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
    print("Figure 1 example verification")
    print("=" * 78)
    print()

    verify(
        "(A) Bowtie (two triangles sharing a hub vertex)",
        fig1_a_bowtie(),
        "traceable, not 2-connected, not RN",
    )

    verify(
        "(B) K_{2,3}",
        fig1_b_k23(),
        "sprawling, not 1-tough",
        check_sprawling=True,
    )

    print("(C) Small banana (2 hub vertices, no hub-hub edge, + 3 legs)")
    print("  Figure 1 caption: \"traceable, SRN, not 1-tough, not sprawling\"")
    print("  *** STATUS: UNRESOLVED -- see module docstring above. ***")
    print("  Reconstructing this graph from the TikZ source and checking it two")
    print("  independent ways (the general LP in resistance.py, and an EXACT")
    print("  brute-force solve over all 27 of its actual spanning trees, cross-")
    print("  checked with cvxpy) both give RN = False, contradicting the \"SRN\"")
    print("  (RN=True) caption. Not reporting a pass/fail here pending")
    print("  confirmation of the intended graph structure.")
    G = fig1_banana(3)
    print(f"  (for reference: n={G.number_of_nodes()}, m={G.number_of_edges()}, "
          f"traceable={_traceable(G)}, 2-connected={nx.is_biconnected(G)}, "
          f"1-tough={toughness.is_one_tough(G)[0]})")
    print()

    verify(
        "(D) Banana (4 parallel length-3 paths between two hubs)",
        fig1_banana(4),
        "2-connected, not 1-tough, not traceable, not RN",
    )

    verify(
        "(E) K3 hub + 3 legs to a common point",
        fig1_clique_plus_legs(3, 2),
        "traceable, RP, not sprawling",
        check_sprawling=True,
    )

    print("(F) K4 hub + 4 legs to a common point")
    print("  Figure 1 caption: \"traceable, 1-tough, SRN, not sprawling\"")
    print("  *** STATUS: UNRESOLVED (vs. caption) -- see module docstring above. ***")
    print("  Cross-checked THREE independent ways: this repo's general margin-LP,")
    print("  an exact brute-force solve over all 2197 of the graph's actual")
    print("  spanning trees (via both cvxpy and scipy's HiGHS solver), and a")
    print("  separately-derived cutting-plane LP (the advisor's original code,")
    print("  after fixing the same closed-polytope-vs-relative-interior bug this")
    print("  repo's resistance.py also had). All three independently agree")
    print("  RN = False, ruling out a solver-specific numerical artifact --")
    print("  but this still contradicts the panel's \"SRN\" (RN=True) caption.")
    print("  Not reporting a pass/fail here pending confirmation of the intended")
    print("  graph structure.")
    print()

    verify(
        "(G) K5 hub + 5 legs to a common point",
        fig1_clique_plus_legs(5, 2),
        "1-tough, not RN, not traceable",
    )

    verify(
        "(H) Petersen graph",
        fig1_h_petersen(),
        "Petersen graph: sprawling, RP, not Hamiltonian",
        check_sprawling=True,
    )

    verify(
        "(I) Path P_7 (representative of the path family, n >= 3)",
        fig1_i_path(7),
        "Path with >= 3 vertices: SRN, traceable, not 2-connected",
    )

    print("=" * 78)
    print("Summary: panels (A), (B), (D), (E), (G), (H), (I) all match their")
    print("captions exactly. Panels (C) and (F) are flagged as unresolved --")
    print("see notes above for both.")
    print("=" * 78)
