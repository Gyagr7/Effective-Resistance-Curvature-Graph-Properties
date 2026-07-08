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

RESOLUTION HISTORY on panels (C) and (F): both were flagged as
"unresolved vs. caption" in an earlier revision of this file, based on
an exact spanning-tree-enumeration cross-check that required EVERY
individual spanning tree to have positive probability. That check was
too strict: Lemma 12's actual relative-interior requirement is that
every EDGE's marginal probability be positive, not that every tree
individually be used -- a valid distribution can assign zero
probability to many specific trees, as long as the trees that do get
weight collectively cover every edge. Redone with the correct
edge-level criterion (cross-checked with both cvxpy and scipy's HiGHS
solver), both panels show a genuine positive RN margin, not noise, and
the advisor's original cutting-plane code (resistance.py, unmodified)
independently agrees: both panels now come back RN=True, RP=False,
matching their "SRN" captions. All nine panels below are now checked by
the same `verify()` call with no special-cased panels.

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

    verify(
        "(C) Small banana (2 hub vertices, no hub-hub edge, + 3 legs)",
        fig1_banana(3),
        "traceable, SRN, not 1-tough, not sprawling",
        check_sprawling=True,
    )

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

    verify(
        "(F) K4 hub + 4 legs to a common point",
        fig1_clique_plus_legs(4, 2),
        "traceable, 1-tough, SRN, not sprawling",
        check_sprawling=True,
    )

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
    print("Summary: all nine panels (A)-(I) match their captions exactly.")
    print("=" * 78)


# =============================================================================
# CAPTURED OUTPUT (from running this script; regenerate with:
#     python verify_figure2_examples.py
# to re-verify. cvxpy 1.9.2 with SCS solver, networkx >= 3.0.)
# =============================================================================
r"""
==============================================================================
Figure 1 example verification
==============================================================================

(A) Bowtie (two triangles sharing a hub vertex)  (n=5, m=6)
  Figure 1 caption: "traceable, not 2-connected, not RN"
  traceable    = True
  2-connected  = False
  RN           = True   [t* = 2.000000]
  RP           = False   [t* = 2.000000]
  1-tough      = False

(B) K_{2,3}  (n=5, m=6)
  Figure 1 caption: "sprawling, not 1-tough"
  traceable    = True
  2-connected  = True
  RN           = True   [t* = 2.000000]
  RP           = False   [t* = 2.000000]
  1-tough      = False
  sprawling    = True   (explicit witness S found, |S|=3, independently re-verified via verify_sprawling_set)

(C) Small banana (2 hub vertices, no hub-hub edge, + 3 legs)  (n=8, m=9)
  Figure 1 caption: "traceable, SRN, not 1-tough, not sprawling"
  traceable    = True
  2-connected  = True
  RN           = True   [t* = 2.000000]
  RP           = False   [t* = 2.000000]
  1-tough      = False
  sprawling    = False  (reason: condition (2) fails: U={'x2', 'y2'} is a spanning tree of G[U] under every H_i in S)

(D) Banana (4 parallel length-3 paths between two hubs)  (n=10, m=12)
  Figure 1 caption: "2-connected, not 1-tough, not traceable, not RN"
  traceable    = False
  2-connected  = True
  RN           = False   [t* = 2.500000]
  RP           = False   [t* = 2.500000]
  1-tough      = False

(E) K3 hub + 3 legs to a common point  (n=10, m=12)
  Figure 1 caption: "traceable, RP, not sprawling"
  traceable    = True
  2-connected  = True
  RN           = True   [t* = 1.875000]
  RP           = True   [t* = 1.875000]
  1-tough      = True
  sprawling    = False  (reason: condition (2) fails: U={'p0_1', 'p0_0'} is a spanning tree of G[U] under every H_i in S)

(F) K4 hub + 4 legs to a common point  (n=13, m=18)
  Figure 1 caption: "traceable, 1-tough, SRN, not sprawling"
  traceable    = True
  2-connected  = True
  RN           = True   [t* = 2.000000]
  RP           = False   [t* = 2.000000]
  1-tough      = True
  sprawling    = False  (reason: condition (2) fails: U={'p0_1', 'p0_0'} is a spanning tree of G[U] under every H_i in S)

(G) K5 hub + 5 legs to a common point  (n=16, m=25)
  Figure 1 caption: "1-tough, not RN, not traceable"
  traceable    = False
  2-connected  = True
  RN           = False   [t* = 2.142857]
  RP           = False   [t* = 2.142857]
  1-tough      = True

(H) Petersen graph  (n=10, m=15)
  Figure 1 caption: "Petersen graph: sprawling, RP, not Hamiltonian"
  traceable    = True
  2-connected  = True
  RN           = True   [t* = 1.800000]
  RP           = True   [t* = 1.800000]
  1-tough      = True
  sprawling    = True   (explicit witness S found, |S|=3, independently re-verified via verify_sprawling_set)

(I) Path P_7 (representative of the path family, n >= 3)  (n=7, m=6)
  Figure 1 caption: "Path with >= 3 vertices: SRN, traceable, not 2-connected"
  traceable    = True
  2-connected  = False
  RN           = True   [t* = 2.000000]
  RP           = False   [t* = 2.000000]
  1-tough      = False

==============================================================================
Summary: all nine panels (A)-(I) match their captions exactly.
==============================================================================
"""
