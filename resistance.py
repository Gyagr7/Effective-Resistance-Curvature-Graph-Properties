"""
resistance.py
=============

Decide whether a graph is resistance nonnegative (RN) or resistance
positive (RP), in the sense of Devriendt's discrete resistance
curvature (Devriendt, "Graphs with nonnegative resistance curvature",
Ann. Combin. 2025; Theorem 1 of the accompanying paper).

Method
------
By Theorem 1, letting

    t*(G) = min_{x in P(G)} max_{v in V(G)} x(E(v))

over the spanning tree polytope P(G), G is RN if t* <= 2 and RP if
t* < 2 (each up to the tolerance `tol_rp`).

The spanning tree polytope is described by
    x_e >= 0                     for every edge e
    sum_e x_e = |V(G)| - 1
    x(E[S]) <= |S| - 1           for every subtour set S (2 <= |S| <= n-1)

The subtour constraints are exponential in number, so t* is found by a
cutting-plane method: solve the LP with only the trivial constraints,
use a min s-t cut to find the most-violated subtour constraint (the
construction described in Svensson's notes on the spanning tree
polytope), add it, and repeat until no violation remains.
"""

from __future__ import annotations

from typing import Dict, List, Optional, Tuple

import cvxpy as cp
import networkx as nx


# ---------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------

def _canon_edge(u, v):
    """Canonical (order-independent) representation of an undirected edge."""
    return (u, v) if u <= v else (v, u)


def _build_edge_list(G: nx.Graph) -> List[Tuple]:
    """Canonical edge list for G, in a fixed order."""
    return [_canon_edge(u, v) for (u, v) in G.edges()]


def _incident_sums(nodes, edges, x_vals) -> Dict:
    """For each vertex v, the sum of x_e over edges e incident to v."""
    inc = {v: 0.0 for v in nodes}
    for (u, v), xe in zip(edges, x_vals):
        xe = float(xe)
        inc[u] += xe
        inc[v] += xe
    return inc


# ---------------------------------------------------------------------
# Subtour separation
# ---------------------------------------------------------------------

def _separate_subtour_via_mincut(
    G: nx.Graph,
    nodes: List,
    edges: List[Tuple],
    x_vals,
    sep_eps: float = 1e-8,
) -> Optional[Tuple[set, float, float]]:
    """
    Separation oracle for the subtour constraints x(E(S)) <= |S| - 1,
    over all nonempty proper subsets S.

    Uses the s-t min-cut construction described in Svensson's notes:
    add source s and sink t; arcs s -> v with capacity x(delta(v))/2;
    arcs v -> t with capacity 1; each undirected edge {u, v} becomes two
    directed arcs u -> v and v -> u, each with capacity x_e / 2. A
    minimum s-t cut then corresponds to a most-violated subtour
    inequality.

    To keep the empty set from trivially dominating every cut, a root
    vertex is forced into the source side (by setting its s -> root
    capacity to a very large number); repeating this for every choice
    of root guarantees any violated subtour set will be found by some
    iteration of the loop.

    Returns (S, violation_amount, cut_value) for the most-violated set
    found, or None if no subtour constraint is violated by more than
    `sep_eps`.
    """
    x_dict = {e: float(xe) for e, xe in zip(edges, x_vals)}
    x_total = sum(x_dict.values())
    n = len(nodes)

    inc = _incident_sums(nodes, edges, x_vals)

    # A safe "infinity": cut values are O(n + x_total), so this comfortably
    # dominates any finite cut.
    BIGM = 10.0 * (n + x_total + 1.0)

    best = None  # (S, violation_amount, cut_value)

    for root in nodes:
        H = nx.DiGraph()
        H.add_node("s")
        H.add_node("t")

        for v in nodes:
            cap_sv = BIGM if v == root else inc[v] / 2.0
            H.add_edge("s", v, capacity=cap_sv)
            H.add_edge(v, "t", capacity=1.0)

        for (u, v), xe in x_dict.items():
            cap = xe / 2.0
            H.add_edge(u, v, capacity=cap)
            H.add_edge(v, u, capacity=cap)

        cut_value, (S_side, _T_side) = nx.minimum_cut(
            H, "s", "t", capacity="capacity",
            flow_func=nx.algorithms.flow.preflow_push,
        )

        S = set(S_side)
        if "s" not in S:
            continue
        S.discard("s")

        if len(S) == 0 or len(S) == n:
            continue  # only proper, nonempty subsets are valid subtour sets

        x_inside = sum(
            xe for (u, v), xe in x_dict.items() if u in S and v in S
        )
        violation = x_inside - (len(S) - 1)

        if violation > sep_eps and (best is None or violation > best[1]):
            best = (S, violation, cut_value)

    return best


# ---------------------------------------------------------------------
# Main decision procedure
# ---------------------------------------------------------------------

def resistance_positive_decision(
    G: nx.Graph,
    solver: str = "SCS",
    scs_eps: float = 1e-7,
    scs_max_iters: int = 200_000,
    max_outer_iters: int = 200,
    sep_eps: float = 1e-15,
    tol_rp: float = 1e-6,
    verbose: bool = True,
):
    """
    Decide the RN / RP status of a connected graph G by minimizing
    t = max_v x(E(v)) over the spanning tree polytope via a
    cutting-plane LP.

    Parameters
    ----------
    solver, scs_eps, scs_max_iters : passed through to cvxpy's SCS solver.
    max_outer_iters : cap on the number of cutting-plane rounds.
    sep_eps : minimum violation the separator will act on.
    tol_rp : tolerance used when comparing t* to 2 for the RN/RP decision.
    verbose : print per-iteration diagnostics.

    Returns
    -------
    rp : bool -- True if t* < 2 (within tol_rp)
    rn : bool -- True if t* <= 2 (within tol_rp)
    t_star : float -- the optimal (closed-polytope) value found
    x_dict : dict -- edge -> x_e at the solution
    """
    if not nx.is_connected(G):
        raise ValueError("Graph must be connected (no spanning tree otherwise).")

    nodes = list(G.nodes())
    n = len(nodes)

    edges = _build_edge_list(G)
    m = len(edges)

    # Edge incidence list, for the degree constraints below.
    incident = {v: [] for v in nodes}
    for i, (u, v) in enumerate(edges):
        incident[u].append(i)
        incident[v].append(i)

    # Variables: x_e for each edge, and t = max expected degree.
    x = cp.Variable(m)
    t = cp.Variable()

    constraints = [
        x >= 0,
        cp.sum(x) == n - 1,
        t >= 0,
    ]
    constraints += [cp.sum(x[incident[v]]) <= t for v in nodes]

    subtour_constraints = []

    for it in range(max_outer_iters):
        prob = cp.Problem(cp.Minimize(t), constraints + subtour_constraints)
        prob.solve(solver=solver, verbose=False, eps=scs_eps, max_iters=scs_max_iters)

        if prob.status not in ("optimal", "optimal_inaccurate"):
            raise RuntimeError(f"LP solve failed: status={prob.status}")

        x_vals = x.value
        t_star = float(t.value)

        viol = _separate_subtour_via_mincut(G, nodes, edges, x_vals, sep_eps=sep_eps)

        if verbose:
            if viol is None:
                print(f"iter={it}, t={t_star:.12f}, subtour_violation=None")
            else:
                S, vamt, _cutv = viol
                print(f"iter={it}, t={t_star:.12f}, subtour_violation={vamt:.3e}, |S|={len(S)}")

        if viol is None:
            break

        S, _, _ = viol
        idxs = [i for i, (u, v) in enumerate(edges) if u in S and v in S]
        subtour_constraints.append(cp.sum(x[idxs]) <= len(S) - 1)

    x_dict = {edges[i]: float(x.value[i]) for i in range(m)}

    # RN / RP decision, with a tolerance around t* = 2.
    rn = t_star <= 2.0 + tol_rp
    rp = t_star < 2.0 - tol_rp

    # Hard certificate: an imbalanced bipartite graph can never be RP
    # (this guards against solver noise producing a false RP=True on
    # such graphs).
    if nx.is_bipartite(G):
        color = nx.algorithms.bipartite.color(G)
        a = sum(1 for v in color if color[v] == 0)
        b = len(color) - a
        if a != b:
            rp = False

    return rp, rn, t_star, x_dict


if __name__ == "__main__":

    G = nx.Graph()
    G.add_edges_from([
        # K4 hub clique
        (0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3),

        # Leg from hub 0: 0 -- 4 -- 5 -- 12
        (0, 4), (4, 5), (5, 12),

        # Leg from hub 1: 1 -- 6 -- 7 -- 12
        (1, 6), (6, 7), (7, 12),

        # Leg from hub 2: 2 -- 8 -- 9 -- 12
        (2, 8), (8, 9), (9, 12),

        # Leg from hub 3: 3 -- 10 -- 11 -- 12
        (3, 10), (10, 11), (11, 12),
    ])

    print("n =", G.number_of_nodes())
    print("m =", G.number_of_edges())
    print(list(G.edges()))

    rp, rn, t_star, x = resistance_positive_decision(G, verbose=True)

    print("\nRESULT")
    print("t* =", t_star)
    print("RN?", rn)
    print("RP?", rp)
