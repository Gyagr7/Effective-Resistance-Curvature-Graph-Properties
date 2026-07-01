"""
resistance.py
=============

Decide whether a graph G is resistance nonnegative (RN) or resistance
positive (RP), in the sense of Devriendt's discrete resistance curvature
(see Devriendt, "Graphs with nonnegative resistance curvature", Ann.
Combin. 2025, and Theorem 1 of the accompanying paper).

Background
----------
By Theorem 1 of the paper, G is RN (resp. RP) if and only if there is a
point x in the spanning tree polytope P(G) with x(E(v)) <= 2 (resp. < 2)
for every vertex v. Equivalently, letting

    t*(G) := min_{x in P(G)} max_{v in V(G)} x(E(v))

we have:
    * G is RN  iff  t*(G) <= 2
    * G is RP  iff  t*(G) <  2

The spanning tree polytope P(G) is described by
    x_e >= 0                              for every edge e
    sum_e x_e = |V(G)| - 1
    x(E[S]) <= |S| - 1                    for every proper nonempty S

The subtour constraints are exponential in number, so t*(G) is computed
by a cutting-plane method: solve the LP relaxation with only the trivial
constraints, find (via a min s-t cut computation) the most violated
subtour constraint, add it, and repeat until no violation remains.

This module is a cleaned-up version of the exploratory code used while
writing the paper.
"""

from __future__ import annotations

from itertools import combinations
from typing import Dict, List, Tuple

import cvxpy as cp
import networkx as nx


def _canon_edge(u, v):
    return (u, v) if str(u) <= str(v) else (v, u)


def _edge_list(G: nx.Graph) -> List[Tuple]:
    return [_canon_edge(u, v) for u, v in G.edges()]


def _incident_sums(nodes, edges, x_vals) -> Dict:
    inc = {v: 0.0 for v in nodes}
    for (u, v), xe in zip(edges, x_vals):
        xe = float(xe)
        inc[u] += xe
        inc[v] += xe
    return inc


def _separate_subtour(G: nx.Graph, nodes, edges, x_vals, eps: float = 1e-7):
    """
    Separation oracle for the subtour constraints x(E[S]) <= |S| - 1.

    Builds the standard s-t min-cut network (Svensson-style construction):
    arcs s->v with capacity x(E(v))/2, arcs v->t with capacity 1, and each
    undirected edge {u, v} becomes two directed arcs u->v, v->u each with
    capacity x_e/2. A most-violated subtour set corresponds to a min s-t
    cut, found by forcing each vertex in turn to lie on the source side.

    Returns (S, violation) for the most violated set found, or None if no
    subtour constraint is violated by more than eps.
    """
    x_dict = {e: float(xe) for e, xe in zip(edges, x_vals)}
    inc = _incident_sums(nodes, edges, x_vals)
    n = len(nodes)
    x_total = sum(x_dict.values())
    BIGM = 10.0 * (n + x_total + 1.0)

    best = None  # (S, violation)

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

        _, (S_side, _) = nx.minimum_cut(
            H, "s", "t", capacity="capacity",
            flow_func=nx.algorithms.flow.preflow_push,
        )

        S = set(S_side)
        S.discard("s")
        if len(S) == 0 or len(S) == n:
            continue

        x_inside = sum(
            xe for (u, v), xe in x_dict.items() if u in S and v in S
        )
        violation = x_inside - (len(S) - 1)

        if violation > eps and (best is None or violation > best[1]):
            best = (S, violation)

    return best


def resistance_threshold(
    G: nx.Graph, max_iters: int = 200, eps: float = 1e-7, verbose: bool = False
):
    """
    Compute t*(G) = min_{x in P(G)} max_v x(E(v)) via cutting planes.

    Returns (t_star, x) where x is a dict edge -> weight attaining (or
    numerically close to attaining) the optimum.
    """
    nodes = list(G.nodes())
    edges = _edge_list(G)
    n, m = len(nodes), len(edges)
    if n < 2:
        raise ValueError("Graph must have at least 2 vertices")

    x = cp.Variable(m, nonneg=True)
    t = cp.Variable()

    edge_index = {e: i for i, e in enumerate(edges)}
    incident = {v: [] for v in nodes}
    for (u, v), i in edge_index.items():
        incident[u].append(i)
        incident[v].append(i)

    base_constraints = [cp.sum(x) == n - 1]
    base_constraints += [
        cp.sum(x[incident[v]]) <= t for v in nodes
    ]

    subtour_constraints = []
    x_vals = None

    for it in range(max_iters):
        problem = cp.Problem(
            cp.Minimize(t), base_constraints + subtour_constraints
        )
        problem.solve()

        if x.value is None:
            raise RuntimeError("LP failed to solve; check graph connectivity")

        x_vals = x.value
        violated = _separate_subtour(G, nodes, edges, x_vals, eps=eps)

        if violated is None:
            if verbose:
                print(f"Converged after {it + 1} LP solve(s). t* = {t.value:.6f}")
            break

        S, violation = violated
        S_idx = [
            i for (u, v), i in edge_index.items() if u in S and v in S
        ]
        subtour_constraints.append(cp.sum(x[S_idx]) <= len(S) - 1)

        if verbose:
            print(f"iter {it}: added subtour cut on |S|={len(S)}, violation={violation:.4f}")
    else:
        raise RuntimeError(f"Cutting-plane method did not converge in {max_iters} iterations")

    x_dict = {e: float(x_vals[i]) for e, i in edge_index.items()}
    return float(t.value), x_dict


def resistance_positive_decision(G: nx.Graph, verbose: bool = False, tol: float = 1e-6):
    """
    Decide RN / RP status of a connected graph G.

    Returns (is_rp, is_rn, t_star, x) where:
      * is_rn  = True iff t_star <= 2 + tol
      * is_rp  = True iff t_star <  2 - tol
      * x is the (approximately) optimal point of the spanning tree
        polytope attaining t_star.
    """
    if not nx.is_connected(G):
        raise ValueError("G must be connected")

    t_star, x = resistance_threshold(G, verbose=verbose)

    is_rn = t_star <= 2 + tol
    is_rp = t_star < 2 - tol

    if verbose:
        print(f"t* = {t_star:.6f}  ->  RN={is_rn}, RP={is_rp}")

    return is_rp, is_rn, t_star, x


if __name__ == "__main__":
    # Sanity check: the Petersen graph is known to be RP but not
    # Hamiltonian (Devriendt 2025; see Theorem 1 discussion in the paper).
    G = nx.petersen_graph()
    rp, rn, t_star, x = resistance_positive_decision(G, verbose=True)
    print(f"Petersen graph: RN={rn}, RP={rp}, t*={t_star:.4f}  (expected RP=True)")
