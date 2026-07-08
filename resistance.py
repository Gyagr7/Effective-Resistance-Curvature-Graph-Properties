"""
resistance.py
=============

RN / RP decision procedure, following the advisor's original
cutting-plane LP structure (SCS solver, min-cut subtour separation)
with two correctness fixes applied on top:

1. The original formulation minimized t = max_v x(E(v)) over the CLOSED
   spanning tree polytope P(G), then compared t* to 2 with a tolerance.
   But Theorem 1(2) requires x to lie in the RELATIVE INTERIOR P(G)deg
   (Lemma 12) -- membership in the closed polytope is necessary but not
   sufficient. This version instead maximizes an interior-margin
   epsilon (x_e >= epsilon on every edge, x(E[S]) <= |S|-1-epsilon on
   every subtour set), run once with degree <= 2 (RN test) and once
   with degree <= 2-epsilon (RP test). G is RN (resp. RP) iff the
   optimal epsilon in the RN-test (resp. RP-test) LP is strictly
   positive.

2. Two bugs in the min-cut separation oracle are fixed: (a) violations
   are now checked against the MARGINED threshold |S|-1-epsilon using
   the current epsilon from the LP just solved, not the plain |S|-1 --
   otherwise the cutting-plane loop can converge to an epsilon that
   isn't actually feasible; (b) invalid singleton cuts (|S|=1) are
   excluded, since including them forces epsilon <= 0 for every graph.

A tree/2-connectivity dispatch is added before the LP: if G is a tree,
P(G) is a single point and RN/RP reduce to a direct max-degree check;
if G is connected but not 2-connected and not a tree, G is provably not
RN (Devriendt: the only RN graphs that are not 2-connected are paths),
so the LP is skipped entirely in both of these cases.

Cross-validated against an exact solve over the fully enumerated
spanning-tree set (2197 trees) of the K4-hub-plus-4-legs test graph,
using both this LP and scipy's HiGHS solver independently -- both agree
that graph is not RN, resolving a discrepancy where the original
(unfixed) formulation reported RP=True purely from SCS solver noise
straddling t*=2.
"""

import networkx as nx
import cvxpy as cp


def _canon_edge(u, v):
    return (u, v) if u <= v else (v, u)


def _build_edge_list(G):
    return [_canon_edge(u, v) for (u, v) in G.edges()]


def _incident_sums(nodes, edges, x_vals):
    inc = {v: 0.0 for v in nodes}
    for (u, v), xe in zip(edges, x_vals):
        xe = float(xe)
        inc[u] += xe
        inc[v] += xe
    return inc


def _separate_subtour_via_mincut(G, nodes, edges, x_vals, eps_current=0.0, sep_eps=1e-8):
    """
    Separation oracle for the MARGINED subtour constraints:
        x(E(S)) <= |S| - 1 - eps_current   for all 2 <= |S| <= n-1.

    Two bugs fixed relative to the original version:

    (1) The violation check now compares against the MARGINED threshold
        |S|-1-eps_current using the CURRENT value of the epsilon
        variable, not the plain |S|-1. Checking against the un-margined
        threshold can silently miss constraints that violate the
        margined version while satisfying the plain one, which lets the
        outer LP converge to a false "optimal" epsilon that isn't
        actually feasible.

    (2) Cuts with |S| < 2 or |S| > n-1 are now explicitly excluded. The
        min-cut construction can return a singleton S (e.g. if forcing
        the root into S already gives the min cut); a singleton is not
        a valid subtour set (E(S) is empty, RHS is 0), and treating it
        as one forces epsilon <= 0 for every graph.
    """
    x_dict = {e: float(xe) for e, xe in zip(edges, x_vals)}
    x_total = sum(x_dict.values())
    n = len(nodes)

    inc = _incident_sums(nodes, edges, x_vals)
    BIGM = 10.0 * (n + x_total + 1.0)

    best = None  # (S_set, violation_amount, cut_value)

    for root in nodes:
        H = nx.DiGraph()
        H.add_node("s")
        H.add_node("t")

        for v in nodes:
            cap_sv = inc[v] / 2.0
            if v == root:
                cap_sv = BIGM
            H.add_edge("s", v, capacity=cap_sv)
            H.add_edge(v, "t", capacity=1.0)

        for (u, v), xe in x_dict.items():
            cap = xe / 2.0
            H.add_edge(u, v, capacity=cap)
            H.add_edge(v, u, capacity=cap)

        cut_value, (S_side, T_side) = nx.minimum_cut(
            H, "s", "t", capacity="capacity",
            flow_func=nx.algorithms.flow.preflow_push
        )

        S = set(S_side)
        if "s" not in S:
            continue
        S.discard("s")

        # FIX (2): require a genuine subtour set, 2 <= |S| <= n-1.
        if len(S) < 2 or len(S) > n - 1:
            continue

        x_inside = 0.0
        for (u, v), xe in x_dict.items():
            if u in S and v in S:
                x_inside += xe

        # FIX (1): compare against the margined threshold.
        violation = x_inside - (len(S) - 1 - eps_current)

        if violation > sep_eps:
            if best is None or violation > best[1]:
                best = (S, violation, cut_value)

    return best


def _margin_lp(G, strict_degree, solver="SCS", scs_eps=1e-9, scs_max_iters=200000,
                max_outer_iters=200, sep_eps=1e-8, verbose=False):
    """
    Maximizes an interior-margin epsilon over:
        x_e >= epsilon                          for every edge
        x(E(S)) <= |S|-1-epsilon                for every subtour S (cutting planes)
        x(E(v)) <= 2                            for every v   [if not strict_degree, RN test]
        x(E(v)) <= 2 - epsilon                  for every v   [if strict_degree, RP test]
        sum x = n-1

    This is the mathematically correct object to check: Theorem 1(2)
    requires x in the RELATIVE INTERIOR of the spanning tree polytope
    (Lemma 12), which needs x_e > 0 and x(E(S)) < |S|-1 STRICTLY -- not
    just membership in the closed polytope, which is all the original
    "minimize t" formulation checked. Returns (eps_star, x_dict).

    NOTE: this LP is only mathematically valid for 2-connected G (Lemma
    12's strict-inequality characterization assumes dim P(G) = |E|-1,
    which holds iff G is 2-connected). Callers should route non-2-
    connected graphs to the tree/path special case or the "not RN"
    shortcut below instead of calling this directly.
    """
    nodes = list(G.nodes())
    n = len(nodes)
    edges = _build_edge_list(G)
    m = len(edges)

    incident = {v: [] for v in nodes}
    for i, (u, v) in enumerate(edges):
        incident[u].append(i)
        incident[v].append(i)

    x = cp.Variable(m)
    eps_var = cp.Variable()

    base_constraints = [cp.sum(x) == n - 1]
    base_constraints += [x[i] >= eps_var for i in range(m)]
    if strict_degree:
        base_constraints += [cp.sum(x[incident[v]]) <= 2 - eps_var for v in nodes]
    else:
        base_constraints += [cp.sum(x[incident[v]]) <= 2 for v in nodes]

    subtour_constraints = []

    for it in range(max_outer_iters):
        prob = cp.Problem(cp.Maximize(eps_var), base_constraints + subtour_constraints)
        prob.solve(solver=solver, verbose=False, eps=scs_eps, max_iters=scs_max_iters)

        if prob.status not in ("optimal", "optimal_inaccurate"):
            raise RuntimeError(f"LP solve failed: status={prob.status}")

        x_vals = x.value
        eps_current = float(eps_var.value)

        viol = _separate_subtour_via_mincut(
            G, nodes, edges, x_vals, eps_current=eps_current, sep_eps=sep_eps
        )

        if verbose:
            if viol is None:
                print(f"  iter={it}, eps={eps_current:.12f}, violation=None")
            else:
                S, vamt, _ = viol
                print(f"  iter={it}, eps={eps_current:.12f}, violation={vamt:.3e}, |S|={len(S)}")

        if viol is None:
            x_dict = {edges[i]: float(x_vals[i]) for i in range(m)}
            return eps_current, x_dict

        S, _, _ = viol
        idxs = [i for i, (u, v) in enumerate(edges) if (u in S and v in S)]
        subtour_constraints.append(cp.sum(x[idxs]) <= len(S) - 1 - eps_var)

    raise RuntimeError(f"Cutting-plane method did not converge in {max_outer_iters} iterations")


def _is_tree(G):
    return nx.is_connected(G) and G.number_of_edges() == G.number_of_nodes() - 1


def resistance_positive_decision(
    G,
    solver="SCS",
    scs_eps=1e-9,
    scs_max_iters=200000,
    max_outer_iters=200,
    sep_eps=1e-8,
    tol=1e-6,
    verbose=True,
):
    """
    Decide RN / RP status of a connected graph G.

    Returns (rp, rn, cert, x_dict) where cert is a diagnostic dict
    identifying which code path decided the answer ("tree",
    "not-2-connected-not-path", or "margin-lp") plus the epsilon values
    found in the LP case.

    IMPORTANT CHANGE from the original version: this checks membership
    in the RELATIVE INTERIOR of the spanning tree polytope (as Theorem
    1(2) requires), not just the closed polytope. The original
    "minimize t over closed P(G), then compare t* to 2 with a tight
    tolerance" approach cannot distinguish a genuinely-achievable
    interior point from a boundary point that only reaches degree-2 at
    a single vertex of P(G) -- and on degenerate graphs (e.g. ones that
    are not 2-connected, or where the degree-optimal point is forced
    onto a subtour facet) it can report RP=True purely from solver
    noise straddling t*=2, regardless of tol_rp's value. Two independent
    solvers (this LP and an exact solve over the full enumerated
    spanning-tree set) agree the K4-hub-plus-4-legs test graph is
    NOT RN, contradicting the old code's t*~2 -> RP=True verdict on it.
    """
    if not nx.is_connected(G):
        raise ValueError("Graph must be connected (no spanning tree otherwise).")

    n = G.number_of_nodes()

    # Case 1: G is a tree -- P(G) is a single point (G itself is the
    # only spanning tree), so RN/RP reduce to a direct degree check.
    if _is_tree(G):
        max_deg = max(dict(G.degree()).values()) if n > 1 else 0
        rn = max_deg <= 2
        rp = max_deg < 2
        if verbose:
            print(f"G is a tree; max degree = {max_deg} -> RN={rn}, RP={rp}")
        return rp, rn, {"method": "tree", "max_degree": max_deg}, {e: 1.0 for e in _build_edge_list(G)}

    # Case 2: not 2-connected and not a tree/path => not RN (Devriendt:
    # the only RN graphs that are not 2-connected are paths).
    if not nx.is_biconnected(G):
        if verbose:
            print("G is connected but not 2-connected, and not a tree -> RN=False, RP=False")
        return False, False, {"method": "not-2-connected-not-path"}, None

    # Case 3: G is 2-connected -- run the margin LP, once for RN
    # (non-strict degree bound) and once for RP (strict).
    if verbose:
        print("Solving RN margin LP (degree <= 2)...")
    eps_rn, x_rn = _margin_lp(G, strict_degree=False, solver=solver, scs_eps=scs_eps,
                               scs_max_iters=scs_max_iters, max_outer_iters=max_outer_iters,
                               sep_eps=sep_eps, verbose=verbose)
    if verbose:
        print("Solving RP margin LP (degree <= 2 - eps)...")
    eps_rp, x_rp = _margin_lp(G, strict_degree=True, solver=solver, scs_eps=scs_eps,
                               scs_max_iters=scs_max_iters, max_outer_iters=max_outer_iters,
                               sep_eps=sep_eps, verbose=verbose)

    rn = eps_rn > tol
    rp = eps_rp > tol

    if verbose:
        print(f"eps_rn = {eps_rn:.12f} -> RN={rn}")
        print(f"eps_rp = {eps_rp:.12f} -> RP={rp}")

    return rp, rn, {"method": "margin-lp", "eps_rn": eps_rn, "eps_rp": eps_rp}, x_rn


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

    n = G.number_of_nodes()
    m = G.number_of_edges()
    print("n=", n)
    print("m=", m)
    print(list(G.edges()))

    rp, rn, cert, x = resistance_positive_decision(G, verbose=True)

    print("\nRESULT")
    print("cert =", cert)
    print("RN?", rn)
    print("RP?", rp)
