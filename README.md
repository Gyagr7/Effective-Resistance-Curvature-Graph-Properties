# Code for "On Some Structural Properties of Graphs with Non-Negative Resistance Curvature"

This repository contains the code used to check the RN/RP status,
sprawling property, and toughness of graphs discussed in the paper:

> G. Agrahari, C. Bibby, S. Boros, H. Garcia, F. Heiderscheidt, Z. Wang.
> *On some structural properties of graphs with non-negative resistance
> curvature.* [arXiv link / venue, once available]

## Contents

| File | Purpose |
|---|---|
| `resistance.py` | Decide whether a graph is resistance nonnegative (RN) or resistance positive (RP) (Theorem 1), via a cutting-plane LP over the spanning tree polytope. |
| `sprawling.py` | Decide whether a graph is *sprawling* (Section 5), a sufficient condition for RN (Theorem 8); every "sprawling" verdict comes with an explicit, independently re-verified witness collection. |
| `toughness.py` | Compute exact (vertex) toughness and check 1-toughness, by brute force. |
| `examples.py` | Graph constructions referenced in the paper (Petersen graph, grid graphs, the `G_t(s_1,...,s_t)` family from Theorem 4, etc.). |
| `verify_figure2_examples.py` | Runs all three checkers against every example graph in Figure 2, labeled by panel letter, with real captured output baked into the file so a reviewer can read the results without running anything. |

## Installation

```bash
pip install -r requirements.txt
```

## Usage

```python
import networkx as nx
import resistance, sprawling, toughness

G = nx.petersen_graph()

# RN / RP decision (Theorem 1)
is_rp, is_rn, cert, x = resistance.resistance_positive_decision(G, verbose=False)
# cert is a dict identifying which code path decided the answer:
#   {"method": "tree", ...}                     -- G is a tree, decided by a direct degree check
#   {"method": "not-2-connected-not-path", ...}  -- G is provably not RN (Devriendt)
#   {"method": "margin-lp", "eps_rn": ..., "eps_rp": ...}  -- decided by the LP below

# Sprawling decision (Section 5) -- takes an adjacency-dict graph
is_sprawl, info = sprawling.is_sprawling(sprawling.from_networkx(G))
# info is the explicit witness collection S if is_sprawl=True,
# or the specific failing condition/set if False

# Toughness (used in the proof of Theorem 4)
tau, cut_set = toughness.toughness(G)
is_1_tough, witness = toughness.is_one_tough(G)
```

### How `resistance.py` decides RN / RP

Theorem 1(2) requires a point in the **relative interior** of the
spanning tree polytope $P(G)$ (Lemma 12) -- not just the closed polytope
-- so `resistance_positive_decision` dispatches on structure before
reaching for an LP:

1. **G is a tree** (unique spanning tree = G itself): $P(G)$ is a single
   point, so RN/RP reduce to a direct max-degree check.
2. **G is connected but not 2-connected, and not a tree**: G is
   provably *not* RN (Devriendt: the only RN graphs that are not
   2-connected are paths), decided with no LP at all.
3. **G is 2-connected**: solves a cutting-plane LP that *maximizes an
   interior-margin* $\varepsilon$ (not the closed-polytope minimum
   $t^*$), once with degree $\le 2$ (RN test) and once with degree
   $\le 2-\varepsilon$ (RP test). G is RN (resp. RP) iff the optimal
   $\varepsilon$ in the RN-test (resp. RP-test) LP is strictly positive.

This distinction matters: an earlier, simpler version of this code
minimized $t^* = \max_v x(E(v))$ over the *closed* polytope and compared
$t^*$ to 2 with a tolerance. That's only a necessary condition, and it
can report a false RP=True purely from solver noise straddling
$t^*=2$ on graphs where the true optimum sits exactly on $P(G)$'s
boundary. The current version was cross-validated by an exact solve
over a fully enumerated spanning-tree set (2197 trees, via both cvxpy
and scipy's HiGHS solver independently) on one such test case.

Each module can also be run directly (`python resistance.py`, etc.) to
execute a few built-in sanity checks against known examples from the
paper. `python examples.py` reproduces the key computational claims
end-to-end, including:

- the Petersen graph is RP;
- grid graphs `P_m x P_n` are sprawling (Theorem 17);
- the `G_5(1,1,1,1,1)` construction from Theorem 4 is 1-tough but not RN,
  disproving Fiedler's conjecture that every 1-tough graph is RP.

`python verify_figure2_examples.py` checks all nine panels of Figure 2:
panels (A)-(D), (F), (G), (I) all match their captions exactly. Panel
(E) is flagged as **unresolved** -- reconstructing it from the TikZ
source and checking it two independent ways gives RN=False, which
contradicts its "SRN" caption; see that file's module docstring and the
printed note for panel (E) for details, pending confirmation of the
intended graph structure. Panel (H) (the Thomassen 34-graph) is
intentionally out of scope for this file -- its RP status is already
established directly in the paper by Theorem 5's hand-constructed
rational edge weighting, independent of any code here.

## Scope and caveats

- `resistance.py` solves an LP via cutting planes and is numerically
  exact up to solver tolerance; it works comfortably on graphs with
  dozens of vertices.
- `sprawling.py` and `toughness.py` are brute-force (they enumerate
  Hamiltonian paths / vertex subsets respectively) and are only
  practical for small graphs -- exactly the regime used for the
  examples in the paper (roughly n <= 12). Both `toughness_family` and
  `build_minimal_tough_graph` in `examples.py` build the same
  underlying construction from Theorem 4 / Lemma 15 (a hub connected to
  a clique via subdivided spokes) -- the former with per-branch lengths
  and string labels, the latter with equal branch lengths and integer
  labels -- kept as two entry points since both conventions have been
  used across the project.
- The verification that the Thomassen 34-graph is RP (Theorem 5) uses a
  hand-constructed rational weighting rather than any code in this
  repo; see the proof of Theorem 5 in the paper.
- Panel (E) of Figure 2 (the K4-hub-plus-4-legs example) has an
  open discrepancy between its "SRN" caption and this code's RN=False
  verdict -- see `verify_figure2_examples.py` for the full account.

## Requirements

See `requirements.txt`. Core dependencies are `networkx` and `cvxpy`
(with a standard open-source LP solver; `resistance.py` uses SCS by
default, which ships with cvxpy).

## Citation

If you use this code, please cite the paper:

```bibtex
@article{agrahari2026resistance,
  title   = {On some structural properties of graphs with non-negative resistance curvature},
  author  = {Agrahari, Gyaneshwar and Bibby, Christin and Boros, Sean and Garcia, Hailey and Heiderscheidt, Fernando and Wang, Zhiyu},
  year    = {2026},
  note    = {arXiv preprint}
}
```
