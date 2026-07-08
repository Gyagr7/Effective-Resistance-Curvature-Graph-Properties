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
| `verify_figure2_examples.py` | Runs all three checkers against every example graph in Figure 1 (the paper's main examples figure), labeled by panel letter, with real captured output baked into the file so a reviewer can read the results without running anything. |

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
rp, rn, t_star, x = resistance.resistance_positive_decision(G, verbose=False)
# t_star is the closed-polytope value t* = min_x max_v x(E(v)) for
# 2-connected G; it is None for the tree / not-2-connected structural
# pre-check cases (see below), which decide RN/RP without an LP at all.
# rn = (t_star <= 2 + tol_rp), rp = (t_star < 2 - tol_rp) when t_star is used

# Sprawling decision (Section 5) -- takes an adjacency-dict graph
is_sprawl, info = sprawling.is_sprawling(sprawling.from_networkx(G))
# info is the explicit witness collection S if is_sprawl=True,
# or the specific failing condition/set if False

# Toughness (used in the proof of Theorem 4)
tau, cut_set = toughness.toughness(G)
is_1_tough, witness = toughness.is_one_tough(G)
```

### How `resistance.py` decides RN / RP

`resistance_positive_decision` follows the advisor's original
cutting-plane LP structure for 2-connected graphs: it minimizes
$t^* = \max_v x(E(v))$ over the closed spanning tree polytope $P(G)$ via
a min-cut-based subtour separator (SCS solver), then classifies:

$$\text{RN} \iff t^* \le 2 + \texttt{tol\_rp}, \qquad \text{RP} \iff t^* < 2 - \texttt{tol\_rp}$$

with a hard bipartite-imbalance certificate (`RP=False` whenever $G$ is
bipartite with unequal parts) as an extra guard against solver noise.
Defaults: `sep_eps=1e-15`, `tol_rp=1e-6`.

Two structural cases are checked **before** the LP, because $t^*$
cannot correctly decide RN on them at any tolerance -- not a
solver-precision issue, but a case where $t^*$ is the wrong quantity to
compare:

1. **G is a tree**: $P(G)$ is a single point, decided directly from
   G's own max degree.
2. **G is connected but not 2-connected, and not a tree**: G is
   provably not RN (Devriendt: the only RN graphs that are not
   2-connected are paths). The bowtie graph (two triangles sharing a
   hub vertex) is the case that surfaced this: it has $t^*=2$ exactly,
   which the plain `t* <= 2 + tol` rule reads as RN=True, but the
   bowtie is not 2-connected and not a path, so it must be RN=False.
   $P(G)$ is lower-dimensional for such graphs (some subtour constraint
   is a forced equality across the *whole* polytope, not just at the
   optimum), which the closed-$t^*$ check has no way to detect.

Each module can also be run directly (`python resistance.py`, etc.) to
execute a few built-in sanity checks against known examples from the
paper. `python examples.py` reproduces the key computational claims
end-to-end, including:

- the Petersen graph is RP;
- grid graphs `P_m x P_n` are sprawling (Theorem 17);
- the `G_5(1,1,1,1,1)` construction from Theorem 4 is 1-tough but not RN,
  disproving Fiedler's conjecture that every 1-tough graph is RP.

`python verify_figure2_examples.py` checks all nine panels of Figure 1
(the paper's main examples figure; formerly labeled "Figure 2"). Panel
lettering, current revision: (A) bowtie, (B) K_{2,3}, (C) small
2-hub/3-leg banana, (D) 2-hub/4-leg banana, (E) K3-hub+legs, (F)
K4-hub+legs, (G) K5-hub+legs, (H) Petersen, (I) path family. (A
previous figure revision included a Thomassen 34-graph panel and did
not have panel (C); that panel no longer exists in the figure, and this
file has been updated to match the current panel letters.)

All nine panels (A)-(I) match their captions exactly. Panels (C) and
(F) went through a genuine back-and-forth before landing there: both
are built from the same "hub(s) connected via parallel 2-vertex legs to
a common point" family and both claim "SRN" (RN=True, RP=False). An
earlier exact spanning-tree-enumeration cross-check flagged both as
contradicting their captions -- but that check required every
individual spanning tree to have positive probability, which is
stricter than Lemma 12 actually requires (only every edge's marginal
probability needs to be positive, not every tree). Redone with the
correct edge-level criterion and cross-validated with two solvers, both
panels show a genuine positive RN margin, matching their captions; see
`verify_figure2_examples.py`'s module docstring for the full account.

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
  repo; see the proof of Theorem 5 in the paper. (Note: the current
  revision of Figure 1 no longer includes a Thomassen 34-graph panel.)

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
