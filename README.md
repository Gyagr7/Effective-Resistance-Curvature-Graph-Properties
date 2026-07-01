# Code for "On Some Structural Properties of Graphs with Non-Negative Resistance Curvature"

This repository contains the code used to check the RN/RP status,
sprawling property, and toughness of graphs discussed in the paper:

> G. Agrahari, C. Bibby, S. Boros, H. Garcia, F. Heiderscheidt, Z. Wang.
> *On some structural properties of graphs with non-negative resistance
> curvature.* [arXiv link / venue, once available]

## Contents

| File | Purpose |
|---|---|
| `resistance.py` | Decide whether a graph is resistance nonnegative (RN) or resistance positive (RP), via a cutting-plane LP over the spanning tree polytope (Theorem 1). |
| `sprawling.py` | Decide whether a graph is *sprawling* (Section 5), a sufficient condition for RN (Theorem 8). |
| `toughness.py` | Compute exact (vertex) toughness and check 1-toughness, by brute force. |
| `examples.py` | Graph constructions referenced in the paper (Petersen graph, grid graphs, the `G_t(s_1,...,s_t)` family from Theorem 4, etc.). |

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
is_rp, is_rn, t_star, x = resistance.resistance_positive_decision(G)

# Sprawling decision (Section 5) -- takes an adjacency-dict graph
is_sprawl, witness = sprawling.is_sprawling(sprawling.from_networkx(G))

# Toughness (used in the proof of Theorem 4)
tau, cut_set = toughness.toughness(G)
is_1_tough, witness = toughness.is_one_tough(G)
```

Each module can also be run directly (`python resistance.py`, etc.) to
execute a few built-in sanity checks against known examples from the
paper. `python examples.py` reproduces the key computational claims
end-to-end, including:

- the Petersen graph is RP (t\* = 1.8 < 2);
- grid graphs `P_m x P_n` are sprawling (Theorem 17);
- the `G_5(1,1,1,1,1)` construction from Theorem 4 is 1-tough but not RN,
  disproving Fiedler's conjecture that every 1-tough graph is RP.

## Scope and caveats

- `resistance.py` solves an LP via cutting planes and is numerically
  exact up to solver tolerance; it works comfortably on graphs with
  dozens of vertices.
- `sprawling.py` and `toughness.py` are brute-force (they enumerate
  Hamiltonian paths / vertex subsets respectively) and are only
  practical for small graphs -- exactly the regime used for the
  examples in the paper (roughly n <= 12).
- The verification that the Thomassen 34-graph is RP (Theorem 5) uses a
  hand-constructed rational weighting rather than this code; see the
  proof of Theorem 5 in the paper.

## Requirements

See `requirements.txt`. Core dependencies are `networkx` and `cvxpy`
(with a standard open-source LP solver, e.g. ECOS or SCS, which ships
with cvxpy by default).

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
