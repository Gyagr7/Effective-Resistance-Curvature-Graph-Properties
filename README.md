# Code for "On Some Structural Properties of Graphs with Non-Negative Resistance Curvature"

This repository contains the code used to check the RN/RP status, sprawling property, and toughness of graphs discussed in the paper:

> G. Agrahari, C. Bibby, S. Boros, H. J. Garcia, F. Heidercheidt, Z. Wang.
> *On some structural properties of graphs with non-negative resistance curvature.*
> arXiv:2607.13169, <https://doi.org/10.48550/arXiv.2607.13169>

Citation keys below are those of `references.bib`.

## AI disclosure

Claude Opus 5 (Anthropic) wrote `certify/rn.c` and `certify/sweep.py` in their entirety, including the formulation of the linear program, the structural shortcuts and their lemmas, the separation routine, the exact rational certification, the I/O, and the comments.
It also drafted parts of the prose and the docstrings elsewhere in this repository, and the current form of `resistance.py`.
Hailey Jay Garcia directed and supervised that work, reviewed the result, and is responsible for the mathematics it rests on and for the answers it reports.
The mathematical content of the paper is the authors' own.

## Contents

| File | Purpose |
|---|---|
| `certify/rn.c` | The certifier: decide RN, RP, or SRN by linear programming over the spanning tree polytope, with optional exact rational certification of the answer. Written in C against GLPK and GMP. |
| `certify/sweep.py` | Differential test for the certifier: run two builds of `certify/rn` over every connected graph of a given order, from nauty's `geng`, and compare their verdicts. |
| `resistance.py` | Python wrapper around `certify/rn`. Serialize a graph, run the certifier, and read the verdict and witness back under the graph's own labels. |
| `sprawling.py` | Decide whether a graph is *sprawling* (Section 5), a sufficient condition for RN (Theorem 8). Every "sprawling" verdict comes with an explicit, independently re-verified witness collection. |
| `toughness.py` | Compute exact (vertex) toughness and check 1-toughness, by brute force. |
| `examples.py` | Graph constructions referenced in the paper (Petersen graph, grid graphs, the `G_t(s_1,...,s_t)` family from Theorem 4, etc.). |
| `verify_figure2_examples.py` | Run all three checkers against every example graph in Figure 1 (the paper's main examples figure), labeled by panel letter. |
| `references.bib` | The bibliography cited by key from the source files. |

## Installation

```bash
pip install -r requirements.txt
make -C certify
```

The build needs GLPK and GMP, packaged as `glpk` and `gmp` on Arch and as `libglpk-dev` and `libgmp-dev` on Debian.
`resistance.py` runs that build itself on first use, if the binary is missing and a compiler is available.
Set the environment variable `RN_BIN` to use a binary from elsewhere.

## Usage

```python
import networkx as nx
import resistance, sprawling, toughness

G = nx.petersen_graph()

# RN / RP decision (Theorem 1)
rp, rn, t_star, x = resistance.resistance_positive_decision(G, verbose=False)
# x is the witness point in the spanning tree polytope, keyed by edge.
# t_star is max_v d_v(x) at that witness, and is None when G is not RN.
# Pass exact=False to accept the floating point verdict without a certificate.

# The certifier's full report, including the exact certificate
report = resistance.rn_report(G)

# Sprawling decision (Section 5), which takes an adjacency-dict graph
is_sprawl, info = sprawling.is_sprawling(sprawling.from_networkx(G))
# info is the explicit witness collection S if is_sprawl=True,
# or the specific failing condition/set if False

# Toughness (used in the proof of Theorem 4)
tau, cut_set = toughness.toughness(G)
is_1_tough, witness = toughness.is_one_tough(G)
```

### How `resistance.py` decides RN / RP

Write P(G) for the spanning tree polytope, P(G)^o for its relative interior, and d_v(x) for the sum of x_e over the edges at v.
Theorem 1 (devriendt2025, agraharietal2026) says

```
    G is RN  <=>  P(G)^o INTERSECT { x : d_v(x) <= 2 for all v }  !=  empty
    G is RP  <=>  P(G)^o INTERSECT { x : d_v(x) <  2 for all v }  !=  empty
```

The certifier decides each by the linear program of (guo2026lp, Theorem 2.1),

```
    max t   s.t.   x(E) = n - 1
                   x_e >= t                        (every edge e)
                   x(E[S]) + (|S|-1) t <= |S|-1    (S a proper nonempty set of vertices)
                   d_v(x) + s t <= 2               (every v)
```

with s = 1 for RP and s = 0 for RN, so the property holds exactly when the optimum is positive.
The program is only ever run on a 2-connected graph, which is this implementation's own departure from (guo2026lp) and is how the relative interior is handled.
A connected graph that is not 2-connected is RN exactly when it is a path, by a short lemma proved at the head of `certify/rn.c`, so a graph with a cut vertex is answered structurally: a path is SRN (RP when it has at most one edge) with witness x = 1, and nothing else with a cut vertex is RN.
A second counting bound settles more of them.
In a bipartite graph with parts A, B and |A| <= |B|, every edge has exactly one end in A, so sum over v in A of d_v(x) equals x(E) = n - 1 for every x in P(G), and hence max_v d_v(x) >= (n-1)/|A|.
A gap |B| - |A| >= 2 therefore puts every point of P(G) above 2 somewhere and G is not RN, again with no program at all; a gap of exactly 1 rules out RP, so only the RN program is run.
Since the parts differ in parity with n, the first case is the even orders and the second the odd ones.
The report names the route taken in its `shortcut` field (`path`, `cut-vertex`, `bipartite-gap`, `bipartite-unbalanced`, or `null` when the programs ran), the cut vertex in `cut_vertex`, and the part sizes in `parts`.
On what reaches the program, 2-connected graphs, P(G) carries the single equality x(E) = n - 1 and the relative interior really is the strict system above.
The rank inequalities are separated on demand, each as a minimum cut.

Under `exact=True`, the default, the certifier then re-solves the rows tight at the floating point optimum in rational arithmetic, forwards for the primal vertex and transposed for the dual bound.
Those two solutions bound the optimum from each side, which certifies its sign.
`resistance.py` raises `CertificationError` when no certificate is obtained, rather than returning an unchecked verdict.
There are no solver tolerances to set.

Each module can also be run directly (`python resistance.py`, etc.) to execute a few built-in sanity checks against known examples from the paper.
`python examples.py` reproduces the key computational claims end-to-end, including:

- the Petersen graph is RP;
- grid graphs `P_m x P_n` are sprawling (Theorem 17);
- the `G_5(1,1,1,1,1)` construction from Theorem 4 is 1-tough but not RN, disproving Fiedler's conjecture that every 1-tough graph is RP.

### Figure 1 panels

`python verify_figure2_examples.py` checks all nine panels of Figure 1 (the paper's main examples figure; formerly labeled "Figure 2").
Panel lettering, current revision: (A) bowtie, (B) K_{2,3}, (C) small 2-hub/3-leg banana, (D) 2-hub/4-leg banana, (E) K3-hub+legs, (F) K4-hub+legs, (G) K5-hub+legs, (H) Petersen, (I) path family.
A previous figure revision included a Thomassen 34-graph panel and did not have panel (C).

Seven panels match their captions.
**Panels (C) and (F) do not, under the current certifier.**
Both are built from the same "hub(s) connected via parallel legs to a common point" family, and both captions claim SRN (RN=True, RP=False).
The exact certifier in `certify/rn.c` reports both as not RN.
Earlier floating point runs under cvxpy with the SCS solver reported a positive RN margin for both, which is where the captions come from.
The captured output block at the bottom of `verify_figure2_examples.py` is from those earlier runs, and its summary line still claims all nine panels agree.
Resolve this before the figure captions are taken as verified.

## Scope and caveats

- `resistance.py` is exact by default: the sign of the optimum is certified in rational arithmetic, so a returned verdict does not depend on solver precision.
  It works comfortably on graphs with dozens of vertices.
- `sprawling.py` and `toughness.py` are brute force, since they enumerate Hamiltonian paths and vertex subsets respectively, so they are only practical for small graphs.
  That is the regime used for the examples in the paper, roughly n <= 12.
- Both `toughness_family` and `build_minimal_tough_graph` in `examples.py` build the same underlying construction from Theorem 4 / Lemma 15, a hub connected to a clique via subdivided spokes.
  The former takes per-branch lengths and uses string labels, the latter takes equal branch lengths and uses integer labels.
  Both are kept as entry points, since both conventions have been used across the project.
- The verification that the Thomassen 34-graph is RP (Theorem 5) uses a hand-constructed rational weighting rather than any code in this repository.
  See the proof of Theorem 5 in the paper.
  The current revision of Figure 1 no longer includes a Thomassen 34-graph panel.

## Requirements

`networkx` (see `requirements.txt`), plus GLPK and GMP and a C compiler for `certify/rn`.
The earlier cvxpy and SCS dependency is gone from `resistance.py`.
`verify_figure2_examples.py` still names cvxpy in its comments, which describe how its captured output was produced.

## Citation

If you use this code, please cite the paper:

```bibtex
@misc{agraharietal2026,
  title         = {On some structural properties of graphs with non-negative resistance curvature},
  author        = {Agrahari, Gyaneshwar and Bibby, Christin and Boros, Sean and Garcia, Hailey Jay and Heidercheidt, Fernando and Wang, Zhiyu},
  year          = {2026},
  eprint        = {2607.13169},
  archivePrefix = {arXiv},
  primaryClass  = {math.CO},
  doi           = {10.48550/arXiv.2607.13169}
}
```
