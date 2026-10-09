# Code for "On Some Structural Properties of Graphs with Non-Negative Resistance Curvature"

This repository contains the code used to check the RN/RP status, sprawling property, and toughness of graphs discussed in the paper:

> G. Agrahari, C. Bibby, S. Boros, H. J. Garcia, F. Heidercheidt, Z. Wang.
> *On some structural properties of graphs with non-negative resistance curvature.*
> arXiv:2607.13169, <https://doi.org/10.48550/arXiv.2607.13169>

Citation keys below are those of `references.bib`.

## AI disclosure

Claude Opus 5 (Anthropic) wrote `certify/rn.c` and `certify/sweep.py` in their entirety, including the formulation of the linear program, the structural shortcuts and their lemmas, the separation routine, the exact rational certification, the I/O, and the comments.
It also drafted parts of the prose and the docstrings elsewhere in this repository, and the current form of `resistance.py` and `examples.py`.
Hailey Jay Garcia directed and supervised that work, reviewed the result, and is responsible for the mathematics it rests on and for the answers it reports.
The mathematical content of the paper is the authors' own.

## Contents

| File | Purpose |
|---|---|
| `certify/rn.c` | The certifier: decide RN, RP, or SRN by linear programming over the spanning tree polytope, with optional exact rational certification of the answer. Written in C against GLPK and GMP. |
| `certify/sweep.py` | Differential test for the certifier: run two builds of `certify/rn` over every connected graph of a given order, from nauty's `geng`, and compare their verdicts. |
| `resistance.py` | Python wrapper around `certify/rn`. Serialize a graph, run the certifier, and read the verdict and witness back under the graph's own labels. |
| `sprawling.py` | Decide whether a graph is *sprawling*, a sufficient condition for RN. Every "sprawling" verdict comes with an explicit, independently re-verified witness collection. |
| `toughness.py` | Compute exact (vertex) toughness and check 1-toughness, by brute force. |
| `examples.py` | The witnesses of the containment diagram in the appendix, one per cell. Run it to rebuild all fourteen and check each against the properties its cell asserts. |
| `references.bib` | The bibliography cited by key from the source files. |

## Installation

```bash
pip install -r requirements.txt
make -C certify
```

The build needs GLPK and GMP, packaged as `glpk` and `gmp` on Arch and as `libglpk-dev` and `libgmp-dev` on Debian.
`resistance.py` runs that build itself on first use, if the binary is missing and a compiler is available.
Set the environment variable `RN_BIN` to use a binary from elsewhere.

`make -C certify check` runs the certifier over `certify/examples.json`, which holds the same fourteen graphs as `examples.py`.
`certify/sweep.py` needs nauty's `geng` on PATH as well (`nauty` on Arch, `nauty` on Debian, where the binary is named `nauty-geng`), and takes a second build of `rn` to compare against:

```bash
./certify/sweep.py --new certify/rn --old /path/to/other/rn 4 8
```

## Usage

```python
import networkx as nx
import resistance, sprawling, toughness

G = nx.petersen_graph()

# RN / RP decision
rp, rn, t_star, x = resistance.resistance_positive_decision(G, verbose=False)
# x is the witness point in the spanning tree polytope, keyed by edge.
# t_star is max_v d_v(x) at that witness, and is None when G is not RN.
# Pass exact=False to accept the floating point verdict without a certificate.

# The certifier's full report, including the exact certificate
report = resistance.rn_report(G)

# Sprawling decision, which takes an adjacency-dict graph
is_sprawl, info = sprawling.is_sprawling(sprawling.from_networkx(G))
# info is the explicit witness collection S if is_sprawl=True,
# or the specific failing condition/set if False

# Toughness
tau, cut_set = toughness.toughness(G)
is_1_tough, witness = toughness.is_one_tough(G)

# The witnesses of the containment diagram, by name
import examples
G = examples.x37()                      # X(Gamma, 3, 4): strictly RN, not traceable
G = examples.toughness_family([1, 2, 1])  # G_3(1,2,1): RP, traceable, not sprawling
```

### How `resistance.py` decides RN / RP

Write $P(G)$ for the spanning tree polytope, $P(G)^\circ$ for its relative interior, and $d_v(x)$ for the sum of $x_e$ over the edges at $v$.
Theorem 1 [[Corollary 3.9, Dev26]](#devriendt2026), [[ABB+26]](#agraharietal2026), says

```math
\begin{aligned}
G \text{ is RN} &\iff P(G)^\circ \cap \{\, x : d_v(x) \le 2 \text{ for all } v \,\} \ne \emptyset, \\
G \text{ is RP} &\iff P(G)^\circ \cap \{\, x : d_v(x) < 2 \text{ for all } v \,\} \ne \emptyset.
\end{aligned}
```

[[Remark 3.10, Dev26]](#devriendt2026) observes that this makes the decision a linear program over $P(G)$, at the cost of a tolerance to turn the strict inequalities into non-strict ones, and asks for the complexity of the decision problem in Question 3.11.
The certifier takes that route with the tolerance promoted to a scalar $t$ to maximize, which is also the shape of the two programs in [[Theorem 2.1, GSY26]](#guo2026lp):

```math
\begin{aligned}
\max\ t \quad \text{s.t.} \quad
  & x(E) = n - 1 \\
  & x_e \ge t && (\text{every edge } e) \\
  & x(E[S]) + (\lvert S\rvert - 1)\,t \le \lvert S\rvert - 1 && (S \text{ a proper nonempty set of vertices}) \\
  & d_v(x) + s\,t \le 2 && (\text{every } v)
\end{aligned}
```

with $s = 1$ for RP and $s = 0$ for RN, so the property holds exactly when the optimum is positive.
A connected graph that is not 2-connected is RN exactly when it is a path [[Proposition 3.7, Dev26]](#devriendt2026), which follows there from the cut-vertex curvature computation of its Example 3.6, so a graph with a cut vertex is answered structurally: a path is SRN (RP when it has at most one edge) with witness $x = 1$, and nothing else with a cut vertex is RN.
A counting bound over vertex cuts settles more of them.
For a nonempty proper set of vertices $S$ leaving $c$ components in $G - S$, the rank inequalities give

```math
\max_{v \in S} d_v(x) \ge 1 + \frac{c - 1}{\lvert S\rvert} \qquad \text{for every } x \in P(G),
```

so $c \ge \lvert S\rvert + 1$ rules out RP and $c \ge \lvert S\rvert + 2$ rules out RN, the latter with no program at all.
Both halves are known: the RP half is the 1-toughness of resistance positive graphs [[Theorem 3.4.18, Fie11]](#fiedler2011), and the RN half is the toughness an RN graph can attain below 1 [[Gar26]](#garcia2026srn).
Taking $S$ to be either part of a bipartite graph is the special case where all the components are single vertices, so parts differing by one rule out RP and parts differing by two or more rule out RN, which is [[Proposition 3.5, Dev26]](#devriendt2026).
Finding the worst cut is NP-hard, but a cut is a certificate by itself, so `rn` tries a short deterministic list of candidates (both colour classes when bipartite, the complement of a greedy independent set by ascending and by descending degree, and $N(v)$ for every $v$) and is sound wherever one fires.
Counted over the 194066 two-connected graphs on 9 vertices, that list answers 8.1% of the not-RN graphs outright and rules out RP for 94.0% of the graphs that are not RP, with no false positives; the graphs it answers are cheap ones, so on such an ensemble it saves about 2% of the floating point runtime, and about a quarter of the runtime under `--exact`, where a skipped program is a skipped rational solve.
It pays on larger graphs with a big independent set: $K_{3,120}$ plus an edge in the small part, for instance, is answered instantly in place of a 0.18s program.
The `cut-tight` route catches every SRN graph: [[Gar26]](#garcia2026srn) shows that a 2-connected strictly RN graph is bipartite with parts differing by one, so its smaller part is such a cut.
That classification also implies that a 1-tough RN graph is RP, which would let the RN program be skipped in place of the RP one; `rn` does not do that, since the RN program is the second program only where the RP optimum came out exactly 0, on 953 of the 194066 two-connected graphs on 9 vertices, and everywhere else it runs it is the only program `rn` ran at all.

The report names the route taken in its `shortcut` field (`path`, `cut-vertex`, `cut-gap`, `cut-tight`, or `null` when the programs ran), the cut vertex in `cut_vertex`, the witnessing cut in `cut_set`, `cut_size`, and `cut_components`, and the part sizes in `parts`.
On what reaches the program, 2-connected graphs, $P(G)$ carries the single equality $x(E) = n - 1$ and the relative interior really is the strict system above.
The rank inequalities are separated on demand, each as a minimum cut.

Under `exact=True`, the default, the certifier then re-solves the rows tight at the floating point optimum in rational arithmetic, forwards for the primal vertex and transposed for the dual bound.
Those two solutions bound the optimum from each side, which certifies its sign.
`resistance.py` raises `CertificationError` when no certificate is obtained, rather than returning an unchecked verdict.
There are no solver tolerances to set.

### The containment diagram

The appendix draws a containment diagram of seven properties: 2-connected, traceable, 1-tough, RN, RP, sprawling, and Hamiltonian.
It has fourteen nonempty cells and draws a witness in each, at the smallest order known.
`examples.py` builds all fourteen.
Running it recomputes each witness's class and, where that is in reach, the rest of its property vector, then compares both against what the cell asserts.
It takes under a second and exits nonzero on any disagreement.

| Cell | Witness | $n$ | Class | Optimum |
|---|---|---|---|---|
| outside everything | $K_{1,3}$ | 4 | not RN | shortcut |
| 2-connected, nothing else | $K_{2,4}$ | 6 | not RN | shortcut |
| 1-tough, not RN, not traceable | $G_5(1,1,1,1,1)$ | 11 | not RN | $-1/12$ |
| 1-tough, traceable, not RN | $G_4(1,1,1,1)$ | 9 | not RN | $0$ |
| RP, not traceable | $T_{34}$ | 34 | RP | $1/34$ |
| RP, traceable, not sprawling | $G_3(1,2,1)$ | 8 | RP | $1/11$ |
| RP, sprawling, not Hamiltonian | $G_3(1,1,1)$ | 7 | RP | $2/19$ |
| Hamiltonian | $K_3$ | 3 | RP | $1/3$ |
| RN, not 2-connected | $P_3$ | 3 | SRN | $1$ |
| RN, sprawling, not 1-tough | $K_{2,3}$ | 5 | SRN | $1/9$ |
| strictly RN, not traceable | $X_{37}$ | 37 | SRN | $1/63$ |
| strictly RN, traceable, not sprawling | $X_{29}$ | 29 | SRN | $1/36$ |
| traceable, not 1-tough, not RN | $K_{2,3}+e$ | 5 | not RN | $0$ |
| traceable, not 2-connected, not RN | $K_{1,3}+e$ | 4 | not RN | shortcut |

The optimum is the slack $t$ of the program above.
It is positive exactly when the graph is RN, and positive with the degree rows tightened exactly when the graph is RP.
The rows marked `shortcut` are answered structurally, with no program run and so no optimum to report, and `rn` gives the route in its `shortcut` field.
`certify/examples.json` holds the same graphs, so `make -C certify check` certifies the whole column directly.

The two $X_n$ come from `fragment_product` and `fragment_product_swapped`, the general construction $X(\Gamma,n,m)$ of the appendix.
Take $n$ copies of a bipartite fragment $\Gamma$, and join every vertex of every copy of the larger part to each of $m$ further vertices.
For $X'$, swap one copy for $K_{1,2}$.
Both witnesses take $\Gamma$ to be the $11$-vertex fragment, with $n = 3$ and $m = 4$.
`examples.py` checks each constructed graph against the graph6 string it was certified under, so a constructor cannot drift unnoticed.

The three big witnesses are past the reach of the brute force in `sprawling.py` and `toughness.py`, so a cheap witness settles each of their remaining properties instead.
The output names the route taken for every one:

- the fragment $\Gamma$ is not traceable, which `examples.py` checks by enumeration, and the two fragment lemmas carry that to $X_{37}$ not traceable and $X_{29}$ not sprawling;
- neither $X_{37}$ nor $X_{29}$ is 1-tough, by the vertex cut `examples.py` exhibits and checks: both are bipartite with parts differing by one, so removing the smaller part leaves one component per vertex of the larger;
- an explicit Hamiltonian path, the five-piece decomposition the appendix describes, makes $X_{29}$ traceable, and `examples.py` checks it;
- RP implies 1-tough [[Theorem 3.4.18, Fie11]](#fiedler2011), which settles $T_{34}$.

The non-traceability of $T_{34}$ rests on a hypotraceability search, and `examples.py` does not run it.
The output names that gap rather than passing over it.

`examples.py` also keeps the Petersen graph and the grid graphs $P_m \times P_n$, which no longer hold a cell, and checks that each is sprawling.

Each module can also be run directly (`python resistance.py`, etc.) to execute a few built-in sanity checks against known examples.

## Scope and caveats

- `resistance.py` is exact by default: the sign of the optimum is certified in rational arithmetic, so a returned verdict does not depend on solver precision.
  It works comfortably on graphs with dozens of vertices.
- `sprawling.py` and `toughness.py` are brute force, since they enumerate Hamiltonian paths and vertex subsets respectively, so they are only practical for small graphs, roughly $n \le 12$.
  The three big witnesses fall outside it, and are handled as above.
- `examples.py` builds the whole $G_t(s_1,\dots,s_t)$ family from one constructor, `toughness_family`, which takes the per-branch lengths and labels the vertices `v0..vt` and `xi_k`.
  The advisor's `build_minimal_tough_graph` is gone.
  It built the same graph, a hub joined to a clique by subdivided spokes, with equal branch lengths and integer labels; call `toughness_family([l] * n)` for it.
- The three big witnesses travel as graph6 strings or as constructors.
  `certify/rn --exact` certifies each of them in rational arithmetic in under 40 ms, so none needs a hand-built weighting.

## Requirements

`networkx` (see `requirements.txt`), plus GLPK and GMP and a C compiler for `certify/rn`.

## References

Keyed as in `references.bib`, which the source files cite by the same keys.
The bracketed labels are the alpha-style ones the citations above use.

#### agraharietal2026

**[ABB+26]**
G. Agrahari, C. Bibby, S. Boros, H. J. Garcia, F. Heidercheidt, Z. Wang.
*On some structural properties of graphs with non-negative resistance curvature.*
[arXiv:2607.13169](https://arxiv.org/abs/2607.13169), [doi:10.48550/arXiv.2607.13169](https://doi.org/10.48550/arXiv.2607.13169)

#### devriendt2026

**[Dev26]**
K. Devriendt.
*Graphs with nonnegative resistance curvature.*
Annals of Combinatorics 30 (2026), 415-438.
[doi:10.1007/s00026-025-00774-x](https://doi.org/10.1007/s00026-025-00774-x), [arXiv:2410.07756](https://arxiv.org/abs/2410.07756)

#### fiedler2011

**[Fie11]**
M. Fiedler.
*Matrices and Graphs in Geometry.*
Encyclopedia of Mathematics and its Applications 139, Cambridge University Press, 2011.

#### garcia2026srn

**[Gar26]**
H. J. Garcia.
*Classification of strictly resistance nonnegative graphs.*
[arXiv:2609.23394](https://arxiv.org/abs/2609.23394), [doi:10.48550/arXiv.2609.23394](https://doi.org/10.48550/arXiv.2609.23394)

#### guo2026lp

**[GSY26]**
M. Guo, W. Sun, Y. Yang.
*Resistance curvature: recognition, polyhedral structure, and graph products.*
[arXiv:2608.20778](https://arxiv.org/abs/2608.20778), [doi:10.48550/arXiv.2608.20778](https://doi.org/10.48550/arXiv.2608.20778)

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

## License

Copyright (C) 2026 Gyaneshwar Agrahari and Hailey Jay Garcia.

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
See [`LICENSE`](LICENSE) for the full text.
