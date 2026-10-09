"""
examples.py
===========

The witnesses of the containment diagram in the appendix: one graph per cell, at the smallest order known.

The diagram compares seven properties: 2-connected, traceable, 1-tough, RN, RP, sprawling, and Hamiltonian.
It has fourteen nonempty cells and draws a witness in each.
This module builds all fourteen, and running it

    python examples.py

recomputes each witness's class with `resistance.py`, recomputes the rest of its property vector where that is in reach, and compares both against what the cell asserts.

The three big witnesses are past the reach of the brute force in `sprawling.py` and `toughness.py`, so a cheap witness settles each of their remaining properties instead:

  * the fragment `Gamma` is not traceable, which this module checks by enumeration, and the two fragment lemmas of the appendix carry that to `X_37` not traceable and `X_29` not sprawling;
  * neither `X_37` nor `X_29` is 1-tough, by the vertex cut this module exhibits and checks;
  * an explicit Hamiltonian path, the decomposition the appendix describes, makes `X_29` traceable, and this module checks it;
  * RP implies 1-tough (fiedler2011, Theorem 3.4.18), which settles `T_34`.

The non-traceability of `T_34` rests on a hypotraceability search, and this module does not run it.
The output names that gap rather than passing over it.

Citation keys are those of `references.bib`.
"""

from __future__ import annotations

from typing import Dict, List, Optional, Set, Tuple

import networkx as nx


# ---------------------------------------------------------------------
# The seven properties, in the order the cells of the diagram nest.
# ---------------------------------------------------------------------

PROPS = ("2-connected", "traceable", "1-tough", "RN", "RP", "sprawling",
         "Hamiltonian")

# Above this order the subset and Hamiltonian path enumerations behind
# `toughness.py` and `sprawling.py` stop being usable. It separates the
# eleven small witnesses from the three big ones.
SMALL = 12


# ---------------------------------------------------------------------
# The small witnesses
# ---------------------------------------------------------------------

def claw() -> nx.Graph:
    """K_{1,3}, the witness outside every one of the seven properties."""
    return nx.star_graph(3)


def k13_plus_edge() -> nx.Graph:
    """
    K_{1,3} plus an edge in the independent set, a triangle with one pendant
    vertex. It is traceable, not 2-connected, and not RN.
    The appendix also names it the paw.
    """
    G = nx.complete_graph(3)
    G.add_edge(2, 3)
    return G


def path(n: int = 3) -> nx.Graph:
    """
    P_n, the path on n vertices.
    The paths are the RN graphs that are not 2-connected (devriendt2026,
    Proposition 3.7), and P_3 is the witness for that cell.
    """
    return nx.path_graph(n)


def triangle() -> nx.Graph:
    """K_3, the smallest Hamiltonian graph."""
    return nx.complete_graph(3)


def complete_bipartite(a: int, b: int) -> nx.Graph:
    """K_{a,b}, relabelled so the two parts read `a0..` and `b0..`."""
    G = nx.complete_bipartite_graph(a, b)
    return nx.relabel_nodes(
        G, {i: (f"a{i}" if i < a else f"b{i - a}") for i in range(a + b)})


def k23() -> nx.Graph:
    """K_{2,3}: RN, sprawling, and not 1-tough."""
    return complete_bipartite(2, 3)


def k23_plus_edge() -> nx.Graph:
    """K_{2,3} plus the edge joining the two vertices of the smaller part."""
    G = complete_bipartite(2, 3)
    G.add_edge("a0", "a1")
    return G


def k24() -> nx.Graph:
    """K_{2,4}: 2-connected, and nothing else."""
    return complete_bipartite(2, 4)


def toughness_family(s: List[int]) -> nx.Graph:
    """
    The graph G_t(s_1, ..., s_t): take K_{t+1} on {v0, v1, ..., vt} and
    subdivide the edge v0-vi exactly s_i times, so that v0 keeps degree t
    and the other t vertices keep the clique among themselves.

    `s` has length t. Four members of the family are witnesses:

        G_5(1,1,1,1,1)  1-tough, not RN, not traceable
        G_4(1,1,1,1)    1-tough, traceable, not RN
        G_3(1,2,1)      RP, traceable, not sprawling
        G_3(1,1,1)      RP, sprawling, not Hamiltonian

    Every member is 1-tough.
    A member with t >= 5 is neither RN nor traceable: v0 has t neighbours of
    degree 2, so it is interior to every Hamiltonian path, and t >= 5 leaves
    more such neighbours needing to be endpoints than a path has ends.
    """
    t = len(s)
    G = nx.complete_graph(t + 1)
    G = nx.relabel_nodes(G, {i: f"v{i}" for i in range(t + 1)})

    for i in range(1, t + 1):
        G.remove_edge("v0", f"v{i}")
        prev = "v0"
        for k in range(s[i - 1]):
            node = f"x{i}_{k}"
            G.add_edge(prev, node)
            prev = node
        G.add_edge(prev, f"v{i}")

    return G


# ---------------------------------------------------------------------
# The X_n family: a bipartite fragment, repeated and hubbed
# ---------------------------------------------------------------------

# Gamma, the 11-vertex bipartite fragment the two big SRN witnesses are
# built from: larger part A of size 6, smaller part B of size 5, and no
# Hamiltonian path. Its two vertices of degree 1 must be the ends of any
# Hamiltonian path, so every other vertex is interior to it. That forces
# both edges at each of its five vertices of degree 2 onto the path, and
# six of the eight edges so forced make a 6-cycle. No path contains one.
GAMMA_EDGES = [(0, 6), (0, 7), (0, 10), (1, 6), (1, 9), (2, 7), (2, 9),
               (3, 8), (3, 9), (3, 10), (4, 8), (5, 10)]
GAMMA_A = (0, 1, 2, 3, 4, 5)


def gamma_fragment() -> Tuple[nx.Graph, Set]:
    """Gamma, and its larger part A."""
    G = nx.Graph()
    G.add_nodes_from(range(11))
    G.add_edges_from(GAMMA_EDGES)
    return G, set(GAMMA_A)


def fragment_product(Gamma: nx.Graph, A: Set, n: int, m: int) -> nx.Graph:
    """
    X(Gamma, n, m): take n disjoint copies of the bipartite graph Gamma,
    whose larger part is A, add m further vertices making up a set B', and
    join every vertex of every copy of A to every vertex of B'.

    The two parts are the n copies of A, and B' together with the n copies
    of B.

    Suppose Gamma is not traceable and m + 1 < 2n. Then X(Gamma, n, m) is
    not traceable. The m vertices of B' cut a Hamiltonian path into m + 1
    pieces, and no piece straddles two copies. Fewer than 2n pieces across
    n copies leave some copy covered by a single piece, which is then a
    Hamiltonian path of Gamma.

    Copy i of a vertex v is labelled `Gi:v`, and the vertices of B' are
    labelled `y1..ym`.
    """
    return _hub(Gamma, A, n, m, swap_last=False)


def fragment_product_swapped(Gamma: nx.Graph, A: Set, n: int,
                             m: int) -> nx.Graph:
    """
    X'(Gamma, n, m): X(Gamma, n, m) with one copy of Gamma replaced by
    K_{1,2}, whose two leaves take the role of A.

    Suppose Gamma is not traceable and m = 2n - 2. The counting of
    `fragment_product` then leaves the copy of K_{1,2} covered by exactly
    one piece of any Hamiltonian path. That copy is therefore a contiguous
    interval of every Hamiltonian path, so X'(Gamma, n, 2n - 2) is not
    sprawling.
    """
    return _hub(Gamma, A, n, m, swap_last=True)


def _hub(Gamma: nx.Graph, A: Set, n: int, m: int, swap_last: bool) -> nx.Graph:
    G = nx.Graph()
    for i in range(1, n + 1):
        if swap_last and i == n:
            # K_{1,2}, written a - b - a' with a and a' in the role of A.
            G.add_edges_from([(f"G{i}:a", f"G{i}:b"), (f"G{i}:b", f"G{i}:a'")])
            boundary = [f"G{i}:a", f"G{i}:a'"]
        else:
            G.add_nodes_from(f"G{i}:{v}" for v in Gamma.nodes())
            G.add_edges_from((f"G{i}:{u}", f"G{i}:{v}") for u, v in Gamma.edges())
            boundary = [f"G{i}:{v}" for v in A]
        for v in boundary:
            for j in range(1, m + 1):
                G.add_edge(v, f"y{j}")
    return G


def x37() -> nx.Graph:
    """
    X_37 = X(Gamma, 3, 4): 37 vertices, 108 edges, parts of size 18 and 19.
    Strictly RN and not traceable.
    """
    Gamma, A = gamma_fragment()
    return fragment_product(Gamma, A, 3, 4)


def x29() -> nx.Graph:
    """
    X_29 = X'(Gamma, 3, 4): 29 vertices, 82 edges, parts of size 14 and 15.
    Strictly RN, traceable, and not sprawling.
    """
    Gamma, A = gamma_fragment()
    return fragment_product_swapped(Gamma, A, 3, 4)


# The two graphs as they were certified, for `check_construction` to match
# the constructors above against.
X37_GRAPH6 = ("dOG@EAOBAo?????@????A????@??K??H???W??k???????????A????????O???"
              "?????_???K????Q????@_???D_v`nB]E{LwZoZov`n?v`nB]?")
X29_GRAPH6 = ("\\OG@EAOBAo?????@????A????@??K??H???W??k???????????Bv`nBZov`e"
              "{LwWv`nB?")


def from_graph6(s: str) -> nx.Graph:
    """The graph a graph6 string names."""
    return nx.from_graph6_bytes(s.encode("ascii"))


def check_construction(G: nx.Graph, graph6: str) -> bool:
    """True when the constructed graph is the one the graph6 string names."""
    return nx.is_isomorphic(G, from_graph6(graph6))


# The Hamiltonian path of X_29 that the appendix describes: the two copies
# of Gamma each give up a single edge and a 9-vertex path, the copy of
# K_{1,2} gives its whole self, and the four hub vertices join the five
# pieces. Every piece ends in A but for the two ends of the whole path.
X29_HAMILTONIAN_PATH = (
    ["G1:8", "G1:4"] + ["y1"]
    + ["G1:5", "G1:10", "G1:3", "G1:9", "G1:1", "G1:6", "G1:0", "G1:7", "G1:2"]
    + ["y2"] + ["G3:a", "G3:b", "G3:a'"] + ["y3"]
    + ["G2:2", "G2:7", "G2:0", "G2:6", "G2:1", "G2:9", "G2:3", "G2:10", "G2:5"]
    + ["y4"] + ["G2:4", "G2:8"])


# ---------------------------------------------------------------------
# T_34
# ---------------------------------------------------------------------

# Thomassen's hypotraceable graph on 34 vertices and 52 edges, the smallest
# known RP graph that is not traceable. It has no short constructor, so it
# travels as the graph6 string it was certified under.
T34_GRAPH6 = ("aUYA@CPAG??A?A?B??o?O?C_?@G?@G???C??@??O??B???G???H???C_O?@G?"
              "?????C????O?O???Go??C?_????H???@@G")


def thomassen34() -> nx.Graph:
    """T_34: RP and not traceable."""
    return from_graph6(T34_GRAPH6)


# ---------------------------------------------------------------------
# Two named families that are not cells
# ---------------------------------------------------------------------

def petersen_graph() -> nx.Graph:
    """
    The Petersen graph: sprawling, RP, and not Hamiltonian. It is displaced
    from its cell of the diagram by G_3(1,1,1), on seven vertices.
    """
    return nx.petersen_graph()


def grid_graph(m: int, n: int) -> nx.Graph:
    """The grid graph P_m x P_n, which is sprawling for all m, n >= 2."""
    return nx.grid_2d_graph(m, n)


# ---------------------------------------------------------------------
# Witness checks, for the properties no enumeration can reach
# ---------------------------------------------------------------------

def is_hamiltonian_path(G: nx.Graph, walk: List) -> bool:
    """True when `walk` lists every vertex of G once and consecutively."""
    if sorted(map(str, walk)) != sorted(map(str, G.nodes())):
        return False
    return all(G.has_edge(u, v) for u, v in zip(walk, walk[1:]))


def cut_breaks_toughness(G: nx.Graph, S: Set) -> bool:
    """True when S witnesses that G is not 1-tough, i.e. c(G - S) > |S|."""
    H = G.copy()
    H.remove_nodes_from(S)
    return nx.number_connected_components(H) > len(S)


def smaller_part(G: nx.Graph) -> Set:
    """The smaller of the two parts of a connected bipartite graph."""
    left, right = nx.bipartite.sets(G)
    return left if len(left) <= len(right) else right


# ---------------------------------------------------------------------
# The fourteen cells
# ---------------------------------------------------------------------

# Each row is the cell, the witness's name, how to build it, its order, the
# class and optimum `certify/rn --exact` settles it with, and the property
# vector the cell asserts, flagged in the order of PROPS. An optimum of
# None is a graph `rn` answers by a structural shortcut, with no program
# run and so no optimum to report.
CELLS = (
    ("outside everything",
     "K_{1,3}", claw, 4, "not RN", None, "nnnnnnn"),
    ("2-connected, nothing else",
     "K_{2,4}", k24, 6, "not RN", None, "ynnnnnn"),
    ("1-tough, not RN, not traceable",
     "G_5(1,1,1,1,1)", lambda: toughness_family([1, 1, 1, 1, 1]), 11,
     "not RN", "-1/12", "ynynnnn"),
    ("1-tough, traceable, not RN",
     "G_4(1,1,1,1)", lambda: toughness_family([1, 1, 1, 1]), 9,
     "not RN", "0", "yyynnnn"),
    ("RP, not traceable",
     "T_34", thomassen34, 34, "RP", "1/34", "ynyyynn"),
    ("RP, traceable, not sprawling",
     "G_3(1,2,1)", lambda: toughness_family([1, 2, 1]), 8, "RP", "1/11",
     "yyyyynn"),
    ("RP, sprawling, not Hamiltonian",
     "G_3(1,1,1)", lambda: toughness_family([1, 1, 1]), 7, "RP", "2/19",
     "yyyyyyn"),
    ("Hamiltonian",
     "K_3", triangle, 3, "RP", "1/3", "yyyyyyy"),
    ("RN, not 2-connected",
     "P_3", lambda: path(3), 3, "SRN", "1", "nynynnn"),
    ("RN, sprawling, not 1-tough",
     "K_{2,3}", k23, 5, "SRN", "1/9", "yynynyn"),
    ("strictly RN, not traceable",
     "X_37", x37, 37, "SRN", "1/63", "ynnynnn"),
    ("strictly RN, traceable, not sprawling",
     "X_29", x29, 29, "SRN", "1/36", "yynynnn"),
    ("traceable, not 1-tough, not RN",
     "K_{2,3}+e", k23_plus_edge, 5, "not RN", "0", "yynnnnn"),
    ("traceable, not 2-connected, not RN",
     "K_{1,3}+e", k13_plus_edge, 4, "not RN", None, "nynnnnn"),
)


def asserted(flags: str) -> Dict[str, bool]:
    """The property vector a flag string such as `yynynyn` stands for."""
    return dict(zip(PROPS, [c == "y" for c in flags]))


# ---------------------------------------------------------------------
# Recomputing a property vector
# ---------------------------------------------------------------------

def _small_properties(G: nx.Graph) -> Dict[str, Tuple[Optional[bool], str]]:
    """Every property of a small witness, each by its own enumeration."""
    import sprawling
    import toughness

    H = sprawling.from_networkx(G)
    sprawl, _ = sprawling.is_sprawling(H)
    tough, _ = toughness.is_one_tough(G)
    n = G.number_of_nodes()
    return {
        "2-connected": (nx.is_biconnected(G), "direct"),
        "traceable": (len(sprawling.all_hamiltonian_paths(H)) > 0, "direct"),
        "1-tough": (tough, "direct"),
        "sprawling": (sprawl, "direct"),
        "Hamiltonian": (n >= 3 and _hamiltonian(G), "direct"),
    }


def _hamiltonian(G: nx.Graph) -> bool:
    """Brute force: a Hamiltonian path of G whose ends are adjacent."""
    import sprawling
    H = sprawling.from_networkx(G)
    return any(G.has_edge(p[0], p[-1])
               for p in sprawling.all_hamiltonian_paths(H))


def _big_properties(name: str, G: nx.Graph,
                    rp: bool) -> Dict[str, Tuple[Optional[bool], str]]:
    """
    What is known about a witness too big to enumerate, and how.

    A value of None is a property left to a computation this module does
    not run; the caller prints which one.
    """
    out: Dict[str, Tuple[Optional[bool], str]] = {
        "2-connected": (nx.is_biconnected(G), "direct"),
    }

    if name == "T_34":
        out["1-tough"] = (rp, "from RP (fiedler2011, Thm 3.4.18)")
        for p in ("traceable", "sprawling", "Hamiltonian"):
            out[p] = (None, "not run: the hypotraceability search")
        return out

    # X_37 and X_29 are bipartite with parts differing by one, and removing
    # the smaller part leaves every vertex of the larger one isolated, so
    # that part is a cut of k vertices leaving k + 1 components.
    S = smaller_part(G)
    out["1-tough"] = (not cut_breaks_toughness(G, S),
                      f"from the cut of {len(S)} vertices leaving "
                      f"{G.number_of_nodes() - len(S)} components")
    out["Hamiltonian"] = (False, "from the unequal bipartition parts")

    # The fragment lemmas carry only when Gamma is not traceable, so the
    # two properties they settle are left unsettled if that check fails.
    Gamma, _ = gamma_fragment()
    by_lemma = False if not _traceable(Gamma) else None
    how = "from Gamma not traceable, checked, and the fragment lemma"

    if name == "X_37":
        out["traceable"] = (by_lemma, how)
        out["sprawling"] = (False, "from not traceable")
    else:
        ok = is_hamiltonian_path(G, X29_HAMILTONIAN_PATH)
        out["traceable"] = (ok, "from the explicit Hamiltonian path, "
                                "checked")
        out["sprawling"] = (by_lemma, how)
    return out


def _traceable(G: nx.Graph) -> bool:
    import sprawling
    return len(sprawling.all_hamiltonian_paths(sprawling.from_networkx(G))) > 0


# ---------------------------------------------------------------------
# The driver
# ---------------------------------------------------------------------

def _flag(v: Optional[bool]) -> str:
    return "?" if v is None else ("y" if v else "n")


def check_cell(cell: str, name: str, build, n: int, klass: str,
               optimum: Optional[str], flags: str):
    """
    Build one witness, recompute what can be recomputed, and compare it
    against what its cell asserts.

    Returns (row, problems, notes): the table line, the disagreements found,
    and the properties not recomputed directly, each with the route taken.
    """
    import resistance

    G = build()
    want = asserted(flags)
    problems: List[str] = []
    notes: List[str] = []

    if G.number_of_nodes() != n:
        problems.append(f"{name}: built on {G.number_of_nodes()} vertices, "
                        f"not {n}")

    report = resistance.rn_report(G)
    got_class = report["class"]
    got_opt = report.get("exact", {}).get("optimum")

    if got_class != klass:
        problems.append(f"{name}: class {got_class}, not {klass}")
    if got_opt != optimum:
        problems.append(f"{name}: optimum {got_opt}, not {optimum}")

    have: Dict[str, Tuple[Optional[bool], str]] = {
        "RN": (bool(report["rn"]), "direct"),
        "RP": (bool(report["rp"]), "direct"),
    }
    if n <= SMALL:
        have.update(_small_properties(G))
    else:
        have.update(_big_properties(name, G, bool(report["rp"])))

    for p in PROPS:
        value, how = have[p]
        if how != "direct":
            notes.append((name, p, how))
        if value is not None and value != want[p]:
            problems.append(f"{name}: {p} is {_flag(value)}, the cell "
                            f"asserts {_flag(want[p])}")

    opt = optimum if optimum else report.get("shortcut", "shortcut")
    vec = "".join(_flag(have[p][0]) for p in PROPS)
    row = f"{name:<15} {n:>3}  {got_class:<7} {opt:<10} {vec}  {cell}"
    return row, problems, notes


def main() -> int:
    problems: List[str] = []
    notes: List[Tuple[str, str, str]] = []

    print("properties: " + " ".join(PROPS))
    print(f"{'witness':<15} {'n':>3}  {'class':<7} {'optimum':<10} "
          f"{'props':<7}  cell")
    print("-" * 78)
    for cell in CELLS:
        row, p, nt = check_cell(*cell)
        print(row)
        problems += p
        notes += nt
    print()

    built = []
    for name, G, g6 in (("X_37", x37(), X37_GRAPH6),
                        ("X_29", x29(), X29_GRAPH6)):
        if check_construction(G, g6):
            built.append(name)
        else:
            problems.append(f"{name}: the constructor has drifted off the "
                            f"certified graph")
    print("built and matched against the certified graph: "
          + ", ".join(built))

    print("not by enumeration:")
    for name, prop, how in notes:
        print(f"  {name:<6} {prop:<12} {how}")

    import sprawling
    named = []
    for label, G in (("Petersen", petersen_graph()),
                     ("P_3 x P_3", grid_graph(3, 3)),
                     ("P_2 x P_4", grid_graph(2, 4))):
        sprawl, info = sprawling.is_sprawling(sprawling.from_networkx(G))
        if sprawl:
            named.append(label)
        else:
            problems.append(f"{label}: sprawling is n, expected y ({info})")
    print("sprawling, and outside the cells: " + ", ".join(named))

    print()
    if problems:
        print(f"{len(problems)} disagreement(s):")
        for line in problems:
            print(f"  {line}")
    else:
        print(f"all {len(CELLS)} cells agree with their witnesses")
    return 1 if problems else 0


if __name__ == "__main__":
    raise SystemExit(main())
