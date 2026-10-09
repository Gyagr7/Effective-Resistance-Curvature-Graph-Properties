"""
resistance.py
=============

Decide whether a graph is resistance nonnegative (RN), resistance positive (RP), or strictly resistance nonnegative (SRN), in the sense of Devriendt's discrete resistance curvature (devriendt2026 Corollary 3.9, agraharietal2026 Theorem 1).
Citation keys are those of `references.bib`.

Method
------
This module only wraps the certifier `certify/rn`: it serializes a graph, runs the certifier on it, and reads the verdict back.
The certifier decides by linear programming over the spanning tree polytope, in exact rational arithmetic.
Write P(G) for the spanning tree polytope, P(G)^o for its relative interior, and d_v(x) for the sum of x_e over the edges at v.
Theorem 1 then says

    G is RN  <=>  P(G)^o INTERSECT { x : d_v(x) <= 2 for all v }  !=  empty
    G is RP  <=>  P(G)^o INTERSECT { x : d_v(x) <  2 for all v }  !=  empty

and (devriendt2026, Remark 3.10) notes that this makes the decision a linear program over P(G), once a tolerance turns the strict inequalities into non-strict ones.
Promoting that tolerance to a scalar t to maximize gives the program `rn` runs, which is also the shape of the two programs of (guo2026lp, Theorem 2.1),

    max t   s.t.   x(E) = n - 1
                   x_e >= t                        (every edge e)
                   x(E[S]) + (|S|-1) t <= |S|-1    (S a proper nonempty set of vertices)
                   d_v(x) + s t <= 2               (every v)

with s = 1 for RP and s = 0 for RN, so the property holds exactly when the optimum is positive.
The program is only ever run on a 2-connected graph, which is how `rn` handles the relative interior and is described in the source of `certify/rn.c`: a connected graph that is not 2-connected is RN exactly when it is a path (devriendt2026, Proposition 3.7), so those graphs are answered structurally and the relative interior of P(G) is the strict system on everything that reaches the program.
A counting bound over vertex cuts then shortens the rest.
For a nonempty proper set of vertices S leaving c components in G - S, the rank inequalities give max_{v in S} d_v(x) >= 1 + (c-1)/|S| for every x in P(G), so c >= |S| + 1 rules out RP and c >= |S| + 2 rules out RN with no program at all.
The RP half is the 1-toughness of RP graphs (fiedler2011, Theorem 3.4.18) and the RN half is the toughness an RN graph can attain below 1 (garcia2026srn); `rn` hunts for the witnessing cut over a short deterministic candidate list.
The report names the route in `shortcut`, the cut in `cut_set`, `cut_size`, and `cut_components`, and the parts in `parts`.
The certifier separates the rank inequalities on demand, each as a minimum cut.
Under `--exact` it then re-solves the rows tight at the floating point optimum in rational arithmetic, both forwards for the vertex and transposed for the dual.
Those two solutions bound the optimum from each side, which certifies its sign.

Requirements
------------
Build `certify/rn` once:

    make -C certify

The build needs GLPK and GMP, packaged as `glpk` and `gmp` on Arch and as `libglpk-dev` and `libgmp-dev` on Debian.
This module builds `rn` on first use, if the binary is missing and a compiler is available.
Set the environment variable `RN_BIN` to use a binary from elsewhere.
"""

from __future__ import annotations

import json
import os
import shutil
import subprocess
from fractions import Fraction
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import networkx as nx

_HERE = Path(__file__).resolve().parent
_CERTIFY_DIR = _HERE / "certify"
_VENDORED_BIN = _CERTIFY_DIR / "rn"


class CertificationError(RuntimeError):
    """Raised when `rn` ran but did not come back with a certified answer."""


# ---------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------

def _canon_edge(u, v):
    """Return the edge {u, v} as an ordered pair, independent of how it was given."""
    return (u, v) if u <= v else (v, u)


# ---------------------------------------------------------------------
# Locating and driving the certifier
# ---------------------------------------------------------------------

def _rn_binary(build: bool = True) -> str:
    """
    Return the path to the `rn` executable.

    Check $RN_BIN first, then certify/rn, then $PATH, and finally run `make -C certify` if `build` is set.
    """
    env = os.environ.get("RN_BIN")
    if env:
        if not os.access(env, os.X_OK):
            raise FileNotFoundError(f"RN_BIN is set to {env!r}, which is not executable")
        return env

    if os.access(_VENDORED_BIN, os.X_OK):
        return str(_VENDORED_BIN)

    found = shutil.which("rn")
    if found:
        return found

    if build and (_CERTIFY_DIR / "rn.c").is_file():
        make = shutil.which("make")
        if make:
            proc = subprocess.run([make, "-s", "-C", str(_CERTIFY_DIR)],
                                  capture_output=True, text=True)
            if proc.returncode == 0 and os.access(_VENDORED_BIN, os.X_OK):
                return str(_VENDORED_BIN)
            detail = (proc.stderr or proc.stdout).strip()
            raise FileNotFoundError(
                "could not build the certifier in certify/.\n"
                f"`make -C certify` said:\n{detail}\n\n"
                "It needs GLPK and GMP: `glpk gmp` on Arch, "
                "`libglpk-dev libgmp-dev` on Debian.")

    raise FileNotFoundError(
        "the `rn` certifier was not found. Build it with `make -C certify` "
        "(needs GLPK and GMP), or point $RN_BIN at an existing binary.")


def _to_rn_json(G: nx.Graph) -> Tuple[dict, List]:
    """
    Serialize G for `rn`.

    Vertices go out as indices, so this also returns the table taking an index back to its label.
    Reading the certificate needs that table.
    """
    nodes = list(G.nodes())
    index = {v: i for i, v in enumerate(nodes)}
    payload = {
        "n": len(nodes),
        "edges": [[index[u], index[v]] for u, v in G.edges()],
    }
    return payload, nodes


def _run_rn(payload: dict, exact: bool, build: bool = True) -> dict:
    """Pipe one graph through `rn` and return its parsed report."""
    binary = _rn_binary(build=build)
    args = [binary, "--jsonl"]
    if exact:
        args.append("--exact")
    proc = subprocess.run(args, input=json.dumps(payload),
                          capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(
            f"{binary} exited {proc.returncode}: {proc.stderr.strip()}")
    out = proc.stdout.strip()
    if not out:
        raise RuntimeError(f"{binary} produced no output: {proc.stderr.strip()}")
    return json.loads(out.splitlines()[-1])


def rn_report(G: nx.Graph, exact: bool = True, build: bool = True) -> dict:
    """
    Return the certifier's full report on G.

    See `certify/rn.c` for the field list.
    Vertices and edges come back under G's own labels rather than the indices sent to `rn`.
    """
    if not nx.is_connected(G):
        raise ValueError("Graph must be connected (no spanning tree otherwise).")

    payload, nodes = _to_rn_json(G)
    report = _run_rn(payload, exact=exact, build=build)

    if not report.get("connected", True):
        raise ValueError(report.get("error", "rn reports G as disconnected"))
    if report.get("class") is None:
        raise RuntimeError(f"rn returned no verdict: {report.get('error')}")
    if exact and report.get("exact", {}).get("status") != "certified":
        raise CertificationError(
            f"rn could not certify its answer for this graph "
            f"(status={report.get('exact', {}).get('status')!r}). "
            f"Re-run with exact=False to accept the floating point verdict.")

    # Translate rn's indices back to this graph's own labels.
    report["vertices"] = [nodes[i] for i in report["vertices"]]
    report["edges"] = [_canon_edge(nodes[u], nodes[v]) for u, v in report["edges"]]
    exact_block = report.get("exact")
    witnesses = [report.get("witness"),
                 exact_block.get("witness") if isinstance(exact_block, dict) else None]
    for w in witnesses:
        if isinstance(w, dict) and w.get("zero_curvature"):
            w["zero_curvature"] = [nodes[i] for i in w["zero_curvature"]]
    return report


# ---------------------------------------------------------------------
# The decision
# ---------------------------------------------------------------------

def _witness_of(report: dict) -> Optional[dict]:
    """Return the exact witness if there is one, and the floating point one otherwise."""
    exact = report.get("exact")
    if isinstance(exact, dict) and isinstance(exact.get("witness"), dict):
        return exact["witness"]
    w = report.get("witness")
    return w if isinstance(w, dict) else None


def resistance_positive_decision(
    G: nx.Graph,
    exact: bool = True,
    build: bool = True,
    verbose: bool = True,
):
    """
    Decide the RN / RP / SRN status of a connected graph G.

    This reads the verdict of `rn_report` back as plain Python.
    A graph that is not 2-connected is settled without a program: a path comes back SRN (RP for at most one edge), and every other graph with a cut vertex comes back not RN.
    So is a graph carrying a vertex cut S whose removal leaves |S| + 2 components or more, which comes back not RN.

    Parameters
    ----------
    exact : certify the sign of the optimum in rational arithmetic.
        Raise `CertificationError` if no certificate is found.
        With `exact=False` the floating point verdict is returned unchecked.
        That runs roughly 20% faster and is unsound on edge cases.
    build : build `certify/rn` if it is missing, as described in `_rn_binary`.
    verbose : print the verdict and the certificate.

    Returns
    -------
    rp : bool
    rn : bool
    t_star : float or None
        max_v d_v(x) at the returned witness x.
    x_dict : dict or None
        edge -> x_e for the witness point, labelled by G.

    Raises
    ------
    ValueError : G is not connected.
    CertificationError : `exact` was set and no certificate was obtained.
    FileNotFoundError : the certifier was missing and could not be built.
    """
    report = rn_report(G, exact=exact, build=build)

    rp = bool(report["rp"])
    rn = bool(report["rn"])
    witness = _witness_of(report)

    t_star: Optional[float] = None
    x_dict: Optional[Dict] = None
    if witness is not None:
        # Exact mode sends rationals as strings such as "11/6".
        # Fraction reads those and the floats of the inexact mode alike.
        degrees = [Fraction(str(d)) for d in witness["degree"]]
        t_star = float(max(degrees)) if degrees else 0.0
        x_dict = {e: float(Fraction(str(xe)))
                  for e, xe in zip(report["edges"], witness["point"])}

    if verbose:
        n, m = report["n"], report["m"]
        parts = report.get("parts")
        print(f"n={n}, m={m}, 2-connected={report.get('two_connected')}"
              + (f", parts={parts[0]}+{parts[1]}" if parts else ""))
        if report.get("shortcut"):
            where = report.get("cut_vertex")
            at = f" at {where}" if where is not None else ""
            if report.get("cut_size") is not None:
                at = (f": a cut of {report['cut_size']} vertices leaving "
                      f"{report['cut_components']} components")
            print(f"shortcut: {report['shortcut']}{at}")
        print(f"separation: {report['rank_cuts']} rank cuts over "
              f"{report['rounds']} rounds, {report['seconds']:.4f}s")
        if exact:
            opt = report["exact"].get("optimum")
            lo = report["exact"].get("bound_lower")
            hi = report["exact"].get("bound_upper")
            why = report["exact"].get("reason")
            if opt is not None:
                print(f"certified: optimum = {opt}")
            elif lo is not None or hi is not None:
                print(f"certified: optimum in [{lo}, {hi}]")
            else:
                print(f"certified: {why or 'structurally'}")
        print(f"class={report['class']}  ->  RN={rn}, RP={rp}")
        if t_star is not None:
            print(f"witness: max_v d_v(x) = {t_star:.12f}")

    return rp, rn, t_star, x_dict


def _clique_with_legs() -> nx.Graph:
    """
    K_4 on 0..3, with a leg of length three from each of its vertices to a
    common vertex 12. The optimum comes out exactly 0, so it is RN only if
    the certification is exact.
    """
    G = nx.Graph()
    G.add_edges_from([(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)])
    for hub, (a, b) in zip(range(4), ((4, 5), (6, 7), (8, 9), (10, 11))):
        G.add_edges_from([(hub, a), (a, b), (b, 12)])
    return G


# Each check is a graph and the class the certifier should return. The
# optimum is reported alongside, and is None where a structural shortcut
# answered the graph with no program run.
CHECKS = (
    ("K_3", nx.complete_graph(3), "RP"),
    ("P_3", nx.path_graph(3), "SRN"),
    ("K_{1,3}", nx.star_graph(3), "not RN"),
    ("K_{2,3}", nx.complete_bipartite_graph(2, 3), "SRN"),
    ("Petersen", nx.petersen_graph(), "RP"),
    ("K_4 with legs", _clique_with_legs(), "not RN"),
)


def main() -> int:
    """Run the checks, and return 1 if any graph disagrees with its row."""
    problems = []
    print(f"{'graph':<14} {'n':>3} {'class':<7} {'optimum':<8} route")
    for name, G, want in CHECKS:
        report = rn_report(G)
        got = report["class"]
        optimum = report.get("exact", {}).get("optimum")
        route = report.get("shortcut") or "program"
        print(f"{name:<14} {report['n']:>3} {got:<7} "
              f"{str(optimum if optimum is not None else '-'):<8} {route}")
        if got != want:
            problems.append(f"{name}: class is {got!r}, expected {want!r}")

    print()
    if problems:
        print(f"{len(problems)} disagreement(s):")
        for line in problems:
            print(f"  {line}")
        return 1
    print(f"all {len(CHECKS)} graphs agree with their known class")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
