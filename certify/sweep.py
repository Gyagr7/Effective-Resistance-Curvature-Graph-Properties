#!/usr/bin/env python3
"""
certify/sweep.py -- compare two builds of `rn` over every connected graph of a given order.

Graphs come from nauty's `geng -c`, and the two binaries are run on the same graph6 input.
The comparison is on the verdict fields only (`class`, `rn`, `rp`, `srn`), since the
reported optima and separation counts legitimately differ between builds that take
different routes to the same answer.

    ./sweep.py --new certify/rn --old /path/to/other/rn 4 9
    ./sweep.py --new certify/rn --old build/rn --exact 4 8

Written to check the structural shortcut of `rn.c` (a graph that is not 2-connected is
RN exactly when it is a path) against the block decomposition it replaced.
"""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
import time

FIELDS = ("class", "rn", "rp", "srn")


def geng(n: int, extra: list[str]) -> list[str]:
    exe = shutil.which("geng") or shutil.which("nauty-geng")
    if exe is None:
        sys.exit("sweep: nauty's geng is not on PATH")
    out = subprocess.run([exe, "-cq", *extra, str(n)],
                         capture_output=True, text=True, check=True)
    return [ln for ln in out.stdout.splitlines() if ln.strip()]


def run(binary: str, g6: list[str], exact: bool) -> list[dict]:
    cmd = [binary, "--jsonl", "--no-echo", "--no-witness"]
    if exact:
        cmd.append("--exact")
    out = subprocess.run(cmd, input="\n".join(g6) + "\n",
                         capture_output=True, text=True, check=True)
    return [json.loads(ln) for ln in out.stdout.splitlines() if ln.strip()]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("lo", type=int)
    ap.add_argument("hi", type=int, nargs="?")
    ap.add_argument("--new", default="certify/rn")
    ap.add_argument("--old", required=True)
    ap.add_argument("--exact", action="store_true")
    ap.add_argument("--geng", default="", help="extra geng arguments, space separated")
    a = ap.parse_args()
    hi = a.hi if a.hi is not None else a.lo

    bad = 0
    for n in range(a.lo, hi + 1):
        g6 = geng(n, a.geng.split())
        t0 = time.monotonic()
        new = run(a.new, g6, a.exact)
        t1 = time.monotonic()
        old = run(a.old, g6, a.exact)
        t2 = time.monotonic()
        if len(new) != len(g6) or len(old) != len(g6):
            sys.exit(f"sweep: n={n}: {len(g6)} graphs in, {len(new)}/{len(old)} out")

        mism = 0
        for s, rn_new, rn_old in zip(g6, new, old):
            if any(rn_new.get(k) != rn_old.get(k) for k in FIELDS):
                mism += 1
                bad += 1
                if mism <= 5:
                    print(f"  MISMATCH {s}: "
                          f"new={ {k: rn_new.get(k) for k in FIELDS} } "
                          f"old={ {k: rn_old.get(k) for k in FIELDS} }")
        cls = {}
        for r in new:
            cls[r["class"]] = cls.get(r["class"], 0) + 1
        short = sum(1 for r in new if r.get("shortcut"))
        print(f"n={n}: {len(g6)} graphs, {mism} mismatches, "
              f"new {t1 - t0:.2f}s, old {t2 - t1:.2f}s, "
              f"shortcut {short}, {cls}")

    print("all verdicts agree" if not bad else f"{bad} mismatches")
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main())
