#!/usr/bin/env python3
"""Cyclomatic-complexity ratchet for pyrestoolbox.

Fails when a function's McCabe complexity exceeds MAX_COMPLEXITY unless that
function is already recorded in the baseline at its recorded value or lower.
The point is to stop new branching accreting, not to force a rewrite of the
functions that are already large - those are listed in the baseline and can only
get simpler.

Usage:
    python3 .github/scripts/check_complexity.py            # check
    python3 .github/scripts/check_complexity.py --update   # rewrite the baseline
"""

import json
import subprocess
import sys
from pathlib import Path

MAX_COMPLEXITY = 15
PACKAGE = "pyrestoolbox"
BASELINE = Path(__file__).resolve().parent.parent / "complexity-baseline.json"
SKIP_DIRS = ("/tests/", "\\tests\\")


def measure():
    """Return {"path::function": complexity} for every function in the package."""
    out = subprocess.run(
        [sys.executable, "-m", "radon", "cc", "-j", PACKAGE],
        capture_output=True, text=True, check=True,
    ).stdout
    scores = {}
    for path, blocks in json.loads(out).items():
        if any(d in "/" + path.replace("\\", "/") + "/" for d in ("/tests/",)):
            continue
        if isinstance(blocks, dict):  # radon reports parse errors as a dict
            continue
        for b in blocks:
            name = f"{b['classname']}.{b['name']}" if b.get("classname") else b["name"]
            key = f"{path.replace(chr(92), '/')}::{name}"
            scores[key] = max(scores.get(key, 0), b["complexity"])
    return scores


def main():
    scores = measure()
    over = {k: v for k, v in scores.items() if v > MAX_COMPLEXITY}

    if "--update" in sys.argv:
        BASELINE.write_text(json.dumps(dict(sorted(over.items())), indent=2) + "\n")
        print(f"Wrote {len(over)} entries to {BASELINE.name}")
        return 0

    baseline = json.loads(BASELINE.read_text()) if BASELINE.exists() else {}

    new = {k: v for k, v in over.items() if k not in baseline}
    worse = {k: (baseline[k], v) for k, v in over.items() if k in baseline and v > baseline[k]}
    better = {k: (baseline[k], scores.get(k, 0)) for k in baseline
              if scores.get(k, 0) < baseline[k]}

    for key, (was, now) in sorted(better.items()):
        print(f"improved  {key}: {was} -> {now}")

    if not new and not worse:
        if better:
            print("\nRun with --update to lock in the improvements above.")
        print(f"\nOK - no function above {MAX_COMPLEXITY} outside the baseline.")
        return 0

    print(f"\nCyclomatic complexity above {MAX_COMPLEXITY}:")
    for key, val in sorted(new.items()):
        print(f"  NEW     {key}: {val}")
    for key, (was, now) in sorted(worse.items()):
        print(f"  WORSE   {key}: {was} -> {now}")
    print(
        "\nSplit the function, or - if the branching is irreducible - run\n"
        "  python3 .github/scripts/check_complexity.py --update\n"
        "and say why in the commit message."
    )
    return 1


if __name__ == "__main__":
    sys.exit(main())
