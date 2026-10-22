#!/usr/bin/env python3
"""Strip (or restore) the @examples blocks from the roxygen blocks in R/.

Used to split the work into two commits that are each coherent on their own:

  1. the behaviour changes (seed discipline, the presto network install, the
     duplicate sclet_copilot, the dependency declarations, the CI matrix)
     together with a roxygenise that does NOT yet carry the examples;
  2. the examples, together with the roxygenise that publishes them.

Without the split, the first commit would have roxygen blocks that document a
`seed` parameter the Rd does not mention, and `R CMD check`'s codoc check
would fail on it.

    python3 .dev/tools/strip-examples.py strip <dir>
    python3 .dev/tools/strip-examples.py restore <dir> <backup-dir>
"""
import os
import re
import shutil
import sys

EXAMPLES_TAG = re.compile(r"^\s*#'\s*@examples\b")
ANY_TAG = re.compile(r"^\s*#'\s*@")


def strip_file(path):
    lines = open(path, encoding="utf-8").read().split("\n")
    out = []
    i = 0
    n = len(lines)
    while i < n:
        if EXAMPLES_TAG.match(lines[i]):
            i += 1
            # drop the continuation lines of the examples block
            while i < n and re.match(r"^\s*#'", lines[i]) and not ANY_TAG.match(lines[i]):
                i += 1
            continue
        out.append(lines[i])
        i += 1
    return "\n".join(out)


def main():
    mode = sys.argv[1]
    if mode == "strip":
        root = sys.argv[2]
        changed = 0
        for f in sorted(os.listdir(os.path.join(root, "R"))):
            if not f.endswith(".R"):
                continue
            p = os.path.join(root, "R", f)
            before = open(p, encoding="utf-8").read()
            after = strip_file(p)
            if before != after:
                open(p, "w", encoding="utf-8").write(after)
                changed += 1
        print(f"stripped @examples from {changed} file(s)")
    elif mode == "restore":
        root, backup = sys.argv[2], sys.argv[3]
        for f in sorted(os.listdir(os.path.join(backup, "R"))):
            if f.endswith(".R"):
                shutil.copy2(os.path.join(backup, "R", f), os.path.join(root, "R", f))
        print(f"restored R/ from {backup}")
    else:
        print(__doc__)
        sys.exit(2)


if __name__ == "__main__":
    main()
