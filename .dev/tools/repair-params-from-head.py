#!/usr/bin/env python3
"""Restore the @param entries of roxygen blocks to their committed (HEAD) form.

Why this exists: the first version of .dev/tools/fix-roxygen-drift.R read the
committed Rd's \\arguments with `as.character()` on the parsed Rd fragment. On
an Rd fragment that deparses -- `\\code{colnames(expr)}` comes back as the
string `list("colnames(expr)")` -- so the repaired roxygen blocks were filled
with deparsed junk, and the fixer's name matching did not understand roxygen's
comma-separated `@param a,b` tags, so it ADDED parameters that were already
documented. R CMD check then reported "Duplicated \\argument entries" and
"Undocumented arguments".

The fix here is deliberately blunt and does not try to be clever: for every
roxygen block, the set of @param entries is taken from HEAD, and any @param
entry that HEAD does not have is dropped. Everything else in the block --
notably the @examples the agents wrote -- is left alone.
"""
import re
import subprocess
import sys
from collections import OrderedDict

TAG = re.compile(r"^(\s*)#'\s*@(\w+)\s*(.*)$")
PARAM = re.compile(r"^(\s*)#'\s*@param\s+([^ ]+)(.*)$")


def blocks(lines):
    """Yield (start, end_exclusive, owner, entries) for each roxygen block.

    `entries` is an OrderedDict of tag -> list of raw lines, so a multi-line
    tag keeps its continuation lines.
    """
    i = 0
    n = len(lines)
    while i < n:
        if not re.match(r"^\s*#'", lines[i]):
            i += 1
            continue
        start = i
        while i < n and re.match(r"^\s*#", lines[i]):
            i += 1
        block = lines[start:i]
        # the owner: the next non-comment, non-blank line
        j = i
        while j < n and (not lines[j].strip() or re.match(r"^\s*#", lines[j])):
            j += 1
        owner = None
        if j < n:
            m = re.match(r"^\s*([A-Za-z_.][A-Za-z0-9._]*)\s*<-", lines[j])
            if m:
                owner = m.group(1)
        entries = OrderedDict()
        k = 0
        while k < len(block):
            m = TAG.match(block[k])
            if m:
                tag = m.group(2)
                group = [block[k]]
                k2 = k + 1
                while k2 < len(block) and re.match(r"^\s*#'", block[k2]) \
                        and not TAG.match(block[k2]):
                    group.append(block[k2])
                    k2 += 1
                entries.setdefault(tag, []).extend(group)
                k = k2
            else:
                k += 1
        yield start, i, owner, entries


def param_names(lines):
    """The set of parameter names a @param entry documents."""
    out = set()
    for l in lines:
        m = PARAM.match(l)
        if m:
            for nm in m.group(2).split(","):
                nm = nm.strip()
                if nm:
                    out.add(nm)
    return out


def main(paths):
    for path in paths:
        head = subprocess.run(
            ["git", "show", "HEAD:" + path],
            capture_output=True, text=True, check=True).stdout.splitlines()
        cur = open(path, encoding="utf-8").read().splitlines()

        head_blocks = list(blocks(head))
        head_by_owner = OrderedDict()
        for s, e, owner, ent in head_blocks:
            if owner:
                head_by_owner[owner] = ent

        # rebuild the current file's blocks with HEAD's @param entries
        out = []
        i = 0
        replaced = dropped = 0
        for s, e, owner, ent in blocks(cur):
            out.extend(cur[i:s])
            new_block = list(cur[s:e])
            if owner and owner in head_by_owner:
                head_params = head_by_owner[owner].get("param", [])
                head_names = param_names(head_params)
                # find the @param group(s) in the current block
                idx = [k for k, l in enumerate(new_block) if PARAM.match(l)]
                if idx:
                    first, last = idx[0], idx[-1]
                    cur_names = set()
                    for k in idx:
                        cur_names |= param_names([new_block[k]])
                    extra = cur_names - head_names
                    if extra or any("list(\"" in l for l in new_block[first:last + 1]):
                        # replace the whole @param group with HEAD's, and drop
                        # any parameter HEAD does not document
                        new_block = (new_block[:first] + head_params +
                                     new_block[last + 1:])
                        replaced += 1
                        dropped += len(extra)
            out.extend(new_block)
            i = e
        out.extend(cur[i:])
        new_text = "\n".join(out) + "\n"
        if new_text != open(path, encoding="utf-8").read():
            open(path, "w", encoding="utf-8").write(new_text)
        print(f"{path}: blocks repaired={replaced} extra params dropped={dropped}")


if __name__ == "__main__":
    main(sys.argv[1:])
