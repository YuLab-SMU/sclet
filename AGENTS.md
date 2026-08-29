# AGENTS.md — sclet workspace layout

This is the single entry point for the whole `sclet` workspace. The workspace
contains **separate, nested `git` checkouts** of *different branches* of the
`YuLab-SMU/sclet` repository. The rule is simple:

> **Each directory is its own checkout of its own branch.**
> Entering a directory is equivalent to being on that branch. Do not assume a
> file you see in one directory is the same version as the one in another.

## Directory ↔ branch mapping

| Directory (`<workspace>/`)  | Git branch | What it contains                                                    | Remote |
|-----------------------------|------------|---------------------------------------------------------------------|--------|
| `.` (this top directory)    | `devel`    | **`sclet` R package** source (`R/`, `man/`, `NAMESPACE`, `DESCRIPTION`, `tests/`) | `git@github.com:YuLab-SMU/sclet.git` |
| `bookdown/`                 | `docs`     | Bookdown **source** (`*.Rmd`, `_bookdown.yml`, `index.Rmd`, `references.bib`, ...) | `git@github.com:YuLab-SMU/sclet.git` |
| `bookdown/docs/`            | `gh-pages` | **Published site** (built `.html`, `search_index.json`, `libs/`, `figures/`, `.nojekyll`) | `git@github.com:YuLab-SMU/sclet.git` |

`bookdown` is already git-ignored at this top level, so the nested checkouts do
not pollute the R package repository's index.

## What goes where

- **R package work** (`sclet` the package): this top-level checkout, branch `devel`.
- **Edits to the book**: change files under `bookdown/` (the `docs` branch).
  Add/edit `*.Rmd`, `_bookdown.yml`, `references.bib`, figures, etc. here.
- **Rebuilding / deploying the site**: rebuild the book, and the generated
  `docs/` output is the `bookdown/docs/` checkout (the `gh-pages` branch).
  Publish `gh-pages` to update the live book site.

## Pitfalls to avoid

- Do **not** mix `.Rmd` source and built `.html` output. They live in different
  branches and different directories for a reason.
- Each directory is an **independent repository**. Moving/renaming files across
  them does not move them across branches — each has its own `.git`.
- A change committed in `bookdown/` is **not** visible in `bookdown/docs/` (and
  vice-versa) until the book is actually rebuilt and the `gh-pages` branch is
  re-pushed.
- Never run a command in `bookdown/docs/` expecting it to affect the `docs`
  source branch, and never run bookdown builds from the wrong directory.
- When you `cd` into a directory, keep in mind it is now a **different branch** —
  this top directory is `devel`, not `docs`.

## Quick reference

```bash
# Confirm which branch you are on before editing
git -C .                 branch --show-current   # → devel
git -C bookdown          branch --show-current   # → docs
git -C bookdown/docs     branch --show-current   # → gh-pages

# All checkouts point at the same repository remote
git -C .                 remote get-url origin
git -C bookdown          remote get-url origin
git -C bookdown/docs     remote get-url origin
```
