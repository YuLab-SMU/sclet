#!/usr/bin/env python
"""Download and cache a small public velocity h5ad for the RegVelo smoke test."""

from __future__ import annotations

import argparse
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--out",
        default="regvelo_smoke_source/pancreas.h5ad",
        help="Output h5ad path.",
    )
    args = parser.parse_args()

    out = Path(args.out).expanduser().resolve()
    out.parent.mkdir(parents=True, exist_ok=True)

    if out.exists():
        print(f"Existing h5ad found: {out}")
        return

    import scvelo as scv

    print("Downloading scvelo.datasets.pancreas() through Python scvelo...")
    adata = scv.datasets.pancreas()

    if "spliced" not in adata.layers or "unspliced" not in adata.layers:
        raise RuntimeError("Downloaded AnnData lacks required spliced/unspliced layers.")

    adata.write_h5ad(out)
    print(f"Wrote h5ad: {out}")
    print(f"Shape: {adata.n_obs} cells x {adata.n_vars} genes")


if __name__ == "__main__":
    main()
