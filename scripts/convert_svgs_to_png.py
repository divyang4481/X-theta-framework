#!/usr/bin/env python3
"""
Regenerate all .png assets in paper/figs from their corresponding .svg files
using CairoSVG. Overwrites existing PNGs (useful if some are invalid/corrupted).

Usage:
  python scripts/convert_svgs_to_png.py
"""
from __future__ import annotations

import sys
import os
from pathlib import Path

PNG_SIGNATURE = b"\x89PNG\r\n\x1a\n"

def is_valid_png(path: Path) -> bool:
    try:
        with path.open('rb') as f:
            sig = f.read(8)
        return sig == PNG_SIGNATURE
    except Exception:
        return False

def main() -> int:
    repo_root = Path(__file__).resolve().parents[1]
    figs_dir = repo_root / 'paper' / 'figs'
    if not figs_dir.exists():
        print(f"ERROR: figs directory not found: {figs_dir}")
        return 1

    try:
        from cairosvg import svg2png
    except Exception as e:
        print("ERROR: CairoSVG not installed or failed to import:", e)
        print("Install with: pip install cairosvg pillow tinycss2 cssselect2")
        return 1

    svg_files = sorted(figs_dir.glob('*.svg'))
    if not svg_files:
        print(f"No SVG files found in {figs_dir}")
        return 0

    ok = 0
    failed = 0
    for svg in svg_files:
        png = svg.with_suffix('.png')
        try:
            svg2png(url=str(svg), write_to=str(png))
            if is_valid_png(png):
                print(f"Converted: {svg.name} -> {png.name}")
                ok += 1
            else:
                print(f"WARN: Output is not a valid PNG signature: {png}")
                failed += 1
        except Exception as e:
            print(f"FAILED: {svg.name} -> {png.name}: {e}")
            failed += 1

    print(f"Done. Success: {ok}, Failed: {failed}")
    return 0 if failed == 0 else 2

if __name__ == '__main__':
    sys.exit(main())
