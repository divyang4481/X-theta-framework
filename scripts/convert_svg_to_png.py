from __future__ import annotations

import sys
from pathlib import Path


def main() -> int:
    try:
        import cairosvg  # type: ignore
    except Exception as e:  # pragma: no cover
        print("ERROR: CairoSVG not installed in this environment:", e, file=sys.stderr)
        return 2

    # Defaults to the project figure if no args provided
    svg = (
        sys.argv[1]
        if len(sys.argv) > 1
        else r"c:\workspace\Physics\X-theta-framework\paper\figs\drone_nav_overview.svg"
    )
    png = (
        sys.argv[2]
        if len(sys.argv) > 2
        else r"c:\workspace\Physics\X-theta-framework\paper\figs\drone_nav_overview.png"
    )
    dpi = int(sys.argv[3]) if len(sys.argv) > 3 else 300

    svg_path = Path(svg)
    png_path = Path(png)

    if not svg_path.exists():
        print(f"ERROR: SVG not found: {svg_path}", file=sys.stderr)
        return 3

    try:
        cairosvg.svg2png(url=str(svg_path), write_to=str(png_path), dpi=dpi)
    except Exception as e:
        print("ERROR: Conversion failed:", e, file=sys.stderr)
        return 4

    print(f"OK: wrote {png_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
"""
Convert the overview SVG to a PNG using CairoSVG.

Requirements:
  pip install --user cairosvg==2.7.1 pillow==10.4.0 tinycss2==1.3.0 cssselect2==0.7.0 cairocffi

Output:
  Writes paper/figs/drone_nav_overview.png (300 DPI)
"""
from pathlib import Path
import sys

SVG = Path(r"c:/workspace/Physics/X-theta-framework/paper/figs/drone_nav_overview.svg")
PNG = Path(r"c:/workspace/Physics/X-theta-framework/paper/figs/drone_nav_overview.png")

def main() -> int:
    if not SVG.exists():
        print(f"ERROR: SVG not found: {SVG}")
        return 2
    try:
        import cairosvg  # type: ignore
    except Exception as e:
        print("CairoSVG is not installed or failed to import.")
        print("Install with: python -m pip install --user cairosvg==2.7.1 pillow==10.4.0 tinycss2==1.3.0 cssselect2==0.7.0 cairocffi")
        print(f"Import error: {e}")
        return 3

    try:
        PNG.parent.mkdir(parents=True, exist_ok=True)
        cairosvg.svg2png(url=str(SVG), write_to=str(PNG), dpi=300)
        print(f"OK: wrote {PNG}")
        return 0
    except Exception as e:
        print("ERROR: conversion failed.")
        print(e)
        return 4

if __name__ == "__main__":
    raise SystemExit(main())
