import json
from pathlib import Path


def normalize_notebook(path: Path) -> bool:
    """Ensure an .ipynb is nbformat 4 compatible.

    Returns True if modified.
    """
    raw = path.read_text(encoding="utf-8")
    nb = json.loads(raw)

    modified = False

    if "cells" not in nb or not isinstance(nb["cells"], list):
        raise ValueError(f"{path}: missing/invalid 'cells'")

    # Add top-level nbformat fields if missing
    if nb.get("nbformat") != 4:
        nb["nbformat"] = 4
        modified = True
    if nb.get("nbformat_minor") is None:
        nb["nbformat_minor"] = 5
        modified = True

    # Add reasonable metadata if missing
    metadata = nb.get("metadata")
    if not isinstance(metadata, dict):
        metadata = {}
        nb["metadata"] = metadata
        modified = True

    if "kernelspec" not in metadata:
        metadata["kernelspec"] = {
            "display_name": "Python 3",
            "language": "python",
            "name": "python3",
        }
        modified = True

    if "language_info" not in metadata:
        metadata["language_info"] = {"name": "python", "version": "3"}
        modified = True

    # Normalize cells
    for cell in nb["cells"]:
        if not isinstance(cell, dict):
            raise ValueError(f"{path}: non-dict cell")

        cell_type = cell.get("cell_type")
        if cell_type not in ("code", "markdown", "raw"):
            # Keep unknown cell types but don't touch.
            continue

        # Normalize source into list-of-lines
        src = cell.get("source", [])
        if isinstance(src, str):
            cell["source"] = src.splitlines()
            modified = True
        elif isinstance(src, list):
            # Ensure all elements are strings
            if any(not isinstance(x, str) for x in src):
                cell["source"] = [str(x) for x in src]
                modified = True
        else:
            cell["source"] = [str(src)]
            modified = True

        if cell_type == "code":
            if "execution_count" not in cell:
                cell["execution_count"] = None
                modified = True
            if "outputs" not in cell or not isinstance(cell["outputs"], list):
                cell["outputs"] = []
                modified = True

            # Make sure metadata exists
            if "metadata" not in cell or not isinstance(cell["metadata"], dict):
                cell["metadata"] = {}
                modified = True

    if modified:
        path.write_text(json.dumps(nb, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    return modified


def main() -> None:
    nb_dir = Path(__file__).resolve().parents[1] / "notebook2"
    if not nb_dir.exists():
        raise SystemExit(f"not found: {nb_dir}")

    changed = 0
    checked = 0
    for nb_path in sorted(nb_dir.glob("*.ipynb")):
        checked += 1
        try:
            if normalize_notebook(nb_path):
                changed += 1
                print(f"normalized: {nb_path.name}")
        except Exception as exc:
            print(f"ERROR: {nb_path.name}: {exc}")

    print(f"Done. checked={checked} changed={changed}")


if __name__ == "__main__":
    main()
