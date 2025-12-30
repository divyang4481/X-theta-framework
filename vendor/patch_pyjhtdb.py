import os, pathlib

ROOT = pathlib.Path(r"vendor\pyJHTDB-20200909.0")
if not ROOT.exists():
    raise SystemExit(f"Missing folder: {ROOT}. Did extraction work?")

repl = {
    "np.int": "int",
    "np.float": "float",
    "np.bool": "bool",
    "np.object": "object",
}

patched_files = 0
patched_hits = 0

for p in ROOT.rglob("*.py"):
    s = p.read_text(encoding="utf-8", errors="ignore")
    s2 = s
    for a,b in repl.items():
        if a in s2:
            s2 = s2.replace(a,b)
    if s2 != s:
        p.write_text(s2, encoding="utf-8")
        patched_files += 1
        patched_hits += 1

print("Patched files:", patched_files)

pyproject = ROOT / "pyproject.toml"
if pyproject.exists():
    pyproject.rename(ROOT / "pyproject.toml.bak")
    print("Renamed pyproject.toml -> pyproject.toml.bak (disable PEP517)")
else:
    print("No pyproject.toml found (ok)")
