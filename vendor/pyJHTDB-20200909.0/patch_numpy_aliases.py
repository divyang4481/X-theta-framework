import re
from pathlib import Path

ROOT = Path("pyJHTDB")

# Fix common NumPy alias removals (NumPy 2.x)
REPL = [
    (r"\bnp\.int\b", "int"),
    (r"\bnumpy\.int\b", "int"),
    (r"\bnp\.float\b", "float"),
    (r"\bnumpy\.float\b", "float"),
    (r"\bnp\.bool\b", "bool"),
    (r"\bnumpy\.bool\b", "bool"),
    (r"\bnp\.object\b", "object"),
    (r"\bnumpy\.object\b", "object"),
    # Clean accidental duplication from earlier brute patches
    (r"\bnp\.np\.", "np."),
]

patched_files = 0
for p in ROOT.rglob("*.py"):
    s = p.read_text(encoding="utf-8", errors="ignore")
    s2 = s
    for pat, rep in REPL:
        s2 = re.sub(pat, rep, s2)

    # Keep your earlier fixes safe too
    s2 = s2.replace("astype(float32)", "astype(np.float32)")
    s2 = s2.replace("dtype = float32", "dtype = np.float32").replace("dtype=float32", "dtype=np.float32")

    if s2 != s:
        p.write_text(s2, encoding="utf-8")
        patched_files += 1
        print("patched:", p)

print("\nTotal patched files:", patched_files)

# Make tests optional: don't import test.py during package import
initp = ROOT / "__init__.py"
s = initp.read_text(encoding="utf-8", errors="ignore")
if "from .test import test_plain" in s:
    s = s.replace(
        "from .test import test_plain",
        "try:\n    from .test import test_plain\nexcept Exception:\n    test_plain = None"
    )
    initp.write_text(s, encoding="utf-8")
    print("Updated __init__.py to make test_plain optional.")
else:
    print("__init__.py already patched (or different).")
