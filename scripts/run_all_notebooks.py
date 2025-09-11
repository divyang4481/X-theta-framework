import sys, subprocess, pathlib

NB_DIR = pathlib.Path("notebooks")
for nb in sorted(NB_DIR.glob("*.ipynb")):
    cmd = [
        sys.executable, "-m", "papermill",
        str(nb), str(nb)  # in-place execute
    ]
    print(">>", " ".join(cmd))
    subprocess.check_call(cmd)
print("All notebooks executed.")
