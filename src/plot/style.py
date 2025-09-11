import os
import matplotlib.pyplot as plt

def savefig(path, dpi=300, bbox_inches="tight"):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    plt.savefig(path, dpi=dpi, bbox_inches=bbox_inches)
    plt.close()
