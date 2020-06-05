import yaml
import numpy as np
from collections import OrderedDict
import matplotlib.pyplot as plt

with open("data_sets/bsse.yaml") as f:
    data = yaml.load(f, Loader=yaml.Loader)

basis_sets = OrderedDict({"def2qzvpp": "def2-qzvpp",
                          "def2tzvp": "def2-tzvp",
                          "def2svp": "def2-svp",
                          "6311gdp": "6-311g(d,p),",
                          "631g": "6-31g"})

fig, axes = plt.subplots(nrows=2, figsize=(10, 5), dpi=100, sharex="all")
n_rxn = len(data["bp86"]["def2svp"])
n_basis = len(data["bp86"])
width = 1/(n_basis+2)
colors = ["crimson", "lightblue", "blue", "lightgreen", "silver", "goldenrod"]
lw = 1.25
ec = "black"
fs = 16
pos = np.arange(n_rxn)
cntr = 0

for ax, func in zip(axes, data.keys()):
    i = 0
    for baskey, basval in basis_sets.items():
        rxn, bas = zip(*data[func][baskey])
        ax.bar(pos+i*width, bas, width, color=colors[i], edgecolor=ec, linewidth=lw, label=basval)
        i += 1

    ax.set_xticks(pos + n_basis/2*width)
    ax.set_xticklabels(rxn)

    ax.set_ylabel(r'BSSE $\left[ \frac{kcal}{mol} \right]$', fontsize=fs)
    ax.set_title(f"Functional: {func}", fontsize=fs)
    [ax.tick_params(dim, labelsize=fs-3) for dim in ("x", "y")]
    ax.grid(axis="y")

    if cntr == 0:
        ax.legend(fontsize=fs-3, ncol=3)
    cntr += 1

plt.tight_layout()
plt.savefig("figs/bsse.png")
# plt.show()
