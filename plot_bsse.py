import yaml
import numpy as np
from collections import OrderedDict
import matplotlib.pyplot as plt

with open("data_sets/counterpoise.yaml") as f:
    data = yaml.load(f, Loader=yaml.Loader)

basis_sets = OrderedDict({"def2qzvpp": "def2-qzvpp",
                          "def2tzvp": "def2-tzvp",
                          "aug-6311gdp": "6-311+g(d,p)",
                          "6311gdp": "6-311g(d,p)"})

fig, axes = plt.subplots(nrows=2, figsize=(10, 7.5), dpi=100, sharex="all")
no_reactions = 16
width = 0.2
colors = ["crimson", "lightblue", "blue", "lightgreen"]
lw = 1.25
ec = "black"
fs = 16
ind = np.arange(no_reactions)
cntr = 0
for ax, func in zip(axes, data.keys()):
    i = 0
    for baskey, basval in basis_sets.items():
        rxn, bas = zip(*data[func][baskey])
        ax.bar(ind+i*width, bas, width, color=colors[i], edgecolor=ec, linewidth=lw, label=basval)
        i += 1

    ax.set_xticks(ind + 3/2*width)
    ax.set_xticklabels(rxn)

    ax.set_ylabel(r'BSSE $\left[ \frac{kcal}{mol} \right]$', fontsize=fs)
    ax.set_title(f"Functional: {func}", fontsize=fs)
    [ax.tick_params(dim, labelsize=fs-2) for dim in ("x", "y")]
    ax.grid(axis="y")

    if cntr == 0:
        ax.legend(fontsize=fs, ncol=1)
    cntr += 1

plt.tight_layout()
plt.savefig("figs/bsse.png")
# plt.show()
