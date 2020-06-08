from functions import reactions, basis_sets_ordered, functionals, load_data, get_reaction_energies
import matplotlib.pyplot as plt
import numpy as np

__reactions = load_data("rxn")
gto, mw = get_reaction_energies()

# Compupte relative errors
data = {cp: {func: {bas: [] for bas in basis_sets_ordered} for func in functionals} for cp in ["cp", "noncp"]}

for cp in ["cp", "noncp"]:
    for rxn in gto.keys():
        for func in functionals:
            e_mw  = mw[rxn][func]
            for bas in basis_sets_ordered:
                e_gto = gto[rxn][func][bas][cp]
                re = (e_gto - e_mw) / e_mw * 100

                data[cp][func][bas].append((rxn, re))

fig, axes = plt.subplots(ncols=2, nrows=2, figsize=(6,3), dpi=300, sharex=True, sharey=True)
n_rxn = len(mw.keys())
n_basis = len(basis_sets_ordered) + 1
width = 1/(n_basis+2)
colors = ["crimson", "lightblue", "blue", "lightgreen", "silver", "goldenrod"]
lw = 1.0
ec = "black"
fs = 9
pos = np.arange(n_rxn)

cntr = 0
for col, func in zip(axes, functionals):
    i = 0
    for baskey, basval in basis_sets_ordered.items():
        rxn, re = zip(*data["cp"][func][baskey])
        col[0].bar(pos + i*width, re, width, color=colors[i], edgecolor=ec, linewidth=lw, label=basval)

        _, re = zip(*data["noncp"][func][baskey])
        col[1].bar(pos + i * width, re, width, color=colors[i], edgecolor=ec, linewidth=lw)

        i += 1

    for ax in col:
        ax.set_xticks(pos + (n_basis - 1) / 2 * width)
        ax.set_xticklabels(rxn)

    if cntr == 0:
        col[0].legend(fancybox=True, fontsize=fs - 3, ncol=2)
    cntr += 1

axes[0,0].set_title("CP")
fig.text(0.04, 0.5, "Relative Error [%]", va="center", rotation="vertical", fontsize=fs)
axes[0,1].set_title("Non-CP")
fig.text(0.95, 0.3, "PBE0", ha="center", va="center")
fig.text(0.95, 0.7, "BP86", ha="center", va="center")
[ax.tick_params(axis="both", labelsize=fs) for col in axes for ax in col]
[ax.grid(axis="y", lw=0.5) for col in axes for ax in col]

plt.subplots_adjust(hspace=0, wspace=0)
plt.savefig("figs/relative_errors.png")