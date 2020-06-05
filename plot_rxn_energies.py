from functions import reactions, basis_sets_ordered, functionals, load_data, get_reaction_energies
import matplotlib.pyplot as plt
import numpy as np

__reactions = load_data("rxn")
gto, mw = get_reaction_energies()

# Convert data to format easier to plot
data_cp = {func: {bas: [] for bas in basis_sets_ordered.keys()} for func in functionals}
data_noncp = {func: {bas: [] for bas in basis_sets_ordered.keys()} for func in functionals}

data_mw = {func: [] for func in functionals}

for rxn in gto.keys():
    for bas in basis_sets_ordered:
        for func in functionals:
            data_cp[func][bas].append((rxn, gto[rxn][func][bas]["cp"]))
            data_noncp[func][bas].append((rxn, gto[rxn][func][bas]["noncp"]))

for func in functionals:
    for rxn in mw.keys():
        data_mw[func].append((rxn, mw[rxn][func]))

fig, axes = plt.subplots(nrows=2, ncols=2, dpi=300, figsize=(6, 3), sharex=True, sharey=True)
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
    i = 1
    _, e = zip(*data_mw[func])
    col[0].bar(pos, e, width, color=colors[0], edgecolor=ec, linewidth=lw, label="MW6" if cntr == 0 else None)
    col[1].bar(pos, e, width, color=colors[0], edgecolor=ec, linewidth=lw, label="MW6" if cntr == 0 else None)

    for baskey, basval in basis_sets_ordered.items():

        rxn, e = zip(*data_cp[func][baskey])
        col[0].bar(pos + i * width, e, width, color=colors[i], edgecolor=ec, linewidth=lw, label=basval if cntr == 0 else None)

        _, e = zip(*data_noncp[func][baskey])
        col[1].bar(pos + i * width, e, width, color=colors[i], edgecolor=ec, linewidth=lw)

        i += 1

    for ax in col:
        ax.set_xticks(pos + (n_basis - 1) / 2 * width)
        ax.set_xticklabels(rxn)

    col[0].legend(fancybox=True, fontsize=fs-3, ncol=2)
    cntr += 1

axes[0,0].set_title("CP")
fig.text(0.04, 0.5, "Reaction energy [kcal/mol]", va="center", rotation="vertical", fontsize=fs)
axes[0,1].set_title("Non-CP")
fig.text(0.95, 0.3, "PBE0", ha="center", va="center")
fig.text(0.95, 0.7, "BP86", ha="center", va="center")
[ax.tick_params(axis="both", labelsize=fs) for col in axes for ax in col]

plt.subplots_adjust(hspace=0, wspace=0)
plt.savefig("figs/reaction_energies.png")