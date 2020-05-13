import yaml
import numpy as np
from collections import OrderedDict
from pprint import pprint
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from functions import get_old_data, get_old_mwref


def plot_re_rxn_energy(cp=True):
    GTO = get_old_data(d3=False, zpe=False)
    REF = get_old_mwref(d3=False)
    basis_sets = OrderedDict({"def2qzvpp": "def2-qzvpp",
                              "def2tzvp": "def2-tzvp",
                              "aug-6311gdp": "6-311+g(d,p)"})

    fig, axes = plt.subplots(nrows=2, figsize=(10, 7.5), dpi=100, sharex="all")
    no_reactions = 7
    no_basis = 3
    width = 1/(no_basis+1)
    colors = ["crimson", "lightblue", "blue"]
    lw = 1.25
    ec = "black"
    fs = 16
    ind = np.arange(no_reactions)
    cntr = 0

    for ax, func in zip(axes, GTO.keys()):
        i = 0
        ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
        for baskey, basval in basis_sets.items():
            rxn, gto, gto_cp = zip(*GTO[func][baskey])
            rxn, mw = zip(*REF[func])

            re = (np.asarray(gto) - np.asarray(mw)) / np.asarray(mw) * 100
            re_cp = (np.asarray(gto_cp) - np.asarray(mw)) / np.asarray(mw) * 100

            if cp:
                ax.bar(ind + i * width, re_cp, width, color=colors[i], edgecolor=ec, linewidth=lw, label=basval)
                filename = "figs/old_gto_vs_mw_cp.png"
            else:
                ax.bar(ind + i * width, re, width, color=colors[i], edgecolor=ec, linewidth=lw, label=basval)
                filename = "figs/old_gto_vs_mw_noncp.png"
            i += 1

        ax.set_xticks(ind + (no_basis-1) / 2 * width)
        ax.set_xticklabels(rxn)

        ax.set_ylabel('Relative Error [%]', fontsize=fs)
        ax.set_title(f"Functional: {func}", fontsize=fs)
        [ax.tick_params(dim, labelsize=fs - 2) for dim in ("x", "y")]
        ax.grid(axis="y")

        if cntr == 0:
            ax.legend(fontsize=fs, ncol=1)
        cntr += 1

    plt.tight_layout()
    plt.savefig(filename)
    #plt.show()


plot_re_rxn_energy(cp=True)
plot_re_rxn_energy(cp=False)