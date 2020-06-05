import matplotlib.pyplot as plt
from tqdm import tqdm
from functions import basis_sets, functionals, reactions, molecules, load_rxn
import yaml

rawdata = yaml.load(open("data_sets/raw_data.yaml"), Loader=yaml.SafeLoader)
basis_types = ["mw", "gto"]

# Extract data from rawdata to a more convenient format
nels = set([nel for mol in molecules for rxn in reactions if (nel := rawdata[rxn]["bp86"]["def2qzvpp"][mol]["no_electrons"]) is not None])
data = {basis: {nel: {func: [] for func in functionals} for nel in nels} for basis in basis_types}
for basis in tqdm(basis_types, desc="Restructuring data..."):
    for nel in nels:
        for func in functionals:
            for rxn in reactions:
                for mol in molecules:
                    if rawdata[rxn][func]["def2qzvpp"][mol]["no_electrons"] == nel:
                        walltime = rawdata[rxn][func]["def2qzvpp"][mol][f"{basis}_walltime"]
                        ncores = rawdata[rxn][func]["def2qzvpp"][mol][f"{basis}_no_cores"]
                        nscfs = rawdata[rxn][func]["def2qzvpp"][mol][f"{basis}_no_scf_cycles"]

                        if not any([walltime is None, ncores is None, nscfs is None]):
                            data[basis][nel][func].append(walltime / 3600 * ncores / nscfs)

# Set up figure
FS = 8
fig, ax = plt.subplots(figsize=(6, 3), dpi=300)
ax.set_yscale("log")
ax.set_ylabel("CPU hours per SCF cycle", fontsize=FS)
ax.set_xlabel("Number of electrons", fontsize=FS)
[ax.tick_params(axis, labelsize=FS) for axis in "xy"]

# Plot data
for i, basis in enumerate(tqdm(basis_types, desc="Plotting data...")):
    for j, nel in enumerate(nels):
        for func in functionals:
            ys = list(set(data[basis][nel][func]))  # Use set to ignore duplicates
            xs = [nel for _ in ys]

            if basis == "mw":
                ax.scatter(xs, ys, c="skyblue" if func == "bp86" else "crimson",
                           s=12, edgecolors="black", linewidths=1,
                           label=f"{basis.upper()}6 / {func.upper()}" if j == 0 else None)
            else:
                ax.scatter(xs, ys, c="lightgreen" if func == "bp86" else "goldenrod",
                           s=12, edgecolors="black", linewidths=1,
                           label=f"def2-QZVPP / {func.upper()}" if j == 0 else None)

plt.grid()
plt.legend(fancybox=True, ncol=2, fontsize=FS)
plt.tight_layout()
plt.savefig("figs/scf_timings.png")