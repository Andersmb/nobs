import yaml
from glob import glob
from pprint import pprint
import sys
import os
sys.path.append("/Users/abr121/Documents/dev/QueueGui3/output_parsers")
from orca import OrcaOut
from mrchem import MrchemOut

AU2KCAL = 627.509
functionals = ["bp86", "pbe0"]
basis_sets = ["6311gdp", "aug-6311gdp", "def2tzvp", "def2qzvpp", "def2svp"]


def stem(s):
    return s.split(".")[0]


def get_bsse_data(skip=None) -> dict:
    if skip is None:
        skip = []
    root = "/Volumes/external/phd/rsync-project-stallo/nobs/calcs/cp"
    outputs = sorted(glob(os.path.join(root, "*.out")))

    data = {func: {bas: []
                for bas in basis_sets}
            for func in functionals}

    for out in outputs:
        jobname = stem(os.path.basename(out))
        rxn = jobname.split("_")[1]
        func = jobname.split("_")[2]
        basis = jobname.split("_")[3]
        if any([rxn in skip, func in skip, basis in skip]):
            continue

        output = OrcaOut(out)
        if not output.normaltermination():
            print(f"Skipped: {jobname} due to bad termination")
            continue
        elif func not in functionals:
            print(f"Skipped: {jobname} due to incorrect functional")
            continue
        bsse = output.bsse("counterpoise_correction_kcalmol")
        data[func][basis].append((rxn, bsse))

    # Sort the reaction,bsse tuples
    for func in data.keys():
        for basis in data[func].keys():
            data[func][basis] = sorted(data[func][basis], key=lambda x: int(x[0][1:]))

    with open("data_sets/bsse.yaml", "w") as f:
        yaml.dump(data, f, default_flow_style=False, indent=4)
    return data


def get_old_data(d3=True, zpe=True) -> dict:
    """Return the old GTO reaction energies."""
    with open("data_sets/old_data.yaml") as f:
        data = yaml.load(f, Loader=yaml.Loader)

    # Get relevant data in convenient format
    skip_basis= ["def2svp", "6311gdp"]
    d = {func: {basis: []
                for basis in basis_sets if basis not in skip_basis}
         for func in functionals}

    for rxn in data.keys():
        for func in functionals:
            for basis in basis_sets:
                if basis in skip_basis:
                    continue
                delta_e = data[rxn][func][basis]["complex"]["energy"] \
                - data[rxn][func][basis]["fragment1"]["energy"] \
                - data[rxn][func][basis]["fragment2"]["energy"]

                if zpe:
                    delta_zpe = data[rxn][func][basis]["complex"]["zpe"] \
                    - data[rxn][func][basis]["fragment1"]["zpe"] \
                    - data[rxn][func][basis]["fragment2"]["zpe"]
                else:
                    delta_zpe = 0

                if d3:
                    delta_d3 = data[rxn][func][basis]["complex"]["d3"] \
                    - data[rxn][func][basis]["fragment1"]["d3"] \
                    - data[rxn][func][basis]["fragment2"]["d3"]
                else:
                    delta_d3 = 0

                if "mw" not in basis:
                    bsse = data[rxn][func][basis]["bsse"] * AU2KCAL
                    DeltaE = (delta_e + delta_zpe +delta_d3) * AU2KCAL
                    DeltaE_CP = DeltaE + bsse
                    d[func][basis].append((rxn, DeltaE, DeltaE_CP))
    return d


def get_old_mwref(d3=False) -> dict:
    """Return the MW6 reference reaction energies, optionally with D3 correction."""
    with open("data_sets/old_data.yaml") as f:
        data = yaml.load(f, Loader=yaml.Loader)

    ref = {func: [] for func in functionals}
    for rxn in data.keys():
        for func in functionals:
            DeltaE = data[rxn][func]["mw6"]["complex"]["energy"] - data[rxn][func]["mw6"]["fragment1"]["energy"] - data[rxn][func]["mw6"]["fragment2"]["energy"]
            delta_d3 = data[rxn][func]["def2svp"]["complex"]["d3"] - data[rxn][func]["def2svp"]["fragment1"]["d3"] - data[rxn][func]["def2svp"]["fragment2"]["d3"] if d3 else 0

            ref[func].append((rxn, (DeltaE + delta_d3)*AU2KCAL))

    return ref
