import yaml
from glob import glob
from pprint import pprint
import sys
import os
sys.path.append("/Users/abr121/Documents/dev/QueueGui3/output_parsers")
from orca import OrcaOut
from mrchem import MrchemOut


def stem(s):
    return s.split(".")[0]


def get_bsse_data(skip=None) -> dict:
    if skip is None:
        skip = []
    root = "/Volumes/external/phd/rsync-project-stallo/nobs/calcs/cp"
    outputs = sorted(glob(os.path.join(root, "*.out")))

    data = {func: {bas: []
                for bas in ["6311gdp", "aug-6311gdp", "def2tzvp", "def2qzvpp"]}
            for func in ["bp86", "pbe0"]}

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
        elif func not in ["pbe0", "bp86"]:
            print(f"Skipped: {jobname} due to incorrect functional")
            continue
        bsse = output.bsse("counterpoise_correction_kcalmol")

        data[func][basis].append((rxn, bsse))

    # Sort the reaction,bsse tuples
    for func in data.keys():
        for basis in data[func].keys():
            data[func][basis] = sorted(data[func][basis], key=lambda x: int(x[0][1:]))

    with open("data_sets/counterpoise.yaml", "w") as f:
        yaml.dump(data, f, default_flow_style=False, indent=4)
    return data

get_bsse_data(skip=[f"r{i}" for i in range(16, 20)])