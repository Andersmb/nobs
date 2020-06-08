import yaml
from glob import glob
from tqdm import tqdm, trange
from collections import OrderedDict
from pprint import pprint
import sys
import os
sys.path.append("/Users/abr121/Documents/dev/QueueGui3/output_parsers")
from orca import OrcaOut
from mrchem import MrchemOut
from collections import OrderedDict

AU2KCAL = 627.509
NUM_RXN = 28
functionals = ["bp86", "pbe0"]
basis_sets = ["6311gdp", "def2tzvp", "def2qzvpp", "def2svp", "631g"]  # We skip 6311+G(d,p) for now. Need to converge more calculations.
reactions = [f"r{i}" for i in range(NUM_RXN)]
molecules = ["complex", "frag1", "frag2"]

basis_sets_ordered = OrderedDict({"def2qzvpp": "def2-qzvpp",
                                  "def2tzvp": "def2-tzvp",
                                  "def2svp": "def2-svp",
                                  "6311gdp": "6-311g(d,p)",
                                  "631g": "6-31g"})

files = {"raw": "data_sets/raw_data.yaml",
         "bsse": "data_sets/bsse.yaml",
         "rxn": "/Volumes/external/phd/rsync-project-stallo/nobs/__reactions__.yaml",
         "old": "data_sets/old_data.yaml"}


def load_data(d) -> dict:
    assert d in files.keys(), f"File not recognized! Must be any of {', '.join(files.keys())}"
    with open(files[d]) as f:
        return yaml.load(f, Loader=yaml.Loader)


def stem(s) -> str:
    return s.split(".")[0]


def get_bsse_data(skip=None) -> dict:
    if skip is None:
        skip = []
    root = "/Volumes/external/phd/rsync-project-stallo/nobs/calcs/cp"
    outputs = sorted(glob(os.path.join(root, "*.out")))

    data = {func: {bas: []
                for bas in basis_sets}
            for func in functionals}

    for out in tqdm(outputs, desc="Collecting BSSE data"):
        jobname = stem(os.path.basename(out))
        rxn = jobname.split("_")[0]
        func = jobname.split("_")[1]
        basis = jobname.split("_")[2]

        if any([rxn in skip, func in skip, basis in skip]):
            continue

        output = OrcaOut(out)
        if not output.normaltermination():
            tqdm.write(f"Skipped: {jobname} due to bad termination")
            continue
        elif func not in functionals:
            tqdm.write(f"Skipped: {jobname} due to incorrect functional")
            continue
        bsse = output.cmp_var("bsse_kcalmol")
        data[func][basis].append((rxn, bsse))

    # Sort the reaction,bsse tuples
    for func in data.keys():
        for basis in data[func].keys():
            data[func][basis] = sorted(data[func][basis], key=lambda x: int(x[0][1:]))

    with open(files["bsse"], "w") as f:
        yaml.dump(data, f, default_flow_style=False, indent=4)
    return data


def get_old_data(d3=False, zpe=False) -> dict:
    """Return the old GTO reaction energies."""
    data = load_data("old")

    # Get relevant data in convenient format
    skip_basis= ["def2svp", "6311gdp", "631g"]
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

                delta_zpe = data[rxn][func]["def2svp"]["complex"]["zpe"] \
                - data[rxn][func]["def2svp"]["fragment1"]["zpe"] \
                - data[rxn][func]["def2svp"]["fragment2"]["zpe"] if zpe else 0.0

                delta_d3 = data[rxn][func]["def2svp"]["complex"]["d3"] \
                - data[rxn][func]["def2svp"]["fragment1"]["d3"] \
                - data[rxn][func]["def2svp"]["fragment2"]["d3"] if d3 else 0.0

                if "mw" not in basis:
                    bsse = data[rxn][func][basis]["bsse"] * AU2KCAL
                    DeltaE = (delta_e + delta_zpe + delta_d3) * AU2KCAL
                    DeltaE_CP = DeltaE + bsse
                    d[func][basis].append((rxn, DeltaE, DeltaE_CP))
    return d


def get_old_mwref(d3=False) -> dict:
    """Return the MW6 reference reaction energies, optionally with D3 correction."""
    data = load_data("old")

    ref = {func: [] for func in functionals}
    for rxn in data.keys():
        for func in functionals:
            DeltaE = data[rxn][func]["mw6"]["complex"]["energy"] - data[rxn][func]["mw6"]["fragment1"]["energy"] - data[rxn][func]["mw6"]["fragment2"]["energy"]
            delta_d3 = data[rxn][func]["def2svp"]["complex"]["d3"] - data[rxn][func]["def2svp"]["fragment1"]["d3"] - data[rxn][func]["def2svp"]["fragment2"]["d3"] if d3 else 0

            ref[func].append((rxn, (DeltaE + delta_d3)*AU2KCAL))

    return ref


def get_raw_data() -> dict:
    """"""
    root_mw = "/Volumes/external/phd/rsync-project-saga/nobs/calcs/mw"
    root_go = "/Volumes/external/phd/rsync-project-stallo/nobs/calcs/go"
    root_sp = "/Volumes/external/phd/rsync-project-stallo/nobs/calcs/sp"
    root_cp = "/Volumes/external/phd/rsync-project-stallo/nobs/calcs/cp"

    # Load __reactions file
    __reactions = load_data("rxn")

    # Initialize data structure
    data = {rxn: {func: {bas: {mol: {
        "gto_energy": None,
        "gto_d3": None,
        "gto_zpe": None,
        "gto_no_cores": None,
        "gto_no_scf_cycles": None,
        "gto_walltime": None,
        "gto_scf_timings": None,

        "mw_energy": None,
        "mw_no_cores": None,
        "mw_no_scf_cycles": None,
        "mw_walltime": None,

        "no_atoms": None,
        "no_electrons": None
    } for mol in molecules} for bas in basis_sets} for func in functionals} for rxn in reactions}

    # Add field for reference energy, if any
    for rxn in reactions:
        data[rxn]["ref"] = None

    # Add field for BSSE
    for rxn in reactions:
        for func in functionals:
            for bas in basis_sets:
                data[rxn][func][bas]["bsse"] = {"au": None, "kcalmol": None}

    # Fill up with actual data
    for rxn in tqdm(reactions, desc="Collecting data"):
        if not __reactions[rxn]["COMPLETE"]:
            continue  # We skip all reactions where all data is not yet gathered

        # Determine from which reaction to grab the fragments
        base_complex = __reactions[rxn]["base_complex"]  # Only relevant for reaction 16/17
        base_frag1 = __reactions[rxn]["base_frag1"]
        base_frag2 = __reactions[rxn]["base_frag2"]
        bases = [base_complex, base_frag1, base_frag2]

        for func in functionals:
            for bas in basis_sets:
                for base, mol in zip(bases, molecules):
                    # We will use the 'base_mol' syntax to obtain properties from the correct files

                    # Define files
                    f_e = os.path.join(root_sp, base, func, bas, mol+".out")
                    f_zpe = os.path.join(root_go, base, func, mol+".out")
                    f_cp = os.path.join(root_cp, f"{rxn}_{func}_{bas}.out")
                    f_mw = os.path.join(root_mw, base, func, mol+".out")

                    # Generate output file instances
                    try:
                        output_e = OrcaOut(f_e)
                        output_zpe = OrcaOut(f_zpe)
                        output_cp = OrcaOut(f_cp)
                        output_mw = MrchemOut(f_mw)
                    except FileNotFoundError:
                        tqdm.write("File not found...")

                    # Get data
                    data[rxn][func][bas][mol]["gto_energy"] = output_e.final_total_energy()
                    data[rxn][func][bas][mol]["gto_d3"] = output_e.dispersion_correction()
                    data[rxn][func][bas][mol]["gto_zpe"] = output_zpe.zero_point_energy_correction()
                    data[rxn][func][bas][mol]["gto_walltime"] = output_e.walltime()
                    data[rxn][func][bas][mol]["gto_scf_timings"] = output_e.scf_timings()[0]
                    data[rxn][func][bas][mol]["gto_no_cores"] = output_e.no_cores()
                    data[rxn][func][bas][mol]["gto_no_scf_cycles"] = output_e.no_scfcycles()[0]

                    data[rxn][func][bas][mol]["mw_energy"] = output_mw.final_energy_pot()
                    data[rxn][func][bas][mol]["mw_walltime"] = output_mw.walltime()
                    data[rxn][func][bas][mol]["mw_no_cores"] = output_mw.no_cores()
                    data[rxn][func][bas][mol]["mw_no_scf_cycles"] = output_mw.no_scfcycles()

                    data[rxn][func][bas]["bsse"]["au"] = output_cp.cmp_var("bsse_au")
                    data[rxn][func][bas]["bsse"]["kcalmol"] = output_cp.cmp_var("bsse_kcalmol")

                    data[rxn][func][bas][mol]["no_atoms"] = __reactions[rxn][mol]["natoms"]
                    data[rxn][func][bas][mol]["no_electrons"] = __reactions[rxn][mol]["nel"]

    with open(files["raw"], "w") as f:
        yaml.dump(dict(data), f)

    return data


def get_reaction_energies(zpe=False, d3=False) -> tuple:
    """Return tuple of dicts with reaction energies for completed reactions from raw data."""
    rawdata = load_data("raw")
    __reactions = load_data("rxn")
    gto = {rxn: {func: {bas: {cp: None for cp in ["cp", "noncp"]} for bas in basis_sets} for func in functionals} for rxn in reactions if __reactions[rxn]["COMPLETE"]}
    mw = {rxn: {func: None for func in functionals} for rxn in reactions if __reactions[rxn]["COMPLETE"]}

    for rxn in reactions:
        if not __reactions[rxn]["COMPLETE"]:
            continue
        for func in functionals:
            e_c = rawdata[rxn][func]["def2qzvpp"]["complex"]["mw_energy"]
            e_f1 = rawdata[rxn][func]["def2qzvpp"]["frag1"]["mw_energy"]
            e_f2 = rawdata[rxn][func]["def2qzvpp"]["frag2"]["mw_energy"]

            d3_c = rawdata[rxn][func]["def2qzvpp"]["complex"]["gto_d3"]
            d3_f1 = rawdata[rxn][func]["def2qzvpp"]["frag1"]["gto_d3"]
            d3_f2 = rawdata[rxn][func]["def2qzvpp"]["frag2"]["gto_d3"]

            zpe_c = rawdata[rxn][func]["def2qzvpp"]["complex"]["gto_zpe"]
            zpe_f1 = rawdata[rxn][func]["def2qzvpp"]["frag1"]["gto_zpe"]
            zpe_f2 = rawdata[rxn][func]["def2qzvpp"]["frag2"]["gto_zpe"]

            delta_e = e_c - e_f1 - e_f2
            delta_d3 = d3_c - d3_f1 - d3_f2 if d3 else 0
            delta_zpe = zpe_c - zpe_f1 - zpe_f2 if zpe else 0

            mw[rxn][func] = (delta_e + delta_d3 + delta_zpe) * AU2KCAL

            for bas in basis_sets:
                e_c = rawdata[rxn][func][bas]["complex"]["gto_energy"]
                e_f1 = rawdata[rxn][func][bas]["frag1"]["gto_energy"]
                e_f2 = rawdata[rxn][func][bas]["frag2"]["gto_energy"]

                d3_c = rawdata[rxn][func][bas]["complex"]["gto_d3"]
                d3_f1 = rawdata[rxn][func][bas]["frag1"]["gto_d3"]
                d3_f2 = rawdata[rxn][func][bas]["frag2"]["gto_d3"]

                zpe_c = rawdata[rxn][func][bas]["complex"]["gto_zpe"]
                zpe_f1 = rawdata[rxn][func][bas]["frag1"]["gto_zpe"]
                zpe_f2 = rawdata[rxn][func][bas]["frag2"]["gto_zpe"]

                delta_e = e_c - e_f1 - e_f2
                delta_d3 = d3_c - d3_f1 - d3_f2 if d3 else 0
                delta_zpe = zpe_c - zpe_f1 - zpe_f2 if zpe else 0

                cp = rawdata[rxn][func][bas]["bsse"]["au"]

                gto[rxn][func][bas]["cp"] = (delta_e + delta_d3 + delta_zpe + cp) * AU2KCAL
                gto[rxn][func][bas]["noncp"] = (delta_e + delta_d3 + delta_zpe) * AU2KCAL

    return gto, mw





if __name__ == "__main__":
    # get_raw_data()
    # get_bsse_data(skip=["aug-6311gdp", "r27"])
    pass