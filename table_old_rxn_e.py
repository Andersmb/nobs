from functions import get_old_data, get_old_mwref, basis_sets

GTO = get_old_data(d3=True, zpe=False)
REF = get_old_mwref(d3=True)

skip = ["6311gdp", "def2svp"]

for func in GTO.keys():
    for basis in basis_sets:
        if basis in skip:
            continue
        print(f"-------{func}_{basis}-------")
        for r, e, ecp in GTO[func][basis]:
            print(f"{r}\t&\t{e:.2f} ({ecp:.2f})")

print("MW REF")
for func in REF.keys():
    print(f"-----{func}--------")
    for r, e in REF[func]:
        print(f"{r}\t&\t{e:.2f}")