"""
Script to iterate through the emolcules database looking for molecules that match the
"[#6X3H1:1]~[#6X3:2](~[#6X3H1])-[#6X3:3](~[#6X3H1])~[#6X3H1:4]" smarts pattern.
Usage:
    python doublering_filter.py --
"""

import argparse
from openeye import oechem, oeomega
import pickle

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--smiles_database",
                        type=str,
                        required=True,
                        help=("Name of a SMILES file containing molecules "
                              "to search through."))
    args = parser.parse_args()
    ifs = oechem.oemolistream(args.smiles_database)
    ofs = oechem.oemolostream()

    emol_count = 0
    dr_count = 0
    group_num = 0
    total_count = 0
    mols = []

    SUBS = oechem.OESubSearch("[#6X3H1:1]~[#6X3:2](~[#6X3H1])-[#6X3:3](~[#6X3H1])~[#6X3H1:4]")
    for mol in ifs.GetOEGraphMols():
        try:
            emol_count += 1
            match = SUBS.Match(mol, True)
            if match.IsValid():
                mols.append(mol)
                dr_count += 1
            if dr_count == 1000:
                with open(f"group{group_num}_mols.pkl", "wb") as file:
                    pickle.dump(mols, file)
                total_count += dr_count
                print(f"total count at {total_count}")
                dr_count = 0
                group_num += 1
                mols = []
            if emol_count % 100 == 0:
                print(f"emol count at {emol_count}")
        except KeyboardInterrupt:
            break
        except:
            print(f"mol failed: {oechem.OEMolToSmiles(mol)}")
        
    with open(f"group{group_num}_mols.pkl", "wb") as file:
        pickle.dump(mols, file)
        total_count += dr_count

    print()
    print(f"Total emolecules: {emol_count}")
    print(f"Total molecules matching smarts: {total_count}")

if __name__ == "__main__":
    main()
