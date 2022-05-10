"""
Script to iterate through the emolcules database looking for molecules that match the
"[#6X3H1:1]~[#6X3:2](~[#6X3H1])-[#6X3:3](~[#6X3H1])~[#6X3H1:4]" smarts pattern.

Usage:
    python emolecules_filter.py --sdf_database database.sdf
"""

import argparse
from openeye import oechem, oeomega
import sys

def filter():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sdf_database",
                        type=str,
                        required=True,
                        help=("Name of a SDF file containing molecules "
                              "to search through."))
    args = parser.parse_args()
    
    total_count = 0
    match_count = 0
    fail_count = 0

    ifs = oechem.oemolistream(args.sdf_database)
    ifs.SetFormat(oechem.OEFormat_SDF)
    ofs = oechem.oemolostream(f"pubchem_filtered/{args.sdf_database.split('.')[0]}.oeb")
    ofs.SetFormat(oechem.OEFormat_OEB)

    SUBS = oechem.OESubSearch("[#6X3H1:1]~[#6X3:2](~[#6X3H1])-[#6X3:3](~[#6X3H1])~[#6X3H1:4]")
    
    for mol in ifs.GetOEGraphMols():
        total_count += 1
        print(total_count)
    
        try:
            match = SUBS.Match(mol, True)
            if match.IsValid():
                oechem.OEWriteMolecule(ofs, mol)
                match_count += 1
        except KeyboardInterrupt:
            break
        except Exception as e:
            fail_count += 1
            continue

    print(f"Total emolecules: {total_count}")
    print(f"Molecules found: {match_count}")
    print(f"Failed molecules: {fail_count}")
    
if __name__ == "__main__":
    print(f"Starting program")
    filter()
