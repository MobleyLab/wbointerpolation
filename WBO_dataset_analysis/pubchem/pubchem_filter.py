"""
Script to iterate through the emolcules database looking for molecules that match the
"[#6X3H1:1]~[#6X3:2](~[#6X3H1])-[#6X3:3](~[#6X3H1])~[#6X3H1:4]" smarts pattern.

Usage:
    python pubchem_filter.py --sdf_database database.sdf
"""

from openeye import oechem, oeomega
import os
import sys

# Pattern to search
SUBS = oechem.OESubSearch("[#6X3H1:1]~[#6X3:2](~[#6X3H1])-[#6X3:3](~[#6X3H1])~[#6X3H1:4]")

def filter(filename):
    print(f"Filtering {filename}")
    total_count = 0
    match_count = 0
    fail_count = 0

    # Set up input and output streams
    ifs = oechem.oemolistream(f"pubchem_full_dataset/{filename}")
    ifs.SetFormat(oechem.OEFormat_SDF)
    ofs = oechem.oemolostream(f"pubchem_filtered/{filename.split('.')[0]}.oeb")
    ofs.SetFormat(oechem.OEFormat_OEB)
    
    # Iterate through every molecule, writing it to the output stream if it matches the substructure pattern
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

    print(f"Total emolecules for {filename}: {total_count}")
    print(f"Molecules found for {filename}: {match_count}")
    print(f"Failed molecules for {filename}: {fail_count}")
    print()
    
if __name__ == "__main__":
    print(f"Starting program")
    # Loop through every .sdf file in the pubchem_full_dataset directory and filter it
    for file in os.listdir("pubchem_full_dataset"):
        filename = os.fsdecode(file)
        if filename.endswith(".sdf"):
            filter(filename)
