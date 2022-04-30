"""
Script to iterate through the emolcules database looking for molecules that match the
"[#6X3H1:1]~[#6X3:2](~[#6X3H1])-[#6X3:3](~[#6X3H1])~[#6X3H1:4]" smarts pattern.

Usage:
    python emolecules_filter.py --smiles_database database.smi --group number
"""

import argparse
from openeye import oechem, oeomega
from openff.toolkit.topology import Molecule
import sys

def smiles2oemol(smiles):
    """
    Conforms molecules from a smiles string to an OEMol with all the traits
    necessary for OpenEye functions to return accurate data
    """
    mol = oechem.OEMol()
    oechem.OESmilesToMol(mol, smiles)
    #initialize omega
    omega = oeomega.OEOmega()
    omega.SetMaxConfs(1)
    omega.SetIncludeInput(True)
    omega.SetCanonOrder(True)
    omega.SetSampleHydrogens(True)
    omega.SetStrictStereo(False)
    omega.SetStrictAtomTypes(True)
    omega.SetIncludeInput(False)
    status = omega(mol)
    
    return (mol, status)

def filter():
    print("Made it to filter")
    parser = argparse.ArgumentParser()
    parser.add_argument("--smiles_database",
                        type=str,
                        required=True,
                        help=("Name of a SMILES file containing molecules "
                              "to search through."))
    parser.add_argument("--group",
                        type=str,
                        required=True,
                        help=("Group number of process"))
    args = parser.parse_args()
    emol_count = 0
    match_count = 0
    fail_count = 0
    group_num = args.group
    SUBS = oechem.OESubSearch("[#6X3H1:1]~[#6X3:2](~[#6X3H1])-[#6X3:3](~[#6X3H1])~[#6X3H1:4]")
    ifs = oechem.oemolistream(args.smiles_database)
    ifs.SetFormat(oechem.OEFormat_SMI)
    ofs = oechem.oemolostream(f"oe_results/mols_group{group_num}.oeb")
    ofs.SetFormat(oechem.OEFormat_OEB)
    ofs_safety = oechem.oemolostream(f"results/mols_group{group_num}.smi")
    ofs_safety.SetFormat(oechem.OEFormat_SMI)
    ofs_fails = oechem.oemolostream(f"oe_fails/failed_group{group_num}.oeb")
    ofs_fails.SetFormat(oechem.OEFormat_OEB)
    print(f"starting group {group_num}")
    
    for mol in ifs.GetOEGraphMols():
        emol_count += 1
        print(emol_count)
        oechem.OEWriteMolecule(ofs_safety, mol)
        smiles = oechem.OEMolToSmiles(mol)
        try:
            molecule, status = smiles2oemol(smiles)
            match = SUBS.Match(molecule, True)
            if match.IsValid():
                oechem.OEWriteMolecule(ofs, molecule)
                match_count += 1
        except KeyboardInterrupt:
            break
        except Exception as e:
            oechem_fails.OEWriteMolecule(ofs_fails, molecule)
            fail_count += 1
            continue

    print(f"Total emolecules: {emol_count}")
    print(f"Molecules found: {match_count}")
    print(f"Failed molecules: {fail_count}")
    
if __name__ == "__main__":
    print(f"Starting program")
    filter()
