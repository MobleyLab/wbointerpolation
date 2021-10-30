"""
Script to iterate through the emolcules database looking for molecules that match the
"[#6X3H1:1]~[#6X3:2](~[#6X3H1])-[#6X3:3](~[#6X3H1])~[#6X3H1:4]" smarts pattern.

Usage:
    python doublering_filter.py --smiles_database database.smi
"""

import argparse
from openeye import oechem, oeomega
from openff.toolkit.topology import Molecule
import pickle
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
    parser = argparse.ArgumentParser()
    parser.add_argument("--smiles_database",
                        type=str,
                        required=True,
                        help=("Name of a SMILES file containing molecules "
                              "to search through."))
    args = parser.parse_args()
    
    emol_count = 0
    match_count = 0
    total_match = 0
    fail_count = 0
    total_fail = 0
    group_num = 0
    failed_group_num = 0
    molecules_match = {}
    fails_openff_mol = []
    SUBS = oechem.OESubSearch("[#6X3H1:1]~[#6X3:2](~[#6X3H1])-[#6X3:3](~[#6X3H1])~[#6X3H1:4]")
    
    with open(args.smiles_database, "r") as file:
        for line in file:
            emol_count += 1
            print(emol_count)
            
            smiles = line.split()[0]
            try:
                molecule, status = smiles2oemol(smiles)
                match = SUBS.Match(molecule, True)
                if match.IsValid():
                    molecules_match[smiles] = molecule
                    match_count += 1
                if match_count == 10000:
                    with open(f"oe_results/mols_group{group_num}.pkl", "wb") as handle:
                        pickle.dump(molecules_match, handle, protocol=pickle.HIGHEST_PROTOCOL)
                    group_num += 1
                    molecules_match = {}
                    total_match += match_count
                    match_count = 0
            except KeyboardInterrupt:
                break
            except Exception as e:
                fails_openff_mol.append( (smiles, str(e)) )
                fail_count += 1
                if fail_count == 10000:
                    with open(f"oe_fails/failed_mols_group{failed_group_num}.pkl", 'wb') as handle:
                        pickle.dump(fails_openff_mol, handle, protocol=pickle.HIGHEST_PROTOCOL)
                    failed_group_num += 1
                    fails_openff_mol = {}
                    total_fail += fail_count
                    fail_count = 0
                continue
                
    with open(f"oe_results/mols_group{group_num}.pkl", "wb") as handle:
        pickle.dump(molecules_match, handle, protocol=pickle.HIGHEST_PROTOCOL)
     
    
    with open(f"oe_fails/failed_mols_group{failed_group_num}.pkl", "wb") as handle:
        pickle.dump(fails_openff_mol, handle, protocol=pickle.HIGHEST_PROTOCOL)
    
    total_match += match_count
    total_fail += fail_count
    
    print(f"Total emolecules: {emol_count}")
    print(f"Molecules found: {total_match}")
    print(f"Failed molecules: {total_fail}")
    
if __name__ == "__main__":
    filter()
