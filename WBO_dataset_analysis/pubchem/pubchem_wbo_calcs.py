"""
Script to calculate the Wiberg Bond Order values for all filtered molecules in pubchem_filtered

Usage:
    python pubchem_wbo_calcs.py --oeb_file filename
"""

from openeye import oechem, oeomega, oequacpac
import argparse
import json
import pickle

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

def wbo_calc(mol):
    # Set up tools for calculation
    AM1_CALCULATOR = oequacpac.OEAM1()
    SUBS = oechem.OESubSearch("[#6X3H1:1]~[#6X3:2](~[#6X3H1])-[#6X3:3](~[#6X3H1])~[#6X3H1:4]")

    # Create a match object using the SUBS substructure and list to hold WBO values
    match = SUBS.Match(mol, True)
    wbos = []
    
    # Loop through every instance of a match
    for match_base in match:
        bonds = []
        
        # Get all of the bonds within the match, adding them to a list of bonds
        try:
            for bond in match_base.GetTargetBonds():
                bonds.append(bond)

            # Store the central bond and its indexes
            central_bond = bonds[2]
            central_bond_idxs = (central_bond.GetBgn().GetIdx(), central_bond.GetEnd().GetIdx())

            # Calculate the result using the central bond indexes and the AM1 Calculator
            results = oequacpac.OEAM1Results()
            AM1_CALCULATOR.CalcAM1(results, mol)
            wbos.append(results.GetBondOrder(central_bond_idxs[0], central_bond_idxs[1]))

        except:
            continue
    
    print(wbos)

    # Retun the list of WBO values for the molecule, 
    # done as a list because there can be more than one instance of a match for each molecule
    return wbos

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--oeb_file",
                        type=str,
                        required=True,
                        help=("Name of an OEB file containing molecules "
                              "to search through."))
    args = parser.parse_args()

    wbos = {}
    count = 0

    ifs = oechem.oemolistream(args.oeb_file)
    ifs.SetFormat(oechem.OEFormat_OEB)

    outfile = args.oeb_file.split("/")[-1][:-4]

    # For every mol, calculate the wbo values and save it to the dictionary of {smiles: [wbo values]}
    with open(f"wbo_results_pubchem_full/{outfile}_wbos.pkl", "wb") as file:
        for mol in ifs.GetOEGraphMols():
            smiles = oechem.OEMolToSmiles(mol)
            mol, _ = smiles2oemol(smiles)
            wbo_values = wbo_calc(mol)
            if wbo_values != []:
                wbos[smiles] = wbo_values
                count += 1
                print(count)

        wbos["mol_count"] = count
        
        # Store the data using pickle dump into the set outfile
        pickle.dump(wbos, file)

if __name__ == "__main__":
    main()
