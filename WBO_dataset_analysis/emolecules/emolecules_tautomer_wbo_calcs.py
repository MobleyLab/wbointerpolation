"""
Script to calculate the Wiberg Bond Order values for all tautomers of the filtered molecules in oe_results/emolecules_mols.oeb.

Usage:
    python emolecules_tautomer_wbo_calcs.py
"""

from openeye import oechem, oeomega, oequacpac
import pickle

def tautomer_wbo_calc(mol):
    # Set up tools for calculation
    AM1_CALCULATOR = oequacpac.OEAM1()
    SUBS = oechem.OESubSearch("[#6X3H1:1]~[#6X3:2](~[#6X3H1])-[#6X3:3](~[#6X3H1])~[#6X3H1:4]")

    # Get the tautomers for the molecule
    tautomers = oequacpac.OEGetReasonableTautomers(mol)
    wbos = []

    # Loop through the tautomers
    for tautomer in tautomers:
        try:
            # Necessary setup to get data from the tautomer
            tautomer.Sweep()
            # Create a match object using the SUBS substructure and list to hold WBO values
            match = SUBS.Match(tautomer, True)
            
            # Loop through every instance of a match
            for match_base in match:
                bonds = []
                # Get all of the bonds within the match, adding them to a list of bonds
                for bond in match_base.GetTargetBonds():
                    bonds.append(bond)
                
                # Store the central bond and its indexes
                central_bond = bonds[2]
                central_bond_idxs = (central_bond.GetBgn().GetIdx(), central_bond.GetEnd().GetIdx())

            # Calculate the result using the central bond indexes and the AM1 Calculator
            results = oequacpac.OEAM1Results()
            AM1_CALCULATOR.CalcAM1(results, tautomer)
            wbos.append(results.GetBondOrder(central_bond_idxs[0], central_bond_idxs[1]))
        except:
            continue

    # Retun the list of WBO values for the molecule, 
    # done as a list because there can be more than one instance of a match for each molecule
    return wbos

def main():
    # Set up the input stream
    ifs = oechem.oemolistream(f"oe_results/emolecules_mols.oeb")
    ifs.SetFormat(oechem.OEFormat_OEB)

    # Initialize variables to store wbos, a molecule count, and count for all the conformers
    wbos = {}
    count = 0
    conformer_count = 0

    # For every mol, calculate the wbo values and save it to the dictionary of {smiles: [wbo values]}
    with open("wbo_results/emolecules_tautomer_wbos.pkl", "wb") as file:
        for mol in ifs.GetOEGraphMols():
            print(conformer_count)
            smiles = oechem.OEMolToSmiles(mol)
            tautomer_wbos = tautomer_wbo_calc(mol)
            wbos[smiles] = tautomer_wbos
            conformer_count += len(tautomer_wbos) 
            count += 1
            
        wbos["mol_count"] = count
        wbos["conformer_count"] = conformer_count
        # Store the data using pickle dump into the set outfile
        pickle.dump(wbos, file)

if __name__ == "__main__":
    main()

