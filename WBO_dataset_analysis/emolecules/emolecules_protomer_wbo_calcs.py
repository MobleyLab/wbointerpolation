"""
Script to calculate the Wiberg Bond Order values for all protomers of the filtered molecules in oe_results/emolecules_filtered.oeb.

Usage:
    python emolecules_protomer_wbo_calcs.py
"""

from openeye import oechem, oeomega, oequacpac
import pickle

def protomer_wbo_calc(mol):
    # Set up tools for calculation
    AM1_CALCULATOR = oequacpac.OEAM1()
    SUBS = oechem.OESubSearch("[#6X3H1:1]~[#6X3:2](~[#6X3H1])-[#6X3:3](~[#6X3H1])~[#6X3H1:4]")

    # Get the protomers for the molecule
    protomers = oequacpac.OEGetReasonableProtomers(mol)
    wbos = []

    # Loop through the protomers
    for protomer in protomers:
        try:
            # Necessary setup to get data from the protomer
            protomer.Sweep()
            # Create a match object using the SUBS substructure and list to hold WBO values
            match = SUBS.Match(protomer, True)

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
            AM1_CALCULATOR.CalcAM1(results, protomer)
            wbos.append(results.GetBondOrder(central_bond_idxs[0], central_bond_idxs[1]))
        except:
            continue

    # Retun the list of WBO values for the molecule, 
    # done as a list because there can be more than one instance of a match for each molecule
    return wbos


def main():
    # Set up the input stream
    ifs = oechem.oemolistream(f"oe_results/emolecules_filtered.oeb")
    ifs.SetFormat(oechem.OEFormat_OEB)

    # Initialize variables to store wbos, a molecule count, and count for all the conformers
    wbos = {}
    count = 0
    conformer_count = 0

    # For every mol, calculate the wbo values and save it to the dictionary of {smiles: [wbo values]}
    with open("wbo_results/emolecules_protomer_wbos.pkl", "wb") as file:
        for mol in ifs.GetOEGraphMols():
            #print(count)
            print(conformer_count)
            smiles = oechem.OEMolToSmiles(mol)
            protomer_wbos = protomer_wbo_calc(mol)
            wbos[smiles] = protomer_wbos
            conformer_count += len(protomer_wbos)
            count += 1

        wbos["mol_count"] = count
        wbos["conformer_count"] = conformer_count
        # Store the data using pickle dump into the set outfile
        pickle.dump(wbos, file)

if __name__ == "__main__":
    main()

