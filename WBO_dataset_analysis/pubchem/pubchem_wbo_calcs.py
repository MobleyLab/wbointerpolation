"""
Script to calculate the Wiberg Bond Order values for all filtered molecules in pubchem_results

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
    AM1_CALCULATOR = oequacpac.OEAM1()
    SUBS = oechem.OESubSearch("[#6X3H1:1]~[#6X3:2](~[#6X3H1])-[#6X3:3](~[#6X3H1])~[#6X3H1:4]")

    match = SUBS.Match(mol, True)
    wbos = []
    
    for match_base in match:
        bonds = []
        
        try:
            for bond in match_base.GetTargetBonds():
                bonds.append(bond)

            central_bond = bonds[2]
            central_bond_idxs = (central_bond.GetBgn().GetIdx(), central_bond.GetEnd().GetIdx())

            results = oequacpac.OEAM1Results()
            AM1_CALCULATOR.CalcAM1(results, mol)
            wbos.append(results.GetBondOrder(central_bond_idxs[0], central_bond_idxs[1]))

        except:
            continue
    
    print(wbos)

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
    with open(f"wbo_results/{outfile}_wbos.pkl", "wb") as file:
        for mol in ifs.GetOEGraphMols():
            smiles = oechem.OEMolToSmiles(mol)
            mol, _ = smiles2oemol(smiles)
            wbo_values = wbo_calc(mol)
            if wbo_values != []:
                wbos[smiles] = wbo_values
                count += 1
                print(count)

        wbos["mol_count"] = count
        pickle.dump(wbos, file)

if __name__ == "__main__":
    main()
