"""
Script to calculate the Wiberg Bond Order values for all protomers of the filtered molecules in oe_results/emolecules_mols.oeb.

Usage:
    python emolecules_wbo_calcs.py
"""


from openeye import oechem, oeomega, oequacpac
import pickle

def wbo_calc(mol):
    AM1_CALCULATOR = oequacpac.OEAM1()
    SUBS = oechem.OESubSearch("[#6X3H1:1]~[#6X3:2](~[#6X3H1])-[#6X3:3](~[#6X3H1])~[#6X3H1:4]")

    match = SUBS.Match(mol, True)
    wbos = []
    for match_base in match:
        bonds = []
        for bond in match_base.GetTargetBonds():
            bonds.append(bond)

        central_bond = bonds[2]
        central_bond_idxs = (central_bond.GetBgn().GetIdx(), central_bond.GetEnd().GetIdx())

    results = oequacpac.OEAM1Results()
    AM1_CALCULATOR.CalcAM1(results, mol)
    wbos.append(results.GetBondOrder(central_bond_idxs[0], central_bond_idxs[1]))

    return wbos

def main():
    ifs = oechem.oemolistream(f"oe_results/emolecules_mols.oeb")
    ifs.SetFormat(oechem.OEFormat_OEB)

    wbos = {}
    count = 0

    with open("emolecules_wbos.pkl", "wb") as file:
        for mol in ifs.GetOEGraphMols():
            smiles = oechem.OEMolToSmiles(mol)
            wbos[smiles] = wbo_calc(mol)
            print(count)
            count += 1

        wbos["mol_count"] = count
        pickle.dump(wbos, file)

if __name__ == "__main__":
    main()
