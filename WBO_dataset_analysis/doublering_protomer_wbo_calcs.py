from openeye import oechem, oeomega, oequacpac
import pickle

def protomer_wbo_calc(mol):
    AM1_CALCULATOR = oequacpac.OEAM1()
    SUBS = oechem.OESubSearch("[#6X3H1:1]~[#6X3:2](~[#6X3H1])-[#6X3:3](~[#6X3H1])~[#6X3H1:4]")

    protomers = oequacpac.OEGetReasonableProtomers(mol)
    wbos = []

    for pro in protomers:
        match = SUBS.Match(pro, True)
        for match_base in match:
            bonds = []
            for bond in match_base.GetTargetBonds():
                bonds.append(bond)
            central_bond = bonds[2]
            central_bond_idxs = (central_bond.GetBgn().GetIdx(), central_bond.GetEnd().GetIdx())
        
        results = oequacpac.OEAM1Results()
        AM1_CALCULATOR.CalcAM1(results, pro)
        wbos.append(results.GetBondOrder(central_bond_idxs[0], central_bond_idxs[1]))

    return wbos

def main():
    ifs = oechem.oemolistream(f"oe_results/doublering_mols.oeb")
    ifs.SetFormat(oechem.OEFormat_OEB)

    ofs = oechem.oemolostream(f"segfault_mol.mol2")
    ofs.SetFormat(oechem.OEFormat_MOL2)
    wbos = {}
    count = 0

    with open("doublering_protomer_wbos.pkl", "wb") as file:
        for mol in ifs.GetOEGraphMols():
            smiles = oechem.OEMolToSmiles(mol)
            if smiles == "c1cc(cnc1)c2ccc(cc2)CC3(CCNCC3)C(=O)O":
                oechem.OEWriteMolecule(ofs, mol)   
            wbos[smiles] = protomer_wbo_calc(mol)
            print(count)
            count += 1
        pickle.dump(wbos, file)

if __name__ == "__main__":
    main()

