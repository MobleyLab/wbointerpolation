from openeye import oechem, oeomega, oequacpac
import pickle

def tautomer_wbo_calc(mol):
    AM1_CALCULATOR = oequacpac.OEAM1()
    SUBS = oechem.OESubSearch("[#6X3H1:1]~[#6X3:2](~[#6X3H1])-[#6X3:3](~[#6X3H1])~[#6X3H1:4]")

    tautomers = oequacpac.OEGetReasonableTautomers(mol)
    wbos = []

    for tautomer in tautomers:
        try:
            tautomer.Sweep()
            match = SUBS.Match(tautomer, True)
            for match_base in match:
                bonds = []
                for bond in match_base.GetTargetBonds():
                    bonds.append(bond)
                central_bond = bonds[2]
                central_bond_idxs = (central_bond.GetBgn().GetIdx(), central_bond.GetEnd().GetIdx())

            results = oequacpac.OEAM1Results()
            AM1_CALCULATOR.CalcAM1(results, tautomer)
            wbos.append(results.GetBondOrder(central_bond_idxs[0], central_bond_idxs[1]))
        except:
            continue

    return wbos

def main():
    ifs = oechem.oemolistream(f"oe_results/doublering_mols.oeb")
    ifs.SetFormat(oechem.OEFormat_OEB)

    wbos = {}
    count = 0
    conformer_count = 0

    with open("doublering_tautomer_wbos.pkl", "wb") as file:
        for mol in ifs.GetOEGraphMols():
            #print(count)
            print(conformer_count)
            smiles = oechem.OEMolToSmiles(mol)
            tautomer_wbos = tautomer_wbo_calc(mol)
            wbos[smiles] = tautomer_wbos
            conformer_count += len(tautomer_wbos) 
            count += 1
            
        wbos["mol_count"] = count
        wbos["conformer_count"] = conformer_count
        pickle.dump(wbos, file)

if __name__ == "__main__":
    main()

