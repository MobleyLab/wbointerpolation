"""
Script to compute the wbo calculations for molecules created by doublering_filter.py. The script
iterates through every file in the openff_results directory, creating two new .pkl files for
each group of molecules. For each group, molecules along with their wbo values are placed into the
openff_wbo_calcs directory and failures are tracked in the openff_wbo_fails directory.

Usage:
    wbo_calc.py
"""

import os
import pickle

def calc():
    group_num = 0
    smarts = "[#6X3H1:1]~[#6X3:2](~[#6X3H1])-[#6X3:3](~[#6X3H1])~[#6X3H1:4]"

    for subdir, dirs, files in os.walk("openff_results"):
        for file in files:
            if ".pkl" in file:
                wbos_dict = {}
                failed_mols = {}

                with open(f"openff_results/{file}", "rb") as file:
                    print(file)
                    smiles_dict = pickle.load(file)

                count = 0

                for smiles, off_molecule in smiles_dict.items():
                    count += 1
                    print(count)
                    try:
                        off_molecule.assign_fractional_bond_orders(bond_order_model="am1-wiberg-elf10")
                        matching_indices = off_molecule.chemical_environment_matches(smarts)
                    
                        wbos = []
                        ind_passed = []

                        for dihedrals in matching_indices:
                            if [dihedrals[1], dihedrals[2]] in ind_passed:
                                continue
                            elif [dihedrals[2], dihedrals[1]] in ind_passed:
                                continue
                            else:
                                bond = off_molecule.get_bond_between(dihedrals[1], dihedrals[2])
                                wbos.append(bond.fractional_bond_order)
                                wbos_dict[smiles] = wbos
                                ind_passed.append([dihedrals[1], dihedrals[2]])
                    except KeyboardInterrupt:
                        break
                    except Exception as e:
                        failed_mols[smiles] = e
                        continue
                        
                with open(f"openff_wbo_calcs/wbo_calc_group{group_num}.pkl", "wb") as file:
                    pickle.dump(wbos_dict, file)
                    
                with open(f"openff_wbo_fails/wbo_fails_group{group_num}.pkl", "wb") as file:
                    pickle.dump(failed_mols, file)
                
                group_num += 1

if __name__ == "__main__":
    calc()
