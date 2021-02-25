"""
Script to calculate the WBOs for central torsions of an Ambertoolkit
and OpenEye molecule to compare the difference between the conformers

Usage:
    python central_torsion_wbo.py
    
Ouput:
    A .pkl file containing a list of dictionaries of every mol's smiles string mapped
    to a tuple of the Ambertools WBO and OpenEye WBO
    
    The list is grouped into dictionaries of length 25 to more easily create plots
    to represent the data
"""

import cmiles
import collections
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os
import pickle
import qcportal as ptl
import time

from single_conformer import conform_molecules
from openeye import oechem, oeomega

#Setup
client = ptl.FractalClient()
torsion_datasets = client.list_collections("TorsionDriveDataset")
datasets = []
for i in range(len(torsion_datasets)):
    datasets.append(torsion_datasets.index[i][1])

def get_data(dataset_name):
    """Loads data from a QC archive dataset"""
    data_dict={}
    #for dataset_name in datasets:
    count = 0
    while True:
        try:
            ds = client.get_collection("TorsionDriveDataset", dataset_name)
            ds.status("default", status="COMPLETE")
            break
        except:
            time.sleep(20)
            print(f"failed to get dataset {dataset_name}")
            return data_dict
#            count += 1
#            if count < 2:
#                continue
#            else:
#                break
    params = []
    for index in ds.df.index:
        # get the dihedral indices
        dihedral_indices = ds.df.loc[index].default.keywords.dihedrals[0]
        cmiles=ds.get_entry(index).attributes['canonical_isomeric_explicit_hydrogen_mapped_smiles']
        data_dict[smiles2oemol(cmiles)] = dihedral_indices[1:3]

    counter = collections.Counter(params)
    print(dataset_name, counter)
    print(" ")
    
    return data_dict

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
    omega(mol)
    
    return mol

def wiberg_bond_order(mol, bond_idxs):
    """Calculates the Wiberg Bond Order for a specified bond"""
    #Construction for WBO calculation
    AM1_CALCULATOR = oequacpac.OEAM1()
    results = oequacpac.OEAM1Results()
    
    AM1_CALCULATOR.CalcAM1(results, mol)
    
    return results.GetBondOrder(bond_idxs[0], bond_idxs[1])

def main():
    """
    Creates the conformers for the molecules of a given dataset and the dataset
    to compare the difference between the Ambertools WBO and OpenEye WBO
    """
    
    for dataset_name in datasets:
        data = get_data(dataset_name)
        dataset_file_name = dataset_name.replace(" ", "")
        
        conformed = False
        for file_name in os.listdir("conformer_results"):
            if dataset_file_name in file_name:
                conformed = True
                break
        if not conformed:
            print(f"Conforming dataset: {dataset_name}")
            conform_molecules(data, dataset_file_name)

        benchmark_data = []
        wbo_values = {}
        with open(f"conformer_results/{dataset_file_name}-ambertools.pkl", "rb") as amber_file, open(f"conformer_results/{dataset_file_name}-openeye.pkl", "rb") as openeye_file:
            amber_data = pickle.load(amber_file)
            openeye_data = pickle.load(openeye_file)

            count = 0
            #Iterate through the amber and openeye version of each molecule
            for amber, openeye in zip(amber_data, openeye_data):
                amber_mol = amber[0]
                #amber_smiles = amber[0].to_smiles()
                amber_torsions = amber[1]

                openeye_mol = openeye[0]
                #openeye_smiles = openeye[0].to_smiles()
                openeye_torsions = openeye[1]

                smiles = amber[0].to_smiles()

                #Using WBO as provided by the fractional bond order between central torsion indices
                wbo_values[smiles] = ( amber_mol.get_bond_between(amber_torsions[0], amber_torsions[1]).fractional_bond_order, openeye_mol.get_bond_between(openeye_torsions[0], openeye_torsions[1]).fractional_bond_order )

                #Groups the data into sets of 25 for better visualization
                count += 1
                if count % 25 == 0 or count == len(data):
                    benchmark_data.append(wbo_values)
                    wbo_values = {}

        benchmark = (dataset_name, benchmark_data)
        with open(f"benchmark_results/{dataset_file_name}.pkl", "wb") as file:
            pickle.dump(benchmark, file)
        
        length = 0
        for bd in benchmark_data:
            length += len(bd)
            
        print(f"Completed benchmark data for: {dataset_name} - length {length}")
    
if __name__ == "__main__":
    main()
