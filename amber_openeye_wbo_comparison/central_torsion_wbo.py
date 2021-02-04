
import cmiles
import collections
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pickle
import qcportal as ptl
import time
 
from openeye import oechem, oeomega, oequacpac
from openforcefield.topology import Molecule, Topology
from openforcefield.typing.engines.smirnoff import ForceField
from single_conformer import conform_molecules

#Setup
client = ptl.FractalClient()
torsion_datasets = client.list_collections("TorsionDriveDataset")
datasets = []
for i in range(len(torsion_datasets)):
    datasets.append(torsion_datasets.index[i][1])

def get_data(datasets):
    """Loads data from a QC archive dataset"""
    data_dict={}
    for dataset_name in datasets:
        count = 0
        while True:
            try:
                ds = client.get_collection("TorsionDriveDataset", dataset_name)
                ds.status("default", status="COMPLETE")
                break
            except:
                time.sleep(20)
                count += 1
                if count < 2:
                    continue
                else:
                    break

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
    
def visualize_benchmark(group_num, wbo_values):
    """Creates a double bar graph to visualize the wbo value comparisons"""
    width = .35
    
    fig, ax = plt.subplots()
    openff_wbos = []
    openeye_wbos = []
    
    for mol, wbos in wbo_values.items():
        openff_wbos.append(wbos[0])
        openeye_wbos.append(wbos[1])
    
    x = np.arange(len(wbo_values.keys()))
    
    ax.bar(x - width/2, openff_wbos, color = "#236AB9", width = width, label="OpenFF wbo calc")
    ax.bar(x + width/2, openeye_wbos, color = "#64A9F7", width = width, label="Openeye wbo calc") #64A9F7
    
    ax.set_ylabel("Wiberg Bond Order", fontweight = "bold")
    ax.set_xlabel("Molecules", fontweight = "bold")
    ax.set_xticks(x)
    #ax.set_xticklabels(wbo_values.keys(), rotation = 90) #Label for each molecule
    ax.set_title(f"WBO Benchmark Group Number {group_num}")

    ax.legend(["OpenFF wbo calc", "Openeye wbo calc"])
    plt.savefig(f"QCA_WBO_figures/QCA_WBO_benchmark Group Number{group_num}.png", bbox_inches = "tight")
    
def main():
    dataset_substitutedphenyl = ['OpenFF Substituted Phenyl Set 1']
    data = get_data(dataset_substitutedphenyl)
    
    #conform_molecules(data, dataset_substitutedphenyl[0].replace(" ", "")) #can be commented out to save time if results/ files alreadys exist
    
    benchmark_data = []
    wbo_values = {}
    with open("results/OpenFFSubstitutedPhenylSet1-ambertools.pkl", "rb") as amber_file:
        with open("results/OpenFFSubstitutedPhenylSet1-openeye.pkl", "rb") as openeye_file:
            amber_data = pickle.load(amber_file)
            openeye_data = pickle.load(openeye_file)
            
            count = 0
            #Iterate through the amber and openeye version of each molecule
            for amber, openeye in zip(amber_data, openeye_data):
                #amber_mol = smiles2oemol(amber[0].to_smiles())
                #Using WBO as provided by the fractional bond order between central torsion indices
                #amber_wbo = amber[1]
                #Using central torsion indices to calculate WBO through OpenEye
                #amber_torsions = amber[2]
                
                openeye_mol = smiles2oemol(openeye[0].to_smiles())
                #Using WBO as provided by the fractional bond order between central torsion indices from openFF
                openeye_wbo = openeye[1]
                #Using central torsion indices to calculate WBO through OpenEye
                openeye_torsions = openeye[2]
                #Using WBO as provided by the fractional bond order between central torsion indices and WBO caclulated through openeye to compare openFF vs openeye
                wbo_values[ openeye_mol ] = ( openeye_wbo, wiberg_bond_order(openeye_mol, openeye_torsions) )
                
                #Using central torsion indices to calculate WBO through OpenEye
                #wbo_values[ (amber_mol, openeye_mol) ] = ( wiberg_bond_order(amber_mol, amber_torsions), wiberg_bond_order(openeye_mol, openeye_torsions) )
                
                #Groups the data into sets of 25 for better visualization
                count += 1
                if count % 25 == 0 or count == len(data):
                    benchmark_data.append(wbo_values)
                    wbo_values = {}
                    group_count = 0
    
    #Creates a visual of the wbo values for comparison
    for group_num, wbo_values in enumerate(benchmark_data):
        visualize_benchmark(group_num+1, wbo_values)
    
if __name__ == "__main__":
    main()
