
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
    
def visualize_bargraphs(benchmark_data):
    for group_num, wbo_values in enumerate(benchmark_data):
        create_bargraph(group_num+1, wbo_values)

def visualize_scatterplot(benchmark_data):
    all_wbo_values = {}
    for wbo_values in benchmark_data:
        all_wbo_values.update(wbo_values)

    create_scatterplot(all_wbo_values)

def create_bargraph(group_num, wbo_values):
    """Creates a double bar graph to visualize the wbo value comparisons"""
    width = .35
    
    fig, ax = plt.subplots()
    amber_wbos = []
    openeye_wbos = []
    
    for mol, wbos in wbo_values.items():
        amber_wbos.append(wbos[0])
        openeye_wbos.append(wbos[1])
    
    x = np.arange(len(wbo_values.keys()))
    
    ax.bar(x - width/2, amber_wbos, color = "#236AB9", width = width, label="Ambertools")
    ax.bar(x + width/2, openeye_wbos, color = "#64A9F7", width = width, label="OpenEye") #64A9F7
    
    ax.set_ylabel("Wiberg Bond Order", fontweight = "bold")
    ax.set_xlabel("Molecules", fontweight = "bold")
    ax.set_xticks(x)
    #Label for each molecule
    ax.set_xticklabels(wbo_values.keys(), rotation = 90)
    ax.set_title(f"WBO Benchmark Group Number {group_num}")

    ax.legend(["Ambertools", "OpenEye"])
    plt.savefig(f"QCA_WBO_bargraphs/QCA_WBO_benchmark Group Number {group_num}.png", bbox_inches = "tight")
    
def create_scatterplot(wbo_values):
    """Creates a scatterplot to visualize the wbo value comparisons"""
    fig, ax = plt.subplots()
    amber_wbos = []
    openeye_wbos = []
    
    for mol, wbos in wbo_values.items():
        amber_wbos.append(wbos[0])
        openeye_wbos.append(wbos[1])
    
    ax.scatter(amber_wbos, openeye_wbos, color = "#236AB9")
    
    ax.set_xlabel("Ambertools", fontweight = "bold")
    ax.set_ylabel("OpenEye", fontweight = "bold")
    
    ax.set_title("WBO Benchmark Scatterplot")

    plt.savefig("QCA_WBO_scatterplot/QCA_WBO_benchmark.png", bbox_inches = "tight")

def find_notable_differences(benchmark_data):
    all_wbo_values = {}
    for wbo_values in benchmark_data:
        all_wbo_values.update(wbo_values)
        
    with open("analysis/QCA_WBO_noteworthy_differences.txt", "w") as file:
        for smiles, wbos in sorted(all_wbo_values.items(),
                                   key = lambda x: abs(x[1][0]-x[1][1]))[round(.75*len(all_wbo_values)):]:
            file.write(f"Smiles: {smiles}\n")
            file.write(f"Amber wbo: {wbos[0]}\n")
            file.write(f"OpenEye wbo: {wbos[1]}\n")
            file.write(f"Difference: {abs(wbos[0]-wbos[1])}\n")
            file.write("\n")

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
                amber_mol = amber[0]
                #amber_smiles = amber[0].to_smiles()
                amber_torsions = amber[1]
                
                openeye_mol = openeye[0]
                #openeye_smiles = openeye[0].to_smiles()
                openeye_torsions = openeye[1]
                
                smiles = amber[0].to_smiles()
                
                #Using WBO as provided by the fractional bond order between central torsion indices
                wbo_values[smiles] = ( amber_mol.get_bond_between(amber_torsions[0], amber_torsions[1]).fractional_bond_order, openeye_mol.get_bond_between(openeye_torsions[0], openeye_torsions[1]).fractional_bond_order )
                
                #Using central torsion indices to calculate WBO through OpenEye
                #wbo_values[ (amber_smiles, amber_smiles) ] = ( wiberg_bond_order(amber_mol, amber_torsions), wiberg_bond_order(openeye_mol, openeye_torsions) )
                
                #####IGNORE: SEPARATE TEST ON OPENFF VS OPENEYE WBO CALCULATION#####
                #Using WBO as provided by the fractional bond order between central torsion indices and WBO caclulated through openeye to compare openFF vs openeye
                #wbo_values[ openeye_mol ] = ( openeye_wbo, wiberg_bond_order(openeye_mol, openeye_torsions) )
                
                #Groups the data into sets of 25 for better visualization
                count += 1
                if count % 25 == 0 or count == len(data):
                    benchmark_data.append(wbo_values)
                    wbo_values = {}
    
    #Creates a file that includes the top 25% of differences between Ambertools and OpenEye wbos
    find_notable_differences(benchmark_data)
    
    #Create bargraphs in groups of 25 comparing the Ambertools and OpenEye wbos
    #visualize_bargraphs(benchmark_data)
    
    #Create a scatterplot comparing the Ambertools and OpenEye wbos
    visualize_scatterplot(benchmark_data)
    
if __name__ == "__main__":
    main()
