"""
Script for creating a histogram based on the WBO values calculated in wbo_calcs.py. The script
iterates through the openff_results folder, creating a large dataset of all the wbo values
from every group and plotting them all on a histogram. The result is saved as
emolecules_wbo_calcs.pdf.

Usage:
    wbo_histogram.py --file filename
    
Saving:
    Use the far right icon at the bottom of the pop up window to save the graph
"""

import argparse
import matplotlib.pyplot as plt
import os
import pickle

def create_histogram(filename):
    all_wbos = []
    
    with open(filename, "rb") as file:
        data = pickle.load(file)
        for smiles, wbos in data.items():
            if isinstance(wbos, list):
                all_wbos += wbos
    
    all_wbos[:] = [wbo for wbo in all_wbos if wbo != 0.0]
    
    plt.hist(all_wbos, log=True)
    plt.xlabel("AM1-WIBERG-ELF10 with Openeye")
    plt.ylabel("log(# of molecules)")
    
    if "protomer" in filename:
        plt.title("[#6X3H1:1]~[#6X3:2](~[#6X3H1])-[#6X3:3](~[#6X3H1])~[#6X3H1:4]\nmol protomer WBOs")
    elif "tautomer" in filename:
        plt.title("[#6X3H1:1]~[#6X3:2](~[#6X3H1])-[#6X3:3](~[#6X3H1])~[#6X3H1:4]\nmol tautomer WBOs")
    else:
        plt.title("[#6X3H1:1]~[#6X3:2](~[#6X3H1])-[#6X3:3](~[#6X3H1])~[#6X3H1:4]\nmol WBOs")
    
    if "conformer_count" in data and "protomer" in filename:
        plt.legend([f"{data['mol_count']} molecules\n{data['conformer_count']} protomers"])
    elif "conformer_count" in data and "tautomer" in filename:
        plt.legend([f"{data['mol_count']} molecules\n{data['conformer_count']} tautomers"])
    else:
        plt.legend([f"{data['mol_count']} molecules"])
    
    plt.show()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--file", type=str, required=True,
	help="name of the .pkl file to read data from")
    args = parser.parse_args()

    create_histogram(args.file)

if __name__ == "__main__":
    main()
