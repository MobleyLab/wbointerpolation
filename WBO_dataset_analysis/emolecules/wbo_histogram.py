"""
Script for creating a histogram based on a .pkl from wbo_results. The script iterates through the file, creating a large dataset of all the wbo values from every different molecule and plots them on a graph.

Usage:
    wbo_histogram.py --file filename --title title [--log True]
    
Saving:
    Use the far right icon at the bottom of the pop up window to save the graph
"""

import argparse
import matplotlib.pyplot as plt
import os
import pickle

def create_histogram(args):
    # Initialize variables to store all of the wbos in one list and count for the molecules
    all_wbos = []
    count = 0
    
    # Load the pickled data
    with open(args.file, "rb") as file:
        data = pickle.load(file)
        
        # Loop through the data, adding the wbo values form each {smiles: [wbo values]} to the all_wbos container
        for smiles, wbos in data.items():
            if isinstance(wbos, list):
                all_wbos += wbos
                count += 1    

    # Remove all 0.0 wbo values
    all_wbos[:] = [wbo for wbo in all_wbos if wbo != 0.0]

    # Plot the histogram with relevant labels, using log if the user specified --log True
    plt.hist(all_wbos, log=args.log)
    plt.xlabel("AM1-WIBERG-ELF10 with Openeye")

    # Sets the ylabel depending if it is a log graph or not
    if (args.log):
        plt.ylabel("log(# of molecules)")
    else:
        plt.ylabel("# of molecules")
    
    # Set the legend to use an additional conformer count if it is a protomers/tautomers file
    if "conformer_count" in data and "protomer" in args.file:
        plt.legend([f"{count} molecules\n{data['conformer_count']} protomers"])
    elif "conformer_count" in data and "tautomer" in args.file:
        plt.legend([f"{count} molecules\n{data['conformer_count']} tautomers"])
    else:
        plt.legend([f"{count} molecules"])

    # Set the title and show the plot so that it can be saved
    plt.title(args.title)    
    plt.show()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--file", type=str, required=True,
	    help="name of the .pkl file to read data from")
    
    parser.add_argument("--title", type=str, required=False, default="",
        help="title of the plot")
    parser.add_argument("--log", type=bool, required=False, default=False,
        help="Optional parameter to use log of data while constructing histogram")

    args = parser.parse_args()

    # Create a histogram using the supplied file and title
    create_histogram(args)

if __name__ == "__main__":
    main()
