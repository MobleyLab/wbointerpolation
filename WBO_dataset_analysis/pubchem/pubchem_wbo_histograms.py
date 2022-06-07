"""
Script for creating a histogram based on the .pkl files from wbo_results in the pubchem directory. 

The script has two options: 
    1. Iterate through the entire directory, creating an individual graph for the data in each .pkl file
    2. Iterate through the entire directory, creating a collective graph based on the data in all of the .pkl files combined

The option is selected using the following command line parameters. The --one_plot parameter is default False which does option 1 for the script, 
users can toggle one_plot to True which does option 2 for the script. The --log parameter uses the log of the results when creating the histogram, 
this parameter is also default value False.

Usage:
    python pubchem_wbo_histogram.py --in_dir dir --out_dir dir [--log True] [--one_plot True] [--title title]
"""

import argparse
import matplotlib.pyplot as plt
import os
import pickle

def create_database_histogram(args):
    # Initialize variables to store all of the wbos in one list and counts for molecules and each file
    all_wbos = []
    count = 0
    file_count = 0
    
    # Loop through each file in the desired input directory
    for file in os.listdir(args.in_dir):
        filename = os.fsdecode(file)

        if filename.endswith(".pkl"):
            # Load the last pickle dump in the file since data was saved multiple times while running with slurm to prevent 
            # loss of data in case the job could not complete
            with open(f"{args.in_dir}/{filename}", "rb") as file:
                try:
                    data = pickle.load(file)
                    while (data != None):
                        try:
                            data = pickle.load(file)
                        except Exception as e:
                            break
                except Exception as e:
                    continue
                
                # Loop through the data, adding the wbo values form each {smiles: [wbo values]} to the all_wbos container
                for smiles, wbos in data.items():
                    if isinstance(wbos, list):
                        all_wbos += wbos
                        count += 1  
            
            file_count += 1
            print(file_count)  

    # Remove all 0.0 wbo values
    all_wbos[:] = [wbo for wbo in all_wbos if wbo != 0.0]

    # Plot the histogram with relevant labels, using log if the user specified --log True
    plt.hist(all_wbos, log=args.log)
    plt.xlabel("AM1-WIBERG-ELF10 with Openeye")
    plt.legend([f"{count} molecules"])
    
    # Sets the title depending if the user suplied one with --title
    if args.title:
        plt.title(args.title)
    else:
        plt.title(f"All WBOs for {args.in_dir}")

    # Sets the ylabel depending if it is a log graph or not
    if (args.log):
        plt.ylabel("log(# of molecules)")
    else:
        plt.ylabel("# of molecules")
    
    # Save the figure as all_wbos to show it represents all the wbos for the given input directory  
    plt.savefig(f"{args.out_dir}/all_wbos.png")

def create_single_histograms(args):
    # Initialize variables to store all of the wbos in one list and counts for molecules and each file
    all_wbos = []
    count = 0
    file_count = 0
    
    # Loop through each file in the desired input directory
    for file in os.listdir(args.in_dir):
        filename = os.fsdecode(file)
        # Truncate the filename by removing the .pkl extension and "_wbos" at the end so the name can be used elsewhere
        filename_as_title = filename.split(".")[0][:-5]

        if filename.endswith(".pkl"):
            
            # Load the last pickle dump in the file since data was saved multiple times while running with slurm to prevent 
            # loss of data in case the job could not complete
            with open(f"{args.in_dir}/{filename}", "rb") as file:
                try:
                    data = pickle.load(file)
                    while (data != None):
                        try:
                            data = pickle.load(file)
                        except Exception as e:
                            break
                except Exception as e:
                    continue

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
            plt.legend([f"{count} molecules"])
            plt.title(f"{filename_as_title} WBOs")

            # Sets the ylabel depending if it is a log graph or not
            if (args.log):
                plt.ylabel("log(# of molecules)")
            else:
                plt.ylabel("# of molecules")
            
            # Save the figure using the formatted filename from above
            plt.savefig(f"{args.out_dir}/{filename_as_title}_wbos.png")

            # Clear the plot and all_wbos list so that the next file will have a clean start
            plt.clf()
            all_wbos = []
            
            # File count for debugging purposes
            file_count += 1
            print(file_count)

def main():
    parser = argparse.ArgumentParser()
    
    parser.add_argument("--in_dir", type=str, required=True,
        help="Directory of .pkl files to create histograms based on")
    parser.add_argument("--out_dir", type=str, required=True,
        help="Directory to store saved histograms")
    
    parser.add_argument("--log", type=bool, required=False, default=False,
        help="Optional parameter to use log of data while constructing histogram")
    parser.add_argument("--one_plot", type=bool, required=False, default=False,
        help="Optional parameter to create a histogram using all files in wbo_results in one plot")
    parser.add_argument("--title", type=str, required=False,
        help="Optional parameter to set the title for a histogram, should only be used with --one_plot set True")
    
    args = parser.parse_args()

    # Call option 1 or option 2 from the script description depending on the one_plot parameter being true or false
    if (args.one_plot):
        create_database_histogram(args)
    else:
        create_single_histograms(args)

if __name__ == "__main__":
    main()
