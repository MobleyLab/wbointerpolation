# Wiberg Bond Order Analysis on "[#6X3H1:1]~[#6X3:2](~[#6X3H1])-[#6X3:3](~[#6X3H1])~[#6X3H1:4]" Molecules

## Environment details

Exact details about the environment can be found in env_details.txt

## Data

### oe_results

Contains a .zip file of the combined results of every group from doublering_filter.py in a .oeb file

### wbo_results

Contains the results of doublering scripts in .pkl that calculate the wbos for protomers, tautomers, and regular versions of all emolecules from in oe_results/doublering_mols.oeb

### wbo_visualizatons 

Contains visualizations of the wbo values in wbo_results

## Scripts

### doublering_filter.py

Script to filter a smiles database for all molecules that contain the "[#6X3H1:1]~[#6X3:2](~[#6X3H1])-[#6X3:3](~[#6X3H1])~[#6X3H1:4]" smarts string. The script takes the file and a group number for the molecules as parameters so that large databases can be divided into multiple parts and multiple groups can be ran at once. For other scripts to work, all files should be compiled into one file called oe_results/doublering_mols.oeb.

### doublering_protomer_wbo_calcs.py

Script to calculate the Wiberg Bond Order values for all protomers of the filtered molecules in oe_results/doublering_mols.oeb. Results are stored in wbo_results/doublering_protomer_wbos.pkl.

### doublering_protomer_wbo_calcs.slurm

Slurm script to run the doublering_protomer_wbo_calcs.py script on a green planet machine

### doublering_tautomer_wbo_calcs.py

Script to calculate the Wiberg Bond Order values for all tautomers of the filtered molecules in oe_results/doublering_mols.oeb. Results are stored in wbo_results/doublering_tautomer_wbos.pkl.

### doublering_tautomer_wbo_calcs.slurm

Slurm script to run the doublering_tautomer_wbo_calcs.py script on a green plane
t machine

### doublering_wbo_calcs.py

Script to calculate the Wiberg Bond Order values for all filtered molecules in oe_results/doublering_mols.oeb. Results are stored in wbo_results/doublering_wbos.pkl.

### doublering_wbo_calcs.slurm

Slurm script to run the doublering_wbo_calcs.py script on a green plane
t machine

### wbo_histogram.py

Script for creating a histogram based on a .pkl from wbo_results. The script iterates through the file, creating a large dataset of all the wbo values from every different molecule and plots them on a graph.

