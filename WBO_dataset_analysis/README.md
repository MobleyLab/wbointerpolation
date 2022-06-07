# Wiberg Bond Order Analysis on \[\#6X3H1:1]\~\[\#6X3:2\]\(\~\[\#6X3H1\]\)\-\[\#6X3:3\]\(\~\[\#6X3H1\]\)\~\[\#6X3H1:4\] Molecules

## Generic info

### Environment details

Exact details about the environment can be found in env_details.txt.

## eMolecules experiment

This analysis was done using eMolecules, which can be found in the version.smi.gz file at https://downloads.emolecules.com/free/2022-03-01/.

### Data

#### oe_results

Contains a .zip file of the combined results of every group from emolecules_filter.py in a .oeb file

#### wbo_results

Contains the results of emolecules scripts in .pkl that calculate the wbos for protomers, tautomers, and regular versions of all emolecules from in oe_results/emolecules_mols.oeb

#### wbo_visualizatons 

Contains visualizations of the wbo values in wbo_results

### eMolecules Scripts

#### emolecules_filter.py

Script to filter a smiles database for all molecules that contain the \[\#6X3H1:1]\~\[\#6X3:2\]\(\~\[\#6X3H1\]\)\-\[\#6X3:3\]\(\~\[\#6X3H1\]\)\~\[\#6X3H1:4\] smarts string. The script takes the file and a group number for the molecules as parameters so that large databases can be divided into multiple parts and multiple groups can be ran at once. For other scripts to work, all files should be compiled into one file called oe_results/emolecules_mols.oeb, or the provided file can be used.

#### emolecules_protomer_wbo_calcs.py

Script to calculate the Wiberg Bond Order values for all protomers of the filtered molecules in oe_results/emolecules_mols.oeb. Results are stored in wbo_results/emolecules_protomer_wbos.pkl.

#### emolecules_tautomer_wbo_calcs.py

Script to calculate the Wiberg Bond Order values for all tautomers of the filtered molecules in oe_results/emolecules_mols.oeb. Results are stored in wbo_results/emolecules_tautomer_wbos.pkl.

#### emolecules_wbo_calcs.py

Script to calculate the Wiberg Bond Order values for all filtered molecules in oe_results/emolecules_mols.oeb. Results are stored in wbo_results/emolecules_wbos.pkl.

#### eMolecules Visualization Script

#### wbo_histogram.py

Script for creating a histogram based on a .pkl from wbo_results directory. The script iterates through an input file, creating a large dataset of all the wbo values from every different molecule and plots them on a graph.

## PubChem experiment

This analysis was done using PubChem, the entirety of the database can be found at https://ftp.ncbi.nlm.nih.gov/pubchem/Compound/CURRENT-Full/SDF/, or with the specific search on \[\#6X3H1:1]\~\[\#6X3:2\]\(\~\[\#6X3H1\]\)\-\[\#6X3:3\]\(\~\[\#6X3H1\]\)\~\[\#6X3H1:4\] molecules at https://pubchem.ncbi.nlm.nih.gov/#query=%5B%236X3H1%3A1%5D~%5B%236X3%3A2%5D(~%5B%236X3H1%5D)-%5B%236X3%3A3%5D(~%5B%236X3H1%5D)~%5B%236X3H1%3A4%5D.

### Heavy Atoms Data

The heavy atoms data is taken from the specific search on \[\#6X3H1:1]\~\[\#6X3:2\]\(\~\[\#6X3H1\]\)\-\[\#6X3:3\]\(\~\[\#6X3H1\]\)\~\[\#6X3H1:4\] in the pubchem database. Unfortunately, these searches do not search the entire database so only a small section is taken with each search. Our solution to get more of the database using this means was to toggle the "Search All" box under the "Substructure Search" heading and filter by heavy atom count, going in groups by 10 (ex: 10-20, 30-40, 50-60, ...). 

#### heavy_atoms_datasets

Contains the raw datasets from the PubChem website for each filtered result as .json files, which can be found at https://pubchem.ncbi.nlm.nih.gov/#query=%5B%236X3H1%3A1%5D~%5B%236X3%3A2%5D(~%5B%236X3H1%5D)-%5B%236X3%3A3%5D(~%5B%236X3H1%5D)~%5B%236X3H1%3A4%5D. Based on the names you can see the amount of heavy atoms filtered for in the groups of 10. There were an odd number of groups, so the final group is atoms with a heavy count in the range 150-159.

#### wbo_results_heavy_atoms

Contains .pkl files that contain the calculated wbo values for all of the molecules in the heavy_atoms_datasets directory. There is a corresponding file in this directory for every file in heavy_atoms_datasets.

#### wbo_visualizations_heavy_atoms

Contains the visualizations for the wbo values in wbo_results_heavy_atoms in a non log graph.

#### wbo_visualizations_heavy_atoms_log

Contains the visualizations for the wbo values in wbo_results_heavy_atoms in a log graph.

### Heavy Atoms Scripts

#### pubchem_wbo_calcs_json.py

Takes a .json file from heavy_atoms_datasets as input and create the .pkl file of wbos that ends up in wbo_results_heavy_atoms.

### PubChem Full Database

#### pubchem_full_dataset

Contains the raw datasets as .sdf files from the PubChem website for the entire database at https://ftp.ncbi.nlm.nih.gov/pubchem/Compound/CURRENT-Full/SDF/.

#### pubchem_filtered

Contains filtered datasets from pubchem_full_dataset which are stored as .oeb files. Files in this directory only include molecules that match the \[\#6X3H1:1]\~\[\#6X3:2\]\(\~\[\#6X3H1\]\)\-\[\#6X3:3\]\(\~\[\#6X3H1\]\)\~\[\#6X3H1:4\] string.

#### wbo_results_pubchem_full

Contains .pkl files that contain the calculated wbo values for all of the molecules in the pubchem_filtered directory. There is a corresponding file in this directory for every file in pubchem_filtered.

#### wbo_visualizations_pubchem_full

Contains the visualizations for the wbo values in wbo_results_pubchem_full in a non log graph.

#### wbo_visualizations_pubchem_full_log

Contains the visualizations for the wbo values in wbo_results_pubchem_full in a log graph.

### PubChem Full Database Scripts

#### pubchem_filter.py

Script to iterate through the emolcules database looking for molecules that match the \[\#6X3H1:1]\~\[\#6X3:2\]\(\~\[\#6X3H1\]\)\-\[\#6X3:3\]\(\~\[\#6X3H1\]\)\~\[\#6X3H1:4\] smarts pattern.

#### pubchem_wbo_calcs.py

Script to calculate the Wiberg Bond Order values for all filtered molecules in pubchem_filtered

### PubChem Visualization Script

#### pubchem_wbo_histograms.py

Script for creating a histogram based on the .pkl files from wbo_results in the pubchem directory. 

The script has two options: 
    1. Iterate through the entire directory, creating an individual graph for the data in each .pkl file
    2. Iterate through the entire directory, creating a collective graph based on the data in all of the .pkl files combined

The option is selected using the following command line parameters. The --one_plot parameter is default False which does option 1 for the script, 
users can toggle one_plot to True which does option 2 for the script. The --log parameter uses the log of the results when creating the histogram, 
this parameter is also default value False.

## Greenplanet details

#### emolecules_protomer_wbo_calcs.slurm

Slurm script to run the emolecules_protomer_wbo_calcs.py script on a green planet machine

#### emolecules_tautomer_wbo_calcs.slurm

Slurm script to run the emolecules_tautomer_wbo_calcs.py script on a green planet machine

#### emolecules_wbo_calcs.slurm

Slurm script to run the emolecules_wbo_calcs.py script on a green planet machine