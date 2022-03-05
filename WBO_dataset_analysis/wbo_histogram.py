"""
Script for creating a histogram based on the WBO values calculated in wbo_calcs.py. The script
iterates through the openff_results folder, creating a large dataset of all the wbo values
from every group and plotting them all on a histogram. The result is saved as
emolecules_wbo_calcs.pdf.

Usage:
    wbo_histogram.py
"""


import matplotlib.pyplot as plt
import os
import pickle

all_wbos = []

with open("wbo_results/doublering_wbos.pkl", "rb") as file:
    data = pickle.load(file)
    for smiles, wbos in data.items():
        all_wbos += wbos

all_wbos[:] = [wbo for wbo in all_wbos if wbo != 0.0]

print(max(all_wbos))
print(min(all_wbos))

print(f"Total: {len(all_wbos)}")

plt.hist(all_wbos, log=True)
plt.xlabel("AM1-WIBERG-ELF10 with Openeye")
plt.ylabel("log(# of molecules)")
plt.title("[#6X3H1:1]~[#6X3:2](~[#6X3H1])-[#6X3:3](~[#6X3H1])~[#6X3H1:4]\nmol WBOs")
plt.savefig("wbo_visualizations/doublering_wbos_graph_log.pdf")
