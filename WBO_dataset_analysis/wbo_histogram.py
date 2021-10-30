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

wbos = []

for subdir, dirs, files in os.walk("oe_wbo_calcs"):
    for file in files:
        if ".pkl" in file:
            with open(f"oe_wbo_calcs/{file}", "rb") as file:
                data = pickle.load(file)
                for values in data.values():
                    wbos += values

plt.hist(wbos)
plt.xlabel("AM1-WIBERG-ELF10 with Openeye")
plt.ylabel("# of Molecules")
plt.savefig("oe_emolecules_wbo_calcs.pdf")

