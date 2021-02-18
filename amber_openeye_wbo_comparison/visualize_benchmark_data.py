"""
Script to visualize the benchmark data for the WBO comparison
between Ambertools and OpenEye

Notes:
    This script must be ran using Python version 3.8+

Usage:
    visualize_benchmark_data.py
"""

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pickle
import plotmol

from bokeh.palettes import Spectral4
from bokeh.plotting import Figure, output_file, save, show
from plotmol.plotmol import default_tooltip_template
from scipy import stats

def read_data():
    """
    Reads the dictionary genereated by central_torsion_wbo.py
    """
    
    with open("benchmark_results/OpenFFSubstitutedPhenylSet1_benchmark.pkl", "rb") as file:
        return pickle.load(file)
        
def find_notable_differences(benchmark_data):
    """
    Finds and writes the top 25% of differences between the Ambertools WBO
    and OpenEye WBO of a molecule
    Writes to analysis/QCA_WBO_noteworthy_differences.txt
    """
    
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
            
def visualize_bargraphs(benchmark_data):
    """
    Creates bargraphs for the groups of benchmark data,
    labelling each with a group number
    """
    
    for group_num, wbo_values in enumerate(benchmark_data):
        create_bargraph(group_num+1, wbo_values)

def visualize_scatterplot(benchmark_data):
    """
    Creates a scatterplot for the benchmark data
    """
    
    all_wbo_values = {}
    for wbo_values in benchmark_data:
        all_wbo_values.update(wbo_values)

    create_scatterplot(all_wbo_values)

def create_bargraph(group_num, wbo_values):
    """
    Creates a double bar graph to visualize the WBO value comparisons
    """
    
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
    """
    Creates a scatterplot to visualize the WBO value comparisons
    """
    
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

def plot_interactive(benchmark_data):
    """
    Creates an interactive scatter plot to visualize the WBO value
    comparisons and their 2D structures
    """
    
    all_wbo_values = {}
    for wbo_values in benchmark_data:
        all_wbo_values.update(wbo_values)
        
    mols = []
    amber_wbos = []
    openeye_wbos = []
    for smiles, wbos in all_wbo_values.items():
        mols.append(smiles)
        amber_wbos.append(wbos[0])
        openeye_wbos.append(wbos[1])
    
    figure = Figure(
                    tooltips = default_tooltip_template(),
                    title = "QCA WBO Benchmark",
                    x_axis_label = "Ambertools",
                    y_axis_label = "OpenEye",
                    plot_width = 800,
                    plot_height = 800,
                    x_range = [.8,1.6],
                    y_range = [.8,1.6]
                    )
    
    plotmol.scatter(figure,
                    x = amber_wbos,
                    y = openeye_wbos,
                    smiles = mols,
                    marker = "o",
                    marker_size = 10,
                    marker_color = "black",
                    fill_color = "blue",
                    legend_label = "QCA WBO Benchmark"
                    )
                    
    figure.line(x = [.8,1.6],
                y = [.8, 1.6])
                
    slope, intercept, r_value, p_value, std_err = stats.linregress(amber_wbos, openeye_wbos)
    y = [slope * wbo + intercept for wbo in amber_wbos]
    figure.line(x = amber_wbos,
            y = y)
    
    #show(figure)
    output_file("QCA_WBO_interactiveplot/QCA_WBO_benchmark_interactive.html")
    save(figure)

if __name__ == "__main__":
    benchmark_data = read_data()
    
    #Creates a file that includes the top 25% of differences between Ambertools and OpenEye wbos
    #find_notable_differences(benchmark_data)
    
    #Create bargraphs in groups of 25 comparing the Ambertools and OpenEye wbos
    #visualize_bargraphs(benchmark_data)
    
    #Create a scatterplot comparing the Ambertools and OpenEye wbos
    #visualize_scatterplot(benchmark_data)
    
    plot_interactive(benchmark_data)
