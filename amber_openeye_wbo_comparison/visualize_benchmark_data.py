"""
Script to visualize the benchmark data for the WBO comparison
between Ambertools and OpenEye

Usage:
    visualize_benchmark_data.py
"""

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os
import pickle
import plotmol

from bokeh.palettes import Category20b, Category20c
from bokeh.plotting import Figure, output_file, save, show
from plotmol.plotmol import default_tooltip_template
from scipy import stats

def read_data(dataset_file_name):
    """
    Reads the dictionary genereated by central_torsion_wbo.py
    """
    
    with open(f"benchmark_results/{dataset_file_name}.pkl", "rb") as file:
        return pickle.load(file)
        
def find_notable_differences(benchmark_data, dataset_file_name):
    """
    Finds and writes the top 25% of differences between the Ambertools WBO
    and OpenEye WBO of a molecule
    Writes to analysis/QCA_WBO_noteworthy_differences.txt
    """
        
    all_wbo_values = {}
    for wbo_values in benchmark_data:
        all_wbo_values.update(wbo_values)
        
    with open(f"noteworthy_differences/{dataset_file_name}_noteworthy_differences.txt", "w") as file:
        for smiles, data in sorted(all_wbo_values.items(),
                                   key = lambda x: abs(x[1][0][0]-x[1][0][1]))[round(.75*len(all_wbo_values)):]:
            file.write(f"Smiles: {smiles}\n")
            file.write(f"Amber wbo: {data[0][0]}\n")
            file.write(f"OpenEye wbo: {data[0][1]}\n")
            file.write(f"Difference: {abs(data[0][0]-data[0][1])}\n")
            file.write("\n")

def plot_interactive(benchmark_data, dataset_name):
    """
    Creates an interactive scatter plot to visualize the WBO value
    comparisons and their 2D structures
    """
    
    dataset_file_name = dataset_name.replace(" ", "")
    
    all_wbo_values = {}
    for wbo_values in benchmark_data:
        all_wbo_values.update(wbo_values)
    
    if len(all_wbo_values) > 0:
        mols = {}
        amber_wbos = []
        openeye_wbos = []
        for smiles, data in all_wbo_values.items():
            mols[smiles] = data[1] #Torsion indices
            amber_wbos.append(data[0][0]) #Amber wbo
            openeye_wbos.append(data[0][1]) #Openeye wbo
        
        figure = Figure(
                        tooltips = default_tooltip_template(),
                        title = f"{dataset_name} Benchmark",
                        x_axis_label = "Ambertools",
                        y_axis_label = "OpenEye",
                        plot_width = 800,
                        plot_height = 800,
                        x_range = [.8,1.6],
                        y_range = [.8,1.6]
                        )
                        
        figure.line(x = [.8,1.6],
                    y = [.8, 1.6])
        
        plotmol.scatter(figure,
                        x = amber_wbos,
                        y = openeye_wbos,
                        smiles = mols,
                        marker = "o",
                        marker_size = 10,
                        marker_color = "black",
                        fill_color = "blue",
                        legend_label = f"{dataset_name} Benchmark"
                        )
        
        #show(figure)
        output_file(f"QCA_WBO_interactiveplot/{dataset_file_name}.html")
        save(figure)
    else:
        print(f"No molecules in {dataset_name}")
    
def plot_all_interactive():
    benchmarks = {}
    colors = Category20b[20] + Category20c[20]
    
    for file in os.listdir("benchmark_results"):
        if ".pkl" in file:
            dataset_file_name = file[:-4]
            dataset_name, benchmark_data = read_data(dataset_file_name)
            dataset_file_name = dataset_name.replace(" ", "")
            
            benchmark_wbo_values = {}
            for wbo_values in benchmark_data:
                benchmark_wbo_values.update(wbo_values)
            
            benchmarks[dataset_name] = benchmark_wbo_values
            
    
    figure = Figure(
                tooltips = default_tooltip_template(),
                title = f"Amber-Openeye WBO Benchmark",
                x_axis_label = "Ambertools",
                y_axis_label = "OpenEye",
                plot_width = 1600,
                plot_height = 800,
                x_range = [.8,1.6],
                y_range = [.8,1.6]
                )
                
    figure.line(x = [.8,1.6],
            y = [.8, 1.6])
            
    color_index = 0
    for benchmark_name, benchmark_data in benchmarks.items():
        print(benchmark_name)
        mols = {}
        amber_wbos = []
        openeye_wbos = []
        for smiles, data in benchmark_data.items():
            mols[smiles] = data[1] #Torsion indices
            amber_wbos.append(data[0][0]) #Amber wbo
            openeye_wbos.append(data[0][1]) #Openeye wbo
        
        plotmol.scatter(figure,
                x = amber_wbos,
                y = openeye_wbos,
                smiles = mols,
                marker = "o",
                marker_size = 10,
                marker_color = "black",
                fill_color = colors[color_index],
                legend_label = f"{benchmark_name} Benchmark ({len(benchmark_data)} mols)"
                )
        color_index += 1
                
    figure.legend.location = "top_left"
    figure.legend.click_policy = "hide"
    figure.add_layout(figure.legend[0], "right")
   
    #show(figure)
    output_file(f"QCA_WBO_interactiveplot/All_Datasets_interactive.html")
    save(figure)

def main():
    for file in os.listdir("benchmark_results"):
        if ".pkl" in file:
            dataset_file_name = file[:-4]
            dataset_name, benchmark_data = read_data(dataset_file_name)
            try:
                plot_interactive(benchmark_data, dataset_name)
            except Exception as e:
                print(f"{dataset_name} failed with error:")
                print(f"{e}")
                print()

            #Creates a file that includes the top 25% of differences between Ambertools and OpenEye wbos
            find_notable_differences(benchmark_data, dataset_file_name)

            #Create bargraphs in groups of 25 comparing the Ambertools and OpenEye wbos
            #visualize_bargraphs(benchmark_data, dataset_name)

            #Create a scatterplot comparing the Ambertools and OpenEye wbos
            #visualize_scatterplot(benchmark_data, dataset_name)

#    plot_all_interactive()

if __name__ == "__main__":
    main()
