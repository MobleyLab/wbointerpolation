#wbointerpolation
Analysis scripts for the development of wiberg bond order interpolated parameters in the OpenForceField.

`torsion_analysis.py` contains \
	-  torsion_barrier_for_molecule \
	-  checkTorsion \
	-  makeOEB \
	-  compute_r_ci \
	-  oeb2oemol \
 	-  genPlots \
	-  genData \
	-  visualize_wbo_correlation_compare \
	-  plot_tid_for_datasets


`smarts_torsions.py` contains \
	-  get_torsiondrives_matching_smarts \
	-  get_assigned_torsion_param


`Analyze_WBO_Torsion.ipynb` is a jupyter notebook making use of above functions to cache datasets locally and plot wbo versus torsion barriers


`interactive_plotting.ipynb` is a jupyter notebook using plotly and dash to drag and select data, and view the chemical structures. Inspired from [Pal Walters' blogpost](http://practicalcheminformatics.blogspot.com/2019/11/interactive-plots-with-chemical.html)

`doublering_filter.py` is a filter script to match molecules with the "[#6X3H1:1]~[#6X3:2](~[#6X3H1])-[#6X3:3](~[#6X3H1])~[#6X3H1:4]" smarts pattern. The filter currently treats the molecules as OpenEye molecules and creates .pkl files of the molecules that were found to match the smarts pattern.
