from plot_td_energies import plot_td_targets_data
import pickle

with open('plot_data.pk', 'rb') as f:
	data = pickle.load(f)
	plot_td_targets_data(data)