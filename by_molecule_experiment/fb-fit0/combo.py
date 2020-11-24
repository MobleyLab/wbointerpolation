import os
import json
import tempfile
from openforcefield.topology import Molecule, Topology
from openforcefield.typing.engines.smirnoff import (ForceField,
UnassignedValenceParameterException, BondHandler, AngleHandler,
ProperTorsionHandler, ImproperTorsionHandler,
vdWHandler)
from plot_td_energies import collect_td_targets_data, plot_td_targets_data
import pickle

import sys
import shutil
import numpy as np
from collections import Counter

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

#directory = '/Users/jessica/Downloads/release_1.2.0/fb-fit/targets/'
#ff = '/Users/jessica/Documents/Grad_research/WBO_Torsions_ChayaPaper/release_1.3.0_2/fb-fit/fb-fit0/forcefield/test.offxml'
ff='/Users/jessica/Documents/Grad_research/wbointerpolation/by_molecule_experiment/fb-fit0/forcefield/test.offxml'
directory='/Users/jessica/Documents/Grad_research/wbointerpolation/by_molecule_experiment/fb-fit0/targets/'

#sub_dir = [x[0] for x in os.walk(directory) if os.path.basename(x[0])[:2] == 'td']
sub_dir = [x[0] for x in os.walk(directory) if os.path.basename(x[0]).startswith('td')]

force_balance_file = 'optimize.in'


def collect_td_targets_data(tmp_folder, targets_folder):
    """ Collect td targets QM vs MM data from tmp folder.
    Returns
    -------
    data: dict
        {
            'td_SMIRNOFF_Coverage_Torsion_Set_1_000_C3H8O3': {
                'metadata': {
                    'dihedrals': [[6, 10, 12, 11]],
                    'grid_spacing': [15],
                    'dihedral_ranges': None,
                    'energy_decrease_thresh': None,
                    'energy_upper_limit': 0.05,
                    'dataset_name': 'OpenFF Group1 Torsions',
                    'entry_label': 'c1c[cH:1][c:2](cc1)[CH2:3][c:4]2ccccc2',
                    'canonical_smiles': 'c1ccc(cc1)Cc2ccccc2',
                    'torsion_grid_ids': [[-165], [-150], [-135], [-120], ...],
                    'smirks': ['[*:1]~[#6X3:2]-[#6X4:3]-[*:4]'],
                    'smirks_ids': ['t17'],
                },
                'iterdata': {
                    0: e_compare_data0,
                    1: e_compare_data1,
                }

            }
        }
        where e_compare_data is a 2D numpy array [[eqm, emm, diff, weight], ..]
    """
    #td_target_folders = [f for f in os.listdir(tmp_folder) if f.startswith('td_') and os.path.isdir(os.path.join(tmp_folder, f))]
    with open('tig.pk','rb') as f:
        match_dirs = pickle.load(f)
        match_dirs = [os.path.basename(i) for i in match_dirs]
    td_target_folders = [f for f in os.listdir(tmp_folder) if os.path.basename(f) in match_dirs]
    print('match',len(match_dirs))
    print('target',len(td_target_folders))
    #exit()
    td_target_folders.sort()
    print(f'Collecting td data from {len(td_target_folders)} target folders')
    data = {}
    for i, tgt_name in enumerate(td_target_folders):
        try:
            print(i, tgt_name)
            # load metadata from targets folder
            with open(os.path.join(targets_folder, tgt_name, 'metadata.json')) as jsonfile:
                metadata = json.load(jsonfile)
            data[tgt_name] = {'metadata': metadata, 'iterdata': {}}
            tgt_folder = os.path.join(tmp_folder, tgt_name)
            iter_folders = [f for f in os.listdir(tgt_folder) if f.startswith('iter_') and os.path.isdir(os.path.join(tgt_folder, f))]
            for iter_name in iter_folders:
                iter_idx = int(iter_name[5:])
                data_file = os.path.join(tgt_folder, iter_name, 'EnergyCompare.txt')
                data[tgt_name]['iterdata'][iter_idx] = np.loadtxt(data_file, ndmin=2)
        except:
            print('skip')
            pass
    return data

def plot_td_targets_data(data, folder_name='td_targets_plots', compare_first=True, iteration=None):
    ''' plot the target data, compare first and last iteration '''
    # aggregate metadata for counts
    smirks_counter = Counter(sid for tgt_data in data.values() for sid in tgt_data['metadata']['smirks_ids'])
    # create new folder
    if os.path.exists(folder_name):
        pass
    else:
        os.mkdir(folder_name)
    print(f"Generating plots for {len(data)} targets in {folder_name}")
    for tgt_name, tgt_data in data.items():
        print(tgt_name)
        #print('data', data.keys())
        try:
            metadata = tgt_data['metadata']
            iterdata = tgt_data['iterdata']
            iter_list = sorted(iterdata.keys())
            # create an new figure
            plt.Figure()
            last_iter = iteration if iteration != None else iter_list[-1]
            last_iter_data = iterdata[last_iter]
            # check qm - mm for last iter
            max_diff = np.max(np.abs(last_iter_data[:, 2]))
            if max_diff > 50:
                print(f"Warning, {tgt_name} iter {last_iter} max |qm-mm| > 50 kJ/mol")
            # plot qm
            plt.plot(last_iter_data[:,0], label='QM Relative Energies', color='C0')
            # plot mm
            plt.plot(last_iter_data[:,1], label=f'MM Iter {last_iter}', color='C2')
            # plot first iter
            if compare_first and len(iter_list) > 0:
                first_iter_data = iterdata[iter_list[0]]
                plt.plot(first_iter_data[:, 1], label=f'MM Iter {iter_list[0]}', color='C1')
            plt.legend()
            # use strings as x axis to support >1D torsion scans
            torsion_labels = [str(gid)[1:-1] for gid in metadata['torsion_grid_ids']]
            tick_x_locs = list(range(len(last_iter_data)))
            tick_stripe = max(int(len(last_iter_data)/10), 1) # reduce tick density
            plt.xticks(ticks=tick_x_locs[::tick_stripe], labels=torsion_labels[::tick_stripe])
            # print metadata as footnote
            footnotes = {
                "Dataset Name": metadata['dataset_name'],
                "Entry Label": metadata['entry_label'],
                "Canonical SMILES": metadata['canonical_smiles'],
                "Torsion Atom Indices": ', '.join(map(str, metadata['dihedrals'])),
                "Torsion SMIRKs": ', '.join(metadata['smirks']),
                "Torsion SMIRKs ID": ', '.join(metadata['smirks_ids']),
                "SMIRKs Total Count": ', '.join(str(smirks_counter[sid]) for sid in metadata['smirks_ids']),
            }
            footnote_lines = [f'{key:<20s} {value:>59s}' for key, value in footnotes.items()]
            plt.xlabel('Torsion Angles\n\n' + '\n'.join(footnote_lines), fontdict={'family':'monospace', 'size': 8})
            # plt.xlabel('Torsion Angles')
            plt.ylabel('Relative Energies (kcal/mol)')
            plt.tight_layout()
            filename = os.path.join(folder_name, tgt_name + '.pdf')
            plt.savefig(filename)
            plt.close()
            print('FOUND ONE')
            exit()

            # copy the figure to subfolders of the target torsions SMIRKs
            for sid in metadata['smirks_ids']:
                subfolder = os.path.join(folder_name, sid)
                if not os.path.exists(subfolder):
                    os.mkdir(subfolder)
                # move file if it's the last one
                if sid == metadata['smirks_ids'][-1]:
                    shutil.move(filename, os.path.join(subfolder, tgt_name + '.pdf'))
                else:
                    shutil.copyfile(filename, os.path.join(subfolder, tgt_name + '.pdf'))
        except:
            print('no data for', tgt_name)
            pass





def main():
    import argparse
    parser = argparse.ArgumentParser('collect td targets data from fitting tmp folder and plot them')
    parser.add_argument('-f', '--tmp_folder', default='debug.tmp')
    parser.add_argument('-t', '--targets_folder', default='targets')
    parser.add_argument('-l', '--load_pickle', help='Load data directly from pickle file')
    parser.add_argument('-i', '--iteration', type=int, default=None, help='iteration number to read rmsd')
    args = parser.parse_args()

    if args.load_pickle:
        with open(args.load_pickle, 'rb') as pfile:
            data = pickle.load(pfile)
    else:
        data = collect_td_targets_data('debug.tmp', 'targets')
        # save data as pickle file
        with open('td_plot_data.pickle', 'wb') as pfile:
            pickle.dump(data, pfile)
    plot_td_targets_data(data, iteration=args.iteration)



# Find which line of optimize.in the targets specificaiton begins at
with open(force_balance_file, 'r') as f:
    replace_index = None
    force_balance_lines = f.readlines()
    for index,line in enumerate(force_balance_lines):
        try:
            #if 'name' in line[:4]:
            if line.startswith('name'):
                replace_index = index
                break
        except:
            pass
    #print(replace_index)

count = 0
for asub in sub_dir:
    #do somehting with metadata
    metadata = os.path.join(asub,'metadata.json')
    #do something with mol2 file
    mol2 = os.path.join(asub, 'input.mol2')
    #get indices from json
    with open(metadata, 'r') as f:
        data = json.load(f)
        #print(data['attributes']['canonical_isomeric_explicit_hydrogen_mapped_smiles'])
        cmiles=data['attributes']['canonical_isomeric_explicit_hydrogen_mapped_smiles']
        dihedral=tuple(data['dihedrals'][0])
        #print(dihedral)

    #molecules =  Molecule.from_file(mol2, allow_undefined_stereo=True)
    molecules=Molecule.from_mapped_smiles(cmiles)
    topology = Topology.from_molecules([molecules])
    molecules.visualize()
    # Let's label using the Parsley force field
    forcefield = ForceField(ff, allow_cosmetic_attributes=True)

    # Run the molecule labeling
    molecule_force_list = forcefield.label_molecules(topology)
    #print(dict(molecule_force_list[0]['ProperTorsions']))
    # Print out a formatted description of the torsion parameters applied to this molecule
    plot_dict = {}
    for mol_idx, mol_forces in enumerate(molecule_force_list):
        for force_tag, force_dict in mol_forces.items():
            if force_tag == 'ProperTorsions':
                for (atom_indices, parameter) in force_dict.items():
                    #print('param id', parameter.id)
                    if (parameter.id == "TIG-fit0"):
                        #print('check1', parameter.id)
                        #check dihedral in forward or reverse order
                        if (atom_indices == dihedral) or (atom_indices == dihedral[::-1]):
                            print(atom_indices, dihedral)
                            print(asub.split('/')[-1])

                            print('check2')
                            #run forcebalance
                            tmp_file = tempfile.NamedTemporaryFile(suffix='.in')
                            if 0:
                                with open(tmp_file.name,'w+t') as f1:
                                    for index,line in enumerate(force_balance_lines):
                                        if index != replace_index:
                                            print(index, line)
                                            f1.write(line)
                                        else:
                                            print('original', line)
                                            print('replace ', line.split()[0] + ' ' + os.path.basename(asub))
                                            f1.write(line.split()[0] + ' ' + os.path.basename(asub)+"\n")
                                    f1.flush()
                                    f1.seek(0)
                                    print('ForceBalance ' + f1.name)
                                    print(os.path.isfile(f1.name))
                                    print(os.system('ForceBalance ' +f1.name))
                                    print('done')

                                    import time
                            if 1:
                                with open('debug.in','w') as f1:
                                    for index,line in enumerate(force_balance_lines):
                                        if index != replace_index:
                                            f1.write(line)
                                        else:
                                            f1.write(line.split()[0] + ' ' + os.path.basename(asub)+"\n")
                                print('ForceBalance ' + 'debug.in')
                                print(os.system('ForceBalance ' + 'debug.in'))
                                print('done')
    #makePickle()
    data = collect_td_targets_data('debug.tmp', 'targets')
    # save data as pickle file
    with open('td_plot_data.pickle', 'wb') as pfile:
        pickle.dump(data, pfile)
    plot_td_targets_data(data)



