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
directory = '/Users/jessica/Downloads/release_1.2.0/fb-fit/targets/'
ff = '/Users/jessica/Documents/Grad_research/WBO_Torsions_ChayaPaper/release_1.3.0_2/fb-fit/fb-fit0/forcefield/test.offxml'
#sub_dir = [x[0] for x in os.walk(directory) if os.path.basename(x[0])[:2] == 'td']
sub_dir = [x[0] for x in os.walk(directory) if os.path.basename(x[0]).startswith('td')]

force_balance_file = 'optimize.in'


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
output = []
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
                            output.append(asub)

print(output)
with open('tig.pk', 'wb') as f:
    pickle.dump(output, f)



