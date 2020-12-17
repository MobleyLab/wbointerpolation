# by jessica maat
# minimal example of checking if bond parameter matches to cmiles with custom ff

from openforcefield.topology import Molecule, Topology
from openforcefield.typing.engines.smirnoff import (ForceField,
UnassignedValenceParameterException, BondHandler, AngleHandler,
ProperTorsionHandler, ImproperTorsionHandler,
vdWHandler)



def checkParam(cmiles, ff2):

    molecules=Molecule.from_mapped_smiles(cmiles)
    topology = Topology.from_molecules([molecules])


    #added
    # Let's label using the Parsley force field
    forcefield2 = ForceField(ff2, allow_cosmetic_attributes=True)
    # Run the molecule labeling
    molecule_force_list = forcefield2.label_molecules(topology)
    #print(dict(molecule_force_list[0]['ProperTorsions']))
    # Print out a formatted description of the torsion parameters applied to this molecule
    #plot_dict = {}
    for mol_idx, mol_forces in enumerate(molecule_force_list):
        for force_tag, force_dict in mol_forces.items():
            print(force_tag)
            if force_tag == 'Bonds':
                for (atom_indices, parameter) in force_dict.items():
                    if parameter.id == 'b1':
                        print('match')
                        return cmiles


cmiles='[H:29][c:1]1[c:2]([c:5]([c:6]([c:3]([c:4]1[O:26][C:17]([H:50])([H:51])[H:52])[H:31])[S:28](=[O:24])(=[O:25])[N:23]2[C:13]([C:9]([C:16]([C:10]([C:14]2([H:46])[H:47])([H:38])[H:39])([H:49])[C:20]([H:58])([H:59])[C:21]([H:60])([H:61])[C:19]([H:56])([H:57])[C:15]3([C:7]([C:11]([N:22]([C:12]([C:8]3([H:34])[H:35])([H:42])[H:43])[H:62])([H:40])[H:41])([H:32])[H:33])[H:48])([H:36])[H:37])([H:44])[H:45])[O:27][C:18]([H:53])([H:54])[H:55])[H:30]'
checkParam(cmiles, 'openff-1.3.0.offxml')
