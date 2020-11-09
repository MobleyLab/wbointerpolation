#!/usr/bin/env python

"""
Utility functions for SMARTS-based querying of TorsionDrive datasets on QCArchive
Created jointly by:
- @trevorgokey
- @pavankum
- @dotsdl
"""

from openforcefield.topology import Molecule, Topology
from openforcefield.typing.engines.smirnoff import ForceField

# NOTES: 
# 1. Perhaps a default for datasets grabbing *all* OpenFF TorsionDrives?


def get_torsiondrives_matching_smarts(smarts, datasets, client):
    """Get all TorsionDrive entries from the given QCArchive instance datasets that
    match the given SMARTS pattern.

    Parameters
    ----------
    smarts : str
        SMARTS to query TorsionDrives against;
        if the SMARTS has 4 indexed atoms, then these must match the dihedral atoms exercised in the TorsionDrive;
        if the SMARTS has 2 indexed atoms, then these must match the atoms of the central, rotated bond in the TorsionDrive.
    datasets : iterable of strings
        TorsionDriveDataset names to sample from.
    client : qcportal.FractalClient
        Fractal client to use for database queries.

    Returns
    -------
    tdentries : list
        List of TDEntries matching given SMARTS.

    Examples
    --------
    Get back TDEntries corresponding to a C-C bond torsiondrive:
    
    >>> from qcportal import FractalClient
    >>> client = FractalClient()
    >>> smarts = "[*:1]~[#6:2]-[#6:3]~[*:4]"
    >>> datasets = ["OpenFF Substituted Phenyl Set 1"]
    >>> tdentries = get_torsiondrives_matching_smarts(smarts, dataset, client)

    Equivalent to the above due to wildcard matching, but specifying only central bond:

    >>> smarts = "[#6:1]-[#6:2]"
    >>> tdentries = get_torsiondrives_matching_smarts(smarts, dataset, client)


    """


    # first, we want to grab all datasets into memory
    tdrs = []
    for dataset in datasets:
        ds = client.get_collection("TorsionDriveDataset", dataset)

        for entry in ds.data.records.values():

            # need to build molecule from smiles so we can query against it
            mol_smiles = entry.attributes["canonical_isomeric_explicit_hydrogen_mapped_smiles"]
            offmol = Molecule.from_mapped_smiles(mol_smiles)

            # apply SMARTS to each TorsionDrive object
            matching_indices = offmol.chemical_environment_matches(smarts)

            # if we get back nothing, move on
            if len(matching_indices) == 0:
                continue

            # by convention, we only have one driven torsion
            # would need to revisit if we are working with 2D torsions
            dihedral_indices = entry.td_keywords.dihedrals[0]

            # assemble matches
            for indices in matching_indices:
                if len(indices) == 2:
                    if sorted(indices) == sorted(dihedral_indices[1:3]):
                        tdrs.append(entry)
                        break
                elif len(indices) == 4:
                    if sorted(indices) == sorted(dihedral_indices):
                        tdrs.append(entry)
                        break
                else:
                    raise ValueError("Number of indices returned not 2 or 4;" 
                                     " check number of tagged atoms in SMARTS")

    return tdrs


def get_assigned_torsion_param(tdentry, forcefield):
    """Get the OpenFF forcefield torsion parameter ultimately assigned to the
    given TorsionDrive entry's torsion dihedral.

    Parameters
    ----------
    tdentry : TDEntry
        TDEntry (TorsionDrive entry) to operate on;
        will be used to generate molecule, extract dihedral indices driven.
    forcefield : str, ForceField
        OpenFF forcefield to apply.

    Returns
    -------
    torsion_params : ProperTorsion
        Dict-like object with attributes giving the applied torsion parameters

    Examples
    --------
    Starting with TDEntries from usage of `get_torsiondrives_matching_smarts`
    (see its Example), we can get back the parameter assigned to this by, say
    `"openff-1.0.0.offxml"`:
    
    >>> from openforcefield.typing.engines.smirnoff import ForceField
    >>> tdentries = get_torsiondrives_matching_smarts(smarts, dataset, client)
    >>> ff = ForceField('openff-1.0.0.offxml')
    >>> assigned = [smarts_torsions.get_assigned_torsion_param(tdentry, ff)
                    for tdentry in tdentries]

    >>> print([t.id for t in assigned])
        ['t47', 't47', 't47', 't47', ...]

    """
    mol_smiles = tdentry.attributes["canonical_isomeric_explicit_hydrogen_mapped_smiles"]
    offmol = Molecule.from_mapped_smiles(mol_smiles)

    if isinstance(forcefield, str):
        forcefield = ForceField(forcefield)

    # apply forcefield parameters
    topology = Topology.from_molecules(offmol)
    
    # we only have one molecule by definition here, so extracting 0th
    molecule_forces = forcefield.label_molecules(topology)[0]

    # by convention, we only have one driven torsion
    # would need to revisit if we are working with 2D torsions
    dihedral_indices = tdentry.td_keywords.dihedrals[0]

    # get torsion parameters corresponding to dihedral indices
    torsions = molecule_forces["ProperTorsions"]
    torsion_params = torsions.get(dihedral_indices)

    # if None, try reversing it
    if torsion_params is None:
        torsion_params = torsions[dihedral_indices[::-1]]

    return torsion_params