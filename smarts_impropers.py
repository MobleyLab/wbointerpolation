#!/usr/bin/env python

"""
Utility functions for SMARTS-based querying of GeoOpt datasets on QCArchive
Created by:
- @pavankum
"""

from openff.toolkit.topology import Molecule, Topology
from openff.toolkit.typing.engines.smirnoff import ForceField
from rdkit.Chem import rdMolTransforms
from simtk.unit import bohrs

def get_gopt_matching_improper(smarts, datasets, client, hess=None):
    """Get all Optimization entries from the given QCArchive instance datasets that
    match the given SMARTS pattern.

    Parameters
    ----------
    smarts : str
        SMARTS to query optimization datasets against optimzation datasets
    datasets : iterable of strings
        OptimizationDataset names to sample from.
    client : qcportal.FractalClient
        Fractal client to use for database queries.

    Returns
    -------
    opt_recs:
        List of optimization records matching the given smarts
    ind_set
        Dict of a list of matching indices per record id
    wbos
        Dict of wbos corresponding to the improper 2nd and 4th atoms

    Examples
    --------
    Get back TDEntries corresponding to a C-C bond torsiondrive:

    >>> from qcportal import FractalClient
    >>> client = FractalClient()
    >>> smarts = "[*:1]~[#6:2]-[#6:3]-[*:4]"
    >>> datasets = ["OpenFF Gen 2 Opt Set 1 Roche"]
    >>> opts, ind_set, wbos = get_gopt_matching_dihedral(smarts, dataset, client)

    """
    from collections import defaultdict

    # first, we want to grab all datasets into memory
    opt_recs = []
    ind_set = defaultdict(list)
    dihs_set = []
    wbos = defaultdict(list)
    failed_smiles = []
    for dataset in datasets:
        if hess:
            ds = client.get_collection("Dataset", dataset)
        else:
            ds = client.get_collection("OptimizationDataset", dataset)

        for entry in ds.data.records.values():
            # need to build molecule from smiles so we can query against it
            mol_smiles = entry.attributes[
                "canonical_isomeric_explicit_hydrogen_mapped_smiles"
            ]
            try:
                offmol = Molecule.from_mapped_smiles(mol_smiles, allow_undefined_stereo=True)
            except:
                failed_smiles.append(mol_smiles)
                print("Failed: ", mol_smiles)
                continue

            # apply SMARTS to each Optimization object
            matching_indices = offmol.chemical_environment_matches(smarts)

            # if we get back nothing, move on
            if len(matching_indices) == 0:
                continue
            else:
                optrec = ds.get_record(name=entry.name, specification="default")
                if optrec.status == "ERROR":
                    continue
                else:
                    offmol.add_conformer(optrec.get_final_molecule().geometry * bohrs)
                    try:
                        offmol.assign_fractional_bond_orders(bond_order_model="am1-wiberg-elf10")
                    except:
                        failed_smiles.append(mol_smiles)
                        print("Failed conformer generation for ELF10 charges: ", mol_smiles)
                        continue
                    ind_passed = []
                    for impropers in matching_indices:
                        if [impropers[1], impropers[3]] in ind_passed:
                            continue
                        elif [impropers[3], impropers[1]] in ind_passed:
                            continue
                        else:
                            bond = offmol.get_bond_between(impropers[1], impropers[3])
                            wbos[optrec.id].append(bond.fractional_bond_order)
                            ind_passed.append([impropers[1], impropers[3]])
                            ind_set[optrec.id].append(impropers)

                opt_recs.append(entry)
                #ind_set.append(matching_indices)
    print(failed_smiles)

    return opt_recs, ind_set, wbos

def show_oemol_struc(oemol, torsions=False, atom_indices=[], width=500, height=300):
    from IPython.display import Image
    from openeye import oechem, oedepict

    # Highlight element of interest
    class NoAtom(oechem.OEUnaryAtomPred):
        def __call__(self, atom):
            return False

    class AtomInTorsion(oechem.OEUnaryAtomPred):
        def __call__(self, atom):
            return atom.GetIdx() in atom_indices

    class NoBond(oechem.OEUnaryBondPred):
        def __call__(self, bond):
            return False

    class BondInTorsion(oechem.OEUnaryBondPred):
        def __call__(self, bond):
            return (bond.GetBgn().GetIdx() in atom_indices) and (
                bond.GetEnd().GetIdx() in atom_indices
            )

    class CentralBondInTorsion(oechem.OEUnaryBondPred):
        def __call__(self, bond):
            return (bond.GetBgn().GetIdx() in atom_indices[1:3]) and (
                bond.GetEnd().GetIdx() in atom_indices[1:3]
            )

    opts = oedepict.OE2DMolDisplayOptions(width, height, oedepict.OEScale_AutoScale)
    opts.SetAtomPropertyFunctor(oedepict.OEDisplayAtomIdx())

    oedepict.OEPrepareDepiction(oemol)
    img = oedepict.OEImage(width, height)
    display = oedepict.OE2DMolDisplay(oemol, opts)
    if torsions:
        atoms = oemol.GetAtoms(AtomInTorsion())
        bonds = oemol.GetBonds(NoBond())
        abset = oechem.OEAtomBondSet(atoms, bonds)
        oedepict.OEAddHighlighting(
            display,
            oechem.OEColor(oechem.OEYellow),
            oedepict.OEHighlightStyle_BallAndStick,
            abset,
        )

    oedepict.OERenderMolecule(img, display)
    png = oedepict.OEWriteImageToString("png", img)
    return Image(png)


