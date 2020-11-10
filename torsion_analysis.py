#Created by: 
#    @jmaat (Jessica Maat)
#    @pavankum (Pavan Behara)


# imports
import os
import sys
import time

import cmiles
import fragmenter
import matplotlib.pyplot as plt
import numpy as np
import qcportal as ptl
from fragmenter import chemi
from fragmenter.utils import HARTREE_2_KJMOL
from openeye import oechem
from openforcefield.topology import Molecule, Topology
from openforcefield.typing.engines.smirnoff import ForceField
from scipy import stats
from simtk import unit


# functions
def torsion_barrier_for_molecule(tdr_object, mapped_smiles, show_plots=False):
    """
    Takes in a single torsion drive record that has energies from multiple conformers (at different torsion angles),
    evaluates the torsion barrier

    Parameters
    ----------
    tdr_object : object
        torsion drive record from QC archive for a molecule

    Returns
    -------
    mol: oemol object
        oemol from the smiles in dataframe index that contains datatags with the following:
        tdr_object.id : int (id of the TD record) with datatag "TDid"
        dihedral_indices: list (list of atom indices for which torsion is driven in this record) datatag "TDindices"
        torsion_barrier: float (torsion barrier energy in KJ/mol, maximum of all the barriers) datatag "TB"
        cmiles: str (string for the cmiles of the molecule in canonical_isomeric_explicit_hydrogen_mapped_smiles)
        datatag "cmiles"
    """
    energies = list(tdr_object.get_final_energies().values())
    tmp = list(tdr_object.get_final_energies().keys())
    angles = [i[0] * np.pi / 180 for i in tmp]
    angles, energies = zip(*sorted(zip(angles, energies)))
    angles = np.array(angles)
    energies = np.array(energies)
    angles = np.append(
        angles[-3:] - 2 * np.pi, np.append(angles, angles[:3] + 2 * np.pi)
    )
    energies = np.append(energies[-3:], np.append(energies, energies[:3]))

    idx = []
    for i in range(len(angles) - 2):
        m1 = (energies[i + 1] - energies[i]) / (angles[i + 1] - angles[i])
        m2 = (energies[i + 2] - energies[i + 1]) / (angles[i + 2] - angles[i + 1])
        if np.sign(m1) == np.sign(m2):
            continue
        else:
            idx.append(i + 1)

    if show_plots:
        min_ener = min(energies)
        energies_y = (energies - min_ener) * HARTREE_2_KJMOL
        fontsize = 14
        plt.figure()
        plt.plot(
            angles * 180 / np.pi,
            energies_y,
            "b-X",
            angles[idx] * 180 / np.pi,
            energies_y[idx],
            "ro",
        )
        plt.legend(["QM data", "Max, min"], bbox_to_anchor=(1, 1), fontsize=fontsize)
        plt.title("Torsion drive interpolation", fontsize=fontsize)
        plt.xlabel("Dihedral Angle [Degrees]", fontsize=fontsize)
        plt.ylabel("Relative energy [KJ / mol]", fontsize=fontsize)
        plt.xticks(fontsize=fontsize)
        plt.yticks(fontsize=fontsize)
        fig_name = "plot_" + tdr_object.id + ".png"
        plt.savefig(fig_name)
        plt.show()

    torsion_barriers = []
    for i in range(int(len(idx) - 1)):
        torsion_barriers.append(
            abs(HARTREE_2_KJMOL * abs(energies[idx[i]] - energies[idx[i + 1]]))
        )
    torsion_barriers = np.array(torsion_barriers)

    # get dihedral indices and pass on to get_wbo function
    dihedral_indices = tdr_object.dict()["keywords"]["dihedrals"][0]
    offmol = Molecule.from_mapped_smiles(mapped_smiles)
    offmol.assign_fractional_bond_orders()
    bond = offmol.get_bond_between(dihedral_indices[1], dihedral_indices[2])
    mol = chemi.smiles_to_oemol(mapped_smiles)
    mol.SetData("WBO", bond.fractional_bond_order)
    mol.SetData("TB", max(torsion_barriers))
    mol.SetData("TDindices", dihedral_indices)
    mol.SetData("TDid", tdr_object.id)
    mol.SetData("cmiles", mapped_smiles)

    return mol


def loadDataset_low(datasetName, specification, benchmark_smiles, qca_overlapped_entries):
    """
    Low level call to load each torsion drive dataset and return a list of molecules

        Parameters
        ----------
        datasetName : str
            torsion drive dataset name.
        specification : str
            specification in the dataset. Example: "B3LYP-D3", "default", "UFF"

        Returns
        -------
        molList : list of objects
            each row contains the tdr_object.id, dihedral_indices, torsion_barrier, oemol_object
    """
    while True:
        try:
            assert datasetName
            break
        except AssertionError:
            print("datasetName is empty. Check input list of dataset tuples")
            raise
    while True:
        try:
            assert specification
            break
        except AssertionError:
            print("specification is empty. Check input list of dataset tuples")
            raise

    # initiate qc portal instance
    client = ptl.FractalClient()
    # from the TorsionDriveDataset collection picking up given datasetName
    ds = client.get_collection("TorsionDriveDataset", datasetName)
    ds.status([specification], status="COMPLETE")

    # Serial implementation

    # Hardcoding benchmark molecules from the lim_mobley_parsely_benchmark
    # https://openforcefield.org/force-fields/force-fields/
    # https://github.com/MobleyLab/benchmarkff/blob/91476147f35579bc52bf984839fd20c72a61d76d/molecules/set_v03_non_redundant/trim3_full_qcarchive.smi
    
    with open(benchmark_smiles) as f:
        bm_smiles = f.readlines()
    bm_mols = [Molecule.from_smiles(smiles) for smiles in bm_smiles]
    
    tb = []
    overlaps = 0
    qca_entries = []
    
    for i in range(ds.df.size):
        if ds.df.iloc[i, 0].status == "COMPLETE":
            smiles = ds.df.index[i]
            mapped_smiles = ds.get_entry(smiles).attributes[
                "canonical_isomeric_explicit_hydrogen_mapped_smiles"
            ]
            mol1 = Molecule.from_mapped_smiles(mapped_smiles)
            not_identical = True
            for mol in bm_mols:
                isomorphic,atom_map = Molecule.are_isomorphic(mol1, 
                                                  mol,
                                                  return_atom_map=False,
                                                  aromatic_matching=False,
                                                  formal_charge_matching=False,
                                                  bond_order_matching=False,
                                                  atom_stereochemistry_matching=False,
                                                  bond_stereochemistry_matching=False,
                                                          )
                if(isomorphic):
                    not_identical = False
                    overlaps += 1
                    entry = ds.get_entry(smiles)
                    qca_entries.extend(list(entry.initial_molecules))
                    break
            if(not_identical): 
                tb.append(torsion_barrier_for_molecule(ds.df.iloc[i, 0], mapped_smiles))
    
    # overlaps_qca_ids.txt is also a hardcoded file
    with open(qca_overlapped_entries, "a") as f:
        for item in qca_entries:
            f.write("%s\n" % item)
        
    print("No. of overlaps with benchmark set, qca entries added to overlaps_qca_ids.txt: ", overlaps)
    print("No. of COMPLETE and not overlapping with benchmark in this dataset:", len(tb), "out of ", len(ds.df))
    return tb


def checkTorsion(molList, ff_name):
    """
    Take mollist and check if the molecules in a list match a specific torsion id

        Parameters
        ----------
        molList : List of objects
            List of oemols with datatags generated in genData function

        Returns
        -------
        molList : list of objects
            List of oemol objects that have a datatag "IDMatch" that contain the torsion id
            involved in the QCA torsion drive
    """

    matches = []
    count = 0
    mols = []
    for mol in molList:
        molecule = Molecule.from_mapped_smiles(mol.GetData("cmiles"))
        topology = Topology.from_molecules(molecule)
        # Let's label using the Parsley force field
        forcefield = ForceField(ff_name)
        # Run the molecule labeling
        molecule_force_list = forcefield.label_molecules(topology)
        params = []
        # Print out a formatted description of the torsion parameters applied to this molecule
        for mol_idx, mol_forces in enumerate(molecule_force_list):
            # print(f'Forces for molecule {mol_idx}')
            for force_tag, force_dict in mol_forces.items():
                if force_tag == "ProperTorsions":
                    for (atom_indices, parameter) in force_dict.items():
                        params.append(parameter.id)
                        if atom_indices == mol.GetData("TDindices") or tuple(
                            reversed(atom_indices)
                        ) == mol.GetData("TDindices"):
                            count += 1
                            mol.SetData("IDMatch", parameter.id)
                            mols.append(mol)
    print(
        "Out of "
        + str(len(molList))
        + " molecules, "
        + str(count)
        + " were processed with checkTorsion()"
    )

    return mols


def makeOEB(oemolList, tag):
    """
    Take mollist and create oeb file using the tag as the .oeb file name

        Parameters
        ----------
        molList : List of objects
            List of oemols with datatags generated in genData function
        tag : String
            Title of the oeb file

        Returns
        -------
    """
    ofile = oechem.oemolostream(tag + ".oeb")
    for mol in oemolList:
        oechem.OEWriteConstMolecule(ofile, mol)
    ofile.close()
    return


def compute_r_ci(wbos, max_energies):
    return (stats.linregress(wbos, max_energies)[2]) ** 2


def oeb2oemol(oebfile):
    """
    Takes in oebfile and generates oemolList
        Parameters
        ----------
        oebfile : String
            Title of an oeb file
        Returns
        -------
        mollist : List of objects
            List of OEMols in the .oeb file

    """
    ifs = oechem.oemolistream(oebfile)
    mollist = []

    for mol in ifs.GetOEGraphMols():
        mollist.append(oechem.OEGraphMol(mol))

    return mollist


# plotting functions
def genPlots(fileName, fname):
    """
    Generates a .pdf plot from a .oeb file of wbo versus torsion barrier height
        Parameters
        ----------
        fileName : String
            .oeb file name for the molecules that contain datatags with plotting information
        fname : String
            The output .pdf file name for the resulting files from the plotting function

        Returns
        -------
    """

    molList = oeb2oemol(fileName)

    torsionDict = {}

    for m in molList:
        tid = m.GetData("IDMatch")
        torsionDict[tid] = {}
        torsionDict[tid]["tb"] = []
        torsionDict[tid]["wbo"] = []

    for m in molList:
        tid = m.GetData("IDMatch")
        if np.isfinite(m.GetData("TB")) and np.isfinite(m.GetData("WBO")):
            torsionDict[tid]["tb"].append(m.GetData("TB"))
            torsionDict[tid]["wbo"].append(m.GetData("WBO"))

    colors = fragmenter.chemi._KELLYS_COLORS
    colors = colors * 8

    fig, ax = plt.subplots()

    markers = ["o", "v", "<", "P", ">", "^", "d", "X"]
    for i, (key, item) in enumerate(torsionDict.items()):
        ax.scatter(
            torsionDict[key]["wbo"],
            torsionDict[key]["tb"],
            c=colors[i],
            marker=markers[int(i / 20)],
            label=key,
        )
        if len(torsionDict[key]["wbo"]) < 2:
            continue
        else:
            slope, intercept, r_value, p_value, std_err = stats.linregress(
                torsionDict[key]["wbo"], torsionDict[key]["tb"]
            )
            plt.plot(
                np.unique(torsionDict[key]["wbo"]),
                np.poly1d([slope, intercept])(np.unique(torsionDict[key]["wbo"])),
                c=colors[i],
            )

    l = ax.legend(bbox_to_anchor=(1, 1), fontsize=14)

    plt.xlabel("Wiberg Bond Order", fontsize=14)
    plt.ylabel("Torsion barrier (kJ/mol)", fontsize=14)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.savefig(fname, bbox_inches="tight")
    plt.show()


def visualize_wbo_correlation_compare(fileName, fileName2, fname):
    """
    Generates a .pdf plot from a .oeb file of wbo versus torsion barrier height for two datasets

        Parameters
        ----------
        fileName : String
            .oeb file name for the molecules that contain datatags with plotting information
        fileName2 : String
            .oeb file name for the molecules that contain datatags with plotting information. This dataset
            must have less torsion parametrs than fileName2
        fname : String
            The output .pdf file name for the resulting files from the plotting function

        Returns
        -------
    """
    molList = oeb2oemol(fileName)
    torsionDict = {}
    torsionDict2 = {}

    for m in molList:
        tid = m.GetData("IDMatch")
        torsionDict[tid] = {}
        torsionDict[tid]["tb"] = []
        torsionDict[tid]["wbo"] = []
        torsionDict2[tid] = {}
        torsionDict2[tid]["tb"] = []
        torsionDict2[tid]["wbo"] = []

    for m in molList:
        tid = m.GetData("IDMatch")
        torsionDict[tid]["tb"].append(m.GetData("TB"))
        torsionDict[tid]["wbo"].append(m.GetData("WBO"))

    molList2 = oeb2oemol(fileName2)

    for m in molList2:
        tid = m.GetData("IDMatch")
        torsionDict2[tid] = {}
        torsionDict2[tid]["tb"] = []
        torsionDict2[tid]["wbo"] = []

    for m in molList2:
        tid = m.GetData("IDMatch")
        torsionDict2[tid]["tb"].append(m.GetData("TB"))
        torsionDict2[tid]["wbo"].append(m.GetData("WBO"))

    colors = fragmenter.chemi._KELLYS_COLORS

    fig, ax = plt.subplots()

    for i, (key, item) in enumerate(torsionDict.items()):
        ax.scatter(
            torsionDict[key]["wbo"],
            torsionDict[key]["tb"],
            c=colors[i],
            label=key,
            marker="o",
        )
        ax.scatter(
            torsionDict2[key]["wbo"], torsionDict2[key]["tb"], c=colors[i], marker="x"
        )
        slope, intercept, r_value, p_value, std_err = stats.linregress(
            torsionDict[key]["wbo"] + torsionDict2[key]["wbo"],
            torsionDict[key]["tb"] + torsionDict2[key]["tb"],
        )
        print(key)
        print(r_value)
        print(slope)
        print(intercept)
        ci_r_value = arch.bootstrap.IIDBootstrap(
            np.asarray(torsionDict[key]["wbo"] + torsionDict2[key]["wbo"]),
            np.asarray(torsionDict[key]["tb"] + torsionDict2[key]["tb"]),
        ).conf_int(compute_r_ci, 1000)
        CI_95 = 1.96 * std_err
        plt.plot(
            np.unique(torsionDict[key]["wbo"] + torsionDict2[key]["wbo"]),
            np.poly1d([slope, intercept])(
                np.unique(torsionDict[key]["wbo"] + torsionDict2[key]["wbo"])
            ),
            c=colors[i],
        )

    l = ax.legend(bbox_to_anchor=(1, 1), fontsize=14)
    plt.xlabel("Wiberg Bond Order", fontsize=14)
    plt.ylabel("Energy Barrier height (kJ/mol)", fontsize=14)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.savefig(fname, bbox_inches="tight")
    plt.show()


# main dataset generation function
def genData(dsName, fileName, ff_name, benchmark_smiles="lim_mobley_parsley_benchmark.smi", qca_overlapped_entries="qca_overlapped_entries.txt"):
    """
    Generates oeb files for the QCA datasets and plots that analyze WBO versus torsion barrier
    """
    molList = loadDataset_low(dsName, "default", benchmark_smiles, qca_overlapped_entries)
    mols1 = checkTorsion(molList, ff_name)
    makeOEB(mols1, fileName)


def plot_tid_for_datasets(fileList, t_id):
    """
    Takes in a list of oeb files and plots wbo vs torsion barrier, combining all the datasets
    and plotting by each tid in the combined dataset

    Parameters
    ----------
    fileList: list of strings
    each string is a oeb file name
    Eg. ['rowley.oeb'] or ['rowley.oeb', 'phenyl.oeb']

    t_id: str
    torsion id, eg., 't43'
    """
    import base64
    import ntpath
    from io import BytesIO

    import pandas as pd
    import plotly.express as px
    import plotly.graph_objects as go
    from plotly.validators.scatter.marker import SymbolValidator
    from rdkit import Chem
    from rdkit.Chem.Draw import MolsToGridImage

    df = pd.DataFrame(columns = ['tid', 'tb', 'wbo', 'cmiles', 'TDindices', 'filename']) 
    fig = go.Figure({'layout' : go.Layout(height=900, width=1000,
            xaxis={'title': 'Wiberg Bond Order'},
            yaxis={'title': 'Torsion barrier (kJ/mol)'},
            #paper_bgcolor='white',
            plot_bgcolor='rgba(0,0,0,0)',
            margin={'l': 40, 'b': 40, 't': 10, 'r': 10},
            legend={'orientation': 'h', 'y': -0.2}, 
            legend_font=dict(family='Arial', color='black', size=15),
            hovermode=False,
            dragmode='select')})
    fig.update_xaxes(title_font=dict(size=26, family='Arial', color='black'), 
                     ticks="outside", tickwidth=2, tickcolor='black', ticklen=10,
                     tickfont=dict(family='Arial', color='black', size=20),
                     showgrid=False, gridwidth=1, gridcolor='black', 
                     mirror=True, linewidth=2, linecolor='black', showline=True)
    fig.update_yaxes(title_font=dict(size=26, family='Arial', color='black'),
                    ticks="outside", tickwidth=2, tickcolor='black', ticklen=10,
                    tickfont=dict(family='Arial', color='black', size=20),
                    showgrid=False, gridwidth=1, gridcolor='black', 
                    mirror=True, linewidth=2, linecolor='black', showline=True)
    

    colors = fragmenter.chemi._KELLYS_COLORS
    colors = colors * 2
    raw_symbols = SymbolValidator().values
    symbols = []
    for i in range(0,len(raw_symbols),8):
        symbols.append(raw_symbols[i])
    count = 0
    
    fname = []
    for fileName in fileList:
        molList = []
        fname = fileName
        molList = oeb2oemol(fname)

        for m in molList:
            tid = m.GetData("IDMatch")
            fname = ntpath.basename(fileName)
            df = df.append({'tid': tid, 
                            'tb': m.GetData("TB"),
                            'wbo' : m.GetData("WBO"),
                            'cmiles' : m.GetData("cmiles"),
                            'TDindices' : m.GetData("TDindices"),
                            'filename' : fname}, 
                            ignore_index = True)
        
        x = df[(df.filename == fname) & (df.tid == t_id)].wbo
        y = df.loc[x.index].tb
        fig.add_scatter(x=x,
                        y=y,
                        mode="markers", 
                        name=fname, 
                        marker_color=colors[count],
                        marker_symbol=count,
                        marker_size=13)
        count += 1
    
    x = df[df.tid == t_id].wbo
    y = df.loc[x.index].tb
    slope, intercept, r_value, p_value, std_err =    stats.linregress(x, y)
    print("tid: ", t_id, "r_value: ", r_value, 
          "slope: ", slope, "intercept: ", intercept)

    fig.add_traces(go.Scatter(
        x=np.unique(x), 
        y=np.poly1d([slope, intercept])(np.unique(x)), 
        showlegend=False, mode ='lines'))
    slope_text = 'slope: '+str('%.2f' % slope)
    r_value = 'r_val: '+str('%.2f' % r_value)
    fig_text = slope_text + ', '+ r_value
    fig.add_annotation(text=fig_text, 
                       font = {'family': "Arial", 'size': 22, 'color': 'black'},
                       xref="paper", yref="paper", x=1, y=1,
                       showarrow=False)

    fig.update_layout(
        title={"text": t_id, "x": 0.5, "xanchor": "center", "yanchor": "top"}
    )
    fig.show()
    return fig
