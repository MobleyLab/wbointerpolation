import base64
import json
import os
from io import BytesIO
from typing import Tuple

import matplotlib.pyplot as plt
import mdtraj
import numpy as np
import pandas as pd
from forcebalance.molecule import Molecule as fb_molecule
from forcebalance.openmmio import MTSVVVRIntegrator
from openff.toolkit.topology import Molecule, Topology
from openff.toolkit.typing.engines.smirnoff.forcefield import ForceField
from scipy import optimize
from simtk import openmm, unit


def get_MM_energy(system, positions):
    """
    Return the single point potential energy.

    Parameters
    ----------
    system : simtk.openmm.System
        The system to check
    positions : simtk.unit.Quantity of dimension (natoms,3) with units of length
        The positions to use

    Returns
    ---------
    energy : simtk.unit.Quantity
        sum of total energy
    """

    integrator = MTSVVVRIntegrator(
        300 * unit.kelvin,
        1 / unit.picoseconds,
        1.0 * unit.femtoseconds,
        system,
        ninnersteps=int(1 / 0.25),
    )
    # openmm.VerletIntegrator(1.0 * unit.femtoseconds)
    context = openmm.Context(system, integrator)
    context.setPositions(positions)
    state = context.getState(getEnergy=True)
    energy = state.getPotentialEnergy().in_units_of(unit.kilocalorie_per_mole)
    return energy


def show_oemol_struc(oemol, torsions=False, atom_indices=[]):
    """
    Returns the oedepict image with or without the torsion highlighted

    Parameters
    ----------
    oemol: openeye oechem mol object
    torsions: boolean, to highlight dihedrals
    atom_indices: dihedral atom indices to highlight

    Returns
    -------
    Image: image in png format
    """
    from IPython.display import Image
    from openeye import oechem, oedepict

    width = 400
    height = 300

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


def get_wbo(offmol, dihedrals):
    """
    Returns the specific wbo of the dihedral bonds
    Parameters
    ----------
    offmol: openforcefield molecule
    dihedrals: list of atom indices in dihedral

    Returns
    -------
    bond.fractional_bond_order: wiberg bond order calculated using openforcefield toolkit for the specific dihedral central bond
    """
    offmol.assign_fractional_bond_orders(bond_order_model="am1-wiberg-elf10")
    bond = offmol.get_bond_between(dihedrals[1], dihedrals[2])
    return bond.fractional_bond_order


def get_torsion_barrier(energies, angles):
    """
    gets the torsion barrier as the maximum barrier among the barriers present in the torsion profile
    Parameters
    ----------
    energies: list of energies along a torsion profile, can be MM or QM (y-axis values)
    angles: list of corresponding angles (x-axis values)

    Returns
    -------
    torsion_barriers: torsion barrier in the provided energies units
    """
    idx = []
    flat_list = [item for sublist in angles for item in sublist]
    angles = np.array(flat_list)
    angles = angles * np.pi / 180
    # appending the first three and last three energies to the end and beginning energies so as to have a continuous curve and all barriers at the edge are also covered
    angles = np.append(
        angles[-3:] - 2 * np.pi, np.append(angles, angles[:3] + 2 * np.pi)
    )
    energies = np.append(energies[-3:], np.append(energies, energies[:3]))

    for i in range(len(angles) - 2):
        m1 = (energies[i + 1] - energies[i]) / (angles[i + 1] - angles[i])
        m2 = (energies[i + 2] - energies[i + 1]) / (angles[i + 2] - angles[i + 1])
        if np.sign(m1) == np.sign(m2):
            continue
        else:
            idx.append(i + 1)
    torsion_barriers = []
    for i in range(int(len(idx) - 1)):
        torsion_barriers.append(abs(energies[idx[i]] - energies[idx[i + 1]]))
        # 4.184 KJ/mol = 1 kcal/mol

    torsion_barriers = np.array(torsion_barriers)
    return max(torsion_barriers)


def zero_dihedral_contribution(openmm_system: openmm.System, dihedral_indices: Tuple[int, int]):
    """
    Author Simon Boothroyd
    Link gist.github.com/SimonBoothroyd/667b50314c628aabe5064f0defb6ad8e
    Zeroes out the contributions of a particular dihedral to the total potential
    energy for a given OpenMM system object.

    Parameters
    ----------
    openmm_system:
        The OpenMM system object which will be modified so the specified dihedral will
        not contribute to the total potential energy of the system.
    dihedral_indices:
        The indices of the two central atoms in the dihedral whose contributions should
        be zeroed.
    """
    torsion_forces = [
        force
        for force in openmm_system.getForces()
        if isinstance(force, openmm.PeriodicTorsionForce)
    ]
    for torsion_force in torsion_forces:
        for torsion_index in range(torsion_force.getNumTorsions()):
            (
                index_i,
                index_j,
                index_k,
                index_l,
                periodicity,
                phase,
                k,
            ) = torsion_force.getTorsionParameters(torsion_index)
            if (
                    index_j not in [dihedral_indices[0], dihedral_indices[1]]
                    or index_k not in [dihedral_indices[0], dihedral_indices[1]]
            ):
                continue
            torsion_force.setTorsionParameters(
                torsion_index,
                index_i,
                index_j,
                index_k,
                index_l,
                periodicity,
                phase,
                0.0
            )


def build_restrained_context(offxml, molfile, traj, dihedrals):
    """
    Builds openmm context for a full MM calculation keeping the dihedral frozen
    and applying a restraint of 1 kcal/mol on atoms not involved in torsion
    Parameters
    ----------
    offxml: forcefield file
    molfile: molecule file in sdf format
    positions: list of positions of all conformers
    dihedrals: list of atom indices in the dihedral

    Returns
    -------
    list of openmm contexts
    """
    contexts = []
    forcefield = ForceField(offxml, allow_cosmetic_attributes=True)
    molecule = Molecule.from_file(molfile, allow_undefined_stereo=True)

    for pos in traj:
        system = forcefield.create_openmm_system(molecule.to_topology())
        add_restraints_to_system(system, pos * unit.angstrom, dihedrals)
        for force in system.getForces():
            if force.__class__.__name__ == "CustomExternalForce":
                force.setForceGroup(11)
        for i in dihedrals:
            system.setParticleMass(i, 0.0)
        integrator = openmm.LangevinIntegrator(
            300.0 * unit.kelvin, 1.0 / unit.picosecond, 2.0 * unit.femtosecond
        )
        platform = openmm.Platform.getPlatformByName("Reference")
        contexts.append(openmm.Context(system, integrator, platform))

    return contexts


def build_context_residual(offxml, molfile, dihedrals):
    """Builds context for the openmm calculation keeping the dihedral frozen and the
    specific energy contributions from the dihedral zeroed out
    Parameters
    ----------
    offxml: forcefield file
    molfile: molecule file in sdf format
    dihedrals: list of atom indices in the dihedral
    Returns
    -------
    """
    forcefield = ForceField(offxml, allow_cosmetic_attributes=True)
    molecule = Molecule.from_file(molfile, allow_undefined_stereo=True)
    system = forcefield.create_openmm_system(molecule.to_topology())

    # freeze the dihedral atoms so that the dihedral angle stays fixed
    for i in dihedrals:
        system.setParticleMass(i, 0.0)
    # zero out the energy contributions from the dihedral atoms
    zero_dihedral_contribution(system, (dihedrals[1], dihedrals[2]))
    integrator = openmm.LangevinIntegrator(
        300.0 * unit.kelvin, 1.0 / unit.picosecond, 2.0 * unit.femtosecond
    )
    platform = openmm.Platform.getPlatformByName("Reference")
    context = openmm.Context(system, integrator, platform)
    return context


def add_restraints_to_system(system, positions, dihedrals):
    k = 1.0 * unit.kilocalories_per_mole / unit.angstroms ** 2

    for i in range(system.getNumParticles()):
        if i not in dihedrals:
            positional_restraint = openmm.CustomExternalForce(
                "0.5 * k * ((x-x0)^2 + (y-y0)^2 + (z-z0)^2)"
            )
            positional_restraint.addPerParticleParameter("k")
            positional_restraint.addPerParticleParameter("x0")
            positional_restraint.addPerParticleParameter("y0")
            positional_restraint.addPerParticleParameter("z0")

            position = positions[i]
            positional_restraint.addParticle(
                i,
                [
                    k,
                    position[0].value_in_unit(unit.nanometers),
                    position[1].value_in_unit(unit.nanometers),
                    position[2].value_in_unit(unit.nanometers),
                ],
            )
            system.addForce(positional_restraint)


def evaluate_minimized_energies(context, trajectory, dihedrals):
    """
    Evaluate minimized energies for all the geometries in the trajectory in given openmm context
    Parameters
    ----------
    context: openmm context
    trajectory: list of (list of xyz coordinates) for different points on the trajectory, usually for different points on the torsion profile

    Returns
    -------
    list of energies in kcal/mol
    """
    energies = []
    min_pos = []
    rmsd_pos = []

    for frame_geo in trajectory:
        initial_positions = frame_geo * unit.angstrom
        context.setPositions(initial_positions)
        # apply restraint of 1 kcal/mol per ang^2 to non-torsion atoms
        #         add_restraints_to_system(context.getSystem(), frame_geo * unit.angstrom, dihedrals)
        openmm.LocalEnergyMinimizer_minimize(
            context, tolerance=5.0e-09, maxIterations=15000
        )
        energy = context.getState(getEnergy=True).getPotentialEnergy()
        energy = energy.value_in_unit(unit.kilocalories_per_mole)
        minimized_positions = (
            context.getState(getPositions=True)
            .getPositions(asNumpy=True)
            .value_in_unit(unit.angstroms)
        )
        rmsd_pos.append(
            get_rmsd(
                initial_positions.value_in_unit(unit.angstroms), minimized_positions
            )
        )
        min_pos.append(minimized_positions)
        energies.append(energy)

    return np.array(energies), min_pos, rmsd_pos


def evaluate_restrained_energies(contexts, trajectory, dihedrals):
    """
    Evaluate minimized energies for all the geometries in the trajectory in given openmm context
    Parameters
    ----------
    context: openmm context
    trajectory: list of (list of xyz coordinates) for different points on the trajectory, usually for different points on the torsion profile

    Returns
    -------
    list of energies in kcal/mol
    """
    energies = []
    min_pos = []
    rmsd_pos = []

    for context, frame_geo in zip(contexts, trajectory):
        initial_positions = frame_geo * unit.angstrom
        context.setPositions(initial_positions)
        openmm.LocalEnergyMinimizer_minimize(
            context, tolerance=5.0e-09, maxIterations=1500
        )
        # Remove the restraint energy from the total energy
        groups = set(range(32))
        #         frc = context.getSystem().getForce(11)
        groups.remove(11)
        energy = context.getState(getEnergy=True, groups=groups).getPotentialEnergy()
        energy = energy.value_in_unit(unit.kilocalories_per_mole)
        minimized_positions = (
            context.getState(getPositions=True)
            .getPositions(asNumpy=True)
            .value_in_unit(unit.angstroms)
        )
        rmsd_pos.append(
            get_rmsd(
                initial_positions.value_in_unit(unit.angstroms), minimized_positions
            )
        )
        min_pos.append(minimized_positions)
        energies.append(energy)

    return np.array(energies), min_pos, rmsd_pos


def plot_energies_data(energies_data_dict, title, ylab):
    """
    Build the torsion profile plots from the energy dict for all the forcefields or energy data tagged per key in the dict
    Parameters
    ----------
    energies_data_dict: dict, energy key and the associated list of energies, eg., 'QM': [e1, e2, e3,...]
    title: string, plot title
    ylab: string, label on the y-axis

    Returns
    -------
    plot in iobytes/string format
    """
    CB_color_cycle = ["#0072B2", "#009E73", "#D55E00", "#F0E442", "#CC79A7", "#56B4E9"]
    marks = ["P", "o", "^", "X"]
    plt.Figure(figsize=(3, 3))
    x_axis = energies_data_dict.pop("td_angles")
    count = 0
    for dataname, datavalues in energies_data_dict.items():
        if dataname == "QM":
            plt.plot(
                x_axis,
                datavalues,
                color=CB_color_cycle[count],
                linestyle="solid",
                marker="D",
                label=dataname,
                markersize=7,
            )
        else:
            if ylab == "res":
                dataname = "QM - " + dataname + "_intrinsic"
            plt.plot(
                x_axis,
                datavalues,
                color=CB_color_cycle[count],
                linestyle="dashed",
                marker=marks[count - 1],
                label=dataname,
                markersize=(14 - 3 * count),
            )
        count = count + 1
    plt.legend()
    if ylab == "res":
        plt.title(title + "  Residuals")
        plt.ylabel("Residuals (QM - MM_intrinsic) [ kcal/mol ]")
        plt.xlabel("Torsion Angle [ degree ]")
    elif ylab == "rel":
        plt.title(title + "  Relative energies")
        plt.ylabel("Relative energies (QM Vs MM_full) [ kcal/mol ]")
        plt.xlabel("Torsion Angle [ degree ]")
    plot_iobytes = BytesIO()
    plt.savefig(plot_iobytes, format="png", dpi=75)
    plt.close()
    return plot_iobytes


def fig2inlinehtml(fig):
    """
    small hack to convert image data to image string for html display
    Parameters
    ----------
    fig: image to be converted

    Returns
    -------
    imgestr: image in string format
    """
    figfile = fig.data
    # for python 3.x:
    figdata_png = base64.b64encode(figfile).decode()
    imgstr = '<img src="data:image/png;base64,{}" />'.format(figdata_png)
    return imgstr


def get_assigned_torsion_param(molecule, ff, dihedrals):
    """
    for a molecule and specific dihedral check the assigned torsion parameter
    Parameters
    ----------
    molecule: openforcefield molecule object
    ff: ForceField offxml file
    dihedrals: list of atom indices in the dihedral

    Returns
    -------
    parameter.id: str of the torsion parameter associated with the dihedral
    """
    topology = Topology.from_molecules([molecule])
    # Run the molecule labeling
    forcefield = ForceField(ff, allow_cosmetic_attributes=True)
    molecule_force_list = forcefield.label_molecules(topology)

    # Print out a formatted description of the parameters applied to this molecule
    for mol_idx, mol_forces in enumerate(molecule_force_list):
        for force_tag, force_dict in mol_forces.items():
            if force_tag == "ProperTorsions":
                for (atom_indices, parameter) in force_dict.items():
                    if atom_indices == tuple(dihedrals) or tuple(
                        reversed(atom_indices)
                    ) == tuple(dihedrals):
                        return parameter.id


def create_energy_dataframe(subdirs, offxml_list, molfile_format, qdata_filename):
    """
    Creates a dataframe with the entries listed under cols below

    """
    cols = [
        "Torsion ID",
        "assigned params",
        "wbo",
        "Max - min for (QM - [MM_SF3, MM_1.3.0]) kcal/mol",
        "Chemical Structure",
        "QM-MM_intrinsic relative energies",
        "QM Vs MM_full",
        "rmsd",
    ]
    df = pd.DataFrame(columns=cols)
    ff_names = [os.path.splitext(os.path.basename(f))[0] for f in offxml_list]

    count = 0
    for mol_folder in subdirs:
        count += 1
        #     if(count > 10):
        #         break
        mol_folder = mol_folder.rstrip()
        if count % 10 == 0:
            print("processed ", count)
        count += 1
        mol_file = mol_folder + "/" + molfile_format
        mol = Molecule.from_file(
            mol_folder + "/" + molfile_format, allow_undefined_stereo=True
        )
        # build openmm contexts for each force field file

        # find scan trajectories
        traj_files = [
            f for f in os.listdir(mol_folder) if os.path.splitext(f)[-1] == ".xyz"
        ]
        # create output subfolder
        mol_name = os.path.basename(mol_folder)
        with open(mol_folder + "/metadata.json") as json_data:
            metadata = json.load(json_data)
        dihedrals = metadata["dihedrals"][0]
        assigned_torsion_params = dict(
            (
                os.path.splitext(os.path.basename(key))[0],
                get_assigned_torsion_param(mol, key, dihedrals),
            )
            for key in offxml_list
        )
        #     {key: value for (key, value) in iterable}[ for x in offxml_list]

        #   zero out the assigned torsion parameter for this particular dihedral to get the intrinsic torsion potential
        #   "(full QM energy) - (full MM energy, locally optimized, with that specific torsion zeroed out)"
        context_res = []
        for i in range(len(offxml_list)):
            key = os.path.splitext(os.path.basename(offxml_list[i]))[0]
            context_res.append(
                build_context_residual(offxml_list[i], mol_file, dihedrals)
            )

        #         contexts_full = [build_context_full(offxml, mol_file, dihedrals) for offxml in offxml_list]

        wbo = get_wbo(mol, dihedrals)
        for f in traj_files:
            #         print(f"- {f}")
            # hold energy data for this trajectory
            energies_data_dict = {}
            energies_full_dict = {}
            tb_dict = {}
            # use ForceBalance.molecule to read the xyz file
            fb_mol = fb_molecule(os.path.join(mol_folder, f))
            # read torsion angles
            with open(mol_folder + "/metadata.json") as json_data:
                metadata = json.load(json_data)
            energies_data_dict["td_angles"] = metadata["torsion_grid_ids"]
            # read QM energies
            qdata = np.genfromtxt(
                str(mol_folder) + "/" + str(qdata_filename),
                delimiter="ENERGY ",
                dtype=None,
            )
            # print(qdata.shape)
            eqm = qdata[:, 1]
            # record the index of the minimum energy structure
            ground_idx = np.argmin(eqm)
            eqm = [627.509 * (x - eqm[ground_idx]) for x in eqm]
            energies_data_dict["QM"] = eqm
            energies_full_dict["QM"] = eqm
            tb_dict["QM"] = get_torsion_barrier(
                energies_data_dict["QM"], energies_data_dict["td_angles"]
            )
            #             print("WBO of conformers")

            for i in range(len(fb_mol.xyzs)):
                g_frame = fb_mol.xyzs[i]
                pos = np.array(g_frame) * unit.angstrom
                mol.assign_fractional_bond_orders(
                    bond_order_model="am1-wiberg-elf10"
                )  # , use_conformers=[pos])
                bond = mol.get_bond_between(dihedrals[1], dihedrals[2])
            #                 print(metadata["torsion_grid_ids"][i], bond.fractional_bond_order)

            ### FULL
            # evalute full mm energies for each force field
            rmsd_dict = {}
            energies_full_dict = energies_data_dict.copy()
            for offxml, ffname in zip(offxml_list, ff_names):
                contextf = build_restrained_context(
                    offxml, mol_file, fb_mol.xyzs, dihedrals
                )
                mm_energies, minimized_pos, rmsd = evaluate_restrained_energies(
                    contextf, fb_mol.xyzs, dihedrals
                )
                # mm_energies -= mm_energies[ground_idx]
                mm_energies -= mm_energies.min()
                key_name = ffname + "_full"
                energies_full_dict[ffname] = mm_energies
                rmsd_dict[ffname] = rmsd
                tb_dict[key_name] = get_torsion_barrier(
                    energies_full_dict[ffname], energies_full_dict["td_angles"]
                )

            ### INTRINSIC
            # evalute intrinsic mm energies for each force field
            for context, ffname in zip(context_res, ff_names):
                mm_energies, _, _ = evaluate_minimized_energies(
                    context, fb_mol.xyzs, dihedrals
                )
                # mm_energies -= mm_energies[ground_idx]
                mm_energies -= mm_energies.min()
                energies_data_dict[ffname] = mm_energies
                energies_data_dict[ffname] = (
                    energies_data_dict["QM"] - energies_data_dict[ffname]
                )
                key_name = ffname + "_residual"
                tb_dict[key_name] = get_torsion_barrier(
                    energies_data_dict[ffname], energies_data_dict["td_angles"]
                )

            plot_iobytes_full = plot_energies_data(energies_full_dict, mol_name, "rel")
            figdata_png_full = base64.b64encode(plot_iobytes_full.getvalue()).decode()
            plot_str_full = '<img src="data:image/png;base64,{}" />'.format(
                figdata_png_full
            )

            plot_iobytes = plot_energies_data(energies_data_dict, mol_name, "res")
            figdata_png = base64.b64encode(plot_iobytes.getvalue()).decode()
            plot_str = '<img src="data:image/png;base64,{}" />'.format(figdata_png)
            oemol = mol.to_openeye()
            dihedral_indices = metadata["dihedrals"][0]
            struc_str = fig2inlinehtml(
                show_oemol_struc(oemol, torsions=True, atom_indices=dihedral_indices)
            )
            info_dict = {}
            tid = mol_name[4:]
            info_dict["assigned_params"] = assigned_torsion_params

            df = df.append(
                {
                    cols[0]: tid,
                    cols[1]: info_dict["assigned_params"],
                    cols[2]: wbo,
                    cols[3]: tb_dict,
                    cols[4]: struc_str,
                    cols[5]: plot_str,
                    cols[6]: plot_str_full,
                    cols[7]: rmsd_dict,
                },
                ignore_index=True,
            )
    return df


def fix_dihedral(molecule, system, dihedrals):
    """
    Author: Simon Boothroyd

    impose a dihedral angle constraint so that bonds and angles can change during optimization
    """
    mdtraj_trajectory = mdtraj.Trajectory(
        xyz=molecule.conformers[0],
        topology=mdtraj.Topology.from_openmm(molecule.topology.to_openmm()),
    )
    dihedral_angle = mdtraj.compute_dihedrals(mdtraj_trajectory, np.array([dihedrals]))[
        0
    ][0].item()
    dihedral_restraint = openmm.CustomTorsionForce(
        f"k * min(min(abs(theta - theta_0), abs(theta - theta_0 + 2 * "
        f"{np.pi})), abs(theta - theta_0 - 2 * {np.pi}))^2"
    )
    dihedral_restraint.addPerTorsionParameter("k")
    dihedral_restraint.addPerTorsionParameter("theta_0")
    theta_0 = dihedral_angle
    k = 1.0 * unit.kilocalories_per_mole / unit.radian ** 2
    dihedral_restraint.addTorsion(
        dihedrals[0],
        dihedrals[1],
        dihedrals[2],
        dihedrals[3],
        [k, theta_0],
    )
    system.addForce(dihedral_restraint)


def target_function(xyz, context):
    """
    for using different optimizers from scipy defining the target function as change in positions
    """
    context.setPositions(xyz.reshape(-1, 3))
    state = context.getState(getEnergy=True, getForces=True)
    frc = state.getForces(asNumpy=True)
    ene = state.getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole)
    frc = frc.value_in_unit(unit.kilocalorie_per_mole / unit.nanometer)
    return ene, -frc.flatten()


def use_scipy_opt(pos, context):
    """
    using a scipy optimization method for openmm optimizations
    """
    pos = pos.value_in_unit(unit.nanometer)
    results = optimize.minimize(target_function, pos, context, method="SLSQP")
    # L-BFGS-B',                            jac=True, options=dict(maxiter=20000, disp=True))
    if results.success:
        print(results.nit, ",   ", results.fun)
    else:
        print("didn't converge for mol")
    return results.fun * unit.kilocalorie_per_mole


def get_rmsd(initial, final):
    """
    Evaluate the RMSD between two arrays
    """
    assert len(initial) == len(final)
    n = len(initial)
    if n == 0:
        return 0.0
    diff = np.subtract(initial, final)
    rmsd = np.sqrt(np.sum(diff ** 2) / n)
    return rmsd
