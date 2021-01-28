"""
Script to compare bond order values between for AmberToolkit and OpenEyeToolit for mols
in a .smi database

Usage:
    python single_conformer.py --smiles_database [DATABASE.SMI] --n [NUM MOLS]
"""

import argparse
import copy
import pickle
import re
from multiprocessing.pool import Pool
from typing import Optional, Tuple, Dict

import numpy
from openeye import oechem, oeomega
from openforcefield.topology import Molecule
from openforcefield.utils import OpenEyeToolkitWrapper, AmberToolsToolkitWrapper
from tqdm import tqdm


def generate_conformer(oe_molecule: oechem.OEMol) -> oechem.OEMol:
    """Generates a single conformer for a given molecule"""

    oe_molecule = oechem.OEMol(oe_molecule)

    # Guess the stereochemistry if it is not already perceived
    unspecified_stereochemistry = any(
        entity.IsChiral() and not entity.HasStereoSpecified()
        for entity in [*oe_molecule.GetAtoms(), *oe_molecule.GetBonds()]
    )

    if unspecified_stereochemistry:
        stereoisomer = next(iter(oeomega.OEFlipper(oe_molecule.GetActive(), 12, True)))
        oe_molecule = oechem.OEMol(stereoisomer)

    # Otherwise, generate the conformer.
    omega_options = oeomega.OEOmegaOptions(oeomega.OEOmegaSampling_Sparse)

    omega = oeomega.OEOmega(omega_options)
    omega.SetIncludeInput(False)
    omega.SetCanonOrder(False)

    output_stream = oechem.oeosstream()

    oechem.OEThrow.SetOutputStream(output_stream)
    oechem.OEThrow.Clear()

    status = omega(oe_molecule)

    oechem.OEThrow.SetOutputStream(oechem.oeerr)

    output_string = output_stream.str().decode("UTF-8")

    output_string = output_string.replace("Warning: ", "")
    output_string = re.sub("^: +", "", output_string, flags=re.MULTILINE)
    output_string = re.sub("\n$", "", output_string)

    if not status:
        raise RuntimeError(f"Could not generate conformer: {output_string}")

    oe_conformer = [*oe_molecule.GetConfs()][0]

    conformer = numpy.zeros((oe_molecule.NumAtoms(), 3))

    for atom_index, coordinates in oe_conformer.GetCoords().items():
        conformer[atom_index, :] = coordinates

    oe_molecule.DeleteConfs()
    oe_molecule.NewConf(oechem.OEFloatArray(conformer.flatten()))

    return oe_molecule


def process_molecule(
    oe_molecule: oechem.OEMol
) -> Tuple[Dict[str, Molecule], Optional[str]]:
    """
    Processes every molecule, generating bond orders using AmberToolkit and OpenEyeToolkit Wrappers
    """
    
    error = None
    smiles = None

    try:

        smiles = oechem.OEMolToSmiles(oe_molecule)

        # Generate a single conformer for the molecule.
        oe_molecule = generate_conformer(oe_molecule)

        # Map to an OpenFF molecule object.
        molecule: Molecule = Molecule.from_openeye(oe_molecule)

        # Compute the AM1 partial charges and WBOs
        molecules = {
            "openeye": copy.deepcopy(molecule),
            "ambertools": copy.deepcopy(molecule),
        }

        for charge_backend, molecule in molecules.items():

            toolkit_wrapper = (
                OpenEyeToolkitWrapper() if charge_backend == "openeye"
                else AmberToolsToolkitWrapper()
            )

            print(f"Generating {charge_backend.upper()} Charges.")
            molecule.compute_partial_charges_am1bcc(
                use_conformers=molecule.conformers,
                toolkit_registry=toolkit_wrapper
            )

            print(f"Generating {charge_backend.upper()} WBO.")
            molecule.assign_fractional_bond_orders(
                "am1-wiberg",
                use_conformers=molecule.conformers,
                toolkit_registry=toolkit_wrapper
            ) #something goin on here

    except (BaseException, Exception) as e:
        molecules = {}
        error = f"Failed to process {smiles}: {str(e)}"

    return molecules, error


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--smiles_database",
                        type=str,
                        required=True,
                        help=("Name of a SMILES file containing molecules"))
    parser.add_argument("--n",
                        type=int,
                        required=True,
                        help=("Number of molecules in the .smi file"))
    args = parser.parse_args()
    
    n_processes = 32
    molecule_set = args.smiles_database[:-4]

    input_molecule_stream = oechem.oemolistream()
    input_molecule_stream.open(f"{molecule_set}.smi")

    with Pool(processes=n_processes) as pool:

        processed_molecules = list(
            tqdm(
                pool.imap(process_molecule, input_molecule_stream.GetOEMols()),
                total=args.n,
            )
        )

    # Retain only the molecules which could be processed and print errors
    # to a log file.
    molecules = {"openeye": [], "ambertools": []}

    with open(f"results/{molecule_set}-errors.log", "w") as file:

        for charged_molecules, error in processed_molecules:

            if error is not None:

                file.write("=".join(["="] * 40) + "\n")
                file.write(error + "\n")
                file.write("=".join(["="] * 40) + "\n")

                continue

            for charge_backend in charged_molecules:
                molecules[charge_backend].append(charged_molecules[charge_backend])

    for charge_backend in molecules:
        with open(f"results/{molecule_set}-{charge_backend}.pkl", "wb") as file:
            pickle.dump(molecules[charge_backend], file)

    input_molecule_stream.close()


if __name__ == "__main__":
    main()
