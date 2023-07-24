""" Utilities for interaction with the ASE library. """

from ase import Atoms

from excitingtools import ExcitingStructure
from excitingtools.constants.units import bohr_to_angstrom


def exciting_structure_to_ase(structure: ExcitingStructure) -> Atoms:
    """ Function to extract the physical structure from an exciting structure object
    and transforms it into an ase.atoms.Atoms object.

    :param structure: input exciting structure object
    :returns: ASE Atoms object
    """
    lattice = structure.get_lattice(convert_to_angstrom=True)
    if structure.structure_attributes.get("cartesian"):
        positions = structure.positions * bohr_to_angstrom
        return Atoms(symbols=structure.species, positions=positions, cell=lattice, pbc=True)
    else:
        return Atoms(symbols=structure.species, scaled_positions=structure.positions, cell=lattice, pbc=True)
