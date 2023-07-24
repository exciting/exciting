""" Utilities to interact with the pymatgen library. """

from pymatgen.core import Structure

from excitingtools import ExcitingStructure
from excitingtools.constants.units import angstrom_to_bohr, bohr_to_angstrom


def exciting_structure_to_pymatgen(structure: ExcitingStructure) -> Structure:
    """ Function to extract the physical structure from an exciting structure object
    and transforms it into a pymatgen.core.structure.Structure object.

    :param structure: input exciting structure object
    :returns: pymatgen Structure object
    """
    lattice = structure.get_lattice(convert_to_angstrom=True)
    cartesian = structure.structure_attributes.get("cartesian", False)
    positions = structure.positions * bohr_to_angstrom if cartesian else structure.positions
    return Structure(lattice=lattice, species=structure.species, coords=positions, coords_are_cartesian=cartesian)


def pymatgen_to_exciting_structure(structure: Structure) -> ExcitingStructure:
    """ Initialise lattice, species and positions from a pymatgen Structure Object.
    Note: pymatgen works in Angstrom, whereas exciting expects atomic units

    :param structure: pymatgen Structure object.
    :return exciting structure object
    """
    lattice = structure.lattice.matrix * angstrom_to_bohr
    species = [x.symbol.capitalize() for x in structure.species]
    positions = structure.frac_coords
    atoms = [{"species": atom, "position": positions[i]} for i, atom in enumerate(species)]
    return ExcitingStructure(atoms, lattice)
