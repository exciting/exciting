"""Structure class, mirroring that of exciting's structure XML sub-tree.
https://exciting.wikidot.com/ref:structure
"""
from __future__ import annotations

import copy
from pathlib import Path
from typing import Optional, Union, List, Dict
from xml.etree import ElementTree

import numpy as np

from excitingtools.constants.units import angstrom_to_bohr, bohr_to_angstrom
from excitingtools.input.base_class import ExcitingXMLInput, AbstractExcitingInput
from excitingtools.input.input_classes import ExcitingCrystalInput, ExcitingSpeciesInput
from excitingtools.structure.lattice import check_lattice, check_lattice_vector_norms
from excitingtools.utils import valid_attributes
from excitingtools.utils.utils import list_to_str


class ExcitingStructure(ExcitingXMLInput):
    """ Class allowing exciting XML structure to be written from python data.
    """
    name = "structure"
    # Path type
    path_type = Union[str, Path]

    # Mandatory attribute "coord" taken out because it's specified inside the atoms
    _valid_atom_attributes = set(valid_attributes.atom_valid_attributes) - {"coord"}

    def __init__(self,
                 atoms,
                 lattice: Optional[list | np.ndarray] = None,
                 species_path: path_type = './',
                 crystal_properties: Optional[dict | ExcitingCrystalInput] = None,
                 species_properties: Optional[Dict[str, Union[dict, ExcitingSpeciesInput]]] = None,
                 **kwargs):
        """ Initialise instance of ExcitingStructure.

        TODO(Alex) Issue 117. Create our own class with a subset of methods common to ASE' Atom()
          Then we can have a single API for this init. If ASE is used, xAtom() is just a wrapper of
          Atom(), else we have some light methods.

        All valid attributes can be found in the module valid_attributes.py

        :param atoms: Atoms object of type ase.atoms.Atoms or of the form List[dict], for example:
         atoms = [{'species': 'X', 'position': [x, y, z]}, ...].
         Each dict can also optionally contain the _valid_atom_attributes:
         {'species': 'X', 'position': [x, y, z],
           'bfcmt': [bx, by, bz], 'lockxyz': [lx, ly, lz], 'mommtfix': [mx, my, mz]}.
        If atoms are defined with ASE, optional atomic_properties cannot be specified.
        Eventually, the list of atoms will be replaced with our custom class, which will extend ase.Atoms()
        with the additional, optional attributes.
        Species value can be a file_name without the suffix '.xml', which will be added automatically.

        :param lattice [a, b, c], where a, b and c are lattice vectors with 3 components.
         For example, a = [ax, ay, az]. Only required if one does not pass an ase Atoms object.
        :param species_path: Optional path to the location of species file/s.
        :param crystal_properties: Optional crystal properties.
        :param species_properties: Optional species properties, defined as:
        {'species1': {'rmt': rmt_value}, 'species2': {'rmt': rmt_value}}
        and with subtrees as:
        {'species1': {'rmt': rmt_value, 'LDAplusU': {'J': J, 'U': U, 'l': l}}, species2: ... }
        :param kwargs: Optional structure properties. Passed as kwargs.
        """
        if isinstance(species_path, Path):
            species_path = species_path.as_posix()
        super().__init__(speciespath=species_path, **kwargs)
        self.structure_attributes = {"speciespath": species_path, **kwargs}

        if isinstance(atoms, list) and lattice is None:
            raise ValueError("If atoms is a list, lattice must be passed as a separate argument.")

        # Simple container for atoms, as documented in __init__.
        if isinstance(atoms, list):
            check_lattice(lattice)
            check_lattice_vector_norms(lattice)
            self.lattice = np.asarray(lattice)
            self.species = [atom['species'].capitalize() for atom in atoms]
            self.positions = [atom['position'] for atom in atoms]
            self.atom_properties = self._init_atom_properties(atoms)
        else:
            self.lattice, self.species, self.positions = self._init_lattice_species_positions_from_ase_atoms(atoms)
            self.atom_properties = [{}] * len(self.species)

        self.unique_species = sorted(set(self.species))

        # Optional properties
        self.crystal_properties = self._initialise_subelement_attribute(ExcitingCrystalInput,
                                                                        crystal_properties or {})
        self.species_properties = dict(self._init_species_properties(species_properties))

    def __setattr__(self, name: str, value):
        """ Overload the attribute setting from the base class, since here we use different attribute names than
        defined in the schema.

        :param name: name of the attribute
        :param value: new value, can be anything
        """
        AbstractExcitingInput.__setattr__(self, name, value)

    def _init_lattice_species_positions_from_ase_atoms(self, atoms) -> tuple:
        """ Initialise lattice, species and positions from an ASE Atoms Object.

        Duck typing for atoms, such that ASE is not a hard dependency.

        :param atoms: ASE Atoms object.
        :return  Lattice, species, positions: Lattice, species and positions
        """
        try:
            cell = atoms.get_cell()
            # ASE works in Angstrom, whereas exciting expects atomic units
            lattice = np.asarray(cell) * angstrom_to_bohr
            species = [x.capitalize() for x in atoms.get_chemical_symbols()]
            if self.structure_attributes.get("cartesian"):
                positions = angstrom_to_bohr * atoms.get_positions()
            else:
                positions = atoms.get_scaled_positions()
            return lattice, species, positions
        except AttributeError:
            message = "atoms must either be an ase.atoms.Atoms object or List[dict], of the form" \
                      "[{'species': 'X', 'position': [x, y, z]}, ...]."
            raise AttributeError(message)

    def _init_atom_properties(self, atoms: List[dict]) -> List[dict]:
        """ Initialise atom_properties.

        For atoms that contain optional atomic properties, store them as
        dicts in a list of len(n_atoms). Atoms with none of these properties
        will be defined as empty dicts.

        For each element of atoms, one must have  {'species': 'X', 'position': [x, y, z]}  and
        may have the additional attributes: {'bfcmt': [bx, by, bz], 'lockxyz': [lx, ly, lz], 'mommtfix': [mx, my, mz]}.
        Extract the optional attributes and return in `atom_properties`, with string values.

        :param atoms: List container.
        :return atom_properties: List of atom properties. Each element is a dict.
        and the dict value has been converted to string - ready for XML usage.
        """
        atom_properties: List[dict] = []

        for atom in atoms:
            optional_property_keys = set(atom.keys()) & self._valid_atom_attributes
            optional_atom = {key: atom[key] for key in optional_property_keys}
            optional_properties = {}
            for key, value in optional_atom.items():
                optional_properties[key] = self._attributes_to_input_str[type(value)](value)
            atom_properties.append(optional_properties)

        return atom_properties

    def _init_species_properties(self, species_properties: dict) -> Dict[str, ExcitingXMLInput]:
        """ Initialise species_properties.

        For species without properties, return empty_properties: {'S': {}, 'Al': {}}.

        :param species_properties: Species properties
        :return Dict of ExitingXMLInput-species_properties.
        """
        if species_properties is None:
            species_properties = {}

        for species in self.unique_species:
            props = species_properties.get(species) or {}
            props["speciesfile"] = species + '.xml'
            yield species, self._initialise_subelement_attribute(ExcitingSpeciesInput, props)

    def get_lattice(self, convert_to_angstrom: bool = False) -> np.ndarray:
        """ Get the full lattice, meaning after the application of scale and stretch values to the stored
        lattice vectors.

        :param convert_to_angstrom: if True returns lattice in angstrom, else in bohr
        :return: full lattice vectors, stored row-wise in a matrix
        """
        lattice = copy.deepcopy(self.lattice)
        unit_conversion = bohr_to_angstrom if convert_to_angstrom else 1
        scale = getattr(self.crystal_properties, "scale", 1)
        lattice *= scale * unit_conversion

        stretch = np.array(getattr(self.crystal_properties, "stretch", [1, 1, 1]))
        lattice *= stretch[:, None]
        return lattice

    def _group_atoms_by_species(self) -> dict:
        """Get the atomic indices for atoms of each species.

        For example, for:
          species = ['Cd', 'S', 'Cd]
        return:
          indices = {'Cd': [1, 3], 'S' : [2]}

        :return dict indices: Indices of atoms in species and positions
        """
        indices = {}
        for x in self.unique_species:
            indices[x] = [i for i, element in enumerate(self.species) if element == x]
        return indices

    def _xml_atomic_subtree(self, x: str, species, atomic_indices):
        """ Add the required atomic positions and any optional attributes, per species.

        :param x: Species
        :param species: Empty SubElement for species x, which gets filled
        """
        for iatom in atomic_indices[x]:
            coord_str = list_to_str(self.positions[iatom])
            ElementTree.SubElement(species, "atom", coord=coord_str, **self.atom_properties[iatom]).text = ' '

    def to_xml(self) -> ElementTree.Element:
        """Convert structure attributes to XML ElementTree
        Makes use of the to_xml() function of the ExitingXMLInput class to convert values to string.

        Expect to return an XML structure which looks like:
          <structure speciespath="./">

           <crystal scale="1.00" scale="1.00" >
             <basevect>1.0 1.0 0.0</basevect>
             <basevect>1.0 0.0 1.0</basevect>
             <basevect>0.0 1.0 1.0</basevect>
           </crystal>

           <species speciesfile="Al.xml">
             <atom coord="0.0  0.0  0.0"> </atom>
           </species>

          </structure>

        :return ET structure: Element tree containing structure attributes.
        """
        structure_attributes = {
            key: self._attributes_to_input_str[type(value)](value) for key, value in self.structure_attributes.items()
        }
        structure = ElementTree.Element(self.name, **structure_attributes)
        structure.text = ' '

        # Lattice vectors
        crystal = self.crystal_properties.to_xml()
        structure.append(crystal)
        for vector in self.lattice:
            ElementTree.SubElement(crystal, "basevect").text = list_to_str(vector)

        # Species tags
        atomic_indices = self._group_atoms_by_species()
        for x in self.unique_species:
            species = self.species_properties[x].to_xml()
            structure.append(species)
            self._xml_atomic_subtree(x, species, atomic_indices)

        return structure
