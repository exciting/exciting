"""Structure class, mirroring that of exciting's structure XML sub-tree.
http://exciting.wikidot.com/ref:structure
"""
from typing import Optional, Union, List
import pathlib
import xml.etree.ElementTree as ET

from excitingtools.utils.utils import list_to_str
from excitingtools.utils.dict_utils import check_valid_keys
from excitingtools.structure.lattice import check_lattice, check_lattice_vector_norms
from excitingtools.input.base_class import ExcitingXMLInput


# Set of all elements
all_species = {'Ni', 'La', 'K', 'Xe', 'Ag', 'Bk', 'Co', 'Md', 'Lu', 'Ar',
               'Bi', 'Cm', 'H', 'Yb', 'Zn', 'Te', 'I', 'Cl', 'As', 'Mg',
               'No', 'Ta', 'N', 'Ac', 'Y', 'At', 'Tb', 'Tc', 'Au', 'O',
               'Lr', 'In', 'Ge', 'Re', 'Pm', 'Gd', 'Kr', 'Po', 'Sc', 'Rf',
               'Sb', 'Rb', 'Ru', 'Dy', 'Ho', 'Ra', 'Se', 'Sr', 'Fr', 'Ga',
               'Fe', 'Es', 'Si', 'Pr', 'Pd', 'Er', 'Rn', 'Ir', 'He', 'Eu',
               'Pt', 'Pu', 'Sn', 'Pb', 'Hf', 'Fm', 'Rh', 'Sm', 'Pa', 'Hg',
               'Os', 'B', 'U', 'Zr', 'Cf', 'C', 'Na', 'Li', 'Mo', 'Cs',
               'Al', 'V', 'Cd', 'Tm', 'Tl', 'Ba', 'Ce', 'W', 'Am', 'Cr',
               'Nb', 'Mn', 'S', 'Ca', 'Be', 'Br', 'Th', 'Ti', 'Np', 'Ne',
               'P', 'Cu', 'F', 'Nd'}


def check_muffin_tin_radii():
    """Ensure no two MT spheres are touching.

    Muffin tin radii cannot overlap. If MT radii have been explicitly specifed,
    check that none of the atom-centred MT spheres overlap (which will cause exciting to crash).

    TODO(Fab) Issue 117. Implement check that MT spheres do not overlap, and uncomment the method call above
      Construct distance matrix with scipy, using the unit cell
      Iterate through the d matrix and apply minimum image convention
      Find nearest neighbours (NN) for each atom - build a list of pairwise terms
      For each NN pair, add the MT radii along the bonding axis and evaluate.
    """
    raise NotImplementedError('Check of MT radii requires implementing. See exciting gitlab issue 117')


class ExcitingStructure(ExcitingXMLInput):
    """ Class allowing exciting XML structure to be written from python data.

    TODO(Fabian/Alex) 117. Implement all remaining attributes:
     All elements are species-specific. They should be passed like:
     species_properties = {'S': {'LDAplusU':{'J': J, 'U': U, 'l': l}} }
     Element: LDAplusU: J, U, l
     Element: dfthalfparam: ampl, cut, exponent
     Element: shell: ionization, number
    """
    # Path type
    path_type = Union[str, pathlib.Path]

    # Mandatory attributes not specified
    _valid_structure_attributes = {'autormt', 'cartesian', 'epslat', 'primcell', 'tshift'}
    _valid_crystal_attributes = {'scale', 'stretch'}
    _valid_species_attributes = {'rmt'}
    _valid_atom_attributes = {'bfcmt', 'lockxyz', 'mommtfix'}

    def __init__(self,
                 atoms,
                 lattice: Optional[list] = None,
                 species_path: Optional[path_type] = './',
                 structure_properties: Optional[dict] = None,
                 crystal_properties: Optional[dict] = None,
                 species_properties: Optional[dict] = None):
        """ Initialise instance of ExcitingStructure.

        TODO(Alex) Issue 117. Create our own class with a subset of methods common to ASE' Atom()
          Then we can have a single API for this init. If ASE is used, xAtom() is just a wrapper of
          Atom(), else we have some light methods.
        TODO(Alex/Fabian) Issue 117.
          structure_attributes and crystal_attributes could equally be kwargs.
          Consider changing or extending before the first major version.

        :param atoms: Atoms object of type ase.atoms.Atoms or of the form List[dict], for example:
         atoms = [{'species': 'X', 'position': [x, y, z]}, ...].
         Each dict can also optionally contain the _valid_atom_attributes:
         {'species': 'X', 'position': [x, y, z],
           'bfcmt': [bx, by, bz], 'lockxyz': [lx, ly, lz], 'mommtfix': [mx, my, mz]}.
        If atoms are defined with ASE, optional atomic_properties cannot be specified.
        Eventually, the list of atoms will be replaced with our custom class, which will extend ase.Atoms()
        with the additional, optional attributes.

        :param lattice [a, b, c], where a, b and c are lattice vectors with 3 components.
         For example, a = [ax, ay, az]. Only required if one does not pass an ase Atoms object.
        :param species_path: Optional path to the location of species file/s.
        :param structure_properties: Optional structure properties. See _valid_structure_attributes.
        :param crystal_properties: Optional crystal properties. See _valid_crystal_attributes
        :param species_properties: Optional species properties, defined as:
        {'species1': {'rmt': rmt_value}, 'species2': {'rmt': rmt_value}}
        """
        if isinstance(species_path, pathlib.Path):
            species_path = species_path.as_posix()

        if isinstance(atoms, list) and lattice is None:
            raise ValueError("If atoms is a list, lattice must be passed as a separate argument.")

        # Simple container for atoms, as documented in __init__.
        if isinstance(atoms, list):
            check_lattice(lattice)
            check_lattice_vector_norms(lattice)
            self.lattice = lattice
            self.species = [atom['species'].capitalize() for atom in atoms]
            self.positions = [atom['position'] for atom in atoms]
            self.atom_properties = self._init_atom_properties(atoms)
        else:
            self.lattice, self.species, self.positions = self._init_lattice_species_positions_from_ase_atoms(atoms)
            self.atom_properties = [{}] * len(self.species)

        # TODO(Fab) 117. Implement check that MT spheres do not overlap. check_muffin_tin_radii()
        self.species_path = species_path
        self.unique_species = sorted(set(self.species))

        # Catch symbols that are not valid elements
        check_valid_keys({x.lower() for x in self.unique_species},
                         {x.lower() for x in all_species}, name='Species input')

        # Optional properties
        self.structure_properties = self._init_structure_properties(structure_properties)
        self.crystal_properties = self._init_crystal_properties(crystal_properties)
        self.species_properties = self._init_species_properties(species_properties)

    @staticmethod
    def _init_lattice_species_positions_from_ase_atoms(atoms) -> tuple:
        """ Initialise lattice, species and positions from an ASE Atoms Object.

        Duck typing for atoms, such that ASE is not a hard dependency.

        :param atoms: ASE Atoms object.
        :return  Lattice, species, positions: Lattice, species and positions
        """
        try:
            cell = atoms.get_cell()
            # Convert to consistent form, [a, b, c], where a = [ax, ay, az]
            lattice = [list(cell[i, :]) for i in range(0, 3)]
            species = [x.capitalize() for x in atoms.get_chemical_symbols()]
            positions = atoms.get_positions()
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

    def _init_structure_properties(self, structure_properties: dict) -> dict:
        """ Initialise structure_properties.

        :param structure_properties: Dict of optional structure properties.
        :return Dict of structure_properties, with string values.
        """
        if structure_properties is None:
            return {}

        check_valid_keys(structure_properties.keys(), self._valid_structure_attributes, "structure_properties")

        return {key: str(value).lower() for key, value in structure_properties.items()}

    def _init_crystal_properties(self, crystal_properties: dict) -> dict:
        """ Initialise crystal_properties.

        :param crystal_properties: Dict of optional structure properties.
        :return Dict of crystal_properties, with string values.
        """
        if crystal_properties is None:
            return {}

        check_valid_keys(crystal_properties.keys(), self._valid_crystal_attributes, "crystal_properties")

        return {key: str(value) for key, value in crystal_properties.items()}

    def _init_species_properties(self, species_properties: dict) -> dict:
        """ Initialise species_properties.

        For species without properties, return empty_properties: {'S': {}, 'Al': {}}.

        :param species_properties: Species properties
        :return Dict of species_properties, with string values.
        """
        if species_properties is None:
            empty_properties = {s: {} for s in self.unique_species}
            return empty_properties

        new_species_properties = {}
        for species in self.unique_species:
            try:
                properties = species_properties[species]
                check_valid_keys(properties.keys(),
                                 self._valid_species_attributes,
                                 f"{species} element's species_properties")
                new_species_properties[species] = {key: str(value) for key, value in properties.items()}
            except KeyError:
                new_species_properties[species] = {}

        return new_species_properties

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
        :param species: Empty SubElement for species x
        :return species: species SubElement with all atomic positions included
        """
        for iatom in atomic_indices[x]:
            coord_str = list_to_str(self.positions[iatom])
            ET.SubElement(species, "atom", coord=coord_str, **self.atom_properties[iatom]).text = ' '
        return species

    def to_xml(self) -> ET.Element:
        """Convert structure attributes to XML ElementTree

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
        structure = ET.Element("structure", speciespath=self.species_path, **self.structure_properties)

        # Lattice vectors
        crystal = ET.SubElement(structure, "crystal", **self.crystal_properties)
        for vector in self.lattice:
            ET.SubElement(crystal, "basevect").text = list_to_str(vector)

        # Species tags
        atomic_indices = self._group_atoms_by_species()
        for x in self.unique_species:
            species = ET.SubElement(structure, "species", speciesfile=x + '.xml', **self.species_properties[x])
            species = self._xml_atomic_subtree(x, species, atomic_indices)

        return structure
