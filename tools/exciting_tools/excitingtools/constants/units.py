""" Units.

Physical constants, defined according to [CODATA 2018])(http://physics.nist.gov/constants)
One could also consider importing them from scipy
"""
import enum

bohr_to_angstrom = 0.529177210903
angstrom_to_bohr = 1. / bohr_to_angstrom
Hartree_to_eV = 27.211396641308


class Unit(enum.Enum):
    """
    Enum class for exciting units. All names are defined with
    as lowercase, for consistency.

    This could be replaced with [PINT](https://pint.readthedocs.io/en/stable/), as used by NOMAD, however it's
    currently not required. If/when we wish to do unit manipulation, it should be reconsidered.
    """
    hartree = enum.auto()
    inv_hartree = enum.auto()
    ev = enum.auto()
    inv_ev = enum.auto()
    kelvin = enum.auto()
    bohr = enum.auto()
    bohr_pow_3 = enum.auto()
    inv_bohr = enum.auto()
    inv_bohr_pow_3 = enum.auto()
    au = enum.auto()
    degrees = enum.auto()
    GK_max = enum.auto()
    electron_rest_mass = enum.auto()
    bohr_velocity_over_bohr_radius = enum.auto()
    bohr_velocity = enum.auto()
    force = enum.auto()
    null = enum.auto()


# Map Unit enums to strings.
# Required because JSON cannot dump objects to file.
enum_to_string = {
    Unit.hartree: 'Hartree',
    Unit.inv_hartree: '1/Hartree',
    Unit.ev: 'eV',
    Unit.inv_ev: 'eV^-1',
    Unit.kelvin: 'K',
    Unit.bohr: 'Bohr',
    Unit.bohr_pow_3: 'Bohr^3',
    Unit.inv_bohr: 'Bohr^-1',
    Unit.inv_bohr_pow_3: 'Bohr^-3',
    Unit.au: 'a.u.',
    Unit.degrees: 'degrees',
    Unit.GK_max: 'GK_max',
    Unit.electron_rest_mass: 'm_electron',
    Unit.bohr_velocity_over_bohr_radius: 'v_Bohr/r_Bohr',
    Unit.bohr_velocity: 'Bohr/t_Bohr',
    Unit.force: 'Hartree/Bohr', 
    Unit.null: 'null'
}
