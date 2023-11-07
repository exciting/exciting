"""
Module containing Tol class and constants
"""
from collections import namedtuple
from typing import Optional

from excitingtools.constants.units import Unit, enum_to_string


class Tol:
    """
    Class for initialising a tolerance with a corresponding unit
    """
    def __init__(self, tol, unit: Optional[Unit] = Unit.null):
        self.tol = tol
        self.unit = unit

    def to_dict(self):
        """
        Convert object to dictionary.

        Method required as it performs enum to string conversion for units,
        which is required when writing to JSON
        """
        return {'tol': self.tol, 'unit': enum_to_string[self.unit]}


class TolWithMessage(Tol):
    """
    Class for initialising a tolerance with a corresponding unit and optional message
    """
    def __init__(self, tol, unit: Optional[Unit] = Unit.null, message: Optional[str] = ''):
        if isinstance(tol, Tol):
            super().__init__(tol.tol, tol.unit)
        elif isinstance(tol, (int, float, str)):
            super().__init__(tol, unit)
        else:
            raise ValueError('First argument of TolWithMessage must be an int, float, str, or instance of class Tol')
        self.message = message

    def to_dict(self):
        dictionary = super().to_dict()
        dictionary['message'] = self.message
        return dictionary


# Named tuple for default tolerances.
# Tolerances should be quantities (not units)
# Should be instantiated once per output file.
DefaultTolerances = namedtuple('DefaultTolerances',
                               ['integer',
                                'float',
                                'str',
                                'total_energy',
                                'energy',
                                'length',
                                'volume',
                                'inv_length',
                                'inv_volume',
                                'frequency',
                                'angle',
                                'temperature',
                                'oscillator_strength',
                                'loss_function',
                                'structure_factor',
                                'time',
                                'current_density',
                                'velocity',
                                'force',
                                'au',
                                'dos',
                                'effective_mass'
                                ])

# Set all defaults to None, such that tolerance templates are not required to define
# all namedtuple values.
DefaultTolerances.__new__.__defaults__ = (None, ) * len(DefaultTolerances._fields)


# Method to tolerance file name.
tol_file_name = {'groundstate': 'tolerance_ground_state.json',
                 'gw': 'tolerance_gw.json',
                 'hybrid': 'tolerance_hybrid.json',
                 'bse': 'tolerance_bse.json',
                 'xanes': 'tolerance_xanes.json',
                 'tddft': 'tolerance_tddft.json',
                 'rt_tddft': 'tolerance_rt_tddft.json',
                 'dos': 'tolerance_dos.json',
                 'band_structure': 'tolerance_bandstructure.json',
                 'plot': 'tolerance_plotting.json',
                 'wannier': 'tolerance_wannier.json',
                 'transport': 'tolerance_transport.json',
                 'optical_properties': 'tolerance_optical.json',
                 'electric_properties': 'tolerance_electric.json',
                 'core_properties': 'tolerance_core.json',
                 'spin_properties': 'tolerance_spin.json',
                 'wfplot': 'tolerance_wfplot.json'
                 # TODO(Bene #160): 
                 # Fix fastBSE test      
                 #'fastBSE': 'tolerance_fastBSE.json'
                 }

methods = list(tol_file_name.keys())

tol_file_to_method = {value: key for key, value in tol_file_name.items()}
