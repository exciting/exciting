# Units
from excitingtools.constants.units import Unit

# User-level objects
from excitingtools.dataclasses import BandData, EigenValues

# old deprecated API
# Parsers returning to dicts
# Questionable whether one should expose this - required for test framework recursive comparisons
# Typically not the API one wants to expose to the user, as parsed dict keys are subject to change
from excitingtools.exciting_dict_parsers.parser_factory import parse, parser_chooser

# Parsers returning to objects
from excitingtools.exciting_obj_parsers import parse_band_structure

# Input objects
from excitingtools.input import (
    ExcitingEPHInput,
    ExcitingGroundStateInput,
    ExcitingGWInput,
    ExcitingInputXML,
    ExcitingMDInput,
    ExcitingPhononsInput,
    ExcitingPropertiesInput,
    ExcitingRelaxInput,
    ExcitingStructure,
    ExcitingXSInput,
)

try:
    from importlib import metadata
except ImportError:
    import importlib_metadata as metadata

__version__ = metadata.version("excitingtools")

__all__ = [
    "Unit",
    "BandData",
    "EigenValues",
    "parse",
    "parser_chooser",
    "parse_band_structure",
    "ExcitingGroundStateInput",
    "ExcitingXSInput",
    "ExcitingPropertiesInput",
    "ExcitingRelaxInput",
    "ExcitingPhononsInput",
    "ExcitingGWInput",
    "ExcitingMDInput",
    "ExcitingEPHInput",
    "ExcitingInputXML",
    "ExcitingStructure",
]
