"""
Parsers for real-time TDDFT output files
"""
from xml.etree.ElementTree import ParseError
import numpy as np


def parse_nexc(name, skiprows=1):
    """
    Parser for NEXC.OUT
    """
    try:
        data = np.genfromtxt(name, skip_header=skiprows)
    except:
        raise ParseError
    out = {"Time": list(data[:, 0]), "number_electrons_GroundState": list(data[:, 1]),
           "number_electrons_ExcitedState": list(data[:, 2]), "sum": list(data[:, 3])}

    return out


def parse_jind(name, skiprows=0):
    """
    Parser for JIND.OUT
    """
    try:
        data = np.genfromtxt(name, skip_header=skiprows)
    except:
        raise ParseError
    out = {"time": list(data[:, 0]), "Jx": list(data[:, 1]), "Jy": list(data[:, 2]), "Jz": list(data[:, 3])}

    return out


def parse_etot(name):
    """
    Parser for ETOT_RTTDDFT.OUT
    """
    try:
        data = np.genfromtxt(name, skip_header=1)
    except:
        raise ParseError
    out = {"Time": list(data[:, 0]), "ETOT": list(data[:, 1]), "Madelung": list(data[:, 2]),
           "Eigenvalues-Core": list(data[:, 3]), "Eigenvalues-Valence": list(data[:, 4]), "Exchange": list(data[:, 5]),
           "Correlation": list(data[:, 6]), "XC-potential": list(data[:, 7]), "Coulomb pot. energy": list(data[:, 8])}

    return out


def parse_eigval_screenshots(name:str):
    """
    Parser for EIGVAL_*.OUT
    """
    try:
        data = np.genfromtxt(name, skip_header=1, comments='ik')
    except:
        raise ParseError
    return {"0": list(data[:,1])}


def parse_proj_screenshots(name:str):
    """
    Parser for PROJ_*.OUT

    Returned values are flattened lists
    """
    try:
        data = np.genfromtxt(name, skip_header=1, comments='ik')
    except:
        raise ParseError
    return {"0": list(data[:,:].flatten())}
