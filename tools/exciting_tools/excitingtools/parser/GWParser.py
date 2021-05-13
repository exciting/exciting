"""
Module containing parsers for GW output files
"""
from xml.etree.ElementTree import ParseError
import numpy as np


def parse_efermi_gw(name: str) -> dict:
    """
    Parser for EFERMI_GW.OUT
    """
    data = {}
    try:
        data["EFERMI_GW"] = np.genfromtxt(name)
    except:
        raise ParseError
    return data


def parse_evalqp(name: str) -> dict:
    """
    Parser for EVALQP.DAT
    """
    try:
        data = np.genfromtxt(name, skip_header=2)
    except:
        raise ParseError
    out = {"0": list(data[:, 0]), "1": list(data[:, 1]), "2": list(data[:, 2]), "3": list(data[:, 3]),
           "4": list(data[:, 4]), "5": list(data[:, 5]), "6": list(data[:, 6]), "7": list(data[:, 7]),
           "8": list(data[:, 8]), "9": list(data[:, 9]), "10": list(data[:, 10])}

    return out


def parse_vxcnn(name: str) -> dict:
    """
    Parser for VXCNN.DAT
    """
    try:
        data = np.genfromtxt(name, skip_header=1)
    except:
        raise ParseError
    out = {"0": list(data[:, 0]), "1": list(data[:, 1]), "2": list(data[:, 2])}

    return out


def parse_eps00_gw(name: str) -> dict:
    """
    Parser for EPS00_GW.OUT
    """
    try:
        data = np.loadtxt(name, skiprows=3, comments=["r", "f"])
    except:
        raise ParseError
    out = {"0": list(data[:, 0]), "1": list(data[:, 1]), "2": list(data[:, 2]), "3": list(data[:, 3]),
           "4": list(data[:, 4]), "5": list(data[:, 5])}

    return out
