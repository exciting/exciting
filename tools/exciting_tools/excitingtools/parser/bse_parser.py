"""
Module containing parsers for BSE output files
"""
from xml.etree.ElementTree import ParseError
import numpy as np


def parse_EPSILON_NAR(name):
    """
    Parser for:
        EPSILON_NAR_BSE-singlet-TDA-BAR_SCR-full_OC.OUT.xml,
        EPSILON_NAR_FXCMB1_OC_QMT001.OUT.xml,
        EPSILON_NAR_NLF_FXCMB1_OC_QMT001.OUT.xml,
        LOSS_NAR_FXCMB1_OC_QMT001.OUT.xml
    """
    try:
        data = np.genfromtxt(name, skip_header=14)
    except:
        raise ParseError
    out = {
        "frequency": data[:, 0],
        "real_oscillator_strength": data[:, 1],
        "imag_oscillator_strength": data[:, 2],
        "real_oscillator_strength_kkt": data[:, 3]
    }

    return out


def parse_LOSS_NAR(name):
    """
    Parser for:
     LOSS_NAR_FXCMB1_OC_QMT001.OUT.xml,
     LOSS_NAR_NLF_FXCMB1_OC_QMT001.OUT.xml
    """
    try:
        data = np.genfromtxt(name, skip_header=14)
    except:
        raise ParseError
    out = {
        "frequency": data[:, 0],
        "real_oscillator_strength": data[:, 1],
        "imag_oscillator_strength": data[:, 2]
    }

    return out


def parse_EXCITON_NAR_BSE(name):
    """
    Parser for EXCITON_NAR_BSE-singlet-TDA-BAR_SCR-full_OC.OUT
    """
    try:
        data = np.genfromtxt(name, skip_header=14)
    except:
        raise ParseError
    out = {}
    out["state"] = data[:, 0]
    out["energy"] = data[:, 1]
    out["energy_shifted"] = data[:, 2]
    out["abs_oscillator_strength"] = data[:, 3]
    out["real_oscillator_strength"] = data[:, 4]
    out["imaginary_oscillator_strength"] = data[:, 5]

    return out
