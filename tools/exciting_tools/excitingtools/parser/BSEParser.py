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
    out = {"variable1": list(data[:, 0]), "function1": list(data[:, 1]), "function2": list(data[:, 2]),
           "function3": list(data[:, 3])}

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
    out = {"variable1": list(data[:, 0]), "function1": list(data[:, 1]), "function2": list(data[:, 2])}

    return out
    

def parse_EXCITON_NAR_BSE(name):
    """
    Parser for EXCITON_NAR_BSE-singlet-TDA-BAR_SCR-full_OC.OUT
    """
    try:
        data = np.genfromtxt(name, skip_header=14)
    except:
        raise ParseError
    # TODO. Could probably be refactored to:
    # out = {str(i): data[:, i] for i in range(0, 6)}
    out = {}
    out["0"] = list(data[:,0])
    out["1"] = list(data[:,1])
    out["2"] = list(data[:,2])
    out["3"] = list(data[:,3])
    out["4"] = list(data[:,4])
    out["5"] = list(data[:,5])

    return out
