import os
import sys
import xml.etree.cElementTree as ET
from xml.etree.ElementTree import ParseError
import re
import numpy as np

from ErrornousFileError import ErrornousFileError


#Parser for NEXC.OUT
def parse_nexc(name,skiprows=1):
    try:
        data = np.genfromtxt(name, skip_header=skiprows)
    except:
        raise ParseError
    out = {}
    out["Time"] = list(data[:,0])
    out["number_electrons_GroundState"] = list(data[:,1])
    out["number_electrons_ExcitedState"] = list(data[:,2])
    out["sum"] = list(data[:,3])

    return out

#Parser for JIND.OUT
def parse_jind(name,skiprows=0):
    try:
        data = np.genfromtxt(name, skip_header=skiprows)
    except:
        raise ParseError
    out = {}
    out["time"] = list(data[:,0])
    out["Jx"] = list(data[:,1])
    out["Jy"] = list(data[:,2])
    out["Jz"] = list(data[:,3])

    return out

#Parser for ETOT_RTTDDFT.OUT
def parse_etot(name):
    try:
        data = np.genfromtxt(name, skip_header=1)
    except:
        raise ParseError
    out = {}
    out["Time"] = list(data[:,0])
    out["ETOT"] = list(data[:,1])
    out["Madelung"] = list(data[:,2])
    out["Eigenvalues-Core"] = list(data[:,3])
    out["Eigenvalues-Valence"] = list(data[:,4])
    out["Exchange"] = list(data[:,5])
    out["Correlation"] = list(data[:,6])
    out["XC-potential"] = list(data[:,7])
    out["Coulomb pot. energy"] = list(data[:,8])

    return out

#Parser for EIGVAL_*.OUT and PROJ_*.OUT
def parse_screenshots( name, convertTo1D = False ):
    try:
        data = np.genfromtxt(name, skip_header=1, comments='ik')
    except:
        raise ParseError

    out = {}
    if ( convertTo1D ):
        # convert numpy.array to a 1D array, and then to list
        out["0"] = list( data[:,:].flatten() )
    else:
        out["0"] = list( data[:,1] )
    return out

class excitingJIND():

    def __init__(self, filePath, skiprows=0):
        self.path  = filePath
        if not os.path.exists(self.path):          #check if the path exists
            raise OSError('Path not valid.')
        try:
            self.data  = parse_jind(self.path, skiprows) #data in form of a python dictionary
        except ParseError:                        #checks if the parser works
            raise ErrornousFileError

class excitingNEXC():

    def __init__(self, filePath, skiprows=1):
        self.path  = filePath
        if not os.path.exists(self.path):          #check if the path exists
            raise OSError('Path not valid.')
        try:
            self.data  = parse_nexc(self.path, skiprows) #data in form of a python dictionary
        except ParseError:                        #checks if the parser works
            raise ErrornousFileError


class excitingETOT():

    def __init__( self, filePath, skiprows=0 ):
        self.path  = filePath                         
        if not os.path.exists(self.path):          #check if the path exists
            raise OSError('Path not valid.')
        try:
            self.data  = parse_etot(self.path) #data in form of a python dictionary
        except ParseError:                        #checks if the parser works
            raise ErrornousFileError


class excitingScreenshots():

    def __init__( self, filePath, convertTo1D = False ):
        self.path  = filePath
        if not os.path.exists(self.path):          #check if the path exists
            raise OSError('Path not valid.')
        try:
            self.data  = parse_screenshots( self.path, convertTo1D ) #data in form of a python dictionary
        except ParseError:                        #checks if the parser works
            raise ErrornousFileError
