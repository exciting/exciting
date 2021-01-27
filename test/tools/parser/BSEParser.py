import xml.etree.cElementTree as ET
from xml.etree.ElementTree import ParseError
import re
import os
import sys
import numpy as np

from ErrornousFileError import ErrornousFileError

#Parser for EPSILON_NAR_BSE-singlet-TDA-BAR_SCR-full_OC  .OUT.xml, EPSILON_NAR_FXCMB1_OC  _QMT001.OUT.xml, EPSILON_NAR_NLF_FXCMB1_OC  _QMT001.OUT.xml, LOSS_NAR_FXCMB1_OC  _QMT001.OUT.xml
def parse_EPSILON_NAR(name):
    try:
        data = np.genfromtxt(name, skip_header=14)
    except:
        raise ParseError
    out = {}
    out["variabe1"]=list(data[:,0])
    out["function1"]=list(data[:,1])
    out["function2"]=list(data[:,2])
    out["function3"]=list(data[:,3])
    
    return out

#Parser for LOSS_NAR_FXCMB1_OC  _QMT001.OUT.xml, LOSS_NAR_NLF_FXCMB1_OC  _QMT001.OUT.xml                          
def parse_LOSS_NAR(name):
    try:
        data = np.genfromtxt(name, skip_header=14)
    except:
        raise ParseError
    out = {}
    out["variable1"]=list(data[:,0])
    out["function1"]=list(data[:,1])
    out["function2"]=list(data[:,2])
    
    return out
    
#Parser for EXCITON_NAR_BSE-singlet-TDA-BAR_SCR-full_OC  .OUT
def parse_EXCITON_NAR_BSE(name):
    try:
        data = np.genfromtxt(name, skip_header=14)
    except:
        raise ParseError
    out = {}
    out["0"] = list(data[:,0])
    out["1"] = list(data[:,1])
    out["2"] = list(data[:,2])
    out["3"] = list(data[:,3])
    out["4"] = list(data[:,4])
    out["5"] = list(data[:,5])

    return out

class excitingEPSILON_NAR():
    
    def __init__(self, filePath):
        self.path  = filePath                         #path of the EPSILON_NAR_...                                   
        if not os.path.exists(self.path):          #check if the path exists                                                                                         
            raise OSError('Path not valid.')
        try:
            self.data  = parse_EPSILON_NAR(self.path)         #data in form of a python dictionary                                                                       
        except ParseError:                        #checks if the parser works                                                                                        
            raise ErrornousFileError

    def returnData(self, dataPath):
        '''                                                                                                                                                          
        returns the data, saved at the dataPath.                                                                                                                     
        In:                                                                                                                                                          
        dataPath     string     path of the data                                                                                                                 
        '''
        keys = dataPath.split('/') #datapath splitted in to single keys                                                                                              
        l = len(keys)
        data = parse_EPSILON_NAR(self.path)
        for i in range (0,l):
            try:
                data = data[keys[i]] #value of the given datapath                                                                                                     
            except KeyError:      #checks if the datapath exists                                                                                                      
                print ("invalid key: ")
                print (keys[i])
                sys.exit()
        return data

class excitingEXCITON_NAR_BSE():

    def __init__(self, filePath):
        self.path  = filePath                         #path of the EXCITON_NAR_BSE-singlet-TDA-BAR_SCR-full_OC  .OUT                                   
        if not os.path.exists(self.path):          #check if the path exists                                                                                         
            raise OSError('Path not valid.')
        try:
            self.data  = parse_EXCITON_NAR_BSE(self.path)         #data in form of a python dictionary                                                     
        except ParseError:                        #checks if the parser works                                                                                        
            raise ErrornousFileError

    def returnData(self, dataPath):
        '''                                                                                                                                                          
        returns the data, saved at the dataPath.                                                                                                                     
        In:                                                                                                                                                          
        dataPath     string     path of the data                                                                                                                 
        '''
        keys = dataPath.split('/') #datapath splitted in to single keys                                                                                              
        l = len(keys)
        data = parse_EXCITON_NAR_BSE(self.path)
        for i in range (0,l):
            try:
                data = data[keys[i]] #value of the given datapath                                                                                                     
            except KeyError:      #checks if the datapath exists                                                                                                      
                print ("invalid key: ")
                print (keys[i])
                sys.exit()
        return data

class excitingLOSS_NAR():

    def __init__(self, filePath):
        self.path  = filePath                         #path of the LOSS_NAR_...                                                                                     
        if not os.path.exists(self.path):          #check if the path exists                                                                                           
            raise OSError('Path not valid.')
        try:
            self.data  = parse_LOSS_NAR(self.path)         #data in form of a python dictionary                                                                    \
                                                                                                                                                                       
        except ParseError:                        #checks if the parser works                                                                                          
            raise ErrornousFileError

    def returnData(self, dataPath):
        '''                                                                                                                                                            
        returns the data, saved at the dataPath.                                                                                                                       
        In:                                                                                                                                                            
        dataPath     string     path of the data                                                                                                                       
        '''
        keys = dataPath.split('/') #datapath splitted in to single keys                                                                                                
        l = len(keys)
        data = parse_LOSS_NAR(self.path)
        for i in range (0,l):
            try:
                data = data[keys[i]] #value of the given datapath                                                                                                     
            except KeyError:      #checks if the datapath exists                                                                                                      
                print ("invalid key: ")
                print (keys[i])
                sys.exit()
        return data
