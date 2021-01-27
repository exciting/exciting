import xml.etree.cElementTree as ET
from xml.etree.ElementTree import ParseError
import re
import numpy as np
import os
import sys

from ErrornousFileError import ErrornousFileError

#Parser for EFERMI_GW.OUT
def parse_efermi_gw(name):
    data = {}
    try:
        data["EFERMI_GW"] = np.genfromtxt(name)
    except:
        raise ParseError
    return data

#Parser for EVALQP.DAT
def parse_evalqp(name):
    try:
        data = np.genfromtxt(name, skip_header=2)
    except:
        raise ParseError
    out = {}
    out["0"] = list(data[:,0])
    out["1"] = list(data[:,1])
    out["2"] = list(data[:,2])
    out["3"] = list(data[:,3])
    out["4"] = list(data[:,4])
    out["5"] = list(data[:,5])
    out["6"] = list(data[:,6])
    out["7"] = list(data[:,7])
    out["8"] = list(data[:,8])
    out["9"] = list(data[:,9])
    out["10"] = list(data[:,10])
    
    return out

#Parser for VXCNN.DAT
def parse_vxcnn(name):
    try:
        data = np.genfromtxt(name, skip_header=1)
    except:
        raise ParseError
    out = {}
    out["0"] = list(data[:,0])
    out["1"] = list(data[:,1])
    out["2"] = list(data[:,2])

    return out

#Parser for EPS00_GW.OUT
def parse_eps00_gw(name):
    try:
        data = np.loadtxt(name, skiprows=3, comments=["r", "f"])
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

class excitingEFERMI_GW():

    def __init__(self, filePath):
        self.path  = filePath                         #path of the EFERMI_GW.OUT                                       
        if not os.path.exists(self.path):          #check if the path exists                                                                          
            raise OSError('Path not valid.')
        try:
            self.data  = parse_efermi_gw(self.path)         #data in form of a python dictionary                                                                     
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
        data = parse_efermi_gw(self.path)
        for i in range (0,l):
            try:
                data = data[keys[i]] #value of the given datapath                                                                                                          
            except KeyError:      #checks if the datapath exists                                                                                                           
                print ("invalid key: ")
                print (keys[i])
                sys.exit()
        return data

class excitingEVALQP():

    def __init__(self, filePath):
        self.path  = filePath                         #path of the EVALQP.DAT         
        if not os.path.exists(self.path):          #check if the path exists                                                                                               
            raise OSError('Path not valid.')
        try:
            self.data  = parse_evalqp(self.path)         #data in form of a python dictionary                                                                     
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
        data = parse_evalqp(self.path)
        for i in range (0,l):
            try:
                data = data[keys[i]] #value of the given datapath                                                                                                          
            except KeyError:      #checks if the datapath exists                                                                                                           
                print ("invalid key: ")
                print (keys[i])
                sys.exit()
        return data

class excitingVXCNN():

    def __init__(self, filePath):
        self.path  = filePath                         #path of the VXCNN.DAT
        if not os.path.exists(self.path):          #check if the path exists                                                                                               
            raise OSError('Path not valid.')
        try:
            self.data  = parse_vxcnn(self.path)         #data in form of a python dictionary                                                                     
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
        data = parse_vxcnn(self.path)
        for i in range (0,l):
            try:
                data = data[keys[i]] #value of the given datapath                                                                                                          
            except KeyError:      #checks if the datapath exists                                                                                                           
                print ("invalid key: ")
                print (keys[i])
                sys.exit()
        return data

class excitingEPS00_GW():

    def __init__(self, filePath):
        self.path  = filePath                         #path of the EPS00_GW.OUT                                                   
        if not os.path.exists(self.path):          #check if the path exists                                                                                           
            raise OSError('Path not valid.')
        try:
            self.data  = parse_eps00_gw(self.path)         #data in form of a python dictionary                                                                 
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
        data = parse_eps00_gw(self.path)
        for i in range (0,l):
            try:
                data = data[keys[i]] #value of the given datapath                                                                                                     
            except KeyError:      #checks if the datapath exists                                                                                                      
                print ("invalid key: ")
                print (keys[i])
                sys.exit()
        return data
