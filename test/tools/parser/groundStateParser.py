import os
import sys
import xml.etree.ElementTree as ET
from xml.etree.ElementTree import ParseError
import re
import numpy as np

from ErrornousFileError import ErrornousFileError

def parse_INFO(name):
    ''' 
    Parser exciting INFO.OUT into a python dictionary.
    In:
        name     string     path of the file to parse
    Out: 
        info     dict       contains the content of the file to parse
    '''
    file = open(name)
    l = len(open(name).readlines())
    lines = []
    nscl = []
    nini = []

    # parses the complete file into a list
    for i in range(l):
        line = next(file)
        # stores the number of the first and last line of every iteration into a list
        if ('SCF iteration number' in line) or ('Hybrids iteration number' in line):
            nscl.append(i)
        if ('Convergency criteria checked for the last 2 iterations' in line) or ('Self-consistent loop stopped' in line):
            nscl.append(i)
        # stores the number of the first and last line of the initialization into a list
        if 'Starting initialization' in line:
            nini.append(i)
        if 'Ending initialization' in line:
            nini.append(i)
        lines.append(line)
    
    errorousFile = True
    for line in lines:
        if ("EXCITING" in line) and ("stopped" in line):
            errorousFile = False
    
    if errorousFile:
        raise ErrornousFileError()

    INFO = {}

    INFO['initialization'] = {}
    ini = []
    inits = {}
    k = 0

    # loops through all lines of the initialization
    for i in range (nini[0], nini[1]):
        # stores the lines, which have the format "variable : value" into a list
        if (':' in lines[i]) and ('+' not in lines[i]):
            lines[i] = lines[i].split(':')
            ini.append(lines[i])
            ini[k][0] = ini[k][0].strip()
            ini[k][1] = ini[k][1].strip()
            if ('Lattice vectors' in ini[k][0]) or ('Reciprocal lattice vectors' in ini[k][0]):
                ini[k][1] = []
                lines[i+1]=(lines[i+1].split())
                lines[i+2]=(lines[i+2].split())
                lines[i+3]=(lines[i+3].split())
                for j in range(3):
                    ini[k][1].append(lines[i+1][j])
                    ini[k][1].append(lines[i+2][j])
                    ini[k][1].append(lines[i+3][j])
                if ' ' in ini[k][1]:
                    ini[k][1] = ini[k][1].split()
            # stores variable-value pairs in a dictionary 
            inits.update({ini[k][0]:ini[k][1]})
            k = k+1
        # type of mixing is stored in the dictionary too
        if 'mixing' in lines[i]:
            lines[i] = lines[i].strip()
            inits.update({'mixing':lines[i]})
            
    INFO['initialization'] = inits
    
    scl1 = []
    INFO['scl'] = {}

    # loops through all scl's
    for j in range(len(nscl)-1):
        scls = {}
        scl = []
        k = 0
        # loops through all lines of the scl
        for i in range (nscl[j],nscl[j+1]):
            # stores the lines, which have the format "variable : value" into a list
            if (':' in lines[i]) and ('+' not in lines[i]) and ('(target)' not in lines[i]):
                lines[i] = lines[i].split(':')
                scl.append(lines[i])
                scl[k][0] = scl[k][0].strip()
                scl[k][1] = scl[k][1].strip()
                if ' ' in scl[k][1]:
                    scl[k][1] = scl[k][1].split()
                # stores variable-value pairs in a dictionary
                scls.update({scl[k][0]:scl[k][1]})
                k = k+1
        INFO['scl'][str(j+1)] = scls
    file.close()
    return(INFO)

def parse_info(name):
    '''
    Parser exciting info.xml into a python dictionary.
    In:
        name     string     path of the file to parse
    Out:
        info     dict       contains the content of the file to parse
    '''
    try:
        root = ET.parse(name)
    except AttributeError:
        raise ErrornousFileError
    i=0

    info = root.find('groundstate').attrib

    excitingRun=[]
    i = 0
    for node in root.find('groundstate').find('scl').iter('iter'):
        excitingRun.append(node.attrib)
        excitingRun[i]['energies']=node.find('energies').attrib
        excitingRun[i]['charges']=node.find('charges').attrib
        atom_nr=0
        atomic_charge=[]
        species=[]
        for atoms in node.find('charges').iter('atom'):
            if atom_nr==0 : species_old=atoms.get('species')
            atom_nr=atom_nr+1
            if atoms.get('species') == species_old:
                species.append({'muffin-tin':atoms.get('muffin-tin')})
            else:
                species_old=atoms.get('species')
                atomic_charge.append(species)
                species=[{'muffin-tin':atoms.get('muffin-tin')}]
                atomic_charge.append(species)
                excitingRun[i]['charges']['atomic']=atomic_charge
                excitingRun[i]['timing']=node.find('timing').attrib
        if node.find('moments') is not None:
            moments={}
            moments['momtot']=node.find('moments').find('momtot').attrib
            moments['interstitial']=node.find('moments').find('momtot').attrib
            moments['mommttot']=node.find('moments').find('interstitial').attrib
            excitingRun[i]['moments']=moments
            atom_nr=0
            atomic_moment=[]
            species=[]
            for atoms in node.find('moments').iter('atom'):
                if atom_nr==0 : species_old=atoms.get('species')
                atom_nr=atom_nr+1
                if atoms.get('species') == species_old:
                    species.append(atoms.find('mommt').attrib)
                else:
                    species_old=atoms.get('species')
                    atomic_moment.append(species)
                    species=[atoms.find('mommt').attrib]
                    atomic_moment.append(species)
                    excitingRun[i]['moments']['atomic']=atomic_moment
        i=i+1                
    info['scl'] = {}
    for item in excitingRun:        #converts list of scl-iterations into a dictionary  
        name = item['iteration']
        info['scl'][name] = item
        
    return info


def parse_atoms(name):
    '''                                                                                                  
    Parser exciting atoms.xml into a python dictionary.                                                   
    In:                                                                                                  
        name     string     path of the file to parse                                                    
    Out:                                                                                                 
        info     dict       contains the content of the file to parse                                    
    '''

    root = ET.parse(name)
    atoms = {}
    atoms['Hamiltonian']=root.find('Hamiltonian').attrib
    atom = []
    i = 0
    for node in root.findall('atom'):
        atom.append(node.attrib)
        
        spectrum = []
        states = node.find('spectrum')
        for state in states.findall('state'):
            spectrum.append(state.attrib)
            
        atom[i]['NumericalSetup']=node.find('NumericalSetup').attrib
        atom[i]['spectrum'] = {}
        j = 0
        for item in spectrum: #converts list of states into a dictionary
            name = str(j)
            atom[i]['spectrum'][name] = item
            j = j+1
        i = i+1
    atoms['atom'] = {}
    for item in atom:   #converts list of atoms into a dictionary
        name = item['chemicalSymbol']
        atoms['atom'][name] = item

    return atoms


def parse_eigval(name):
    '''                                                                                                  
    Parser exciting eigval.xml into a python dictionary.                                                   
    In:                                                                                                  
        name     string     path of the file to parse                                                    
    Out:                                                                                                 
        info     dict       contains the content of the file to parse                                    
    '''

    root = ET.parse(name).getroot()
    eigval = root.attrib

    kpts = []
    for node in root.findall('kpt'):
        kpt = node.attrib
        state = []
        for subnode in node:
            state.append(subnode.attrib)
            kpt['state']= {}   #converts list of states into a dictionary
        for item in state:
            name = item['ist']
            kpt['state'][name] = item
            kpts.append(kpt)
            eigval['kpt'] = {}
    for item in kpts:   #converts list of kpts into a dictionary
        name = item['ik']
        eigval['kpt'][name] = item
        
    return eigval

def parse_evalcore(name):
    '''                                                                                                  
    Parser exciting evalcore.xml into a python dictionary.                                                   
    In:                                                                                                  
        name     string     path of the file to parse                                                    
    Out:                                                                                                 
        info     dict       contains the content of the file to parse                                    
    '''

    root = ET.parse(name).getroot()
    evalcore = root.attrib

    speciess = []
    for node in root.findall('species'):
        species = node.attrib
        atoms = []
        for subnode in node:
            atom = subnode.attrib
            states = []
            for subnode1 in subnode:
                state = subnode1.attrib
                states.append(state)
            atom['state']= {}
            for item in states:    #converts list of states into a dictionary
                name = item['ist']
                atom['state'][name] = item
            atoms.append(atom)
            species['atom'] = {}
            for item in atoms:    #converts list of atoms into a dictionary
                name = item['ia']
                species['atom'][name] = item
        speciess.append(species)
    evalcore['species'] = {}
    for item in speciess:   #converts list of species into a dictionary
        name = item['chemicalSymbol']
        evalcore['species'][name] = item
        
    return evalcore

def parse_geometry(name):
    '''                                                                                                  
    Parser exciting geometry.xml into a python dictionary.                                                   
    In:                                                                                                  
        name     string     path of the file to parse                                                    
    Out:                                                                                                 
        info     dict       contains the content of the file to parse                                    
    '''

    root = ET.parse(name).getroot()
    structure = root.find('structure').attrib
    crystal = root.find('structure').find('crystal').attrib
    geometry = {}
    geometry['structure'] = structure
    structure['crystal'] = crystal
    speciess = []
    for node in root.find('structure').findall('species'):
        species = node.attrib
        atoms = []
        for subnode in node:
            atom = subnode.attrib
            atoms.append(atom)
            species['atom'] = {}
        i = 1
        for item in atoms:
            name = str(i)
            species['atom'][name] = item['coord'].split()
            i = i+1
        speciess.append(species)
    structure['species'] = {}
    j = 1
    for item in speciess:
        name = str(j)
        structure['species'][name] = item
        j = j+1

    basevects = []
    for node in root.find('structure').find('crystal').findall('basevect'):
        basevect = node.text
        basevects.append(basevect)
        structure['crystal']['basevect'] = {}
        k = 1
    for item in basevects:
        name = str(k)
        structure['crystal']['basevect'][name] = item
        k = k+1
    return(geometry)

    
    

class excitingINFO():

    def __init__(self, filePath):
        self.path  = filePath                         #path of the INFO.OUT
        if not os.path.exists(self.path):          #check if the path exists
            raise OSError('Path not valid.')
        
        try:
            self.data  = parse_INFO(self.path)         #data in form of a python dictionary 
        except (ParseError, AttributeError):
            raise ErrornousFileError                   #checks if the parser works 
            
    def returnData(self, dataPath):
        ''' 
        returns the data, saved at the dataPath.
        In: 
            dataPath     string     path of the data
        '''
        keys = dataPath.split('/') #datapath splitted in to single keys
        l = len(keys)
        data = parse_INFO(self.path)
        for i in range (0,l):
            try:
                data = data[keys[i]] #value of the given datapath
            except KeyError:      #checks if the datapath exists
                print ("invalid key: ")
                print (keys[i])
                sys.exit()
        return data


class excitingInfo():
    
    def __init__(self, filePath):
        self.path  = filePath                         #path of the info.xml
        if not os.path.exists(self.path):          #check if the path exists
            raise OSError('Path not valid.')

        try:
            self.data  = parse_info(self.path)         #data in form of a python dictionary
        except (ParseError, AttributeError):                        #checks if the parser works
            raise ErrornousFileError

    def returnData(self, dataPath):
        '''
        returns the data, saved at the dataPath.
        In:
            dataPath     string     path of the data
        '''
        keys = dataPath.split('/') #datapath splitted in to single keys
        l = len(keys)
        data = parse_info(self.path)
        for i in range (0,l):
            try:
                data = data[keys[i]] #value of the given datapath
            except KeyError:      #checks if the datapath exists 
                print ("invalid key: ")
                print (keys[i])
                sys.exit()
        return data

    
class excitingAtoms():

    def __init__(self, filePath):
        self.path  = filePath                         #path of the atoms.xml                            
        if not os.path.exists(self.path):          #check if the path exists                           
            raise OSError('Path not valid.')

        try:
            self.data  = parse_atoms(self.path)         #data in form of a python dictionary            
        except (ParseError, AttributeError):                        #checks if the parser works                           
            raise ErrornousFileError

    def returnData(self, dataPath):
        '''                                                                                             
        returns the data, saved at the dataPath.                                                        
        In:                                                                                             
            dataPath     string     path of the data                                                    
        '''
        keys = dataPath.split('/') #datapath splitted in to single keys                                 
        l = len(keys)
        data = parse_atoms(self.path)
        for i in range (0,l):
           try:
               data = data[keys[i]] #value of the given datapath                                        
           except KeyError:      #checks if the datapath exists                                         
               print ("invalid key: ")
               print (keys[i])
               sys.exit()
        return data
    

class excitingEigval():

    def __init__(self, filePath):
        self.path  = filePath                         #path of the eigval.xml             
        if not os.path.exists(self.path):          #check if the path exists              
            raise OSError('Path not valid.')
        try:
            self.data  = parse_eigval(self.path)         #data in form of a python dictionary
        except (ParseError, AttributeError):                        #checks if the parser works             
            raise ErrornousFileError

    def returnData(self, dataPath):
            '''                                                                                             
            returns the data, saved at the dataPath.                                                        
            In:                                                                                             
            dataPath     string     path of the data                                                 
            '''
            keys = dataPath.split('/') #datapath splitted in to single keys                   
            l = len(keys)
            data = parse_eigval(self.path)
            for i in range (0,l):
                try:
                    data = data[keys[i]] #value of the given datapath                 
                except KeyError:      #checks if the datapath exists                            
                    print ("invalid key: ")
                    print (keys[i])
                    sys.exit()
            return data


class excitingEvalcore():

    def __init__(self, filePath):
        self.path  = filePath                         #path of the evalcore.xml                         
        if not os.path.exists(self.path):          #check if the path exists                            
            raise OSError('Path not valid.')
        try:
            self.data  = parse_evalcore(self.path)         #data in form of a python dictionary         
        except (ParseError, AttributeError):                        #checks if the parser works                           
            raise ErrornousFileError

    def returnData(self, dataPath):
        '''                                                                                             
        returns the data, saved at the dataPath.                                                        
        In:                                                                                             
            dataPath     string     path of the data                                                    
        '''
        keys = dataPath.split('/') #datapath splitted in to single keys                                 
        l = len(keys)
        data = parse_evalcore(self.path)
        for i in range (0,l):
           try:
               data = data[keys[i]] #value of the given datapath                                        
           except KeyError:      #checks if the datapath exists
               print ("invalid key: ")
               print (keys[i])
               sys.exit()
        return data


class excitingGeometry():

    def __init__(self, filePath):
        self.path  = filePath                         #path of the geometry.xml                        
        if not os.path.exists(self.path):          #check if the path exists                           
            raise OSError('Path not valid.')
        try:
            self.data  = parse_geometry(self.path)         #data in form of a python dictionary         
        except (ParseError, AttributeError):                        #checks if the parser works                           
            raise ErrornousFileError

    def returnData(self, dataPath):
        '''                                                                                             
        returns the data, saved at the dataPath.                                                        
        In:                                                                                             
            dataPath     string     path of the data                                                   
        '''
        keys = dataPath.split('/') #datapath splitted in to single keys                                
        l = len(keys)
        data = parse_geometry(self.path)
        for i in range (0,l):
           try:
               data = data[keys[i]] #value of the given datapath                                        
           except KeyError:      #checks if the datapath exists                                         
               print ("invalid key: ")
               print (keys[i])
               sys.exit()
        return data
