import os
import sys
import xml.etree.cElementTree as ET
from xml.etree.ElementTree import ParseError
import re
import numpy as np

#Parser for RHO3D.xml, VCL3D.xml, VXC3D.xml, WF3D.xml, ELF3D.xml, EF3D.xml
def parse_plot_3D(name):
    root = ET.parse(name)
    plot_3D = {}
    plot_3D["title"] = root.find("title").text
    grid = root.find("grid").attrib
    axis = []
    for ax in root.find("grid").findall("axis"):
        axis.append(ax.attrib)
    grid["axis"] = {}
    for item in axis:
        name = item["name"]
        grid["axis"][name] = item
    grid["value"] = root.find("grid").find("value").attrib
    plot_3D["grid"] = grid

    func_node = root.find("function")
    function = root.find("function").attrib
    row0 = []
    for node in func_node.findall("row"):
        row = node.attrib
        row1 = []
        for nod in node:
            row2 = nod.attrib
            row2["data"] = nod.text.split()
            row1.append(row2)
        row["row"] = {}
        for item in row1:
            name = item["index"]
            row["row"][name] = item
        row0.append(row)
    function["row"] = {}
    for item in row0:
        name = item["index"]
        function["row"][name] = item
    plot_3D["function"] = function

    return plot_3D

#Parser for LSJ.xml
def parse_LSJ(name):
    root = ET.parse(name).getroot()
    LSJ = {}
    species = []
    for node in root.findall("species"):
        spec = node.attrib
        atom = []
        for nod in node.findall("atom"):
            at = nod.attrib
            at["L"] = nod.find("L").text.split()
            at["S"] = nod.find("S").text.split()
            at["J"] = nod.find("J").text.split()
            atom.append(at)
        spec["atom"] = {}
        for item in atom:
            name = item["n"]
            spec["atom"][name] = item
        species.append(spec)
    LSJ["species"] = {}
    for item in species:
        name = item['n']
        LSJ['species'][name] = item

    return LSJ

#Parser for EFG.xml
def parse_EFG(name):
    root = ET.parse(name).getroot()
    EFG = {}
    species = []
    for node in root.findall("species"):
        spec = node.attrib
        atom = []
        for nod in node.findall("atom"):
            at = nod.attrib
            ef = []
            for no in nod.findall("EFG-tensor"):
                ef = no.attrib
                line = []
                for n in no.findall("line"):
                    li = n.text.split(" ")
                    line.append(li)
                ef["matrix"] = line
            at["efg"] = ef
            atom.append(at)
        spec["atom"] = {}
        for item in atom:
            name = item["n"]
            spec["atom"][name] = item
        species.append(spec)
    EFG["species"] = {}
    for item in species:
        name = item["n"]
        EFG["species"][name] = item

    return EFG

#Parser for mossbauer.xml
def parse_mossbauer(name):
    root = ET.parse(name).getroot()
    mossbauer = {}
    species = []
    for node in root.findall("species"):
        spec = node.attrib
        atom = []
        for nod in node:
            atom.append(nod.attrib)
        spec["atom"] = {}
        for item in atom:
            name = item["n"]
            spec["atom"][name] = item
        species.append(spec)
    mossbauer["species"] = {}
    for item in species:
        name = item["n"]
        mossbauer["species"][name] = item
    
    return mossbauer

#parser for expiqr.xml
def parse_expiqr(name):
    root = ET.parse(name).getroot()
    expiqr = {}
    expiqr["q-vector"] = root.find("q-vector").attrib
    kgrid = {}
    for k in root.find("k-grid"):
        kgrid = k.attrib
        states = []
        for st in k.findall("state"):
            state = st.attrib
            states.append(state)
            statesj = []
            for s in st.findall("state"):
                statej = s.attrib
                statesj.append(statej)
            state["state"] = {}
            for item in statesj:
                name = item["j"]
                state["state"][name] = item
        kgrid["state"] = {}
        for item in states:
            name = item["i"]
            kgrid["state"][name] = item
    expiqr["k-grid"] = kgrid
    return expiqr

#parser for effmass.xml
def parse_effmass(name):
    root = ET.parse(name).getroot()
    effmass = {}
    effmass["k-point"] = root.find("k-point").attrib
    state = []
    for node in root.findall("state"):
        st = node.attrib
        evdk = node.find("evdk").attrib
        matrix1 = []
        for line in node.find("evdk").findall("line"):
            matrix1.append(line.text.split())
        evdk["matrix"] = matrix1
        emt = node.find("emt").attrib
        matrix2 = []
        for line in node.find("emt").findall("line"):
            matrix2.append(line.text.split())
        emt["matrix"] = matrix2
        st["evdk"] = evdk
        st["emt"] = emt
        state.append(st)
    effmass["state"] = {}
    for item in state:
        name = item["n"]
        effmass["state"][name] = item

    return effmass

#parser for bandstructure.xml
def parse_bandstructure(name):
    root = ET.parse(name).getroot()
    bandstructure = {}
    bandstructure["title"] = root.find("title").text
    bands = []
    for node in root.findall("band"):
        band = node.attrib
        points = []
        for nod in node.findall("point"):
            point = nod.attrib
            points.append(point)
        bands.append(band)
        band["point"] = {}
        i = 1
        for item in points:
            name = str(i)
            band["point"][name] = item
            i = i+1
    bandstructure["band"] = {}
    j = 1
    for item in bands:
        name = str(j)
        bandstructure["band"][name] = item
        j = j+1
    return bandstructure

#parser for dos.xml
def parse_dos(name):
    root = ET.parse(name).getroot()
    dos = {}
    dos["title"] = root.find("title").text
    totaldos = root.find("totaldos").attrib
    dos["totaldos"] = totaldos
    diagram  = root.find("totaldos").find("diagram").attrib
    dos["totaldos"]["diagram"] = diagram
    points = []
    for node in root.find("totaldos").find("diagram").findall("point"):
        point = node.attrib
        points.append(point)
    dos["totaldos"]["diagram"]["point"] = {}
    i = 1
    for item in points:
        name = str(i)
        dos["totaldos"]["diagram"]["point"][name] = item
        i = i+1
    return dos
    
#parser for KERR.OUT
def parse_kerr(name):
    try:
        data = np.genfromtxt(name)
    except:
        raise ParseError
    x = data[:,0]
    y = data[:,1]

    k=0
    startval = x[0]
    for i in range(1,len(x)):
        if x[i]==0.0:
            k=i
            break
    
    out = {}
    out['x']  = list(x[:k])
    out['re'] = list(y[:k])
    out['im'] = list(y[k:])

    return out

#parser for epsilon
def parse_epsilon(name):
    try:
        data = np.genfromtxt(name)
    except:
        raise ParseError
    out = {}
    out['e']  = list(data[:,0])
    out['re'] = list(data[:,1])
    out['im'] = list(data[:,2])

    return out

#parser for CHI_111.OUT
def parse_chi(name):
    try:
        data = np.genfromtxt(name)
    except:
        raise ParseError
    out = {}
    out['0'] = list(data[:,0])
    out['1'] = list(data[:,1])
    out['2'] = list(data[:,2])
    out['3'] = list(data[:,3])

    return out

#parser for ELNES.OUT
def parse_elnes(name):
    try:
        data = np.genfromtxt(name)
    except:
        raise ParseError
    out = {}
    out['0'] = list(data[:,0])
    out['1'] = list(data[:,1])

    return out

#parser for SEEBECK_11.OUT
def parse_seebeck(name): 
    try:
        data = np.genfromtxt(name)
    except:
        raise ParseError
    out = {}
    out['t']  = list(data[:,0])
    out['mu'] = list(data[:,1])
    out['realpart'] = list(data[:,2])
    out['imaginarypart'] = list(data[:,2])

    return out

#parser for ldos.out
def parse_ldos(name):
    try:
        data = np.genfromtxt(name)
    except:
        raise ParseError
    out = {}
    out['0'] = list(data[:,0])
    out['1'] = list(data[:,1])

    return out
    
#parser for band_edges.out
def parse_band_edges(name):
    try:
        data = np.genfromtxt(name)
    except:
        raise ParseError
    out = {}
    out['0'] = list(data[:,0])
    out['1'] = list(data[:,1])
    out['2'] = list(data[:,2])

    return out

#parser for spintext.xml
def parse_spintext(name):
    try:
        treespin = ET.parse(name)
    except:
        raise ParseError
    root = treespin.getroot()
    spintext = {}
    bands = []
    for node in root.findall('band'):
        band = node.attrib
        bands.append(band)
        kpoints = []
        for subnode in node.findall('k-point'):
            kpoint = subnode.attrib
            kpoints.append(kpoint)
        band['kpoint'] = {}
        i = 1
        for item in kpoints:
            band['kpoint'][str(i)] = item
            band['kpoint'][str(i)]['vec'] = item['vec'].split()
            band['kpoint'][str(i)]['spin'] = item['spin'].split()
            i += 1
    spintext['band'] = {}
    for item in bands:
        spintext['band'][item['ist']] = item
    return spintext

#parser for POLARIZATION.OUT
def parse_polarization(name):
    file = open(name)
    lines = []
    for i in range(len(open(name).readlines())):
        line = next(file)
        if '#' not in line:
            lines.append(line.split())
    polarization = {}
    polarization['total'] = lines[0]
    polarization['electronic'] = lines[1]
    polarization['ionic'] = lines[2]
    return polarization

#parser for TDOS_WANNIER.OUT
def parse_tdos_wannier(name):
    try:
        data = np.genfromtxt(name)
    except:
        raise ParseError
    out = {}
    out['0'] = list(data[:,0])
    out['1'] = list(data[:,1])

    return out

#parser for WANNIER_INFO.OUT
def parse_wannier_info(name):
    file = open(name)
    wannier = {}
    lines = []
    list = []
    total = []
    start = False
    for i in range(len(open(name).readlines())):
        line = next(file)
        if "* Wannier functions" in line:
            start = True
        if start:
            lines.append(line)
    for i in range(len(lines)):
        if lines[i].strip().startswith( '1'):
            for j in range (4):
                list.append(lines[i+j].split())
        if lines[i].strip().startswith( '5'):
            for j in range (4):
                list.append(lines[i+j].split())
        if lines[i].strip().startswith( 'total'):
            total.append(lines[i].split())
            
    wannier['vec'] = {}
    wannier['Omega'] = {}
    wannier['Omega_I'] = {}
    wannier['Omega_D'] = {}
    wannier['Omega_OD'] = {}
    i = 0
    for item in list:
        wannier['vec'][str(i)] = item[1:4]
        wannier['Omega'][str(i)] = item[4]
        wannier['Omega_I'][str(i)] = item[5]
        wannier['Omega_D'][str(i)] = item[6]
        wannier['Omega_OD'][str(i)] = item[7]
        i += 1

    wannier['vec']['total'] = {}
    wannier['Omega']['total'] = {}
    wannier['Omega_I']['total'] = {}
    wannier['Omega_D']['total'] = {}
    wannier['Omega_OD']['total'] = {}

    j = 0
    for item in total:
        wannier['Omega']['total'][str(j)] = item[1]
        wannier['Omega_I']['total'][str(j)] = item[2]
        wannier['Omega_D']['total'][str(j)] = item[3]
        wannier['Omega_OD']['total'][str(j)] = item[4]
        j += 1

    return(wannier)
    
    # parser for coreoverlap.xml
    def parse_core_overlap(name):
        try:
            tree= ET.parse(name)
        except:
            raise ParseError
        
        root=tree.getroot()
        core_overlap = {}
        core_overlap['nkpt'] = root.attrib['nkpt']
        core_overlap['nstfv'] = root.attrib['nstfv']
        core_overlap['ncg'] = root.attrib['ncg']

        kpoints = []
        for kpoint in root:
            kpt = {}
            kpt["index"] = kpoint.attrib['index']
            pairs = []
            for pair_xml in kpoint:
                pair = pair_xml.attrib
                pair["overlap"] = pair["overlap"].split()
                pairs.append(pair)
            kpt["pairs"] = pairs
            kpoints.append
        core_overlap["kpoints"] = kpoints

        return core_overlap

class excitingPlot3D():

    def __init__(self, filePath):
        self.path  = filePath                         #path of the RHO3D.xml, VCL3D.xml, VXC3D.xml, WF3D.xml, ELF3D.xml, EF3D.xml                                   
        if not os.path.exists(self.path):          #check if the path exists                                      
            raise OSError('Path not valid.')
        try:
            self.data  = parse_plot_3D(self.path)         #data in form of a python dictionary                   
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
        data = parse_plot_3D(self.path)
        for i in range (0,l):
           try:
               data = data[keys[i]] #value of the given datapath                                                  
           except KeyError:      #checks if the datapath exists                                                   
               print ("invalid key: ")
               print (keys[i])
               sys.exit()
        return data

class excitingLSJ():

    def __init__(self, filePath):
        self.path  = filePath                         #path of the LSJ.xml                                   
        if not os.path.exists(self.path):          #check if the path exists                                      
            raise OSError('Path not valid.')
        try:
            self.data  = parse_LSJ(self.path)         #data in form of a python dictionary                   
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
        data = parse_LSJ(self.path)
        for i in range (0,l):
           try:
               data = data[keys[i]] #value of the given datapath                                                  
           except KeyError:      #checks if the datapath exists                                                   
               print ("invalid key: ")
               print (keys[i])
               sys.exit()
        return data

class excitingEFG():

    def __init__(self, filePath):
        self.path  = filePath                         #path of the EFG.xml                                   
        if not os.path.exists(self.path):          #check if the path exists                                      
            raise OSError('Path not valid.')
        try:
            self.data  = parse_EFG(self.path)         #data in form of a python dictionary                   
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
        data = parse_EFG(self.path)
        for i in range (0,l):
           try:
               data = data[keys[i]] #value of the given datapath                                                  
           except KeyError:      #checks if the datapath exists                                                   
               print ("invalid key: ")
               print (keys[i])
               sys.exit()
        return data

class excitingMossbauer():

    def __init__(self, filePath):
        self.path  = filePath                         #path of the mossbauer.xml                                   
        if not os.path.exists(self.path):          #check if the path exists                                      
            raise OSError('Path not valid.')
        try:
            self.data  = parse_mossbauer(self.path)         #data in form of a python dictionary                   
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
        data = parse_mossbauer(self.path)
        for i in range (0,l):
           try:
               data = data[keys[i]] #value of the given datapath                                                  
           except KeyError:      #checks if the datapath exists                                                   
               print ("invalid key: ")
               print (keys[i])
               sys.exit()
        return data

class excitingExpiqr():

    def __init__(self, filePath):
        self.path  = filePath                         #path of the expiqr.xml                                   
        if not os.path.exists(self.path):          #check if the path exists                                      
            raise OSError('Path not valid.')
        try:
            self.data  = parse_expiqr(self.path)         #data in form of a python dictionary                   
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
        data = parse_expiqr(self.path)
        for i in range (0,l):
           try:
               data = data[keys[i]] #value of the given datapath                                                  
           except KeyError:      #checks if the datapath exists                                                   
               print ("invalid key: ")
               print (keys[i])
               sys.exit()
        return data

class excitingEffmass():

    def __init__(self, filePath):
        self.path  = filePath                         #path of the effmass.xml                                   
        if not os.path.exists(self.path):          #check if the path exists                                      
            raise OSError('Path not valid.')
        try:
            self.data  = parse_effmass(self.path)         #data in form of a python dictionary                   
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
        data = parse_effmass(self.path)
        for i in range (0,l):
           try:
               data = data[keys[i]] #value of the given datapath                                                  
           except KeyError:      #checks if the datapath exists                                                   
               print ("invalid key: ")
               print (keys[i])
               sys.exit()
        return data

class excitingBandstructure():

    def __init__(self, filePath):
        self.path  = filePath                         #path of the bandstructure.xml                                   
        if not os.path.exists(self.path):          #check if the path exists                                      
            raise OSError('Path not valid.')
        try:
            self.data  = parse_bandstructure(self.path)         #data in form of a python dictionary                   
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
        data = parse_bandstructure(self.path)
        for i in range (0,l):
           try:
               data = data[keys[i]] #value of the given datapath                                                  
           except KeyError:      #checks if the datapath exists                                                   
               print ("invalid key: ")
               print (keys[i])
               sys.exit()
        return data

class excitingDos():

    def __init__(self, filePath):
        self.path  = filePath                         #path of the dos.xml                                   
        if not os.path.exists(self.path):          #check if the path exists                                      
            raise OSError('Path not valid.')
        try:
            self.data  = parse_dos(self.path)         #data in form of a python dictionary                   
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
        data = parse_dos(self.path)
        for i in range (0,l):
           try:
               data = data[keys[i]] #value of the given datapath                                                  
           except KeyError:      #checks if the datapath exists                                                   
               print ("invalid key: ")
               print (keys[i])
               sys.exit()
        return data

class excitingKerr():

    def __init__(self, filePath):
        self.path  = filePath                         #path of the KERR.OUT                                   
        if not os.path.exists(self.path):          #check if the path exists                                      
            raise OSError('Path not valid.')
        try:
            self.data  = parse_kerr(self.path)         #data in form of a python dictionary                   
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
        data = parse_kerr(self.path)
        for i in range (0,l):
           try:
               data = data[keys[i]] #value of the given datapath                                                  
           except KeyError:      #checks if the datapath exists                                                   
               print ("invalid key: ")
               print (keys[i])
               sys.exit()
        return data

class excitingEpsilon():

    def __init__(self, filePath):
        self.path  = filePath                         #path of the EPSILON.OUT                           
        if not os.path.exists(self.path):          #check if the path exists                             
            raise OSError('Path not valid.')
        try:
            self.data  = parse_epsilon(self.path)         #data in form of a python dictionary         
        except ParseError:                       #checks if the parser works                           
            raise ErrornousFileError

    def returnData(self, dataPath):
        '''                                                                                             
        returns the data, saved at the dataPath.                                                        
        In:                                                                                             
        dataPath     string     path of the data                                                     
        '''
        keys = dataPath.split('/') #datapath splitted in to single keys                                  
        l = len(keys)
        data = parse_epsilon(self.path)
        for i in range (0,l):
            try:
                data = data[keys[i]] #value of the given datapath                                        
            except KeyError:      #checks if the datapath exists                                         
                print ("invalid key: ")
                print (keys[i])
                sys.exit()
        return data

class excitingChi():

    def __init__(self, filePath):
        self.path  = filePath                         #path of the CHI_111.OUT                          
        if not os.path.exists(self.path):          #check if the path exists                             
            raise OSError('Path not valid.')
        try:
            self.data  = parse_chi(self.path)         #data in form of a python dictionary
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
        data = parse_chi(self.path)
        for i in range (0,l):
            try:
                data = data[keys[i]] #value of the given datapath                                        
            except KeyError:      #checks if the datapath exists                                         
                print ("invalid key: ")
                print (keys[i])
                sys.exit()
        return data

class excitingElnes():

    def __init__(self, filePath):
        self.path  = filePath                         #path of the ELNES.OUT                          
        if not os.path.exists(self.path):          #check if the path exists                             
            raise OSError('Path not valid.')
        try:
            self.data  = parse_elnes(self.path)         #data in form of a python dictionary
        except ParseError:                        #checks if the parser works                           
            raise ErrornousFileError

    def returnData(self, dataPath):
        '''                                                                                             
        returns the data, saved at the dataPath.                                                        
        In:                                                                                             
        dataPath     string     path of the data                                                     
        '''
        keys = dataPath.split('/') #datapath splitted in to single keys                                  
        l = len(skeys)
        data = parse_elnes(self.path)
        for i in range (0,l):
            try:
                data = data[keys[i]] #value of the given datapath                                        
            except KeyError:      #checks if the datapath exists                                         
                print ("invalid key: ")
                print (keys[i])
                sys.exit()
        return data

class excitingSeebeck():

    def __init__(self, filePath):
        self.path  = filePath 
        if not os.path.exists(self.path):  
            raise OSError('Path not valid.')
        try:
            self.data  = parse_seebeck(self.path)     
        except ParseError:                                              
            raise ErrornousFileError

    def returnData(self, dataPath):
        '''                                                                                             
        returns the data, saved at the dataPath.                                                        
        In:                                                                                             
        dataPath     string     path of the data                                                     
        '''
        keys = dataPath.split('/') #datapath splitted in to single keys                                  
        l = len(keys)
        data = parse_seebeck(self.path)
        for i in range (0,l):
            try:
                data = data[keys[i]] #value of the given datapath                                        
            except KeyError:      #checks if the datapath exists                                         
                print ("invalid key: ")
                print (keys[i])
                sys.exit()
        return data
    
class excitingLdos():

    def __init__(self, filePath):
        self.path  = filePath                         #path of the ldos.out   
        if not os.path.exists(self.path):          #check if the path exists        
            raise OSError('Path not valid.')
        try:
            self.data  = parse_ldos(self.path)         #data in form of a python dictionary     
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
        data = parse_ldos(self.path)
        for i in range (0,l):
            try:
                data = data[keys[i]] #value of the given datapath       
            except KeyError:      #checks if the datapath exists           
                print ("invalid key: ")
                print (keys[i])
                sys.exit()
        return data

class excitingBand_edges():

    def __init__(self, filePath):
        self.path  = filePath                         #path of the band_edges.out                              
        if not os.path.exists(self.path):          #check if the path exists                             
            raise OSError('Path not valid.')
        try:
            self.data  = parse_band_edges(self.path)         #data in form of a python dictionary              
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
        data = parse_band_edges(self.path)
        for i in range (0,l):
            try:
                data = data[keys[i]] #value of the given datapath                                        
            except KeyError:      #checks if the datapath exists                                         
                print ("invalid key: ")
                print (keys[i])
                sys.exit()
        return data

class excitingSpintext():

    def __init__(self, filePath):
        self.path  = filePath                         #path of the spintext.xml            
        if not os.path.exists(self.path):          #check if the path exists   
            raise OSError('Path not valid.')
        try:
            self.data  = parse_spintext(self.path)         #data in form of a python dictionary
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
        data = parse_spintext(self.path)
        for i in range (0,l):
            try:
                data = data[keys[i]] #value of the given datapath      
            except KeyError:      #checks if the datapath exists         
                print ("invalid key: ")
                print (keys[i])
                sys.exit()
        return data

class excitingPolarization():

    def __init__(self, filePath):
        self.path  = filePath                         #path of the POLARIZATION.OUT    
        if not os.path.exists(self.path):          #check if the path exists    
            raise OSError('Path not valid.')
        try:
            self.data  = parse_polarization(self.path)         #data in form of a python dictionary
        except ParseError:                        #checks if the parser works           
            raise ErrornousFileError
    def returnData(self, dataPath):
        '''                                                                                                                                                                                                                returns the data, saved at the dataPath.                                                                                                                                                              
        In:                                                                                                                                                              
        dataPath     string     path of the data                                                                
        '''
        keys = dataPath.split('/') #datapath splitted in to single keys                      
        l = len(keys)
        data = parse_polarization(self.path)
        for i in range (0,l):
            try:
                data = data[keys[i]] #value of the given datapath        
            except KeyError:      #checks if the datapath exists                       
                print ("invalid key: ")
                print (keys[i])
                sys.exit()
        return data

class excitingWannierInfo():

    def __init__(self, filePath):
        self.path  = filePath                         #path of the WANNIER_INFO.OUT   
        if not os.path.exists(self.path):          #check if the path exists          
            raise OSError('Path not valid.')
        try:
            self.data  = parse_wannier_info(self.path)         #data in form of a python dictionary
        except:                        #checks if the parser works                      
            raise ErrornousFileError
    def returnData(self, dataPath):
        '''                                                                                        
        returns the data, saved at the dataPath.                                                   
        In:                                                                                        
        dataPath     string     path of the data                                                   
        '''
        keys = dataPath.split('/') #datapath splitted in to single keys                            
        l = len(keys)
        data = parse_wannier_info(self.path)
        for i in range (0,l):
            try:
                data = data[keys[i]] #value of the given datapath                                  
            except KeyError:      #checks if the datapath exists                                   
                print ("invalid key: ")
                print (keys[i])
                sys.exit()
        return data

class excitingTdosWannier():

    def __init__(self, filePath):
        self.path  = filePath                         #path of the TDOS_WANNIER.OUT
        if not os.path.exists(self.path):          #check if the path exists   
            raise OSError('Path not valid.')
        try:
            self.data  = parse_tdos_wannier(self.path)         #data in form of a python dictionary
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
        data = parse_tdos_wannier(self.path)
        for i in range (0,l):
            try:
                data = data[keys[i]] #value of the given datapath                                    
            except KeyError:      #checks if the datapath exists                                     
                print ("invalid key: ")
                print (keys[i])
                sys.exit()
        return data
        
class coreOverlap():
    def __init__(self, filePath):
        self.path  = filePath                      #path of the RHO3D.xml, VCL3D.xml, VXC3D.xml, WF3D.xml, ELF3D.xml, EF3D.xml                                   
        if not os.path.exists(self.path):          #check if the path exists                                      
            raise OSError('Path not valid.')
        try:
            self.data  = parse_core_overlap(self.path)     #data in form of a python dictionary                   
        except ParseError:                                 #checks if the parser works                                     
            raise ErrornousFileError
