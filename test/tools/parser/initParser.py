import os
import xml.etree.ElementTree as ET
from xml.etree.ElementTree import ParseError

def parseInit(fileName):
    root = ET.parse(fileName).getroot()

    init =  root.attrib
    init['tests'] = []
    for child in root:
        test = child.attrib
        try:
            test['tolValuePair'] = child.find('tolValuePair').attrib
        except:
            test['tolValuePair'] = {}
            test['tolValuePair']['tol'] = '0'
            test['tolValuePair']['value'] = ''
        try:
            test['eigval'] = child.find('eigval').attrib
        except:
            test['eigval'] = {}
        init['tests'].append(test)

    for i in range(0, len(init['tests'])):
        if 'tolIter' not in init['tests'][i]:
            init['tests'][i]['tolIter'] = 0
        if 'maxIter' not in init['tests'][i]:
            init['tests'][i]['maxIter'] = 0
    return init

def getInitFile(directory):
    path = '.'
    files = [i for i in os.listdir(path) if os.path.isfile(os.path.join(path,i)) and 'init' in i]
    if len(files)==0:
        print(' No init file in %s. '%directory)
        raise FileNotFoundError
    elif len(files)>1:
        print(' More than one init file in %s. Skip test.'%directory)
        raise OSError 
    else:
        return os.path.join(path,files[0])

    
