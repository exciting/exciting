"""
Module containing parsers for exciting properties
"""
import xml.etree.cElementTree as ET
from xml.etree.ElementTree import ParseError
import numpy as np


def parse_plot_3D(name):
    """
    Parser for RHO3D.xml, VCL3D.xml, VXC3D.xml, WF3D.xml, ELF3D.xml, EF3D.xmlit
    """
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


def parse_LSJ(name):
    """
    Parser for LSJ.xml
    """
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


def parse_EFG(name):
    """
    Parser for EFG.xml
    """
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


def parse_mossbauer(name):
    """
    Parser for mossbauer.xml
    """
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


def parse_expiqr(name):
    """
    parser for expiqr.xml
    """
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


def parse_effmass(name):
    """
    parser for effmass.xml
    """
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


def parse_bandstructure(name):
    """
    parser for bandstructure.xml
    """
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


def parse_dos(name):
    """
    parser for dos.xml
    """
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
    

def parse_kerr(name):
    """
    parser for KERR.OUT
    """
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


def parse_epsilon(name):
    """
    parser for epsilon
    """
    try:
        data = np.genfromtxt(name)
    except:
        raise ParseError
    out = {'e': list(data[:, 0]), 're': list(data[:, 1]), 'im': list(data[:, 2])}
    return out


def parse_chi(name):
    """
    parser for CHI_111.OUT
    """
    try:
        data = np.genfromtxt(name)
    except:
        raise ParseError
    out = {'0': list(data[:, 0]), '1': list(data[:, 1]), '2': list(data[:, 2]), '3': list(data[:, 3])}
    return out


def parse_elnes(name):
    """
    parser for ELNES.OUT
    """
    try:
        data = np.genfromtxt(name)
    except:
        raise ParseError
    out = {'0': list(data[:, 0]), '1': list(data[:, 1])}
    return out


def parse_seebeck(name):
    """
    parser for SEEBECK_11.OUT
    """
    try:
        data = np.genfromtxt(name)
    except:
        raise ParseError
    out = {'t': list(data[:, 0]), 'mu': list(data[:, 1]), 'realpart': list(data[:, 2]),
           'imaginarypart': list(data[:, 2])}

    return out


def parse_ldos(name):
    """
    parser for ldos.out
    """
    try:
        data = np.genfromtxt(name)
    except:
        raise ParseError
    out = {'0': list(data[:, 0]), '1': list(data[:, 1])}

    return out
    

def parse_band_edges(name):
    """
    parser for band_edges.out
    """
    try:
        data = np.genfromtxt(name)
    except:
        raise ParseError
    out = {'0': list(data[:, 0]), '1': list(data[:, 1]), '2': list(data[:, 2])}

    return out


def parse_spintext(name):
    """
    parser for spintext.xml
    """
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


def parse_polarization(name):
    """
    parser for POLARIZATION.OUT
    """
    file = open(name)
    lines = []
    for i in range(len(open(name).readlines())):
        line = next(file)
        if '#' not in line:
            lines.append(line.split())
    polarization = {'total': lines[0], 'electronic': lines[1], 'ionic': lines[2]}
    return polarization


def parse_tdos_wannier(name):
    """
    parser for TDOS_WANNIER.OUT
    """
    try:
        data = np.genfromtxt(name)
    except:
        raise ParseError
    out = {'0': list(data[:, 0]), '1': list(data[:, 1])}

    return out


def parse_wannier_info(name):
    """
    parser for WANNIER_INFO.OUT
    """
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

    return wannier


def parse_core_overlap(name):
    """
    parser for coreoverlap.xml
    """
    try:
        tree = ET.parse(name)
    except:
        raise ParseError

    root = tree.getroot()
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
            pair["de"] = float(pair["de"])
            pair["overlap"] = float(pair["overlap"].split()[0])**2+float(pair["overlap"].split()[0])**2
            pairs.append(pair)
        kpt["pairs"] = pairs
        kpoints.append(kpt)
    core_overlap["kpoints"] = kpoints

    return core_overlap
