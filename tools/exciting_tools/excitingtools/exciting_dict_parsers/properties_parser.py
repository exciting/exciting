"""Parsers for exciting properties.
"""
import xml.etree.cElementTree as ET
from xml.etree.ElementTree import ParseError
import numpy as np
import os
from typing import Dict

from excitingtools.parser_utils.parser_decorators import xml_root


def parse_plot_3d(name: str) -> dict:
    """
    Parser for RHO3D.xml, VCL3D.xml, VXC3D.xml, WF3D.xml, ELF3D.xml, EF3D.xmlit

    :param str name: File name
    :return dict output: Parsed data
    """
    root = ET.parse(name)
    plot_3d = {"title": root.find("title").text}
    grid = root.find("grid").attrib
    axis = []
    for ax in root.find("grid").findall("axis"):
        axis.append(ax.attrib)
    grid["axis"] = {}
    for item in axis:
        name = item["name"]
        grid["axis"][name] = item
    grid["value"] = root.find("grid").find("value").attrib
    plot_3d["grid"] = grid

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
    plot_3d["function"] = function

    return plot_3d


def parse_lsj(name: str) -> dict:
    """
    Parser for LSJ.xml

    :param str name: File name
    :return dict output: Parsed data
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


def parse_efg(name: str) -> dict:
    """
    Parser for EFG.xml

    Returns a dictionary of the form:
      data = {'species1': {'chemicalSymbol': chemicalSymbol,
                          'atom1': { 'trace': trace,
                                     'efg: efg,
                                     'eigenvalues': eigenvalues
                                   },
                          'atom2': {},
                          }
              'species2':...
            }

    :param str name: File name
    :return dict output: Parsed data
    """
    root = ET.parse(name).getroot()
    data = {}

    for species in root.findall('species'):
        species_key = species.tag + str(species.attrib['n'])
        data[species_key] = {
            'chemicalSymbol': species.attrib['chemicalSymbol']
        }

        for atom in species.findall('atom'):
            atom_key = atom.tag + atom.attrib['n']

            for efg_tensor in atom.findall("EFG-tensor"):
                efg = np.empty(shape=(3, 3))

                for i, line in enumerate(efg_tensor.findall("line")):
                    efg[i, :] = [float(x) for x in line.text.split()]

                for eigenvalues in efg_tensor.findall("eigenvalues"):
                    eig = [float(e) for e in eigenvalues.text.split()]

            data[species_key][atom_key] = {
                'trace': float(efg_tensor.attrib['trace']),
                'efg': efg,
                'eigenvalues': eig
            }
    return data


def parse_mossbauer(name: str) -> dict:
    """
    Parser for mossbauer.xml

    :param str name: File name
    :return dict output: Parsed data
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


def parse_expiqr(name: str) -> dict:
    """
    Parser for expiqr.xml

    :param str name: File name
    :return dict output: Parsed data
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


def parse_effmass(name: str) -> dict:
    """
    Parser for effmass.xml

    :param str name: File name
    :return dict output: Parsed data
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
        evdk["evdk_matrix"] = matrix1

        emt = node.find("emt").attrib
        matrix2 = []
        for line in node.find("emt").findall("line"):
            matrix2.append(line.text.split())
        emt["emt_matrix"] = matrix2

        st["evdk"] = evdk
        st["emt"] = emt
        state.append(st)

    effmass["state"] = {}
    for item in state:
        name = item["n"]
        effmass["state"][name] = item

    return effmass


# TODO(Hannah). Issue 138. Ensure test cases work with `parse_bandstructure` and remove `parse_bandstructure_depreciated`
# This parser is depreciated. Please do not use.
def parse_bandstructure_depreciated(name: str) -> dict:
    """
    Parser for bandstructure.xml.

    Used for parsing in the test framework, as returns a dict.

    :param str name: File name
    :return dict output: Parsed data
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
            i = i + 1
    bandstructure["band"] = {}
    j = 1
    for item in bands:
        name = str(j)
        bandstructure["band"][name] = item
        j = j + 1
    return bandstructure


@xml_root
def parse_band_structure_xml(root) -> dict:
    """ Parse KS band structure from bandstructure.xml.

    :param root: Band structure XML file name, XML string or ElementTree.Element as input.
    :return: Band data
    """
    # Split band structure file contents: title, bands and vertices
    bs_xml: Dict[str, list] = {'title': [], 'band': [], 'vertex': []}

    for item in list(root):
        try:
            bs_xml[item.tag].append(item)
        except KeyError:
            raise KeyError(f'Element tag {item.tag} requires implementing in band structure parser')

    n_bands = len(bs_xml['band'])
    first_band = bs_xml['band'][0]
    n_kpts = len(list(first_band))

    # Same set of flattened k-points, per band - so parse once
    k_points_along_band = np.array([point.get('distance') for point in list(first_band)], dtype=float)

    # Read E(k), per band
    band_energies = np.empty(shape=(n_kpts, n_bands))
    for ib, band in enumerate(bs_xml['band']):
        for ik, point in enumerate(list(band)):
            band_energies[ik, ib] = point.get('eval')

    vertices = []
    for element in bs_xml['vertex']:
        vertices.append({'distance': float(element.get('distance')),
                         'label': element.get('label'),
                         'coord':  [float(x) for x in element.get('coord').split()]})

    return {'title': bs_xml['title'], 'n_kpts': n_kpts, 'n_bands': n_bands,
            'k_points_along_band': k_points_along_band,
            'band_energies': band_energies, 'vertices': vertices}


def parse_band_structure_dat(name: str) -> dict:
    """Parser for bandstructure.dat

    :param str name: File name
    :return dict output: Parsed data
    """
    bs_dat = np.loadtxt(name)
    with open(name) as f:
        header = f.readline()

    n_kpts = int(header.split()[3])
    n_bands = int(header.split()[2])
    dimensions = 3

    k_points = np.empty(shape=(n_kpts, dimensions))
    flattened_k_points = np.empty(n_kpts)
    for i in range(n_kpts):
        k_points[i] = np.array([k for k in bs_dat[i, 2:5]])
        flattened_k_points[i] = bs_dat[i, 5]

    band_energies = np.reshape(bs_dat[:, 6], (n_kpts, n_bands), order='F')

    return {
            'n_kpts': n_kpts,
            'n_bands': n_bands,
            'k_points': k_points,
            'flattened_k_points': flattened_k_points,
            'band_energies': band_energies
        }


def parse_dos(name: str) -> dict:
    """
    Parser for dos.xml

    :param str name: File name
    :return dict output: Parsed data
    """
    root = ET.parse(name).getroot()
    dos = {}
    dos["title"] = root.find("title").text
    totaldos = root.find("totaldos").attrib
    dos["totaldos"] = totaldos
    diagram = root.find("totaldos").find("diagram").attrib
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
        i = i + 1
    return dos


@xml_root
def parse_charge_density(root) -> np.ndarray:
    """ Parse charge density from RHO1D.xml file.

    `axis` and `vertex` sub-trees ignored in the parsing.

    :param root: XML file name, XML string or ElementTree.Element as input.
    :return: Numpy array containing rho[:, 1] = distance and rho[:, 2] = density.
    """
    function_points = root.find('grid').find('function')
    rho = np.empty(shape=(len(function_points), 2))
    for i, point in enumerate(function_points):
        rho[i, :] = [point.attrib['distance'], float(point.attrib['value'])]
    return rho


def parse_kerr(name: str) -> dict:
    """
    Parser for KERR.OUT

    :param str name: File name
    :return dict output: Parsed data
    """
    try:
        data = np.genfromtxt(name, skip_header=1)
    except:
        raise ParseError

    out = {'energy': data[:, 0], 're': data[:, 1], 'im': data[:, 2]}

    return out


def parse_epsilon(name: str) -> dict:
    """
    Parser for EPSILON_ij.OUT

    :param str name: File name
    :return dict output: Parsed data
    """
    try:
        data = np.genfromtxt(name, skip_header=1)
    except:
        raise ParseError
    out = {'energy': data[:, 0], 're': data[:, 1], 'im': data[:, 2]}
    return out


def parse_chi(name: str) -> dict:
    """
    Parser for CHI_111.OUT

    :param str name: File name
    :return dict output: Parsed data
    """
    try:
        data = np.genfromtxt(name, skip_header=1)
    except:
        raise ParseError
    out = {
        'energy': data[:, 0],
        're': data[:, 1],
        'im': data[:, 2],
        'modulus': data[:, 3]
    }
    return out


def parse_elnes(name: str) -> dict:
    """
    Parser for ELNES.OUT

    :param str name: File name
    :return dict output: Parsed data
    """
    try:
        data = np.genfromtxt(name)
    except:
        raise ParseError
    out = {'energy': data[:, 0], 'elnes': data[:, 1]}
    return out


def parse_seebeck(name: str) -> dict:
    """
    Parser for SEEBECK_11.OUT

    :param str name: File name
    :return dict output: Parsed data
    """
    try:
        data = np.genfromtxt(name)
    except:
        raise ParseError
    out = {
        'temperature': data[:, 0],
        'mu': data[:, 1],
        're': data[:, 2],
        'im': data[:, 3]
    }

    return out


def parse_ldos(name: str) -> dict:
    """
    Parser for ldos.out

    :param str name: File name
    :return dict output: Parsed data
    """
    try:
        data = np.genfromtxt(name)
    except:
        raise ParseError
    out = {'energy': data[:, 0], 'ldos': data[:, 1]}

    return out


def parse_band_edges(name: str) -> dict:
    """
    Parser for band_edges.out

    Keys
    * c_axis corresponds to the linear grid along the magnitude
      of the c vector of the unit cell.
    * VBM = Valence band maximum
    * CBm = Conduction band minimum

    :param str name: File name
    :return dict output: Parsed data
    """
    try:
        data = np.genfromtxt(name)
    except:
        raise ParseError
    out = {'c_axis': data[:, 0], 'VBM': data[:, 1], 'CBm': data[:, 2]}

    return out


def parse_spintext(name: str) -> dict:
    """
    Parse spintext.xml

    TODO(Bene) Issue 87 Refactor to return a dict

    Each element of the list contains a dict with keys:
     ['ist', 'k-point', 'spin', 'energy']

    :param str name: Path to the spintext.xml that will be parsed
    :return dict spintext: List that holds the parsed spintexture.xml
    """
    # parse file
    file_name = 'spintext.xml'
    if name.split('/')[-1] != file_name:
        name = os.path.join(name, file_name)

    tree_spin = ET.parse(name)
    root_spin = tree_spin.getroot()

    spintext = {}
    i = 0
    for band in root_spin.findall("band"):
        k_point = []
        spin = []
        energy = []

        for val in band.findall("k-point"):
            k_point.append([float(k) for k in val.attrib["vec"].split()])
            spin.append([float(s) for s in val.attrib["spin"].split()])
            energy.append(float(val.attrib["energy"]))

        spintext[str(i)] = {
            "ist": int(band.attrib["ist"]),
            "k-point": k_point,
            "spin": spin,
            "energy": energy
        }
        i += 1

    return spintext


def parse_polarization(name: str) -> dict:
    """
    Parser for POLARIZATION.OUT

    :param str name: File name
    :return dict output: Parsed data
    """
    file = open(name)
    lines = []
    for i in range(len(open(name).readlines())):
        line = next(file)
        if '#' not in line:
            lines.append(line.split())
    polarization = {
        'total': lines[0],
        'electronic': lines[1],
        'ionic': lines[2]
    }
    return polarization


def parse_tdos_wannier(name: str) -> dict:
    """
    Parser for TDOS_WANNIER.OUT

    :param str name: File name
    :return dict output: Parsed data
    """
    try:
        data = np.genfromtxt(name)
    except:
        raise ParseError
    out = {'energy': data[:, 0], 'dos': data[:, 1]}

    return out


def parse_wannier_info(name: str) -> dict:
    """
    Parser for WANNIER_INFO.OUT

    :param str name: File name
    :return dict output: Parsed data
    """
    file = open(name)

    # Extract data
    lines = []
    data = []
    total = []
    start = False
    for i, line in enumerate(file.readlines()):
        if "* Wannier functions" in line:
            start = True
        if start:
            lines.append(line)
    for i in range(len(lines)):
        if lines[i].strip().startswith('1'):
            for j in range(4):
                data.append(lines[i + j].split())
        if lines[i].strip().startswith('5'):
            for j in range(4):
                data.append(lines[i + j].split())
        if lines[i].strip().startswith('total'):
            total.append(lines[i].split())
    file.close()

    # Package data into dictionary
    n_wannier = len(data)
    localisation_center = np.empty(shape=(n_wannier, 3))
    wannier = {
        'n_wannier': n_wannier,
        'Omega': [],
        'Omega_I': [],
        'Omega_D': [],
        'Omega_OD': []
    }

    for i, item in enumerate(data):
        localisation_center[i, :] = [float(x) for x in item[1:4]]
        wannier['Omega'].append(float(item[4]))
        wannier['Omega_I'].append(float(item[5]))
        wannier['Omega_D'].append(float(item[6]))
        wannier['Omega_OD'].append(float(item[7]))

    wannier['localisation_center'] = localisation_center

    totals = {'Omega': [], 'Omega_I': [], 'Omega_D': [], 'Omega_OD': []}
    for j, item in enumerate(total):
        totals['Omega'].append(float(item[1]))
        totals['Omega_I'].append(float(item[2]))
        totals['Omega_D'].append(float(item[3]))
        totals['Omega_OD'].append(float(item[4]))

    wannier['total'] = totals

    return wannier


def parse_core_overlap(name: str) -> dict:
    """
    Parser for coreoverlap.xml

    Parsed dictionary has the structure:

        output = {'nkpt':  nkpt
                  'nstfv': nstfv
                  'ncg':   ncg
                  'kpoints': [{'index': index, 'pairs': pairs},
                              {'index': index, 'pairs': pairs},
                              ...]
                  }

    where output['kpoints'][ik]['pairs'] =
      [{'ist1': '1', 'ist2': '1', 'de': 12.97849772, 'overlap': 3.35753859e-07},
       {'ist1': '1', 'ist2': '2', 'de': 12.97849772, 'overlap': 3.35753859e-07},
       ...
       n_pairs]

    :param str name: File name
    :return dict output: Parsed data
    """
    try:
        tree = ET.parse(name)
    except:
        raise ParseError

    root = tree.getroot()
    core_overlap = {
        'nkpt': int(root.attrib['nkpt']),
        'nstfv': int(root.attrib['nstfv']),
        'ncg': int(root.attrib['ncg'])
    }

    k_points = []
    for k_point in root:
        kpt = {"index": int(k_point.attrib['index'])}
        pairs = []
        for pair_xml in k_point:
            pair = pair_xml.attrib
            pair['ist1'] = int(pair['ist1'])
            pair['ist2'] = int(pair['ist2'])
            pair["de"] = float(pair["de"])
            pair["overlap"] = float(pair["overlap"].split()[0]) ** 2 + float(
                pair["overlap"].split()[1]) ** 2
            pairs.append(pair)
        kpt["pairs"] = pairs
        k_points.append(kpt)
    core_overlap["kpoints"] = k_points

    return core_overlap


def parse_lossfunction(fname: str) -> tuple:
    """
    Parses files containing loss function
    e.g. LOSS_FXCRPA_OC11_QMT001.OUT

    :param str fname: name of the file
    """
    xdata = []
    ydata = []
    file = open(fname, 'r')
    for lines in file:
        if 'Frequency' in lines:
            break
    for lines in file:
        data = lines.split()
        xdata.append(float(data[0]))
        ydata.append(float(data[1]))
    file.close()
    return xdata, ydata


def parse_wf1d(fname: str) -> dict:
  """
  Parses files containing one dimensional wave function plot as saved in
  the files, _e.g._, wf1d-0001-0001.dat, where the first 0001 indicates the k-point
  and the second the state index.
  
  :param str fname: name of the file
  """

  data = np.genfromtxt(fname)

  output = {}
  output["path"] = data[:, 0]
  output["|psi|^2"] = data[:, 1]
  output["Re(psi)"] = data[:, 2]
  output["Im(psi)"] = data[:, 3]

  return output 


def parse_wf2d(fname: str) -> dict:
  """
  Parses files containing the two dimensional wave function plot as saved in the files,
  _e.g._, wf2d-0001-0001.xsf, where the first 0001 indicates the k-point
  and the second the state index.

  Does not parse PRIMVEC and PRIMCOORD. These can be parsed with ase.io.read(fname).

  :param str fname: name of the file
  """

  file = open(fname, 'r')
  lines = file.readlines()
  output = {}
  
  i = 0
  while i < len(lines):
    if 'BEGIN_BLOCK_DATAGRID_2D' in lines[i]:
      i = i + 1
      data_name = lines[i].replace('\n', '').lstrip()
      
      i = i + 2
      grid = np.fromstring(lines[i], dtype=int, sep=' ')

      i = i + 1
      origin = np.fromstring(lines[i], dtype=np.double, sep = ' ')
      
      i = i + 1
      point1 = np.fromstring(lines[i], dtype=np.double, sep = ' ')
      
      i = i + 1
      point2 = np.fromstring(lines[i], dtype=np.double, sep = ' ')

      i = i + 1
      data = np.fromstring(lines[i], dtype=np.double, sep=' ')
      while 'END_DATAGRID_2D' not in lines[i]:
        i = i + 1
        new_data = np.fromstring(lines[i], dtype=np.double, sep=' ')
        data = np.concatenate((data, new_data), axis=None)
    
      output['grid'] = grid
      output['origin'] = origin
      output['point1'] = point1
      output['point2'] = point2
      output[data_name] = data # data_name: "module squared", "real", "imaginary"

    i = i + 1

  return output


def parse_wf3d(fname: str) -> dict:
  """
  Parses files containing the two dimensional wave function plot as saved in the files,
  _e.g._, wf3d-0001-0001.xsf, where the first 0001 indicates the k-point
  and the second the state index.

  Does not parse PRIMVEC and PRIMCOORD. These can be parsed with ase.io.read(fname).

  :param str fname: name of the file
  """

  file = open(fname, 'r')
  lines = file.readlines()
  output = {}
  
  i = 0
  while i < len(lines):
    if 'BEGIN_BLOCK_DATAGRID_3D' in lines[i]:
      i = i + 1
      data_name = lines[i].replace('\n', '').lstrip()
      
      i = i + 2
      grid = np.fromstring(lines[i], dtype=int, sep=' ')

      i = i + 1
      origin = np.fromstring(lines[i], dtype=np.double, sep = ' ')
      
      i = i + 1
      point1 = np.fromstring(lines[i], dtype=np.double, sep = ' ')
      
      i = i + 1
      point2 = np.fromstring(lines[i], dtype=np.double, sep = ' ')

      i = i + 1
      point3 = np.fromstring(lines[i], dtype=np.double, sep = ' ')

      i = i + 1
      data = np.fromstring(lines[i], dtype=np.double, sep=' ')
      while 'END_DATAGRID_3D' not in lines[i]:
        i = i + 1
        new_data = np.fromstring(lines[i], dtype=np.double, sep=' ')
        data = np.concatenate((data, new_data), axis=None)
    
      output['grid'] = grid
      output['origin'] = origin
      output['point1'] = point1
      output['point2'] = point2
      output['point3'] = point3
      output[data_name] = data # data_name: "module squared", "real", "imaginary"

    i = i + 1

  return output  

def parse_cube(fname: str) -> dict:
  """
  Parses .cube files. These files contain data calculated on a regular grid in a box.
  All vectors are given in cartesian coordinates. `output['cube_data']` contains the data. 
  It is described by `output['description']`.

  :param str fname: name of the file
  """

  file = open(fname, 'r')
  lines = file.readlines()
  output = {}

  output['title'] = str(lines[0])
  output['description'] = str(lines[1])

  n_atoms = int(lines[2].split()[0])
  output['n_atoms'] = n_atoms
  output['origin'] = np.array(lines[2].split()[1:], dtype=np.double)

  output['n_1'] = int(lines[3].split()[0])
  output['v_increment_1'] = np.array(lines[3].split()[1:], dtype=np.double)

  output['n_2'] = int(lines[4].split()[0])
  output['v_increment_2'] = np.array(lines[4].split()[1:], dtype=np.double)

  output['n_3'] = int(lines[5].split()[0])
  output['v_increment_3'] = np.array(lines[5].split()[1:], dtype=np.double)

  line_offset = 6
  output['atoms'] = []
  for line in lines[line_offset : line_offset + n_atoms]:
      atom_number = int(line.split()[0])
      charge = np.double(line.split()[1])
      coordinate = np.array(line.split()[2:], dtype=np.double)
      output['atoms'].append(dict(atom_number=atom_number, charge=charge, coordinate=coordinate))

  line_offset = line_offset + n_atoms

  cube_data = []
  for line in lines[line_offset:]:
      cube_data = cube_data + line.split()
  
  output['cube_data'] = np.array(cube_data, dtype=np.double)

  return output
  
