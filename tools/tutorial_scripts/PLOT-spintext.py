#!/usr/bin/python2

"""
Plot square spin texture. This script is designed to work with SETUP-spintexture-plane.py together.
"""

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import gridspec
import os
import argparse as ap
import sys
import xml.etree.cElementTree as ET

# TODO(Bene) Issue #3: Remove all the function available in excitingtools from this script and import them from there


def triple_product(a, b, c):
    """
    Vector triple product, defined as 
      \mathbf{a} \cdot (\mathbf{b} \wedge \mathbf{c})

    :param a: Vector a 
    :param b: Vector b
    :param c: Vector c
    :return triple product
    """
    return np.dot(a, np.cross(b, c))

def plane_transformation(rec_lat_vec, plot_vec):
    """
    Take reciprocal lattice vectors and ONS of a plane in rec. lat. coordinates where the first two vectors span the plane and the third is normal to them 
    and calculate a matrix that transforms points in the plane to the xy plane in cartesian coordinates.
    input:
    :param rec_lat_vec:            reciprocal lattice vectors
    :param plot_vec:               ONS of the plotting plane 
    :return transformation_matrix: matrix that transforms k and spin vectors to the plot plane
    """
    norm = np.linalg.norm
    # transform plot vec in cartesian coordinates
    plot_vec = (rec_lat_vec.dot(plot_vec)).transpose()
    # extend plot vec to an orthogonal system
    plot_vec = np.array([(plot_vec[1] - plot_vec[0]) / norm(plot_vec[1] - plot_vec[0]), 
                         (plot_vec[2] - plot_vec[0]) / norm(plot_vec[2] - plot_vec[0]), 
                         np.cross(plot_vec[1] - plot_vec[0], plot_vec[2] - plot_vec[0]) \
                           / norm(np.cross(plot_vec[1] - plot_vec[0], plot_vec[2] - plot_vec[0])) ])
    transformation_matrix = np.linalg.inv(plot_vec)
    for v in transformation_matrix:
        v = v / norm(v)
    transformation_matrix = np.transpose(transformation_matrix)

    return transformation_matrix


def reciprocal_lattice_vectors(a):
    """
    Get the reciprocal lattice vectors of real-space lattice vectors \{\mathbf{a}\}:

      \mathbf{b}_0 = 2 \pi \frac{\mathbf{a}_1 \wedge \mathbf{a}_2} {\mathbf{a}_0 \cdot (\mathbf{a}_1 \wedge \mathbf{a}_2)}
      \mathbf{b}_1 = 2 \pi \frac{\mathbf{a}_2 \wedge \mathbf{a}_3} {\mathbf{a}_0 \cdot (\mathbf{a}_1 \wedge \mathbf{a}_2)}
      \mathbf{b}_2 = 2 \pi \frac{\mathbf{a}_0 \wedge \mathbf{a}_1} {\mathbf{a}_0 \cdot (\mathbf{a}_1 \wedge \mathbf{a}_2)}

    :param np.ndarray a: Lattice vectors, stored column-wise
    :return: np.ndarray b: Reciprocal lattice vectors, stored column-wise
    """
    volume = triple_product(a[:, 0], a[:, 1], a[:, 2])
    b = np.empty(shape=(3, 3))
    b[:, 0] = 2 * np.pi * np.cross(a[:, 1], a[:, 2]) / volume
    b[:, 1] = 2 * np.pi * np.cross(a[:, 2], a[:, 0]) / volume
    b[:, 2] = 2 * np.pi * np.cross(a[:, 0], a[:, 1]) / volume
    return b


def parse_lattice_vectors(file_path):
    """
    Parse the lattice coordinate from the input.xml. Units are in bohr.
    :param file_path:    relative path to the input.xml
    :return lattvec:    matrix that holds the lattice vectors.
    """
    file_name = 'input.xml'
    treeinput = ET.parse(file_path)
    if file_path.split('/')[-1] != file_name:
        file_path = os.path.join(file_path, fileName)

    treeinput = ET.parse(file_path)
    root_input = treeinput.getroot()

    try:
        scale = float(root_input.find("structure").find("crystal").attrib["scale"])
    except Exception:
        scale = 1.0
    lat_vec = [np.array(val.text.split(), dtype=float)*scale for val in root_input.find("structure").find("crystal").findall("basevect")]
    lat_vec = np.array(lat_vec).transpose()
    
    return lat_vec


def parse_spintext(name):
    """
    Parse spintext.xml
    :param name:    path to the spintext.xml that will be parsed
    :return spintext:      dictionary that holds the parsed spintexture.xml
    """
    file_name = 'spintext.xml'
    if name.split('/')[-1] != file_name:
        name = os.path.join(name, file_name)
    
    tree_spin = ET.parse(name)
    root_spin = tree_spin.getroot()
    spintext = []
    for band in root_spin.findall("band"):
        b = {}
        b["ist"] = int(band.attrib["ist"])
        b["k-point"] = [val.attrib["vec"].split() for val in band.findall("k-point")]
        b["spin"] = [val.attrib["spin"].split() for val in band.findall("k-point")]
        b["energy"] = [float(val.attrib["energy"]) for val in band.findall("k-point")]
        spintext.append(b)
    return spintext

def parse_spintext_grid(file_path):
    """
    Parse the grid from the spintext element in input.xml
    :param file_path:    relative path to the input.xml that is parsed for spin texture grid
    :return grid:        List with the number of grid points per directions
    """
    file_name = 'input.xml'
    if file_path.split('/')[-1] != file_name:
        file_path = os.path.join(file_path, file_name)

    tree_input = ET.parse(file_path)
    root_input = tree_input.getroot()
    grid = np.array(root_input.find("properties").find("spintext").find("plot2d").\
                     find("parallelogram").attrib["grid"].split(), dtype=int)
    
    return grid

def parse_spintext_plane(file_path):
    """
    Parse the vectors that span the k-plane for the spin texture calculation and extend it to an ortho normal system.
    :param file_path:    relative path to the input.xml
    :return plot_vec:    matrix that holds the ONS
    """
    file_name = 'input.xml'
    if file_path.split('/')[-1] != file_name:
        file_path = os.path.join(file_path, file_name)
    
    tree_input = ET.parse(file_path)
    root_input = tree_input.getroot()
    
    plot_vec = np.array([val.attrib["coord"].split() for val in root_input.find("properties").find("spintext").find("plot2d").find("parallelogram")], dtype=float).transpose()

    return plot_vec

def option_parser():
    """
    Parse command line inputs 

    Parse: 
        directory
        band
        contour
        contour_threshhold

    :return input_options: Dictionary of parsed command line arguments 
    """
    p = ap.ArgumentParser(description=\
        'Usage: SETUP-spintexture-plane.py -d <directory> -b <band> -c <contour> -cthr <contour_threshhold>')

    help_directory = 'Defines the relative path to the directory \
                      where input.xml and spintext.xml is stored.'

    help_band = 'Defines the band index for the plot.'

    help_contour = 'Defines the context of the contour plot. \
                     Choises are energy and spin_z'
    
    help_contourthr = 'Defines a threshhold for the contour in case of spin_z.'

    p.add_argument('-d',
                   metavar = 'directory',
                   help = help_directory,
                   type = str,
                   default = '.')
    
    p.add_argument('-b',
                   metavar = 'band',
                   help = help_band,
                   type = int,
                   default = 0)
    
    p.add_argument('-c',
                   metavar = 'contour',
                   help = help_contour,
                   type = str,
                   choices = ['energy', 'spin_z'],
                   default = 'spin_z')
    
    p.add_argument('-cthr',
                   metavar = 'contour_threshhold',
                   help = help_contourthr,
                   type = str,
                   default = 'max')

    args = p.parse_args()

    input_options = {'directory' : args.d,
                     'band' : args.b,
                     'contour' : args.c}

    try:
        cthr = float(args.cthr)
    except:
        cthr = args.cthr
        if cthr != 'max':
            raise p.error('contour_threshhold needs to be float or "max"')
    input_options['contour_threshhold'] = cthr
    
    return input_options


def main(directory, band, contour, contour_threshhold):
    """
    Plot the spin texture for a given band. 
    :param directory: Directory of the exciting calculation.
    :param band:      Number of the band for the plot.
    :param contour:   Variable that will be plotted as contour. Can be either energy or spin_z.
    :param contourthreshold: Threshhold for the contour plit. Can be either max or float. If max, the threshhold is the absolute maximum value of the contour.
    """
    file_path_input = os.path.join(directory, 'input.xml')
    file_path_spintext = os.path.join(directory, 'spintext.xml')

    
    lat_vec = parse_lattice_vectors(file_path_input)
    plot_vec = parse_spintext_plane(file_path_input)
    grid = parse_spintext_grid(file_path_input)
    
    spin_text_data = parse_spintext(file_path_spintext)
    
    if band <= spin_text_data[-1]['ist'] and band >= spin_text_data[0]['ist']:
        band_index = band - spin_text_data[0]['ist']
    else:
        sys.exit('Band must be in the range of the bands considered for the spin texture.')
    

    rec_lat_vec = reciprocal_lattice_vectors(lat_vec)
    trans_mat = plane_transformation(rec_lat_vec, plot_vec)
    spin_text = spin_text_data[band_index]
    k_point = [trans_mat.dot(np.array(k, dtype=np.float)) for k in spin_text["k-point"]] 
    spin = [trans_mat.dot(np.array(s, dtype=np.float)) for s in spin_text["spin"]]
    energy = np.array(spin_text["energy"], dtype=np.float)
    bohr_to_angstrom = 0.529177

    k_x = np.array([k[0] / bohr_to_angstrom for k in k_point]).reshape(grid)
    k_y = np.array([k[1] / bohr_to_angstrom for k in k_point]).reshape(grid)
    k_z = np.array([k[2] / bohr_to_angstrom for k in k_point]).reshape(grid)

    s_x = np.array([s[0] for s in spin]).reshape(grid)
    s_y = np.array([s[1] for s in spin]).reshape(grid)
    s_z = np.array([s[2] for s in spin])

    mpl.rcParams['grid.linewidth']  = 3
    mpl.rcParams['xtick.labelsize'] = 25
    mpl.rcParams['ytick.labelsize'] = 25
    mpl.rcParams['xtick.major.width'] = 2
    mpl.rcParams['ytick.major.width'] = 2
    mpl.rcParams['xtick.major.size'] = 5
    mpl.rcParams['ytick.major.size'] = 5
    mpl.rcParams['axes.edgecolor']  = 'black'
    mpl.rcParams['axes.labelsize']  = 33      # fontsize of the x any y labels
    mpl.rcParams['axes.labelcolor'] = 'black'
    mpl.rcParams['axes.linewidth']  = 2.0     # set the value globally
    plt.rcParams['xtick.major.pad'] = 10
    plt.rcParams['ytick.major.pad'] = 10
    plt.rcParams.update({'mathtext.default':'regular'})
    
    fig = plt.figure(figsize=(12,10),dpi=300)
    gs = gridspec.GridSpec(1, 1, figure=fig)
    ax = fig.add_subplot(gs[0])
    
    # plot colour map
    if (contour == "spin_z"):
        if (contour_threshhold == "max"):
            thr = round(max(np.abs(s_z)), 2)
        else:
            thr = float(contour_threshhold)
            
        if thr==0.0:
            raise ValueError('Threshhold of contour plot for s_z is zero. You can set it manually with the -cthr option.')
        
        ticks = np.linspace(-1, 1, 5) * thr
        s_z = s_z.reshape(grid)
        
        cp = ax.contourf(k_x, k_y, s_z, 100, cmap="bwr", vmin= -thr, vmax = thr)
        
        cbar=fig.colorbar(cp, orientation="vertical", ticks=ticks)
        cbar.set_label('$s_z$')
        
    elif (contour == "energy"):
        thr_min, thr_max = round(min(energy), 1), round(max(energy),1)
        ticks = thr_min + np.linspace(0.0, 1.0, 8) * (thr_max - thr_min)
        energy = np.array(energy).reshape(grid)
        cp = ax.contourf(k_x, k_y, energy, 100, cmap="autumn", vmin= thr_min, vmax = thr_max)
        cbar=fig.colorbar(cp, orientation="vertical", ticks=ticks)
        cbar.set_label('E [eV]', labelpad=50, rotation=-90)   
             
    # plot spin texture
    ax.quiver(k_x, k_y, s_x, s_y)
    ax.set_xlabel("$k_1$ [a. u.]")
    ax.locator_params(nbins=3, axis='x')
    ax.set_ylabel("$k_2$ [a. u.]")
    ax.locator_params(nbins=3, axis='y')
    

    plt.savefig("PLOT-spintext.png",  bbox_inches='tight')

if __name__ == "__main__":
    input_options = option_parser()
    main(**input_options)
