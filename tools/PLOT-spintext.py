#!/usr/bin/python3

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

from excitingtools.parser.propertiesParser import parse_spintext
from excitingtools.lattice import parse_lattice_vectors, reciprocal_lattice_vectors, plane_transformation

def parse_spintext_grid(file_path:str) -> np.array:
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

def parse_spintext_plane(file_path:str) -> np.array:
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



    fig = plt.figure(figsize=(12,10),dpi=300)
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
    gs = gridspec.GridSpec(1, 1, figure=fig)
    

    # plot the colour map
    ax = fig.add_subplot(gs[0])
    if (contour == "spin_z"):
        if (contour_threshhold == "max"):
            thr = max(np.abs(s_z))
        else:
            thr = float(contour_threshhold)
        s_z = s_z.reshape(grid)
        cp = ax.contourf(k_x, k_y, s_z, 100, cmap="bwr", vmin= -thr, vmax = thr)
        cbar=fig.colorbar(cp, orientation="vertical", ticks=[-thr, -thr/2,0, thr/2, thr])
        cbar.set_label('$s_z$ [b$^{-1}$]')
    elif (contour == "energy"):
        thr_min, thr_max = min(energy), max(energy)
        energy = np.array(energy).reshape(grid)
        cp = ax.contourf(k_x, k_y, energy, 100, cmap="autumn", vmin= thr_min, vmax = thr_max)
        cbar=fig.colorbar(cp, orientation="vertical", ticks=[thr_min, thr_min+(thr_max-thr_min)/2, thr_max])
        cbar.set_label('E [eV]')
    ax.quiver(k_x, k_y, s_x, s_y)
    ax.set_xlabel("$k_1$ [Å $^{-1}$]")
    ax.set_ylabel("$k_2$ [Å $^{-1}$]")

    plt.savefig("PLOT-spintext.png",  bbox_inches='tight')
    plt.show()

if __name__ == "__main__":
    input_options = option_parser()
    main(**input_options)
