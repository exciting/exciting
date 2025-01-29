#!/usr/bin/env python
# -*- coding: utf-8 -*-

import matplotlib
matplotlib.use('Agg')  # Use 'Agg' backend for non-interactive plotting

import numpy as np
import matplotlib.pyplot as plt
import argparse
import sys

# Define plotting styles globally
line_styles = ['-', '--', ':', '-', '--', ':', '-.', '-']  # Ensure 8 styles
colors = [
    '#1f77b4',  # Dataset 1: Blue
    '#ff7f0e',  # Dataset 2: Orange
    '#006400',  # Dataset 3: Green
    '#ff0000',  # Dataset 4: Black (강렬한 색상으로 변경)
    '#9467bd',  # Dataset 5: Purple
    '#8c564b',  # Dataset 6: Brown
    '#e377c2',  # Dataset 7: Pink
    '#7f7f7f'   # Dataset 8: Gray
]
markers = ['o', 's', '^', 'D', 'v', '<', '>', 'X']  # 다양한 마커

def read_epsilon_file(filename):
    """
    Reads an EPSILON file and extracts frequencies and Im(epsilon_M) values.

    Parameters:
        filename (str): Path to the EPSILON file.

    Returns:
        tuple: Two NumPy arrays containing frequencies and Im(epsilon_M) values.
    """
    frequencies = []
    im_epsm = []
    try:
        with open(filename, 'r') as f:
            for line in f:
                # Skip empty lines and comments
                if line.strip() == '' or line.strip().startswith('#'):
                    continue
                tokens = line.strip().split()
                if len(tokens) < 3:
                    continue
                try:
                    freq = float(tokens[0])
                    im_eps = float(tokens[2])
                    frequencies.append(freq)
                    im_epsm.append(im_eps)
                except ValueError:
                    # Skip lines with non-numeric data
                    continue
    except IOError:
        print("Error: File '{}' not found or unreadable.".format(filename))
        sys.exit(1)
    return np.array(frequencies), np.array(im_epsm)

def plot_spectra(filenames, labels=None, x_range=None, title_prefix=None):
    """
    Plots the spectra from multiple EPSILON files.

    Parameters:
        filenames (list): List of EPSILON filenames to plot.
        labels (list, optional): List of labels for each file.
        x_range (tuple, optional): Tuple specifying the x-axis range (min, max).
        title_prefix (str, optional): Prefix to add to the plot title.
    """
    plt.figure(figsize=(14, 8))  # Increased width to accommodate legend outside
    for i, filename in enumerate(filenames):
        freq, im_epsm = read_epsilon_file(filename)
        label = labels[i] if labels else filename
        plt.plot(freq, im_epsm, label=label,
                 linestyle=line_styles[i % len(line_styles)],
                 color=colors[i % len(colors)],
                 linewidth=2.5)  # Adjusted linewidth
    plt.xlabel('Energy [eV]', fontsize=22, fontweight='bold')
    plt.ylabel(r'Im $\epsilon_{\mathrm{M}}$', fontsize=22, fontweight='bold', labelpad=15)  # Added labelpad=15
    # Construct the title with optional prefix
    if title_prefix:
        title = '{} '.format(title_prefix)
    else:
        title = 'Macroscopic Dielectric Function'
    plt.title(title, fontsize=26, fontweight='bold', y=1.02)  # Adjusted y
    # Place legend outside the plot
    plt.legend(loc='upper left', bbox_to_anchor=(1, 1), fontsize=16)
    plt.grid(True)
    if x_range:
        plt.xlim(x_range)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.tight_layout()
    plt.subplots_adjust(right=0.75)  # Adjust right to make space for legend
    plt.savefig('spectra.png', dpi=300)
    plt.close()

def plot_difference(filenames, labels=None, x_range=None, title_prefix=None):
    """
    Plots the absolute differences in Im(epsilon_M) with respect to the first EPSILON file.

    Parameters:
        filenames (list): List of EPSILON filenames to plot differences.
        labels (list, optional): List of labels for each file.
        x_range (tuple, optional): Tuple specifying the x-axis range (min, max).
        title_prefix (str, optional): Prefix to add to the plot title.
    """
    plt.figure(figsize=(14, 8))  # Increased width to accommodate legend outside
    freq_list = []
    im_epsm_list = []
    for filename in filenames:
        freq, im_epsm = read_epsilon_file(filename)
        freq_list.append(freq)
        im_epsm_list.append(im_epsm)
    # Use the first file as the base
    base_freq = freq_list[0]
    base_im_epsm = im_epsm_list[0]
    reference_label = labels[0] if labels else filenames[0]
    for i in range(1, len(im_epsm_list)):
        # Check if frequencies match, if not, interpolate
        if not np.array_equal(freq_list[i], base_freq):
            print("Frequencies of file '{}' do not match the base file. Interpolating.".format(filenames[i]))
            interp_im_epsm = np.interp(base_freq, freq_list[i], im_epsm_list[i])
            diff = np.abs(interp_im_epsm - base_im_epsm)
        else:
            diff = np.abs(im_epsm_list[i] - base_im_epsm)
        # Set label as |base - current|
        current_label = labels[i] if labels else filenames[i]
        label = '|{} - {}|'.format(reference_label, current_label)
        plt.plot(base_freq, diff, label=label,
                 linestyle=line_styles[(i-1) % len(line_styles)],
                 color=colors[(i-1) % len(colors)],
                 linewidth=2.5)  # Adjusted linewidth
    # Construct the title with optional prefix
    if title_prefix:
        title = '{} Absolute Difference (with respect to {})'.format(title_prefix, reference_label)
    else:
        title = 'Absolute Difference (with respect to {})'.format(reference_label)
    plt.xlabel('Energy [eV]', fontsize=22, fontweight='bold')
    plt.ylabel(r'Absolute Difference in Im $\epsilon_{\mathrm{M}}$', fontsize=22, fontweight='bold', labelpad=15)  # Added labelpad=15
    plt.title(title, fontsize=18, fontweight='bold', y=1.02)  # Adjusted y and fontsize
    # Place legend outside the plot
    plt.legend(loc='upper left', bbox_to_anchor=(1, 1), fontsize=16)
    plt.grid(True)
    if x_range:
        plt.xlim(x_range)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.tight_layout()
    plt.subplots_adjust(right=0.75)  # Adjust right to make space for legend
    plt.savefig('difference.png', dpi=300)
    plt.close()

def main():
    parser = argparse.ArgumentParser(description='Plot EPSILON files and their absolute differences.')
    parser.add_argument('filenames', nargs='+', help='EPSILON filenames to plot (up to 8 files)')
    parser.add_argument('-lx', '--label_x', default='Energy [eV]', help='Label for x-axis')
    parser.add_argument('-ly', '--label_y', default='Im $\epsilon_{\mathrm{M}}$', help='Label for y-axis')
    parser.add_argument('-t', '--title', default='', help='Prefix for the plot title')
    parser.add_argument('-x', '--x_range', nargs=2, type=float, help='Range for x-axis (min max)')
    parser.add_argument('-l', '--labels', nargs='+', help='Labels for each file')
    parser.add_argument('-d', '--difference', action='store_true', help='Plot absolute difference between spectra')
    args = parser.parse_args()

    max_files = len(colors)
    if len(args.filenames) > max_files:
        print("Maximum of {} files can be plotted.".format(max_files))
        sys.exit(1)

    if args.labels:
        if len(args.labels) != len(args.filenames):
            print("Number of labels must match number of filenames.")
            sys.exit(1)

    x_range = tuple(args.x_range) if args.x_range else None

    if args.difference and len(args.filenames) < 2:
        print("At least two filenames are required to plot differences.")
        sys.exit(1)

    title_prefix = args.title.strip() if args.title else None

    if args.difference:
        plot_difference(args.filenames, labels=args.labels, x_range=x_range, title_prefix=title_prefix)
    else:
        plot_spectra(args.filenames, labels=args.labels, x_range=x_range, title_prefix=title_prefix)

if __name__ == '__main__':
    main()

