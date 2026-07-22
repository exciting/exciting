""" Plot (renormalized) electron / phonon band structure and spectral functions.

Requires the following files:
    - eph_el_disp.json
    - eph_ph_disp.json
    - eph_el_qp_disp.json
    - eph_else+sfun_P###_T####.json

Located at `excitingscripts/plot/eph_dispersion.py`.

Call as:

```bash
python3 -m excitingscripts.plot.eph_dispersion
```
"""

import argparse as ap
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import json
import os
import re
from scipy.interpolate import pchip_interpolate as pchip
from excitingtools.exciting_dict_parsers.eph_parser import parse_bandgap_renormalization
from excitingtools.dataclasses import BandData
from excitingtools.constants.units import Hartree_to_eV

def option_parser():
    """
    Parse command line inputs 

    Parse: 
        directory
        prefix
        vector-valence-band
        vector-conduction-band

    :return input_options: Dictionary of parsed command line arguments 
    """
    # define valid arguments
    p = ap.ArgumentParser(description=\
            'Parse and plot electron-phonon renormalized band gap for different temperatures.')
    # dispersion
    help_disp = 'JSON file containing the electron / phonon dispersion, e.g. "eph_el_disp.json" or "eph_ph_disp.json".'
    p.add_argument('-disp', '--dispersion', default=None, type=str, help=help_disp)
    # QP energies
    help_qp = 'JSON file containing the quasi-particle dispersion, e.g. "eph_el_qp_disp.json".'
    p.add_argument('-qp', '--quasi-particles', default=None, type=str, help=help_qp)
    # spectral function
    help_sfun = 'List of JSON files containing self-energy and spectral function, e.g. "eph_else+sfun_P*_T0300.json" for all files with spectral function at 300 K.'
    p.add_argument('-sfun', '--spectral-function', nargs="*", default=None, type=str, help=help_sfun)
    # unit
    help_unit = 'Energy unit for output. Choose from "eV" (electronvolt, default), "meV" (millielectronvolt), "icm" (inverse centimeter), "THz" (terrahertz).'
    p.add_argument('-u', '--unit', default="eV", type=str, help=help_unit)
    # energy range
    help_erange = 'Lower and upper energy boundary to display. Use "nan" for automatic determination.'
    p.add_argument('-e', '--energy-range', nargs=2, default=[None, None], type=float, help=help_erange)
    # line color
    help_lc = 'Pair of predefined named colors. The first color applies to the dispersion specified via -disp, the second to the quasi-particle dispersion specified via -qp. Default is "teal tomato".'
    p.add_argument('-lc', '--linecolor', nargs=2, default=['teal', 'tomato'], type=str, help=help_lc)
    # maximal spectral function
    help_smax = 'Upper limit for spectral function plot. If not specified, a reasonable default value depending on the specified unit applies.'
    p.add_argument('-smax', '--spectral-function-maximum', default=None, type=float, help=help_smax)
    # number of frequencies
    help_nfreq = 'Number of frequency points to sample spectral function at. Default is 500. Might to be increased, if large energy range is displayed.'
    p.add_argument('-nf', '--num-frequencies', default=500, type=int, help=help_nfreq)
    # colormap
    help_cmap = 'Name of one of the predefined matplotlib colormaps. Default is "terrain".'
    p.add_argument('-cmap', '--colormap', default='terrain', type=str, help=help_cmap)

    # parse arguments
    args = p.parse_args()
    input_options = {}
    input_options['dispersion'] = args.dispersion
    input_options['qp'] = args.quasi_particles
    input_options['sfun'] = args.spectral_function
    input_options['unit'] = args.unit
    input_options['erange'] = args.energy_range
    input_options['lc'] = args.linecolor
    input_options['smax'] = args.spectral_function_maximum
    input_options['nfreq'] = args.num_frequencies
    input_options['cmap'] = args.colormap

    return input_options
#_______________________________________________________________________________

def refine_sfun(xc, yc):
    ny, nxc, nz = yc.shape
    nxf = 2*nxc - 1

    xf = np.zeros(nxf)
    xf[::2] = xc
    xf[1::2] = (xc[:-1] + xc[1:]) / 2

    yf = np.zeros((ny, nxf, nz))
    for iz in range(nz):
        y1 = yc[:, 0, iz]
        yf[:, 0, iz] = y1
        for ixc in range(1, nxc):
            ixf = 2*ixc
            y2 = yc[:, ixc, iz]
            c = np.convolve(y1, np.flip(y2), mode='full')
            d = ny - np.argmax(c)
            s1 = int(d/2)
            s2 = s1 - d
            yf[:, ixf-1, iz] = (np.roll(y1, s1) + np.roll(y2, s2) + np.roll(y1, -s2) + np.roll(y2, -s1)) / 4
            yf[:, ixf, iz] = y2
            y1 = y2

    return xf, yf
#_______________________________________________________________________________

def main(input_options) -> None:
    '''
    input:
    :input_options: dictionary that holds the input options parsed from the command line arguments
    '''
    Ha2eV = Hartree_to_eV
    Ha2icm = Ha2eV * 1e3 * 8.06573
    Ha2THz = Ha2icm * 0.0299793

    # set input variables
    match input_options['unit']:
        case 'Ha':
            scale = 1
            unit = 'Ha'
        case 'meV':
            scale = Ha2eV * 1e3
            unit = 'meV'
        case 'icm':
            scale = Ha2icm
            unit = r'cm$^-1$'
        case 'THz':
            scale = Ha2THz
            unit = 'THz'
        case _:
            scale = Ha2eV
            unit = 'eV'
    disp = input_options['dispersion']
    qp = input_options['qp']
    sfun_files = input_options['sfun']
    erange = input_options['erange']
    lc1, lc2 = input_options['lc']
    smax = input_options['smax'] if input_options['smax'] else 200 / scale
    nfreq = input_options['nfreq']
    cmap = input_options['cmap']

    # parse band structure
    path = None
    bands = None
    qpbands = None
    pdist = None
    if disp:
        with open(disp) as file:
            data = json.load(file)
            path = data['path']
            bands = np.array(data['bands']) * scale
            erange[0] = np.min(bands) if not erange[0] or np.isnan(erange[0]) else erange[0]
            erange[1] = np.max(bands) if not erange[1] or np.isnan(erange[1]) else erange[1]
    if qp:
        with open(qp) as file:
            data = json.load(file)
            path = data['path']
            qpbands = np.array(data['bands']) * scale
            qplinewidth = np.array(data['imaginary part']) * scale
            erange[0] = np.min(qpbands) if not erange[0] or np.isnan(erange[0]) else erange[0]
            erange[1] = np.max(qpbands) if not erange[1] or np.isnan(erange[1]) else erange[1]
    if sfun_files:
        with open(sfun_files[0]) as file:
            data = json.load(file)
            nband = len(data['frequencies'])
            erange[0] = np.min(data['frequencies'][0]) * scale if not erange[0] or np.isnan(erange[0]) else erange[0]
            erange[1] = np.max(data['frequencies'][0]) * scale if not erange[1] or np.isnan(erange[1]) else erange[1]
    if path:
        points = np.array([p['coord_lat'] for p in path['points']])
        pdist = np.array([p['distance'] for p in path['points']])
        vertices = [{'distance': v['distance'], 'label': v['label'], 'coord': v['coord_lat']} for v in path['vertices']]
        bands = BandData(bands=bands, k_points=points, e_fermi=0, flattened_k_points=pdist, vertices=vertices) if type(bands) is np.ndarray else None
        qpbands = BandData(bands=qpbands, k_points=points, e_fermi=0, flattened_k_points=pdist, vertices=vertices) if type(qpbands) is np.ndarray else None

    # parse and interpolate spectral function
    if sfun_files:
        idx = [re.search(r'P\d{3}', f) for f in sfun_files]
        idx = [int(i[0][1:]) for i in idx if i]
        idx = np.argsort(idx)
        sfun_files = [sfun_files[i] for i in idx]
        npt = len(idx)
        pdist_sfun = pdist if type(pdist) is np.ndarray else np.linspace(0, 1, npt)
        freqs = np.linspace(1.1*erange[0]-0.1*erange[1], 1.1*erange[1]-0.1*erange[0], nfreq+1)
        freqs = (freqs[:nfreq] + freqs[1:])/2
        sfun = np.zeros((nfreq, npt, nband))
        for ip, sfun_file in enumerate(sfun_files):
            with open(sfun_file) as file:
                data = json.load(file)
                for ib in range(nband):
                    e0 = data['electron energy'][ib] * scale
                    freqs0 = np.array(data['frequencies'][ib]) * scale
                    selfen = np.squeeze(np.array(data['self-energy'][ib]).view(np.complex128), axis=-1) * scale
                    selfen = pchip(freqs0, np.real(selfen), freqs) + 1j * pchip(freqs0, np.imag(selfen), freqs)
                    sfun[:,ip,ib] = -1/np.pi * np.imag(selfen) / ((freqs - (e0 + np.real(selfen)))**2 + np.imag(selfen)**2)
        depth = int(np.round(max(0, np.log(nfreq/npt * 16/9) / np.log(2))))
        for i in range(depth):
            pdist_sfun, sfun = refine_sfun(pdist_sfun, sfun)
        sfun = np.sum(sfun, axis=2)

    plt.rcParams.update({
        'xtick.major.width': 2,
        'ytick.major.width': 2,
        'xtick.labelsize': 30,
        'ytick.labelsize': 30,
        'axes.linewidth': 2,
        'lines.linewidth': 3,
        'axes.labelsize': 30,
    })

    fig, ax = plt.subplots(figsize=(16,9), layout='constrained')
    if disp:
        vertices, labels = bands.band_path()
    if qp:
        vertices, labels = qpbands.band_path()
    if qp or disp:
        ax.set(xticks=vertices, xticklabels=labels, xlim=(vertices[0], vertices[-1]))
        for x in vertices[1:-1]:
            ax.axvline(x, linestyle='--', color='black')
    ax.set(ylabel=f'Energy [{unit}]', ylim=erange)
    ax.tick_params(axis='both', which='major')
    if disp:
        for ib in range(0, bands.n_bands):
            ax.plot(bands.flattened_k_points, bands.bands[:, ib], color=lc1)
    if qp:
        for ib in range(0, qpbands.n_bands):
            ax.plot(qpbands.flattened_k_points, qpbands.bands[:, ib], color=lc2, ls='--')
            ax.fill_between(qpbands.flattened_k_points, qpbands.bands[:,ib]-qplinewidth[:,ib], qpbands.bands[:,ib]+qplinewidth[:,ib], color=lc2, alpha=0.3), 

    if sfun_files:
        norm = mpl.colors.Normalize(0, smax, clip=True)
        levels = smax*np.linspace(0, 1, 200)**3
        im = ax.contourf(pdist_sfun, freqs, sfun, levels=levels, cmap=cmap, norm=norm, alpha=1, zorder=-1000, extend='max')
        fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=plt.get_cmap(name=cmap)), ax=ax, pad=0.0, label=r'Spectral function ['+unit+r'$^{-1}$]')

    fig.savefig('PLOT.png', format='png', dpi=100, bbox_inches='tight')

if __name__ == "__main__":
    input_options = option_parser()
    main(input_options)
