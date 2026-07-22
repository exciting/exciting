""" Parse and plot electron-phonon renormalized band gap for different temperatures.

Requires the following files:
    - EVALQP_T####.dat

Located at `excitingscripts/plot/eph_gap_renormalization.py`.

Call as:

```bash
python3 -m excitingscripts.plot.eph_gap_renormalization
```
"""

import argparse as ap
import numpy as np
import matplotlib.pyplot as plt
from excitingtools.exciting_dict_parsers.eph_parser import parse_bandgap_renormalization
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
    # directory
    help_directory = 'Directories in which the data to be plotted have to be found. If only one or no directory is specified, the data for the plots are taken from the same directory. Default value is the current directory.'
    p.add_argument('-d', '--directory', nargs='*', default=["."], type=str, help=help_directory)
    # file prefix
    help_prefix = 'List of allowed prefixes with which the quasi-particle energy files might start. Default is "EVALQP".'
    p.add_argument('--prefix', nargs='*', default=["EVALQP"], type=str, help=help_prefix)
    # VBM vector
    help_vvb = 'Vector (in lattice coordinates) of point to take for valence band. If no vector is specified, all valence band is reported for all vectors found in the directories.'
    p.add_argument('-vvb', '--vector-valence-band', nargs=3, default=None, type=float, help=help_vvb)
    # CBM vector
    help_vcb = 'Vector (in lattice coordinates) of point to take for conduction band. If no vector is specified, all conduction band is reported for all vectors found in the directories.'
    p.add_argument('-vcb', '--vector-conduction-band', nargs=3, default=None, type=float, help=help_vcb)
    # plot
    help_plot = 'Plot band gap (renormalization) vs temperature.'
    p.add_argument('-p', '--plot', action='store_true', help=help_plot)
    # renormalization
    help_renorm = 'Plot energy renormalization instead of absolute energy.'
    p.add_argument('--renorm', action='store_true', help=help_renorm)

    # parse arguments
    args = p.parse_args()
    input_options = {}
    input_options['directory'] = args.directory
    input_options['prefix'] = args.prefix
    input_options['vvb'] = args.vector_valence_band
    input_options['vcb'] = args.vector_conduction_band
    input_options['plot'] = args.plot
    input_options['renorm'] = args.renorm

    return input_options
#_______________________________________________________________________________

def main(input_options) -> None:
    '''
    input:
    :input_options: dictionary that holds the input options parsed from the command line arguments
    '''
    # set input variables
    directory = input_options['directory']
    prefix = input_options['prefix']
    key_vb = ', '.join([f'{x:20.12g}' for x in input_options['vvb']]) if input_options['vvb'] else None
    key_cb = ', '.join([f'{x:20.12g}' for x in input_options['vcb']]) if input_options['vcb'] else None
    renorm = input_options['renorm']

    # get parsed data
    data = parse_bandgap_renormalization(directory, prefixes=prefix)

    # select points for valence and conduction band
    ipvb = [i for i, point in enumerate(data['points']) if not key_vb or ', '.join([f'{x:20.12g}' for x in point['vector']]) == key_vb]
    ipvb = ipvb * len(data['points']) if key_vb and not key_cb else ipvb
    if len(ipvb) == 0:
        raise ValueError('No data could be found for given valence band vector.')
    ipcb = [i for i, point in enumerate(data['points']) if not key_cb or ', '.join([f'{x:20.12g}' for x in point['vector']]) == key_cb]
    ipcb = ipcb * len(data['points']) if key_cb and not key_vb else ipcb
    if len(ipcb) == 0:
        raise ValueError('No data could be found for given conduction band vector.')

    # select data
    selected_data = []
    for ipv, ipc in zip(ipvb, ipcb):
        pvb = data['points'][ipv]
        pcb = data['points'][ipc]
        selected_data.append({ \
                'valence band vector': pvb['vector'], \
                'conduction band vector': pcb['vector'] \
                })
        selected_data[-1]['valence band index'] = pvb['valence band index'] if 'valence band index' in pvb else None
        selected_data[-1]['conduction band index'] = pcb['conduction band index'] if 'conduction band index' in pcb else None
        selected_data[-1]['valence band KS energy'] = pvb['valence band KS energy'] * Hartree_to_eV if 'valence band KS energy' in pvb else np.nan
        selected_data[-1]['conduction band KS energy'] = pcb['conduction band KS energy'] * Hartree_to_eV if 'conduction band KS energy' in pcb else np.nan
        selected_data[-1]['temperatures'] = sorted(list(set(pvb['temperatures'].keys()) | set(pcb['temperatures'].keys())))
        selected_data[-1]['valence band QP energy'] = []
        selected_data[-1]['conduction band QP energy'] = []
        for temp in selected_data[-1]['temperatures']:
            vb = pvb['temperatures'][temp] if temp in pvb['temperatures'].keys() else {}
            cb = pcb['temperatures'][temp] if temp in pcb['temperatures'].keys() else {}
            selected_data[-1]['valence band QP energy'].append(vb['valence band QP energy'].real * Hartree_to_eV if 'valence band QP energy' in vb else np.nan)
            selected_data[-1]['conduction band QP energy'].append(cb['conduction band QP energy'].real * Hartree_to_eV if 'conduction band QP energy' in cb else np.nan)

    # pretty print selected data
    for data in selected_data:
        ivb = f'n = {data['valence band index']:d}' if data['valence band index'] else 'not found'
        icb = f'n = {data['conduction band index']:d}' if data['conduction band index'] else 'not found'
        vvb = f'[{data['valence band vector'][0]:.3f}, {data['valence band vector'][1]:.3f}, {data['valence band vector'][2]:.3f}]'
        vcb = f'[{data['conduction band vector'][0]:.3f}, {data['conduction band vector'][1]:.3f}, {data['conduction band vector'][2]:.3f}]'
        print('-' * 89)
        print(f'{"":5s} | ' \
              f'{f"VALENCE BAND ({ivb})":^25s} | ' \
              f'{f"CONDUCTION BAND ({icb})":^25s} | ' \
              f'{"BAND GAP":^25s}')
        print(f'{"":^5s} | ' \
              f'{f"{vvb}":^25s} | ' \
              f'{f"{vcb}":^25s} | ' \
              f'{"":^25s}')
        print(f'{"":^5s} | ' \
              f'{f"E_0 = {data['valence band KS energy']:.3f} eV":^25s} | ' \
              f'{f"E_0 = {data['conduction band KS energy']:.3f} eV":^25s} | ' \
              f'{f"E_g_0 = {data['conduction band KS energy'] - data['valence band KS energy']:.3f} eV":^25s}')
        print(f'{"T [K]":>5s} | ' \
              f'{"E_EPH [eV]":>12s} {"dE [meV]":>12s} | ' \
              f'{"E_EPH [eV]":>12s} {"dE [meV]":>12s} | ' \
              f'{"E_g_EPH [eV]":>12s} {"dE_g [meV]":>12s}')
        print('-' * 89)
        for temp, evb_eph, ecb_eph in zip(data['temperatures'], data['valence band QP energy'], data['conduction band QP energy']):
            print(f'{temp:5d} | ' \
                  f'{evb_eph:12.3f} {(evb_eph - data['valence band KS energy']) * 1e3:+12.1f} | ' \
                  f'{ecb_eph:12.3f} {(ecb_eph - data['conduction band KS energy']) * 1e3:+12.1f} | ' \
                  f'{ecb_eph - evb_eph:12.3f} {((ecb_eph - data['conduction band KS energy'])- (evb_eph - data['valence band KS energy'])) * 1e3:+12.1f}')
        print('-' * 89)
        print()

    # plot selected data
    if input_options['plot']:
        plt.rcParams.update({
            'xtick.major.width': 2,
            'ytick.major.width': 2,
            'xtick.labelsize': 20,
            'ytick.labelsize': 20,
            'axes.linewidth': 2,
            'lines.linewidth': 3,
            'axes.labelsize': 20,
            'axes.grid': True,
        })

        fig, axs = plt.subplots(figsize=(16,9), ncols=2, nrows=2, sharex=True, layout='constrained')
        ax_vb = axs[1,1]
        ax_cb = axs[0,1]
        gs = axs[0,0].get_gridspec()
        for ax in axs[0:,0]:
          ax.remove()
        ax_gap = fig.add_subplot(gs[0:,0])

        prefix = r'$\Delta E$' if renorm else r'$E$'
        unit = 'meV' if renorm else 'eV'
        scale = 1e3 if renorm else 1
        ax_gap.set(xlabel='Temperature [K]', ylabel=f'{prefix} Gap [{unit}]')
        ax_vb.set(xlabel='Temperature [K]', ylabel=f'{prefix} valence [{unit}]')
        ax_cb.set(ylabel=f'{prefix} conduction [{unit}]')

        for data in selected_data:
            if len(data['temperatures']) <= 1:
                continue
            temps = np.array(data['temperatures'])
            evb = np.array(data['valence band QP energy']) * scale
            ecb = np.array(data['conduction band QP energy']) * scale
            if renorm:
                evb -= data['valence band KS energy'] * scale
                ecb -= data['conduction band KS energy'] * scale
            ax_gap.plot(temps, ecb-evb) 
            ax_vb.plot(temps, evb)
            ax_cb.plot(temps, ecb)

        ylimgap = ax_gap.get_ylim()
        ylimvb = ax_vb.get_ylim()
        ylimcb = ax_cb.get_ylim()
        if renorm:
            ylimgap = (min(ylimgap[0], 0), max(ylimgap[1], 0))
            ax_gap.set(ylim=ylimgap)
            ylimvb = (min(ylimvb[0], 0), max(ylimvb[1], 0))
            ylimcb = (min(ylimcb[0], 0), max(ylimcb[1], 0))
            span = max(ylimvb[1]-ylimvb[0], ylimcb[1]-ylimcb[0])
            ax_vb.set(ylim=(ylimvb[0], ylimvb[0]+span))
            ax_cb.set(ylim=(ylimcb[1]-span, ylimcb[1]))
        else:
            span = max(ylimvb[1]-ylimvb[0], ylimcb[1]-ylimcb[0])
            ax_vb.set(ylim=((np.sum(ylimvb)-span)/2, (np.sum(ylimvb)+span)/2))
            ax_cb.set(ylim=((np.sum(ylimcb)-span)/2, (np.sum(ylimcb)+span)/2))

        fig.savefig('PLOT.png', format='png', dpi=300, bbox_inches='tight')

if __name__ == "__main__":
    input_options = option_parser()
    main(input_options)
