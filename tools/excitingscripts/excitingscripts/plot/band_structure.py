"""Plot single and multiple electronic and phonon band structures (BS)

Require the following files:
- for electronic BS: input.xml, BAND.OUT, BAND-QP.OUT, BAND_WANNIER.OUT, BANDLINES.OUT
- for phonon BS: input.xml, PHDISP.OUT, PHLINES.OUT

More details can be found **[here](https://www.exciting-code.org/home/the-python-script-plot.band_structure)**.

Located at `excitingscripts/plot/band_structure.py`.

Call as:

```bash
python3 -m excitingscripts.plot.band_structure
```
"""

import argparse as ap
import os
import sys
from xml.etree import ElementTree as ET

import matplotlib.pyplot as plt
import matplotlib.style
import matplotlib.ticker as ticker
import numpy as np

if matplotlib.__version__.split(".")[0] == "2": matplotlib.style.use('classic')

#START_DEF++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def option_parser():
    """
    Parse command line inputs 

    Parse: 
        directory
        eboundary
        assign_type
        phonon
        no_legend
        ezero
        scale_box
        kboundary
        kpoint_boundary
        eunit
        funit
        title
        no_title
        legend_position
        scale_box
        invert_colors
        invert_plots
        max_ticks_y
        legend_label
        show

    :return input_options: Dictionary of parsed command line arguments 
    """
    p = ap.ArgumentParser(description=\
                'Plot single and multiple electronic/phonon band structures.')
    
    help_directory = 'List of the directories in which the data to be plotted have to be found.    If only one or no directory is specified, the data for the plots are taken from the same directory. Default value is the current directory.'
    
    help_eboundary = 'One or two floats corresponding to the minimum and maximum energy in the plot of the electronic band-structure in the units specified by --eunit. If the argument --phonon is present, the option values specify the minimum and maximum frequency in the phonon-dispersion plot in the units specified by --funit. If either the argument --eboundary is not present or no floats are given, the energy (frequency) boundaries are chosen correspondingly to the maximum and minimum of the full data to be plotted.'
    
    help_assign_type = "List of the description keys for each plot of the electronic band-structure. Possible choices are 'KS' (standard Kohn-Sham calculation), 'GW' (G0W0 calculation), and 'WA' (data are interpolated by using Wannier functions). If not present, the option 'KS 'is assumed for all plots. Not used if the argument --phonon is present."
    
    help_ezero = "Set the energy zero in the plot to either 'Ef' (Fermi energy, default) or 'vbM' (valence-band maximum)."
    
    help_eunit = "Set the units of the energy appearing in the plot of the electronic band-structure. Possible choices are 'eV' (electronvolt, default) and 'Ha' (Hartree)."

    help_funit = "Set the units of the frequency appearing in the plot of the phonon band-structure. Possible choices are 'icm' (inverse centimeter, cm^-1, default), 'meV' (millielectronvolt), and 'THz' (terahertz)."

    help_legend_position = "The location of the legend. The strings 'upper left', 'upper right', 'lower left', 'lower right' place the legend at the corresponding corner of the axes/figure. The strings 'upper center', 'lower center', 'center left', 'center right' place the legend at the center of the corresponding edge of the axes/figure. The string 'center' places the legend at the center of the axes/figure. The string 'best' places the legend at the location, among the nine locations defined so far, with the minimum overlap with other drawn artists. This option can be quite slow for plots with large amounts of data; your plotting speed may benefit from providing a specific location. For back-compatibility, 'center right' (but no other location) can also be spelled 'right', and each string locations can also be given as the corresponding numeric value."
                  
    help_phonon = 'If present, it tags the plotting of the phonon-dispersion curves. If absent, the electronic band-structure is plotted.'
    
    help_no_legend = 'If present, it disables the plotting of the legend.'
    
    help_title = "Used as --title 'String as a title' assign a title to the plot."
    
    help_no_title = 'If present, it disables the writing of the title.'
    
    help_scale_box = "One or two floats corresponding to the scaling factor in the horizontal and vertical size of the plot appearence, respectively."
    
    help_kboundary = "One or two floats in the interval [0,1] corresponding to the minimum and maximum relative coordinate along the path in reciprocal space, respectively. The value 0 correspond to the initial point in the path, the value 1 to the last one."
    
    help_kpoint_boundary = "One or two integers indicating the initial and final special k-point in the horizontal path."
    
    help_reverse_colors = "If present, the order of the sequence of colors of the plots is reversed."
    
    help_reverse_plots = "If present, the order of appearance of the plots is reversed."

    help_max_ticks_y = "Specifies the maximum number of ticks along the y-axis in the plot."
    
    help_legend_label = "Specifies the labels to appear in the legend for each plot."

    help_show = "Opens plot in new window"

    #---------------------------------------------------------------------------
    
    p.add_argument('-d','--directory',
                   nargs = '*', default = ["."],
                   type = str, help = help_directory)
    
    p.add_argument('-e','--eboundary',
                   nargs = '*', default = [None, None],
                   type = float, help = help_eboundary)
    
    p.add_argument('-s','--scale_box',
                   nargs = '*', default = [1.0, 1,0],
                   type = float, help = help_scale_box)
    
    p.add_argument('-k','--kboundary',
                   nargs = '*', default = [0.0, 1,0],
                   type = float, help = help_kboundary)
    
    p.add_argument('-kp','--kpoint_boundary',
                   nargs = '*', default = [None, None],
                   type = int, help = help_kpoint_boundary)

    p.add_argument('-a','--assign_type',
                   nargs = '*', default = ['KS'],
                   choices = ['KS', 'GW', 'WA'],
                   type = str, help = help_assign_type)
    
    p.add_argument('-p','--phonon', action='store_true', help = help_phonon)
    
    p.add_argument('-rc','--reverse_colors', action='store_true', help = help_reverse_colors)

    p.add_argument('-rp','--reverse_plots', action='store_true', help = help_reverse_plots)

    p.add_argument('-nl','--no_legend', action='store_true', help = help_no_legend)
     
    p.add_argument('-z','--ezero',
                   type = str, help = help_ezero,
                   choices = ['Ef', 'vbM'], default = 'Ef')

    p.add_argument('-eu','--eunit',
                   type = str, help = help_eunit,
                   choices = ['eV', 'Ha'], default = 'eV')
    
    p.add_argument('-fu','--funit',
                   type = str, help = help_funit,
                   choices = ['icm', 'meV', 'THz'], default = 'icm')
                   
    p.add_argument('-t','--title',
                   type = str, default = None, help = help_title)
    
    p.add_argument('-nt','--no_title', action='store_true', help = help_no_title)
    
    p.add_argument('-mty','--max_ticks_y',
                   type = int, default = None, help = help_max_ticks_y)
    
    p.add_argument('-ll','--legend_label',
                   nargs = '*', default = [],
                   type = str, help = help_legend_label)    

    p.add_argument('-l','--legend_position',
                   type = str, help = help_legend_position,
                   choices = ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10',
                              'best', 'upper right', 'upper left', 'lower left', 'lower right', 'right', 'center left', 'center right', 'lower center', 'upper center', 'center'],
                   default = 'best')

    p.add_argument('-sh', '--show',
                   action='store_true',
                   help=help_show)
                       
    #---------------------------------------------------------------------------
    
    args = p.parse_args()
    input_options = {}
    
    input_options['directory'] = args.directory
    input_options['atype'] = args.assign_type

    input_options['emin'] = None
    if len(args.eboundary) >= 1: input_options['emin'] = args.eboundary[0]
    input_options['emax'] = None
    if len(args.eboundary) >= 2: input_options['emax'] = args.eboundary[1]
    
    input_options['sx'] = 1.0    
    if len(args.scale_box) >= 1: input_options['sx'] = args.scale_box[0]
    input_options['sy'] = 1.0
    if len(args.scale_box) >= 2: input_options['sy'] = args.scale_box[1]
    
    input_options['kmin'] = 0.0    
    if len(args.kboundary) >= 1: input_options['kmin'] = args.kboundary[0]
    input_options['kmax'] = 1.0
    if len(args.kboundary) >= 2: input_options['kmax'] = args.kboundary[1]
    
    input_options['kp'] = args.kpoint_boundary
   
    if ( input_options['kmin'] < 0.0 or
         input_options['kmax'] > 1.0 or
         input_options['kmin']>=input_options['kmax'] ):
        print("\n ERROR: It must be  0 ≤ kmin < kmax ≤ 1 !\n")
        print("        kmin = "+str(input_options['kmin']))
        sys.exit("        kmax = "+str(input_options['kmax'])+"\n")
        
    input_options['phonon'] = args.phonon
    
    input_options['title'] = args.title
    input_options['no_title'] = args.no_title
    
    input_options['ezero'] = args.ezero
    input_options['eunit'] = args.eunit
    input_options['funit'] = args.funit
    
    input_options['no_legend'] = args.no_legend
    input_options['leg_pos'] = args.legend_position
    if len(args.legend_position)<=2:
        input_options['leg_pos'] = int(args.legend_position)
        
    input_options['reverse_plots'] = args.reverse_plots
    input_options['reverse_colors'] = args.reverse_colors

    input_options['maxticksy'] = args.max_ticks_y
    
    input_options['legend'] = args.legend_label

    input_options['show'] = args.show

    return input_options
#_______________________________________________________________________________

def extract_title_text(xmlfile):
    '''
    extract text in title from xmlfile
    '''
    inquire_file(xmlfile)
    text = inquire_element("title",ET.parse(xmlfile).getroot()).text
    return text
#_______________________________________________________________________________

def inquire_spin(xmlfile):
    '''
    check in input.xml if a spin-polarized calculation is performed
    '''
    inquire_file(xmlfile)
    lspin = False
    groundstate = inquire_element("groundstate",ET.parse(xmlfile).getroot())
    for spin in groundstate.findall("spin"): 
        so = str(spin.get("spinorb")).lower()
        if not so == 'true': lspin = True
    return lspin
#_______________________________________________________________________________

def inquire_element(element, xmlpath):
    '''
    check if element is in xmlfile
    '''
    if xmlpath.find(element) is None:
        sys.exit("\n WARNING: Element '"+element+"' not found in "+str(xmlpath)+"\n")
    return xmlpath.find(element)
#_______________________________________________________________________________

def extract_xticks_labels(xmlfile, phonon):
    '''
    extract x-ticks labels for band structure from xmlfile
    '''
    inquire_file(xmlfile)
    bandxml = "bandstructure"
    propxml = "properties"
    if phonon:
        bandxml = "phonondispplot"
        propxml = "phonons"
    properties = inquire_element(propxml,ET.parse(xmlfile).getroot())
    bandstructure = inquire_element(bandxml,properties)
    plot1d = inquire_element("plot1d",bandstructure)
    path = inquire_element("path",plot1d)
    label = []
    buffer_ = None
    for point in path.findall("point"):
        lab = point.attrib["label"]
        if lab.lower() == 'gamma': lab = '\u0393'
        if point.attrib.get("breakafter") == "true":
            buffer_ = lab
            continue
        if buffer_ is not None:
            lab = f"{buffer_},{lab}"
            buffer_ = None
        label.append(lab)
    return label
#_______________________________________________________________________________

def read_xticks_position(infile):
    '''
    read position of ticks on the horizontal axis
    '''
    inquire_file(infile)
    bandlines=list(set(np.genfromtxt(infile)[:,0]))
    bandlines.sort()
    return bandlines
#_______________________________________________________________________________

def set_legend_label(leg_label,spin,leg_spin,local_single,k,j):
    '''
    set legend label for the plot
    '''
    legend_label = leg_label
    if local_single: leg_spin[k] = "spin" + leg_spin[k]
    if spin: legend_label = leg_label + leg_spin[k]
    if j!=0: legend_label = None
    return legend_label
#_______________________________________________________________________________

def set_energy_zero(band, band_spin, spin, ezero):
    '''
    set energy zero, assume in input bands aligned to the Fermi energy 
    '''
    vbm = -1e30  ;  cbm =  1e30
    for j in range(len(band)):
        all_negative = [x for x in band[j] if x < 0]
        if len(all_negative)!=0: vbm = max(vbm, max(all_negative))
        all_positive = [x for x in band[j] if x > 0]
        if len(all_positive)!=0: cbm = min(cbm, min(all_positive))
    eshift = vbm + (cbm-vbm) / 2.
    if ezero== "vbM": eshift = vbm
    band = band-eshift
    if spin: band_spin = band_spin - eshift
    return band, band_spin
#_______________________________________________________________________________
 
def plot_ezero(ezero, ymin, ymax, xmin, xmax, kmin, kmax, sx, phonon, ax1):
    '''
    plot the energy zero
    '''
    ezero_size = "40"
    ezero_label = '$E_{{F}}$'
    xshift = 0.01
    if ezero == "vbM": ezero_label = '$E_{_{VBM}}$'
    dx = xmax-xmin
    dk = kmax-kmin
    ezero_position = xmax+xshift*dx/sx
    if ymin <= 0.0 and 0.0 <= ymax and (not phonon):
        plt.text(ezero_position, 0.0, ezero_label, size=ezero_size, 
                 ha='left', va='center', transform=ax1.transData)
    return
#_______________________________________________________________________________
 
def read_phonon_dispersion(kvec,band,infile,funit):
    '''
    read phonon-dispersion curves frome files infile 
    '''
    au2icm = 2.194746313705e5
    icm2mev = 1.0 / 8.06573
    icm2thz = 0.0299793
    
    conversion_factor = au2icm
    if funit == 'meV': conversion_factor = au2icm * icm2mev
    if funit == 'THz': conversion_factor = au2icm * icm2thz
   
    band.append([])  ;  kvec.append([])
    list1=[]  ;  list2=[]
    
    for line in open(infile):
        i_line=line.split()
        if len(i_line):
            list1.append(float(i_line[0]))
            frequency = (float(i_line[1]))
            if abs(frequency) * au2icm < 1.e-4:
                list2.append(0.0)
            else:
                list2.append(frequency*conversion_factor)
        else:
            kvec[-1].append(list1)
            band[-1].append(list2)
            list1=[]  ;  list2=[]  
    return            
#_______________________________________________________________________________
            
def read_electronic_band(kvec,band,band_spin,spin,infile,eunit):
    '''
    read electronic band-structure frome files infile 
    '''
    ha2ev = 27.211396132 
    
    conversion_factor = ha2ev
    if eunit== "Ha": conversion_factor = 1.0
    
    band.append([])  ;  band_spin.append([])  ;  kvec.append([])
    inquire_file(infile)
    steps=find_steps(infile)
    k_list = np.genfromtxt(infile)[:,0]
    e_list = np.genfromtxt(infile)[:,1] * conversion_factor
    for j in range(0, len(k_list), steps): kvec[-1].append(k_list[j:j + steps])
    for j in range(0, len(e_list) ,steps): band[-1].append(e_list[j:j + steps])
    if spin:
        nbands=int(len(band[-1]) / 2)
        kspin = [] ; band_up = []  ;  band_down = []
        for j in range(nbands): kspin.append(kvec[-1][j])
        for j in range(nbands): band_up.append(band[-1][j])
        for j in range(nbands,2 * nbands): band_down.append(band[-1][j])
        band[-1] = band_up
        band_spin[-1] = band_down
    return
#_______________________________________________________________________________

def find_steps(infile):
    '''
    find number of k values in the electronic energies of a single band
    '''
    inf = open(infile,"r")
    steps = 0
    while True:
       line = inf.readline().strip().split()
       if len(line) == 0: break
       steps += 1
    inf.close()
    return steps
#_______________________________________________________________________________

def inquire_file(infile):
    '''
    inquire file existence
    '''
    if not os.path.exists(infile):
       sys.exit("\n ERROR: File "+infile+" does NOT exist!\n")
    return
#_______________________________________________________________________________

def check_number_of_plots(nop,nmax):
    '''
    print warning if the number_of_plots is larger than 4:
    un update of the colors list could be necessary.
    no_leg = input_options['no_legend']
    '''
    if nop > nmax:
       print("\n WARNING: Number of plots = "+str(nop)+" is larger than "+str(nmax)+" !")
       print("          Updating the list of colors may be necessary.\n")
    return

#END_DEF++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def main(input_options):
    '''
    input:
    :input_options: dictionary that holds the input options parsed 
                    from the command line arguments
    '''
    directory = input_options['directory']
    atype = input_options['atype']
    phonon = input_options['phonon']
    emin = input_options['emin'] 
    emax = input_options['emax'] 
    ezero = input_options['ezero']
    eunit = input_options['eunit']
    funit = input_options['funit']
    sx = input_options['sx']
    sy = input_options['sy']
    kmin = input_options['kmin']
    kmax = input_options['kmax']
    kp = input_options['kp']
    title = input_options['title']
    no_title = input_options['no_title']
    leg_pos = input_options['leg_pos']
    no_leg = input_options['no_legend']
    reverse_colors = input_options['reverse_colors']
    reverse_plots = input_options['reverse_plots']
    maxticksy = input_options['maxticksy']
    legend = input_options['legend']
    show = input_options['show']
    
    #---------------------------------------------------------------------------
    # Initialize cases
    
    if phonon:
        number_of_plots = len(directory)
    else:
        number_of_plots = max(len(directory),len(atype))
        if len(directory) == 1 and len(directory) < len(atype):
            for i in range(1,len(atype)): directory.append(directory[0])
        if len(directory) > len(atype):
            len_atype = len(atype)
            for i in range(len_atype,len(directory)): atype.append("KS")

    check_number_of_plots(number_of_plots,5)

    gw = ( "GW" in atype )
    wannier = ( "WA" in atype )
    
    for i in range(len(directory)): 
        if directory[i] == "./": directory[i] = "."

    local_single = ( number_of_plots==1 and directory[0]=="." )
    local_only = ( all(x == directory[0] for x in directory) and directory[0] == "." )
    
    #---------------------------------------------------------------------------
    # Initialize labels and filenames
    
    leg_spin = ""

    bar = "/"
    if local_only: bar = ""
    
    bandfileroot = "BAND"
    bandlinesfile = "BANDLINES.OUT"
    if atype[0]== "WA": bandlinesfile = "BANDLINES_WANNIER.OUT"
    if phonon:
        bandfileroot = "PHDISP"
        bandlinesfile = "PHDLINES.OUT"

    infile = []
    leg_label = []
    
    for i in range(number_of_plots):
        leg_label.append(directory[i])
        filename = directory[i]+"/"+bandfileroot
        if directory[i] == ".": leg_label[i] = ""
        if not phonon:
            if atype[i] == "GW":
                leg_label[i] = leg_label[i]+bar+'$G_0W_0$'
                filename = filename+"-QP"
            if atype[i] == "WA":
                leg_label[i] = leg_label[i]+bar+'WA'
                filename = filename+"_WANNIER"
            if atype[i] == "KS" and gw: leg_label[i] = leg_label[i] + bar + 'KS'
            if atype[i] == "KS" and wannier: leg_label[i] = leg_label[i] + bar + 'KS'
        infile.append(filename+".OUT")
        if i<len(legend): leg_label[i] = legend[i]
        
    #-------------------------------------------------------------------------------
    # Plot defaults

    xpos_title = 1
    ypos_title = 1.0 + 0.05 / sy
    size_title = "40"
    
    elab = 'Energy $-$ $E_{{F}}$ ['+eunit+']'
    if ezero== "vbM": elab = 'Energy $-$ $E_{_{VBM}}$ [' + eunit + ']'
    if phonon:
        elab = 'Frequency [cm$^{-1}$]'
        if funit!= 'icm': elab = 'Frequency [' + funit + ']'

    xplot_size = 16 * sx
    yplot_size = 9 * sy
    
    line_thickness = "3.0"
    sline_thickness = "1.4"
    axes_thickness = "4.0"
    
    if phonon:
        line_thickness = "4.0"
        axes_thickness = "5.0"
    
    leg_size = 30
    if number_of_plots>=3: leg_size = 27
    
    dpi = 300
    
    figcolor = 'white'
    
    line_color = ["mediumblue", "firebrick", 
                  "green", "darkgoldenrod"]
    length_line_color = len(line_color)
    if number_of_plots > length_line_color:
        for i in range(length_line_color,number_of_plots): 
            line_color.append("darkslategrey")

    if phonon:
        line_color[0] = "firebrick"
        line_color[1] = "mediumblue"

    spin_color = line_color
    
    #-------------------------------------------------------------------------------
    # Check if the calculations are spin-polarized
 
    spin = []  ;  xml_label = [] ; global_spin = False
    for i in range(number_of_plots):
        spin.append(inquire_spin(directory[i] + "/" + "input.xml"))
        xml_label.append(extract_title_text(directory[i] + "/" + "input.xml"))
        global_spin = ( global_spin or spin[i] )
    
    if global_spin: leg_size = 27

    #-------------------------------------------------------------------------------
    # Set title
       
    if local_single and title is None: title = xml_label[0]
        
    #-------------------------------------------------------------------------------
    # Read band-structure data 

    band = []  ;  kvec = []  ;  band_spin = []
  
    for i in range(number_of_plots):
        if phonon:  read_phonon_dispersion(kvec, band, infile[i], funit)
        else: read_electronic_band(kvec,band,band_spin,spin[i],infile[i],eunit)
    
    #-------------------------------------------------------------------------------
    # Set energy zero

    if not phonon:
        for i in range(number_of_plots):
            band[i], band_spin[i] = set_energy_zero(band[i],band_spin[i],spin[i],ezero)

    #-------------------------------------------------------------------------------
    # Read x-ticks position and labels (assuming all plots with the same k-path)

    bandlines = read_xticks_position(directory[0] + "/" + bandlinesfile)
    label = extract_xticks_labels(directory[0] + "/" + "input.xml", phonon)

    #-------------------------------------------------------------------------------
    # Settings for the plot 
    
    fig = plt.figure(figsize=(xplot_size,yplot_size),dpi=dpi)
    fig.patch.set_edgecolor(figcolor)
    fig.patch.set_facecolor(figcolor)

    plt.rcParams['axes.linewidth'] = axes_thickness # set the value globally
    plt.rcParams['grid.linewidth'] = 1.5
    plt.rcParams['xtick.labelsize'] = 40
    plt.rcParams['ytick.labelsize'] = 30
    plt.rcParams['axes.edgecolor'] = 'black'
    plt.rcParams['axes.labelsize'] = 40      # fontsize of the x any y labels
    plt.rcParams['axes.labelcolor'] = 'black'
    plt.rcParams['axes.axisbelow'] = 'True'  # axis gridlines and ticks are below
                                              # the axes elements (lines, text, etc)
    plt.rcParams['legend.fontsize'] = leg_size
    plt.rcParams['xtick.major.pad'] = 10
    plt.rcParams['ytick.major.pad'] = 10

    plt.rcParams.update({'mathtext.default':'regular'})

    #-------------------------------------------------------------------------------
    # Band structure plot 
    
    ax1 = fig.add_subplot(111)
    ax1.xaxis.grid(True,which='major', color='k', linestyle='-', linewidth=line_thickness)
    ax1.axhline(y=0, linestyle="dashed", linewidth=line_thickness, color="black")

    ax1.xaxis.set_label_position('bottom')
    ax1.set_xticks(bandlines)
    ax1.set_xticklabels(label)
    ax1.set_ylabel(elab, labelpad=18 * sx)

    for line in ax1.get_xticklines() + ax1.get_yticklines():
        line.set_markersize(10)
        line.set_markeredgewidth(line_thickness)
        
    if maxticksy is not None:
        ax1.yaxis.set_major_locator(ticker.MaxNLocator(maxticksy))    

    plot_range = list(range(number_of_plots))
    if reverse_plots: plot_range = reversed(list(range(number_of_plots)))
    
    for i in plot_range:
        npspin = not phonon and spin[i]
        lc = line_color[i]
        sc = spin_color[i]
        if reverse_colors:
            lc = line_color[number_of_plots-1-i]
            sc = spin_color[number_of_plots-1-i]
        leg_spin = ["$\\uparrow$","$\\downarrow$"]
        for j in range(len(band[i])):
            llab = set_legend_label(leg_label[i],npspin,leg_spin,local_single,0,j)
            ax1.plot(kvec[i][j], band[i][j], color=lc, lw=line_thickness, label=llab)
            if npspin:
                llab = set_legend_label(leg_label[i],npspin,leg_spin,local_single,1,j)
                ax1.plot(kvec[i][j], band_spin[i][j],
                         color=sc, lw=sline_thickness, label=llab)
 
    if title is not None:  
        if not no_title:
            ax1.text(xpos_title, ypos_title, title, size=size_title, 
                     transform=ax1.transAxes, ha='right', va='center', rotation=0)
    
    if ( ((not local_single or (local_single and global_spin)) and not no_leg) or
         (len(legend)>0) ):       
        leg=ax1.legend(loc=leg_pos,borderaxespad=0.5,
                       framealpha=0.9,fancybox=True)
        leg.get_frame().set_linewidth(float(axes_thickness))
        leg.get_frame().set_edgecolor("grey")
        leg.draw_frame(True)

    if emin is not None: plt.ylim(ymin=emin)
    if emax is not None: plt.ylim(ymax=emax)

    ymin, ymax = ax1.get_ylim()
    
    xmin, xmax = (bandlines[0], bandlines[-1])
    
    for i in range(2):
        if (kp[i] is not None) and (kp[i] not in list(range(1, len(bandlines) + 1))):
            sys.exit("\n ERROR: k-point index = '" + str(kp[i])+
                     "' out of range [1," + str(len(bandlines))+"] !\n")
    
    if kp[0] is not None: xmin = bandlines[kp[0]-1]
    if kp[1] is not None: xmax = bandlines[kp[1]-1]
    
    x0 = xmin
    dx = xmax-xmin
    
    xmin = x0 + kmin*dx
    xmax = x0 + kmax*dx
    
    plt.xlim(xmin, xmax)
    
    plot_ezero(ezero, ymin, ymax, xmin, xmax, kmin, kmax, sx, phonon, ax1)
        
    fig.savefig('PLOT.png', format='png', dpi=300, bbox_inches='tight')

    if show:
       plt.show()

    sys.exit()    

    #-------------------------------------------------------------------------------
    
if __name__ == "__main__":
    input_options = option_parser()
    main(input_options)
