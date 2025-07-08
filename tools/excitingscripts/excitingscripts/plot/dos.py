"""Visualize desisty of states.

More details can be found **[here](https://www.exciting-code.org/home/the-python-script-plot.dos)**.

Located at `excitingscripts/plot/dos.py`.

Call as:

```bash
python3 -m excitingscripts.plot.dos
```
"""


#_______________________________________________________________________________
import argparse as ap
import os
import sys
from xml.etree import ElementTree as ET

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import matplotlib.style
import matplotlib.ticker as ticker
import numpy as np
import pylab as pyl

if matplotlib.__version__.split(".")[0]=="2": matplotlib.style.use('classic')

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
        scale_box
        dos_boundary
        eunit
        funit
        title
        no_title
        legend_position
        scale_box
        reverse_colors
        reverse_plots
        no_fill
        no_reverse_spin
        max_ticks_x
        max_ticks_y
        legend_label
        grid
        show
         
    :return input_options: Dictionary of parsed command line arguments 
    """
    p = ap.ArgumentParser(description=\
                'Plot single and multiple electronic/phonon density of sataes.')
    
    help_directory = 'List of the directories in which the data to be plotted have to be found.    If only one or no directory is specified, the data for the plots are taken from the same directory. Default value is the current directory.'
    
    help_eboundary = 'One or two floats corresponding to the minimum and maximum energy in the plot of the electronic DOS in the units specified by --eunit, respectively. If the argument --phonon is present, the option values specify the minimum and maximum frequency in cm-1 appearing in the phonon-dispersion plot. If either the argument --eboundary is not present or no floats are given, the energy (frequency) boundaries are chosen correspondingly to the maximum and minimum of the full data to be plotted.'
    
    help_dos_boundary = 'One or two floats corresponding to the minimum and maximum density in the plot of the density of states (DOS) according to the units specified by --eunit, respectively. If the argument --phonon is present, the option values specify the minimum and maximum DOS, according to the units specified by --funit. If either the argument --dos_boundary is not present or no floats are given, the DOS boundaries are chosen correspondingly to the maximum and minimum of the full data to be plotted. Please notice that in the case of spin-polarized system the electronic spin-down DOS is represented by default by negative values.'
    
    help_assign_type = "List of the description keys for each plot of the electronic DOS. Possible choices are 'KS' (standard Kohn-Sham calculation), 'GW' (G0W0 calculation), and 'WA' (data are interpolated by using Wannier functions). If not present, the option 'KS' is assumed for all plots. Not used if the argument --phonon is present."
        
    help_eunit = "Set the units of the energy appearing in the plot of the electronic density of states. Possible choices are 'eV' (electronvolt, default) and 'Ha' (Hartree)."

    help_funit = "Set the units of the frequency appearing in the plot of the phonon density of states. Possible choices are 'icm' (inverse centimeter, cm^-1, default), 'meV' (millielectronvolt), and 'THz' (terahertz)."
    
    help_legend_position = "The location of the legend. The strings 'upper left', 'upper right', 'lower left', 'lower right' place the legend at the corresponding corner of the axes/figure. The strings 'upper center', 'lower center', 'center left', 'center right' place the legend at the center of the corresponding edge of the axes/figure. The string 'center' places the legend at the center of the axes/figure. The string 'best' places the legend at the location, among the nine locations defined so far, with the minimum overlap with other drawn artists. This option can be quite slow for plots with large amounts of data; your plotting speed may benefit from providing a specific location. For back-compatibility, 'center right' (but no other location) can also be spelled 'right', and each string locations can also be given as the corresponding numeric value."
                  
    help_phonon = 'If present, it tags the plotting of the phonon DOS. If absent, the electronic DOS is plotted.'
    
    help_no_legend = 'If present, it disables the plotting of the legend.'
    
    help_no_fill = 'If present, plots are not filled.'
     
    help_title = "Used as --title 'String as a title' assign a title to the plot."
   
    help_no_title = 'If present, it disables the writing of the title.'
    
    help_scale_box = "One or two floats corresponding to the scaling factor in the horizontal and vertical size of the plot appearence, respectively."
    
    help_reverse_colors = "If present, the order of the sequence of colors of the plots is reversed."
    
    help_reverse_plots = "If present, the order of appearance of the plots is reversed."
    
    help_no_reverse_spin = "If present and a spin-polarized density of state is plotted, the spin-down DOS is represented by positive values instead of negative ones."
    
    help_max_ticks_x = "Specifies the maximum number of ticks along the x-axis in the plot."

    help_max_ticks_y = "Specifies the maximum number of ticks along the y-axis in the plot."
    
    help_legend_label = "Specifies the labels to appear in the legend for each plot."

    help_grid = 'If present, a grid is plotted in correspondence to the position of the major ticks.'

    help_show = "Opens plot in new window"
    
    #---------------------------------------------------------------------------
    
    p.add_argument('-d','--directory',
                   nargs = '*', default = ["."],
                   type = str, help = help_directory)
    
    p.add_argument('-e','--eboundary',
                   nargs = '*', default = [None, None],
                   type = float, help = help_eboundary)
    
    p.add_argument('-db','--dos_boundary',
                   nargs = '*', default = [None, None],
                   type = float, help = help_dos_boundary)
    
    p.add_argument('-s','--scale_box',
                   nargs = '*', default = [1.0, 1,0],
                   type = float, help = help_scale_box)

    p.add_argument('-a','--assign_type',
                   nargs = '*', default = ['KS'],
                   choices = ['KS', 'GW', 'WA'],
                   type = str, help = help_assign_type)
    
    p.add_argument('-p','--phonon', action='store_true', help = help_phonon)
    
    p.add_argument('-nf','--no_fill', action='store_true', help = help_no_fill)
    
    p.add_argument('-rc','--reverse_colors', action='store_true', help = help_reverse_colors)

    p.add_argument('-rp','--reverse_plots', action='store_true', help = help_reverse_plots)

    p.add_argument('-nl','--no_legend', action='store_true', help = help_no_legend)
    
    p.add_argument('-nrs','--no_reverse_spin', action='store_true', help = help_no_reverse_spin)
        
    p.add_argument('-eu','--eunit',
                   type = str, help = help_eunit,
                   choices = ['eV', 'Ha'], default = 'eV')
    
    p.add_argument('-fu','--funit',
                   type = str, help = help_funit,
                   choices = ['icm', 'meV', 'THz'], default = 'icm')
    
    p.add_argument('-t','--title',
                   type = str, default = None, help = help_title)
    
    p.add_argument('-nt','--no_title', action='store_true', help = help_no_title)
    
    p.add_argument('-mtx','--max_ticks_x',
                   type = int, default = None, help = help_max_ticks_x)
    
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
                   
    p.add_argument('-g','--grid', action='store_true', help = help_grid)

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
    
    input_options['dmin'] = None
    if len(args.dos_boundary) >= 1: input_options['dmin'] = args.dos_boundary[0]
    input_options['dmax'] = None
    if len(args.dos_boundary) >= 2: input_options['dmax'] = args.dos_boundary[1]

    input_options['sx'] = 1.0    
    if len(args.scale_box) >= 1: input_options['sx'] = args.scale_box[0]
    input_options['sy'] = 1.0
    if len(args.scale_box) >= 2: input_options['sy'] = args.scale_box[1]

    input_options['phonon'] = args.phonon
    
    input_options['title'] = args.title
    input_options['no_title'] = args.no_title
    
    input_options['eunit'] = args.eunit
    input_options['funit'] = args.funit
    
    input_options['maxticksx'] = args.max_ticks_x
    input_options['maxticksy'] = args.max_ticks_y
    
    input_options['legend'] = args.legend_label

    input_options['no_legend'] = args.no_legend
    input_options['leg_pos'] = args.legend_position
    if len(args.legend_position)<=2:
        input_options['leg_pos'] = int(args.legend_position)
        
    input_options['no_reverse_spin'] = args.no_reverse_spin
    input_options['reverse_plots'] = args.reverse_plots
    input_options['reverse_colors'] = args.reverse_colors
    input_options['no_fill'] = args.no_fill
    input_options['grid'] = args.grid

    input_options['show'] = args.show
    
    return input_options
#_______________________________________________________________________________

def extract_title_text(xmlfile):
    '''
    extract text in title from xmlfile
    '''
    inquire_file(xmlfile)
    text = ET.parse(xmlfile).getroot().find("title").text
    return text
#_______________________________________________________________________________

def inquire_spin(xmlfile):
    '''
    check in input.xml if a spin-polarized calculation is performed
    '''
    inquire_file(xmlfile)
    lspin = False
    groundstate = ET.parse(xmlfile).getroot().find("groundstate")
    for spin in groundstate.findall("spin"): 
        so = str(spin.get("spinorb")).lower()
        if not so == 'true': lspin = True
    return lspin
#_______________________________________________________________________________

def set_legend_label(leg_label,spin,leg_spin,local_single,k):
    '''
    set legend label for the plot
    '''
    legend_label = leg_label
    if local_single: leg_spin[k] = "spin" + leg_spin[k]
    if spin: legend_label = leg_label + leg_spin[k]
    return legend_label
#_______________________________________________________________________________

def inquire_element(element,xmlpath):
    '''
    check if element is in xmlfile
    '''
    if xmlpath.find(element) is None:
        sys.exit("\n WARNING: Element '"+element+"' not found in "+str(xmlpath)+"\n")
    return xmlpath.find(element)
#_______________________________________________________________________________            

def find_steps(infile):
    inf = open(infile,"r")
    steps = 0
    while True:
       line=inf.readline().strip().split()
       if len(line)==0: break
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
    if nop>nmax:
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
    dmin = input_options['dmin'] 
    dmax = input_options['dmax'] 
    eunit = input_options['eunit']
    funit = input_options['funit']
    sx = input_options['sx']
    sy = input_options['sy']
    title = input_options['title']
    no_title = input_options['no_title']
    leg_pos = input_options['leg_pos']
    no_leg = input_options['no_legend']
    no_fill = input_options['no_fill']
    reverse_colors = input_options['reverse_colors']
    reverse_plots = input_options['reverse_plots']
    no_reverse_spin = input_options['no_reverse_spin']
    maxticksx = input_options['maxticksx']
    maxticksy = input_options['maxticksy']
    legend = input_options['legend']
    grid = input_options['grid']
    show = input_options['show']
    
    #---------------------------------------------------------------------------
    # Initialize cases
    
    if phonon:
        number_of_plots = len(directory)
    else:
        number_of_plots = max(len(directory),len(atype))
        if len(directory)==1 and len(directory)<len(atype):
            for i in range(1,len(atype)): directory.append(directory[0])
        if len(directory)>len(atype):
            len_atype = len(atype)
            for i in range(len_atype,len(directory)): atype.append("KS")

    check_number_of_plots(number_of_plots,5)

    gw = ( "GW" in atype )
    wannier = ( "WA" in atype )
    
    for i in range(len(directory)): 
        if directory[i]== "./": directory[i] = "."

    local_single = ( number_of_plots==1 and directory[0]=="." )
    local_only = ( all(x==directory[0] for x in directory) and directory[0]=="." )
            
    #---------------------------------------------------------------------------
    # Initialize labels and filenames
    
    leg_spin = ""

    bar = "/"
    if local_only: bar = ""
    
    bandfileroot = "TDOS"
    if phonon:
        bandfileroot = "PHDOS"
 
    infile = []
    leg_label = []
    
    for i in range(number_of_plots):
        leg_label.append(directory[i])
        filename = directory[i]+"/"+bandfileroot
        if directory[i]== ".": leg_label[i] = ""
        if not phonon:
            if atype[i]== "GW":
                leg_label[i] = leg_label[i]+bar+'$G_0W_0$'
                filename = filename+"-QP"
            if atype[i]== "WA":
                leg_label[i] = leg_label[i]+bar+'WA'
                filename = filename+"_WANNIER"
            if atype[i]== "KS" and gw: leg_label[i] = leg_label[i] + bar + 'KS'
            if atype[i]== "KS" and wannier: leg_label[i] = leg_label[i] + bar + 'KS'
        infile.append(filename+".OUT")
        if i<len(legend): leg_label[i] = legend[i]

    #-------------------------------------------------------------------------------
    # Plot defaults
    
    xpos_title = 1
    ypos_title = 1.0+0.05/sy
    size_title = "40"  
    
    elab = 'Energy$-E_F$ ['+eunit+']'
    if phonon:
        elab = 'Frequency [cm$^{-1}$]'
        if funit!= 'icm': elab = 'Frequency [' + funit + ']'
    
    dlab = 'DOS [states/'+eunit+'/unit cell]'
    if phonon:
        dlab = 'Phonon DOS [states/cm$^{-1}$]'
        if funit!= 'icm':
            dlab = 'Phonon DOS [states/'+funit+']'
            if funit== 'THz': dlab = r'Phonon DOS [states/$\,$' + funit + ']'

    xplot_size = 16*sx
    yplot_size = 9*sy
    
    line_thickness = "3.0"
    sline_thickness = "2.0"
    axes_thickness = "4.0"
    leg_size = 30
    if number_of_plots>=3: leg_size = 27
    
    dpi = 300
    
    figcolor = 'white'
    
    line_color = ["mediumblue", "firebrick", "green", "darkgoldenrod"]
    fill_color = ["cornflowerblue", "lightsalmon", "lightgreen", "moccasin"]
    length_line_color = len(line_color)
    if number_of_plots>length_line_color:
        for i in range(length_line_color,number_of_plots): 
            line_color.append("darkslategrey")
            fill_color.append("white")

    if phonon:
        line_color[0] = "firebrick"
        line_color[1] = "mediumblue"
        fill_color[0] = "lightsalmon"
        fill_color[1] = "cornflowerblue"
    
    spin_color = line_color
    sfill_color = fill_color
    
    #-------------------------------------------------------------------------------
    # Check if the calculations are spin-polarized
 
    spin = []  ;  xml_label = [] ; global_spin = False
    for i in range(number_of_plots):
        spin.append(inquire_spin(directory[i]+"/"+"input.xml"))
        xml_label.append(extract_title_text(directory[i]+"/"+"input.xml"))
        global_spin = ( global_spin or spin[i] )
    
    if global_spin: leg_size = 27

    #-------------------------------------------------------------------------------
    # Set title
       
    if local_single and title is None: title = xml_label[0]
        
    #-------------------------------------------------------------------------------
    # Read DOS data

    ene = []  ;  dos = []  ;  ene_spin = []  ;  dos_spin = []
    
    au2icm = 2.194746313705e5
    icm2mev = 1.0/8.06573
    icm2thz = 0.0299793
    ha2ev = 27.211396132 

    if phonon:
        conversion_factor = au2icm
        if funit== 'meV': conversion_factor = au2icm * icm2mev
        if funit== 'THz': conversion_factor = au2icm * icm2thz
    else:
        conversion_factor = ha2ev
        if eunit== "Ha": conversion_factor = 1.0

    for i in range(number_of_plots):
        inquire_file(infile[i])
        elist = np.genfromtxt(infile[i])[:,0]*conversion_factor
        dlist = np.genfromtxt(infile[i])[:,1]/conversion_factor
        if ( spin[i] ):
            ndos = int(len(elist)/2)
            edown = [elist[j] for j in range(ndos,2*ndos)]
            ddown = [dlist[j] for j in range(ndos,2*ndos)]
            elist = [elist[j] for j in range(ndos)] 
            dlist = [dlist[j] for j in range(ndos)]
        len_e = len(elist)
        #elist = [elist[j] for j in range(len_e) if dlist[j]>=1e-6]
        #dlist = [dlist[j] for j in range(len_e) if dlist[j]>=1e-6]
        ene.append(elist)   
        dos.append(dlist)
        if spin[i]:
            sfactor = 1.0
            if no_reverse_spin: sfactor = -1.0
            len_e = len(edown)
            edown = [edown[j] for j in range(len_e)]# if ddown[j]<=sfactor*1e-6]
            ddown = [sfactor*ddown[j] for j in range(len_e)]# if ddown[j]<=sfactor*1e-6]
            ene_spin.append(edown)
            dos_spin.append(ddown)
 
    #-------------------------------------------------------------------------------
    # Settings for the plot 
    
    fig = plt.figure(figsize=(xplot_size,yplot_size),dpi=dpi)
    fig.patch.set_edgecolor(figcolor)
    fig.patch.set_facecolor(figcolor)

    plt.rcParams['axes.linewidth']  = axes_thickness # set the value globally
    plt.rcParams['grid.linewidth']  = 1.5
    plt.rcParams['xtick.labelsize'] = 30
    plt.rcParams['ytick.labelsize'] = 30
    plt.rcParams['axes.edgecolor']  = 'black'
    plt.rcParams['axes.labelsize']  = 40      # fontsize of the x any y labels
    plt.rcParams['axes.labelcolor'] = 'black'
    plt.rcParams['axes.axisbelow']  = 'True'  # axis gridlines and ticks are below
                                              # the axes elements (lines, text, etc)
    plt.rcParams['legend.fontsize'] = leg_size
    plt.rcParams['xtick.major.pad'] = 10
    plt.rcParams['ytick.major.pad'] = 10

    plt.rcParams.update({'mathtext.default':'regular'})

    #-------------------------------------------------------------------------------
    # DOS plot
    
    ax1 = fig.add_subplot(111)

    ax1.xaxis.set_label_position('bottom')
    ax1.set_xlabel(elab)
    ax1.set_ylabel(dlab, labelpad=20)

    if maxticksx is not None:
        ax1.xaxis.set_major_locator(ticker.MaxNLocator(maxticksx))
        
    if maxticksy is not None:
        ax1.yaxis.set_major_locator(ticker.MaxNLocator(maxticksy))
       
    transparency = 1.0
    my_legend = []
       
    plot_range = list(range(number_of_plots))
    if reverse_plots: plot_range = reversed(list(range(number_of_plots)))
    
    for i in plot_range:
        npspin = not phonon and spin[i]
        lc = line_color[i]
        sc = spin_color[i]
        fc = fill_color[i]
        if reverse_colors:
            lc = line_color[number_of_plots-1-i]
            sc = spin_color[number_of_plots-1-i]
            fc = fill_color[number_of_plots-1-i]
        leg_spin = ["",""]
        if no_reverse_spin and not phonon: leg_spin = ["$\\uparrow$", "$\\downarrow$"]
        llab = set_legend_label(leg_label[i],npspin,leg_spin,local_single,0)
        if not no_fill: ax1.fill_between(ene[i], dos[i], color=fc, alpha=transparency)
        ax1.plot(ene[i], dos[i], color=lc, lw=line_thickness, label=llab)
        thislegend = mpatches.Patch(facecolor=fc, edgecolor=lc, lw=line_thickness, label=llab)
        my_legend.append(thislegend)
        if spin[i]:
            llab = set_legend_label(leg_label[i],npspin,leg_spin,local_single,1)
            if not no_fill and not no_reverse_spin:
                ax1.fill_between(ene_spin[i], dos_spin[i], color=fc, alpha=transparency/2.0)
            ax1.plot(ene_spin[i], dos_spin[i], color=sc, lw=sline_thickness, label=llab) 
            if no_reverse_spin:
                thislegend = mpatches.Patch(facecolor="white", edgecolor=sc, 
                                            lw=sline_thickness, label=llab)
                my_legend.append(thislegend)
            
    for line in ax1.get_xticklines() + ax1.get_yticklines():
        line.set_markersize(10)
        line.set_markeredgewidth(line_thickness)
                     
    plt.axvline(x=0, lw=3.0, color="black", linestyle='dashed') 
    plt.axhline(y=0, lw=3.0, color="black", linestyle='-') 
                
    if title is not None:  
        if not no_title:
            ax1.text(xpos_title, ypos_title, title, size=size_title, 
                     transform=ax1.transAxes, ha='right', va='center', rotation=0)

    xpos_spin = 0.03/sx
    ha_spin = 'left'
    if local_single or no_leg:
        xpos_spin = 1-0.03/sx  
        ha_spin = 'right'
    
    if ( (not local_single and not no_leg) or 
         (no_reverse_spin and not phonon) or 
         (len(legend)>0) ):       
        leg=ax1.legend(handles=my_legend,loc=leg_pos,borderaxespad=0.7,
                      framealpha=0.9,fancybox=True)
        leg.get_frame().set_linewidth(float(axes_thickness))
        leg.get_frame().set_edgecolor("grey")
        leg.draw_frame(True)

    ymin, ymax = (np.min(dos), ax1.get_ylim()[1])
    xmin, xmax = (np.min(ene), np.max(ene))

    if emin is not None: xmin=emin
    if emax is not None: xmax=emax

    if dmin is not None: ymin=dmin
    if dmax is not None: ymax=dmax
    
    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)
    
    props = dict(boxstyle="Round, pad=0.2", facecolor="beige", 
                 edgecolor='palegoldenrod', lw=5,alpha=0.9)   
  
    if global_spin and not no_reverse_spin:
        if ymin<0.0:
            ax1.text(xpos_spin, 0.08/sy, "spin$\\downarrow$", size="30", #color="darkblue",
                     transform=ax1.transAxes, ha=ha_spin, va='center', 
                     rotation=0, bbox = props)
        if ymax>0.0:
            ax1.text(xpos_spin, 1-0.08/sy, "spin$\\uparrow$", size="30", #color="darkblue",
                     transform=ax1.transAxes, ha=ha_spin, va='center', 
                     rotation=0, bbox = props) 
            
    if grid: pyl.grid(True)

    fig.savefig('PLOT.png',format='png',dpi=300, bbox_inches='tight')

    if show:
       plt.show()

    sys.exit()    

    #-------------------------------------------------------------------------------
    
if __name__ == "__main__":
    input_options = option_parser()
    main(input_options)
