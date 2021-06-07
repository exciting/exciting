#!/usr/bin/env python
# -*- coding: utf-8 -*-
#_______________________________________________________________________________
'''
Plot single and multiple electronic and phonon band structures (BS)

Require the following files:

- for electronic BS: input.xml, BAND.OUT, BAND-QP.OUT, BAND_WANNIER.OUT, BANDLINES.OUT

- for phonon BS: input.xml, PHDISP.OUT, PHLINES.OUT

'''

from   xml.etree import ElementTree as ET
import sys
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import pylab as pyl
import os
import numpy as np
import argparse as ap

import matplotlib.style
if matplotlib.__version__.split(".")[0]=="2": matplotlib.style.use('classic')

#START_DEF++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def option_parser():
    """
    Parse command line inputs 

    Parse: 
        directory
        files
        legend_label
        column_x 
        column_y
        xboundary
        yboundary
        label_x
        label_y
        x_scale
        y_scale
        legend_position
        title
        max_ticks_x
        no_legend
        grid
        scale_box
        reverse_colors
        reverse_plots

    :return input_options: Dictionary of parsed command line arguments 
    """
    p = ap.ArgumentParser(description=\
                'Plot data from single and multiple files/directories from selected columns.')
    
    help_directory = 'List of the directories in which the data to be plotted have to be found. If only one or no directory is specified, the data for the plots are taken from the same directory. Default value is the current directory.'
    
    help_files = 'List of the names of the files in which the data to be plotted have to be found. At least a file must be specified.'
   
    help_legend_label = "Specifies the labels to appear in the legend for each plot."
    
    help_column_x = "List of integer values, each of them specifying the number of the text column in the corresponding file from which the x values must be read. Default value is '1' for each plot. If the value '0' is given, a running natural index is considered as x value." 

    help_column_y = "List of integer values, each of them specifying number of the text column in the file from which the y values must be read. Default value is '1' for each plot. Default value is '2' for each file." 
    
    help_xboundary = "One or two floats corresponding to the minimum and maximum x value in the plot."
    
    help_yboundary = "One or two floats corresponding to the minimum and maximum y value in the plot."

    help_label_x = "A string corresponding to the x label." 
    
    help_label_y = "A string corresponding to the y label." 
    
    help_x_scale = "A float corresponding to a scaling factor for the x values"
    
    help_y_scale = "A float corresponding to a scaling factor for the y values"
   
    help_legend_position = "The location of the legend. The strings 'upper left', 'upper right', 'lower left', 'lower right' place the legend at the corresponding corner of the axes/figure. The strings 'upper center', 'lower center', 'center left', 'center right' place the legend at the center of the corresponding edge of the axes/figure. The string 'center' places the legend at the center of the axes/figure. The string 'best' places the legend at the location, among the nine locations defined so far, with the minimum overlap with other drawn artists. This option can be quite slow for plots with large amounts of data; your plotting speed may benefit from providing a specific location. For back-compatibility, 'center right' (but no other location) can also be spelled 'right', and each string locations can also be given as the corresponding numeric value."             

    help_title = "Used as --title 'String as a title' assign a title to the plot."
    
    help_max_ticks_x = "Specifies the maximum number of ticks along the x-axis in the plot."

    help_no_legend = 'If present, it disables the plotting of the legend.'
    
    help_grid = 'If present, a grid is plotted in correspondence to the position of the major ticks.'

    help_scale_box = "One or two floats corresponding to the scaling factor in the horizontal and vertical size of the plot appearence, respectively."
    
    help_reverse_colors = "If present, the order of the sequence of colors of the plots is reversed."
    
    help_reverse_plots = "If present, the order of appearance of the plots is reversed."

    #---------------------------------------------------------------------------

    p.add_argument('-d','--directory',
                   nargs = '*', default = ["."],
                   type = str, help = help_directory)
    
    p.add_argument('-f','--files',
                   nargs = '*', default = [],
                   type = str, help = help_files)
    
    p.add_argument('-ll','--legend_label',
                   nargs = '*', default = [],
                   type = str, help = help_legend_label)
   
    p.add_argument('-cx','--column_x',
                   nargs = '*', default = [1],
                   type = int, help = help_column_x)
    
    p.add_argument('-cy','--column_y',
                   nargs = '*', default = [2],
                   type = int, help = help_column_y)
    
    p.add_argument('-x','--xboundary',
                   nargs = '*', default = [None, None],
                   type = float, help = help_yboundary)
    
    p.add_argument('-y','--yboundary',
                   nargs = '*', default = [None, None],
                   type = float, help = help_yboundary)
      
    p.add_argument('-lx','--label_x', default = 'x [x]',
                   type = str, help = help_label_x)
    
    p.add_argument('-ly','--label_y', default = 'y [y]',
                   type = str, help = help_label_y)

    p.add_argument('-xs','--x_scale', default = 1.0,
                   type = float, help = help_x_scale)
    
    p.add_argument('-ys','--y_scale', default = 1.0,
                   type = float, help = help_y_scale)
    
    p.add_argument('-mtx','--max_ticks_x',
                   type = int, default = None, help = help_max_ticks_x)

    p.add_argument('-s','--scale_box',
                   nargs = '*', default = [1.0, 1,0],
                   type = float, help = help_scale_box)
                   
    p.add_argument('-t','--title',
                   type = str, default = None, help = help_title)
   
    p.add_argument('-lp','--legend_position',
                   type = str, help = help_legend_position,
                   choices = ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10',
                              'best', 'upper right', 'upper left', 'lower left', 'lower right', 'right', 'center left', 'center right', 'lower center', 'upper center', 'center'],
                   default = 'best')
                   
    p.add_argument('-nl','--no_legend', action='store_true', help = help_no_legend)
    
    p.add_argument('-g','--grid', action='store_true', help = help_grid)
    
    p.add_argument('-rc','--reverse_colors', action='store_true', help = help_reverse_colors)

    p.add_argument('-rp','--reverse_plots', action='store_true', help = help_reverse_plots)
    
    #---------------------------------------------------------------------------

    args = p.parse_args()
    input_options = {}

    input_options['directory'] = args.directory
    
    if ( len(args.files)==0 ):
        sys.exit("\n ERROR: At least a filename must be specified:\n\n" 
                 +"        PLOT-files.py -f [FILES [FILES ...]]\n")
    input_options['files'] = args.files
    
    input_options['legend'] = args.legend_label
    input_options['columnx'] = args.column_x
    input_options['columny'] = args.column_y

    input_options['sx'] = 1.0    
    if ( len(args.scale_box) >= 1 ): input_options['sx'] = args.scale_box[0]
    input_options['sy'] = 1.0
    if ( len(args.scale_box) >= 2 ): input_options['sy'] = args.scale_box[1]

    input_options['xmin'] = None
    if ( len(args.xboundary) >= 1 ): input_options['xmin'] = args.xboundary[0]
    input_options['xmax'] = None
    if ( len(args.xboundary) >= 2 ): input_options['xmax'] = args.xboundary[1]
    
    input_options['ymin'] = None
    if ( len(args.yboundary) >= 1 ): input_options['ymin'] = args.yboundary[0]
    input_options['ymax'] = None
    if ( len(args.yboundary) >= 2 ): input_options['ymax'] = args.yboundary[1]

    input_options['labelx'] = args.label_x
    input_options['labely'] = args.label_y
    
    input_options['xscale'] = args.x_scale
    input_options['yscale'] = args.y_scale
    
    input_options['title'] = args.title
    input_options['maxticksx'] = args.max_ticks_x
    input_options['no_legend'] = args.no_legend
    input_options['grid'] = args.grid
    input_options['leg_pos'] = args.legend_position
    if ( len(args.legend_position)<=2 ):
        input_options['leg_pos'] = int(args.legend_position)
      
    input_options['reverse_plots'] = args.reverse_plots
    input_options['reverse_colors'] = args.reverse_colors
  
    return input_options
#_______________________________________________________________________________
 
def read_data(x,y,cx,cy,nop,sx,sy,infile):
    '''
    read data frome files infile 
    '''
    for i in range(nop):
        inquire_file(infile[i])
        ix = cx[i]-1  ;  iy = cy[i]-1
        listx = []  ;  listy = []
        for line in open(infile[i]).readlines():
            noline = ( line.strip()[0][:1] != "#" and
                       line.strip()[0][:1] != "&" and
                       line.strip()[0][:1] != "@" )
            if (len(line.strip()) != 0 and noline):
                i_line = line.strip().split()
                listy.append(float(i_line[iy])*sy)
                if ( cx[i]!=0 ): 
                    listx.append(float(i_line[ix])*sx)
                else:
                    listx.append(len(listy)*sx)
        x.append(listx)
        y.append(listy)
    return 
#_______________________________________________________________________________

def inquire_file(infile):
    '''
    inquire file existence
    '''
    if ( not os.path.exists(infile) ):
        sys.exit("\n ERROR: File "+infile+" does NOT exist!\n")
    return
#_______________________________________________________________________________

def check_number_of_plots(nop,nmax):
    '''
    inquire file existence
    '''
    if ( nop>nmax ):
        sys.exit("\n ERROR: Number of plots = "+str(nop)
                 +" is larger than "+str(nmax)+" !\n")
    return

#END_DEF++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def main(input_options):
    '''
    input:
    :input_options: dictionary that holds the input options parsed 
                    from the command line arguments
    '''
    directory = input_options['directory']
    files = input_options['files']
    legend = input_options['legend']
    columnx = input_options['columnx']
    columny = input_options['columny']
    
    xxmin = input_options['xmin'] 
    xxmax = input_options['xmax'] 

    yymin = input_options['ymin'] 
    yymax = input_options['ymax']
    
    labelx = input_options['labelx']
    labely = input_options['labely']
    
    xs = input_options['xscale']
    ys = input_options['yscale']
    
    sx = input_options['sx']
    sy = input_options['sy']
    
    maxticksx = input_options['maxticksx']
    title = input_options['title']
    leg_pos = input_options['leg_pos']
    no_leg = input_options['no_legend']
    grid = input_options['grid']

    reverse_colors = input_options['reverse_colors']
    reverse_plots = input_options['reverse_plots']
    
    #---------------------------------------------------------------------------
    # Initialize cases
    
    number_of_plots = max(len(files),len(directory))
    
    only_one_directory = ( len(directory)==1 )
    only_one_file = ( len(files)==1 )
    only_one_plot = ( number_of_plots==1 )
    one_file_for_directory = ( len(files)==len(directory) )
    
    if ( not only_one_file and 
         not only_one_directory and 
         not one_file_for_directory):
        sys.exit("\n ERROR: The number of specified files ("
                 +str(len(files))+") and specified directories ("
                 +str(len(directory))+") must be equal in this case!\n")

    if ( only_one_directory and not one_file_for_directory ):
            for i in range(1,len(files)): directory.append(directory[0])
    if ( only_one_file and not one_file_for_directory ):
            for i in range(1,len(directory)): files.append(files[0])

    same_columnx = ( len(columnx)==1 )
    same_columny = ( len(columny)==1 )

    if ( not same_columnx and (len(columnx)!=number_of_plots) ):
        sys.exit("\n ERROR: The value of the x-column number"
                 +" must be specified in this case for each plot ("
                 +str(number_of_plots)+"!\n")
    
    if ( not same_columny and (len(columny)!=number_of_plots) ):
        sys.exit("\n ERROR: The value of the y-column number"
                 +" must be specified in this case for each plot ("
                 +str(number_of_plots)+"!\n")

    if ( same_columnx and not only_one_plot ):
            for i in range(1,len(files)): columnx.append(columnx[0])
    if ( same_columny and not only_one_plot ):
            for i in range(1,len(files)): columny.append(columny[0])
            
    check_number_of_plots(number_of_plots,8)
    
    for i in range(len(directory)): 
        if ( directory[i]=="./" ): directory[i] = "."
        
    local_only = ( all(x==directory[0] for x in directory) and directory[0]=="." )

    #---------------------------------------------------------------------------
    # Initialize labels, filenames, and title

    bar = "/"
    if ( local_only ): bar = ""
    
    infile = []
    leg_label = []
   
    for i in range(number_of_plots):
        infile.append(directory[i]+"/"+files[i])
        if ( title is None and only_one_plot ): 
            leg_label.append("")
            title = directory[0]+bar+files[0]
            if ( local_only ): title = files[0]
            no_leg = True
        elif ( only_one_file and not only_one_plot):
            leg_label.append(directory[i])
            if ( title is None ): title = files[i]
        elif ( only_one_directory and not only_one_plot):
            leg_label.append(files[i])
            if ( title is None ): 
                title = directory[i]
                if (local_only): title = ""
        else:
            leg_label.append(directory[i]+bar+files[i])
        if ( i<len(legend) ): leg_label[i] = legend[i] 

    #-------------------------------------------------------------------------------
    # Plot defaults

    xplot_size = 14*sx
    yplot_size = 9*sy
    
    line_thickness = "4.0"
    axes_thickness = "4.0"
    
    leg_size = 30
    if ( number_of_plots>=3 ): leg_size = 27
    
    dpi = 300
    
    figcolor = 'white'
    
    line_color = ["mediumblue", "firebrick", "green", "orange",
                  "darkcyan", "purple", "darkgoldenrod", "grey"]
    
    #-------------------------------------------------------------------------------
    # Read data 

    xx = [] ; yy = [] 
  
    read_data(xx,yy,columnx,columny,number_of_plots,xs,ys,infile)
 
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
    # Plot 
    
    ax1 = fig.add_subplot(111)

    ax1.xaxis.set_label_position('bottom')
    ax1.set_xlabel(labelx, labelpad=10)
    ax1.set_ylabel(labely, labelpad=20)

    if ( maxticksx is not None ): 
        ax1.xaxis.set_major_locator(plt.MaxNLocator(maxticksx))

    for line in ax1.get_xticklines() + ax1.get_yticklines():
        line.set_markersize(10)
        line.set_markeredgewidth(line_thickness)

    plot_range = range(number_of_plots)
    if ( reverse_plots ): plot_range = reversed(range(number_of_plots))

    for i in plot_range:
        lc = line_color[i]
        if ( reverse_colors ): 
            lc = line_color[number_of_plots-1-i]
        ax1.plot(xx[i], yy[i], color=lc, lw=line_thickness, label=leg_label[i])
                
    if title is not None:  
        ax1.text(1,1.05,title,size="40",
                 transform=ax1.transAxes,ha='right',va='center',rotation=0)

    if ( not no_leg ):       
        leg=ax1.legend(loc=leg_pos,borderaxespad=0.7,
                       framealpha=0.9,fancybox=True)
        leg.get_frame().set_linewidth(axes_thickness)
        leg.get_frame().set_edgecolor("grey")
        leg.draw_frame(True)

    if xxmin is not None: plt.xlim(xmin=xxmin)
    if xxmax is not None: plt.xlim(xmax=xxmax)
    if yymin is not None: plt.ylim(ymin=yymin)
    if yymax is not None: plt.ylim(ymax=yymax)
    
    if (grid): pyl.grid(True)

    fig.savefig('PLOT.png',format='png',dpi=300, bbox_inches='tight')
    fig.savefig('PLOT.eps',format='eps',bbox_inches='tight')

    sys.exit()    

    #-------------------------------------------------------------------------------
    
if __name__ == "__main__":
    input_options = option_parser()
    main(input_options)
