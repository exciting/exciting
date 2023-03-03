#!/usr/bin/python2
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import os
import sys

import matplotlib.style
if matplotlib.__version__.split(".")[0]=="2": matplotlib.style.use('classic')

def parse_lossfunction(fname):
    """
    Parses files containing loss function
    e.g. LOSS_FXCRPA_OC11_QMT001.OUT

    :param str fname: name of the file
    """
    xdata = []
    ydata = []
    with open(fname, 'r') as file:
        for lines in file:
            if 'Frequency' in lines:
                break
        for lines in file:
            data = lines.split()
            xdata.append(float(data[0]))
            ydata.append(float(data[1]))
    return xdata, ydata

#-------------------------------------------------------------------------------

def shell_value(variable, vlist, default):
    v = default
    e = False
    for i in range(len(vlist)):
        if vlist[i] == variable:
            v = os.environ[variable]
            e = True
            break
    return v, e

#-------------------------------------------------------------------------------

current = os.environ['PWD']
ev_list = os.environ.keys()

rundir = shell_value('EXCITINGRUNDIR',ev_list,current)[0]
rlabel = shell_value('RLABEL',ev_list,"rundir-")[0]
showpyplot = shell_value('SHOWPYPLOT',ev_list,"")[1]
dpipng = int(shell_value('DPIPNG',ev_list,300)[0])

#-------------------------------------------------------------------------------
#Check arguments

nfiles=0
handler=""

print
fnames=[]
for i in sys.argv[1:]:
    if i[:2] != '--':
        nfiles=nfiles+1
        fnames.append(i)
        if not os.path.isfile(i):
           print " Error: file \"%s\" doesn't exist!\n"%(i)
           sys.exit()
        print ' file: '+i
    else:
        handler=i[2:]

if nfiles > 1 and handler == "":
    print "\n Error: The legend handler must be given if more than one file is plotted!"
    print "\n **Usage type 1**:    PLOT-loss-function.py --legend-handler lossfile-1.OUT lossfile-2.OUT"
    print   " **Usage type 2**:    PLOT-loss-function.py lossfile.OUT\n"
    sys.exit()

if nfiles==1: handler="kernel"

if not handler in ["qmt","kernel"]:
    print "\n Error: _\""+handler+"\"_ The legend handler has to be either 'kernel' or 'qmt'!"
    print "\n **Usage type 1**:    PLOT-loss-function.py --legend-handler lossfile-1.OUT lossfile-2.OUT"
    print   " **Usage type 2**:    PLOT-loss-function.py lossfile.OUT\n"
    sys.exit()

if nfiles<1:
    print "\n Error: At least one input file should be given!"
    print "\n** Usage type 1**:    PLOT-loss-function.py --legend-handler lossfile-1.OUT lossfile-2.OUT"
    print   "** Usage type 2**:    PLOT-loss-function.py lossfile.OUT\n"
    sys.exit()

print

#-------------------------------------------------------------------------------
#Parse LOSS function data files

xdata=[] ; ydata=[] ; labels=[]; legends=[]
for i,fname in enumerate(fnames):
    labels.append({})
    sfname = fname.split("/")[-1].split("_")

    #legends for BSE calculations

    if "BSE" in sfname[1].split("-"):
        if handler=='qmt':
            legend=fname.split(".")[-3].split("_")[-1]
            if legend[:2]=="OC": legend='QMT001('+legend[2:]+')'
        if handler=='kernel':
            print "\n Error: BSE Calculations do not have kernels!\n"
            sys.exit()

    #legends for TDDFT calculations
    else:
        if handler=="kernel":
            if "FXCRPA" in sfname: legend="RPA "
            if "FXCALDA" in sfname: legend="ALDA "
            if not "NLF" in sfname: legend=legend+"(LFE) "
            if "NLF" in sfname: legend=legend+"(no-LFE) "
            if "toplot" in sfname: legend=sfname[0]+" "+legend
        if handler=="qmt":
            if sfname[-2][0:2]=="OC": legend="Optical(%s)"%(sfname[-2][2:])
            else:
                legend=fname.split('.')[-2].split('_')[-1]

    legends.append(legend)
    labels[i]["xlabel"] = "Energy"
    labels[i]["ylabel"] = "Loss Function"
    x, y = parse_lossfunction( fname )
    xdata.append(x)
    ydata.append(y)


#-------------------------------------------------------------------------------
#Plot LOSS function/s

colors=['k','r','g','b','y','c','m']
fig=plt.figure(1,figsize=(8,5.5))

params = {'font.size'             : 20,
          'xtick.major.width'     :  2,
          'ytick.major.width'     :  2,
          'xtick.major.size'      :  5,
          'ytick.major.size'      :  5,
          'patch.linewidth'       :  1.5,
          'axes.linewidth'        :  2.,
          'axes.formatter.limits' : (-4, 6),
          'lines.linewidth'       : 1.8,
          'lines.markeredgewidth' : 2.0,
          'lines.markersize'      : 18,
          'legend.fontsize'       : 16,
          'legend.borderaxespad'  :  1,
          'legend.borderpad'      :  0.5,
          'savefig.dpi'           : 80}

plt.rcParams.update(params)

ax = fig.add_axes([0.16,0.14,0.8,0.8])

for i in range(nfiles): ax.plot(xdata[i],ydata[i],colors[np.mod(i,7)],label=legends[i])

ax.legend(loc=2, edgecolor="grey")
ax.set_xlim(min([min(xdata[i]) for i in range(nfiles)]),max([max(xdata[i]) for i in range(nfiles)]))
ax.set_ylim(0)
ax.set_xlabel(str.capitalize(labels[0]["xlabel"])+" [eV]", fontsize=24)
ax.set_ylabel(str.capitalize(labels[0]["ylabel"]), fontsize=24)

fig.subplots_adjust(left=None, bottom=None, right=None, wspace=None, hspace=None)

plt.savefig('PLOT.ps',  orientation='portrait',format='eps')
plt.savefig('PLOT.png', orientation='portrait',format='png',dpi=dpipng)

if (showpyplot): plt.show()
#-------------------------------------------------------------------------------
