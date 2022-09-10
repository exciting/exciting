#!/usr/bin/env python
#%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%#
#%!% ---------------------------------- OPTIMIZE-analyze.py ---------------------------------- %!%#
#%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%#
#
# AUTHOR:
#    Rostam Golesorkhtabar
#    r.golesorkhtabar@gmail.com
# 
# DATE:
#    Wed Jan 01 00:00:00 2014
#
# SYNTAX:
#    python OPTIMIZE-analyze.py
#           OPTIMIZE-analyze.py
# 
# EXPLANATION:
#
#__________________________________________________________________________________________________

from   pylab import *
import os
import sys
import glob
import copy
import math
import os.path
import numpy as np
from   lxml  import etree as ET
import matplotlib.pyplot as plt
import pylab             as pyl
from   scipy.optimize import fmin_powell

import matplotlib.style
if matplotlib.__version__.split(".")[0]=="2": matplotlib.style.use('classic')

#%!%!%--- CONSTANTS ---%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!
_e   =  1.602176565e-19         # elementary charge
Bohr =  5.291772086e-11         # a.u. to meter
Ha2eV= 27.211396132             # Ha to eV
ToGPa= (_e*Ha2eV)/(1e9*Bohr**3) # Ha/[a.u.]^3 to GPa
#__________________________________________________________________________________________________

#%!%!%--- SUBROUTINS AND FUNCTIONS ---%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%
def E_eos(p0, V):
    if (eos=='M'):
        """ Murnaghan Energy"""
        E0, V0, B0, Bp = p0
        E = E0 + (B0*V/Bp*(1/(Bp-1)*(V0/V)**Bp +1)-B0*V0/(Bp-1))
    else:
        """ Birch-Murnaghan Energy"""
        E0, V0, B0, Bp = p0
        E = E0 + (9.*B0*V0/16)*(((((V0/V)**(2./3))-1.)**3.)*Bp \
               + ((((V0/V)**(2/3.))-1.)**2.)*(6.-4.*((V0/V)**(2./3.))))
    return E
#--------------------------------------------------------------------------------------------------

def P_eos(p0, V):
    if (eos=='M'):
        """ Murnaghan Pressure"""
        E0, V0, B0, Bp = p0
        P = B0/Bp*((V0/V)**Bp - 1.)
    else:
        """ Birch-Murnaghan Pressure"""
        E0, V0, B0, Bp = p0
        P = 3./2*B0*((V0/V)**(7./3) - (V0/V)**(5./3))*(1. + 3./4*(Bp-4.)*((V0/V)**(2./3) - 1.))
    return P
#--------------------------------------------------------------------------------------------------

def snr(p0, v, e):
    """ Squared norm of residue vector calculation """
    return np.sum((e - E_eos(p0, v))**2.)
#--------------------------------------------------------------------------------------------------

def sortlist(lst1, lst2):
    temp = copy.copy(lst1)

    lst3 = []
    lst4 = []

    temp.sort()

    for i in range(len(lst1)):
        lst3.append(lst1[lst1.index(temp[i])])
        lst4.append(lst2[lst1.index(temp[i])])

    return lst3, lst4
#--------------------------------------------------------------------------------------------------

def readenergy():
    os.system("grep \"Total energy     \" INFO.OUT > tempfile") 
    tmpfile = open('tempfile', 'r')
    e = float(tmpfile.readlines()[-1].strip().split()[3])
    tmpfile.close()
    os.system("rm -f tempfile")
    return e
#__________________________________________________________________________________________________

#%!%!%--- Reading the INFO file ---%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!
INFO_file = str(glob.glob('INFO_*')[0])

INFO=open(INFO_file,'r')
mod =INFO_file[5:]

l1  = INFO.readline()
l2  = INFO.readline()

l3  = INFO.readline()
mdr = float(l3.split()[-1])

l4  = INFO.readline()
NoP = int(l4.split()[-1])

INFO.close()
#--------------------------------------------------------------------------------------------------

#%!%--- Reading the energies ---%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%
if (mod == 'VOL'):

    volume = []
    energy = []

    vollist= glob.glob('vol_??')
    for vol_num in vollist:
        os.chdir(vol_num)

        if (os.path.exists('INFO.OUT') == False):
            print'\n     ... Oops NOTICE: There is NO "INFO.OUT" file in "'+ vol_num + \
            '" directory !?!?!?    \n'

        for line in open('INFO.OUT','r'):
            if (line.find('Unit cell volume')>=0): 
                vol = float(line.split()[-1])
                break
        volume.append(vol)

        energy.append(readenergy())

        os.chdir('../')

    volume, energy = sortlist(volume, energy)
    
    fvol = open('energy-vs-volume', 'w')
    for i in range(len(energy)):
        print >>fvol, volume[i],'   ', energy[i]
    fvol.close()

    data = np.loadtxt('energy-vs-volume')
    vi, ei = data.T
    if (len(ei) < 3): sys.exit('\n     ... Oops ERROR: EOS fit needs at least 3 points.    \n')

    #%!%!%!%!--- Reading the EOS type ---%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!
    eos = raw_input('\n>>>> Murnaghan or Birch-Murnaghan EOS: [M/B] ').upper()
    if (eos != 'B' and eos != 'M'): sys.exit("\n    ... Oops ERROR: Choose 'B' or 'M' \n")
    if (eos == 'B'): eos = 'BM'
    #----------------------------------------------------------------------------------------------

    #%!%!%!%!--- FIT CALCULATIONS ---%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!
    a2, a1, a0 = np.polyfit(vi, ei, 2)
    V0 = -a1/(2.*a2)
    E0 = a2*V0**2. + a1*V0 + a0
    B0 = a2*V0
    Bp = 2.

    p0 = [E0, V0, B0, Bp]

    viei = sorted([zip(vi, ei)])
    v, e = np.array(viei).T

    p1, fopt, direc, n_iter, n_funcalls, warnflag = \
    fmin_powell(snr, p0, args=(v, e), full_output=True, disp=0)
    E0, V0, B0, Bp = p1

    print\
    '\n =====================================================================',\
    '\n Fit accuracy:',\
    '\n     Log(Final residue in [Ha]): '+str(round(log10(sqrt(fopt)),2)),'\n'\
    '\n Final parameters:'                                \
    '\n     E_min = ' + str(round(E0,7))      +' [Ha]'    \
    '\n     V_min = ' + str(round(V0,4))      +' [Bohr^3]'\
    '\n     B_0   = ' + str(round(B0*ToGPa,3))+' [GPa]'   \
    "\n     B'    = " + str(round(Bp,3))      +'\n'

    str_V = []
    str_de= []
    str_P = []

    for i in range(len(ei)):
        Pi     = P_eos(p1, vi[i])*ToGPa
        ei_eos = E_eos(p1, vi[i])

        str_vi = str(round(vi[i],4))

        if (Pi > 0):
            str_Pi = '+'+str(round(Pi,3))
        else:
            str_Pi =     str(round(Pi,3))

        dei = ei[i] - ei_eos
        if (dei > 0):
            str_dei = '+'+str('%8.8f'%(dei))
        else:
            str_dei =     str('%8.8f'%(dei))

        str_V.append(str_vi)
        str_de.append(str_dei)
        str_P.append(str_Pi)

    sum_Vi = 0
    sum_dei= 0
    sum_Pi = 0
    for i in range(len(ei)):
        sum_Vi = sum_Vi  + len(str_V[i])
        sum_dei= sum_dei + len(str_de[i])
        sum_Pi = sum_Pi  + len(str_P[i])

    len_Vi = int(sum_Vi /len(ei)) + 1
    len_dei= int(sum_dei/len(ei)) + 1
    len_Pi = int(sum_Pi /len(ei)) + 1

    #%!%!%--- WRITING THE OUTPUT FILE ---%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!
    outf = open(eos + '_eos.out', 'w')

    if (eos=='M'):
        print >>outf,' === Murnaghan eos ==============================='
    else:
        print >>outf,' === Birch-Murnaghan eos ========================='

    print >>outf,                                                               \
      ' Fit accuracy:',                                                         \
    '\n     Log(Final residue in [Ha]): '+str(round(log10(sqrt(fopt)),2)),  '\n'\
    '\n Final parameters:'                                                      \
    '\n     E_min = ' + str(round(E0,7))      +' [Ha]'                          \
    '\n     V_min = ' + str(round(V0,4))      +' [Bohr^3]'                      \
    '\n     B_0   = ' + str(round(B0*ToGPa,3))+' [GPa]'                         \
    "\n     B'    = " + str(round(Bp,3))      +                                 \
    '\n =================================================\n'\
    '\n Volume' + ((len_Vi-3)*' ') + 'E_dft-E_eos     Pressure [GPa]'

    for i in range(len(ei)):
        print >>outf, str_V[i]  + ((len_Vi -len(str_V[i]  ))*' ') + '    ' \
                    + str_de[i] + ((len_dei-len(str_de[i]))*' ') + '    ' + str_P[i]

    outf.close()
    #----------------------------------------------------------------------------------------------

    #%!%!%--- Writing the 'eos-optimized.xml' file ---%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%
    INOBJ= open(mod.lower()+'-xml/input.xml', 'r')
    doc  = ET.parse(INOBJ)
    root = doc.getroot()

    scale = map(float,doc.xpath('/input/structure/crystal/@scale'))
    if (scale==[]):
        ascale=1.
    else:
        ascale=scale[0]
        
    stretchstr = doc.xpath('/input/structure/crystal/@stretch')
    if (stretchstr==[]):
        stretch=[1.,1.,1.]
    else:
        stretch=np.array(map(float,stretchstr[0].split()))

    basevectsn = doc.xpath('//basevect/text()')
    bv = []
    for basevect in basevectsn:
        bv.append(map(float,basevect.split()))

    M_old= np.array(bv)
    D    = np.linalg.det(M_old)
    V0_in= abs(stretch[0]*stretch[1]*stretch[2]*ascale**3*D)
    s_min= (V0/V0_in)**(1./3.)-1.

    def_matrix={\
    'VOL'  :[[ 1.+s_min, 0.      , 0.      ],
             [ 0.      , 1.+s_min, 0.      ],
             [ 0.      , 0.      , 1.+s_min]]}

    M_min = np.array(def_matrix[mod])
    M_new = np.dot(M_old, M_min)

    bsvct = doc.xpath('//crystal/basevect')
    for j in range(3):
        bdummy = '%22.16f'%(M_new[j,0]) + '%22.16f'%(M_new[j,1]) + '%22.16f'%(M_new[j,2])+' '
        bsvct[j].text = bdummy

    #---Writing the structure file-----------------------------------------------------------------
    OUTOBJ = open(eos+'-optimized.xml', 'w')
    OUTOBJ.write(ET.tostring(root, method         ='xml',
                                   pretty_print   =True ,
                                   xml_declaration=False ,
                                   encoding       ='UTF-8'))
    OUTOBJ.close()
    print   ' Optimized lattice parameter saved into the file: "' + eos + '-optimized.xml".'\
          '\n =====================================================================\n'
#--------------------------------------------------------------------------------------------------
if (mod != 'VOL'):

    fee  = open('energy-vs-strain', 'w')

    for i in range(1, NoP+1):
        if (0  < i and i <   10): dir_num = mod.lower() + '_0'+ str(i)
        if (9  < i and i <  100): dir_num = mod.lower() + '_' + str(i)

        if (os.path.exists(dir_num) == False):
            print '\n    ... Oops NOTICE: There is NO '+ dir_num +' directory !?!?!?    \n'
            break

        os.chdir(dir_num)

        if (os.path.exists('INFO.OUT') == False):
            print '\n     ... Oops NOTICE: There is NO "INFO.OUT" file in "'+ dir_num + \
            '" directory !?!?!?    \n'

        s = i-(NoP+1)/2
        r = 2*mdr*s/(NoP-1)
        if (s==0): r=0.00001

        if (r>0):
            strain ='+'+str(round(r,10))
        else:
            strain =    str(round(r,10))

        print >>fee, strain,'   ', readenergy()
        os.chdir('../')

    fee.close()

    data = np.loadtxt('energy-vs-strain')
    si, ei = data.T
    vs = sorted([zip(si, ei)])
    s, e = np.array(vs).T
    if (len(e) < 5):
        sys.exit('\n     ... Oops ERROR: 4th order polynomial fit needs at least 5 points.\n')

    #%!%!%!%!--- FIT CALCULATIONS ---%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!
    coeffitions = np.polyfit(si, ei, 4)
    f4 = np.poly1d(coeffitions)

    s_fit = np.linspace(mdr*-1.2, mdr*1.2, 1000)
    e_fit = f4(s_fit)
    s_min = s_fit[e_fit.argmin()]
    e_min = e_fit[e_fit.argmin()]

    #%!%--- Writing the 'mod-optimized.xml' file ---%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%
    INOBJ= open(mod.lower()+'-xml/input.xml', 'r')
    doc  = ET.parse(INOBJ)
    root = doc.getroot()

    stretchstr = doc.xpath('/input/structure/crystal/@stretch')
    if (stretchstr==[]):
        stretch=[1.,1.,1.]
    else:
        stretch=np.array(map(float,stretchstr[0].split()))

    basevectsn = doc.xpath('//basevect/text()')
    bv = []
    for basevect in basevectsn:
        bv.append(map(float,basevect.split()))

    M_old= np.array(bv)

    def_matrix={\
    'BOA'  :[[(1+s_min)**-.5, 0.          , 0.           ],
             [ 0.           , 1.+s_min    , 0.           ],
             [ 0.           , 0.          ,(1+s_min)**-.5]],\

    'COA'  :[[(1+s_min)**-.5, 0.            , 0.         ],
             [ 0.           , (1+s_min)**-.5, 0.         ],
             [ 0.           , 0.            , 1.+s_min   ]],\

    'ALPHA':[[1./(1-s_min**2), 0.           , 0.         ],
             [ 0.            , 1.           ,s_min       ],
             [ 0.            ,s_min         , 1.         ]],\

    'BETA' :[[ 1.          , 0.             ,s_min       ],
             [ 0.          , 1./(1-s_min**2), 0.         ],
             [s_min        , 0.             , 1.         ]],\

    'GAMMA':[[ 1.          ,s_min       , 0.             ],
             [s_min        , 1.         , 0.             ],
             [ 0.          , 0.         , 1./(1-s_min**2)]]}

    M_min = np.array(def_matrix[mod])
    M_new = np.dot(M_old, M_min)

    bsvct = doc.xpath('//crystal/basevect')
    for j in range(3):
        bdummy = '%22.16f'%(M_new[j,0]) + '%22.16f'%(M_new[j,1]) + '%22.16f'%(M_new[j,2])+' '
        bsvct[j].text = bdummy

    #---Writing the structure file-----------------------------------------------------------------
    OUTOBJ   = open(mod.lower()+'-optimized.xml', 'w')
    OUTOBJ.write(ET.tostring(root, method         ='xml',
                                   pretty_print   =True ,
                                   xml_declaration=False ,
                                   encoding       ='UTF-8'))
    OUTOBJ.close()

    print '\n ====================================================================='\
          '\n Optimized lattice parameter saved into the file: "'+ mod.lower() +'-optimized.xml".'\
          '\n =====================================================================\n'
#--------------------------------------------------------------------------------------------------

#%!%!%--- PLOT DEFINITIONS ---%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%
params = {'axes.linewidth'        : 2.,
          'figure.subplot.bottom' : 0.14,
          'figure.subplot.right'  : 0.93, 
          'figure.subplot.left'   : 0.20,
          'xtick.major.pad'       : 8, 
          }

plt.rcParams.update(params)

#%!%!%--- PLOT SECTION ---%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%
fig = plt.figure()
ax  = fig.add_subplot(111)
fig.subplots_adjust(left=0.20)

if (mod == 'VOL'):

    if (eos=='M'):
        fit_label = 'Murnaghan eos'
    else:
        fit_label = 'Birch-Murnaghan eos'

    xlabel = u'Volume [Bohr\u00B3]'
    ylabel = 'Energy [Ha]'

    plt.text(0.3,0.75, 'E$_{min}$ = '+str(round(E0,7))+' [Ha]'         , transform = ax.transAxes)
    plt.text(0.3,0.70, 'V$_{min}$ = '+str(round(V0,4))+u' [Bohr\u00B3]', transform = ax.transAxes)
    plt.text(0.3,0.65, 'B$_0$ = '+str(round(B0*ToGPa,3))+' [GPa]'      , transform = ax.transAxes)
    plt.text(0.3,0.60, 'B$^\prime$ = '+str(round(Bp,3))                , transform = ax.transAxes)

    vmn  = min(min(vi), V0)
    vmx  = max(max(vi), V0)
    dv   = vmx - vmn
    v_eos= np.linspace(vmn-(0.1*dv), vmx+(0.1*dv), 1000)
    e_eos= E_eos(p1, v_eos)

    xx = [] ; xx = v_eos
    yy = [] ; yy = e_eos
    x0 = [] ; x0 = vi
    y0 = [] ; y0 = ei

if (mod != 'VOL'):
    xlabel    = 'Physical strain $\epsilon$'
    ylabel    = 'Energy [Ha]'
    fit_label = '4th order polynomial fit'

    plt.text(0.3,0.75, 'E$_{min}$ = '+str(round(e_min, 7))+' [Ha]', transform = ax.transAxes)
    plt.text(0.3,0.70, '$\epsilon_{min}$  = '+str(round(s_min,5)) , transform = ax.transAxes)

    xx = [] ; xx = s_fit
    yy = [] ; yy = e_fit
    x0 = [] ; x0 = si
    y0 = [] ; y0 = ei    

ax.set_xlabel(xlabel, fontsize = 18)
ax.set_ylabel(ylabel, fontsize = 18)

ax.plot(xx, yy, 'k'               ,
                color     = 'red' ,
                linewidth = 2     ,
                label     = fit_label)

ax.plot(x0, y0, 'o'                     ,
                color          = 'green',
                markersize     = 8      ,
                markeredgecolor= 'black',
                markeredgewidth= 1      ,
                label          = 'DFT Calc.')
ax.legend(numpoints=1,loc=9)

for label in ax.xaxis.get_ticklabels(): label.set_fontsize(15)
for label in ax.yaxis.get_ticklabels(): label.set_fontsize(15)
for line in ax.get_xticklines() + ax.get_yticklines():
    line.set_markersize(6)
    line.set_markeredgewidth(2)

pyl.grid(True)

ax.xaxis.set_major_locator(MaxNLocator(7))

max_y = max(max(yy), max(y0))
min_y = min(min(yy), min(y0))

max_x = max(max(xx), max(x0))
min_x = min(min(xx), min(x0))

dyy = (max_y-min_y)/15
ax.set_ylim(min_y-dyy,max_y+dyy)
dxx = (max_x-min_x)/18
ax.set_xlim(min_x-dxx,max_x+dxx)

if (mod == 'VOL'):
    plt.savefig(eos+'_eos.png', orientation='portrait',format='png',dpi=300)
    plt.savefig(eos+'_eos.eps', orientation='portrait',format='eps')
else:
    plt.savefig(mod.lower()+'.png', orientation='portrait',format='png',dpi=300)
    plt.savefig(mod.lower()+'.eps', orientation='portrait',format='eps')

#--------------------------------------------------------------------------------------------------