#!/usr/bin/python -u
# -*- coding: utf-8 -*-
#_______________________________________________________________________________

from   subprocess import *
from   lxml       import etree as ET
from   scipy.optimize.optimize import fmin
from   scipy.optimize.optimize import fmin_powell
import numpy as np
import time
import sys
import os
       
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def run_exciting(run_label):
    
    # DUNE pas-style 
    #cline = "cd workdir ; $HOME/JOB.sh "+run_label+" ; qsub job.job ; cd ../"
    
    # DUNE general style
    #cline = "cd workdir ; qsub runscript.sh ; cd ../"
    
    # LOCAL workstation
    cline = "cd workdir; $EXCITINGROOT/bin/excitingser ; cd ../"

    os.system(cline)
    
    return
  
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def read_input():
  
    # check if inputfile exists
    lfile = os.path.exists("input_ini.xml")
    if (not(lfile)): sys.exit("\n ERROR: file "+inputfile+" not found!\n")
    tree  = ET.parse(inputfile)

    return tree
  
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def get_minimum(efile):
  
    inpf  = open(efile,"r")
    lines = inpf.readlines()
    imin  = 1
    emin  = 0.
    for iline in range(len(lines)):
        line = float(lines[iline].split()[1])
        lmin = int(lines[iline].split()[0])
        if (line < emin):
            emin = line   
            imin = lmin
    inpf.close()
    
    if (imin < 1000): 
        imt = '%3i'
        dirmin = "rundir-"+(imt%imin)
    if (imin < 100): 
        imt = '%2i'
        dirmin = "rundir-"+(imt%imin)
    if (imin < 10): 
        imt = '%1i'
        dirmin = "rundir-0"+(imt%imin)

    return dirmin

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def generate_input(eps,tree):
 
    root = tree.getroot()
    basevect_txt = tree.xpath('//basevect/text()')
    basevect = []
    for bv in basevect_txt: basevect.append(map(float,bv.split()))
        
    M_old= np.array(basevect)
    
    def_matrix={\
    'VOL'  :[[1.+eps[0]      , 0.          , 0.             ],
             [0.          , 1.+eps[0]      , 0.             ],
             [0.          , 0.          , 1+eps[0]          ]],\

    'BOA'  :[[(1+eps[1])**-.5, 0.          , 0.             ],
             [ 0.         , 1.+eps[1]      , 0.             ],
             [ 0.         , 0.          ,(1+eps[1])**-.5    ]],\

    'COA'  :[[(1+eps[2])**-.5, 0.          , 0.             ],
             [ 0.         , (1+eps[2])**-.5, 0.             ],
             [ 0.         , 0.          , 1.+eps[2]         ]],\

    'ALPHA':[[1./(1-eps[3]**2), 0.           , 0.           ],
             [ 0.          , 1.           ,eps[3]           ],
             [ 0.          ,eps[3]           , 1.           ]],\

    'BETA' :[[ 1.          , 0.           ,eps[4]           ],
             [ 0.          , 1./(1-eps[4]**2), 0.           ],
             [eps[4]          , 0.           , 1.           ]],\

    'GAMMA':[[ 1.          ,eps[5]           , 0.           ],
             [eps[5]          , 1.           , 0.           ],
             [ 0.          , 0.           , 1./(1-eps[5]**2)]]}

    M_VOL   = np.array(def_matrix['VOL'])   ; M_new   = np.dot(M_old, M_VOL)
    M_BOA   = np.array(def_matrix['BOA'])   ; M_new   = np.dot(M_new, M_BOA)
    M_COA   = np.array(def_matrix['COA'])   ; M_new   = np.dot(M_new, M_COA)
    M_ALPHA = np.array(def_matrix['ALPHA']) ; M_new   = np.dot(M_new, M_ALPHA)
    M_BETA  = np.array(def_matrix['BETA'])  ; M_new   = np.dot(M_new, M_BETA)
    M_GAMMA = np.array(def_matrix['GAMMA']) ; M_new   = np.dot(M_new, M_GAMMA)
    
    fmt = '%22.16f'
    basevect_new = tree.xpath('//crystal/basevect')
    for j in range(3):
        basevect_new[j].text =\
	                fmt%(M_new[j,0]) + fmt%(M_new[j,1]) + fmt%(M_new[j,2])+' '

    # write new input.xml
    OUTOBJ   = open('workdir/input.xml', 'w')
    OUTOBJ.write(ET.tostring(root, method         ='xml',
                                   pretty_print   =True ,
                                   xml_declaration=False ,
                                   encoding       ='UTF-8'))
    OUTOBJ.close()

    os.system("INPUT-relaxupdate.py workdir/input.xml geometry_opt.xml")
    os.system('mv -f workdir/input_rel.xml workdir/input.xml')
    
    return
  
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def get_energy(eps):

    # get total energy
    string_one = '''grep "Total energy" workdir/INFO.OUT|\
                    tail -n 1|awk '{print $4}' '''
    string_two = '''grep "Unit cell volume" workdir/INFO.OUT|\
                    awk '{print $5}' '''
    
    tot_energy = float(Popen(string_one, shell=True,\
                             stdout=PIPE).communicate()[0].strip())
    volume     = float(Popen(string_two, shell=True,\
                             stdout=PIPE).communicate()[0].strip())
         
    return tot_energy
  
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
def tot_energy(eps):
  
    os.system("rm -rf workdir")
    os.system("mkdir workdir")
    
    #read scale factor 
    inpf  = open('scale',"r")
    scale = float(inpf.readline().split()[0])
    inpf.close()
    
    #scale eps because first step in fmin is to small
    fac   = [x*scale for x in eps]
    
    e = [0,0,0,0,0,0]

    if len(eps)==1: e[0]=fac[0]
    if len(eps)==2: e[0]=fac[0]; e[2]=fac[1]
    if len(eps)==3: e[0]=fac[0]; e[1]=fac[1]; e[2]=fac[2]
    if len(eps)==4: 
       # check angle for monoclinic structures
       inpf = open('monoclinic_angle',"r")
       mc_angle = int(inpf.readline().split()[0])
       e[0]=fac[0]; e[1]=fac[1]; e[2]=fac[2]; e[mc_angle]=fac[3]
    if len(eps)==6: e[:]=fac[:]

    tree = read_input()
    generate_input(e,tree)

    # run exciting
    ldone = False
    run_exciting("run_label")
     
    # wait until the calculation is completed
    string_zero = '''grep "$(head -n 2 workdir/INFO.OUT|tail -n 1|\
                     awk '{print $2" "$3" stopped"}')" workdir/INFO.OUT|wc -l'''
    while True:
        time.sleep(5)
        if (os.path.exists('workdir/INFO.OUT')):
            ldone = bool(int(Popen(string_zero, shell=True,\
	                     stdout=PIPE).communicate()[0].strip()))
        if ldone: break
    
    # get total energy
    energy = get_energy(e)
    
    gfile = "workdir/geometry_opt.xml"
    lgfile = os.path.exists(gfile)
    if (lgfile): os.system('cp workdir/geometry_opt.xml ./')

    os.system("LATTICE-stepdirectory.sh")
    
    return energy

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def check_monoclinic(tree):

    stretch  = [1.,1.,1.]
    basevect = []
    
    basevect_txt = tree.xpath('//basevect/text()')
    for bv in basevect_txt: basevect.append(map(float,bv.split()))  
    BV= np.array(basevect)

    xml_stretch = tree.xpath('//crystal/@stretch')
    if (len(xml_stretch) > 0): 
        xml_stretch[0] = xml_stretch[0].replace('d', 'e')
        xml_stretch[0] = xml_stretch[0].replace('D', 'e')
        stretch = nu.array(map(float,xml_stretch[0].split()))
 
    epsangle = 1.e-8
    alphascalarproduct = 0.0
    betascalarproduct  = 0.0
    gammascalarproduct = 0.0

    for icar in range(3):
        alphascalarproduct = alphascalarproduct\
                           + BV[1,icar]*BV[2,icar]*stretch[1]*stretch[2]
        betascalarproduct  = betascalarproduct\
                           + BV[0,icar]*BV[2,icar]*stretch[0]*stretch[2]
        gammascalarproduct = gammascalarproduct\
                           + BV[0,icar]*BV[1,icar]*stretch[0]*stretch[1]
    a = alphascalarproduct
    b = betascalarproduct
    g = gammascalarproduct
    e = epsangle

    conventional = False
    all_right    = False
    mc_angle     = 6
    
    if ( (a< e) and (b< e) and (g>=e) ):  
        conventional = True ; mc_angle  = 5  
    if ( (a< e) and (b>=e) and (g< e) ):  
        conventional = True ; mc_angle  = 4
    if ( (a>=e) and (b< e) and (g< e) ):  
        conventional = True ; mc_angle  = 3
    if ( (a< e) and (b< e) and (g< e) ):  
        conventional = True ; all_right = True 

    if (not(conventional)):
        sys.exit(\
        "\n     ... Oops ERROR: Your MONOCLINIC structure is not in the conventional cell."+\
        "\n                     This is not compatible with LATTICE-optimization.py!\n"+\
        "\n                     Please, CHECK your input file!\n")
    
    if (not(all_right)):
        outf = open("monoclinic_angle","w")
        print>>outf, '%2i'%mc_angle
        print>>outf, " 3 => alpha  "
        print>>outf, " 4 => beta   "
        print>>outf, " 5 => gamma  "  
        outf.close()
    
    return all_right

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#START==========================================================================          
# reading input file name

narg  = len(sys.argv)-1
inputfile = "input.xml"
if (narg>0): inputfile = str(sys.argv[1])
linputfile = os.path.exists(inputfile)
if (not(linputfile)): 
    sys.exit('\n     ... Oops ERROR: There is NO '+inputfile+' file !?!?!? \n')
os.system('cp '+inputfile+' input_ini.xml')

eps_tolerance = 1.0e-4
ene_tolerance = 1.0e-5
scale         = 30.e0

if (narg>1): eps_tolerance = float(sys.argv[2])
if (narg>2): ene_tolerance = float(sys.argv[3])
if (narg>3): scale         = float(sys.argv[4])

outf = open("scale","w")  ;  print>>outf, scale  ;  outf.close()

#-------------------------------------------------------------------------------          
# create dictionary of symmetries in the Laue classification

LC_Dic = {              \
          'CI' :'Cubic I'        ,\
          'CII':'Cubic II'       ,\
          'HI' :'Hexagonal I'    ,\
          'HII':'Hexagonal II'   ,\
          'RI' :'Rhombohedral I' ,\
          'RII':'Rhombohedral II',\
          'TI' :'Tetragonal I'   ,\
          'TII':'Tetragonal II'  ,\
          'O'  :'Orthorhombic'   ,\
          'M'  :'Monoclinic'     ,\
          'N'  :'Triclinic'}

#-------------------------------------------------------------------------------            
# calculating space-group number and classifying it

os.system('$EXCITINGTOOLS/exciting2sgroup.py '+inputfile+' sgroup.in')
os.system('sgroup sgroup.in 1>sgroup.out 2>sgroup.err')
os.system('rm -f sgroup.in ')
os.system('touch exciting')
os.system('touch initial-step')

if (os.path.getsize('sgroup.err') != 0):
    fer  = open('sgroup.err', 'r')
    lines= fer.readlines()
    print '\n     ... Oops '+ lines[0]
    for i in range(1, len(lines)):
        print '                 '+ lines[i]
    print
    fer.close()
    sys.exit()
else:
    os.system('rm -f sgroup.err')

SGfile  = open('sgroup.out', 'r')
SGlines = SGfile.readlines()
SGfile.close()

os.system('rm -f sgroup.out ')

for i in range(len(SGlines)):
    if (SGlines[i].find('Number and name of space group:') >= 0):
        SGN = int(float(SGlines[i].split()[6]))
        SGN_comment = SGlines[i].strip()
        break

if  (   1 <= SGN and SGN <=   2 ):   LC = 'N'    # Triclinic
elif(   3 <= SGN and SGN <=  15 ):   LC = 'M'    # Monoclinic
elif(  16 <= SGN and SGN <=  74 ):   LC = 'O'    # Orthorhombic
elif(  75 <= SGN and SGN <=  88 ):   LC = 'TII'  # Tetragonal II
elif(  89 <= SGN and SGN <= 142 ):   LC = 'TI'   # Tetragonal I
elif( 143 <= SGN and SGN <= 148 ):   LC = 'RII'  # Rhombohedral II
elif( 149 <= SGN and SGN <= 167 ):   LC = 'RI'   # Rhombohedral I
elif( 168 <= SGN and SGN <= 176 ):   LC = 'HII'  # Hexagonal II
elif( 177 <= SGN and SGN <= 194 ):   LC = 'HI'   # Hexagonal I
elif( 195 <= SGN and SGN <= 206 ):   LC = 'CII'  # Cubic II
elif( 207 <= SGN and SGN <= 230 ):   LC = 'CI'   # Cubic I
else: sys.exit('\n     ... Oops ERROR: WRONG Space-Group Number !?!?!?\n')

print
print "=============================================================================="
print '\n     '+ SGN_comment +'\
       \n     '+ LC_Dic[LC] +' structure in the Laue classification.\n'
print "=============================================================================="
print
       
#-------------------------------------------------------------------------------          
# read optimization type                          

if (LC=='CI' or\
    LC=='CII'):  opt_type = ['VOL',]

if (LC=='HI' or\
    LC=='HII'or\
    LC=='RI' or\
    LC=='RII'or\
    LC=='TI' or\
    LC=='TII'):  opt_type = ['VOL','COA']
    
if (LC=='O'):    opt_type = ['VOL','BOA','COA']
if (LC=='M'):    opt_type = ['VOL','BOA','COA','GAMMA']
if (LC=='N'):    opt_type = ['VOL','BOA','COA','ALPHA','BETA','GAMMA']

#-------------------------------------------------------------------------------          
# check compatibility for monoclinic structures

if (LC=='M'): 
    tree = read_input()
    if (check_monoclinic(tree)): opt_type = ['VOL','BOA','COA']

#-------------------------------------------------------------------------------          
# performing optimization step                          

print "     a          b          c          alpha      beta       gamma      delta_E"
print "------------------------------------------------------------------------------"

eps0 = [0,] ; eps0 = len(opt_type)*eps0

res = fmin(tot_energy,\
           eps0,\
           xtol=eps_tolerance,\
           ftol=ene_tolerance,\
           retall=0,\
           disp=0,\
           full_output=0)
print
print "=============================================================================="
print "     Convergence has been achieved."
print
#
#
dirmin = get_minimum("energy-vs-step")
os.system("INPUT-relaxupdate.py "+dirmin+"/input.xml "+dirmin+"/geometry_opt.xml")
os.system('cp '+dirmin+'/input_rel.xml ./input_opt_rel.xml')
#
print "=============================================================================="
print "INITIAL LATTICE PARAMETERS:"

os.system("LATTICE-parameters.sh "+inputfile) 
#read scale factor 
inpf  = open('scale',"r")
scale = float(inpf.readline().split()[0])
inpf.close()

print "=============================================================================="
print "OPTIMIZED LATTICE PARAMETERS (E_TOL = "+str(ene_tolerance)+" Ha;"\
                                  +" X_TOL = "+str(eps_tolerance)+";"\
                                  +" X_FAC = "+str(scale)+"):"
os.system("LATTICE-parameters.sh input_opt_rel.xml")    
print "=============================================================================="
print
os.system('rm -f scale')

#END============================================================================          
         

