#!/usr/bin/env python
#%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%#
#%!%!%-------------------------------- STRESS-exciting-setup --------------------------------%!%!%#
#%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%#
#
# AUTHOR:
# Rostam Golesorkhtabar
# r.golesorkhtabar@gmail.com
# 
# DATE:
# Wed Jan 01 00:00:00 2014
#
# SYNTAX:
# python STRESS-exciting-setup
# 
# EXPLANATION:
# 
#__________________________________________________________________________________________________

from lxml  import etree as ET
from sys   import stdin
from numpy import *
import numpy as np
import os.path
import shutil
import glob
import math
import sys
import os

#%!%!%--- DICTIONARIS ---%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!
def_str_dic = {             \
'1':
'\n[ 1+eps  0      0     ]' \
'\n[ 0      1+eps  0     ]' \
'\n[ 0      0      1+eps ]',\

'2':
'\n[(1+eps)^-.5   0           0          ]' \
'\n[ 0           (1+eps)^+1.  0          ]' \
'\n[ 0            0          (1+eps)^-.5 ]',\

'3':
'\n[(1+eps)^-.5   0           0          ]' \
'\n[ 0           (1+eps)^-.5  0          ]' \
'\n[ 0            0          (1+eps)^+1. ]',\

'4':
'\n[ 1/(1-eps^2)  0           0          ]' \
'\n[ 0            1          eps         ]' \
'\n[ 0           eps          1          ]',\

'5':
'\n[ 1           0           eps         ]' \
'\n[ 0           1/(1-eps^2)  0          ]' \
'\n[eps          0            1          ]',\

'6':
'\n[ 1          eps           0          ]' \
'\n[eps          1            0          ]' \
'\n[ 0           0            1/(1-eps^2)]'}

Laue_dic = {            \
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
#--------------------------------------------------------------------------------------------------

#%!%!%--- Checking the *.xml file exist ---%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!
INF = raw_input('\n>>>> Please enter the exciting input file name: ')
if (os.path.exists(INF) == False):
    sys.exit('\n ... Oops ERROR: There is NO file with "'+ INF +'" name !?!?!?    \n')
#--------------------------------------------------------------------------------------------------

#%!%!%--- Calculate Space-Group Number and classify it ---%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%
os.system('xsltproc $EXCITINGSCRIPTS/exciting2sgroup.xsl '+ INF +' > sgroup.in')

SIN = open('sgroup.in', 'r')
l1  = SIN.readline()
l2  = SIN.readline()
gamma = float(l2.split()[-1])
beta  = float(l2.split()[-2])
alpha = float(l2.split()[-3])

os.system('sgroup sgroup.in 1>sgroup.out 2>sgroup.err')
os.system('rm -f sgroup.in ')

if (os.path.getsize('sgroup.err') != 0):
    fer  = open('sgroup.err', 'r')
    lines= fer.readlines()
    print '\n     ... Oops '+ lines[0]
    for i in range(1, len(lines)):
        print '                 '+ lines[i]
    fer.close()
    sys.exit()
else: os.system('rm -f sgroup.err')

SGf   = open('sgroup.out', 'r')
SGlins= SGf.readlines()
SGf.close()

for i in range(len(SGlins)):
    if (SGlins[i].find('Number and name of space group:') >= 0):
        SGN = int(float(SGlins[i].split()[6]))
        SGN_line = SGlins[i].strip()
        break

if (1 <= SGN and SGN <= 2):      # Triclinic
    LC = 'N'
    SCs= 6

elif(3 <= SGN and SGN <= 15):    # Monoclinic
    LC = 'M'
    SCs= 4

    if (SGN ==  5 or \
        SGN ==  8 or \
        SGN ==  9 or \
        SGN == 12 or \
        SGN == 15):

        print '\n ... Oops PRBLM: The monoclinic unique axis cannot be specified '
        unique_axis = raw_input('\n>>>> Please enter the monoclinic uique axis [a, b, or c]: ')
        if (unique_axis != 'a' or \
            unique_axis != 'b' or \
            unique_axis != 'c'):
            sys.exit('\n ... Oops ERROR: The monoclinic uique axis must be a, b, or c. !!!!!!\n')

    elif(alpha != 90):
        unique_axis = 'a'

    elif(beta != 90):
        unique_axis = 'b'

    else:
        unique_axis = 'c'

elif(16 <= SGN and SGN <= 74):   # Orthorhombic
    LC = 'O'
    SCs= 3

elif(75 <= SGN and SGN <= 88):   # Tetragonal II
    LC = 'TII'
    SCs= 2

elif(89 <= SGN and SGN <= 142):  # Tetragonal I
    LC = 'TI'
    SCs= 2

elif(143 <= SGN and SGN <= 148): # Rhombohedral II 
    LC = 'RII'
    SCs= 2

elif(149 <= SGN and SGN <= 167): # Rhombohedral I
    LC = 'RI'
    SCs= 2

elif(168 <= SGN and SGN <= 176): # Hexagonal II
    LC = 'HII'
    SCs= 2

elif(177 <= SGN and SGN <= 194): # Hexagonal I
    LC = 'HI'
    SCs= 2

elif(195 <= SGN and SGN <= 206): # Cubic II
    LC = 'CII'
    SCs= 1

elif(207 <= SGN and SGN <= 230): # Cubic I
    LC = 'CI'
    SCs= 1

else: sys.exit('\n ... Oops ERROR: WRONG Space-Group Number !?!?!?    \n')

if (LC != 'M'):
    print '\n     '+ SGN_line     +'                                      \
           \n     '+ Laue_dic[LC] +' structure in the Laue classification.\
           \n     This structure has '+ str(SCs) +' independent stress components.'

if (LC == 'M'):
    print '\n     '+ SGN_line[0:SGN_line.find('[')] +'                            \
           \n     '+ Laue_dic[LC] +' structure in the Laue classification.        \
           \n     This structure has '+ str(SCs) +' independent stress components.\
           \n     The monoclinic unique axis is "' + unique_axis + '".'

#--------------------------------------------------------------------------------------------------

#%!%!%--- Reading the maximum Physical strain ---%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!
mdr = input('\n>>>> Please enter the maximum amount of strain '\
            '\n     The suggested value is between 0.0010 and 0.0100 : ')

if (1 < mdr or mdr < 0):
    sys.exit('\n ... Oops ERROR: The maximum physical strain is out of range !!!!!!\n')

mdr = round(mdr, 3)
print '     The maximum amount of strain is '+ str(mdr) + '\n'
#--------------------------------------------------------------------------------------------------

#%!%!%--- Reading the number of the distorted structures ---%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%
NoP = input('>>>> Please enter the number of the distorted structures [odd number > 4]: ')
NoP = int(abs(NoP))

if (NoP < 5):
    sys.exit('\n ... Oops ERROR: The NUMBER of the distorted structures < 5 !!!!!!    \n')
if (99 < NoP):
    sys.exit('\n ... Oops ERROR: The NUMBER of the distorted structures > 99 !!!!!!   \n')

if (NoP%2 == 0):
    NoP   += 1
print '     The number of the distorted structures is '+ str(NoP) + '\n'

ptn = int((NoP-1)/2)

if (mdr/ptn <=  0.00001):
    sys.exit(' ... Oops ERROR: The interval strain value is < 0.00001'  \
           '\n                 Choose a larger maximum amount of strain'\
           '\n                 or a less number of distorted structures.\n')
#--------------------------------------------------------------------------------------------------

#%!%!%--- Reading the input file ---%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%
INOBJ= open(INF, 'r')
doc  = ET.parse(INOBJ)
root = doc.getroot()
#--------------------------------------------------------------------------------------------------

#%%-- Reading the scale, stretch, and base vectors from the input file and calculate the volume -%%
scale = map(float,doc.xpath('/input/structure/crystal/@scale'))
if (scale == []): scale = [1.]

stretchstr = doc.xpath('/input/structure/crystal/@stretch')
if (stretchstr == []): stretch = [1.,1.,1.]

basevectsn = doc.xpath('//basevect/text()')
bv = []
for basevect in basevectsn:
    bv.append(map(float,basevect.split()))

M_old= np.array(bv)
D    = np.linalg.det(M_old)
V0   = abs(stretch[0]*stretch[1]*stretch[2]*scale[0]**3*D)
#--------------------------------------------------------------------------------------------------

#%!%!%--- Writing the "INFO_ElaStic" file ---%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!
INFO = open('INFO_Stress','w')
print >>INFO,'Space-group number              =', SGN          ,\
           '\nVolume of equilibrium unit cell =', V0, '[a.u^3]',\
           '\nMaximum physical strain         =', mdr          ,\
           '\nNumber of distorted structures  =', NoP

if (LC == 'M'):
    print >>INFO,\
             'Monoclinic unique axis          =', unique_axis

INFO.close
#--------------------------------------------------------------------------------------------------

#%!%!%--- Directory management ---%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%
OLDlist = glob.glob('Dst??_old')
for Dstn_old in OLDlist:
    shutil.rmtree(Dstn_old)

Dstlist = glob.glob('Dst??')
for Dstn in Dstlist:
    os.rename(Dstn, Dstn+'_old')

if (os.path.exists('Structures_exciting_old')):
    shutil.rmtree( 'Structures_exciting_old')

if (os.path.exists('Structures_exciting')):
    os.rename(     'Structures_exciting', 'Structures_exciting_old')
#--------------------------------------------------------------------------------------------------

if (LC == 'CI' or \
    LC == 'CII'):
    def_list = ['1']

if (LC == 'HI' or \
    LC == 'HII'or \
    LC == 'RI' or \
    LC == 'RII'or \
    LC == 'TI' or \
    LC == 'TII'):
    def_list = ['1','3']

if (LC == 'O'):
    def_list = ['1','2','3']

if (LC == 'M'):
    if (unique_axis == 'a'): def_list = ['1','2','3','4']
    if (unique_axis == 'b'): def_list = ['1','2','3','5']
    if (unique_axis == 'c'): def_list = ['1','2','3','6']

if (LC == 'N'):
    def_list = ['1','2','3','4','5','6']

#%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%#
#%!%! ---------------------------------- Structures Maker ----------------------------------- !%!%#
#%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%#

fdis = open('Distorted_Parameters','w')
cont1= 0
for i in def_list:

    cont1  = cont1 + 1
    if (cont1 < 10):
        Dstn = 'Dst0'+str(cont1)
    else:
        Dstn = 'Dst' +str(cont1)
    
    os.mkdir(Dstn)
    os.chdir(Dstn)

    print>>fdis, Dstn + ', Deformation Matrix = ' + def_str_dic [i] + '\n'

    cont2 = 0
    for s in range(-ptn, ptn+1):
        r = mdr*s/ptn
        if (s==0): r = 0.00001
        eps = r

        def_mtx_dic = {                                       \
        '1' : [[1.+eps      , 0.          , 0.             ],
               [0.          , 1.+eps      , 0.             ],
               [0.          , 0.          , 1+eps          ]],\

        '2' : [[(1+eps)**-.5, 0.          , 0.             ],
               [ 0.         , 1.+eps      , 0.             ],
               [ 0.         , 0.          ,(1+eps)**-.5    ]],\

        '3' : [[(1+eps)**-.5, 0.          , 0.             ],
               [ 0.         , (1+eps)**-.5, 0.             ],
               [ 0.         , 0.          , 1.+eps         ]],\

        '4' : [[1./(1-eps**2), 0.           , 0.           ],
               [ 0.          , 1.           ,eps           ],
               [ 0.          ,eps           , 1.           ]],\

        '5' : [[ 1.          , 0.           ,eps           ],
               [ 0.          , 1./(1-eps**2), 0.           ],
               [eps          , 0.           , 1.           ]],\

        '6' : [[ 1.          ,eps           , 0.           ],
               [eps          , 1.           , 0.           ],
               [ 0.          , 0.           , 1./(1-eps**2)]]}

        M_eps = np.array(def_mtx_dic[i])
        M_new = np.dot(M_old, M_eps)

        cont2 = cont2 + 1
        if (cont2 < 10):
            Dstn_cont2 = Dstn+ '_0'+str(cont2)
        else:
            Dstn_cont2 = Dstn+ '_' +str(cont2)

        print>>fdis, Dstn_cont2 + ',  eps = ' + str(r)

        bsvct = doc.xpath('//crystal/basevect')
        for j in range(3):
            bdummy = '%22.16f'%(M_new[j,0]) + '%22.16f'%(M_new[j,1]) + '%22.16f'%(M_new[j,2])+' '
            print>>fdis, 'V' + str(j+1) + ' --=> ' + bdummy
            bsvct[j].text = bdummy
        print>>fdis

        os.mkdir(Dstn_cont2)
        os.chdir(Dstn_cont2)

        # Writing the structure file --------------------------------------------------------------
        fileName = Dstn_cont2 + '.xml'
        OUTOBJ   = open(fileName, 'w')
        OUTOBJ.write(ET.tostring(root, method         ='xml',
                                       pretty_print   =True ,
                                       xml_declaration=True ,
                                       encoding       ='UTF-8'))
        OUTOBJ.close()
        # -----------------------------------------------------------------------------------------
        os.chdir('../')
    os.chdir('../')
print>>fdis,'   Distorted parameters: END'
fdis.close()
os.system('mkdir Structures_exciting; cp Dst??/Dst??_??/Dst??_??.xml Structures_exciting/')
#--------------------------------------------------------------------------------------------------
