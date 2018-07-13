#!/usr/bin/env python
#%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%#
#%!%!% --------------------------------- OPTIMIZE-setup.py --------------------------------- %!%!%#
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
#    python OPTIMIZE-setup.py
#           OPTIMIZE-setup.py
# 
# EXPLANATION:
# 
#__________________________________________________________________________________________________

import os
import sys
import math
import shutil
import os.path
import numpy as np
from sys  import stdin
from lxml import etree as ET

#%!%!%--- Dictionaries ---%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%
Def_matrix = {              \
'VOL'  :
'\n[ 1+eps  0      0     ]' \
'\n[ 0      1+eps  0     ]' \
'\n[ 0      0      1+eps ]',\

'BOA'  :
'\n[(1+eps)^-.5   0           0          ]' \
'\n[ 0           (1+eps)^+1.  0          ]' \
'\n[ 0            0          (1+eps)^-.5 ]',\

'COA'  :
'\n[(1+eps)^-.5   0           0          ]' \
'\n[ 0           (1+eps)^-.5  0          ]' \
'\n[ 0            0          (1+eps)^+1. ]',\

'ALPHA':
'\n[ 1/(1-eps^2)  0           0          ]' \
'\n[ 0            1          eps         ]' \
'\n[ 0           eps          1          ]',\

'BETA' :
'\n[ 1           0           eps         ]' \
'\n[ 0           1/(1-eps^2)  0          ]' \
'\n[eps          0            1          ]',\

'GAMMA':
'\n[ 1          eps           0          ]' \
'\n[eps          1            0          ]' \
'\n[ 0           0            1/(1-eps^2)]'}

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
#--------------------------------------------------------------------------------------------------

#%!%!%--- Checking the input file exist ---%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!
INF = raw_input('\n>>>> Please enter the exciting input file name: ')
if (os.path.exists(INF) == False):
    sys.exit('\n     ... Oops ERROR: There is NO '+ INF +' file !?!?!?    \n')
#--------------------------------------------------------------------------------------------------

#%!%!%--- Calculating the Space-Group Number and classifying it ---%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!
os.system('$EXCITINGTOOLS/exciting2sgroup.py '+ INF +' sgroup.in')
os.system('sgroup sgroup.in 1>sgroup.out 2>sgroup.err; rm -f sgroup.in ')

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

SGf   = open('sgroup.out', 'r')
SGlins= SGf.readlines()
SGf.close()

for i in range(len(SGlins)):
    if (SGlins[i].find('Number and name of space group:') >= 0):
        SGN = int(float(SGlins[i].split()[6]))
        SGN_explanation=SGlins[i].strip()
        break

if (1 <= SGN and SGN <= 2):      # Triclinic
    LC = 'N'

elif(3 <= SGN and SGN <= 15):    # Monoclinic
    LC = 'M'

elif(16 <= SGN and SGN <= 74):   # Orthorhombic
    LC = 'O'

elif(75 <= SGN and SGN <= 88):   # Tetragonal II
    LC = 'TII'
  
elif(89 <= SGN and SGN <= 142):  # Tetragonal I
    LC = 'TI'

elif(143 <= SGN and SGN <= 148): # Rhombohedral II
    LC = 'RII'

elif(149 <= SGN and SGN <= 167): # Rhombohedral I
    LC = 'RI'

elif(168 <= SGN and SGN <= 176): # Hexagonal II
    LC = 'HII'

elif(177 <= SGN and SGN <= 194): # Hexagonal I
    LC = 'HI'

elif(195 <= SGN and SGN <= 206): # Cubic II
    LC = 'CII'

elif(207 <= SGN and SGN <= 230): # Cubic I
    LC = 'CI'

else: sys.exit('\n     ... Oops ERROR: WRONG Space-Group Number !?!?!?    \n')

print '\n     '+ SGN_explanation +'\
       \n     '+ LC_Dic[LC] +' structure in the Laue classification.'

#--------------------------------------------------------------------------------------------------



























#%!%!%--- Reading the input file ---%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%
INOBJ= open(INF, 'r')
doc  = ET.parse(INOBJ)
root = doc.getroot()
#--------------------------------------------------------------------------------------------------

#%!%!%--- Reading the scale, stretch, and base vectors from the input file and calculate the volume
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
V0   = abs(stretch[0]*stretch[1]*stretch[2]*ascale**3*D)

#--------------------------------------------------------------------------------------------------
# Check compatibility for monoclinic structures

if (LC=='M'):
  
    #print M_old
  
    epsangle = 1.e-8
    gammascalarproduct = 0.0

    for icar in range(3):
        #print M_old[0,icar], M_old[1,icar]
        gammascalarproduct = gammascalarproduct + M_old[0,icar]*M_old[1,icar]*stretch[0]*stretch[1]

    #print gammascalarproduct 

    if (abs(gammascalarproduct) < epsangle):
        sys.exit("\n     ... Oops ERROR: Your MONOCLINIC structure is not compatible"+\
                 "\n                     with the OPTIMIZE internal representation, where"+\
                 "\n                     the angle GAMMA (between bvec_1 and bvec_2) is the"+\
                 "\n                     ONLY non right angle between the crystal basis vectors!\n"+\
                 "\n                     Please, CHECK your input file!\n")

#--------------------------------------------------------------------------------------------------

#%!%!%--- Read optimization type ---%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%
if (LC=='CI' or LC=='CII'):
    dirn = 'VOL'
    print

if (LC=='HI' or\
    LC=='HII'or\
    LC=='RI' or\
    LC=='RII'or\
    LC=='TI' or\
    LC=='TII'):
    print '\n     Which parameter would you like to optimize?'\
          '\n     1 ... volume                               '\
          '\n     2 ... c/a ratio with constant volume       '
    num = input(">>>> Please choose '1' or '2': ")
    if (num != 1 and num != 2):
        sys.exit("\n     ... Oops ERROR: Choose '1' or '2' \n")
    if (num == 1 ): dirn = 'VOL'
    if (num == 2 ): dirn = 'COA'

if (LC=='O'):
    print '\n     Which parameter would you like to optimize?'\
          '\n     1 ... volume                               '\
          '\n     2 ... b/a ratio with constant volume       '\
          '\n     3 ... c/a ratio with constant volume       '
    num = input(">>>> Please choose '1' or '2' or '3': ")
    if (num != 1 and num != 2 and num != 3):
        sys.exit("\n     ... Oops ERROR: Choose '1' or '2' or '3'\n")
    if (num == 1 ): dirn = 'VOL'
    if (num == 2 ): dirn = 'BOA'
    if (num == 3 ): dirn = 'COA'

if (LC=='M'):
    print '\n     Which parameter would you like to optimize?'\
          '\n     1 ... volume                               '\
          '\n     2 ... b/a ratio with constant volume       '\
          '\n     3 ... c/a ratio with constant volume       '\
          '\n     4 ... gamma angle with constant volume     '
    num = input(">>>> Please choose '1' or '2' or '3' or '4': ")
    if (num != 1 and num != 2 and num != 3 and num != 4):
        sys.exit("\n     ... Oops ERROR: Choose '1' or '2' or '3' or '4'\n")
    if (num == 1 ): dirn = 'VOL'
    if (num == 2 ): dirn = 'BOA'
    if (num == 3 ): dirn = 'COA'
    if (num == 4 ): dirn = 'GAMMA'

if (LC=='N'):
    print '\n     Which parameter would you like to optimize?'\
          '\n     1 ... volume                               '\
          '\n     2 ... b/a ratio with constant volume       '\
          '\n     3 ... c/a ratio with constant volume       '\
          '\n     4 ... alpha angle with constant volume     '\
          '\n     5 ... beta  angle with constant volume     '\
          '\n     6 ... gamma angle with constant volume     '
    num = input(">>>> Please choose '1' or '2' or '3' or '4' or '5' or '6': ")
    if (num != 1 and num != 2 and num != 3 and num != 4 and num != 5 and num != 6):
        sys.exit("\n     ... Oops ERROR: Choose '1' or '2' or '3' or '4' or '5' or '6'\n")
    if (num == 1 ): dirn = 'VOL'
    if (num == 2 ): dirn = 'BOA'
    if (num == 3 ): dirn = 'COA'
    if (num == 4 ): dirn = 'ALPHA'
    if (num == 5 ): dirn = 'BETA'
    if (num == 6 ): dirn = 'GAMMA'
#--------------------------------------------------------------------------------------------------

#%!%!%--- Directory Management ---%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%
if (os.path.exists(dirn+'_old')):
    shutil.rmtree( dirn+'_old')

if (os.path.exists(dirn)):
    os.rename(dirn,dirn+'_old')

os.mkdir(dirn)
os.chdir(dirn)
#--------------------------------------------------------------------------------------------------

#%!%!%--- Reading the maximum physical strain and number of distorted structures ---%!%!%!%!%!%!%!%
mdr = input('\n>>>> Please enter the maximum physical strain value  '\
            '\n     The suggested value is between 0.001 and 0.050: ')

if (1 < mdr or mdr < 0):
    sys.exit('\n     ... Oops ERROR: The maximum physical strain is OUT of range !!!!!!\n')

mdr = round(mdr, 4)
print '     The maximum physical strain is '+ str(mdr)

NoP = input('\n>>>> Please enter the number of the distorted structures [odd number > 4]: ')
NoP = int(abs(NoP))

if (NoP < 5):
    sys.exit('\n     ... Oops ERROR: The NUMBER of the distorted structures < 5 !!!!!!    \n')
if (99 < NoP):
    sys.exit('\n     ... Oops ERROR: The NUMBER of the distorted structures > 99 !!!!!!   \n')

if (NoP%2 == 0):
    NoP   += 1
print '     The number of the distorted structures is '+ str(NoP) + '\n'

ptn = int((NoP-1)/2)

if (mdr/ptn <= 0.00001):
    sys.exit('     ... Oops ERROR: The interval of the strain values is < 0.00001'\
           '\n                     Choose a larger maximum physical strain       '\
           '\n                     or a less number of distorted structures.   \n')
#--------------------------------------------------------------------------------------------------

#%!%!%--- Making the INFO file ---%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%
INFO = open('INFO_'+ dirn , 'w')
print >>INFO, 'Space-group number              =', SGN       ,\
            '\nStructure type                  =', LC_Dic[LC],\
            '\nMaximum physical strain         =', mdr       ,\
            '\nNumber of distorted structures  =', NoP
INFO.close
#--------------------------------------------------------------------------------------------------

#%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%#
#%!%! ----------------------------------- Structures maker ---------------------------------- %!%!#
#%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%!%#

fdis = open(dirn.lower()+'-Parameters', 'w')
print >>fdis,dirn+', Deformation Matrix = ' + Def_matrix[dirn]

cont= 0
for s in range(-ptn, ptn+1):
    r = mdr*s/ptn
    if (s==0): r = 0.00001
    eps = r
    def_matrix={\
    'VOL'  :[[1.+eps      , 0.          , 0.             ],
             [0.          , 1.+eps      , 0.             ],
             [0.          , 0.          , 1+eps          ]],\

    'BOA'  :[[(1+eps)**-.5, 0.          , 0.             ],
             [ 0.         , 1.+eps      , 0.             ],
             [ 0.         , 0.          ,(1+eps)**-.5    ]],\

    'COA'  :[[(1+eps)**-.5, 0.          , 0.             ],
             [ 0.         , (1+eps)**-.5, 0.             ],
             [ 0.         , 0.          , 1.+eps         ]],\

    'ALPHA':[[1./(1-eps**2), 0.           , 0.           ],
             [ 0.          , 1.           ,eps           ],
             [ 0.          ,eps           , 1.           ]],\

    'BETA' :[[ 1.          , 0.           ,eps           ],
             [ 0.          , 1./(1-eps**2), 0.           ],
             [eps          , 0.           , 1.           ]],\

    'GAMMA':[[ 1.          ,eps           , 0.           ],
             [eps          , 1.           , 0.           ],
             [ 0.          , 0.           , 1./(1-eps**2)]]}
        
    M_eps = np.array(def_matrix[dirn])
    M_new = np.dot(M_old, M_eps)

    cont += 1
    if (0  < cont and cont <  10): dirn_num = dirn.lower() + '_0'+ str(cont)
    if (9  < cont and cont < 100): dirn_num = dirn.lower() + '_' + str(cont)

    bsvct = doc.xpath('//crystal/basevect')
    print>>fdis, '\n'+dirn_num+', Physical strain = '+ str(eps)
    for j in range(3):
        bdummy = '%22.16f'%(M_new[j,0]) + '%22.16f'%(M_new[j,1]) + '%22.16f'%(M_new[j,2])+' '
        print>>fdis, 'V' + str(j+1) + ' --=> ' + bdummy
        bsvct[j].text = bdummy

    os.mkdir(dirn_num)
    os.chdir(dirn_num)

    # Writing the structure file ------------------------------------------------------------------
    fileName = dirn_num + '.xml'
    OUTOBJ   = open(fileName, 'w')
    OUTOBJ.write(ET.tostring(root, method         ='xml',
                                   pretty_print   =True ,
                                   xml_declaration=False ,
                                   encoding       ='UTF-8'))
    OUTOBJ.close()
    os.chdir('../')
#--------------------------------------------------------------------------------------------------

print>>fdis,'\n   Distorted parameters: END'
fdis.close()

os.mkdir(dirn.lower()+'-xml')
os.system('cp -f */*_??.xml '+dirn.lower()+'-xml/')
os.system('cp -f ../'+INF+' '+dirn.lower()+'-xml/input.xml')
os.chdir('../')
#--------------------------------------------------------------------------------------------------