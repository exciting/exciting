#!/usr/bin/python
import sys
import numpy as np

if len(sys.argv[1:]) < 2:
    print('Error('+sys.argv[0]+'): Problem with input arguments!')
    print('Usage: SpectralFunction.dat <ik list>')
    sys.exit()

fname = sys.argv[1]
iklist = sys.argv[2:]

# get info about k-points
vkc = []
fh = open(fname, 'r')
line = fh.readline()
while line:
    if "# k-point" in line:
        line = line.split()
        vkc.append([float(line[9]), float(line[10]), float(line[11])])
        line = fh.readline().split()
        nst = len(line)-1
        nw = 0
        while len(line):
            nw = nw+1
            line = fh.readline().split()
    line = fh.readline()
fh.close()
nkpt = len(vkc)

# get info abour frequencies
w = []
fh = open(fname, 'r')
line = fh.readline()
while line:
    if "# k-point" in line:
        line = fh.readline().split()
        while len(line):
            w.append(float(line[0]))
            line = fh.readline().split()
nw = len(w)

print('nkpt =', nkpt)
print('nw =', nw)

# collect data from file
data = np.zeros((nkpt,nw))
fh = open(fname, 'r')
line = fh.readline()
ik = -1
while line:
    if "# k-point" in line:
        ik = ik+1
        line = fh.readline().split()
        iw = -1
        while len(line):
            iw = iw+1
            d = np.array([float(i) for i in line[1:]])
            data[ik,iw] = np.sum(d)
            # print ik, iw, data[ik,iw]
            line = fh.readline().split()
    line = fh.readline()

# output
fh = open("heatMap.dat", "w")
kpath = 0.0
v0 = vkc[0]
for i in range(len(iklist)):
    ik = int(iklist[i])-1
    kpath = kpath + np.sqrt( (vkc[ik][0]-v0[0])**2 + (vkc[ik][1]-v0[1])**2 + (vkc[ik][2]-v0[2])**2 )
    v0 = vkc[ik]
    for iw in range(nw):
        line = "{}  {}  {}\n".format(kpath, w[iw], data[ik,iw])
        fh.write(line)
    fh.write("\n")
