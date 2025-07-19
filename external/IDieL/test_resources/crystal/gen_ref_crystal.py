"""
 Copyright (C) 2020-2024 GreenX library
 
 Licensed under the Apache License, Version 2.0 (the "License");
 you may not use this file except in compliance with the License.
 You may obtain a copy of the License at
 
   http://www.apache.org/licenses/LICENSE-2.0
 
 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
 implied. See the License for the specific language governing
 permissions and limitations under the License.
 
  @file
  This script generates the references to compare against for crystal structure and symmetry tests
"""
import numpy as np
import scipy.linalg as la
import spglib 

#Build lattice
a = 3.9224652813161520 * np.array([1.0,  0.0,               0.0])
b = 3.9224652813161520 * np.array([-0.5, np.sqrt(3.0)/2.0, 0.0])
c = 3.9224652813161520 * np.array([0.0,  0.0,               1.6452947797075332])

lattice = np.zeros(shape=(3,3),dtype=np.double)
lattice[0,:] = a
lattice[1,:] = b
lattice[2,:] = c

#Number of atoms
natoms = 4

# Atomic types
types = np.array([30,30,34,34])

# Reduced positions
redpos = np.zeros(shape=(natoms,3),dtype=np.double)
redpos[0,:] = np.array([0.3333333333333357, 0.6666666666666643, 0.0000000000000000])
redpos[1,:] = np.array([0.6666666666666643, 0.3333333333333357, 0.5000000000000000])
redpos[2,:] = np.array([0.3333333333333357, 0.6666666666666643, 0.3736232276372607])
redpos[3,:] = np.array([0.6666666666666643, 0.3333333333333357, 0.8736232276372607])

# Compute vuc
vuc = la.det(lattice)
print("vuc = {:16.12f} nm**3".format(vuc))

#Compute the reciprocal lattice
ilat = 2.0 * np.pi * la.inv(lattice).T
print('a* : {:16.12f} {:16.12f} {:16.12f}'.format(ilat[0,0],ilat[0,1],ilat[0,2]))
print('b* : {:16.12f} {:16.12f} {:16.12f}'.format(ilat[1,0],ilat[1,1],ilat[1,2]))
print('c* : {:16.12f} {:16.12f} {:16.12f}'.format(ilat[2,0],ilat[2,1],ilat[2,2]))

# Now compute the crystal symmetry
cell = (lattice, redpos, types)
sym  = spglib.get_symmetry_dataset(cell)
print('Space group : ', sym['number'])

# Generate symmetrizers
Rs = sym['rotations']
print(np.shape(Rs))
i = 1
for r in Rs:
    print('rot_ref(1,:,{}) = (/{}_i64, {}_i64, {}_i64/)'.format(i,r[0,0],r[0,1],r[0,2]))
    print('rot_ref(2,:,{}) = (/{}_i64, {}_i64, {}_i64/)'.format(i,r[1,0],r[1,1],r[1,2]))
    print('rot_ref(3,:,{}) = (/{}_i64, {}_i64, {}_i64/)'.format(i,r[2,0],r[2,1],r[2,2]))
    i += 1
print()
#Computing cartesian rotation
i = 1
lat = lattice.T
for r0 in Rs:
    t1 = lat.T
    t2 = (lat @ r0).T
    r = (la.inv(t1) @ t2).T
    print('crot_ref(1,:,{}) = (/{}_r64, {}_r64, {}_r64/)'.format(i,r[0,0],r[0,1],r[0,2]))
    print('crot_ref(2,:,{}) = (/{}_r64, {}_r64, {}_r64/)'.format(i,r[1,0],r[1,1],r[1,2]))
    print('crot_ref(3,:,{}) = (/{}_r64, {}_r64, {}_r64/)'.format(i,r[2,0],r[2,1],r[2,2]))
    i += 1


