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
 This script generates the references to compare against for spherical harmonics tests
"""
import numpy as np
from scipy.special import sph_harm

# Read the mesh file
data = np.loadtxt('ang21.txt')

# Extract phi and theta from the data
theta = data[:, 0]
phi   = data[:, 1]

# Compute spherical harmonics
l_max = 5  # Set the maximum angular momentum quantum number
m_values = np.arange(-l_max, l_max + 1)
spherical_harmonics = np.empty((len(phi), (l_max+1)**2 ), dtype=np.complex128)

i = 0
for l in range(0,l_max+1):
    m_values = np.arange(-l, l + 1)
    for m in m_values:
        print(l,m,i)
        spherical_harmonics[:, i] = sph_harm(m, l, phi, theta)
        i += 1

# Separate real and imaginary parts
real_part = np.real(spherical_harmonics)
imag_part = np.imag(spherical_harmonics)

# Save real and imaginary parts to separate files
np.savetxt('real_part.txt', real_part, fmt='%.16e')
np.savetxt('imag_part.txt', imag_part, fmt='%.16e')
