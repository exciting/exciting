#! /usr/bin/python3
#
#  Copyright (C) 2020-2023 GreenX
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#   http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
# implied. See the License for the specific language governing
# permissions and limitations under the License.
#
# @file
# Contains procedures to test Python bindings to IDieL
#
# This module contains the procedures to test Python bindings

#First add IDieL bindings to the path
import sys

#sys.path.append('../IDieL/interfaces/')
sys.path.append('../build/interfaces')

import IDieLPy
import numpy as np

def test_IDieL_python_wrappers():
    # Create an instance of IDieL
    print("[TEST: IDieLPy]")
    try:
        IDieL_object = IDieLPy.IDieLPy()

        # Initialize mesh data
        ngrid = np.array([6, 6, 6], dtype=np.int64)

        # Initialize silicon crystal data
        natoms = np.int64(2)
        latpar = np.float64(0.513)
        lattice = np.array([[0.2565, 0.2565, 0.0],
                            [0.2565, 0.0, 0.2565],
                            [0.0, 0.2565, 0.2565]], dtype=np.float64, order='F')
        redpos = np.array([[0.00, 0.25], [0.00, 0.25], [0.00, 0.25]], dtype=np.float64, order='F')
        types =  np.array([14, 14], dtype=np.int64)

        # Example usage
        print("Initialize")
        IDieL_object.initialize(lattice, natoms, redpos, types, ngrid, 3)
        print("Initialize: Done")

        # Destroy the instance
        IDieL_object.destroy()

        # Return sys exit 0 for CTest
        sys.exit(0)
    except:
        sys.exit(1)

if __name__ == "__main__":
    test_IDieL_python_wrappers()

