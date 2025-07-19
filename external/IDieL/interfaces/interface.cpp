// Copyright (C) 2020-2024 GreenX library
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//   http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
// implied. See the License for the specific language governing
// permissions and limitations under the License.

/// @file
/// This file contains a module for exposing the idiel_t to C++, and python
/// Python bindings are created through pybind11 library

#include <memory>
#include <complex>
#include <pybind11/pybind11.h>
#include <pybind11/complex.h>
#include <pybind11/numpy.h>
#include <iostream>

/// Shortcut to the namespace
namespace py = pybind11;

/// Declaring type aliases
using Complex = std::complex<double>;
using NumpyComplexArray = py::array_t<Complex,py::array::f_style>;
using NumpyRealArray = py::array_t<double,py::array::f_style>;
using NumpyIntArray = py::array_t< int,py::array::f_style>;


/// External functions to call Fortran functions
/// The Fortran type is stored as a C pointer for purposes of this file
extern "C" {
    extern void allocate_idiel_t(void** object_ptr);
    extern void deallocate_idiel_t(void** object_ptr);
    extern void init_common(void** object_ptr, double* lattice,  int natoms, double* redpos,  int* elements,  int* nq,  int dim);
    extern void set_dielectric_blocks(void** object_ptr, Complex* h, Complex* wl, Complex* wu, Complex* ib);
    extern void compute_anisotropic_avg_inversedielectric_3d(void** object_ptr, bool hermitian);
    extern void compute_anisotropic_avg_scrcoulomb_2d(void** object_ptr, bool hermitian);
    extern void invert_body(void** object_ptr, Complex* body);
    extern  int get_n_basis(void** object_ptr);
    extern void* head(void** object_ptr);
    extern void* wing_lower(void** object_ptr);
    extern void* wing_upper(void** object_ptr);
    extern void* body(void** object_ptr);
}

/// Now we define a C++ class to manipulate the object
class idiel_cxx {

private:
    // This is a raw pointer containing the Fortran type
    // TODO: Change to smart pointer (unique_ptr)
    void* idiel_f90;

public:

    // Constructor: we need this call as Fortran to be manipulated pointers
    // must be initialized by Fortran
    idiel_cxx() {
        allocate_idiel_t(&idiel_f90);
    }

    // Destructor: Fortran clean-up
    ~idiel_cxx() {
        deallocate_idiel_t(&idiel_f90);
    }

    // Declare copy and movement constructor deleted
    // This is necessary because of the raw pointer
    idiel_cxx(const idiel_cxx& A) = delete;
    idiel_cxx(idiel_cxx&& A)      = delete;

    /// Only in case, this is called automatically in C++ when out-of-scope
    /// So we call Fortran deallocation routines
    void destroy(){
        deallocate_idiel_t(&idiel_f90);
    }

    /// Initialize common elements, note that I make dim mandatory as it makes the interfacing easier
    void initialize(double* lattice,  int natoms, double* redpos,  int* elements,  int* nq,  int dim) {
        init_common(&idiel_f90, lattice, natoms, redpos, elements, nq, dim);
    }

    /// We support both functions, C++ treat optional as two different functions
    void setDielectricBlocksFull(Complex* h, Complex* wl, Complex* wu, Complex* ib) {
        set_dielectric_blocks(&idiel_f90, h, wl, wu, ib);
    }

    /// We support both functions, C++ treat optional as two different functions
    void setDielectricBlocksPartial(Complex* h, Complex* wl, Complex* wu) {
        set_dielectric_blocks(&idiel_f90, h, wl, wu, nullptr);
    }

    /// Computing the average (3D)
    void ComputeAniAvgInvDielMat3D(bool hermitian) {
        compute_anisotropic_avg_inversedielectric_3d(&idiel_f90, hermitian);
    }

    /// Compute the average (2D)
    void ComputeAniAvgScrCoulomb2D(bool hermitian) {
        compute_anisotropic_avg_scrcoulomb_2d(&idiel_f90, hermitian);
    }

    /// Invert body and save data to Fortran object
    void invertBody(Complex* body) {
        invert_body(&idiel_f90, body);
    }

    /// Return the basis size
     int getNBasis(){
        return get_n_basis(&idiel_f90);
    }

    /// Return C-pointer head
    void* get_head(){
        return head(&idiel_f90);
    }

    /// Return C-pointer wingL
    void* get_wingL(){
        return wing_lower(&idiel_f90);
    }

    /// Return C-pointer wingU
    void* get_wingU(){
        return wing_upper(&idiel_f90);
    }

    /// Return C-pointer body
    void* get_body(){
        return body(&idiel_f90);
    }
};

/// This is how a C++ class translates to Python module
/// TODO: Add documentation to the python module
PYBIND11_MODULE(IDieLPy, m) {

    m.doc() = "Python bindings for IDieL";

    py::class_<idiel_cxx>(m, "IDieLPy")
        .def(py::init<>())
        .def("initialize", [](idiel_cxx &self, NumpyRealArray& lattice, int natoms, NumpyRealArray& redpos, NumpyIntArray& elements, NumpyIntArray& nq, int dim){
            auto lattice_ptr = lattice.mutable_unchecked<2>();
            auto redpos_ptr = redpos.mutable_unchecked<2>();
            auto elements_ptr = elements.mutable_unchecked<1>();
            auto nq_ptr = nq.mutable_unchecked<1>();

            // Get rid of constness
            double* lattice_data = const_cast<double*>(lattice_ptr.data(0, 0));
            double* redpos_data = const_cast<double*>(redpos_ptr.data(0, 0));
             int* elements_data = const_cast< int*>(elements_ptr.data(0));
             int* nq_data = const_cast< int*>(nq_ptr.data(0));

            self.initialize(lattice_data, static_cast< int>(natoms), redpos_data, elements_data, nq_data, static_cast< int>(dim));
        })
        .def("setDielectricBlocksFull", [](idiel_cxx &self, NumpyComplexArray& head, NumpyComplexArray& wingL, NumpyComplexArray& wingU, NumpyComplexArray& Binv){

            auto head_ptr  =  head.mutable_unchecked<2>();
            auto wingL_ptr =  wingL.mutable_unchecked<2>();
            auto wingU_ptr =  wingU.mutable_unchecked<2>();
            auto Binv_ptr  =  Binv.mutable_unchecked<2>();

            // Get rid of constness
            Complex* head_data  = const_cast<Complex*>(head_ptr.data(0,0));
            Complex* wingL_data = const_cast<Complex*>(wingL_ptr.data(0,0));
            Complex* wingU_data = const_cast<Complex*>(wingU_ptr.data(0,0));
            Complex* Binv_data  = const_cast<Complex*>(Binv_ptr.data(0,0));

            self.setDielectricBlocksFull(head_data, wingL_data, wingU_data, Binv_data);
        })
        .def("setDielectricBlocksPartial", [](idiel_cxx &self, NumpyComplexArray& head, NumpyComplexArray& wingL, NumpyComplexArray& wingU){

            auto head_ptr  =  head.mutable_unchecked<2>();
            auto wingL_ptr =  wingL.mutable_unchecked<2>();
            auto wingU_ptr =  wingU.mutable_unchecked<2>();

            // Get rid of constness
            Complex* head_data  = const_cast<Complex*>(head_ptr.data(0,0));
            Complex* wingL_data = const_cast<Complex*>(wingL_ptr.data(0,0));
            Complex* wingU_data = const_cast<Complex*>(wingU_ptr.data(0,0));

            self.setDielectricBlocksPartial(head_data, wingL_data, wingU_data);
        })
        .def("ComputeAniAvgInvDielMat3D",  &idiel_cxx::ComputeAniAvgInvDielMat3D)
        .def("ComputeAniAvgScrCoulomb2D",  &idiel_cxx::ComputeAniAvgScrCoulomb2D)
        .def("invertBody", [](idiel_cxx &self, NumpyComplexArray& body){
            auto body_ptr  =  body.mutable_unchecked<2>();
            Complex* body_data  = const_cast<Complex*>(body_ptr.data(0,0));
            self.invertBody(body_data);
        })
        .def("getNBasis", &idiel_cxx::getNBasis)
        .def("get_head", [](idiel_cxx &self)->Complex {
            auto ptr = self.get_head();
            return *(reinterpret_cast<Complex*>(ptr));
        })
        .def("get_wingL", [](idiel_cxx &self)->NumpyComplexArray {
            auto n   = self.getNBasis();
            auto ptr = self.get_wingL();
            return NumpyComplexArray(n, reinterpret_cast<Complex*>(ptr));
        })
        .def("get_wingU", [](idiel_cxx &self)->NumpyComplexArray {
            auto n   = self.getNBasis();
            auto ptr = self.get_wingU();
            return NumpyComplexArray(n, reinterpret_cast<Complex*>(ptr));
        })
        .def("get_body", [](idiel_cxx &self)->NumpyComplexArray {
            auto n   = self.getNBasis();
            auto ptr = self.get_body();
            NumpyComplexArray body(n, reinterpret_cast<Complex*>(ptr));
            body.reshape({n,n});
            return body;
        })
        .def("destroy", &idiel_cxx::destroy);
};
