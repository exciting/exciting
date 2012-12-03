#!/usr/bin/env python

from ase.visualize.vtk import requirevtk, probe_vtk_kilobyte

requirevtk()
vtk_kilobyte = probe_vtk_kilobyte(1024)

import numpy as np

import sys, unittest, gc
from ase.test import CustomTestCase, CustomTextTestRunner
from ase.utils.memory import MemoryStatistics, MemorySingleton, shapeopt

from vtk import vtkDataArray
from ase.visualize.vtk.data import vtkFloatArrayFromNumPyArray, \
                                   vtkDoubleArrayFromNumPyArray, \
                                   vtkFloatArrayFromNumPyMultiArray, \
                                   vtkDoubleArrayFromNumPyMultiArray

# -------------------------------------------------------------------

class UTConversionDataArrayNumPy(CustomTestCase):
    """
    Abstract class with test cases for VTK/NumPy data conversion.

    Leak tests the six possible permutations of deletion order for the
    objects involved in conversion between VTK and NumPy data formats.

    Objects:

    conv: instance of vtkDataArrayFromNumPyBuffer of subclass thereof
        The object in charge of the conversion
    data: NumPy array
        NumPy data with a specific memory footprint
    vtk_da: instance of vtkDataArray of subclass thereof
        VTK data array with a similar memory footprint

    Permutations:

    Case A: 012 i.e. deletion order is conv, data, vtk_da
    Case B: 021 i.e. deletion order is conv, vtk_da, data
    Case C: 102 i.e. deletion order is data, conv, vtk_da
    Case D: 120 i.e. deletion order is data, vtk_da, conv
    Case E: 201 i.e. deletion order is vtk_da, conv, data
    Case F: 210 i.e. deletion order is vtk_da, data, conv
    """
    footprint = 100*1024**2
    dtype = None
    verbose = 0
    gc_threshold = (300,5,5) #default is (700,10,10)
    gc_flags = gc.DEBUG_LEAK # | gc.DEBUG_STATS
    ctol = -7 #10MB
    etol = -7 #10MB

    def setUp(self):
        self.mem_ini = MemorySingleton(self.verbose-1)
        self.mem_ref = MemoryStatistics(self.verbose-1)
        self.mem_cur = self.mem_ref.copy()

        self.gc_threshold_old = gc.get_threshold()
        self.gc_flags_old = gc.get_debug()
        gc.set_threshold(*self.gc_threshold)
        gc.set_debug(self.gc_flags)

        # Try to obtain a clean slate
        gc.collect()
        self.gc_count = len(gc.garbage)
        del gc.garbage[:]

    def tearDown(self):
        gc.collect()
        self.assertEqual(len(gc.garbage), self.gc_count)
        if len(gc.garbage)>0:
            if self.verbose>1: print gc.get_objects()
            #TODO be pedantic and fail?
        del gc.garbage[:]
        gc.set_threshold(*self.gc_threshold_old)
        gc.set_debug(self.gc_flags_old)

    def assertAlmostConsumed(self, bytes, digits=0, key='VmSize'):
        self.mem_cur.update()
        dm = self.mem_cur-self.mem_ref
        self.assertAlmostEqual(dm[key], bytes, digits)

    def assertAlmostExceeded(self, bytes, digits=0, key='VmPeak'):
        self.mem_cur.update()
        dm = self.mem_cur-self.mem_ini
        #self.assertAlmostEqual(dm[key], bytes, digits) #TODO what really?
        #self.assertAlmostEqual(max(0, dm[key]-bytes), 0, digits) #TODO ???
        #dm = 200 MB, bytes = 100MB     ok
        #dm = 101 MB, bytes = 100MB     ok
        #dm = 0 MB, bytes = 100MB       bad

    def convert_to_vtk_array(self, data):
        """Convert an ndarray to a VTK data array.

         data: NumPy array
            NumPy data with a specific memory footprint
        """

        raise RuntimeError('Virtual member function.')

    def get_leaktest_scenario(self):
        """Construct the necessary conversion objects for leak testing.

        Returns tuple of the form (conv, data, vtk_da,) where:

        conv: instance of vtkDataArrayFromNumPyBuffer of subclass thereof
            The object in charge of the conversion
        data: NumPy array
            NumPy data with a specific memory footprint
        vtk_da: instance of vtkDataArray of subclass thereof
            VTK data array with a similar memory footprint
        """

        raise RuntimeError('Virtual member function.')

    # =================================

    def test_deletion_case_a(self):
        # Case A: 012 i.e. deletion order is conv, data, vtk_da
        (conv, data, vtk_da,) = self.get_leaktest_scenario()

        # Conversion cleanup
        del conv
        self.assertAlmostConsumed(2*self.footprint, self.ctol)
        self.assertAlmostExceeded(2*self.footprint, self.etol)
        if self.verbose>=1: print 'Conversion cleanup=', self.mem_cur-self.mem_ref

        # NumPy cleanup
        del data
        self.assertAlmostConsumed(self.footprint, self.ctol)
        self.assertAlmostExceeded(2*self.footprint, self.etol)
        if self.verbose>=1: print 'NumPy cleanup=', self.mem_cur-self.mem_ref

        # VTK cleanup
        del vtk_da
        self.assertAlmostConsumed(0, self.ctol)
        self.assertAlmostExceeded(2*self.footprint, self.etol)
        if self.verbose>=1: print 'VTK cleanup=', self.mem_cur-self.mem_ref

    def test_deletion_case_b(self):
        # Case B: 021 i.e. deletion order is conv, vtk_da, data
        (conv, data, vtk_da,) = self.get_leaktest_scenario()

        # Conversion cleanup
        del conv
        self.assertAlmostConsumed(2*self.footprint, self.ctol)
        self.assertAlmostExceeded(2*self.footprint, self.etol)
        if self.verbose>=1: print 'Conversion cleanup=', self.mem_cur-self.mem_ref

        # VTK cleanup
        del vtk_da
        self.assertAlmostConsumed(self.footprint, self.ctol)
        self.assertAlmostExceeded(2*self.footprint, self.etol)
        if self.verbose>=1: print 'VTK cleanup=', self.mem_cur-self.mem_ref

        # Numpy cleanup
        del data
        self.assertAlmostConsumed(0, self.ctol)
        self.assertAlmostExceeded(2*self.footprint, self.etol)
        if self.verbose>=1: print 'NumPy cleanup=', self.mem_cur-self.mem_ref

    def test_deletion_case_c(self):
        # Case C: 102 i.e. deletion order is data, conv, vtk_da
        (conv, data, vtk_da,) = self.get_leaktest_scenario()

        # NumPy cleanup
        del data
        self.assertAlmostConsumed(self.footprint, self.ctol)
        self.assertAlmostExceeded(2*self.footprint, self.etol)
        if self.verbose>=1: print 'NumPy cleanup=', self.mem_cur-self.mem_ref

        # Conversion cleanup
        del conv
        self.assertAlmostConsumed(self.footprint, self.ctol)
        self.assertAlmostExceeded(2*self.footprint, self.etol)
        if self.verbose>=1: print 'Conversion cleanup=', self.mem_cur-self.mem_ref

        # VTK cleanup
        del vtk_da
        self.assertAlmostConsumed(0, self.ctol)
        self.assertAlmostExceeded(2*self.footprint, self.etol)
        if self.verbose>=1: print 'VTK cleanup=', self.mem_cur-self.mem_ref

    def test_deletion_case_d(self):
        # Case D: 120 i.e. deletion order is data, vtk_da, conv
        (conv, data, vtk_da,) = self.get_leaktest_scenario()

        # NumPy cleanup
        del data
        self.assertAlmostConsumed(self.footprint, self.ctol)
        self.assertAlmostExceeded(2*self.footprint, self.etol)
        if self.verbose>=1: print 'NumPy cleanup=', self.mem_cur-self.mem_ref

        # VTK cleanup
        del vtk_da
        self.assertAlmostConsumed(self.footprint, self.ctol)
        self.assertAlmostExceeded(2*self.footprint, self.etol)
        if self.verbose>=1: print 'VTK cleanup=', self.mem_cur-self.mem_ref

        # Conversion cleanup
        del conv
        self.assertAlmostConsumed(0, self.ctol)
        self.assertAlmostExceeded(2*self.footprint, self.etol)
        if self.verbose>=1: print 'Conversion cleanup=', self.mem_cur-self.mem_ref

    def test_deletion_case_e(self):
        # Case E: 201 i.e. deletion order is vtk_da, conv, data
        (conv, data, vtk_da,) = self.get_leaktest_scenario()

        # VTK cleanup
        del vtk_da
        self.assertAlmostConsumed(2*self.footprint, self.ctol)
        self.assertAlmostExceeded(2*self.footprint, self.etol)
        if self.verbose>=1: print 'VTK cleanup=', self.mem_cur-self.mem_ref

        # Conversion cleanup
        del conv
        self.assertAlmostConsumed(self.footprint, self.ctol)
        self.assertAlmostExceeded(2*self.footprint, self.etol)
        if self.verbose>=1: print 'Conversion cleanup=', self.mem_cur-self.mem_ref

        # NumPy cleanup
        del data
        self.assertAlmostConsumed(0, self.ctol)
        self.assertAlmostExceeded(2*self.footprint, self.etol)
        if self.verbose>=1: print 'NumPy cleanup=', self.mem_cur-self.mem_ref

    def test_deletion_case_f(self):
        # Case F: 210 i.e. deletion order is vtk_da, data, conv
        (conv, data, vtk_da,) = self.get_leaktest_scenario()

        # VTK cleanup
        del vtk_da
        self.assertAlmostConsumed(2*self.footprint, self.ctol)
        self.assertAlmostExceeded(2*self.footprint, self.etol)
        if self.verbose>=1: print 'VTK cleanup=', self.mem_cur-self.mem_ref

        # NumPy cleanup
        del data
        self.assertAlmostConsumed(self.footprint, self.ctol)
        self.assertAlmostExceeded(2*self.footprint, self.etol)
        if self.verbose>=1: print 'NumPy cleanup=', self.mem_cur-self.mem_ref

        # Conversion cleanup
        del conv
        self.assertAlmostConsumed(0, self.ctol)
        self.assertAlmostExceeded(2*self.footprint, self.etol)
        if self.verbose>=1: print 'Conversion cleanup=', self.mem_cur-self.mem_ref

# -------------------------------------------------------------------

# class UTDataArrayFromNumPyBuffer(...): TODO

# -------------------------------------------------------------------

class UTDataArrayFromNumPyArray_Scalar(UTConversionDataArrayNumPy):
    """
    Test cases for memory leaks during VTK/NumPy data conversion.
    Conversion of 1D NumPy array to VTK data array using buffers."""

    def setUp(self):
        UTConversionDataArrayNumPy.setUp(self)
        self.shape = self.footprint//np.nbytes[self.dtype]

    def get_leaktest_scenario(self):
        self.assertAlmostEqual(np.prod(self.shape)*np.nbytes[self.dtype], \
                               self.footprint, -2) #100B

        # Update memory reference
        self.mem_ref.update()

        # NumPy allocation
        data = np.empty(self.shape, self.dtype)
        self.assertAlmostConsumed(self.footprint, self.ctol)
        self.assertAlmostExceeded(self.footprint, self.etol)
        if self.verbose>=1: print 'NumPy allocation=', self.mem_cur-self.mem_ref

        # NumPy to VTK conversion
        np2da = self.convert_to_vtk_array(data)
        self.assertAlmostConsumed(2*self.footprint, self.ctol)
        self.assertAlmostExceeded(2*self.footprint, self.etol)
        if self.verbose>=1: print 'Conversion=', self.mem_cur-self.mem_ref

        # VTK retrieval
        vtk_da = np2da.get_output()
        self.assertTrue(isinstance(vtk_da, vtkDataArray))
        self.assertAlmostEqual(vtk_da.GetActualMemorySize()*vtk_kilobyte, \
                               self.footprint, -3) #1kB
        if self.verbose>=1: print 'VTK retrieval=', self.mem_cur-self.mem_ref

        return (np2da, data, vtk_da,)

class UTFloatArrayFromNumPyArray_Scalar(UTDataArrayFromNumPyArray_Scalar):
    __doc__ = UTDataArrayFromNumPyArray_Scalar.__doc__
    dtype = np.float32
    convert_to_vtk_array = vtkFloatArrayFromNumPyArray

class UTDoubleArrayFromNumPyArray_Scalar(UTDataArrayFromNumPyArray_Scalar):
    __doc__ = UTDataArrayFromNumPyArray_Scalar.__doc__
    dtype = np.float64
    convert_to_vtk_array = vtkDoubleArrayFromNumPyArray

class UTDataArrayFromNumPyArray_Vector(UTConversionDataArrayNumPy):
    """
    Test cases for memory leaks during VTK/NumPy data conversion.
    Conversion of 2D NumPy array to VTK data array using buffers."""

    def setUp(self):
        UTConversionDataArrayNumPy.setUp(self)
        size = self.footprint//np.nbytes[self.dtype]
        self.shape = (size//3, 3)

    def get_leaktest_scenario(self):
        self.assertAlmostEqual(np.prod(self.shape)*np.nbytes[self.dtype], \
                               self.footprint, -2) #100B

        # Update memory reference
        self.mem_ref.update()

        # NumPy allocation
        data = np.empty(self.shape, self.dtype)
        self.assertAlmostConsumed(self.footprint, self.ctol)
        self.assertAlmostExceeded(self.footprint, self.etol)
        if self.verbose>=1: print 'NumPy allocation=', self.mem_cur-self.mem_ref

        # NumPy to VTK conversion
        np2da = self.convert_to_vtk_array(data)
        self.assertAlmostConsumed(2*self.footprint, self.ctol)
        self.assertAlmostExceeded(2*self.footprint, self.etol)
        if self.verbose>=1: print 'Conversion=', self.mem_cur-self.mem_ref

        # VTK retrieval
        vtk_da = np2da.get_output()
        self.assertTrue(isinstance(vtk_da, vtkDataArray))
        self.assertAlmostEqual(vtk_da.GetActualMemorySize()*vtk_kilobyte, \
                               self.footprint, -3) #1kB
        if self.verbose>=1: print 'VTK retrieval=', self.mem_cur-self.mem_ref

        return (np2da, data, vtk_da,)

class UTFloatArrayFromNumPyArray_Vector(UTDataArrayFromNumPyArray_Vector):
    __doc__ = UTDataArrayFromNumPyArray_Vector.__doc__
    dtype = np.float32
    convert_to_vtk_array = vtkFloatArrayFromNumPyArray

class UTDoubleArrayFromNumPyArray_Vector(UTDataArrayFromNumPyArray_Vector):
    __doc__ = UTDataArrayFromNumPyArray_Vector.__doc__
    dtype = np.float64
    convert_to_vtk_array = vtkDoubleArrayFromNumPyArray

# -------------------------------------------------------------------

class UTDataArrayFromNumPyMultiArray_Scalar(UTConversionDataArrayNumPy):
    """
    Test cases for memory leaks during VTK/NumPy data conversion.
    Conversion of NumPy grid scalars to VTK data array using buffers."""

    def setUp(self):
        UTConversionDataArrayNumPy.setUp(self)
        size = self.footprint//np.nbytes[self.dtype]
        digits, shape = shapeopt(1000, size, ndims=3, ecc=0.3)
        if self.verbose>=1: print 'digits=%8.3f, shape=%s' % (digits,shape)
        self.shape = shape + (1,)
        self.assertAlmostEqual(np.prod(self.shape)*np.nbytes[self.dtype], \
                               self.footprint, -4) #10kB

    def get_leaktest_scenario(self):

        # Update memory reference
        self.mem_ref.update()

        # NumPy allocation
        data = np.empty(self.shape, self.dtype)
        self.assertAlmostConsumed(self.footprint, self.ctol)
        self.assertAlmostExceeded(self.footprint, self.etol)
        if self.verbose>=1: print 'NumPy allocation=', self.mem_cur-self.mem_ref

        # NumPy to VTK conversion
        np2da = self.convert_to_vtk_array(data)
        self.assertAlmostConsumed(2*self.footprint, self.ctol)
        self.assertAlmostExceeded(2*self.footprint, self.etol)
        if self.verbose>=1: print 'Conversion=', self.mem_cur-self.mem_ref

        # VTK retrieval
        vtk_da = np2da.get_output()
        self.assertTrue(isinstance(vtk_da, vtkDataArray))
        self.assertAlmostEqual(vtk_da.GetActualMemorySize()*vtk_kilobyte, \
                               self.footprint, -4) #10kB
        if self.verbose>=1: print 'VTK retrieval=', self.mem_cur-self.mem_ref

        return (np2da, data, vtk_da,)

class UTFloatArrayFromNumPyMultiArray_Scalar(UTDataArrayFromNumPyMultiArray_Scalar):
    __doc__ = UTDataArrayFromNumPyMultiArray_Scalar.__doc__
    dtype = np.float32
    convert_to_vtk_array = vtkFloatArrayFromNumPyMultiArray

class UTDoubleArrayFromNumPyMultiArray_Scalar(UTDataArrayFromNumPyMultiArray_Scalar):
    __doc__ = UTDataArrayFromNumPyMultiArray_Scalar.__doc__
    dtype = np.float64
    convert_to_vtk_array = vtkDoubleArrayFromNumPyMultiArray

class UTDataArrayFromNumPyMultiArray_Vector(UTConversionDataArrayNumPy):
    """
    Test cases for memory leaks during VTK/NumPy data conversion.
    Conversion of NumPy grid vectors to VTK data array using buffers."""

    def setUp(self):
        UTConversionDataArrayNumPy.setUp(self)
        size = self.footprint//np.nbytes[self.dtype]
        digits, shape = shapeopt(1000, size//3, ndims=3, ecc=0.3)
        if self.verbose>=1: print 'digits=%8.3f, shape=%s' % (digits,shape)
        self.shape = shape + (3,)
        self.assertAlmostEqual(np.prod(self.shape)*np.nbytes[self.dtype], \
                               self.footprint, -4) #10kB

    def get_leaktest_scenario(self):

        # Update memory reference
        self.mem_ref.update()

        # NumPy allocation
        data = np.empty(self.shape, self.dtype)
        self.assertAlmostConsumed(self.footprint, self.ctol)
        self.assertAlmostExceeded(self.footprint, self.etol)
        if self.verbose>=1: print 'NumPy allocation=', self.mem_cur-self.mem_ref

        # NumPy to VTK conversion
        np2da = self.convert_to_vtk_array(data)
        self.assertAlmostConsumed(2*self.footprint, self.ctol)
        self.assertAlmostExceeded(2*self.footprint, self.etol)
        if self.verbose>=1: print 'Conversion=', self.mem_cur-self.mem_ref

        # VTK retrieval
        vtk_da = np2da.get_output()
        self.assertTrue(isinstance(vtk_da, vtkDataArray))
        self.assertAlmostEqual(vtk_da.GetActualMemorySize()*vtk_kilobyte, \
                               self.footprint, -4) #10kB
        if self.verbose>=1: print 'VTK retrieval=', self.mem_cur-self.mem_ref

        return (np2da, data, vtk_da,)

class UTFloatArrayFromNumPyMultiArray_Vector(UTDataArrayFromNumPyMultiArray_Vector):
    __doc__ = UTDataArrayFromNumPyMultiArray_Vector.__doc__
    dtype = np.float32
    convert_to_vtk_array = vtkFloatArrayFromNumPyMultiArray

class UTDoubleArrayFromNumPyMultiArray_Vector(UTDataArrayFromNumPyMultiArray_Vector):
    __doc__ = UTDataArrayFromNumPyMultiArray_Vector.__doc__
    dtype = np.float64
    convert_to_vtk_array = vtkDoubleArrayFromNumPyMultiArray

# -------------------------------------------------------------------

if __name__ in ['__main__', '__builtin__']:
    # We may have been imported by test.py, if so we should redirect to logfile
    if __name__ == '__builtin__':
        testrunner = CustomTextTestRunner('vtk_data.log', verbosity=2)
    else:
        testrunner = unittest.TextTestRunner(stream=sys.stdout, verbosity=2)

    testcases = [UTFloatArrayFromNumPyArray_Scalar, \
                 UTDoubleArrayFromNumPyArray_Scalar, \
                 UTFloatArrayFromNumPyArray_Vector, \
                 UTDoubleArrayFromNumPyArray_Vector, \
                 UTFloatArrayFromNumPyMultiArray_Scalar, \
                 UTDoubleArrayFromNumPyMultiArray_Scalar, \
                 UTFloatArrayFromNumPyMultiArray_Vector, \
                 UTDoubleArrayFromNumPyMultiArray_Vector]

    for test in testcases:
        info = '\n' + test.__name__ + '\n' + test.__doc__.strip('\n') + '\n'
        testsuite = unittest.defaultTestLoader.loadTestsFromTestCase(test)
        testrunner.stream.writeln(info)
        testresult = testrunner.run(testsuite)
        # Provide feedback on failed tests if imported by test.py
        if __name__ == '__builtin__' and not testresult.wasSuccessful():
            raise SystemExit('Test failed. Check vtk_data.log for details.')

