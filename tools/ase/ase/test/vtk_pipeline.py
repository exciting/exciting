#!/usr/bin/env python

from ase.visualize.vtk import requirevtk

requirevtk()

import sys, unittest
from ase.test import CustomTestCase, CustomTextTestRunner

from vtk import vtkContourFilter, vtkPolyDataNormals, \
                vtkLinearSubdivisionFilter, vtkPolyDataMapper
from ase.visualize.vtk.pipeline import vtkPolyDataPipeline

# -------------------------------------------------------------------

class UTPipeline(CustomTestCase):
    """
    Abstract test case class - TODO."""

    def assertConnected(self, vtk_one, vtk_two, port=0):
        self.assertEqual(vtk_two.GetNumberOfInputConnections(port), 1)
        self.assertEqual(vtk_one.GetOutputPort(),
                         vtk_two.GetInputConnection(port, 0))

    def assertNotConnected(self, vtk_one, vtk_two, port=0):
        self.assertEqual(vtk_two.GetNumberOfInputConnections(port), 0)


# -------------------------------------------------------------------

class UTPolyDataPipeline(UTPipeline):
    """
    TODO"""

    def setUp(self):
        self.vtk_iso = vtkContourFilter()
        #self.vtk_iso.SetInput(...)

        self.vtk_dnorm = vtkPolyDataNormals()
        self.vtk_subdiv = vtkLinearSubdivisionFilter()
        self.vtk_dmap = vtkPolyDataMapper()

    def tearDown(self):
        del self.vtk_dmap
        del self.vtk_subdiv
        del self.vtk_dnorm
        del self.vtk_iso

# -------------------------------------------------------------------

class UTPolyDataPipeline_PureVTK(UTPolyDataPipeline):
    """
    Consistency tests for pure-VTK objects with the purpose
    of ultimately mapping the polygonal data within the pipeline."""

    # =================================

    def test_consistent_data_oninit(self):
        pipeline = vtkPolyDataPipeline(self.vtk_iso)
        self.assertTrue(pipeline.hasdata())
        self.assertEqual(pipeline.vtkish_data, self.vtk_iso)
        self.assertFalse(pipeline.hasfilters())
        self.assertFalse(pipeline.isclosed())

        pipeline.append(self.vtk_dnorm)
        self.assertConnected(self.vtk_iso, self.vtk_dnorm)

        pipeline.append(self.vtk_subdiv)
        self.assertConnected(self.vtk_dnorm, self.vtk_subdiv)

        pipeline.append(self.vtk_dmap)
        self.assertConnected(self.vtk_subdiv, self.vtk_dmap)

        self.assertTrue(pipeline.hasfilters())
        self.assertTrue(pipeline.isclosed())

    def test_consistent_data_direct(self):
        pipeline = vtkPolyDataPipeline()
        self.assertFalse(pipeline.hasdata())
        self.assertFalse(pipeline.hasfilters())
        self.assertFalse(pipeline.isclosed())

        pipeline.set_data(self.vtk_iso)
        self.assertTrue(pipeline.hasdata())
        self.assertEqual(pipeline.vtkish_data, self.vtk_iso)

        pipeline.append(self.vtk_dnorm)
        self.assertConnected(self.vtk_iso, self.vtk_dnorm)

        pipeline.append(self.vtk_subdiv)
        self.assertConnected(self.vtk_dnorm, self.vtk_subdiv)

        pipeline.append(self.vtk_dmap)
        self.assertConnected(self.vtk_subdiv, self.vtk_dmap)

        self.assertTrue(pipeline.hasfilters())
        self.assertTrue(pipeline.isclosed())

    def test_consistent_data_postponed(self):
        pipeline = vtkPolyDataPipeline()
        self.assertFalse(pipeline.hasdata())
        self.assertFalse(pipeline.hasfilters())
        self.assertFalse(pipeline.isclosed())

        pipeline.append(self.vtk_dnorm)
        self.assertTrue(pipeline.hasfilters())
        self.assertNotConnected(self.vtk_iso, self.vtk_dnorm)

        pipeline.append(self.vtk_subdiv)
        self.assertConnected(self.vtk_dnorm, self.vtk_subdiv)

        pipeline.set_data(self.vtk_iso)
        self.assertTrue(pipeline.hasdata())
        self.assertEqual(pipeline.vtkish_data, self.vtk_iso)
        self.assertConnected(self.vtk_iso, self.vtk_dnorm)

        pipeline.append(self.vtk_dmap)
        self.assertConnected(self.vtk_subdiv, self.vtk_dmap)
        self.assertTrue(pipeline.hasdata())
        self.assertTrue(pipeline.isclosed())

    def test_consistent_data_onclose(self):
        pipeline = vtkPolyDataPipeline()
        self.assertFalse(pipeline.hasdata())
        self.assertFalse(pipeline.hasfilters())
        self.assertFalse(pipeline.isclosed())

        pipeline.append(self.vtk_dnorm)
        self.assertTrue(pipeline.hasfilters())
        self.assertNotConnected(self.vtk_iso, self.vtk_dnorm)

        pipeline.append(self.vtk_subdiv)
        self.assertConnected(self.vtk_dnorm, self.vtk_subdiv)

        pipeline.append(self.vtk_dmap)
        self.assertConnected(self.vtk_subdiv, self.vtk_dmap)
        self.assertTrue(pipeline.isclosed())

        pipeline.set_data(self.vtk_iso)
        self.assertTrue(pipeline.hasdata())
        self.assertEqual(pipeline.vtkish_data, self.vtk_iso)
        self.assertConnected(self.vtk_iso, self.vtk_dnorm)

    # =================================

    def test_failure_duplicate_data(self):
        pipeline = vtkPolyDataPipeline(self.vtk_iso)
        self.assertRaises(ValueError, pipeline.append, self.vtk_iso)

    def test_failure_duplicate_mixed(self):
        pipeline = vtkPolyDataPipeline(self.vtk_iso)
        pipeline.append(self.vtk_dnorm)
        self.assertRaises(ValueError, pipeline.append, self.vtk_iso)

    def test_failure_duplicate_cyclic(self):
        pipeline = vtkPolyDataPipeline(self.vtk_iso)
        pipeline.append(self.vtk_dnorm)
        pipeline.append(self.vtk_subdiv)
        self.assertRaises(ValueError, pipeline.append, self.vtk_dnorm)

    def test_failure_duplicate_filter(self):
        pipeline = vtkPolyDataPipeline(self.vtk_iso)
        pipeline.append(self.vtk_dnorm)
        self.assertRaises(ValueError, pipeline.append, self.vtk_dnorm)

    # =================================

    def test_failure_output_missing(self):
        pipeline = vtkPolyDataPipeline()
        self.assertRaises(RuntimeError, pipeline.get_output_port)

    def test_failure_output_closed(self):
        pipeline = vtkPolyDataPipeline(self.vtk_iso)
        pipeline.append(self.vtk_dmap)
        self.assertRaises(RuntimeError, pipeline.append, self.vtk_dnorm)

# -------------------------------------------------------------------

class UTPolyDataPipeline_PipelineVTK(UTPolyDataPipeline):
    """
    Consistency tests for mixed VTK objects and pipelines with the purpose
    of ultimately mapping the polygonal data through embedded pipelines."""

    # =================================

    def get_consistent_datapipe(self):
        datapipe = vtkPolyDataPipeline(self.vtk_iso)
        self.assertTrue(datapipe.hasdata())
        self.assertTrue(self.vtk_iso in datapipe)
        self.assertEqual(datapipe.vtkish_data, self.vtk_iso)
        self.assertFalse(datapipe.hasfilters())
        self.assertFalse(datapipe.isclosed())

        return datapipe

    def test_consistent_datapipe_oninit(self):
        datapipe = self.get_consistent_datapipe()

        pipeline = vtkPolyDataPipeline(datapipe)
        self.assertTrue(pipeline.hasdata())
        self.assertTrue(datapipe in pipeline)
        self.assertTrue(self.vtk_iso in pipeline)
        self.assertEqual(pipeline.vtkish_data, datapipe)
        self.assertRaises(RuntimeError, datapipe.isclosed)
        self.assertFalse(pipeline.hasfilters())
        self.assertFalse(pipeline.isclosed())

        pipeline.append(self.vtk_dnorm)
        self.assertConnected(self.vtk_iso, self.vtk_dnorm)

        pipeline.append(self.vtk_subdiv)
        self.assertConnected(self.vtk_dnorm, self.vtk_subdiv)

        pipeline.append(self.vtk_dmap)
        self.assertConnected(self.vtk_subdiv, self.vtk_dmap)

        self.assertTrue(pipeline.hasfilters())
        self.assertTrue(pipeline.isclosed())

    def test_consistent_datapipe_direct(self):
        datapipe = self.get_consistent_datapipe()

        pipeline = vtkPolyDataPipeline()
        self.assertFalse(pipeline.hasdata())
        self.assertFalse(datapipe in pipeline)
        self.assertFalse(self.vtk_iso in pipeline)
        self.assertFalse(pipeline.hasfilters())
        self.assertFalse(pipeline.isclosed())

        pipeline.set_data(datapipe)
        self.assertTrue(pipeline.hasdata())
        self.assertTrue(datapipe in pipeline)
        self.assertTrue(self.vtk_iso in pipeline)
        self.assertEqual(pipeline.vtkish_data, datapipe)
        self.assertRaises(RuntimeError, datapipe.isclosed)

        pipeline.append(self.vtk_dnorm)
        self.assertConnected(self.vtk_iso, self.vtk_dnorm)

        pipeline.append(self.vtk_subdiv)
        self.assertConnected(self.vtk_dnorm, self.vtk_subdiv)

        pipeline.append(self.vtk_dmap)
        self.assertConnected(self.vtk_subdiv, self.vtk_dmap)

        self.assertTrue(pipeline.hasfilters())
        self.assertTrue(pipeline.isclosed())

    def test_consistent_datapipe_postponed(self):
        datapipe = self.get_consistent_datapipe()

        pipeline = vtkPolyDataPipeline()
        self.assertFalse(pipeline.hasdata())
        self.assertFalse(datapipe in pipeline)
        self.assertFalse(self.vtk_iso in pipeline)
        self.assertFalse(pipeline.hasfilters())
        self.assertFalse(pipeline.isclosed())

        pipeline.append(self.vtk_dnorm)
        self.assertTrue(pipeline.hasfilters())
        self.assertNotConnected(self.vtk_iso, self.vtk_dnorm)

        pipeline.append(self.vtk_subdiv)
        self.assertConnected(self.vtk_dnorm, self.vtk_subdiv)

        pipeline.set_data(datapipe)
        self.assertTrue(pipeline.hasdata())
        self.assertTrue(datapipe in pipeline)
        self.assertTrue(self.vtk_iso in pipeline)
        self.assertEqual(pipeline.vtkish_data, datapipe)
        self.assertConnected(self.vtk_iso, self.vtk_dnorm)
        self.assertTrue(datapipe.isclosed())

        pipeline.append(self.vtk_dmap)
        self.assertConnected(self.vtk_subdiv, self.vtk_dmap)
        self.assertTrue(pipeline.hasdata())
        self.assertTrue(pipeline.isclosed())

    def test_consistent_datapipe_onclose(self):
        datapipe = self.get_consistent_datapipe()

        pipeline = vtkPolyDataPipeline()
        self.assertFalse(pipeline.hasdata())
        self.assertFalse(datapipe in pipeline)
        self.assertFalse(self.vtk_iso in pipeline)
        self.assertFalse(pipeline.hasfilters())
        self.assertFalse(pipeline.isclosed())

        pipeline.append(self.vtk_dnorm)
        self.assertTrue(pipeline.hasfilters())
        self.assertNotConnected(self.vtk_iso, self.vtk_dnorm)

        pipeline.append(self.vtk_subdiv)
        self.assertConnected(self.vtk_dnorm, self.vtk_subdiv)

        pipeline.append(self.vtk_dmap)
        self.assertConnected(self.vtk_subdiv, self.vtk_dmap)
        self.assertTrue(pipeline.isclosed())

        pipeline.set_data(datapipe)
        self.assertTrue(pipeline.hasdata())
        self.assertEqual(pipeline.vtkish_data, datapipe)
        self.assertConnected(self.vtk_iso, self.vtk_dnorm)
        self.assertTrue(datapipe.isclosed())

    # =================================

    def get_consistent_filterpipe(self):
        filterpipe = vtkPolyDataPipeline()
        self.assertFalse(filterpipe.hasdata())
        self.assertFalse(self.vtk_iso in filterpipe)
        self.assertFalse(filterpipe.hasfilters())
        self.assertFalse(filterpipe.isclosed())

        filterpipe.append(self.vtk_dnorm)
        self.assertTrue(filterpipe.hasfilters())
        self.assertNotConnected(self.vtk_iso, self.vtk_dnorm)

        filterpipe.append(self.vtk_subdiv)
        self.assertConnected(self.vtk_dnorm, self.vtk_subdiv)

        return filterpipe

    def test_consistent_filterpipe_oninit(self):
        filterpipe = self.get_consistent_filterpipe()

        pipeline = vtkPolyDataPipeline(self.vtk_iso)
        self.assertTrue(pipeline.hasdata())
        self.assertEqual(pipeline.vtkish_data, self.vtk_iso)
        self.assertFalse(pipeline.hasfilters())
        self.assertFalse(pipeline.isclosed())

        pipeline.append(filterpipe)
        self.assertRaises(RuntimeError, filterpipe.isclosed)
        self.assertConnected(self.vtk_iso, self.vtk_dnorm)
        self.assertConnected(self.vtk_dnorm, self.vtk_subdiv)
        self.assertFalse(pipeline.isclosed())

        pipeline.append(self.vtk_dmap)
        self.assertTrue(filterpipe.isclosed())
        self.assertConnected(self.vtk_subdiv, self.vtk_dmap)
        self.assertTrue(pipeline.hasfilters())
        self.assertTrue(pipeline.isclosed())

    def test_consistent_filterpipe_direct(self):
        filterpipe = self.get_consistent_filterpipe()

        pipeline = vtkPolyDataPipeline()
        self.assertFalse(pipeline.hasdata())
        self.assertFalse(pipeline.hasfilters())
        self.assertFalse(pipeline.isclosed())

        pipeline.set_data(self.vtk_iso)
        self.assertTrue(pipeline.hasdata())
        self.assertEqual(pipeline.vtkish_data, self.vtk_iso)

        pipeline.append(filterpipe)
        self.assertRaises(RuntimeError, filterpipe.isclosed)
        self.assertConnected(self.vtk_iso, self.vtk_dnorm)
        self.assertConnected(self.vtk_dnorm, self.vtk_subdiv)
        self.assertFalse(pipeline.isclosed())

        pipeline.append(self.vtk_dmap)
        self.assertTrue(filterpipe.isclosed())
        self.assertConnected(self.vtk_subdiv, self.vtk_dmap)
        self.assertTrue(pipeline.hasfilters())
        self.assertTrue(pipeline.isclosed())

    def test_consistent_filterpipe_postponed(self):
        filterpipe = self.get_consistent_filterpipe()

        pipeline = vtkPolyDataPipeline()
        self.assertFalse(pipeline.hasdata())
        self.assertFalse(pipeline.hasfilters())
        self.assertFalse(pipeline.isclosed())

        pipeline.append(filterpipe)
        self.assertRaises(RuntimeError, filterpipe.isclosed)
        self.assertNotConnected(self.vtk_iso, self.vtk_dnorm)
        self.assertConnected(self.vtk_dnorm, self.vtk_subdiv)

        pipeline.set_data(self.vtk_iso)
        self.assertTrue(pipeline.hasdata())
        self.assertEqual(pipeline.vtkish_data, self.vtk_iso)
        self.assertConnected(self.vtk_iso, self.vtk_dnorm)

        pipeline.append(self.vtk_dmap)
        self.assertTrue(filterpipe.isclosed())
        self.assertConnected(self.vtk_subdiv, self.vtk_dmap)
        self.assertTrue(pipeline.hasdata())
        self.assertTrue(pipeline.isclosed())

    def test_consistent_filterdata_onclose(self):
        filterpipe = self.get_consistent_filterpipe()

        pipeline = vtkPolyDataPipeline()
        self.assertFalse(pipeline.hasdata())
        self.assertFalse(pipeline.hasfilters())
        self.assertFalse(pipeline.isclosed())

        pipeline.append(filterpipe)
        self.assertRaises(RuntimeError, filterpipe.isclosed)
        self.assertNotConnected(self.vtk_iso, self.vtk_dnorm)
        self.assertConnected(self.vtk_dnorm, self.vtk_subdiv)

        pipeline.append(self.vtk_dmap)
        self.assertTrue(filterpipe.isclosed())
        self.assertConnected(self.vtk_subdiv, self.vtk_dmap)
        self.assertFalse(pipeline.hasdata())
        self.assertTrue(pipeline.isclosed())

        pipeline.set_data(self.vtk_iso)
        self.assertTrue(pipeline.hasdata())
        self.assertEqual(pipeline.vtkish_data, self.vtk_iso)
        self.assertConnected(self.vtk_iso, self.vtk_dnorm)

    # =================================

    def test_failure_recursive_data(self):
        pipeline = vtkPolyDataPipeline()
        self.assertRaises(ValueError, pipeline.set_data, pipeline)

    def test_failure_recursive_mixed(self):
        pipeline = vtkPolyDataPipeline(self.vtk_iso)
        pipeline.append(self.vtk_dnorm)
        self.assertRaises(ValueError, pipeline.append, pipeline)

    def test_failure_recursive_filter(self):
        pipeline = vtkPolyDataPipeline()
        self.assertRaises(ValueError, pipeline.append, pipeline)


# -------------------------------------------------------------------

if __name__ in ['__main__', '__builtin__']:
    # We may have been imported by test.py, if so we should redirect to logfile
    if __name__ == '__builtin__':
        testrunner = CustomTextTestRunner('vtk_pipeline.log', verbosity=2)
    else:
        testrunner = unittest.TextTestRunner(stream=sys.stdout, verbosity=2)

    testcases = [UTPolyDataPipeline_PureVTK, UTPolyDataPipeline_PipelineVTK]

    for test in testcases:
        info = '\n' + test.__name__ + '\n' + test.__doc__.strip('\n') + '\n'
        testsuite = unittest.defaultTestLoader.loadTestsFromTestCase(test)
        testrunner.stream.writeln(info)
        testresult = testrunner.run(testsuite)
        # Provide feedback on failed tests if imported by test.py
        if __name__ == '__builtin__' and not testresult.wasSuccessful():
            raise SystemExit('Test failed. Check vtk_pipeline.log for details.')

