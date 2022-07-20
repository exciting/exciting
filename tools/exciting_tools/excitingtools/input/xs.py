"""Module for class of exciting xs (excited states).
http://exciting.wikidot.com/ref:xs
"""
from typing import Optional, List, Union
from xml.etree import ElementTree as ET

import numpy as np
from excitingtools.input.base_class import ExcitingXMLInput
from excitingtools.utils.utils import list_to_str
from excitingtools.utils.dict_utils import check_valid_keys


class ExcitingXSInput(ExcitingXMLInput):
    """ Class allowing to write attributes to XML."""

    # TODO(Fabian): Add all the other subelements, see http://exciting.wikidot.com/ref:xs
    # Issue 121: https://git.physik.hu-berlin.de/sol/exciting/-/issues/121
    _valid_xs_attributes = {'bfieldc', 'broad', 'dbglev', 'dfoffdiag', 'dogroundstate', 'emattype', 'emaxdf',
                            'epsdfde', 'fastpmat', 'gqmax', 'gqmaxtype', 'lmaxapwwf', 'lmaxemat', 'maxscl', 'nempty',
                            'ngridk', 'ngridq', 'nosym', 'pwmat', 'reducek', 'reduceq', 'rgkmax', 'scissor', 'skipgnd',
                            'swidth', 'tappinfo', 'tevout', 'vkloff', 'writexsgrids', 'xstype'}

    @staticmethod
    def _initialise_subelement_attribute(XSClass, element):
        """
        Initialize given elements to the ExcitingXSInput constructor. If element is already ExcitingXMLInput class
        object, nothing happens. For None elements None is returned. In any other case, the class constructor of the
        given XSClass is called.
        """
        if isinstance(element, ExcitingXMLInput):
            return element
        elif element is None:
            return None
        else:
            # Assume the element type is valid for the class constructor
            return XSClass(element)

    def __init__(self, xstype: str,
                 xs: Optional[dict] = None,
                 screening: Optional[Union[dict, ExcitingXMLInput]] = None,
                 BSE: Optional[Union[dict, ExcitingXMLInput]] = None,
                 qpointset: Optional[Union[np.ndarray, List[List[float]]]] = None,
                 energywindow: Optional[Union[dict, ExcitingXMLInput]] = None,
                 plan: Optional[List[str]] = None):
        """
        Initialize instance of ExcitingXS
        """
        if xs is None:
            xs = {}
        self.xs = ExcitingXMLInput('xs', self._valid_xs_attributes, xstype=xstype, **xs)
        self.screening = self._initialise_subelement_attribute(ExcitingXSScreeningInput, screening)
        self.BSE = self._initialise_subelement_attribute(ExcitingXSBSEInput, BSE)
        self.energywindow = self._initialise_subelement_attribute(ExcitingXSEnergywindowInput, energywindow)
        self.qpointset = self._initialise_subelement_attribute(ExcitingXSQpointsetInput, qpointset)
        self.plan = self._initialise_subelement_attribute(ExcitingXSPlanInput, plan)

    def to_xml(self) -> ET.Element:
        """Put class attributes into an XML tree, 'xs'.
        """
        xs_tree = self.xs.to_xml()

        attributes = list(filter(lambda x: isinstance(x, ExcitingXMLInput), vars(self).values()))
        attributes = list(filter(lambda x: x.name != 'xs', attributes))
        for attribute in attributes:
            xs_tree.append(attribute.to_xml())

        return xs_tree


class ExcitingXSBSEInput(ExcitingXMLInput):
    """
    Class for exciting BSE Input
    """
    _valid_BSE_attributes = {'aresbse', 'blocks', 'bsedirsing', 'bsetype', 'checkposdef', 'chibar0', 'chibar0comp',
                             'chibarq', 'coupling', 'cuttype', 'distribute', 'econv', 'eecs', 'efind', 'fbzq',
                             'iqmtrange', 'lmaxdielt', 'measure', 'nexc', 'ngridksub', 'nleblaik', 'nosym', 'nstlbse',
                             'nstlxas', 'outputlevel', 'reducek', 'rgkmax', 'sciavbd', 'sciavqbd', 'sciavqhd',
                             'sciavqwg', 'sciavtype', 'scrherm', 'vkloff', 'writehamhdf5', 'writepotential', 'xas',
                             'xasatom', 'xasedge', 'xasspecies', 'xes'}

    def __init__(self, BSE: dict):
        super().__init__('BSE', self._valid_BSE_attributes, **BSE)


class ExcitingXSScreeningInput(ExcitingXMLInput):
    """
    Class for exciting Screening Input
    """
    _valid_screening_attributes = {'do', 'intraband', 'nempty', 'ngridk', 'nosym', 'reducek', 'rgkmax', 'screentype',
                                   'tr', 'vkloff'}

    def __init__(self, screening: dict):
        super().__init__('screening', self._valid_screening_attributes, **screening)


class ExcitingXSEnergywindowInput(ExcitingXMLInput):
    """
    Class for exciting Energywindow Input
    """
    _valid_energywindow_attributes = {'intv', 'points'}

    def __init__(self, energywindow: dict):
        super().__init__('energywindow', **energywindow)


class ExcitingXSQpointsetInput(ExcitingXMLInput):
    """
    Class for exciting Qpointset Input
    """

    def __init__(self, qpointset: Optional[Union[np.ndarray, List[List[float]]]] = np.array([0.0, 0.0, 0.0])):
        """
        Qpointset should be passed either as numpy array or as a list of lists, so either
        np.array([[0., 0., 0.], [0.0, 0.0, 0.01], ...])
        or
        [[0., 0., 0.], [0.0, 0.0, 0.01], ...]
        """
        super().__init__('qpointset')
        self.qpointset = qpointset

    def to_xml(self) -> ET.Element:
        qpointset = ET.Element('qpointset')
        for qpoint in self.qpointset:
            ET.SubElement(qpointset, 'qpoint').text = list_to_str(qpoint)

        return qpointset


class ExcitingXSPlanInput(ExcitingXMLInput):
    """
    Class for exciting Plan Input
    """
    _valid_plan_elements = {'xsgeneigvec', 'tetcalccw', 'writepmatxs', 'writeemat', 'df', 'df2', 'idf', 'scrgeneigvec',
                            'scrtetcalccw', 'scrwritepmat', 'screen', 'scrcoulint', 'exccoulint', 'bse', 'bsegenspec',
                            'writebevec', 'writekpathweights', 'bsesurvey', 'kernxc_bse', 'writebandgapgrid',
                            'write_wfplot', 'write_screen', 'writepmat', 'dielectric', 'writepmatasc', 'pmatxs2orig',
                            'writeoverlapxs', 'writeematasc', 'writepwmat', 'emattest', 'x0toasc', 'x0tobin',
                            'fxc_alda_check', 'kernxc_bse3', 'testxs', 'xsestimate', 'testmain', 'excitonWavefunction',
                            'portstate(1)', 'portstate(2)', 'portstate(-1)', 'portstate(-2)'}

    def __init__(self, plan: List[str]):
        """
        Plan doonly elements are passed as a List of strings in the order exciting shall execute them:
            ['bse', 'xseigval', ...]
        """
        super().__init__('plan')
        check_valid_keys(plan, self._valid_plan_elements, 'Plan')
        self.plan = plan

    def to_xml(self) -> ET.Element:
        plan = ET.Element('plan')
        for task in self.plan:
            ET.SubElement(plan, 'doonly', task=task)

        return plan
