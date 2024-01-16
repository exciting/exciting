""" Automatic generation of all standard input classes plus definiton of exceptions. """
from typing import List, Union, Any
from xml.etree import ElementTree

import numpy as np

from excitingtools.input.base_class import AbstractExcitingInput
from excitingtools.input.dynamic_class import generate_classes_str
from excitingtools.utils.dict_utils import check_valid_keys
from excitingtools.utils.utils import list_to_str
from excitingtools.utils.valid_attributes import valid_plan_entries

# define names of classes which are meant to be available for a user or used directly elsewhere in excitingtools
# type hint as Any to not conflict with static type checkers as the input classes are generated dynamically
ExcitingCrystalInput: Any
ExcitingSpeciesInput: Any
ExcitingGroundStateInput: Any
ExcitingXSInput: Any
ExcitingPropertiesInput: Any
ExcitingPointInput: Any
ExcitingBandStructureInput: Any
ExcitingRelaxInput: Any
ExcitingPhononsInput: Any
ExcitingGWInput: Any
ExcitingMDInput: Any
ExcitingEPHInput: Any

# execute dynamically generated string with all standard class defintions
exec(generate_classes_str())


class ExcitingTitleInput(AbstractExcitingInput):
    """ Holds only the title but for consistency reasons as class. """
    name = "title"

    def __init__(self, title: str):
        self.title = title

    def to_xml(self, **kwargs) -> ElementTree:
        """ Puts title to xml, only the text is title. """
        title_tree = ElementTree.Element(self.name)
        title_tree.text = self.title
        return title_tree


class ExcitingQpointsetInput(AbstractExcitingInput):
    """
    Class for exciting Qpointset Input
    """
    name = "qpointset"

    def __init__(self, qpointset: Union[np.ndarray, List[List[float]]] = np.array([0.0, 0.0, 0.0])):
        """
        Qpointset should be passed either as numpy array or as a list of lists, so either
        np.array([[0., 0., 0.], [0.0, 0.0, 0.01], ...])
        or
        [[0., 0., 0.], [0.0, 0.0, 0.01], ...]
        """
        self.qpointset = qpointset

    def to_xml(self) -> ElementTree.Element:
        """ Special implementation of to_xml for the qpointset element. """
        qpointset = ElementTree.Element(self.name)
        for qpoint in self.qpointset:
            ElementTree.SubElement(qpointset, 'qpoint').text = list_to_str(qpoint)

        return qpointset


class ExcitingPlanInput(AbstractExcitingInput):
    """
    Class for exciting Plan Input
    """
    name = "plan"

    def __init__(self, plan: List[str]):
        """
        Plan doonly elements are passed as a List of strings in the order exciting shall execute them:
            ['bse', 'xseigval', ...]
        """
        check_valid_keys(plan, valid_plan_entries, self.name)
        self.plan = plan

    def to_xml(self) -> ElementTree.Element:
        """ Special implementation of to_xml for the plan element. """
        plan = ElementTree.Element(self.name)
        for task in self.plan:
            ElementTree.SubElement(plan, 'doonly', task=task)

        return plan


class ExcitingKstlistInput(AbstractExcitingInput):
    """
    Class for exciting Kstlist Input
    """
    name = "kstlist"

    def __init__(self, kstlist: Union[np.ndarray, List[List[int]]]):
        """
        Kstlist should be passed either as numpy array or as a list of lists, so either
        np.array([[1, 2], [3, 4], ...])
        or
        [[1, 2], [3, 4], ...]
        """
        self.kstlist = kstlist

    def to_xml(self) -> ElementTree.Element:
        """ Special implementation of to_xml for the kstlist element. """
        kstlist = ElementTree.Element(self.name)
        for pointstatepair in self.kstlist:
            ElementTree.SubElement(kstlist, 'pointstatepair').text = list_to_str(pointstatepair)

        return kstlist
