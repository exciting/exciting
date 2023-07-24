""" Automatic generation of all standard input classes plus definiton of exceptions. """

from pathlib import Path
from typing import Optional, List, Callable
from typing import Union
from xml.etree import ElementTree

import numpy as np

from excitingtools.input.base_class import AbstractExcitingInput, ExcitingXMLInput
from excitingtools.input.dynamic_class import generate_classes_str
from excitingtools.input.xml_utils import xml_tree_to_pretty_str, prettify_tag_attributes
from excitingtools.utils.dict_utils import check_valid_keys
from excitingtools.utils.utils import list_to_str
from excitingtools.utils.valid_attributes import valid_plan_entries

# define names of classes which are meant to be available for a user or used directly elsewhere in excitingtools
ExcitingCrystalInput: Callable
ExcitingSpeciesInput: Callable
ExcitingGroundStateInput: Callable
ExcitingXSInput: Callable
ExcitingPropertiesInput: Callable
ExcitingPointInput: Callable
ExcitingBandStructureInput: Callable

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


class ExcitingInputXML(ExcitingXMLInput):
    """
    Container for a complete input xml file.
    """
    name = "input"
    _default_filename = "input.xml"

    def set_title(self, title: str):
        """ Set a new title. """
        self.__dict__["title"].title = title

    def to_xml_str(self) -> str:
        """Compose XML ElementTrees from exciting input classes to create an input xml string.

        :return: Input XML tree as a string, with pretty formatting.
        """
        xml_tree = self.to_xml()
        tags_to_prettify = ["\t<structure", "\t\t<crystal", "\t\t<species", "\t\t\t<atom", "\t<groundstate", "\t<xs",
                            "\t\t<screening", "\t\t<BSE", "\t\t<energywindow"]
        input_xml_str = prettify_tag_attributes(xml_tree_to_pretty_str(xml_tree), tags_to_prettify)
        return input_xml_str

    def write(self, filename: Union[str, Path] = _default_filename):
        """Writes the xml string to file.

        :param filename: name of the file.
        """
        with open(filename, "w") as fid:
            fid.write(self.to_xml_str())


class ExcitingQpointsetInput(AbstractExcitingInput):
    """
    Class for exciting Qpointset Input
    """
    name = "qpointset"

    def __init__(self, qpointset: Optional[Union[np.ndarray, List[List[float]]]] = np.array([0.0, 0.0, 0.0])):
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


class ExcitingPathInput(ExcitingXMLInput):
    """
    Class for exciting path input.
    
    Note: Not to be confused by the class name, it is NOT directly a band path. Defining what is in the
    exciting input file under the subtree called 'path'.
    """
    name = "path"

    def __init__(self,
                 points: List[Union[dict, ExcitingPointInput]],
                 **kwargs):
        """Generate an object of ExcitingXMLInput for the path attributes."""
        super().__init__(**kwargs)
        self.__dict__["points"] = [self._initialise_subelement_attribute(ExcitingPointInput, point) for point in points]

    def to_xml(self) -> ElementTree:
        """Put class attributes into an XML tree, with the element given by self.name.
        :return ElementTree.Element sub_tree: sub_tree element tree, with class attributes inserted.
        """
        attributes = {key: str(value) for key, value in vars(self).items() if key != "points"}
        xml_tree = ElementTree.Element(self.name, **attributes)

        for point in self.points:
            xml_tree.append(point.to_xml())

        return xml_tree
