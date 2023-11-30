""" Input xml class. """
from pathlib import Path
from typing import Union

from excitingtools.input.base_class import ExcitingXMLInput
from excitingtools.input.input_classes import ExcitingGroundStateInput, ExcitingTitleInput
from excitingtools.input.structure import ExcitingStructure
from excitingtools.input.xml_utils import prettify_tag_attributes, xml_tree_to_pretty_str


class ExcitingInputXML(ExcitingXMLInput):
    """
    Container for a complete input xml file.
    """
    name = "input"
    _default_filename = "input.xml"
    structure: ExcitingStructure
    groundstate: ExcitingGroundStateInput
    title: ExcitingTitleInput

    def set_title(self, title: str):
        """ Set a new title. """
        self.__dict__["title"].title = title

    def to_xml_str(self) -> str:
        """Compose XML ElementTrees from exciting input classes to create an input xml string.

        :return: Input XML tree as a string, with pretty formatting.
        """
        return prettify_tag_attributes(xml_tree_to_pretty_str(self.to_xml()))

    def write(self, filename: Union[str, Path] = _default_filename):
        """Writes the xml string to file.

        :param filename: name of the file.
        """
        with open(filename, "w") as fid:
            fid.write(self.to_xml_str())
