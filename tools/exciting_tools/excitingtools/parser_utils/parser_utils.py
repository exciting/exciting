"""General parser utility functions.
"""
from xml.etree import ElementTree


def find_element(root: ElementTree.Element, tag: str) -> ElementTree.Element:
    """ Finds a given tag in an Element, either return the full ElementTree
    if the tag is correct or find the tag in that ElementTree.
    :param root: Element to find the tag in
    :param tag: tag to search for
    :returns: the (sub)element with the correct tag
    """
    if root.tag == tag:
        return root
    return root.find(tag)
