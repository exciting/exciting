from lxml import etree

def validate(xmlPath, xsdPath):
    xmlSchema_doc = etree.parse(xsdPath)
    xmlSchema = etree.XMLSchema(xmlSchema_doc)

    xml_doc = etree.parse(xmlPath)
    result = xmlSchema.validate(xml_doc)

    return result