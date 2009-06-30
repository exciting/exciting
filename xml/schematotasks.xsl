<?xml version="1.0" encoding="UTF-8" ?>
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform" xmlns:xs="http://www.w3.org/2001/XMLSchema">
  <xsl:output method="text"/>
  <xsl:variable name="newline">
    <xsl:text>
    </xsl:text>
  </xsl:variable>
  <xsl:output method="text"></xsl:output>
  <xsl:template match="/">
  <xsl:for-each select="//xs:element[@name='properties']/*/*/xs:element/@name">
  <xsl:text>if(associated(input%properties%</xsl:text>
  <xsl:value-of select="." />
    <xsl:text>)) then</xsl:text>
  <xsl:value-of select="$newline" /> 
   <xsl:text>endif</xsl:text>  <xsl:value-of select="$newline" /> 
  </xsl:for-each>
  </xsl:template>
  </xsl:stylesheet>