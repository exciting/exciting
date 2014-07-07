<?xml version="1.0" encoding="UTF-8" ?>
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform" xmlns:str="http://exslt.org/strings" >
  <xsl:output method="text" cdata-section-elements="at" />
  <xsl:template match="/">
   
  <xsl:for-each select="/plot1d/grid/point">
    <xsl:value-of select="@distance" />
    <xsl:text> </xsl:text>
    <xsl:value-of select="@function1" />
    <xsl:text>
</xsl:text>
  </xsl:for-each>
  
  </xsl:template>
</xsl:stylesheet>
