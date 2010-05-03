<?xml version="1.0" encoding="UTF-8" ?>
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform" xmlns:str="http://exslt.org/strings">
  <xsl:output method="text" cdata-section-elements="at" />
  <xsl:variable name="da" select="str:tokenize(/plot3d/grid/axis[1]/@delta)" />
  <xsl:variable name="db" select="str:tokenize(/plot3d/grid/axis[2]/@delta)" />
  <xsl:variable name="dc" select="str:tokenize(/plot3d/grid/axis[3]/@delta)" />
  <xsl:variable name="grid" select="str:tokenize(/plot3d/grid/@gridticks)" />
  <xsl:template match="/">
    <xsl:text>
</xsl:text>
    <xsl:for-each select="/plot2d/function/row">
      <xsl:variable name="indexx" select="@index" />
      <xsl:value-of select="." />
      <xsl:text>
</xsl:text>
    </xsl:for-each>
  </xsl:template>
</xsl:stylesheet>