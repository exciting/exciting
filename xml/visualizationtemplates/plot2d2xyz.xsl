<?xml version="1.0" encoding="UTF-8" ?>
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform" xmlns:str="http://exslt.org/strings" >
  <xsl:output method="text" cdata-section-elements="at" />
  <xsl:variable name="da" select="str:tokenize(/plot2d/grid/axis[1]/@deltas)" />
  <xsl:variable name="db" select="str:tokenize(/plot2d/grid/axis[2]/@deltas)" />
  <xsl:variable name="grid" select="str:tokenize(/plot3d/grid/@gridticks)" />
  <xsl:template match="/">
    <xsl:for-each select="/plot2d/function">
      <xsl:text>#grid </xsl:text>
      <xsl:value-of select="/plot2d/grid/@gridticks"/>
      <xsl:text>
</xsl:text>
      <xsl:for-each select="/plot2d/function/row">
        <xsl:variable name="indexx" select="@index"/>
        <xsl:for-each select="str:tokenize(.)">
          <xsl:variable name="indexy" select="position()"/>
          <xsl:value-of select="$da*$indexx"/>
          <xsl:text> </xsl:text>
          <xsl:value-of select="$db*$indexy"/>
          <xsl:text> </xsl:text>
          <xsl:value-of select="."/>
          <xsl:text>
</xsl:text>
        </xsl:for-each>
      </xsl:for-each>
      <xsl:text>
</xsl:text>      
    </xsl:for-each>
    
  </xsl:template>
</xsl:stylesheet>
