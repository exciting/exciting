<?xml version="1.0" encoding="UTF-8" ?>
  <!--
    ##################################################################
     Use: xsltproc [dash dash]path "../../species/" exciting2sgroup.xsl input.xml 
     ###################################################################
  -->
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform" xmlns:str="http://exslt.org/strings" xmlns:math="http://exslt.org/math">
  <xsl:output method="text" />
  <xsl:include href="basevec2abc.xsl" />
  <xsl:template match="/">
    <xsl:text>P
</xsl:text>
    <xsl:value-of select="$a" />
    <xsl:text> </xsl:text>
    <xsl:value-of select="$b" />
    <xsl:text> </xsl:text>
    <xsl:value-of select="$c" />
    <xsl:text> </xsl:text>
    <xsl:value-of select="$alpha" />
    <xsl:text> </xsl:text>
    <xsl:value-of select="$beta" />
    <xsl:text> </xsl:text>
    <xsl:value-of select="$gamma" />
    <xsl:text>
    
</xsl:text>
    <xsl:value-of select="count(//atom)" />
    <xsl:text>
</xsl:text>
    <xsl:for-each select="//species">
      <xsl:variable name="speciesfile">
        <xsl:value-of select="../@speciespath" />
        <xsl:value-of select="@speciesfile" />
      </xsl:variable>
      <xsl:for-each select="atom">
        <xsl:value-of select="@coord" />
        <xsl:text>
</xsl:text>
        <xsl:value-of select="../@speciesfile" />
        <xsl:text>
</xsl:text>
      </xsl:for-each>
    </xsl:for-each>
  </xsl:template>
  
</xsl:stylesheet>