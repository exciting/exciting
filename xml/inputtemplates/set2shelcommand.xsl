<?xml version="1.0" encoding="UTF-8" ?>
<xsl:stylesheet version="1.0"
    xmlns:xsl="http://www.w3.org/1999/XSL/Transform">
<xsl:output method="text" />
<xsl:template match="/">
 
<!-- Authors: chm, jus, sag -->
 
<!--/////////////////////////////////////////////////////////////////////////-->
<!--/// define the name of the input file ///////////////////////////////////-->
<!--/////////////////////////////////////////////////////////////////////////-->
<xsl:variable name="inputfilename"><xsl:text></xsl:text></xsl:variable>
 
<!-- Loop over all elements named "set" from reference xml-file -->
<xsl:for-each select = "/experiment/set">
 
  <!-- Define path here -->
  <xsl:variable name="path">
          <xsl:value-of select="@path"/>
          <xsl:text>/</xsl:text>
  </xsl:variable>
 
 <!-- Write document at Path $path -->
<xsl:text>cd </xsl:text> <xsl:value-of select="$path"/>
<xsl:text> 
excitingser
cd -
</xsl:text>
</xsl:for-each>
</xsl:template>
</xsl:stylesheet>
