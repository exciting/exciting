<?xml version="1.0" encoding="UTF-8" ?>
<xsl:stylesheet version="1.0"
    xmlns:xsl="http://www.w3.org/1999/XSL/Transform">
<xsl:output method="xml" />
<xsl:template match="/">
 
<!-- Authors: chm, jus, sag -->
 
<xsl:variable name="inputfilename"><xsl:text>input.xml</xsl:text></xsl:variable>

<!-- Loop over all elements named "set" from reference xml-file -->
<xsl:for-each select = "/experiment/set">
<xsl:variable name="path"> 
<xsl:value-of select="@path"/>
<xsl:text>/</xsl:text>
<xsl:value-of select="$inputfilename"></xsl:value-of>
</xsl:variable> 

 <!-- Write document at Path $path -->
 <xsl:document href="{$path}" method="xml" indent="yes">
 <xsl:comment>
 This file is generated with XSLTPROC using a template file and a reference file
 All parameters from the set filled in by XSLTPROC are listed below:
<!-- list all attributes in input file -->
<xsl:for-each select="./@*">
<xsl:text/>   
<xsl:value-of select="name()"/><xsl:text/> 
<xsl:value-of select="."/><xsl:text/>

</xsl:for-each><xsl:text/>  
</xsl:comment>
<!-- ////////////////////////////////////////////////////////////////////////-->
<!-- /// The input file begins here /////////////////////////////////////////-->
<!-- ////////////////////////////////////////////////////////////////////////-->
<input xsi:noNamespaceSchemaLocation="excitinginput.xsd"
  xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsltpath="../../../xml/">
  <title>Ag Energy-Volume</title>
  <structure speciespath="../../../species/">
    <crystal >
    <xsl:attribute name="scale">
    <xsl:value-of select="3.86*@unitcellscaling"/>
    </xsl:attribute>

      <basevect>1.0 1.0 0.0</basevect>
      <basevect>1.0 0.0 1.0</basevect>
      <basevect>0.0 1.0 1.0</basevect>
    </crystal>
    <species speciesfile="Ag.xml">
      <atom coord="0.0  0.0  0.0" />
    </species>
  </structure>
  <groundstate ngridk="8 8 8" rgkmax="8" />
</input>

<!-- ////////////////////////////////////////////////////////////////////////-->
<!-- /// The input file ends here ///////////////////////////////////////////-->
<!-- ////////////////////////////////////////////////////////////////////////-->
 
</xsl:document>
</xsl:for-each>
</xsl:template>
</xsl:stylesheet>

