<?xml version="1.0" encoding="UTF-8" ?>
<xsl:stylesheet version="1.0"
  xmlns:xsl="http://www.w3.org/1999/XSL/Transform" xmlns:math="http://exslt.org/math"
  xmlns:str="http://exslt.org/strings">

  <xsl:output method="xml" indent="yes" />
  <xsl:template match="/">
    <graph>
      <xsl:for-each select="/experiment/set">
        <!-- Define path here -->
        <xsl:variable name="path">
          <xsl:value-of select="@path" />
          <xsl:text>/</xsl:text>
        </xsl:variable>
        <!--collect data -->
        <xsl:variable name="infopath">
          <xsl:value-of select="$path" />
          <xsl:text>info.xml</xsl:text>
        </xsl:variable>
        <xsl:variable name="inputpath">
          <xsl:value-of select="$path" />
          <xsl:text>input.xml</xsl:text>
        </xsl:variable>
        <point>
          <xsl:attribute name="volume">
          <!-- to calculate the volume we  -->
<xsl:value-of
            select="math:power(document($inputpath)//crystal/@scale,3) * math:power(math:sqrt(2.0),3)"></xsl:value-of>
 </xsl:attribute>
          <xsl:attribute name="scale">
 <xsl:value-of select="document($inputpath)//crystal/@scale"></xsl:value-of>
 </xsl:attribute>
          <xsl:attribute name="totalEnergy">
  <xsl:value-of
            select="document($infopath)//iter[last()]/energies/@totalEnergy"></xsl:value-of>
 </xsl:attribute>
        </point>
      </xsl:for-each>
    </graph>
  </xsl:template>
</xsl:stylesheet>

