<?xml version="1.0" encoding="UTF-8" ?>
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">
  <xsl:output method="text" cdata-section-elements="at"/>
  <xsl:template match="/">
  <xsl:variable name="newline"><xsl:text>
</xsl:text></xsl:variable>
    <xsl:text>BEGIN_INFO
 # Band-XCRYSDEN-Structure-File for Fermi surface plotting
 # created by xmlfermis2bxsf.xsd
 # Launch as: xcrysden --bxsf filename.bxsf
 # numberg of grid points </xsl:text>
    <xsl:value-of select="count(/fermisurface/point)"/>
    <xsl:text>
   Fermi Energy:    0.000000000    
 END_INFO
 BEGIN_BLOCK_BANDGRID_3D
 band_energies
 BANDGRID_3D_BANDS
   </xsl:text>
    <xsl:value-of select="count(/fermisurface/point[1]/band)"/>
    <xsl:text>
    </xsl:text>
    <xsl:value-of select="/fermisurface/runitcell/@grid"/>
    <xsl:text>
   0.000000000       0.000000000       0.000000000
</xsl:text>
    <xsl:for-each select="/fermisurface/runitcell/bvec">
      <xsl:text>   </xsl:text>
      <xsl:value-of select="."/>
      <xsl:value-of select="$newline"/>
    </xsl:for-each>
    <xsl:for-each select="/fermisurface/point[1]/band">
      <xsl:text> BAND:   </xsl:text>
      <xsl:value-of select="@nr"/>
      <xsl:value-of select="$newline"/>
     
      <xsl:variable name="nr" select="@nr"/>
      <xsl:for-each select="/fermisurface/point/band[@nr=$nr]">
      <xsl:text> </xsl:text>   
      <xsl:value-of select="@eval"/>
      <xsl:value-of select="$newline"/>
      </xsl:for-each>
    </xsl:for-each>
    <xsl:text> END_BANDGRID_3D
 END_BLOCK_BANDGRID_3D</xsl:text>
  </xsl:template>
</xsl:stylesheet>
  