<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
  xmlns:exsl="http://exslt.org/common" xmlns:math="http://exslt.org/math"
  xmlns:set="http://exslt.org/sets" extension-element-prefixes="exsl math set" version="1.0">
  <xsl:output method="text" cdata-section-elements="at"/>
  
  <xsl:template match="/">
    
    <xsl:variable name="newline">
      <xsl:text>
</xsl:text>
    </xsl:variable>
    
    <xsl:text>BEGIN_INFO
 # Band-XCRYSDEN-Structure-File for Fermi surface plotting
 # created by xmlfermis2bxsf.xsl
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
    <xsl:choose>
      <xsl:when test="/fermisurface/point/band">
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
            <xsl:text>  </xsl:text>
            <xsl:value-of select="@eval"/>
            <xsl:value-of select="$newline"/>
          </xsl:for-each>
        </xsl:for-each>
      </xsl:when>
      
      <xsl:otherwise>
        <!-- test="not(/fermisurface/point/band)" -->
        <xsl:value-of select="1"/>
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
        <xsl:choose>
          <xsl:when test="/fermisurface/point/@product">
            <xsl:text> prod:   </xsl:text>
            <xsl:value-of select="1"/>
            <xsl:value-of select="$newline"/>
            <xsl:for-each select="fermisurface/point">
              <xsl:value-of select="@product"/>
              <xsl:value-of select="$newline"/>
            </xsl:for-each>
          </xsl:when>
          <xsl:otherwise>
            <xsl:text> prod:up</xsl:text>
            <xsl:value-of select="$newline"/>
            <xsl:for-each select="/fermisurface/point">
              <xsl:value-of select="@up"/>
              <xsl:value-of select="$newline"/>
            </xsl:for-each>
            <xsl:text> prod:down   </xsl:text>
            <xsl:value-of select="$newline"/>
            <xsl:for-each select="/fermisurface/point">
              <xsl:value-of select="@down"/>
              <xsl:value-of select="$newline"/>
            </xsl:for-each>
          </xsl:otherwise>
        </xsl:choose>
        
      </xsl:otherwise>
    </xsl:choose>
    <xsl:text> END_BANDGRID_3D
 END_BLOCK_BANDGRID_3D</xsl:text>
    
  </xsl:template>
</xsl:stylesheet>
