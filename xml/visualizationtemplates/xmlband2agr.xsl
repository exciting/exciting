<?xml version="1.0" encoding="UTF-8" ?>
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">
  <xsl:output method="text" />
  <xsl:template name="bandstograce" >
  <xsl:param name="data"/>
  <xsl:param name="l" select="false"/>
    <xsl:text>
@version 50122
@ default linewidth 2.0
@ page resize 600, 600
@ with line
@     line on
@     line loctype world
@     line g0
@     line 0, 0,</xsl:text> <xsl:value-of select="//point[last()]/@distance"></xsl:value-of> <xsl:text>, 0
@     line linewidth 2.0
@     line linestyle 3
@ line def
@ with string
@     string on
@     string loctype world
@     string g0
@     string </xsl:text> <xsl:value-of select="//point[last()]/@distance*1.01 "></xsl:value-of> <xsl:text>, -0.2
@     string font 12
@     string char size 1.650000
@     string def "E\sF"
@ r0 off
@ with g0
@     view 0.160000, 0.120000, 0.850000, .85
@     title "</xsl:text>
<xsl:value-of select="/bandstructure/title"/>
<xsl:if test="$l">
<xsl:text> character </xsl:text>
<xsl:value-of select="../@chemicalSymbol"/>
<xsl:text> l: </xsl:text>
<xsl:value-of select="$l"/>
</xsl:if>
<xsl:text>"  
@     world 0, -8, </xsl:text> <xsl:value-of select="//point[last()]/@distance"></xsl:value-of> <xsl:text>, 8  
@     yaxis  label "Energy (eV)"
@     yaxis  label char size 1.800000
@     yaxis  ticklabel char size 1.50000
@     xaxis  ticklabel char size 1.65000
@     xaxis  tick place both
@     xaxis  tick spec type both
@     xaxis  tick major grid on
@     xaxis  tick spec   </xsl:text>
    <xsl:value-of select="count(/bandstructure/vertex)" />
    <xsl:text>
</xsl:text>
    <xsl:for-each select="/bandstructure/vertex">
      <xsl:text>@     xaxis  tick major  </xsl:text>
      <xsl:value-of select="position()-1" />
      <xsl:text> , </xsl:text>
      <xsl:value-of select="@distance" />
      <xsl:text>
</xsl:text>
      <xsl:text>@     xaxis  ticklabel </xsl:text>
      <xsl:value-of select="position()-1" />
      <xsl:text>, "</xsl:text>
      <xsl:choose>
        <xsl:when test="@label='GAMMA'">
          <xsl:text>\xG\f{}</xsl:text>
        </xsl:when>
        <xsl:when test="@label='Gamma'">
          <xsl:text>\xG\f{}</xsl:text>
        </xsl:when>
        <xsl:when test="@label='DELTA'">
          <xsl:text>\xD\f{}</xsl:text>
        </xsl:when>
        <xsl:when test="@label='Delta'">
          <xsl:text>\xD\f{}</xsl:text>
        </xsl:when>
        <xsl:when test="@label='THETA'">
          <xsl:text>\xQ\f{}</xsl:text>
        </xsl:when>
        <xsl:when test="@label='Theta'">
          <xsl:text>\xQ\f{}</xsl:text>
        </xsl:when>
        <xsl:when test="@label='XI'">
          <xsl:text>\xX\f{}</xsl:text>
        </xsl:when>
        <xsl:when test="@label='Xi'">
          <xsl:text>\xX\f{}</xsl:text>
        </xsl:when>
        <xsl:when test="@label='PI'">
          <xsl:text>\xP\f{}</xsl:text>
        </xsl:when>
        <xsl:when test="@label='Pi'">
          <xsl:text>\xP\f{}</xsl:text>
        </xsl:when>
        <xsl:when test="@label='SIGMA'">
          <xsl:text>\xS\f{}</xsl:text>
        </xsl:when>
        <xsl:when test="@label='Sigma'">
          <xsl:text>\xS\f{}</xsl:text>
        </xsl:when>
        <xsl:when test="@label='PHI'">
          <xsl:text>\xF\f{}</xsl:text>
        </xsl:when>
        <xsl:when test="@label='Phi'">
          <xsl:text>\xF\f{}</xsl:text>
        </xsl:when>
        <xsl:when test="@label='PSI'">
          <xsl:text>\xY\f{}</xsl:text>
        </xsl:when>
        <xsl:when test="@label='Psi'">
          <xsl:text>\xY\f{}</xsl:text>
        </xsl:when>
        <xsl:when test="@label='OMEGA'">
          <xsl:text>\xW\f{}</xsl:text>
        </xsl:when>
        <xsl:when test="@label='Omega'">
          <xsl:text>\xW\f{}</xsl:text>
        </xsl:when>
        <xsl:otherwise>
          <xsl:value-of select="@label" />
        </xsl:otherwise>
      </xsl:choose>
      <xsl:text>"
</xsl:text>
    </xsl:for-each>
 <xsl:text>@ autoticks
</xsl:text>
<xsl:if test="not( $l)">
<xsl:for-each select="$data">
<xsl:variable name="setnr" select="position()-1"/>
<xsl:variable name="species" select="../../@chemicalSymbol"/>
<xsl:if test="point/bc">
<xsl:for-each select="point">
<xsl:variable name="x" select="@distance"/>
<xsl:variable name="y" select="@eval*27.21138386"/>
<xsl:variable name="n" select="position()"/>
<xsl:for-each select="bc">
<xsl:if test="@character>0.1 and $n mod 2=@l mod 2">
@with string
@    string on
@    string loctype world
@    string g0
@    string <xsl:value-of select="$x"/> , <xsl:value-of select="$y"/>
@    string color <xsl:value-of select="@l"/>
@    string rot 40
@    string font 12
@    string just 0
@    string char size <xsl:value-of select="@character*2"/>
@    string def "<xsl:value-of select="$species"/> l=<xsl:value-of select="@l"/>"
</xsl:if>
</xsl:for-each>
</xsl:for-each>
</xsl:if>
</xsl:for-each>
</xsl:if>
    <xsl:for-each select="$data">
      <xsl:text>@target G0.S</xsl:text>
      <xsl:value-of select="position()" />
  
<xsl:choose>
<xsl:when test="$l">
<xsl:text>
@type xysize
</xsl:text>
</xsl:when>
<xsl:otherwise>
<xsl:text>
@type xy
</xsl:text>
</xsl:otherwise>
</xsl:choose>

      <xsl:for-each select="./point">
        <xsl:value-of select="@distance" />
        <xsl:text>  </xsl:text>
        <xsl:value-of select="@eval*27.21138386" />
        <xsl:if test="$l">
          <xsl:text>  </xsl:text>
        <xsl:value-of select="bc[@l=$l]/@character*2" />
        </xsl:if>
        <xsl:text>
</xsl:text>
      </xsl:for-each>
      <xsl:text>&amp;
</xsl:text>
    </xsl:for-each>
  </xsl:template>
  <xsl:template match="/">
<xsl:if test="/bandstructure/band" >
<xsl:call-template name="bandstograce">
<xsl:with-param name="data" select="bandstructure/band"/>
</xsl:call-template>
</xsl:if> 
<xsl:if test="/bandstructure/species/atom/band" >
<xsl:for-each select="/bandstructure/species[1]/atom[1]/band[1]/point[1]/bc">
<xsl:variable name="l" select="@l"/>
<xsl:for-each select="/bandstructure/species">
<xsl:for-each select="./atom">

<xsl:variable name="filename">
<xsl:value-of select="/bandstructure/title"/>
<xsl:text>_band_</xsl:text><xsl:value-of select="../@name"/>
<xsl:text>_</xsl:text><xsl:value-of select="position()"/>
<xsl:text>_l</xsl:text><xsl:value-of select="$l"/>
<xsl:text>.agr</xsl:text>
</xsl:variable>
<xsl:document href="{$filename}"  method="text">
<xsl:call-template name="bandstograce">
<xsl:with-param name="data" select="./band"/>
<xsl:with-param name="l" select="$l"/>
</xsl:call-template>
</xsl:document>
</xsl:for-each>
</xsl:for-each>
</xsl:for-each>
</xsl:if> 

  </xsl:template>
</xsl:stylesheet>
