<?xml version="1.0" encoding="UTF-8" ?>
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform" 
xmlns:xs="http://www.w3.org/2001/XMLSchema"
 xmlns:math="http://exslt.org/math"
>
  <xsl:output method="xml" indent='yes'/> 
 <xsl:template match="/">
 <report>
  <test>
    <status>
    <xsl:choose>
    <xsl:when test="/info/groundstate/@status='finished'"><xsl:text>passed</xsl:text></xsl:when>
    <xsl:otherwise><xsl:text>failed</xsl:text></xsl:otherwise>
    </xsl:choose>
    </status>
    <name>  libxc  </name>
    <description>passes if groundstate using libxc works</description>
    <directory>test05/ </directory>
  </test>
  
  <test>
    <status>
    <xsl:choose>
    <xsl:when test="math:abs(/info/groundstate/scl/iter[last()]/energies/@totalEnergy+ 242.346369077)&lt;0.0001"><xsl:text>passed</xsl:text></xsl:when>
    <xsl:otherwise><xsl:text>failed</xsl:text></xsl:otherwise>
    </xsl:choose>
    </status>
    <name>  libxc energy  </name>
    <description>passes if groundstate using libxc yields total energy within range +-0.0001
    <xsl:value-of select="/info/groundstate/scl/iter[last()]/energies/@totalEnergy"/>
    </description>
    <directory>test05/ </directory>
  </test>
  
  <!--test>
  <status>
    <xsl:choose>
    <xsl:when test="document('runelectrstr/dos.xml')/dos"><xsl:text>passed</xsl:text></xsl:when>
    <xsl:otherwise><xsl:text>failed</xsl:text></xsl:otherwise>
    </xsl:choose>
    </status>
    <name>  dos works  </name>
    <description>passes if DOS doesn't fail'
 
    </description>
    <directory>test05/runelectrstr</directory>
  </test-->
  
  <!--test>
    <status>
    <xsl:choose>
    <xsl:when test="document('runelectrstr/bandstructure.xml')/bandstructure"><xsl:text>passed</xsl:text></xsl:when>
    <xsl:otherwise><xsl:text>failed</xsl:text></xsl:otherwise>
    </xsl:choose>
    </status>
    <name>  bandstructure  </name>
    <description>passes if bandstructure doesn't fail
 
    </description>
    <directory>test05/runelectrstr</directory>
  </test-->
  
  <!--test>
  <status>
    <xsl:choose>
    <xsl:when test="document('runelectrstr/fermisurface.xml')/fermisurface"><xsl:text>passed</xsl:text></xsl:when>
    <xsl:otherwise><xsl:text>failed</xsl:text></xsl:otherwise>
    </xsl:choose>
    </status>
    <name>  fermisurfaceplot works  </name>
    <description>pass if  fermisurfaceplot doesn't fail'
 
    </description>
    <directory>test05/runelectrstr</directory>
  </test-->
  
</report>
</xsl:template>
</xsl:stylesheet> 
