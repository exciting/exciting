<?xml version="1.0" encoding="UTF-8" ?>
<!-- 
##################################################################
 Use:
 
 xsltproc --path "../../species/" exciting3wienstruct.xsl input.xml
 
 give the path option when species definitions should be loaded from
 that path
###################################################################
 -->
<xsl:stylesheet version="1.0"
  xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
  xmlns:str="http://exslt.org/strings" 
  xmlns:math="http://exslt.org/math">
  <xsl:output method="text" />
  <xsl:template match="/">
<xsl:value-of select="//title"/>
<xsl:text>
P                           </xsl:text>
<xsl:value-of select="format-number(count(//atom),'000000')"/>
<xsl:text>
             RELA
</xsl:text><xsl:value-of select="format-number($a,' 00.000000')"/>
  <xsl:value-of select="format-number($b,' 00.000000 ')"/>
  <xsl:value-of select="format-number($c,' 00.000000 ')"/>
  <xsl:value-of select="format-number($alpha,' 00.000000 ')"/>
  <xsl:value-of select="format-number($beta,' 00.000000 ')"/>
  <xsl:value-of select="format-number($gamma,' 00.000000 ')"/>
  <xsl:text>
</xsl:text>
<xsl:for-each select="//species">
<xsl:variable name="speciesfile">
<xsl:value-of select="@speciesfile"/>
</xsl:variable>
<xsl:for-each select="atom">
<xsl:text>
ATOM </xsl:text>
<xsl:value-of select ="format-number(position(),'000')"/>
<xsl:text>: X=</xsl:text>
<xsl:value-of select ="format-number(str:tokenize(@coord)[1],'0.00000000')"/>
<xsl:text> Y=</xsl:text>
<xsl:value-of select ="format-number(str:tokenize(@coord)[2],'0.00000000')"/>
<xsl:text> Z=</xsl:text>
<xsl:value-of select ="format-number(str:tokenize(@coord)[3],'0.00000000')"/>
<xsl:text>          MULT= 1          ISPLIT= 8
</xsl:text>
<xsl:value-of select="../@chemicalSymbol|document($speciesfile)/speciesdb/species/@chemicalSymbol"/>
<xsl:value-of select="str:padding(10 - string-length(../@chemicalSymbol|document($speciesfile)/speciesdb/species/@chemicalSymbol),' ')"/>
<xsl:text>NPT=</xsl:text>  
<xsl:choose>
<xsl:when test="document($speciesfile)">
<xsl:value-of select="format-number(document($speciesfile)/speciesdb/species/muffinTin/@radialmeshPoints,'00000')"/>
</xsl:when>
<xsl:otherwise>
<xsl:text>  700</xsl:text>
</xsl:otherwise>
</xsl:choose>
<xsl:text>  R0=</xsl:text>
<xsl:choose>
<xsl:when test="document($speciesfile)">
<xsl:value-of select="format-number(document($speciesfile)/speciesdb/species/muffinTin/@rmin,'0.00000000')"/>
</xsl:when>
<xsl:otherwise><xsl:text>0.00005   </xsl:text> </xsl:otherwise>
</xsl:choose>
<xsl:text>   RMT=</xsl:text>
<xsl:value-of select="format-number(../@rmt|document($speciesfile)/speciesdb/species/muffinTin/@radius,'0000.00000')"/>
<xsl:text>     Z:</xsl:text>
<xsl:value-of select="format-number( math:abs(../@atomicNumber|document($speciesfile)/speciesdb/species/@z),'00.00')"/>
<xsl:text>  
LOCAL ROT MATRIX:    1.0000000 0.0000000 0.0000000
                     0.0000000 1.0000000 0.0000000
                     0.0000000 0.0000000 1.0000000
</xsl:text>
</xsl:for-each>
</xsl:for-each>
<xsl:text>
  0      NUMBER OF SYMMETRY OPERATIONS
</xsl:text>
  </xsl:template>
   <xsl:template name="norm">
    <xsl:param name="vectorstring" />
    <xsl:value-of
      select="math:sqrt(
math:power(str:tokenize($vectorstring)[1]*$scale,2)
+math:power(str:tokenize($vectorstring)[2]*$scale,2)
+math:power(str:tokenize($vectorstring)[3]*$scale,2)
)" />
  </xsl:template>
    <xsl:variable name="bohr2angstr" select="0.529177" />
  <xsl:variable name="scale">
    <xsl:choose>
      <xsl:when test="/input/structure/crystal/@scale">
        <xsl:value-of select="/input/structure/crystal/@scale" />
      </xsl:when>
      <xsl:otherwise>
        <xsl:value-of select="1" />
      </xsl:otherwise>
    </xsl:choose>
  </xsl:variable>
  <xsl:variable name="a">
    <xsl:call-template name="norm">
      <xsl:with-param name="vectorstring" select="/input/structure/crystal/basevect[1]" />
    </xsl:call-template>
  </xsl:variable>
  <xsl:variable name="b">
    <xsl:call-template name="norm">
      <xsl:with-param name="vectorstring" select="/input/structure/crystal/basevect[2]" />
    </xsl:call-template>
  </xsl:variable>
  <xsl:variable name="c">
    <xsl:call-template name="norm">
      <xsl:with-param name="vectorstring" select="/input/structure/crystal/basevect[3]" />
    </xsl:call-template>
  </xsl:variable>
  <xsl:template name="angle">
  <xsl:param name="vecs1"/>
  <xsl:param name="vecs2"/>
  <xsl:variable name="vec1" select="str:tokenize($vecs1)"/>
  <xsl:variable name="vec2" select="str:tokenize($vecs2)"/>
  <xsl:variable name="norm1" >
   <xsl:call-template name="norm">
      <xsl:with-param name="vectorstring" select="$vecs1"/>
   </xsl:call-template>   
  </xsl:variable>
  <xsl:variable name="norm2">
    <xsl:call-template name="norm">
      <xsl:with-param name="vectorstring" select="$vecs2"/>
    </xsl:call-template>
  </xsl:variable>
  <xsl:value-of select="57.2957795130823*math:acos(($vec1[1]*$vec2[1]+$vec1[2]*$vec2[2]+$vec1[3]*$vec2[3])div $norm1 div $norm2)"/>
</xsl:template>
<xsl:variable name="alpha">
<xsl:call-template name="angle">
<xsl:with-param name="vecs1" select="/input/structure/crystal/basevect[2]"/>
<xsl:with-param name="vecs2" select="/input/structure/crystal/basevect[3]"/>
</xsl:call-template>
</xsl:variable>
<xsl:variable name="beta">
<xsl:call-template name="angle">
<xsl:with-param name="vecs1" select="/input/structure/crystal/basevect[1]"/>
<xsl:with-param name="vecs2" select="/input/structure/crystal/basevect[3]"/>
</xsl:call-template>
</xsl:variable>
<xsl:variable name="gamma">
<xsl:call-template name="angle">
<xsl:with-param name="vecs1" select="/input/structure/crystal/basevect[1]"/>
<xsl:with-param name="vecs2" select="/input/structure/crystal/basevect[2]"/>
</xsl:call-template>
</xsl:variable>
</xsl:stylesheet>