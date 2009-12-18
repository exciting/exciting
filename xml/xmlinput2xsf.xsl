<?xml version="1.0" encoding="UTF-8" ?>
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform" xmlns:str="http://exslt.org/strings" xmlns:math="http://exslt.org/math">
  <xsl:output method="text" />
  <!-- usage: xsltproc xmlinput2xcf.xsl input.xml >exciting.in -->
  <xsl:template name="norm">
    <xsl:param name="vectorstring" />
    <xsl:value-of
      select="math:sqrt(
math:power(str:tokenize($vectorstring)[1]*$scale,2)
+math:power(str:tokenize($vectorstring)[2]*$scale,2)
+math:power(str:tokenize($vectorstring)[3]*$scale,2)
)*$bohr2angstr" />
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
   <xsl:variable name="stretch">
    <xsl:choose>
      <xsl:when test="/input/structure/crystal/@stretch">
        <xsl:value-of select="/input/structure/crystal/@stretch" />
      </xsl:when>
      <xsl:otherwise>
        <xsl:value-of select="'1 1 1'" />
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
  <xsl:template match="/">
    <xsl:text>
 INFO
 nunit      1    1    1
 unit   cell
 celltype   convcell
 shape   parapipedal
 END_INFO
 DIM-GROUP
           3           1
 PRIMVEC
 </xsl:text>
    <xsl:for-each select="/input/structure/crystal/basevect">
    <xsl:variable name="basevn" select="position()"/>
      <xsl:for-each select="str:tokenize(.)">
        <xsl:value-of select="$scale * ./.*$bohr2angstr * str:tokenize($stretch)[$basevn]" />
        <xsl:text>   </xsl:text>
      </xsl:for-each>
      <xsl:text>
   </xsl:text>
    </xsl:for-each>
    <xsl:variable name="a1" select=
      "str:tokenize(/input/structure/crystal/basevect[1])[1]*
      $bohr2angstr*$scale*str:tokenize($stretch)[1]" />
    <xsl:variable name="a2" select=
      "str:tokenize(/input/structure/crystal/basevect[1])[2]*
      $bohr2angstr*$scale*str:tokenize($stretch)[1]" />
    <xsl:variable name="a3" select=
      "str:tokenize(/input/structure/crystal/basevect[1])[3]*
      $bohr2angstr*$scale*str:tokenize($stretch)[1]" />
    <xsl:variable name="b1" select=
      "str:tokenize(/input/structure/crystal/basevect[2])[1]*
      $bohr2angstr*$scale*str:tokenize($stretch)[2]" />
    <xsl:variable name="b2" select=
      "str:tokenize(/input/structure/crystal/basevect[2])[2]*
      $bohr2angstr*$scale*str:tokenize($stretch)[2]" />
    <xsl:variable name="b3" select=
      "str:tokenize(/input/structure/crystal/basevect[2])[3]*
      $bohr2angstr*$scale*str:tokenize($stretch)[2]" />
    <xsl:variable name="c1" select=
      "str:tokenize(/input/structure/crystal/basevect[3])[1]*
      $bohr2angstr*$scale*str:tokenize($stretch)[3]" />
    <xsl:variable name="c2" select=
      "str:tokenize(/input/structure/crystal/basevect[3])[2]*
      $bohr2angstr*$scale*str:tokenize($stretch)[3]" />
    <xsl:variable name="c3" select=
      "str:tokenize(/input/structure/crystal/basevect[3])[3]*
      $bohr2angstr*$scale*str:tokenize($stretch)[3]" />
    <xsl:text>  
 PRIMCOORD
    </xsl:text>
    <xsl:value-of select="count(input/structure/species/atom)" />
    <xsl:text> 1</xsl:text>
    <xsl:for-each select="input/structure/species/atom">
      <xsl:text>
</xsl:text>
      <xsl:choose>
        <xsl:when test="../@atomicNumber">
          <xsl:value-of select="../@atomicNumber" />
        </xsl:when>
        <xsl:otherwise>
          <xsl:variable name="speciesfile">
            <xsl:value-of select="/input/structure/@speciespath" />
            <xsl:text>/</xsl:text>
            <xsl:value-of select="../@speciesfile" />
          </xsl:variable>
          <xsl:value-of select="-1*document($speciesfile)/spdb/sp/@z" />
        </xsl:otherwise>
      </xsl:choose>
      <xsl:text>   </xsl:text>
      <xsl:value-of select="str:tokenize(@coord)[1]*$a1 + str:tokenize(@coord)[2]*$b1  + str:tokenize(@coord)[3]*$c1" />
      <xsl:text>  </xsl:text>
      <xsl:value-of select="str:tokenize(@coord)[1]*$a2 + str:tokenize(@coord)[2]*$b2  + str:tokenize(@coord)[3]*$c2" />
      <xsl:text>  </xsl:text>
      <xsl:value-of select="str:tokenize(@coord)[1]*$a3 + str:tokenize(@coord)[2]*$b3  + str:tokenize(@coord)[3]*$c3" />
      <xsl:text>  </xsl:text>
    </xsl:for-each>
    <xsl:text>
   </xsl:text>
  </xsl:template>
</xsl:stylesheet>
