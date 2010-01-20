<?xml version="1.0" encoding="UTF-8" ?>
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform" xmlns:str="http://exslt.org/strings" xmlns:math="http://exslt.org/math">
  <xsl:output method="xml" />
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
      <xsl:when test="/input/structire/crystal/@scale">
        <xsl:value-of select="/input/structire/crystal/@scale" />
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
  <xsl:template match="/">
   <xsl:element name="system">
   
   </xsl:element>
    <xsl:for-each select="/input/structure/crystal/basevect">
      <xsl:for-each select="str:tokenize(.)">
        <xsl:value-of select="$scale * ./.*$bohr2angstr" />
        <xsl:text>   </xsl:text>
      </xsl:for-each>
      <xsl:text>
   </xsl:text>
    </xsl:for-each>
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
      <xsl:value-of select="str:tokenize(@coord)[1]*$a" />
      <xsl:text>  </xsl:text>
      <xsl:value-of select="str:tokenize(@coord)[2]*$b" />
      <xsl:text>  </xsl:text>
      <xsl:value-of select="str:tokenize(@coord)[3]*$c" />
      <xsl:text>  </xsl:text>
    </xsl:for-each>
    <xsl:text>
   </xsl:text>
  </xsl:template>
</xsl:stylesheet>