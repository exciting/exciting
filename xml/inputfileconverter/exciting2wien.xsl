<?xml version="1.0" encoding="UTF-8" ?>
<!-- 
##################################################################
 Use:
 
 xsltproc [dash dash]path "../../species/" exciting3wienstruct.xsl input.xml
 
 give the path option when species definitions should be loaded from
 that path
###################################################################
 -->
<xsl:stylesheet version="1.0"
  xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
  xmlns:str="http://exslt.org/strings" 
  xmlns:math="http://exslt.org/math">
  <xsl:output method="text" />
  <xsl:include href="basevec2abc.xsl" />
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
<xsl:value-of select="../@chemicalSymbol|document($speciesfile)/spdb/sp/@chemicalSymbol"/>
<xsl:value-of select="str:padding(10 - string-length(../@chemicalSymbol|document($speciesfile)/spdb/sp/@chemicalSymbol),' ')"/>
<xsl:text>NPT=</xsl:text>  
<xsl:choose>
<xsl:when test="document($speciesfile)">
<xsl:value-of select="format-number(document($speciesfile)/spdb/sp/muffinTin/@radialmeshPoints,'00000')"/>
</xsl:when>
<xsl:otherwise>
<xsl:text>  700</xsl:text>
</xsl:otherwise>
</xsl:choose>
<xsl:text>  R0=</xsl:text>
<xsl:choose>
<xsl:when test="document($speciesfile)">
<xsl:value-of select="format-number(document($speciesfile)/spdb/sp/muffinTin/@rmin,'0.00000000')"/>
</xsl:when>
<xsl:otherwise><xsl:text>0.00005   </xsl:text> </xsl:otherwise>
</xsl:choose>
<xsl:text>   RMT=</xsl:text>
<xsl:value-of select="format-number(../@rmt|document($speciesfile)/spdb/sp/muffinTin/@radius,'0000.00000')"/>
<xsl:text>     Z:</xsl:text>
<xsl:value-of select="format-number( math:abs(../@atomicNumber|document($speciesfile)/spdb/sp/@z),'00.00')"/>
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
</xsl:stylesheet>