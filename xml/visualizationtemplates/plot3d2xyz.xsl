<?xml version="1.0" encoding="UTF-8" ?>
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform" xmlns:str="http://exslt.org/strings" >
  <xsl:output method="text" cdata-section-elements="at" />
  <xsl:variable name="da" select="str:tokenize(/plot3d/grid/axis[1]/@delta)" />
  <xsl:variable name="db" select="str:tokenize(/plot3d/grid/axis[2]/@delta)" />
  <xsl:variable name="dc" select="str:tokenize(/plot3d/grid/axis[3]/@delta)" />
  <xsl:variable name="grid" select="str:tokenize(/plot3d/grid/@gridticks)" />
  
  <xsl:template match="/">
  <xsl:for-each select="/plot3d/function">
  
<xsl:text>#grid </xsl:text><xsl:value-of select="/plot3d/grid/@gridticks" />
 <xsl:text>
</xsl:text>   
    <xsl:for-each select="row">
    <xsl:variable name="indexx" select="@index"/>
      <xsl:for-each select="row">
       <xsl:variable name="indexy" select="@index"/>
        <xsl:for-each select="str:tokenize(.)">
          <xsl:variable name="indexz" select="position()"/>
          <xsl:value-of select="$indexx*$da[1]+ $indexy*$db[1] + position()*$dc[1]" />
          <xsl:text> </xsl:text>
          <xsl:value-of select="$indexx*$da[2]+ $indexy*$db[2]+position()*$dc[2]" />
          <xsl:text> </xsl:text>
          <xsl:value-of select="$indexx*$da[3]+$indexy*$db[3]+position()*$dc[3]" />
          <xsl:text> </xsl:text>
          <xsl:value-of select="." />
          <xsl:text>
</xsl:text>
        </xsl:for-each>
      </xsl:for-each>
      <xsl:text>
</xsl:text>
    </xsl:for-each>
    <xsl:text>
</xsl:text>
  </xsl:for-each>
  </xsl:template>
</xsl:stylesheet>
