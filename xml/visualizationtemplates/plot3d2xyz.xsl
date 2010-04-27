<?xml version="1.0" encoding="UTF-8" ?>
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">
  <xsl:output method="text" cdata-section-elements="at" />
  <xsl:variable name="da" select="str:tokenize(/plot3d/grid/axis[1]/@delta)" />
  <xsl:variable name="db" select="str:tokenize(/plot3d/grid/axis[2]/@delta)" />
  <xsl:variable name="dc" select="str:tokenize(/plot3d/grid/axis[3]/@delta)" />
  <xsl:template match="/">
<xsl:text>
</xsl:text>
grid <xsl:value-of select="/plot3d/grid/@gridticks" />
 <xsl:text>
</xsl:text>   
    <xsl:for-each select="/function[1]/set">
      <xsl:for-each select="set">
        <xsl:for-each select="str:tokenize(.)">
          <xsl:value-of select="../../@index*$da[1]+ ../@index*$db[1]+position()*$dc[1]" />
          <xsl:text> </xsl:text>
          <xsl:value-of select="../../@index*$da[1]+ ../@index*$db[1]+position()*$dc[1]" />
          <xsl:text> </xsl:text>
          <xsl:value-of select="../../@index*$da[1]+ ../@index*$db[1]+position()*$dc[1]" />
          <xsl:text> </xsl:text>
          <xsl:value-of select="." />
          <xsl:for-each select="/plot3d/function[position()!=1]">
            <xsl:text> </xsl:text>
            <xsl:value-of select="str:tokenize(.)[$count]" />
          </xsl:for-each>
           <xsl:text>
</xsl:text>
        </xsl:for-each>
      </xsl:for-each>
    </xsl:for-each>
  </xsl:template>
</xsl:stylesheet>