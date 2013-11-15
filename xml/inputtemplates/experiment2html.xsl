<?xml version="1.0" encoding="UTF-8" ?>
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform" xmlns:exsl="http://exslt.org/common"
 xmlns:str="http://exslt.org/strings">
  <!--
    documentation :
    http://exciting-code.org/expand-all-parameter-permutations
  -->
  <xsl:output method="xml" encoding="UTF-8" indent="yes" />
  
  <xsl:template match="/">
  <table>
  <xsl:for-each select="//set">
 <tr>
 <td>
 <a>
 <xsl:attribute name="href">
 <xsl:text>/exist/rest/db/calculations/</xsl:text>
 <xsl:value-of select="@hash"/>
 <xsl:text>/input.xml?_xsl=http://xml.exciting-code.org/inputfileconverter/xmlinput2html.xsl</xsl:text>
 </xsl:attribute>
 calc
 </a>
 </td>
 <xsl:for-each select="@*">
 <xsl:if test="not (name(.)='hash' or  name(.)='path')">
 <td>
 <xsl:value-of select="name(.)"/><xsl:text>:</xsl:text>  <xsl:value-of select="."/> 
 </td>
 </xsl:if>
 </xsl:for-each>
 </tr>
  </xsl:for-each>
  </table>
  </xsl:template>
  </xsl:stylesheet>